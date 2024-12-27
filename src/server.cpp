#include "server.h"
#include "external_prod.h"
#include "utils.h"
#include <cassert>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <fstream>

// copy the pir_params and set evaluator equal to the context_. 
// client_galois_keys_, client_gsw_keys_, and db_ are not set yet.
PirServer::PirServer(const PirParams &pir_params)
    : pir_params_(pir_params), context_(pir_params.get_seal_params()),
      DBSize_(pir_params.get_DBSize()), evaluator_(context_), dims_(pir_params.get_dims()),
      hashed_key_width_(pir_params_.get_hashed_key_width()) {
  // delete the raw_db_file if it exists
  std::remove(RAW_DB_FILE);

  // allocate enough space for the database, init with std::nullopt
  db_ = std::make_unique<std::optional<seal::Plaintext>[]>(DBSize_);
  fill_inter_res();

}

PirServer::~PirServer() {
  // delete the raw_db_file
  std::remove(RAW_DB_FILE);
}


Entry generate_entry(const uint64_t id, const size_t entry_size, std::ifstream &random_file) {
  Entry entry(entry_size);
  // use the random file in ther kernal
  entry.resize(entry_size);
  random_file.read(reinterpret_cast<char *>(entry.data()), entry_size);
  return entry;
}


Entry generate_entry_with_key(uint64_t key_id, size_t entry_size, size_t hashed_key_width) {
  if (entry_size < hashed_key_width) {
    throw std::invalid_argument("Entry size is too small for the hashed key width");
  }

  Entry entry;
  entry.reserve(entry_size);
  std::mt19937_64 rng(key_id);
  // generate the entire entry using random numbers for simplicity.
  for (int i = 0; i < entry_size; i++) {
    entry.push_back(rng() % 256);
  }
  return entry;
}


// Fills the database with random data
void PirServer::gen_data() {

  // init the database with std::nullopt
  db_.reset(new std::optional<seal::Plaintext>[DBSize_]);

  auto fst_dim_sz = dims_[0];
  const size_t other_dim_sz = DBSize_ / fst_dim_sz;
  std::ifstream random_file("/dev/urandom", std::ios::binary);
  if (!random_file.is_open()) {
    throw std::invalid_argument("Unable to open /dev/urandom");
  }

  for (size_t j = 0; j < other_dim_sz; ++j) {
    std::vector<Entry> one_chunk(fst_dim_sz, Entry(pir_params_.get_entry_size()));
    for (size_t k = 0; k < fst_dim_sz; ++k) {
      one_chunk[k] = generate_entry(other_dim_sz * k + j, pir_params_.get_entry_size(), random_file);
    }
    write_one_chunk(one_chunk);
    push_database_chunk(one_chunk, j);
    print_progress(j+1, other_dim_sz);
  }
  random_file.close();
  // transform the ntt_db_ from coefficient form to ntt form. db_ is not transformed.
  BENCH_PRINT("Transforming the database to NTT form...");
  preprocess_ntt();
}

std::vector<CuckooSeeds> PirServer::gen_keyword_data(size_t max_iter, uint64_t keyword_seed) {
  // Generate random keywords for the database.
  std::vector<Key> keywords;
  size_t key_num = pir_params_.get_num_entries() / pir_params_.get_blowup_factor(); // TODO: put this as pir params
  keywords.reserve(key_num);
  // We randomly generate a bunch of keywords. Then, we treat each keyword in the key-value pair as a seed.
  // In this the current method, all keyword is generated using the same keyword_seed given by client.
  std::mt19937_64 key_rng(keyword_seed); 
  for (size_t i = 0; i < key_num; ++i) {
    keywords.push_back(key_rng()); 
  }

  DEBUG_PRINT(keywords.size() << " keywords generated");
  // check if the keywords are all unique: 
  std::unordered_set<Key> unique_keywords(keywords.begin(), keywords.end());
  if (unique_keywords.size() != keywords.size()) {
    std::cerr << "Keywords are not unique" << std::endl;
    return {{}, {}};
  } else {
    DEBUG_PRINT("Keywords are unique");
  }

  // Insert data into the database using cuckoo hashing
  std::vector<CuckooSeeds> used_seeds;
  // std::mt19937_64 hash_rng;
  for (size_t i = 0; i < max_iter; i++) {
    uint64_t seed1 = key_rng();
    uint64_t seed2 = key_rng();
    used_seeds.push_back({seed1, seed2});
    std::vector<Key> cuckoo_hash_table = cuckoo_insert(seed1, seed2, 100, keywords, pir_params_.get_blowup_factor());
    // now we have a successful insertion. We create the database using the keywords we have and their corresponding values.
    if (cuckoo_hash_table.size() > 0) {
      std::vector<Entry> data(key_num); 
      
      // we insert key-value pair one by one. Generating the entries on the fly.
      size_t entry_size = pir_params_.get_entry_size();
      size_t hashed_key_width = pir_params_.get_hashed_key_width();
      for (size_t j = 0; j < pir_params_.get_num_entries(); ++j) {
        // Keyword(string) -> hash to fixed size bit string
        Entry entry = generate_entry_with_key(keywords[j], entry_size, hashed_key_width);
        size_t index1 = std::hash<Key>{}(keywords[j] ^ seed1) % cuckoo_hash_table.size();
        size_t index2 = std::hash<Key>{}(keywords[j] ^ seed2) % cuckoo_hash_table.size(); 
        if (cuckoo_hash_table[index1] == keywords[j]) {
          data[index1] = entry;
        } else {
          data[index2] = entry;
        }
      }

      // set the database and return the used seeds and the database to the client. Data is returned for debugging purposes.
      // set_database(data);  // TODO: fix this set_database using new method.
      return used_seeds;
    }
  }
  std::cerr << "Failed to insert data into cuckoo hash table" << std::endl;
  // resize the cuckoo_hash_table to 0
  return used_seeds;
}


// Computes a dot product between the fst_dim_query and the database for the
// first dimension with a delayed modulus optimization. fst_dim_query should
// be transformed to ntt.
std::vector<seal::Ciphertext>
PirServer::evaluate_first_dim(std::vector<seal::Ciphertext> &fst_dim_query) {
  const size_t fst_dim_sz = dims_[0];  // number of entries in the first dimension
  const size_t other_dim_sz = DBSize_ / fst_dim_sz;  // number of entries in the other dimensions
  const auto seal_params = context_.get_context_data(fst_dim_query[0].parms_id())->parms();
  const auto coeff_modulus = seal_params.coeff_modulus();
  const size_t poly_degree = seal_params.poly_modulus_degree();
  const size_t coeff_mod_count = coeff_modulus.size();
  const size_t coeff_val_cnt = poly_degree * coeff_mod_count; // polydegree * RNS moduli count
  constexpr size_t num_poly = 2;  // ciphertext has two polynomials
  
  // transform the selection vector to ntt form
  for (size_t i = 0; i < fst_dim_query.size(); i++) {
    evaluator_.transform_to_ntt_inplace(fst_dim_query[i]);
  }

  // fill the intermediate result with zeros
  std::fill(inter_res.begin(), inter_res.end(), 0);

  /*
  I imagine DB as a (other_dim_sz * fst_dim_sz) matrix, each column is
  other_dim_sz many consecutive entries in the database. We are going to
  multiply the selection_vector with the DB. Then only one row of the result
  is going to be added to the result vector.
  The high level is summing C_{BFV_k} * DB_{N_1 * j + k}
  */

  const size_t one_ct_size = num_poly * coeff_val_cnt;


  // const size_t tile_size = 8;

  // for (size_t k_base = 0; k_base < fst_dim_sz; k_base += tile_size) {
  //   for (size_t j = 0; j < other_dim_sz; ++j) {
  //     for (size_t k = k_base; k < std::min(k_base + tile_size, fst_dim_sz); k++) {
  //       for (size_t poly_id = 0; poly_id < num_poly; poly_id++) {
  //         // j are already filled with the results
  //         auto inter_shift = j * one_ct_size + poly_id * coeff_val_cnt;
  //         utils::multiply_poly_acum(fst_dim_query[k].data(poly_id),
  //                                   (*db_[fst_dim_sz * j + k]).data(),
  //                                   coeff_val_cnt, inter_res.data() + inter_shift);
  //       }
  //     }
  //   }
  // }


  for (size_t j = 0; j < other_dim_sz; ++j) {
    for (size_t k = 0; k < fst_dim_sz; k++) {
      for (size_t poly_id = 0; poly_id < num_poly; poly_id++) {
        // j are already filled with the results
        auto inter_shift = j * one_ct_size + poly_id * coeff_val_cnt;
        utils::multiply_poly_acum(fst_dim_query[k].data(poly_id),
                                  (*db_[fst_dim_sz * j + k]).data(),
                                  coeff_val_cnt, inter_res.data() + inter_shift);
      }
    }
  }

  // ========== transform the intermediate to coefficient form. Delay the modulus operation ==========
  std::vector<seal::Ciphertext> result; // output vector
  result.reserve(other_dim_sz);
  seal::Ciphertext ct_acc;

  for (size_t j = 0; j < other_dim_sz; ++j) {
    ct_acc = fst_dim_query[fst_dim_sz - 1]; // just a quick way to construct a new ciphertext. Will overwrite data in it.
    for (size_t poly_id = 0; poly_id < num_poly; poly_id++) {   // Each ciphertext has two polynomials
      auto mod_acc_ptr = ct_acc.data(poly_id); // pointer to store the modulus of accumulated value
      auto inter_shift = j * one_ct_size + poly_id * coeff_val_cnt;
      auto buff_ptr = inter_res.data() + inter_shift;
      
      for (int mod_id = 0; mod_id < coeff_mod_count; mod_id++) {  // RNS has two moduli
        // Now we calculate the modulus for the accumulated value.
        auto rns_padding = (mod_id * poly_degree);
        for (int coeff_id = 0; coeff_id < poly_degree; coeff_id++) {  // for each coefficient, we mod it with RNS modulus
          // the following is equivalent to: mod_acc_ptr[coeff_id + rns_padding] = buff_ptr[coeff_id + rns_padding] % coeff_modulus[mod_id]
          auto x = buff_ptr[coeff_id + rns_padding];
          uint64_t raw[2] = {static_cast<uint64_t>(x), static_cast<uint64_t>(x >> 64)};
          mod_acc_ptr[coeff_id + rns_padding] = util::barrett_reduce_128(raw, coeff_modulus[mod_id]);
        }
      }
    }
    evaluator_.transform_from_ntt_inplace(ct_acc);  // transform the result back to coefficient form
    result.push_back(ct_acc);
  }

  return result;
}


// NO, THIS IS TOO SLOW.
std::vector<seal::Ciphertext>
PirServer::evaluate_first_dim_direct_mod(std::vector<seal::Ciphertext> &fst_dim_query) {
  const size_t fst_dim_sz = dims_[0];  // number of entries in the first dimension
  const size_t other_dim_sz = DBSize_ / fst_dim_sz;  // number of entries in the other dimensions
  
  // transform the selection vector to ntt form
  for (size_t i = 0; i < fst_dim_query.size(); i++) {
    evaluator_.transform_to_ntt_inplace(fst_dim_query[i]);
  }

  std::vector<seal::Ciphertext> result(other_dim_sz); // output vector
  seal::Ciphertext temp;

  // I imagine DB as a (other_dim_sz * fst_dim_sz) matrix, each column is
  // other_dim_sz many consecutive entries in the database. We are going to
  // multiply the selection_vector with the DB. Then only one row of the result
  // is going to be added to the result vector.
  for (size_t j = 0; j < other_dim_sz; ++j) {
    // summing C_{BFV_k} * DB_{N_1 * j + k}
    evaluator_.multiply_plain(fst_dim_query[0], db_[fst_dim_sz * j].value(), result[j]);
    for (size_t k = 1; k < fst_dim_sz; k++) {
      evaluator_.multiply_plain(fst_dim_query[k], db_[fst_dim_sz * j + k].value(), temp);
      evaluator_.add_inplace(result[j], temp);
    }
  }
  // transform the result to coefficient form
  for (size_t i = 0; i < result.size(); i++) {
    evaluator_.transform_from_ntt_inplace(result[i]);
  }

  return result;
}


void PirServer::evaluate_gsw_product(std::vector<seal::Ciphertext> &result,
                                                              GSWCiphertext &selection_cipher) {
  
  /**
   * Note that we only have a single GSWCiphertext for this selection.
   * Here is the logic:
   * We want to select the correct half of the "result" vector. 
   * Suppose result = [x || y], where x and y are of the same size(block_size).
   * If we have RGSW(0), then we want to set result = x, 
   * If we have RGSW(1), then we want to set result = y.
   * The simple formula is: 
   * result = RGSW(b) * (y - x) + x, where "*" is the external product, "+" and "-" are homomorphic operations.
   */
  auto block_size = result.size() / 2;
  for (int i = 0; i < block_size; i++) {
    auto &x = result[i];
    auto &y = result[i + block_size];
    evaluator_.sub_inplace(y, x);  // y - x
    data_gsw.external_product(selection_cipher, y, y);  // b * (y - x)
    data_gsw.ciphertext_inverse_ntt(y);
    evaluator_.add_inplace(result[i], y);  // x + b * (y - x)
  }
  result.resize(block_size);
}

// This function is using the algorithm 5 in Constant-weight PIR: Single-round Keyword PIR via Constant-weight Equality Operators.
// https://www.usenix.org/conference/usenixsecurity22/presentation/mahdavi. Basically, the algorithm 3 in Onion-Ring ORAM has some typos.
// And we can save one Subs(c_b, k) operation in the algorithm 3. The notations of this function follows the constant-weight PIR paper.
std::vector<seal::Ciphertext> PirServer::expand_query(uint32_t client_id,
                                                      seal::Ciphertext &ciphertext) const {
  seal::EncryptionParameters params = pir_params_.get_seal_params();
  int poly_degree = params.poly_modulus_degree();   // n in paper. The degree of the polynomial

  // This aligns with the number of coeff used by the client.
  int num_cts = dims_[0] + pir_params_.get_l() * (dims_.size() - 1);  

  int log2N = 0; // log2(num_cts) rounds up. This is the same as padding num_cts to the next power of 2 then taking the log2.
  while ((1 << log2N) < num_cts) {
    log2N++;
  }

  // The access pattern to this array looks like this: https://raw.githubusercontent.com/helloboyxxx/images-for-notes/master/uPic/expansion.png
  // It helps me to understand this recursion :)
  std::vector<Ciphertext> cts((size_t)pow(2, log2N));
  cts[0] = ciphertext;   // c_0 = c in paper

  const auto& client_galois_key = client_galois_keys_.at(client_id); // used for substitution

  for (size_t a = 0; a < log2N; a++) {

    int expansion_const = pow(2, a);

    for (size_t b = 0; b < expansion_const; b++) {
      Ciphertext cipher0 = cts[b];   // c_b in paper
      evaluator_.apply_galois_inplace(cipher0, poly_degree / expansion_const + 1,
                                      client_galois_key); // Subs(c_b, n/k + 1)
      Ciphertext temp;
      evaluator_.sub(cts[b], cipher0, temp);
      utils::shift_polynomial(params, temp,
                              cts[b + expansion_const],
                              -expansion_const);
      evaluator_.add_inplace(cts[b], cipher0);
    }
  }

  return cts;
}

void PirServer::set_client_galois_key(const uint32_t client_id, std::stringstream &galois_stream) {
  seal::GaloisKeys client_key;
  client_key.load(context_, galois_stream);
  client_galois_keys_[client_id] = client_key;
}

void PirServer::set_client_gsw_key(const uint32_t client_id, std::stringstream &gsw_stream) {
  std::vector<seal::Ciphertext> temp_gsw;
  // load 2l ciphertexts from the stream
  for (size_t i = 0; i < 2 * pir_params_.get_l_key(); i++) {
    seal::Ciphertext row;
    row.load(context_, gsw_stream);
    temp_gsw.push_back(row);
  }
  GSWCiphertext gsw_key;

  key_gsw.sealGSWVecToGSW(gsw_key, temp_gsw);
  key_gsw.gsw_ntt_negacyclic_harvey(gsw_key); // transform the GSW ciphertext to NTT form

  client_gsw_keys_[client_id] = gsw_key;
}



Entry PirServer::direct_get_entry(const uint64_t abstract_idx) {
  // Calculate the actual index based on the abstract index
  auto actual_idx = entry_idx_to_actual(abstract_idx, dims_[0], DBSize_);

  // read the entry from raw_db_file
  std::ifstream in_file(RAW_DB_FILE, std::ios::binary);
  if (!in_file.is_open()) {
    throw std::invalid_argument("Unable to open file for reading");
  }

  // Seek to the correct position to the plaintext in the file
  in_file.seekg(actual_idx * pir_params_.get_entry_size());

  // Read the entry from the file
  Entry entry(pir_params_.get_entry_size());
  in_file.read(reinterpret_cast<char *>(entry.data()), pir_params_.get_entry_size());
  in_file.close();

  return entry;
}



std::vector<Key> cuckoo_insert(uint64_t seed1, uint64_t seed2, size_t swap_limit,
                                 std::vector<Key> &keywords, float blowup_factor) {
  std::vector<uint64_t> two_tables(keywords.size() * blowup_factor, 0); // cuckoo hash table for keywords
  size_t half_size = two_tables.size();

  // loop and insert each key-value pair into the cuckoo hash table.
  std::hash<Key> hasher;
  for (size_t i = 0; i < keywords.size(); ++i) {
    Key holding = keywords[i]; // initialy holding is the keyword. Also used for swapping.
    // insert the holding value
    bool inserted = false;
    for (size_t j = 0; j < swap_limit; ++j) {
      // hash the holding keyword to indices in the table
      size_t index1 = std::hash<Key>{}(holding ^ seed1) % half_size;
      
      if (two_tables[index1] == 0) {
        two_tables[index1] = holding;
        inserted = true;
        break;
      }
      std::swap(holding, two_tables[index1]); // swap the holding value with the value in the table
      
      // hash the holding keyword to another index in the table
      size_t index2 = (std::hash<Key>{}(holding ^ seed2) % half_size);
      assert(index1 + half_size != index2); // two hash functions should not hash to the same "index".
      if (two_tables[index2] == 0) {
        two_tables[index2] = holding;
        inserted = true;
        break;
      }
      std::swap(holding, two_tables[index2]); // swap the holding value with the value in the table
    }
    if (inserted == false) {
      DEBUG_PRINT("num_inserted: " << i);
      // print the two indices that are causing the problem.
      size_t holding_index1 = std::hash<Key>{}(holding ^ seed1) % half_size;
      size_t holding_index2 = (std::hash<Key>{}(holding ^ seed2) % half_size);

      Key first = two_tables[holding_index1];
      Key second = two_tables[holding_index2];
      DEBUG_PRINT("index1: " << holding_index1 << " index2: " << holding_index2);
      DEBUG_PRINT("first: " << first << " second: " << second << " holding: " << holding);
      
      // the two hashed indices for first
      size_t first_index1 = std::hash<Key>{}(first ^ seed1) % half_size;
      size_t first_index2 = (std::hash<Key>{}(first ^ seed2) % half_size);
      DEBUG_PRINT("first_index1: " << first_index1 << " first_index2: " << first_index2);

      // the two hashed indices for second
      size_t second_index1 = std::hash<Key>{}(second ^ seed1) % half_size;
      size_t second_index2 = (std::hash<Key>{}(second ^ seed2) % half_size);
      DEBUG_PRINT("second_index1: " << second_index1 << " second_index2: " << second_index2 << "\n");


      return {};  // return an empty vector if the insertion is not successful.
    }
  }
  return two_tables;
}

std::vector<seal::Ciphertext> PirServer::make_query(const uint32_t client_id, PirQuery &&query) {

  // ========================== Expansion & conversion ==========================
  // Query expansion
  auto expand_start = CURR_TIME;
  std::vector<seal::Ciphertext> query_vector = expand_query(client_id, query);
  auto expand_end = CURR_TIME;

  // Reconstruct RGSW queries
  auto convert_start = CURR_TIME;
  auto l = pir_params_.get_l();
  std::vector<GSWCiphertext> gsw_vec(dims_.size() - 1); // GSW ciphertexts
  for (int i = 1; i < dims_.size(); i++) {
    std::vector<seal::Ciphertext> lwe_vector; // BFV ciphertext, size l * 2. This vector will be reconstructed as a single RGSW ciphertext.
    for (int k = 0; k < l; k++) {
      auto ptr = dims_[0] + (i - 1) * l + k;
      lwe_vector.push_back(query_vector[ptr]);
    }
    // Converting the BFV ciphertext to GSW ciphertext
    key_gsw.query_to_gsw(lwe_vector, client_gsw_keys_[client_id], gsw_vec[i - 1]);
  }
  auto convert_end = CURR_TIME;

  // ========================== Evaluations ==========================
  // Evaluate the first dimension
  auto first_dim_start = CURR_TIME;
  std::vector<seal::Ciphertext> result = evaluate_first_dim(query_vector);
  auto first_dim_end = CURR_TIME;

  // Evaluate the other dimensions
  auto other_dim_start = CURR_TIME;
  for (int i = 1; i < dims_.size(); i++) {
    evaluate_gsw_product(result, gsw_vec[i - 1]);
  }
  auto other_dim_end = CURR_TIME;


  // ========================== Post-processing ==========================
  // modulus switching so to reduce the response size by half
  if(pir_params_.get_seal_params().coeff_modulus().size() > 2) {
    DEBUG_PRINT("Modulus switching...");
    evaluator_.mod_switch_to_next_inplace(result[0]); // result.size() == 1.
  }

  // ========================== Timing ==========================
  BENCH_PRINT("\t\tExpand time:\t" << TIME_DIFF(expand_start, expand_end) << "ms");
  BENCH_PRINT("\t\tConvert time:\t" << TIME_DIFF(convert_start, convert_end) << "ms");
  BENCH_PRINT("\t\tFirst dim time:\t" << TIME_DIFF(first_dim_start, first_dim_end) << "ms");
  BENCH_PRINT("\t\tOther dim time:\t" << TIME_DIFF(other_dim_start, other_dim_end) << "ms");

  return result;
}


std::vector<seal::Ciphertext> PirServer::make_seeded_query(const uint32_t client_id, std::stringstream &data_stream) {
  // Deserialize the query
  PirQuery query;
  query.load(context_, data_stream);
  return make_query(client_id, std::move(query));
}


void PirServer::push_database_chunk(std::vector<Entry> &chunk_entry, const size_t chunk_idx) {

  // Flattens data into vector of u8s and pads each entry with 0s to entry_size number of bytes.
  // This is actually resizing from entry.size() to pir_params_.get_entry_size()
  // This is redundent if the given entries uses the same pir parameters.
  for (Entry &entry : chunk_entry) {
    if (entry.size() != 0 && entry.size() <= pir_params_.get_entry_size()) {
      entry.resize(pir_params_.get_entry_size(), 0);
    }

    if (entry.size() > pir_params_.get_entry_size()) {
        std::invalid_argument("Entry size is too large");
    }
  }

  const auto bits_per_coeff = pir_params_.get_num_bits_per_coeff();
  const auto num_coeffs = pir_params_.get_seal_params().poly_modulus_degree();
  const auto num_bits_per_plaintext = pir_params_.get_num_bits_per_plaintext();
  const auto num_entries_per_plaintext = pir_params_.get_num_entries_per_plaintext();
  const auto num_plaintexts = chunk_entry.size() / num_entries_per_plaintext;  // number of plaintexts in the new chunk
  const uint128_t coeff_mask = (uint128_t(1) << (bits_per_coeff)) - 1;  // bits_per_coeff many 1s
  
  const auto fst_dim_sz = dims_[0];
  const auto chunk_offset = fst_dim_sz * chunk_idx;  // offset for the current chunk

  // Now we handle plaintexts one by one.
  for (int i = 0; i < num_plaintexts; i++) {
    seal::Plaintext plaintext(num_coeffs);

    // Loop through the entries that corresponds to the current plaintext. 
    // Then calculate the total size (in bytes) of this plaintext.
    // NOTE: it is possible that some entry is empty, which has size 0.
    int additive_sum_size = 0;
    for (int j = num_entries_per_plaintext * i;
         j < std::min(num_entries_per_plaintext * (i + 1), chunk_entry.size()); j++) {
      additive_sum_size += chunk_entry[j].size();
    }

    if (additive_sum_size == 0) {
      continue; // leave std::nullopt in the chunk if the plaintext is empty.
    }

    int index = 0;  // index for the current coefficient to be filled
    uint128_t data_buffer = 0;
    size_t data_offset = 0;
    // For each entry in the current plaintext
    for (int j = num_entries_per_plaintext * i;
         j < std::min(num_entries_per_plaintext * (i + 1), chunk_entry.size()); j++) {
      // For each byte in this entry
      for (int k = 0; k < pir_params_.get_entry_size(); k++) {
        // data_buffer temporarily stores the data from entry bytes
        data_buffer += uint128_t(chunk_entry[j][k]) << data_offset;
        data_offset += 8;
        // When we have enough data to fill a coefficient
        // We will one by one fill the coefficients with the data_buffer.
        while (data_offset >= bits_per_coeff) {
          plaintext[index] = data_buffer & coeff_mask;
          index++;
          data_buffer >>= bits_per_coeff;
          data_offset -= bits_per_coeff;
        }
      }
    }
    // add remaining data to a new coefficient
    if (data_offset > 0) {
      plaintext[index] = data_buffer & coeff_mask;
      index++;
    }
    db_[i + chunk_offset] = std::move(plaintext);
  }
}

void PirServer::preprocess_ntt() {
  // tutorial on Number Theoretic Transform (NTT): https://youtu.be/Pct3rS4Y0IA?si=25VrCwBJuBjtHqoN
  for (size_t i = 0; i < DBSize_; ++i) {
    if (db_[i].has_value()) {
      seal::Plaintext &pt = db_[i].value();
      evaluator_.transform_to_ntt_inplace(pt, context_.first_parms_id());
    }
  }
}


void PirServer::fill_inter_res() {
  // We need to store 1/dim[0] many ciphertexts in the intermediate result.
  // However, in the first dimension, we want to store them in uint128_t.
  // So, we need to calculate the number of uint128_t we need to store.

  auto seal_params = pir_params_.get_seal_params();
  auto poly_degree = seal_params.poly_modulus_degree();
  // number of rns modulus
  auto num_mods = seal_params.coeff_modulus().size();
  // number of coefficients in a ciphertext
  auto coeff_count = poly_degree * num_mods * 2;  // 2 for two polynomials
  auto poly_count = DBSize_ / dims_[0];
  // number of uint128_t we need to store in the intermediate result
  auto elem_cnt = poly_count * coeff_count;
  // allocate memory for the intermediate result
  inter_res.resize(elem_cnt);
}

void PirServer::write_one_chunk(std::vector<Entry> &data) {
  // write the database to a binary file in CACHE_DIR
  std::string filename = std::string(RAW_DB_FILE);
  std::ofstream out_file(filename, std::ios::binary | std::ios::app); // append to the file
  if (out_file.is_open()) {
    for (auto &entry : data) {
      out_file.write(reinterpret_cast<const char *>(entry.data()), entry.size());
    }
    out_file.close();
  } else {
    std::cerr << "Unable to open file for writing" << std::endl;
  }
}


std::vector<uint64_t> PirServer::get_dims() const {
  return dims_;
}