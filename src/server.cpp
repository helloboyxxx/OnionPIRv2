#include "server.h"
#include "external_prod.h"
#include "utils.h"
#include <cassert>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <fstream>


// copy the pir_params and set evaluator equal to the context_. 
// client_galois_keys_, client_gsw_keys_, and db_ are not set yet.
PirServer::PirServer(const PirParams &pir_params)
    : pir_params_(pir_params), context_(pir_params.get_seal_params()),
      num_pt_(pir_params.get_num_pt()), evaluator_(context_), dims_(pir_params.get_dims()) {
  // delete the raw_db_file if it exists
  std::remove(RAW_DB_FILE);

  // allocate enough space for the database, init with std::nullopt
  db_ = std::make_unique<std::optional<seal::Plaintext>[]>(num_pt_);
  fill_inter_res();
}

PirServer::~PirServer() {
  // delete the raw_db_file
  std::remove(RAW_DB_FILE);
}

// Fills the database with random data
void PirServer::gen_data() {

  std::ifstream random_file("/dev/urandom", std::ios::binary);
  if (!random_file.is_open()) {
    throw std::invalid_argument("Unable to open /dev/urandom");
  }

  // init the database with std::nullopt
  db_.reset(new std::optional<seal::Plaintext>[num_pt_]);
  const size_t fst_dim_sz = pir_params_.get_fst_dim_sz();
  const size_t other_dim_sz = pir_params_.get_other_dim_sz();
  const size_t num_en_per_pt = pir_params_.get_num_entries_per_plaintext();
  const size_t entry_size = pir_params_.get_entry_size();

  for (size_t row = 0; row < other_dim_sz; ++row) {
    std::vector<Entry> one_chunk(fst_dim_sz * num_en_per_pt, Entry(entry_size));
    for (size_t col = 0; col < fst_dim_sz; ++col) {
      size_t poly_id = row * fst_dim_sz + col;
      for (size_t local_id = 0; local_id < num_en_per_pt; ++local_id) {
        size_t entry_id = poly_id * num_en_per_pt + local_id;
        one_chunk[col * num_en_per_pt + local_id] = generate_entry(entry_id, entry_size, random_file);
      }
    }
    write_one_chunk(one_chunk);
    push_database_chunk(one_chunk, row);
    print_progress(row+1, other_dim_sz);
  }
  random_file.close();
  // transform the ntt_db_ from coefficient form to ntt form. db_ is not transformed.
  BENCH_PRINT("\nTransforming the database to NTT form...");
  preprocess_ntt();
}


// Computes a dot product between the fst_dim_query and the database for the
// first dimension with a delayed modulus optimization. fst_dim_query should
// be transformed to ntt.
std::vector<seal::Ciphertext>
PirServer::evaluate_first_dim(std::vector<seal::Ciphertext> &fst_dim_query) {
  const size_t fst_dim_sz = pir_params_.get_fst_dim_sz();  // number of plaintexts in the first dimension
  const size_t other_dim_sz = pir_params_.get_other_dim_sz();  // number of plaintexts in the other dimensions
  const auto seal_params = context_.get_context_data(fst_dim_query[0].parms_id())->parms();
  const auto coeff_modulus = seal_params.coeff_modulus();
  const size_t rns_mod_cnt = coeff_modulus.size();
  const size_t coeff_val_cnt = DatabaseConstants::PolyDegree * rns_mod_cnt; // polydegree * RNS moduli count
  const size_t one_ct_sz = 2 * coeff_val_cnt; // Ciphertext has two polynomials

  // transform the selection vector to ntt form
  for (size_t i = 0; i < fst_dim_query.size(); i++) {
    evaluator_.transform_to_ntt_inplace(fst_dim_query[i]);
  }

  // fill the intermediate result with zeros
  std::fill(inter_res.begin(), inter_res.end(), 0);

  /*
  I imagine DB as a (other_dim_sz * fst_dim_sz) matrix, where each element is a
  vector. The goal is to calculate DB * fst_dim_query, but each multiplication is
  a vector-vector elementwise multiplication.
  TODO: I believe this can be optimized. Currently, the throughput is same as
  the vec-vec elementwise multiplication. However, ideally, the first dimension
  is equivalent to poly_degree many mat-vec multiplications.
  Well, obviously, I need help.
  */
  TIME_START(CORE_TIME);

  size_t row, col, elem_id;
  uint64_t *db_ptr;
  uint64_t *q0, *q1;
  uint128_t *inter_ptr0, *inter_ptr1;
  for (row = 0; row < other_dim_sz; ++row) {
    for (col = 0; col < fst_dim_sz; ++col) {
      db_ptr = db_[row * fst_dim_sz + col].value().data();
      q0 = fst_dim_query[col].data(0); // query polynomial 0
      q1 = fst_dim_query[col].data(1); // query polynomial 1
      inter_ptr0 = inter_res.data() + row * one_ct_sz; // intermediate result ptr for polynomial 0
      inter_ptr1 = inter_res.data() + row * one_ct_sz + coeff_val_cnt;
      #pragma GCC unroll 32
      for (elem_id = 0; elem_id < coeff_val_cnt; ++elem_id) {
        inter_ptr0[elem_id] += static_cast<uint128_t>(q0[elem_id]) * db_ptr[elem_id];
        inter_ptr1[elem_id] += static_cast<uint128_t>(q1[elem_id]) * db_ptr[elem_id];
      }
    }
  }

  TIME_END(CORE_TIME);

  // ========== transform the intermediate to coefficient form. Delay the modulus operation ==========
  TIME_START(FST_DELEY_MOD_TIME);
  std::vector<seal::Ciphertext> result; // output vector
  result.reserve(other_dim_sz);
  seal::Ciphertext ct_acc;
  for (size_t j = 0; j < other_dim_sz; ++j) {
    ct_acc = fst_dim_query[fst_dim_sz - 1]; // just a quick way to construct a new ciphertext. Will overwrite data in it.
    for (size_t poly_id = 0; poly_id < 2; poly_id++) {   // Each ciphertext has two polynomials
      auto mod_acc_ptr = ct_acc.data(poly_id); // pointer to store the modulus of accumulated value
      auto inter_shift = j * one_ct_sz + poly_id * coeff_val_cnt;
      auto buff_ptr = inter_res.data() + inter_shift;
      
      for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {  // RNS has two moduli
        // Now we calculate the modulus for the accumulated value.
        auto rns_padding = (mod_id * DatabaseConstants::PolyDegree);
        for (size_t coeff_id = 0; coeff_id < DatabaseConstants::PolyDegree; coeff_id++) {  // for each coefficient, we mod it with RNS modulus
          // the following is equivalent to: mod_acc_ptr[coeff_id + rns_padding] = buff_ptr[coeff_id + rns_padding] % coeff_modulus[mod_id]
          auto x = buff_ptr[coeff_id + rns_padding];
          uint64_t raw[2] = {static_cast<uint64_t>(x), static_cast<uint64_t>(x >> 64)};
          mod_acc_ptr[coeff_id + rns_padding] = util::barrett_reduce_128(raw, coeff_modulus[mod_id]);
        }
      }
    }
    TIME_END(FST_DELEY_MOD_TIME);
    TIME_START(FST_NTT_TIME);
    evaluator_.transform_from_ntt_inplace(ct_acc);  // transform the result back to coefficient form
    TIME_END(FST_NTT_TIME);
    result.push_back(ct_acc);
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
  for (size_t i = 0; i < block_size; i++) {
    auto &x = result[i];
    auto &y = result[i + block_size];
    evaluator_.sub_inplace(y, x);  // y - x
    data_gsw.external_product(selection_cipher, y, y);  // b * (y - x)
    y.is_ntt_form() = true;
    evaluator_.transform_from_ntt_inplace(y);
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
  // This aligns with the number of coeff used by the client.
  size_t num_cts = dims_[0] + pir_params_.get_l() * (dims_.size() - 1);  

  size_t log2N = 0; // log2(num_cts) rounds up. This is the same as padding num_cts to the next power of 2 then taking the log2.
  while ((1 << log2N) < num_cts) {
    log2N++;
  }

  // The access pattern to this array looks like this: https://raw.githubusercontent.com/helloboyxxx/images-for-notes/master/uPic/expansion.png
  // It helps me to understand this recursion :)
  std::vector<Ciphertext> cts((size_t)pow(2, log2N));
  cts[0] = ciphertext;   // c_0 = c in paper

  const auto& client_galois_key = client_galois_keys_.at(client_id); // used for substitution

  for (size_t a = 0; a < log2N; a++) {

    size_t expansion_const = pow(2, a);

    for (size_t b = 0; b < expansion_const; b++) {
      Ciphertext cipher0 = cts[b];   // c_b in paper
      evaluator_.apply_galois_inplace(cipher0, DatabaseConstants::PolyDegree / expansion_const + 1,
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


Entry PirServer::direct_get_entry(const uint64_t entry_idx) {
  // read the entry from raw_db_file
  std::ifstream in_file(RAW_DB_FILE, std::ios::binary);
  if (!in_file.is_open()) {
    throw std::invalid_argument("Unable to open file for reading");
  }
  // Read the entry from the file
  auto entry_size = pir_params_.get_entry_size();
  in_file.seekg(entry_idx * entry_size);
  Entry entry(entry_size);
  in_file.read(reinterpret_cast<char *>(entry.data()), entry_size);
  in_file.close();

  return entry;
}


std::vector<seal::Ciphertext> PirServer::make_query(const uint32_t client_id, PirQuery &&query) {

  // ========================== Expansion & conversion ==========================
  // Query expansion
  TIME_START(EXPAND_TIME);
  std::vector<seal::Ciphertext> query_vector = expand_query(client_id, query);
  TIME_END(EXPAND_TIME);

  // Reconstruct RGSW queries
  TIME_START(CONVERT_TIME);
  std::vector<GSWCiphertext> gsw_vec; // GSW ciphertexts
  if (dims_.size() != 1) {  // if we do need futher dimensions
    gsw_vec.resize(dims_.size() - 1);
    for (size_t i = 1; i < dims_.size(); i++) {
      std::vector<seal::Ciphertext> lwe_vector; // BFV ciphertext, size l * 2. This vector will be reconstructed as a single RGSW ciphertext.
      for (size_t k = 0; k < DatabaseConstants::GSW_L; k++) {
        auto ptr = dims_[0] + (i - 1) * DatabaseConstants::GSW_L + k;
        lwe_vector.push_back(query_vector[ptr]);
      }
      // Converting the BFV ciphertext to GSW ciphertext
      key_gsw.query_to_gsw(lwe_vector, client_gsw_keys_[client_id], gsw_vec[i - 1]);
    }
  }
  TIME_END(CONVERT_TIME);

  // ========================== Evaluations ==========================
  // Evaluate the first dimension
  TIME_START(FST_DIM_TIME);
  std::vector<seal::Ciphertext> result = evaluate_first_dim(query_vector);
  TIME_END(FST_DIM_TIME);

  // Evaluate the other dimensions
  TIME_START(OTHER_DIM_TIME);
  if (dims_.size() != 1) {
    for (size_t i = 1; i < dims_.size(); i++) {
      evaluate_gsw_product(result, gsw_vec[i - 1]);
    }
  }
  TIME_END(OTHER_DIM_TIME);

  // ========================== Post-processing ==========================
  // modulus switching so to reduce the response size by half
  if(pir_params_.get_seal_params().coeff_modulus().size() > 2) {
    DEBUG_PRINT("Modulus switching...");
    evaluator_.mod_switch_to_next_inplace(result[0]); // result.size() == 1.
  }

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
  const auto num_bits_per_plaintext = pir_params_.get_num_bits_per_plaintext();
  const auto num_entries_per_plaintext = pir_params_.get_num_entries_per_plaintext();
  const auto num_plaintexts = chunk_entry.size() / num_entries_per_plaintext;  // number of plaintexts in the new chunk
  const uint128_t coeff_mask = (uint128_t(1) << (bits_per_coeff)) - 1;  // bits_per_coeff many 1s
  
  const auto fst_dim_sz = pir_params_.get_fst_dim_sz();  // number of plaintexts in the first dimension
  const auto chunk_offset = fst_dim_sz * chunk_idx;  // offset for the current chunk

  // Now we handle plaintexts one by one.
  for (size_t i = 0; i < num_plaintexts; i++) {
    seal::Plaintext plaintext(DatabaseConstants::PolyDegree);

    // Loop through the entries that corresponds to the current plaintext. 
    // Then calculate the total size (in bytes) of this plaintext.
    // NOTE: it is possible that some entry is empty, which has size 0.
    size_t additive_sum_size = 0;
    for (size_t j = num_entries_per_plaintext * i;
         j < std::min(num_entries_per_plaintext * (i + 1), chunk_entry.size()); j++) {
      additive_sum_size += chunk_entry[j].size();
    }

    if (additive_sum_size == 0) {
      continue; // leave std::nullopt in the chunk if the plaintext is empty.
    }

    size_t index = 0;  // index for the current coefficient to be filled
    uint128_t data_buffer = 0;
    size_t data_offset = 0;
    // For each entry in the current plaintext
    for (size_t j = num_entries_per_plaintext * i;
         j < std::min(num_entries_per_plaintext * (i + 1), chunk_entry.size()); j++) {
      // For each byte in this entry
      for (size_t k = 0; k < pir_params_.get_entry_size(); k++) {
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
  for (size_t i = 0; i < num_pt_; ++i) {
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
  // number of rns modulus
  auto num_mods = pir_params_.get_ct_coeff_mod_cnt();
  DEBUG_PRINT("Number of RNS moduli: " << num_mods);
  auto other_dim_sz = pir_params_.get_other_dim_sz();
  // number of uint128_t we need to store in the intermediate result
  auto elem_cnt = other_dim_sz * DatabaseConstants::PolyDegree * num_mods * 2;
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


std::vector<size_t> PirServer::get_dims() const {
  return dims_;
}
