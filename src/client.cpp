#include "client.h"
#include "pir.h"
#include "utils.h"
#include "gsw_eval.h"


// constructor
PirClient::PirClient(const PirParams &pir_params)
    : client_id_(0), context_(pir_params.get_seal_params()), keygen_(context_),
      secret_key_(keygen_.secret_key()), decryptor_(context_, secret_key_),
      encryptor_(context_, secret_key_), evaluator_(context_),
      pir_params_(pir_params), dims_(pir_params.get_dims()) {
  // set the client id to 0
  client_id_ = rand();
}

std::vector<Ciphertext> PirClient::generate_gsw_from_key() {
  std::vector<seal::Ciphertext> gsw_enc; // temporary GSW ciphertext using seal::Ciphertext
  const auto sk_ = secret_key_.data();
  const auto ntt_tables = context_.first_context_data()->small_ntt_tables();
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  std::vector<uint64_t> sk_ntt(sk_.data(), sk_.data() + coeff_count * rns_mod_cnt);

  RNSIter secret_key_iter(sk_ntt.data(), coeff_count);
  inverse_ntt_negacyclic_harvey(secret_key_iter, rns_mod_cnt, ntt_tables);

  GSWEval key_gsw(pir_params_, pir_params_.get_l_key(), pir_params_.get_base_log2_key());
  key_gsw.plain_to_half_gsw(sk_ntt, encryptor_, secret_key_, gsw_enc);
  return gsw_enc;
}

size_t PirClient::get_database_plain_index(size_t entry_index) {
  return entry_index / pir_params_.get_num_entries_per_plaintext();
}

std::vector<size_t> PirClient::get_query_indices(size_t plaintext_index) {
  std::vector<size_t> query_indices;
  const size_t col_idx = plaintext_index % dims_[0];  // the first dimension
  size_t row_idx = plaintext_index / dims_[0];  // the rest of the dimensions
  size_t remain_pt_num = pir_params_.get_num_pt() / dims_[0];

  query_indices.push_back(col_idx);
  for (size_t i = 1; i < dims_.size(); i++) {
    size_t dim_size = dims_[i];
    remain_pt_num /= dim_size;
    query_indices.push_back(row_idx / remain_pt_num);
    row_idx = row_idx % remain_pt_num;
  }

  return query_indices;
}



seal::Ciphertext PirClient::generate_query(const size_t entry_index) {

  // ================== Setup parameters ==================
  // Get the corresponding index of the plaintext in the database
  const size_t plaintext_index = get_database_plain_index(entry_index);
  std::vector<size_t> query_indices = get_query_indices(plaintext_index);
  PRINT_INT_ARRAY("\t\tquery_indices", query_indices.data(), query_indices.size());
  const size_t bits_per_ciphertext = 1 << pir_params_.get_expan_height(); // padding msg_size to the next power of 2

  // Algorithm 1 from the OnionPIR Paper

  // empty plaintext
  seal::Plaintext plain_query(DatabaseConstants::PolyDegree); 
  // We set the corresponding coefficient to the inverse so the value of the
  // expanded ciphertext will be 1
  uint64_t inverse = 0;
  const uint64_t plain_modulus = pir_params_.get_plain_mod();
  seal::util::try_invert_uint_mod(bits_per_ciphertext, plain_modulus, inverse);

  // Add the first dimension query vector to the query
  plain_query[ query_indices[0] ] = inverse;
  
  // Encrypt plain_query first. Later we will insert the rest. $\tilde c$ in paper
  seal::Ciphertext query;
  encryptor_.encrypt_symmetric_seeded(plain_query, query);

  // no further dimensions
  if (query_indices.size() == 1) {
    DEBUG_PRINT("No further dimensions");
    return query;
  }

  // ================== Add GSW values to the query ==================
  const size_t l = pir_params_.get_l();
  const size_t base_log2 = pir_params_.get_base_log2();
  const auto coeff_modulus = pir_params_.get_coeff_modulus();
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();

  // The following two for-loops calculates the powers for GSW gadgets.
  std::vector<uint128_t> inv(rns_mod_cnt);
  for (size_t k = 0; k < rns_mod_cnt; k++) {
    uint64_t result;
    seal::util::try_invert_uint_mod(bits_per_ciphertext, coeff_modulus[k], result);
    inv[k] = result;
  }

  // rns_mod_cnt many rows, each row is B^{l-1},, ..., B^0 under different moduli
  std::vector<std::vector<uint64_t>> gadget = utils::gsw_gadget(l, base_log2, rns_mod_cnt, coeff_modulus);

  // This for-loop corresponds to the for-loop in Algorithm 1 from the OnionPIR paper
  auto q_head = query.data(0); // points to the first coefficient of the first ciphertext(c0) 
  for (size_t i = 1; i < query_indices.size(); i++) {  // dimensions
    // we use this if statement to replce the j for loop in Algorithm 1. This is because N_i = 2 for all i > 0
    // When 0 is requested, we use initial encrypted value of seal::Ciphertext query, where the coefficients decrypts to 0. 
    // When 1 is requested, we add special values to the coefficients of the query so that they decrypts to correct GSW(1) values.
    if (query_indices[i] == 1) {
      for (size_t k = 0; k < l; k++) {
        for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {
          const size_t pad = mod_id * DatabaseConstants::PolyDegree;   // We use two moduli for the same gadget value. They are apart by coeff_count.
          const size_t coef_pos = dims_[0] + (i-1) * l + k + pad;  // the position of the coefficient in the query
          uint128_t mod = coeff_modulus[mod_id].value();
          // the coeff is (B^{l-1}, ..., B^0) / bits_per_ciphertext
          uint128_t coef = gadget[mod_id][k] * inv[mod_id] % mod;
          q_head[coef_pos] = (q_head[coef_pos] + coef) % mod;
        }
      }
    }
  }

  return query;
}


seal::Ciphertext PirClient::fast_generate_query(const size_t entry_index) {
  // ================== Setup parameters ==================
  // Get the corresponding index of the plaintext in the database
  const size_t plaintext_index = get_database_plain_index(entry_index);
  std::vector<size_t> query_indices = get_query_indices(plaintext_index);
  PRINT_INT_ARRAY("\t\tquery_indices", query_indices.data(), query_indices.size());
  const size_t expan_height = pir_params_.get_expan_height();
  const size_t bits_per_ciphertext = 1 << expan_height; // padding msg_size to the next power of 2

  // Algorithm 1 from the OnionPIR Paper

  // empty plaintext
  seal::Plaintext plain_query(DatabaseConstants::PolyDegree); // we allow 4096 coefficients in the plaintext polynomial to be set as suggested in the paper.
  // We set the corresponding coefficient to the inverse so the value of the
  // expanded ciphertext will be 1
  uint64_t inverse = 0;
  const uint64_t plain_modulus = pir_params_.get_plain_mod();
  seal::util::try_invert_uint_mod(bits_per_ciphertext, plain_modulus, inverse);

  // Since we are using new expansion method, where index b is split to even
  // part and odd part, even part is stored in index 2b, and odd part is stored
  // in index 2b+1. This results in a bit-reversed order of the indices.
  const size_t reversed_index = utils::bit_reverse(query_indices[0], expan_height);
  plain_query[ reversed_index ] = inverse; // Add the first dimension query vector to the query
  DEBUG_PRINT("reversed_index: " << reversed_index << ", query_indices[0]: " << query_indices[0]);

  // Encrypt plain_query first. Later we will insert the rest. $\tilde c$ in paper
  seal::Ciphertext query;
  encryptor_.encrypt_symmetric_seeded(plain_query, query);

  // add gsw values to the query bfv 
  add_gsw_to_query(query, query_indices);

  return query;
}


void PirClient::add_gsw_to_query(seal::Ciphertext &query, const std::vector<size_t> query_indices) {
  // no further dimensions
  if (query_indices.size() == 1) { return; }
  const size_t expan_height = pir_params_.get_expan_height();
  const size_t bits_per_ciphertext = 1 << expan_height; // padding msg_size to the next power of 2
  const size_t l = pir_params_.get_l();
  const size_t base_log2 = pir_params_.get_base_log2();
  const auto coeff_modulus = pir_params_.get_coeff_modulus();
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();
  const size_t fst_dim_sz = dims_[0];

  // The following two for-loops calculates the powers for GSW gadgets.
  std::vector<uint128_t> inv(rns_mod_cnt);
  for (size_t k = 0; k < rns_mod_cnt; k++) {
    uint64_t result;
    seal::util::try_invert_uint_mod(bits_per_ciphertext, coeff_modulus[k], result);
    inv[k] = result;
  }

  // rns_mod_cnt many rows, each row is B^{l-1},, ..., B^0 under different moduli
  std::vector<std::vector<uint64_t>> gadget = utils::gsw_gadget(l, base_log2, rns_mod_cnt, coeff_modulus);

  // This for-loop corresponds to the for-loop in Algorithm 1 from the OnionPIR paper
  auto q_head = query.data(0); // points to the first coefficient of the first ciphertext(c0) 
  for (size_t i = 1; i < query_indices.size(); i++) {  // dimensions
    // we use this if statement to replce the j for loop in Algorithm 1. This is because N_i = 2 for all i > 0
    // When 0 is requested, we use initial encrypted value of seal::Ciphertext query, where the coefficients decrypts to 0. 
    // When 1 is requested, we add special values to the coefficients of the query so that they decrypts to correct GSW(1) values.
    if (query_indices[i] == 1) {
      for (size_t k = 0; k < l; k++) {
        const size_t coef_pos = fst_dim_sz + (i-1) * l + k;  // the position of the coefficient in the resulting query
        const size_t reversed_idx = utils::bit_reverse(coef_pos, expan_height);  // the position of the coefficient in the query
        DEBUG_PRINT("reversed_idx: " << reversed_idx << ", coef_pos: " << coef_pos);
        for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {
          const size_t pad = mod_id * DatabaseConstants::PolyDegree;   // We use two moduli for the same gadget value. They are apart by coeff_count.
          uint128_t mod = coeff_modulus[mod_id].value();
          // the coeff is (B^{l-1}, ..., B^0) / bits_per_ciphertext
          uint128_t coef = gadget[mod_id][k] * inv[mod_id] % mod;
          q_head[reversed_idx + pad] = (q_head[reversed_idx + pad] + coef) % mod;
        }
      }
    }
  }
}




size_t PirClient::write_query_to_stream(const seal::Ciphertext &query, std::stringstream &data_stream) {
  return query.save(data_stream);
}

size_t PirClient::write_gsw_to_stream(const std::vector<Ciphertext> &gsw, std::stringstream &gsw_stream) {
  size_t total_size = 0;
  for (auto &ct : gsw) {
    size_t size = ct.save(gsw_stream);
    total_size += size;
  }
  return total_size;
}

size_t PirClient::create_galois_keys(std::stringstream &galois_key_stream) {
  std::vector<uint32_t> galois_elts;

  // This is related to the unpacking algorithm.
  // expansion height is the height of the expansion tree such that
  // 2^get_expan_height() is equal to the number of packed values padded to the next power of 2.
  const size_t expan_height = pir_params_.get_expan_height();
  const size_t poly_degree = DatabaseConstants::PolyDegree;
  for (size_t i = 0; i < expan_height; i++) {
    galois_elts.push_back(1 + (poly_degree >> i));
  }
  // PRINT_INT_ARRAY("galois_elts: ", galois_elts, galois_elts.size());
  auto written_size = keygen_.create_galois_keys(galois_elts).save(galois_key_stream);
  return written_size;
}

std::vector<seal::Plaintext> PirClient::decrypt_result(const std::vector<seal::Ciphertext> reply) {
  std::vector<seal::Plaintext> result(reply.size(), seal::Plaintext(DatabaseConstants::PolyDegree));
  for (size_t i = 0; i < reply.size(); i++) {
    decryptor_.decrypt(reply[i], result[i]);
  }

  return result;
}

Entry PirClient::get_entry_from_plaintext(const size_t entry_index, const seal::Plaintext plaintext) const {
  // Offset in the plaintext in bits
  const size_t start_position_in_plaintext = (entry_index % pir_params_.get_num_entries_per_plaintext()) *
                                       pir_params_.get_entry_size() * 8;

  // Offset in the plaintext by coefficient
  const size_t num_bits_per_coeff = pir_params_.get_num_bits_per_coeff();
  size_t coeff_index = start_position_in_plaintext / num_bits_per_coeff;

  // Offset in the coefficient by bits
  const size_t coeff_offset = start_position_in_plaintext % num_bits_per_coeff;

  // Size of entry in bits
  const size_t entry_size = pir_params_.get_entry_size();
  Entry result;

  uint128_t data_buffer = plaintext.data()[coeff_index] >> coeff_offset;
  uint128_t data_offset = num_bits_per_coeff - coeff_offset;

  while (result.size() < entry_size) {
    if (data_offset >= 8) {
      result.push_back(data_buffer & 0xFF);
      data_buffer >>= 8;
      data_offset -= 8;
    } else {
      coeff_index += 1;
      uint128_t next_buffer = plaintext.data()[coeff_index];
      data_buffer |= next_buffer << data_offset;
      data_offset += num_bits_per_coeff;
    }
  }

  return result;
}

