#include "gsw_eval.h"
#include "utils.h"
#include "seal/util/rlwe.h"
#include "logging.h"
#include <cassert>

// Here we compute a cross product between the transpose of the decomposed BFV
// (a 2l vector of polynomials) and the GSW ciphertext (a 2lx2 matrix of
// polynomials) to obtain a size-2 vector of polynomials, which is exactly our
// result ciphertext. We use an NTT multiplication to speed up polynomial
// multiplication, assuming that both the GSWCiphertext and decomposed bfv is in
// polynomial coefficient representation.


void GSWEval::gsw_ntt_negacyclic_harvey(GSWCiphertext &gsw) {
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();
  const auto context = pir_params_.get_context();
  auto ntt_tables = context.first_context_data()->small_ntt_tables();

  for (auto &poly : gsw) {
    seal::util::CoeffIter gsw_poly_ptr(poly.data());
    for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {
      seal::util::ntt_negacyclic_harvey(gsw_poly_ptr + coeff_count * mod_id, *(ntt_tables + mod_id));
    }
    seal::util::CoeffIter gsw_poly_ptr2(poly.data() + coeff_count * rns_mod_cnt);
    for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {
      seal::util::ntt_negacyclic_harvey(gsw_poly_ptr2 + coeff_count * mod_id, *(ntt_tables + mod_id));
    }
  }
}

void GSWEval::external_product(GSWCiphertext const &gsw_enc, seal::Ciphertext const &bfv,
                              seal::Ciphertext &res_ct) {
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const size_t coeff_val_cnt = DatabaseConstants::PolyDegree * rns_mod_cnt; // polydegree * RNS moduli count

  // Decomposing the BFV ciphertext to 2l polynomials. Transform to NTT form.
  std::vector<std::vector<uint64_t>> decomposed_bfv;
  TIME_START(DECOMP_RLWE_TIME);
  if (rns_mod_cnt == 1) {
    decomp_rlwe_single_mod(bfv, decomposed_bfv);
  } else {
    decomp_rlwe(bfv, decomposed_bfv);
  }
  TIME_END(DECOMP_RLWE_TIME);

  std::vector<std::vector<uint128_t>> result(
      2, std::vector<uint128_t>(coeff_val_cnt, 0));

  TIME_START(EXTERN_PROD_MAT_MULT_TIME);
  // matrix multiplication: decomp(bfv) * gsw = (1 x 2l) * (2l x 2) = (1 x 2)
  for (size_t k = 0; k < 2; ++k) {
    for (size_t j = 0; j < 2 * l_; j++) {
      seal::util::ConstCoeffIter encrypted_gsw_ptr(gsw_enc[j].data() + k * coeff_val_cnt);
      seal::util::ConstCoeffIter encrypted_rlwe_ptr(decomposed_bfv[j]);
      #pragma GCC unroll 32
      for (size_t i = 0; i < coeff_val_cnt; i++) {
        result[k][i] += static_cast<uint128_t>(encrypted_rlwe_ptr[i]) * encrypted_gsw_ptr[i];
      }
    }
  }
  TIME_END(EXTERN_PROD_MAT_MULT_TIME);

  // taking mods.
  const auto coeff_modulus = pir_params_.get_coeff_modulus();
  for (size_t poly_id = 0; poly_id < 2; poly_id++) {
    auto ct_ptr = res_ct.data(poly_id);
    auto pt_ptr = result[poly_id];

    for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {
      auto mod_idx = (mod_id * coeff_count);
      for (size_t coeff_id = 0; coeff_id < coeff_count; coeff_id++) {
        auto x = pt_ptr[coeff_id + mod_idx];
        uint64_t raw[2] = {static_cast<uint64_t>(x), static_cast<uint64_t>(x >> 64)};
        ct_ptr[coeff_id + mod_idx] = util::barrett_reduce_128(raw, coeff_modulus[mod_id]);
      }
    }
  }
  res_ct.is_ntt_form() = true;  // the result of two NTT form polynomials is still in NTT form.
}

void GSWEval::decomp_rlwe(seal::Ciphertext const &ct, std::vector<std::vector<uint64_t>> &output) {
  assert(output.size() == 0);
  output.reserve(2 * l_);

  // Setup parameters
  const uint128_t base = uint128_t(1) << base_log2_;
  const uint128_t mask = base - 1;
  const auto &coeff_modulus = pir_params_.get_coeff_modulus();
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();
  const auto &context = pir_params_.get_context();
  const auto context_data = context.first_context_data();
  auto ntt_tables = context_data->small_ntt_tables();
  seal::util::RNSBase *rns_base = context_data->rns_tool()->base_q();
  auto pool = seal::MemoryManager::GetPool();

  std::vector<uint64_t> data(coeff_count * rns_mod_cnt);
  for (size_t poly_id = 0; poly_id < 2; poly_id++) {
    const uint64_t *poly_ptr = ct.data(poly_id);

    memcpy(data.data(), poly_ptr, coeff_count * rns_mod_cnt * sizeof(uint64_t));
    rns_base->compose_array(data.data(), coeff_count, pool);

    for (int p = l_ - 1; p >= 0; p--) {
      std::vector<uint64_t> row = data;
      TIME_START(RIGHT_SHIFT_TIME);
      for (size_t k = 0; k < coeff_count; k++) {
        auto ptr = row.data() + k * rns_mod_cnt;
        // ! This right shift is very time consuming. About 3 times slower than the actual external product multiplication.
        seal::util::right_shift_uint(ptr, p * base_log2_, rns_mod_cnt, ptr); // shift right by p * base_log2
        ptr[0] &= mask;
        for (size_t i = 1; i < rns_mod_cnt; i++) {
          ptr[i] = 0;
        }
      }
      TIME_END(RIGHT_SHIFT_TIME);
      rns_base->decompose_array(row.data(), coeff_count, pool);

      TIME_START(EXTERN_NTT_TIME);
      // transform to NTT form
      for (size_t i = 0; i < rns_mod_cnt; i++) {
        seal::util::ntt_negacyclic_harvey(row.data() + coeff_count * i, ntt_tables[i]);
      }
      TIME_END(EXTERN_NTT_TIME);
      output.emplace_back(std::move(row));
    }
  }
}

void GSWEval::decomp_rlwe_single_mod(seal::Ciphertext const &ct, std::vector<std::vector<uint64_t>> &output) {

  assert(output.size() == 0);
  output.reserve(2 * l_);

  // Get parameters
  const uint64_t base = uint64_t(1) << base_log2_;
  const uint64_t mask = base - 1;

  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const auto &context = pir_params_.get_context();
  const auto &context_data = context.first_context_data();
  const auto ntt_tables = context_data->small_ntt_tables();
  constexpr size_t bits_per_uint64_sz = 64;

  seal::util::RNSBase *rns_base = context_data->rns_tool()->base_q();
  auto pool = seal::MemoryManager::GetPool();

  std::vector<uint64_t> data(coeff_count);

  for (size_t poly_id = 0; poly_id < 2; poly_id++) {
    const uint64_t *poly_ptr = ct.data(poly_id);
    memcpy(data.data(), poly_ptr, coeff_count * sizeof(uint64_t));
    rns_base->compose_array(data.data(), coeff_count, pool);

    for (int p = l_ - 1; p >= 0; p--) {
      std::vector<uint64_t> row = data;
      const size_t shift_amount = p * base_log2_;
      TIME_START(RIGHT_SHIFT_TIME);
      for (size_t k = 0; k < coeff_count; k++) {
        row[k] = (row[k] >> shift_amount) & mask;
      }
      TIME_END(RIGHT_SHIFT_TIME);
      rns_base->decompose_array(row.data(), coeff_count, pool);

      // transform to NTT form
      TIME_START(EXTERN_NTT_TIME);
      seal::util::ntt_negacyclic_harvey(row.data(), ntt_tables[0]);
      TIME_END(EXTERN_NTT_TIME);
      output.emplace_back(std::move(row));
    }
  }
}

void GSWEval::query_to_gsw(std::vector<seal::Ciphertext> query, GSWCiphertext gsw_key,
                           GSWCiphertext &output) {
  const size_t curr_l = query.size();
  output.resize(curr_l);
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();

  for (size_t i = 0; i < curr_l; i++) {
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i].push_back(query[i].data(0)[j]);
    }
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i].push_back(query[i].data(1)[j]);
    }
  }
  gsw_ntt_negacyclic_harvey(output);
  output.resize(2 * curr_l);
  for (size_t i = 0; i < curr_l; i++) {
    external_product(gsw_key, query[i], query[i]);
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i + curr_l].push_back(query[i].data(0)[j]);
    }
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i + curr_l].push_back(query[i].data(1)[j]);
    }
  }
}

void GSWEval::plain_to_gsw(std::vector<uint64_t> const &plaintext,
                  seal::Encryptor const &encryptor, seal::SecretKey const &sk,
                  std::vector<seal::Ciphertext> &output) { 
  output.clear(); 
  // when poly_id = 0, we are working on the first half of the GSWCiphertext
  for (size_t poly_id = 0; poly_id < 2; poly_id++) {
    for (size_t k = 0; k < l_; k++) {
      seal::Ciphertext cipher; 
      plain_to_gsw_one_row(plaintext, encryptor, sk, poly_id, k, cipher);
      output.push_back(cipher);
    }
  }
}

void GSWEval::plain_to_gsw_one_row(std::vector<uint64_t> const &plaintext,
                                   seal::Encryptor const &encryptor,
                                   seal::SecretKey const &sk,
                                   const size_t level, const size_t half,
                                   seal::Ciphertext &output) {

  // Accessing context data within this function instead of passing these parameters
  const auto &context = pir_params_.get_context();
  const auto &params = context.first_context_data()->parms();
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const auto &coeff_modulus = pir_params_.get_coeff_modulus();
  const size_t rns_mod_cnt = coeff_modulus.size();
  assert(plaintext.size() == coeff_count * rns_mod_cnt || plaintext.size() == coeff_count);

  // Create RGSW gadget.
  // We only use rns_mod_cnt many gadgets, we can save some computation here. Not a big deal though. 
  std::vector<std::vector<uint64_t>> gadget = gsw_gadget(l_, base_log2_, rns_mod_cnt, coeff_modulus);

  encryptor.encrypt_zero_symmetric_seeded(output);
  auto ct = output.data(half);
  // Many(2) moduli are used
  for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {
    size_t pad = (mod_id * coeff_count);
    uint128_t mod = coeff_modulus[mod_id].value();
    uint64_t gadget_coef = gadget[mod_id][level];
    auto pt = plaintext.data();
    if (plaintext.size() == coeff_count * rns_mod_cnt) {
      pt = plaintext.data() + pad;
    }
    // Loop through plaintext coefficients
    for (size_t j = 0; j < coeff_count; j++) {
      // TODO: We can use barret reduction here.
      uint128_t val = (uint128_t)pt[j] * gadget_coef % mod;
      ct[j + pad] =
          static_cast<uint64_t>((ct[j + pad] + val) % mod);
    }
  }
}

void GSWEval::plain_to_half_gsw(std::vector<uint64_t> const &plaintext,
                                   seal::Encryptor const &encryptor,
                                   seal::SecretKey const &sk,
                                   std::vector<seal::Ciphertext> &output) {
  output.clear();
  // when poly_id = 0, we are working on the first half of the GSWCiphertext
  for (size_t poly_id = 0; poly_id < 2; poly_id++) {
    for (size_t k = 0; k < l_; k++) {
      seal::Ciphertext cipher; 
      plain_to_half_gsw_one_row(plaintext, encryptor, sk, poly_id, k, cipher);
      output.push_back(cipher);
    }
  }
}

void GSWEval::plain_to_half_gsw_one_row(std::vector<uint64_t> const &plaintext,
                                  seal::Encryptor const &encryptor,
                                  seal::SecretKey const &sk, const size_t half,
                                  const size_t level, seal::Ciphertext &output) {

  // Accessing context data within this function instead of passing these parameters
  const auto &context = pir_params_.get_context();
  const auto &params = context.first_context_data()->parms();
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const auto &coeff_modulus = pir_params_.get_coeff_modulus();
  const size_t rns_mod_cnt = coeff_modulus.size();
  assert(plaintext.size() == coeff_count * rns_mod_cnt || plaintext.size() == coeff_count);

  // Create RGSW gadget.
  std::vector<std::vector<uint64_t>> gadget = gsw_gadget(l_, base_log2_, rns_mod_cnt, coeff_modulus);

  // ================== Second half of the seeded GSW ==================
  if (half == 1) {
    // extract the level column of gadget
    std::vector<uint64_t> col;
    for (size_t i = 0; i < rns_mod_cnt; i++) {
      col.push_back(gadget[i][level]);
    }
    seal::util::prepare_seeded_gsw_key(sk, col, context,
                                       params.parms_id(), false, output);
    return;
  }

  // ================== Other cases ==================
  assert(half == 0);
  // If we are at the first half of the GSW, we are adding new things to the
  // first polynomial (c0) of the given BFV ciphertext. c1 is not touched.
  encryptor.encrypt_zero_symmetric_seeded(output);
  auto ct = output.data(0);
  // Many(2) moduli are used
  for (size_t mod_id = 0; mod_id < rns_mod_cnt; mod_id++) {
    size_t pad = (mod_id * coeff_count);
    uint128_t mod = coeff_modulus[mod_id].value();
    uint64_t gadget_coef = gadget[mod_id][level];
    auto pt = plaintext.data();
    if (plaintext.size() == coeff_count * rns_mod_cnt) {
      pt = plaintext.data() + pad;
    }
    // Loop through plaintext coefficients
    for (size_t j = 0; j < coeff_count; j++) {
      // TODO: We can use barret reduction here.
      uint128_t val = (uint128_t)pt[j] * gadget_coef % mod;
      ct[j + pad] =
          static_cast<uint64_t>((ct[j + pad] + val) % mod);
    }
  }
}

void GSWEval::sealGSWVecToGSW(GSWCiphertext &output, const std::vector<seal::Ciphertext> &gsw_vec) {
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const size_t rns_mod_cnt = pir_params_.get_rns_mod_cnt();

  output.clear();
  for (auto &ct : gsw_vec) {
    std::vector<uint64_t> row;
    for (size_t i = 0; i < coeff_count * rns_mod_cnt; i++) {
      row.push_back(ct.data(0)[i]);
    }
    for (size_t i = 0; i < coeff_count * rns_mod_cnt; i++) {
      row.push_back(ct.data(1)[i]);
    }
    output.push_back(row);
  }
}