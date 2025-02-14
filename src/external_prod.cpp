#include "external_prod.h"
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

GSWEval data_gsw, key_gsw;

void GSWEval::gsw_ntt_negacyclic_harvey(GSWCiphertext &gsw) {
  const auto &context_data = context->first_context_data();
  auto &parms2 = context->first_context_data()->parms();
  auto &coeff_modulus = parms2.coeff_modulus();
  size_t coeff_count = parms2.poly_modulus_degree();
  size_t rns_mod_cnt = coeff_modulus.size();
  auto ntt_tables = context_data->small_ntt_tables();

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

void GSWEval::ciphertext_inverse_ntt(seal::Ciphertext &ct) {
  const auto &context_data = context->first_context_data();
  auto &parms2 = context_data->parms();
  auto &coeff_modulus = parms2.coeff_modulus();
  size_t coeff_count = parms2.poly_modulus_degree();
  size_t rns_mod_cnt = coeff_modulus.size();
  auto ntt_tables = context_data->small_ntt_tables();

  for (size_t i = 0; i < rns_mod_cnt; i++) {
    seal::util::inverse_ntt_negacyclic_harvey(ct.data(0) + coeff_count * i, *(ntt_tables + i));
  }
  for (size_t i = 0; i < rns_mod_cnt; i++) {
    seal::util::inverse_ntt_negacyclic_harvey(ct.data(1) + coeff_count * i, *(ntt_tables + i));
  }
}

void GSWEval::external_product(GSWCiphertext const &gsw_enc, seal::Ciphertext const &bfv,
                              seal::Ciphertext &res_ct) {

  const auto &coeff_modulus = context->first_context_data()->parms().coeff_modulus();
  const size_t rns_mod_cnt = coeff_modulus.size();  // rns mod count
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  const size_t coeff_val_cnt = DatabaseConstants::PolyDegree * rns_mod_cnt; // polydegree * RNS moduli count

  // Decomposing the BFV ciphertext to 2l polynomials. Transform to NTT form.
  std::vector<std::vector<uint64_t>> decomposed_bfv;
  TIME_START(DECOMP_RLWE_TIME);
  decomp_rlwe(bfv, decomposed_bfv);
  TIME_END(DECOMP_RLWE_TIME);

  std::vector<std::vector<uint128_t>> result(
      2, std::vector<uint128_t>(coeff_val_cnt, 0));

  TIME_START(EXTERN_PROD_MAT_MULT_TIME);
  // matrix multiplication: decomp(bfv) * gsw = (1 x 2l) * (2l x 2) = (1 x 2)
  for (size_t k = 0; k < 2; ++k) {
    for (size_t j = 0; j < 2 * l; j++) {
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
}

void GSWEval::decomp_rlwe(seal::Ciphertext const &ct, std::vector<std::vector<uint64_t>> &output) {

  assert(output.size() == 0);
  output.reserve(2 * l);

  // Get parameters
  const uint128_t base = uint128_t(1) << base_log2;
  const uint128_t mask = base - 1;

  const auto &context_data = context->first_context_data();
  const auto &parms = context_data->parms();
  const auto &coeff_modulus = parms.coeff_modulus();
  const size_t coeff_count = parms.poly_modulus_degree();
  const size_t rns_mod_cnt = coeff_modulus.size();
  const auto ntt_tables = context_data->small_ntt_tables();


  seal::util::RNSBase *rns_base = context_data->rns_tool()->base_q();
  auto pool = seal::MemoryManager::GetPool();

  std::vector<uint64_t> data(coeff_count * rns_mod_cnt);

  for (size_t poly_id = 0; poly_id < 2; poly_id++) {
    const uint64_t *poly_ptr = ct.data(poly_id);

    memcpy(data.data(), poly_ptr, coeff_count * rns_mod_cnt * sizeof(uint64_t));
    rns_base->compose_array(data.data(), coeff_count, pool);

    for (size_t p = l - 1; p >= 0; p--) {
      std::vector<uint64_t> row = data;
      for (size_t k = 0; k < coeff_count; k++) {
        auto ptr = row.data() + k * rns_mod_cnt;
        // ! This right shift is very time consuming. About 3 times slower than the actual external product multiplication.
        seal::util::right_shift_uint(ptr, p * base_log2, rns_mod_cnt, ptr); // shift right by p * base_log2
        ptr[0] &= mask;
        for (size_t i = 1; i < rns_mod_cnt; i++) {
          ptr[i] = 0;
        }
      }
      rns_base->decompose_array(row.data(), coeff_count, pool);

      // transform to NTT form
      for (size_t i = 0; i < rns_mod_cnt; i++) {
        seal::util::ntt_negacyclic_harvey(row.data() + coeff_count * i, ntt_tables[i]);
      }

      output.emplace_back(std::move(row));
    }
  }
}

void GSWEval::query_to_gsw(std::vector<seal::Ciphertext> query, GSWCiphertext gsw_key,
                           GSWCiphertext &output) {
  size_t cl = query.size();
  assert(output.size() == 0);
  output.resize(cl);

  const auto &context_data = context->get_context_data(query[0].parms_id());
  auto &parms = context_data->parms();
  auto &coeff_modulus = parms.coeff_modulus();
  size_t coeff_count = parms.poly_modulus_degree();
  size_t rns_mod_cnt = coeff_modulus.size();

  for (size_t i = 0; i < cl; i++) {
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i].push_back(query[i].data(0)[j]);
    }
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i].push_back(query[i].data(1)[j]);
    }
  }
  gsw_ntt_negacyclic_harvey(output);
  output.resize(2 * cl);
  for (size_t i = 0; i < cl; i++) {
    external_product(gsw_key, query[i], query[i]);
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i + cl].push_back(query[i].data(0)[j]);
    }
    for (size_t j = 0; j < coeff_count * rns_mod_cnt; j++) {
      output[i + cl].push_back(query[i].data(1)[j]);
    }
  }
}

void GSWEval::encrypt_plain_to_gsw(std::vector<uint64_t> const &plaintext,
                                   seal::Encryptor const &encryptor,
                                   seal::SecretKey const &sk,
                                   std::vector<seal::Ciphertext> &output) {
  output.clear();
  // when poly_id = 0, we are working on the first half of the GSWCiphertext
  for (size_t poly_id = 0; poly_id < 2; poly_id++) {
    for (size_t k = 0; k < l; k++) {
      seal::Ciphertext cipher =
          enc_plain_to_gsw_one_row(plaintext, encryptor, sk, poly_id, k);
      output.push_back(cipher);
    }
  }
}

seal::Ciphertext
GSWEval::enc_plain_to_gsw_one_row(std::vector<uint64_t> const &plaintext,
                                  seal::Encryptor const &encryptor,
                                  seal::SecretKey const &sk, const size_t half,
                                  const size_t level) {

  // Accessing context data within this function instead of passing these parameters
  const auto &parms = context->first_context_data()->parms();
  size_t coeff_count = parms.poly_modulus_degree();
  size_t rns_mod_cnt = parms.coeff_modulus().size();
  const auto &coeff_modulus = parms.coeff_modulus();
  assert(plaintext.size() == coeff_count * rns_mod_cnt || plaintext.size() == coeff_count);

  // Create RGSW gadget.
  std::vector<std::vector<uint64_t>> gadget = gsw_gadget(l, base_log2, rns_mod_cnt, coeff_modulus);

  // ================== Second half of the seeded GSW ==================
  if (half == 1) {
    seal::Ciphertext cipher;
    // extract the level column of gadget
    std::vector<uint64_t> col;
    for (size_t i = 0; i < rns_mod_cnt; i++) {
      col.push_back(gadget[i][level]);
    }
    seal::util::prepare_seeded_gsw_key(sk, col, *context,
                                       parms.parms_id(), false, cipher);
    return cipher;
  }

  // ================== Other cases ==================

  // If we are at the first half of the GSW, we are adding new things to the
  // first polynomial (c0) of the given BFV ciphertext. c1 is not touched.
  seal::Ciphertext cipher;
  if (half == 0) {
    encryptor.encrypt_zero_symmetric_seeded(cipher);
  } else {
    encryptor.encrypt_zero_symmetric(cipher);
  }
  auto ct = cipher.data(half);
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
      uint128_t val = (uint128_t)pt[j] * gadget_coef % mod;
      ct[j + pad] =
          static_cast<uint64_t>((ct[j + pad] + val) % mod);
    }
  }
  return cipher;
}

void GSWEval::sealGSWVecToGSW(GSWCiphertext &output, const std::vector<seal::Ciphertext> &gsw_vec) {
  const auto &context_data = context->first_context_data();
  auto &parms = context_data->parms();
  auto &coeff_modulus = parms.coeff_modulus();
  size_t coeff_count = parms.poly_modulus_degree();
  size_t rns_mod_cnt = coeff_modulus.size();

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