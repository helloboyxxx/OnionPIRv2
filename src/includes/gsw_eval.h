#pragma once
#include "seal/seal.h"
#include "pir.h"
#include <vector>


// A GSWCiphertext is a flattened 2lx2 matrix of polynomials
typedef std::vector<std::vector<uint64_t>> GSWCiphertext;

class GSWEval {
  private:
    PirParams pir_params_;
    size_t l_;
    size_t base_log2_;
  
  public:
    GSWEval(const PirParams &pir_params, const size_t l, const size_t base_log2)
        : pir_params_(pir_params), l_(l), base_log2_(base_log2) {}
    ~GSWEval() = default;
    GSWEval(const GSWEval &gsw_eval) = default;

    /*!
      Computes the external product between a GSW ciphertext and a decomposed BFV
      ciphertext.
      @param gsw_enc -GSW Ciphertext, should only encrypt 0 or 1 to prevent large
      noise growth
      @param rlwe_expansion - decomposed vector of BFV ciphertext
      @param ct_poly_size - number of ciphertext polynomials
      @param res_ct - output ciphertext
    */

    void external_product(GSWCiphertext const &gsw_enc, seal::Ciphertext const &bfv,
                          seal::Ciphertext &res_ct);

    /*!
      Performs a gadget decomposition of a size 2 BFV ciphertext into 2 sets of
      rows of l polynomials (the 2 sets are concatenated into a single vector of
      vectors). Each polynomial coefficient encodes the value congruent to the
      original ciphertext coefficient modulus the value of base^(l-row).
      @param ct - input BFV ciphertext. Should be of size 2.
      @param output - output to store the decomposed ciphertext as a vector of
      vectors of polynomial coefficients
    */
    void decomp_rlwe(seal::Ciphertext const &ct, std::vector<std::vector<uint64_t>> &output);


    // Similar to decomp_rlwe. Use this when rn_mod_cnt = 1. It uses faster right shift.
    void decomp_rlwe_single_mod(seal::Ciphertext const &ct, std::vector<std::vector<uint64_t>> &output);

    /*!
      Generates a GSW ciphertext from a BFV ciphertext query.

      @param query - input BFV ciphertext. Should be of size l * 2.
      @param gsw_key - GSW encryption of -s
      @param output - output to store the GSW ciphertext as a vector of vectors of
      polynomial coefficients
    */
    void query_to_gsw(std::vector<seal::Ciphertext> query, GSWCiphertext gsw_key,
                      GSWCiphertext &output);


    // The input plaintext is a vector of all the coefficients. Assert that the size of this vector 
    // is equal to the number of coefficients in the Plaintext, including the zero coefficients.
    void encrypt_plain_to_gsw(std::vector<uint64_t> const &plaintext,
                              seal::Encryptor const &encryptor,
                              seal::SecretKey const &sk,
                              std::vector<seal::Ciphertext> &output);

    /**
    @brief Helper function for encrypt_plain_to_gsw. This function encrypts the
    plaintext to a single row of the GSW ciphertext at the given "half" and the
    given l. half = 0 means the first half of the GSW.

    @param plaintext
    @param half 0 denotes the top l rows, 1 denotes the bottom l rows
    @param level level in the given half
    @return seal::Ciphertext
    */
    seal::Ciphertext
    enc_plain_to_gsw_one_row(std::vector<uint64_t> const &plaintext,
                            seal::Encryptor const &encryptor,
                            seal::SecretKey const &sk, const size_t half,
                            const size_t level);

    /**
    * @brief Transform the given GSWCipher text from polynomial representation to NTT representation.
    * 
    */
    void gsw_ntt_negacyclic_harvey(GSWCiphertext &gsw);

    void ciphertext_inverse_ntt(seal::Ciphertext &ct);

    // helper functions
    void sealGSWVecToGSW(GSWCiphertext &output, const std::vector<seal::Ciphertext> &gsw_vec);
};