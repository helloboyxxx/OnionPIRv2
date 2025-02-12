#pragma once

#include "pir.h"
class PirClient {
public:
  PirClient(const PirParams &pirparms);
  ~PirClient();

  /**
  This is the core function for the client.
  High level steps:
  1. Compute the query indices.
  2. Creates a plain_query (pt in paper), add the first dimension, then encrypts it.
  3. For the rest dimensions, calculate required RGSW coefficients and insert
  them into the ciphertext. Result is $\tilde c$ in paper.
  @param entry_index The input to the PIR blackbox.
  @return returns a seal::Ciphertext with a a seed stored in
  c_1, which should not be touched before doing serialization.
  */
  PirQuery generate_query(const std::uint64_t entry_index);
  size_t write_query_to_stream(const PirQuery &query, std::stringstream &data_stream);
  size_t write_gsw_to_stream(const std::vector<Ciphertext> &gsw, std::stringstream &gsw_stream);
  size_t create_galois_keys(std::stringstream &galois_key_stream) const;
  std::vector<seal::Plaintext> decrypt_result(std::vector<seal::Ciphertext> reply);
  seal::Decryptor *get_decryptor();
  // Retrieves an entry from the plaintext containing the entry.
  std::vector<Ciphertext> generate_gsw_from_key();
  // The height of the expansion tree during packing unpacking stages
  inline const size_t get_expan_height() const {
    return std::ceil(std::log2(dims_[0] + pir_params_.get_l() * (dims_.size() - 1)));
  }
  
  Entry get_entry_from_plaintext(size_t entry_index, seal::Plaintext plaintext);
  uint32_t get_client_id() const;

private:
  uint32_t client_id_;
  seal::EncryptionParameters enc_params_;
  PirParams pir_params_;
  std::vector<size_t> dims_;
  seal::KeyGenerator *keygen_;
  seal::Decryptor *decryptor_;
  seal::Encryptor *encryptor_;
  seal::Evaluator *evaluator_;
  seal::SEALContext *context_;
  const seal::SecretKey *secret_key_;
  
  // Gets the corresponding plaintext index in a database for a given entry index
  size_t get_database_plain_index(size_t entry_index);

  // Gets the query indices for a given plaintext
  std::vector<size_t> get_query_indices(size_t plaintext_index);
};
