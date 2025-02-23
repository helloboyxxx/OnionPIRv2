#pragma once

#include "gsw_eval.h"
#include "pir.h"
#include <optional>

#define RAW_DB_FILE "./rawDB.bin"

// typedef std::vector<std::optional<seal::Plaintext>> DatabaseChunk;  // 256 plaintexts
typedef std::unique_ptr<std::optional<seal::Plaintext>[]> DatabaseChunk;  // Heap allocation for N_1 plaintexts
typedef std::unique_ptr<std::optional<seal::Plaintext>[]> Database;       // One consecutive huge vector for the entire database on heap

class PirServer {
public:
  PirServer(const PirParams &pir_params);
  ~PirServer();

  /**
   * Generate random data for the server database and directly set the database.
   * It pushes the data to the database in chunks.
   */
  void gen_data();

  // push one chunk of entry to the given database
  void push_database_chunk(std::vector<Entry> &chunk_entry, const size_t chunk_idx);

  // Given the client id and a packed client query, this function first unpacks the query, then returns the retrieved encrypted result.
  std::vector<seal::Ciphertext> make_query(const size_t client_id, std::stringstream &query_stream);

  /**
   * @brief A clever way to evaluate the external product for second to last dimensions. 
   * 
   * @param result The BFV ciphertexts
   * @param selection_cipher A single RGSW(b) ciphertext, where b \in {0, 1}. 0 to get the first half of the result, 1 to get the second half.
   */
  void evaluate_gsw_product(std::vector<seal::Ciphertext> &result, GSWCiphertext &selection_cipher);
  
  void set_client_galois_key(const size_t client_id, std::stringstream &gsw_stream);
  void set_client_gsw_key(const size_t client_id, std::stringstream &gsw_stream);

  /**
  Asking the server to return the entry at the given (abstract) index.
  This is not doing PIR. So this reveals the index to the server. This is
  only for testing purposes.
  */
  Entry direct_get_entry(const size_t index) const;

private:
  size_t num_pt_;
  seal::SEALContext context_;
  seal::Evaluator evaluator_;
  std::vector<size_t> dims_;
  std::map<size_t, seal::GaloisKeys> client_galois_keys_;
  std::map<size_t, GSWCiphertext> client_gsw_keys_;
  Database db_; // pointer to the entire database vector
  std::vector<uint64_t> db_aligned_; // aligned database for fast first dim
  std::vector<uint128_t> inter_res_; // pointer to the intermediate result vector for fst dim
  PirParams pir_params_;
  GSWEval key_gsw_;
  GSWEval data_gsw_;

  /*!
    Expands the first query ciphertext into a selection vector of ciphertexts
    where the ith ciphertext encodes the ith bit of the first query ciphertext.
  */
  std::vector<seal::Ciphertext> expand_query(size_t client_id, seal::Ciphertext &ciphertext) const;
  /*!
    Performs a cross product between the first selection vector and the
    database.
  */
  std::vector<seal::Ciphertext> evaluate_first_dim(std::vector<seal::Ciphertext> &selection_vector);
  /*!
    Transforms the plaintexts in the database into their NTT representation.
    This speeds up computation but takes up more memory.
  */
  void preprocess_ntt();

  void realign_db();

  // Fill the intermediate_db_ with some ciphertext. We just need to allocate the memory.
  void fill_inter_res();

  // write one chunk of the database to a binary file in CACHE_DIR
  void write_one_chunk(std::vector<Entry> &chunk);
};
