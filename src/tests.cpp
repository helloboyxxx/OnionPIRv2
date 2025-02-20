#include "tests.h"
#include "gsw_eval.h"
#include "pir.h"
#include "server.h"
#include "client.h"
#include "utils.h"
#include <cassert>
#include <iostream>



#define EXPERIMENT_ITERATIONS 3

void print_func_name(std::string func_name) {
  PRINT_BAR;
  #ifdef _DEBUG
    std::cout << "                    "<< func_name << "(Debug build)" << std::endl;
  #endif
  #ifdef _BENCHMARK
    std::cout << "                    "<< func_name << "(Benchmark build)" << std::endl;
  #endif
  PRINT_BAR;
}

void run_tests() {
  test_pir();
  bfv_example();
}


void test_pir() {
  print_func_name(__FUNCTION__);
  auto success_count = 0;
  
  // ============== setting parameters for PIR scheme ==============
  PirParams pir_params;
  pir_params.print_params();
  PirServer server(pir_params); // Initialize the server with the parameters
  
  BENCH_PRINT("Initializing server...");
  // Data to be stored in the database.
  server.gen_data();
  BENCH_PRINT("Server initialized");

  // some global results
  size_t galois_key_size = 0;
  size_t gsw_key_size = 0;
  size_t query_size = 0;
  size_t response_size = 0;

  // Run the query process many times.
  srand(time(0)); // reset the seed for the random number generator
  for (size_t i = 0; i < EXPERIMENT_ITERATIONS; i++) {
    
    // ============= OFFLINE PHASE ==============
    // Initialize the client
    PirClient client(pir_params);
    const size_t client_id = client.get_client_id();
    std::stringstream galois_key_stream, gsw_stream, data_stream;

    // Client create galois keys and gsw keys and writes to the stream (to the server)
    galois_key_size = client.create_galois_keys(galois_key_stream);
    gsw_key_size = client.write_gsw_to_stream(
        client.generate_gsw_from_key(), gsw_stream);
    //--------------------------------------------------------------------------------
    // Server receives the gsw keys and galois keys and loads them when needed
    server.set_client_galois_key(client_id, galois_key_stream);
    server.set_client_gsw_key(client_id, gsw_stream);

    // ===================== ONLINE PHASE =====================
    // Client start generating query
    size_t entry_index = rand() % pir_params.get_num_entries();
  
    // ============= CLIENT ===============
    TIME_START(CLIENT_TOT_TIME);
    PirQuery query = client.generate_query(entry_index);
    query_size = client.write_query_to_stream(query, data_stream);
    TIME_END(CLIENT_TOT_TIME);
    
    // ============= SERVER ===============
    TIME_START(SERVER_TOT_TIME);
    auto result = server.make_query(client_id, data_stream);
    TIME_END(SERVER_TOT_TIME);

    // ============= CLIENT ===============
    TIME_START(CLIENT_TOT_TIME);
    // client gets result from the server and decrypts it
    auto decrypted_result = client.decrypt_result(result);
    Entry result_entry = client.get_entry_from_plaintext(entry_index, decrypted_result[0]);
    TIME_END(CLIENT_TOT_TIME);

    // write the result to the stream to test the size
    std::stringstream result_stream;
    response_size = result[0].save(result_stream);
    result_stream.str(std::string()); // clear the stream

    // Directly get the plaintext from server. Not part of PIR.
    Entry actual_entry = server.direct_get_entry(entry_index);
    // extract and print the actual entry index
    uint64_t actual_entry_idx = get_entry_idx(actual_entry);
    uint64_t result_entry_idx = get_entry_idx(result_entry);
    
    END_EXPERIMENT();
    // ============= PRINTING RESULTS ===============    
    DEBUG_PRINT("\t\tWanted/result/actual idx:\t" << entry_index << " / " << result_entry_idx << " / " << actual_entry_idx);
    #ifdef _DEBUG
    PRINT_RESULTS(i+1);
    #endif

    if (entry_is_equal(result_entry, actual_entry)) {
      // print a green success message
      std::cout << "\033[1;32mSuccess!\033[0m" << std::endl;
      success_count++;
    } else {
      // print a red failure message
      std::cout << "\033[1;31mFailure!\033[0m" << std::endl;
      std::cout << "PIR Result:\t";
      print_entry(result_entry);
      std::cout << "Actual Entry:\t";
      print_entry(actual_entry);
    }
    PRINT_BAR;
  }

  double avg_server_time = GET_AVG_TIME(SERVER_TOT_TIME);
  double throughput = pir_params.get_DBSize_MB() / (avg_server_time / 1000);
  

  // ============= PRINTING FINAL RESULTS ===============]
  BENCH_PRINT("                       FINAL RESULTS                         ")
  PRINT_BAR;
  BENCH_PRINT("Success rate: " << success_count << "/" << EXPERIMENT_ITERATIONS);
  BENCH_PRINT("galois key size: " << galois_key_size << " bytes");
  BENCH_PRINT("gsw key size: " << gsw_key_size << " bytes");
  BENCH_PRINT("total key size: " << static_cast<double>(galois_key_size + gsw_key_size) / 1024 / 1024 << "MB");
  BENCH_PRINT("query size: " << query_size << " bytes");
  BENCH_PRINT("response size: " << response_size << " bytes");
  
  // PRINT_AVERAGE_RESULTS();
  PRETTY_PRINT();
  BENCH_PRINT("Throughput: " << throughput << " MB/s");
}

// TODO: I deleted some small test cases and examples because I am refactoring the code. I will add them back later.


  // This is an example of how to use the BFV scheme in SEAL and in our PIR scheme.
 void bfv_example() {
  print_func_name(__FUNCTION__);

  // You need a a chunk of code to init the seal parameters. Here is the minimum you need:
  seal::EncryptionParameters params(seal::scheme_type::bfv);
  const size_t coeff_count = 4096;  // you can try other powers of two.
  params.set_poly_modulus_degree(coeff_count); // example: a_1 x^4095 + a_2 x^4094 + ...
  const uint64_t pt_mod = generate_prime(49); // 49 bits for the plain modulus, then you can use 48 bits for storing data.
  params.set_plain_modulus(pt_mod);
  std::vector<int> bit_sizes({60, 60,60}); // You can also try our own DatabaseConstants::CoeffMods
  const auto coeff_modulus = CoeffModulus::Create(coeff_count, bit_sizes);
  params.set_coeff_modulus(coeff_modulus);
  // ================== END OF SEAL PARAMS INIT ==================
  // The following are things you need to encrypt, evaluate, and decrypt BFV.
  SEALContext context_(params);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);
  // =============================================================
  BENCH_PRINT("coeff_count: " << coeff_count);
  BENCH_PRINT("Num of coeff mods that SEAL uses: "
              << context_.key_context_data()->parms().coeff_modulus().size());
  BENCH_PRINT("Num of coeff mods used for actual ciphertexts"
              << context_.first_context_data()->parms().coeff_modulus().size());

  // ============= Now let's try some BFV * BFV multiplication in coefficient form ==============
  seal::Plaintext a(coeff_count), b(coeff_count), result;
  a[0] = 1; a[1] = 9;
  b[0] = 3; b[1] = 6;
  BENCH_PRINT("Vector a: " << a.to_string());
  BENCH_PRINT("Vector b: " << b.to_string());

  seal::Ciphertext a_encrypted, b_encrypted, cipher_result;
  encryptor_->encrypt_symmetric(a, a_encrypted);
  encryptor_->encrypt_symmetric(b, b_encrypted);
  
  BENCH_PRINT("Noise budget before: " << decryptor_->invariant_noise_budget(a_encrypted));
  evaluator_.multiply(a_encrypted, b_encrypted, cipher_result);
  decryptor_->decrypt(cipher_result, result);
  // You can see that this direct multiplication consumes a lot of noise budget.
  BENCH_PRINT("Noise budget after: " << decryptor_->invariant_noise_budget(cipher_result));
  BENCH_PRINT("BFV x BFV result: " << result.to_string());
  PRINT_BAR;
  // ============= Now let's try addition in coefficient form ==============
  a.set_zero(); b.set_zero();
  a[0] = 1; a[1] = 9;
  b[0] = 3; b[1] = 6;
  BENCH_PRINT("Vector a: " << a.to_string());
  BENCH_PRINT("Vector b: " << b.to_string());

  encryptor_->encrypt_symmetric(a, a_encrypted);
  encryptor_->encrypt_symmetric(b, b_encrypted);
  BENCH_PRINT("Noise budget before: " << decryptor_->invariant_noise_budget(a_encrypted));
  evaluator_.add(a_encrypted, b_encrypted, cipher_result);
  decryptor_->decrypt(cipher_result, result);
  BENCH_PRINT("Noise budget after: " << decryptor_->invariant_noise_budget(cipher_result));
  BENCH_PRINT("BFV + BFV result: " << result.to_string());
  PRINT_BAR;

  // ============= Now let's try addition and multiplication in ntt form ==============
  a.set_zero(); b.set_zero();
  a[0] = 1; a[1] = 9;
  b[0] = 3; b[1] = 6;
  BENCH_PRINT("Vector a: " << a.to_string());
  BENCH_PRINT("Vector b: " << b.to_string());
  encryptor_->encrypt_symmetric(a, a_encrypted);
  encryptor_->encrypt_symmetric(b, b_encrypted);
  BENCH_PRINT("Noise budget before: " << decryptor_->invariant_noise_budget(a_encrypted));

  evaluator_.transform_to_ntt_inplace(a_encrypted);
  evaluator_.transform_to_ntt_inplace(b_encrypted);
  evaluator_.add(a_encrypted, b_encrypted, cipher_result);
  evaluator_.transform_from_ntt_inplace(cipher_result);
  
  decryptor_->decrypt(cipher_result, result);
  BENCH_PRINT("Noise budget after: " << decryptor_->invariant_noise_budget(cipher_result)); // noise budget is almost the same.
  BENCH_PRINT("NTT + NTT result: " << result.to_string());  // and the result is correct! NTT form polynomial is additive

  // ============= Now let's try BFV multiplied by a constant in ntt form ==============
  seal::Plaintext scalar(coeff_count);
  scalar[0] = 2; scalar[1] = 3;
  evaluator_.transform_to_ntt_inplace(scalar, context_.first_parms_id()); // This happens in preprocess_ntt
  // Now instead of using multiply_plain, I want to demonstrate what happens in the first dimension evaluation. 
  // This is demonstrating what you can do in ntt form, but the actual order of computation in OnionPIRv2 fst dim can be different.
  size_t rns_mod_cnt = coeff_modulus.size() - 1;
  std::vector<uint128_t> res(coeff_count * rns_mod_cnt);
  std::fill(res.begin(), res.end(), 0);
  uint64_t *ct0_ptr = a_encrypted.data(0);
  uint64_t *ct1_ptr = a_encrypted.data(1);
  uint128_t *res0_ptr = res.data();
  uint128_t *res1_ptr = res.data() + coeff_count;
  uint64_t *pt_ptr = scalar.data();
  // element wise vector multiplication.
  for (size_t i = 0; i < coeff_count * rns_mod_cnt; i++) {
    res0_ptr[i] = static_cast<uint128_t>(ct0_ptr[i]) * pt_ptr[i];
    res1_ptr[i] = static_cast<uint128_t>(ct1_ptr[i]) * pt_ptr[i];
  }
  // Another scan on the res to reduce the modulus.
  // Meanwhile we can reconstruct the ciphertext from the res vector and decrypt it.
  seal::Ciphertext scalar_mul_result = a_encrypted; // just copy a random ciphertext with correct format, we will overwrite it.
  uint64_t *scal_mul_ct0_ptr = scalar_mul_result.data(0);
  uint64_t *scal_mul_ct1_ptr = scalar_mul_result.data(1);
  for (size_t i = 0; i < coeff_count; i++) {
    for (size_t j = 0; j < rns_mod_cnt; j++) {
      auto curr_mod = coeff_modulus[j].value();
      scal_mul_ct0_ptr[i] = res0_ptr[i + j * coeff_count] % curr_mod;
      scal_mul_ct0_ptr[i] = res1_ptr[i + j * coeff_count] % curr_mod;
    }
  }

  evaluator_.transform_from_ntt_inplace(scalar_mul_result);
  decryptor_->decrypt(scalar_mul_result, result);
  BENCH_PRINT("NTT x scalar result: " << result.to_string());  // and the result is correct! NTT form polynomial is multiplicative
  PRINT_BAR;
  // ============= Now let's try subtraction in coefficient form ==============
  // This is actually interesting
}