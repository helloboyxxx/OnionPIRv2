#include "tests.h"
#include "external_prod.h"
#include "pir.h"
#include "server.h"
#include "client.h"
#include "utils.h"
#include <cassert>
#include <iostream>
#include <bitset>



#define EXPERIMENT_ITERATIONS 5
#define WARMUP_ITERATIONS     3

void print_func_name(std::string func_name) {
#ifdef _DEBUG
  std::cout << "                    "<< func_name << "(Debug build)" << std::endl;
#endif
#ifdef _BENCHMARK
  std::cout << "                    "<< func_name << "(Benchmark build)" << std::endl;
#endif
}

void run_tests() {
  DEBUG_PRINT("Running tests");
  PRINT_BAR;

  // If we compare the following two examples, we do see that external product increase the noise much slower than BFV x BFV.
  // bfv_example();
  // test_external_product();
  // test_ntt_add();
  // test_ntt_scalar_mul();
  // test_ct_sub();
  // serialization_example();
  // test_plain_to_gsw();

  test_pir();

  // test_prime_gen();
  // test_reading_pt();
  // test_reading_ct();
  // test_matrix_mult();

  PRINT_BAR;
  DEBUG_PRINT("Tests finished");
}

/**
 * @brief This is a BFV x BFV example. The coefficients in example vectors and the result are in hex.
 */
void bfv_example() {
  print_func_name(__FUNCTION__);

  EncryptionParameters parms(scheme_type::bfv);
  size_t poly_degree = 4096;
  parms.set_poly_modulus_degree(poly_degree);
  parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_degree));

  SEALContext context_(parms);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);
  DEBUG_PRINT("poly_degree: " << poly_degree);
  std::cout << "Size f: " << context_.key_context_data()->parms().coeff_modulus().size()
            << std::endl;
  std::cout << "Size f: " << context_.first_context_data()->parms().coeff_modulus().size()
            << std::endl;
  seal::Plaintext a(poly_degree), b(poly_degree), result;
  a[0] = 1;
  a[1] = 9;

  b[0] = 3;
  b[1] = 6;

  DEBUG_PRINT("Vector a: " << a.to_string());
  DEBUG_PRINT("Vector b: " << b.to_string());

  seal::Ciphertext a_encrypted, b_encrypted, cipher_result;
  encryptor_->encrypt_symmetric(a, a_encrypted);
  encryptor_->encrypt_symmetric(b, b_encrypted);
  
  std::cout << "Noise budget before: " << decryptor_->invariant_noise_budget(a_encrypted)
            << std::endl;

  evaluator_.multiply(a_encrypted, b_encrypted, cipher_result);
  decryptor_->decrypt(cipher_result, result);
  std::cout << "Noise budget after: " << decryptor_->invariant_noise_budget(cipher_result) << std::endl;
  std::cout << "BFV x BFV result: " << result.to_string() << std::endl;
}

// This is a BFV x GSW example
void test_external_product() {
  print_func_name(__FUNCTION__);
  PirParams pir_params;
  pir_params.print_params();
  auto parms = pir_params.get_seal_params();    // This parameter is set to be: seal::scheme_type::bfv
  auto context_ = seal::SEALContext(parms);   // Then this context_ knows that it is using BFV scheme
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = seal::Encryptor(context_, secret_key_);
  auto decryptor_ = seal::Decryptor(context_, secret_key_);
  size_t coeff_count = parms.poly_modulus_degree();
  uint64_t poly_degree = pir_params.get_seal_params().poly_modulus_degree();

  DEBUG_PRINT("poly_degree: " << poly_degree);
  // the test data vector a and results are both in BFV scheme.
  seal::Plaintext a(poly_degree), result;
  size_t plain_coeff_count = a.coeff_count();
  seal::Ciphertext a_encrypted(context_), cipher_result(context_);    // encrypted "a" will be stored here.
  auto &context_data = *context_.first_context_data();

  // vector b
  std::vector<uint64_t> b(poly_degree);

  // vector a is in the context of BFV scheme. 
  // 0, 1, 2, 4 are coeff_index of the term x^i, 
  // the index of the coefficient in the plaintext polynomial
  a[0] = 1;
  a[1] = 2;
  a[2] = 3;

  DEBUG_PRINT("Vector a: " << a.to_string());

  // vector b is in the context of GSW scheme.
  // b[0] = 3;
  b[2] = 5;
  
  // print b
  std::string b_result = "Vector b: ";
  for (int i = 0; i < 5; i++) {
    b_result += std::to_string(b[i]) + " ";
  }
  DEBUG_PRINT(b_result);
  
  // Since a_encrypted is in a context of BFV scheme, the following function encrypts "a" using BFV scheme.
  encryptor_.encrypt_symmetric(a, a_encrypted);

  std::cout << "Noise budget before: " << decryptor_.invariant_noise_budget(a_encrypted)
            << std::endl;
  GSWCiphertext b_gsw;
  std::vector<seal::Ciphertext> temp_gsw;
  data_gsw.encrypt_plain_to_gsw(b, encryptor_, secret_key_, temp_gsw);
  data_gsw.sealGSWVecToGSW(b_gsw, temp_gsw);
  data_gsw.gsw_ntt_negacyclic_harvey(b_gsw);  // transform b_gsw to NTT form

  size_t mult_rounds = 1;

  for (int i = 0; i < mult_rounds; i++) {
    data_gsw.external_product(b_gsw, a_encrypted, a_encrypted);
    a_encrypted.is_ntt_form() = true;
    evaluator_.transform_from_ntt_inplace(a_encrypted);
    decryptor_.decrypt(a_encrypted, result);
    std::cout << "Noise budget after: " << decryptor_.invariant_noise_budget(a_encrypted)
              << std::endl;
  
  // output decrypted result
  std::cout << "External product result: " << result.to_string() << std::endl;
  }
}

// This is a simple test for demonstrating that we are expecting a "transparent"
// ciphertext if we have two identical ciphertext(value equal) subtracted. This
// is because it is possible to have two entries in the database that are thesame. 
// Please use -DSEAL_THROW_ON_TRANSPARENT_CIPHERTEXT=OFF flag when compiling SEAL.
void test_ct_sub() {
  print_func_name(__FUNCTION__);
  PirParams pir_params;
  auto parms = pir_params.get_seal_params();    // This parameter is set to be: seal::scheme_type::bfv
  auto context_ = seal::SEALContext(parms);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);

  // Create a ciphertext of 1
  seal::Plaintext c(pir_params.get_seal_params().poly_modulus_degree());
  c[0] = 2;
  seal::Ciphertext c_encrypted(context_);
  encryptor_->encrypt_symmetric(c, c_encrypted);

  // Create two plaintexts of 3. This mimics the two identical entries in the database
  seal::Plaintext pt1(pir_params.get_seal_params().poly_modulus_degree());
  seal::Plaintext pt2(pir_params.get_seal_params().poly_modulus_degree());
  pt1[0] = 3; 
  pt2[0] = 3;

  // Multiplication of a * pt1 and a * pt2
  seal::Ciphertext result_1(context_);
  seal::Ciphertext result_2(context_);
  evaluator_.multiply_plain(c_encrypted, pt1, result_1);
  evaluator_.multiply_plain(c_encrypted, pt2, result_2);

  // Subtraction
  evaluator_.sub_inplace(result_1, result_2);

  // Decrypt the result
  seal::Plaintext result_pt;
  decryptor_->decrypt(result_1, result_pt);
  std::cout << "Result: " << result_pt.to_string() << std::endl;
}

void test_ntt_add() {
  print_func_name(__FUNCTION__);
  PirParams pir_params;
  auto parms = pir_params.get_seal_params();    // This parameter is set to be: seal::scheme_type::bfv
  auto context_ = seal::SEALContext(parms);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);

  // Create two ciphertexts
  seal::Plaintext c1(pir_params.get_seal_params().poly_modulus_degree());
  seal::Plaintext c2(pir_params.get_seal_params().poly_modulus_degree());
  c1[0] = 4;
  c1[1] = 23;
  c2[0] = 2;
  c2[4] = 12;
  seal::Ciphertext c1_encrypted(context_);
  seal::Ciphertext c2_encrypted(context_);
  encryptor_->encrypt_symmetric(c1, c1_encrypted);
  encryptor_->encrypt_symmetric(c2, c2_encrypted);

  // transform both of them to NTT form
  evaluator_.transform_to_ntt_inplace(c1_encrypted);
  evaluator_.transform_to_ntt_inplace(c2_encrypted);


  // Add the two ciphertexts and store them in a new ciphertext
  seal::Ciphertext result_1(context_);
  evaluator_.add(c1_encrypted, c2_encrypted, result_1);

  // transform the result back to normal form
  evaluator_.transform_from_ntt_inplace(result_1);

  // Decrypt the result
  seal::Plaintext result_pt;
  decryptor_->decrypt(result_1, result_pt);
  std::cout << "Result: " << result_pt.to_string() << std::endl;
}

// create a ciphtertext and transform it to ntt form. Then multiply it with a scalar.
void test_ntt_scalar_mul() {
  print_func_name(__FUNCTION__);
  PirParams pir_params;
  auto parms = pir_params.get_seal_params();    // This parameter is set to be: seal::scheme_type::bfv
  auto context_ = seal::SEALContext(parms);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);

  // Create a ciphertext
  seal::Plaintext c(pir_params.get_seal_params().poly_modulus_degree());
  c[1] = 4;
  c[0] = 3;
  seal::Ciphertext c_encrypted(context_);
  encryptor_->encrypt_symmetric(c, c_encrypted);

  // transform it to NTT form
  // evaluator_.transform_to_ntt_inplace(c_encrypted);

  // we use seal::util::left_shift_uint to do the scalar multiplication
  seal::util::left_shift_uint(c_encrypted.data(0), 2, 1, c_encrypted.data(0));
  seal::util::left_shift_uint(c_encrypted.data(1), 2, 1, c_encrypted.data(1));

  // transform the result back to normal form
  // c_encrypted.is_ntt_form() = true;
  // evaluator_.transform_from_ntt_inplace(c_encrypted);

  // Decrypt the result
  seal::Plaintext result_pt;
  decryptor_->decrypt(c_encrypted, result_pt);
  std::cout << "Result: " << result_pt.to_string() << std::endl;
}


void serialization_example() {
  PirParams pir_params;
  const auto params = pir_params.get_seal_params();
  const auto context_ = seal::SEALContext(params);
  const auto evaluator_ = seal::Evaluator(context_);
  const auto keygen_ = seal::KeyGenerator(context_);
  const auto secret_key_ = keygen_.secret_key();
  const auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  const auto decryptor_ = new seal::Decryptor(context_, secret_key_);

  std::stringstream data_stream;

  // ================== Raw Zero ciphertext ==================
  seal::Ciphertext raw_zero;
  encryptor_->encrypt_zero_symmetric(raw_zero);
  auto raw_size = raw_zero.save(data_stream); // store the raw zero in the stream

  // ================== SEAL original method for creating serialized zero ==================
  // Original method for creating a serializable object
  Serializable<Ciphertext> orig_serialized_zero = encryptor_->encrypt_zero_symmetric();
  auto s_size = orig_serialized_zero.save(data_stream);   // ! Storing the original zero

  // ================== New way to create a ciphertext with a seed ==================
  // New way to create a ciphertext with a seed, do some operations and then convert it to a serializable object.
  seal::Ciphertext new_seeded_zero;
  encryptor_->encrypt_zero_symmetric_seeded(new_seeded_zero); // This function allows us to change the ciphertext.data(0).

  // Add something in the third coeeficient of seeded_zero
  DEBUG_PRINT("Size: " << new_seeded_zero.size());
  auto ptr_0 = new_seeded_zero.data(0);
  auto ptr_1 = new_seeded_zero.data(1); // corresponds to the second polynomial (c_1)
  // print the binary value of the first coefficient
  BENCH_PRINT("Indicator:\t" << std::bitset<64>(ptr_1[0]));  // used in has_seed_marker()
  // the seed is stored in here. By the time I write this code, it takes 81
  // bytes to store the prng seed. Notice that they have common headers.
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[1]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[2]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[3]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[4]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[5]));
  
  auto mods = context_.first_context_data()->parms().coeff_modulus();
  auto plain_modulus = params.plain_modulus().value();
  uint128_t mod_0 = mods[0].value();
  uint128_t mod_1 = mods[1].value();
  uint128_t delta = mod_0 * mod_1 / plain_modulus;
  uint128_t message = 15;
  uint128_t to_add = delta * message;
  auto padding = params.poly_modulus_degree();
  ptr_0[0] = (ptr_0[0] + (to_add % mod_0)) % mod_0;
  ptr_0[0 + padding] = (ptr_0[0 + padding] + (to_add % mod_1)) % mod_1;

  // write the serializable object to the stream
  auto s2_size = new_seeded_zero.save(data_stream); // ! Storing new ciphertext with a seed

  // ================== Deserialize and decrypt the ciphertexts ==================
  seal::Ciphertext raw_ct, orig_ct, new_ct;
  raw_ct.load(context_, data_stream);  // ! loading the raw zero
  orig_ct.load(context_, data_stream);  // ! loading the original zero
  new_ct.load(context_, data_stream); // ! loading the new ciphertext with a seed 

  // decrypt the ciphertexts
  seal::Plaintext raw_pt, orig_pt, new_pt;
  decryptor_->decrypt(raw_ct, raw_pt);
  decryptor_->decrypt(orig_ct, orig_pt);
  decryptor_->decrypt(new_ct, new_pt);

  // ================== Print the results ==================
  BENCH_PRINT("Raw zero size: " << raw_size);
  BENCH_PRINT("Serializable size 1: " << s_size);
  BENCH_PRINT("Serializable size 2: " << s2_size); // smaller size, but allow us to work on the ciphertext!

  BENCH_PRINT("Raw plaintext: " << raw_pt.to_string());
  BENCH_PRINT("Original plaintext: " << orig_pt.to_string());
  BENCH_PRINT("New plaintext: " << new_pt.to_string()); // Hopefully, this decrypts to the message.
}


void test_pir() {
  print_func_name(__FUNCTION__);
  auto server_time_sum = 0;
  auto client_time_sum = 0;
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
  for (int i = 0; i < EXPERIMENT_ITERATIONS + WARMUP_ITERATIONS; i++) {
    
    // ============= OFFLINE PHASE ==============
    // Initialize the client
    PirClient client(pir_params);
    const int client_id = rand();
    std::stringstream galois_key_stream, gsw_stream, data_stream;

    // Client create galois keys and gsw keys and writes to the stream (to the server)
    galois_key_size = client.create_galois_keys(galois_key_stream);
    gsw_key_size = client.write_gsw_to_stream(
        client.generate_gsw_from_key(), gsw_stream);
    //--------------------------------------------------------------------------------
    server.decryptor_ = client.get_decryptor();
    // Server receives the gsw keys and galois keys and loads them when needed
    server.set_client_galois_key(client_id, galois_key_stream);
    server.set_client_gsw_key(client_id, gsw_stream);

    // ===================== ONLINE PHASE =====================
    // Client start generating query
    size_t entry_index = rand() % pir_params.get_num_entries();
    BENCH_PRINT("Experiment [" << i+1 << "]");
    BENCH_PRINT("\t\tEntry index:\t" << entry_index);

    // ============= CLIENT ===============
    auto c_start_time = CURR_TIME;  // client start time for the query
    PirQuery query = client.generate_query(entry_index);
    query_size = client.write_query_to_stream(query, data_stream);
    
    // ============= SERVER ===============
    auto s_start_time = CURR_TIME;  // server start time for processing the query
    auto result = server.make_seeded_query(client_id, data_stream);
    auto s_end_time = CURR_TIME;


    // ============= CLIENT ===============
    // client gets result from the server and decrypts it
    auto decrypted_result = client.decrypt_result(result);
    Entry entry = client.get_entry_from_plaintext(entry_index, decrypted_result[0]);
    auto c_end_time = CURR_TIME;

    // write the result to the stream to test the size
    std::stringstream result_stream;
    response_size = result[0].save(result_stream);
    result_stream.str(std::string()); // clear the stream

    // Directly get the plaintext from server. Not part of PIR.
    Entry actual_entry = server.direct_get_entry(entry_index);

    // extract and print the actual entry index
    uint64_t actual_entry_idx = get_entry_idx(actual_entry);

    // ============= PRINTING RESULTS ===============    
    BENCH_PRINT("\t\tActual Index:\t" << actual_entry_idx);
    BENCH_PRINT("\t\tServer time:\t" << TIME_DIFF(s_start_time, s_end_time) << " ms");
    BENCH_PRINT("\t\tClient Time:\t" << TIME_DIFF(c_start_time, c_end_time) - TIME_DIFF(s_start_time, s_end_time) << " ms"); 
    BENCH_PRINT("\t\tNoise budget:\t" << client.get_decryptor()->invariant_noise_budget(result[0]));
    
    if (i < WARMUP_ITERATIONS) {
      PRINT_BAR;
      continue;
    }
    server_time_sum += TIME_DIFF(s_start_time, s_end_time);
    client_time_sum += TIME_DIFF(c_start_time, c_end_time) - TIME_DIFF(s_start_time, s_end_time);
    if (entry_is_equal(entry, actual_entry)) {
      // print a green success message
      std::cout << "\033[1;32mSuccess!\033[0m" << std::endl;
      success_count++;
    } else {
      // print a red failure message
      std::cout << "\033[1;31mFailure!\033[0m" << std::endl;

      std::cout << "PIR Result:\t";
      print_entry(entry);
      std::cout << "Actual Entry:\t";
      print_entry(actual_entry);
    }
    PRINT_BAR;
  }

  auto avg_server_time = server_time_sum / EXPERIMENT_ITERATIONS;
  BENCH_PRINT("Average server time: " << avg_server_time << " ms");
  BENCH_PRINT("Average client time: " << client_time_sum / EXPERIMENT_ITERATIONS << " ms");
  BENCH_PRINT("galois key size: " << galois_key_size << " bytes");
  BENCH_PRINT("gsw key size: " << gsw_key_size << " bytes");
  BENCH_PRINT("total key size: " << static_cast<double>(galois_key_size + gsw_key_size) / 1024 / 1024 << "MB");
  BENCH_PRINT("query size: " << query_size << " bytes");
  BENCH_PRINT("response size: " << response_size << " bytes");
  BENCH_PRINT("Throughput: " << pir_params.get_DBSize_MB() / (static_cast<double>(avg_server_time) / 1000) << " MB/s");
  BENCH_PRINT("Success rate: " << success_count << "/" << EXPERIMENT_ITERATIONS);
}

void test_prime_gen() {
  print_func_name(__FUNCTION__);
  for (size_t i = 2; i < 65; ++i) {
    DEBUG_PRINT(generate_prime(i));
  }
}



void test_reading_pt() {
    print_func_name(__FUNCTION__);
  // This test generate a super simple plaintext database in coefficient form. Then we
  // simply try to read all the uint64_t values from the database and perform a
  // simple operation on them. The goal is to force reading the database once.
  // We want to know the cache performance.

  // For comparison, we create another vector of the same size (same number
  // of uint64_t values) and perform the same operation on them.

  // Setup PIR param for convenience
  PirParams pir_params;
  const size_t num_pt = pir_params.get_num_pt();
  const auto seal_params = pir_params.get_seal_params();
  const auto context = seal::SEALContext(seal_params);
  const size_t coeff_count = seal_params.poly_modulus_degree();
  const auto bits_per_coeff = pir_params.get_num_bits_per_coeff();
  const uint64_t coeff_mask = (uint64_t(1) << (bits_per_coeff)) - 1;  // bits_per_coeff many 1s
  const size_t storage_MB = num_pt * coeff_count * sizeof(uint64_t) / 1024 / 1024;
  auto evaluator = seal::Evaluator(context);

  std::cout << "num_pt: " << num_pt << std::endl;
  std::cout << "coeff_count: " << coeff_count << std::endl;
  std::cout << "storage size in MB: " << storage_MB << std::endl;

  // Plaintext database
  Database db = std::make_unique<std::optional<seal::Plaintext>[]>(num_pt);
  // initialize the database with reserved Plaintexts

  BENCH_PRINT("Generating the database...");
  for (size_t i = 0; i < num_pt; i++) {
    auto curr_pt = seal::Plaintext(coeff_count);
    for (size_t j = 0; j < coeff_count; j++) {
      curr_pt[j] = i * j & coeff_mask;
    }
    // transform the plaintext to NTT form
    evaluator.transform_to_ntt_inplace(curr_pt, context.first_parms_id());
    db[i] = std::move(curr_pt);
  }


  // Create a vector of the same size on the heap
  std::vector<uint64_t> vec;
  vec.reserve(num_pt * coeff_count);
  for (size_t i = 0; i < num_pt; i++) {
    for (size_t j = 0; j < coeff_count; j++) {
      vec.push_back(i * j & coeff_mask);
    }
  }

  // Read the database and perform a simple operation
  auto db_start_time = CURR_TIME;
  for (size_t i = 0; i < num_pt; i++) {
    for (size_t poly_id = 0; poly_id < 2; ++poly_id) {
      seal::Plaintext &pt = db[i].value();
      #pragma GCC unroll 32
      for (size_t j = 0; j < coeff_count; j++) {
        asm volatile("" : : "r"(pt[j]) : "memory");
      }
    }
  }
  auto db_end_time = CURR_TIME;

  // Perform the same operation on the vector
  auto vec_start_time = CURR_TIME;
  for (size_t i = 0; i < num_pt; i++) {
    for (size_t poly_id = 0; poly_id < 2; ++poly_id) {
      #pragma GCC unroll 32
      for (size_t j = 0; j < coeff_count; j++) {
        asm volatile("" : : "r"(vec[i * coeff_count + j]) : "memory");
      }
    }
  }
  auto vec_end_time = CURR_TIME;


  // simple tile
  uint64_t *p = vec.data();
  uint64_t *end = p + vec.size();
  const size_t tile_size = 4096;

  auto vec_simple_start = CURR_TIME;
  while (p < end)
  {
      #pragma GCC unroll 32
      for (size_t i = 0; i < tile_size; i++)
      {
          asm volatile("" : : "r"(p[i]) : "memory");
      }
      #pragma GCC unroll 32
      for (size_t i = 0; i < tile_size; i++)
      {
          asm volatile("" : : "r"(p[i]) : "memory");
      }
      p += tile_size;
  }
  auto vec_simple_end = CURR_TIME;

  BENCH_PRINT(
      "Database read time: "
      << TIME_DIFF(db_start_time, db_end_time) << " ms. Throughput: "
      << (2.0 * storage_MB /
          (static_cast<double>(TIME_DIFF(db_start_time, db_end_time)) / 1000))
      << " MB/s"

  );

  BENCH_PRINT(
      "Vector operation time: "
      << TIME_DIFF(vec_start_time, vec_end_time) << " ms. Throughput: "
      << (2.0 * storage_MB /
          (static_cast<double>(TIME_DIFF(vec_start_time, vec_end_time)) / 1000))
      << " MB/s");


  BENCH_PRINT(
      "Vector simple time: "
      << TIME_DIFF(vec_simple_start, vec_simple_end) << " ms. Throughput: "
      << (2.0 * storage_MB /
          (static_cast<double>(TIME_DIFF(vec_simple_start, vec_simple_end)) / 1000))
      << " MB/s");
}



void test_reading_ct() {
  // In this test case, we use the parameters from the PIR scheme to create a
  // BFV ciphertext query vector. Basically, we need fst_dim_sz many BFV
  // ciphertexts. In total, we need to read the same vector other_dim_sz many
  // times. And since there are two polynomials for each BFV ciphertext, we
  // basically need to read 2 * num_pt many polynomials. For now, lets ignore
  // the fact that we might also have RNS. We simply use ct that has 60 bits
  // coeffs.

  print_func_name(__FUNCTION__);
  PirParams pir_params;
  auto parms = pir_params.get_seal_params();    // This parameter is set to be: seal::scheme_type::bfv
  auto context_ = seal::SEALContext(parms);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);
  const size_t fst_dim_sz = pir_params.get_fst_dim_sz();
  const size_t other_dim_sz = pir_params.get_other_dim_sz();
  const auto coeff_modulus = context_.first_context_data()->parms().coeff_modulus();
  const size_t coeff_mod_count = coeff_modulus.size();
  const size_t coeff_val_cnt = DatabaseConstants::PolyDegree * coeff_mod_count; // polydegree * RNS moduli count
  const size_t one_ct_size = 2 * coeff_val_cnt; // 2 polynomials

  BENCH_PRINT("fst_dim_sz: " << fst_dim_sz << ", other_dim_sz: " << other_dim_sz << ", coeff_val_cnt: " << coeff_val_cnt);
  BENCH_PRINT("uint64_t accessed: " << 2 * fst_dim_sz * other_dim_sz * coeff_val_cnt);


  // create a vector of BFV of zeros. vector size = fst_dim_sz. Then transform to NTT form.
  std::vector<seal::Ciphertext> ct_vec(fst_dim_sz);
  for (size_t i = 0; i < fst_dim_sz; i++) {
    encryptor_->encrypt_zero_symmetric(ct_vec[i]);
    evaluator_.transform_to_ntt_inplace(ct_vec[i]);
  }



  // ================== Normal Matrix vector multiplication order ==================
  auto normal_start = CURR_TIME;
  for (size_t i = 0; i < other_dim_sz; i++) {
    for (size_t j = 0; j < fst_dim_sz; j++) {
      for (size_t k = 0; k < 2; k++) {
        #pragma GCC unroll 32
        for (size_t coeff_id = 0; coeff_id < coeff_val_cnt; coeff_id++) {
          asm volatile("" : : "r"(ct_vec[j].data(k)[coeff_id]) : "memory");
        }
      }
    }
  }
  auto normal_end = CURR_TIME;
  BENCH_PRINT("Normal time: " << TIME_DIFF(normal_start, normal_end) << " ms");


  /// ================== With Tiling ==================
  auto tiling_start = CURR_TIME;
  // const size_t tile_size = DatabaseConstants::TileSz;
  const size_t tile_size = 32;
  for (size_t k_base = 0; k_base < fst_dim_sz; k_base += tile_size) {
    for (size_t j = 0; j < other_dim_sz; ++j) {
      for (size_t k = k_base; k < std::min(k_base + tile_size,fst_dim_sz); ++k) {
        for (size_t poly_id = 0; poly_id < 2; ++poly_id) {
          #pragma GCC unroll 32
          for (size_t coeff_id = 0; coeff_id < coeff_val_cnt; ++coeff_id) {
            asm volatile("" : : "r"(ct_vec[k].data(poly_id)[coeff_id]) : "memory");
          }
        }
      }
    }
  }
  auto tiling_end = CURR_TIME;
  BENCH_PRINT("Tiling time: " << TIME_DIFF(tiling_start, tiling_end) << " ms");


  // ================== with fst_dim_data rearranged and tiling ==================
  std::vector<uint64_t> fst_dim_data(fst_dim_sz * one_ct_size);

  for (size_t k = 0; k < fst_dim_sz; k++) {
    for (size_t poly_id = 0; poly_id < 2; poly_id++) {
      auto shift = k * one_ct_size + poly_id * coeff_val_cnt;
      std::copy(ct_vec[k].data(poly_id),
                ct_vec[k].data(poly_id) + coeff_val_cnt,
                fst_dim_data.begin() + shift);
    }
  }


  auto tiling_rearrange_start = CURR_TIME;
  for (size_t k_base = 0; k_base < fst_dim_sz; k_base += tile_size) {
    for (size_t j = 0; j < other_dim_sz; ++j) {
      for (size_t k = k_base; k < std::min(k_base + tile_size,fst_dim_sz); ++k) {
        for (size_t poly_id = 0; poly_id < 2; poly_id++) {
          #pragma GCC unroll 32
          for (size_t coeff_id = 0; coeff_id < coeff_val_cnt; coeff_id++) {
            asm volatile("" : : "r"(fst_dim_data[k * one_ct_size + poly_id * coeff_val_cnt + coeff_id]) : "memory");
          }
        }
      }
    }
  }
  auto tiling_rearrange_end = CURR_TIME;
  BENCH_PRINT("Tiling and Rearrange time: " << TIME_DIFF(tiling_rearrange_start, tiling_rearrange_end) << " ms");

}


// // Function to create a large dummy vector for padding
// std::vector<uint8_t> create_padding(size_t size_in_bytes) {
//     return std::vector<uint8_t>(size_in_bytes, 0);
// }


// void test_trivial_loop() {
//     print_func_name(__FUNCTION__);
    
//     PirParams pir_params;
//     const size_t num_pt = pir_params.get_num_pt();
//     const auto seal_params = pir_params.get_seal_params();
//     const size_t coeff_count = seal_params.poly_modulus_degree();

//     // Define padding size (e.g., 1 MB)
//     constexpr size_t padding_size = 1 << 20; // 1 MB

//     auto padding1 = create_padding(padding_size);
//     std::vector<uint64_t> vec1(num_pt * coeff_count);

//     auto padding2 = create_padding(padding_size);
//     std::vector<uint64_t> vec2(num_pt * coeff_count);

//     auto padding3 = create_padding(padding_size);
//     std::vector<uint64_t> vec3(num_pt * coeff_count);

//     // Fill vec1 and vec2 with random values
//     std::mt19937_64 rng(12345ULL);
//     std::uniform_int_distribution<uint64_t> dist_matrix(0ULL, std::numeric_limits<uint64_t>::max());
//     for (size_t i = 0; i < num_pt * coeff_count; i++) {
//         vec1[i] = dist_matrix(rng);
//         vec2[i] = dist_matrix(rng);
//         vec3[i] = dist_matrix(rng);
//     }

//     uint64_t sum = 0;
//     // Perform the linear scan
//     auto start = CURR_TIME;
//     for (size_t i = 0; i < num_pt * coeff_count; i++) {
//         // sum ^= vec1[i] ^ vec2[i] ^ vec3[i];
//         asm volatile("" : : "r"(vec1[i]) : "memory");
//         asm volatile("" : : "r"(vec2[i]) : "memory");
//         asm volatile("" : : "r"(vec3[i]) : "memory");
//     }
//     auto end = CURR_TIME;

//     BENCH_PRINT("Time: " << TIME_DIFF(start, end) << " ms");
//     BENCH_PRINT("Sum: " << sum);
// }


void test_matrix_mult () {
  print_func_name(__FUNCTION__);

  // PirParams pir_params;
  // const size_t num_pt = pir_params.get_num_pt();
  // const auto seal_params = pir_params.get_seal_params();
  // const size_t coeff_count = seal_params.poly_modulus_degree();

  const size_t num_pt = 524288;
  const size_t coeff_count = 1024;
  const size_t db_size = num_pt * coeff_count * sizeof(uint64_t) / 1024 / 1024; // in MB
  BENCH_PRINT("db_size: " << db_size << " MB");


  // Create a random generator
  // (Seed with a fixed value or std::random_device{}() for more variability)
  std::mt19937_64 rng(12345ULL);

  // Distribution for the matrix entries: all possible uint64_t
  std::uniform_int_distribution<uint64_t> dist_matrix(
      0ULL, std::numeric_limits<uint64_t>::max()
  );

  // Distribution for the vector entries: only 0 or 1
  std::uniform_int_distribution<uint64_t> dist_vector(0ULL, 1ULL);

  // Allocate and fill the matrix (coeff_count rows, each row has num_pt columns)
  std::vector<std::vector<uint64_t>> matrix(coeff_count, std::vector<uint64_t>(num_pt));
  for (size_t i = 0; i < coeff_count; i++)
  {
      for (size_t j = 0; j < num_pt; j++)
      {
          matrix[i][j] = dist_matrix(rng);
      }
  }

  // Allocate and fill the vector (length = num_pt) with 0s or 1s
  std::vector<uint64_t> vec(num_pt);
  for (size_t j = 0; j < num_pt; j++)
  {
      vec[j] = dist_vector(rng);
  }

  // Prepare a result vector (length = coeff_count)
  std::vector<uint64_t> result(coeff_count, 0ULL);


  auto matrix_mult_start = CURR_TIME;
  // --- Matrix-vector multiplication ---
  // result[i] = sum_{j=0}^{num_pt-1} matrix[i][j] * vec[j]
  for (size_t i = 0; i < coeff_count; i++)
  {
      uint64_t sum = 0ULL;
      for (size_t j = 0; j < num_pt; j++)
      {
          // In 64-bit, watch out for potential overflow if matrix entries are very large,
          // but standard unsigned 64-bit wrap-around rules will apply
          sum += matrix[i][j] * vec[j];
      }
      result[i] = sum;
  }
  auto matrix_mult_end = CURR_TIME; 

  auto elapsed_time_us = std::chrono::duration_cast<std::chrono::microseconds>(matrix_mult_end - matrix_mult_start).count();

  BENCH_PRINT("Matrix-vector multiplication time: " << elapsed_time_us << " us");

  // throughput in MB/s, intermediate result is in us
  double throughput = db_size / (elapsed_time_us / 1000000.0);
  BENCH_PRINT("Throughput: " << throughput << " MB/s");
  

  // sum the result to force reading the result
  uint64_t tot = 0;
  for (size_t i = 0; i < coeff_count; i++) {
    tot += result[i];
  }

  // print the total
  BENCH_PRINT("Total: " << tot);


}
