#include "tests.h"
#include "gsw_eval.h"
#include "pir.h"
#include "server.h"
#include "client.h"
#include "utils.h"
#include <cassert>
#include <iostream>
#include <bitset>

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
  // test_pir();
  bfv_example();
  serialization_example();
  test_external_product();
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
  PRINT_BAR;

  // ============= Now let's try BFV multiplied by a constant in ntt form ==============
  seal::Plaintext scalar(coeff_count);
  scalar[0] = 2;
  scalar[1] = 3;
  BENCH_PRINT("Vector a: " << a.to_string());
  BENCH_PRINT("Scalar: " << scalar.to_string());
  evaluator_.transform_to_ntt_inplace(scalar, context_.first_parms_id()); // This happens in preprocess_ntt
  // Now instead of using multiply_plain, I want to demonstrate what happens in the first dimension evaluation. 
  // This is demonstrating what you can do in ntt form, but the actual order of computation in OnionPIRv2 fst dim can be different.
  size_t rns_mod_cnt = coeff_modulus.size() - 1;
  std::vector<uint128_t> res(coeff_count * rns_mod_cnt);
  std::fill(res.begin(), res.end(), 0);
  uint64_t *ct0_ptr = a_encrypted.data(0);
  uint64_t *ct1_ptr = a_encrypted.data(1);
  uint128_t *res0_ptr = res.data();
  uint128_t *res1_ptr = res.data() +  coeff_count * rns_mod_cnt * 2;
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
      scal_mul_ct0_ptr[i + j * coeff_count] = res0_ptr[i + j * coeff_count] % curr_mod;
      scal_mul_ct1_ptr[i + j * coeff_count] = res1_ptr[i + j * coeff_count] % curr_mod;
    }
  }
  evaluator_.transform_from_ntt_inplace(scalar_mul_result);
  decryptor_->decrypt(scalar_mul_result, result);
  BENCH_PRINT("NTT x scalar result: " << result.to_string());  // and the result is correct! NTT form polynomial is multiplicative
  /*
  Now, in the old OnionPIR, this kind of elementwise multiplication is computed for num_poly many times. That is, the smallest operation
  is this vector-vector elementwise multiplication. However, this is actually bad for the cache. TODO: I need to find a way to optimize this.
  */
  PRINT_BAR;

  // ============= Now let's try BFV multiplied by two identical constants then subtract ==============
  // Actually, this creates something called transparant ciphertext, which is warned in the SEAL documentation.
  seal::Plaintext constant(coeff_count);
  constant[0] = 2;
  seal::Ciphertext fst_mult_result, snd_mult_result;
  evaluator_.multiply_plain(a_encrypted, scalar, fst_mult_result);
  evaluator_.multiply_plain(a_encrypted, scalar, snd_mult_result);
  BENCH_PRINT("If you see an error about 'transparent ciphertext' below, "
              "please make sure you are using "
              "-DSEAL_THROW_ON_TRANSPARENT_CIPHERTEXT=OFF when building SEAL");
  evaluator_.sub_inplace(fst_mult_result, snd_mult_result);
  evaluator_.transform_from_ntt_inplace(fst_mult_result);
  decryptor_->decrypt(fst_mult_result, result);
  BENCH_PRINT("You should see a zero ¬_¬: " << result.to_string()); 
}


void serialization_example() {
  print_func_name(__FUNCTION__);
  PirParams pir_params;
  auto params = pir_params.get_seal_params();
  auto context_ = seal::SEALContext(params);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);

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
  // bytes to store the prng seed. Notice that they have common prefix.
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[1]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[2]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[3]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[4]));
  BENCH_PRINT("Seed: \t\t" << std::bitset<64>(ptr_1[5]));
  
  auto mods = pir_params.get_coeff_modulus();
  auto plain_modulus = pir_params.get_plain_mod();
  uint128_t ct_mod = 1; 
  for (size_t mod_id = 0; mod_id < mods.size(); mod_id++) {
    ct_mod *= mods[mod_id].value();
  }
  uint128_t delta = ct_mod / plain_modulus;  // delta = floor (ciphertext modulus / plaintext modulus)
  uint128_t message = 15;
  uint128_t to_add = delta * message;
  auto padding = params.poly_modulus_degree();
  for (size_t mod_id = 0; mod_id < mods.size(); mod_id++) {
    ptr_0[mod_id * padding] = (ptr_0[mod_id * padding] + (to_add % mods[mod_id].value())) % mods[mod_id].value();
  }

  // write the serializable object to the stream
  auto s2_size = new_seeded_zero.save(data_stream); // Storing new ciphertext with a seed

  BENCH_PRINT("Size of the ciphertexts: " << new_seeded_zero.size());

  // ================== Deserialize and decrypt the ciphertexts ==================
  seal::Ciphertext raw_ct, orig_ct, new_ct;
  raw_ct.load(context_, data_stream);  // loading the raw zero
  orig_ct.load(context_, data_stream);  // loading the original zero
  new_ct.load(context_, data_stream); // loading the new ciphertext with a seed 

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

// This is a BFV x GSW example
void test_external_product() {
  print_func_name(__FUNCTION__);
  PirParams pir_params;
  const auto params = pir_params.get_seal_params();
  auto context_ = seal::SEALContext(params);
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);
  const size_t coeff_count = DatabaseConstants::PolyDegree;

  // the test data vector a and results are both in BFV scheme.
  seal::Plaintext a(coeff_count), result;
  std::vector<uint64_t> b(coeff_count); // vector b is in the context of GSW scheme.
  a[0] = 1; a[1] = 2; a[2] = 3;
  b[0] = 2;
  BENCH_PRINT("Vector a: " << a.to_string());
  std::string b_str = "Vector b: ";
  for (int i = 0; i < 5; i++) {
    b_str += std::to_string(b[i]) + " ";
  }
  BENCH_PRINT(b_str);  
  
  seal::Ciphertext a_encrypted;    // encrypted "a" will be stored here. 
  encryptor_->encrypt_symmetric(a, a_encrypted);

  // encrypt the plaintext b to GSW ciphertext
  std::vector<seal::Ciphertext> temp_gsw;
  GSWEval data_gsw(pir_params, pir_params.get_l(), pir_params.get_base_log2());
  data_gsw.plain_to_gsw(b, *encryptor_, secret_key_, temp_gsw); 
  GSWCiphertext b_gsw;
  data_gsw.sealGSWVecToGSW(b_gsw, temp_gsw);
  data_gsw.gsw_ntt_negacyclic_harvey(b_gsw);

  // actual external product
  std::cout << "Noise budget before: " << decryptor_->invariant_noise_budget(a_encrypted) << std::endl;
  data_gsw.external_product(b_gsw, a_encrypted, a_encrypted);
  evaluator_.transform_from_ntt_inplace(a_encrypted);
  decryptor_->decrypt(a_encrypted, result);
  // output decrypted result
  std::cout << "External product result: " << result.to_string() << std::endl;

  std::cout << "Noise budget after: " << decryptor_->invariant_noise_budget(a_encrypted)
            << std::endl;
}