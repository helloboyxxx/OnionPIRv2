#include "tests.h"
#include "gsw_eval.h"
#include "pir.h"
#include "server.h"
#include "client.h"
#include "utils.h"
#include "logging.h"
#include "matrix.h"
#include <cassert>
#include <iostream>
#include <bitset>

#ifdef HAVE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Core>
#endif

#define EXPERIMENT_ITERATIONS 5

void PirTest::run_tests() {
  test_pir();
  // bfv_example();
  // serialization_example();
  // test_external_product();
  // test_single_mat_mult();
  // test_fst_dim_mult();
  // test_batch_decomp();
  // test_fast_expand_query();
}


void PirTest::test_pir() {
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
    // size_t query_index = rand() % pir_params.get_num_entries();
    size_t query_index = 1; 

    // ============= CLIENT ===============
    TIME_START(CLIENT_TOT_TIME);
    // seal::Ciphertext query = client.generate_query(query_index);
    seal::Ciphertext query = client.fast_generate_query(query_index);
    query_size = client.write_query_to_stream(query, data_stream);
    TIME_END(CLIENT_TOT_TIME);
    
    // ============= SERVER ===============
    TIME_START(SERVER_TOT_TIME);
    std::vector<seal::Ciphertext> result = server.make_query(client_id, data_stream);
    TIME_END(SERVER_TOT_TIME);

    // ============= CLIENT ===============
    TIME_START(CLIENT_TOT_TIME);
    // client gets result from the server and decrypts it
    std::vector<seal::Plaintext> decrypted_result = client.decrypt_result(result);
    Entry response_entry = client.get_entry_from_plaintext(query_index, decrypted_result[0]);
    TIME_END(CLIENT_TOT_TIME);

    // write the result to the stream to test the size
    std::stringstream result_stream;
    response_size = result[0].save(result_stream);
    result_stream.str(std::string()); // clear the stream

    // Directly get the plaintext from server. Not part of PIR.
    Entry actual_entry = server.direct_get_entry(query_index);
    // extract and print the actual entry index
    uint64_t actual_entry_idx = utils::get_entry_idx(actual_entry);
    uint64_t resp_entry_idx = utils::get_entry_idx(response_entry);
    
    END_EXPERIMENT();
    // ============= PRINTING RESULTS ===============    
    DEBUG_PRINT("\t\tWanted/resp/actual idx:\t" << query_index << " / " << resp_entry_idx << " / " << actual_entry_idx);
    #ifdef _DEBUG
    PRINT_RESULTS(i+1);
    #endif

    if (utils::entry_is_equal(response_entry, actual_entry)) {
      // print a green success message
      std::cout << "\033[1;32mSuccess!\033[0m" << std::endl;
      success_count++;
    } else {
      // print a red failure message
      std::cout << "\033[1;31mFailure!\033[0m" << std::endl;
      std::cout << "PIR Result:\t";
      utils::print_entry(response_entry, 20);
      std::cout << "Actual Entry:\t";
      utils::print_entry(actual_entry, 20);
    }
    PRINT_BAR;
  }

  double avg_server_time = GET_AVG_TIME(SERVER_TOT_TIME);
  double throughput = pir_params.get_DBSize_MB() / (avg_server_time / 1000);

  // ============= PRINTING FINAL RESULTS ===============]
  BENCH_PRINT("                                FINAL RESULTS")
  PRINT_BAR;
  BENCH_PRINT("Success rate: " << success_count << "/" << EXPERIMENT_ITERATIONS);
  BENCH_PRINT("galois key size: " << galois_key_size << " bytes");
  BENCH_PRINT("gsw key size: " << gsw_key_size << " bytes");
  BENCH_PRINT("total key size: " << static_cast<double>(galois_key_size + gsw_key_size) / 1024 / 1024 << "MB");
  BENCH_PRINT("query size: " << query_size << " bytes");
  BENCH_PRINT("response size: " << response_size << " bytes");
  
  PRETTY_PRINT();
  BENCH_PRINT("Server throughput: " << throughput << " MB/s");
}

  // This is an example of how to use the BFV scheme in SEAL and in our PIR scheme.
void PirTest::bfv_example() {
  print_func_name(__FUNCTION__);
  // You need a a chunk of code to init the seal parameters. Here is the minimum you need:
  seal::EncryptionParameters params(seal::scheme_type::bfv);
  const size_t coeff_count = 4096;  // you can try other powers of two.
  params.set_poly_modulus_degree(coeff_count); // example: a_1 x^4095 + a_2 x^4094 + ...
  const uint64_t pt_mod = utils::generate_prime(49); // 49 bits for the plain modulus, then you can use 48 bits for storing data.
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


void PirTest::serialization_example() {
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
void PirTest::test_external_product() {
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
  a[0] = 1; a[1] = 2; a[2] = 4;
  b[0] = 2; // You can also try 1, then you can do external product hundreds of times.
  BENCH_PRINT("Vector a: " << a.to_string());
  std::string b_str = "Vector b: ";
  for (int i = 0; i < 5; i++) {
    b_str += std::to_string(b[i]) + " ";
  }
  BENCH_PRINT(b_str);  
  
  seal::Ciphertext a_encrypted;    // encrypted "a" will be stored here. 
  encryptor_->encrypt_symmetric(a, a_encrypted);

  // encrypt the plaintext b to GSW ciphertext
  // You can also try different gsw_l and base_log2. But you need to follow the equation:
  // base_log2 = (bits + l - 1) / l; where bits is the bit width of the ciphertext modulus. 
  const size_t gsw_l = pir_params.get_l();
  const size_t base_log2 = pir_params.get_base_log2();
  GSWEval data_gsw(pir_params, gsw_l, base_log2);
  std::vector<seal::Ciphertext> temp_gsw;
  data_gsw.plain_to_gsw(b, *encryptor_, secret_key_, temp_gsw); // In OnionPIR, client use a similar function to encrypt the secret key. 
  GSWCiphertext b_gsw;
  data_gsw.seal_GSW_vec_to_GSW(b_gsw, temp_gsw);
  data_gsw.gsw_ntt_negacyclic_harvey(b_gsw);  // We need NTT form RGSW.

  // actual external product
  BENCH_PRINT("Noise budget before: " << decryptor_->invariant_noise_budget(a_encrypted));
  const size_t num_iter = 10; // And you can do this external product many times when the data in GSW is small. 
  for (size_t i = 0; i < num_iter; ++i) {
    data_gsw.external_product(b_gsw, a_encrypted, a_encrypted); // The decomposition requires coefficient form BFV
    evaluator_.transform_from_ntt_inplace(a_encrypted);
    decryptor_->decrypt(a_encrypted, result);
    // output decrypted result
    BENCH_PRINT("External product result: " << result.to_string());
  }
  BENCH_PRINT("Noise budget after: " << decryptor_->invariant_noise_budget(a_encrypted));
  PRINT_BAR;
  // ============= Now, let's try profiling the external product ==============
  // I prefer to use samply. It works for both mac and linux. 
  // I will also log the time elapsed in the external product function.
  
  // when poly_degree = 2048, a single BFV is 32KB.
  const size_t num_samples = 10000;
  std::vector<seal::Ciphertext> a_encrypted_vec(num_samples);
  for (size_t i = 0; i < num_samples; i++) {
    encryptor_->encrypt_symmetric(a, a_encrypted_vec[i]);
  }
  CLEAN_TIMER();
  TIME_START(OTHER_DIM_MUX_EXTERN);
  for (size_t i = 0; i < num_samples; i++) {
    data_gsw.external_product(b_gsw, a_encrypted_vec[i], a_encrypted_vec[i]);

    TIME_START("inverse ntt");
    evaluator_.transform_from_ntt_inplace(a_encrypted_vec[i]); // Try uncommenting this line and see the difference.
    TIME_END("inverse ntt");
  }
  TIME_END(OTHER_DIM_MUX_EXTERN);

  // print the timing result
  // roughly the result should be in the structure: 
  /*
    External product
      - Decomposition
        - memcpy
        - compose
        - right shift
        - decompose
        - ntt
      - mat mat mult
      - delayed mod
    Inverse NTT
  */
  END_EXPERIMENT();
  PRINT_RESULTS(); 
}


void PirTest::test_single_mat_mult() {
  print_func_name(__FUNCTION__);
  CLEAN_TIMER();
  // This is testing mat mat multiplication: A x B = C 
  // with a special condition that the width of B is 2 and width of A is DatabaseConstants::MaxFstDimSz.
  // Ideally, this tells the limit of the first dimension throughput.
  constexpr size_t rows = 1 << 20; 
  constexpr size_t cols = DatabaseConstants::MaxFstDimSz; 
  constexpr size_t b_cols = 2; // two polynomials 
  constexpr size_t db_size = rows * cols * sizeof(uint64_t);  // we only care the big matrix
  BENCH_PRINT("Matrix size: " << db_size / 1024 / 1024 << " MB");

  //Allocate memory for A, B, out.
  std::vector<uint64_t> A_data(rows * cols);
  std::vector<uint64_t> B_data(cols * b_cols);
  std::vector<uint64_t> C_data(rows * b_cols);
  std::vector<uint128_t> C_data128(rows * b_cols);

  // Fill A and B with random data
  utils::fill_rand_arr(A_data.data(), rows * cols);
  utils::fill_rand_arr(B_data.data(), cols * b_cols);
  // Wrap them in our matrix_t structures
  matrix_t A_mat { A_data.data(), rows, cols, 1 };
  matrix_t B_mat { B_data.data(), cols, b_cols, 1 };
  matrix_t C_mat { C_data.data(), rows, b_cols, 1 };
  matrix128_t C_mat128 { C_data128.data(), rows, b_cols, 1 };
  uint128_t sum128 = 0;
  size_t sum = 0;


  // ============= naive mat-vec mult ==============
  const std::string NAIVE_MAT_MULT = "Naive matrix multiplication";
  TIME_START(NAIVE_MAT_MULT);
  naive_mat_vec(&A_mat, &B_mat, &C_mat);
  TIME_END(NAIVE_MAT_MULT);

  // ============= naive mat-vec mult 128 bits ==============
  const std::string NAIVE_MAT_MULT_128 = "Naive matrix multiplication 128 bits";
  TIME_START(NAIVE_MAT_MULT_128);
  naive_mat_vec_128(&A_mat, &B_mat, &C_mat128);
  TIME_END(NAIVE_MAT_MULT_128);

  // ============= naive level mat mult ==============
  const std::string NAIVE_LEVEL_MAT_MULT = "Naive level matrix multiplication";
  TIME_START(NAIVE_LEVEL_MAT_MULT);
  naive_level_mat_mat(&A_mat, &B_mat, &C_mat);
  TIME_END(NAIVE_LEVEL_MAT_MULT);

  // ============= naive level mat mult 128 bits ==============
  const std::string NAIVE_LEVEL_MAT_MULT_128 = "Naive level matrix multiplication 128 bits";
  TIME_START(NAIVE_LEVEL_MAT_MULT_128);
  naive_level_mat_mat_128(&A_mat, &B_mat, &C_mat128);
  TIME_END(NAIVE_LEVEL_MAT_MULT_128);

  // ============= level mat mult ==============
  const std::string LV_MAT_MULT = "Matrix multiplication";
  TIME_START(LV_MAT_MULT);
  level_mat_mat(&A_mat, &B_mat, &C_mat);
  TIME_END(LV_MAT_MULT);

  // ============= level mat mult 128 bits ==============
  const std::string LV_MAT_MULT_128 = "Matrix multiplication 128 bits";
  TIME_START(LV_MAT_MULT_128);
  level_mat_mat_128(&A_mat, &B_mat, &C_mat128);
  TIME_END(LV_MAT_MULT_128);

  // ============= Eigen mat mult ==============
  const std::string EIGEN_MULT = "Matrix multiplication Eigen";
  Eigen::setNbThreads(1);  // Force Eigen to use only 1 thread
  TIME_START(EIGEN_MULT);
  level_mat_mult_eigen(&A_mat, &B_mat, &C_mat);
  TIME_END(EIGEN_MULT);

  // ============= avx mat mat mult 128 bits ==============
  const std::string AVX_MAT_MULT_128 = "AVX matrix multiplication 128 bits";
  TIME_START(AVX_MAT_MULT_128);
  avx_mat_mat_mult_128(A_data.data(), B_data.data(), C_data128.data(), rows, cols);
  TIME_END(AVX_MAT_MULT_128);


  // ============= Print the results to avoid over optimization ==============
  for (size_t i = 0; i < rows * b_cols; i++) { sum += C_data[i]; }
  BENCH_PRINT("Sum: " << sum);
  for (size_t i = 0; i < rows * b_cols; i++) { sum128 += C_data128[i]; }
  BENCH_PRINT("Sum: " << utils::uint128_to_string(sum128));

  // ============= Profiling the matrix multiplication ==============
  END_EXPERIMENT();
  // PRINT_RESULTS(); // uncomment this line to see the actual time elapsed in each function.
  PRINT_BAR;
  double naive_time = GET_AVG_TIME(NAIVE_MAT_MULT);
  double naive_time_128 = GET_AVG_TIME(NAIVE_MAT_MULT_128);
  double naive_lvl_time = GET_AVG_TIME(NAIVE_LEVEL_MAT_MULT);
  double naive_lvl_time_128 = GET_AVG_TIME(NAIVE_LEVEL_MAT_MULT_128);
  double lv_time = GET_AVG_TIME(LV_MAT_MULT);
  double lv_time_128 = GET_AVG_TIME(LV_MAT_MULT_128);
  double eigen_time = GET_AVG_TIME(EIGEN_MULT);
  double avx_time_128 = GET_AVG_TIME(AVX_MAT_MULT_128);

  double naive_throughput = db_size / (naive_time * 1000);
  double naive_throughput_128 = db_size / (naive_time_128 * 1000);
  double naive_lvl_throughput = db_size / (naive_lvl_time * 1000);
  double naive_lvl_throughput_128 = db_size / (naive_lvl_time_128 * 1000);
  double lv_throughput = db_size / (lv_time * 1000);
  double lv_throughput_128 = db_size / (lv_time_128 * 1000);
  double eigen_throughput = db_size / (eigen_time * 1000);
  double avx_throughput_128 = db_size / (avx_time_128 * 1000);

  BENCH_PRINT("Naive mat vec throughput: \t\t" << naive_throughput << " MB/s");
  BENCH_PRINT("Naive mat vec 128 throughput: \t\t" << naive_throughput_128 << " MB/s");
  BENCH_PRINT("Naive level mat mat throughput: \t" << naive_lvl_throughput << " MB/s");
  BENCH_PRINT("Naive level mat mat 128 throughput: \t" << naive_lvl_throughput_128 << " MB/s");
  BENCH_PRINT("Level mat mat throughput: \t\t" << lv_throughput << " MB/s");
  BENCH_PRINT("Level mat mat 128 throughput: \t\t" << lv_throughput_128 << " MB/s");
  BENCH_PRINT("Eigen mat mat throughput: \t\t" << eigen_throughput << " MB/s");
  BENCH_PRINT("AVX mat mat 128 throughput: \t\t" << avx_throughput_128 << " MB/s");
}


void PirTest::test_fst_dim_mult() {
  print_func_name(__FUNCTION__);
  CLEAN_TIMER();
  // for this test, I want to know if the matrix multiplication is memory bound
  // or compute bound. If possible, please re-write this test case for GPU as
  // well as it indicates the limit of the first dimension.

  // Let's write the best code we can to compute (m x n) x (n x p) matrix
  // multiplication for k times.
  constexpr size_t m = 1 << 9; // the other_dim_sz
  constexpr size_t n = DatabaseConstants::MaxFstDimSz;
  constexpr size_t p = 2; // coz we have only 2 polynomials in the ciphertext.
  constexpr size_t k = DatabaseConstants::PolyDegree;
  constexpr size_t db_size = m * n * k * sizeof(uint64_t);  // we only care the big matrix
  PirParams pir_params;
  BENCH_PRINT("Matrix size: " << db_size / 1024 / 1024 << " MB");

  // Allocate memory for A, B, out. 
  // We interpret these as stacked (k) matrices.
  std::vector<uint64_t> A_data(m * n * k);
  std::vector<uint64_t> B_data(n * p * k);
  std::vector<uint64_t> C_data(m * p * k);
  std::vector<uint128_t> C_data_128(m * p * k);
  // Fill A and B with random data
  utils::fill_rand_arr(A_data.data(), m * n * k); 
  utils::fill_rand_arr(B_data.data(), n * p * k);
  // Wrap them in our matrix_t structures
  matrix_t A_mat { A_data.data(), m, n, k };
  matrix_t B_mat { B_data.data(), n, p, k };
  matrix_t C_mat { C_data.data(), m, p, k };
  matrix128_t C_mat_128 { C_data_128.data(), m, p, k };
  size_t sum = 0;
  uint128_t sum128 = 0;


  // ============= Naive level mat mat ==============
  const std::string NAIVE_MAT_MULT = "Naive level mat mat";
  TIME_START(NAIVE_MAT_MULT);
  naive_level_mat_mat(&A_mat, &B_mat, &C_mat);
  TIME_END(NAIVE_MAT_MULT);

  // ============= Naive level mat mat 128bits ==============
  const std::string NAIVE_MAT_MULT_128 = "Naive level mat mat 128 bits";
  TIME_START(NAIVE_MAT_MULT_128);
  naive_level_mat_mat_128(&A_mat, &B_mat, &C_mat_128);
  TIME_END(NAIVE_MAT_MULT_128);


  // ===================== Performing matrix multiplication by levels ===================== 
  // So, the idea is that we can do k many matrix matrix
  // multiplications. Instead of doing the component wise multiplication, which
  // I think is not cache friendly, matrix multiplication can benefit from local caching. 
  // Note that these two functions are processing the data in a very different order. 
  const std::string LV_MAT_MULT = "Matrix multiplication";
  TIME_START(LV_MAT_MULT);
  level_mat_mat(&A_mat, &B_mat, &C_mat);
  TIME_END(LV_MAT_MULT);

  // ============= level mat mult 128 bits ==============
  const std::string LV_MAT_MULT_128 = "Matrix multiplication 128 bits";
  TIME_START(LV_MAT_MULT_128);
  level_mat_mat_128(&A_mat, &B_mat, &C_mat_128);
  TIME_END(LV_MAT_MULT_128);

  // ============= Level mat mult direct mod ==============
  const std::string LV_MAT_MULT_DIRECT_MOD = "Matrix multiplication direct mod";
  seal::Modulus mod = pir_params.get_coeff_modulus()[0];
  TIME_START(LV_MAT_MULT_DIRECT_MOD);
  level_mat_mat_direct_mod(&A_mat, &B_mat, &C_mat, mod);
  TIME_END(LV_MAT_MULT_DIRECT_MOD);

  // ============= OnionPIR v1 elementwise multiplication ==============
  const std::string ELEM_MULT = "elementwise multiplication";
  TIME_START(ELEM_MULT);
  component_wise_mult(&A_mat, &B_mat, &C_mat); 
  TIME_END(ELEM_MULT);

  // ============= component wise mult 128 bits ==============
  const std::string ELEM_MULT_128 = "Old elementwise multiplication 128 bits";
  TIME_START(ELEM_MULT_128);
  component_wise_mult_128(&A_mat, &B_mat, &C_mat_128);
  TIME_END(ELEM_MULT_128);

  // ============= component wise mult direct mod using hexl ==============
  const std::string ELEM_MULT_DIRECT_MOD = "elementwise multiplication direct mod";
  uint64_t mod_val = pir_params.get_coeff_modulus()[0].value();
  TIME_START(ELEM_MULT_DIRECT_MOD);
  component_wise_mult_direct_mod(&A_mat, &B_mat, C_data.data(), mod_val);
  TIME_END(ELEM_MULT_DIRECT_MOD);


  // ============= level mat mult using Eigen ==============
  #ifdef HAVE_EIGEN
  const std::string EIGEN_MULT = "Matrix multiplication Eigen";
  Eigen::setNbThreads(1);  // Force Eigen to use only 1 thread
  TIME_START(EIGEN_MULT);
  level_mat_mult_eigen(&A_mat, &B_mat, &C_mat);
  TIME_END(EIGEN_MULT);
  #endif

  // some simple code to make sure it is not optimized out
  sum = 0; 
  for (size_t i = 0; i < m * p * k; i++) { sum += C_data[i]; }
  BENCH_PRINT("Sum: " << sum);
  sum128 = 0;
  for (size_t i = 0; i < m * p * k; i++) { sum128 += C_data_128[i]; }
  BENCH_PRINT("Sum: " << utils::uint128_to_string(sum128));


  // ============= Profiling the matrix multiplication ==============
  END_EXPERIMENT();
  // PRINT_RESULTS(); // uncomment this line to see the actual time elapsed in each function.
  PRINT_BAR;

  // Let's calculate the throughput of the matrix multiplication, express in MB/s
  double naive_mat_mult_time = GET_AVG_TIME(NAIVE_MAT_MULT);
  double naive_mat_mult_128_time = GET_AVG_TIME(NAIVE_MAT_MULT_128);
  double level_mat_mult_time = GET_AVG_TIME(LV_MAT_MULT);
  double level_mat_mult_128_time = GET_AVG_TIME(LV_MAT_MULT_128);
  double level_mat_mult_direct_mod_time = GET_AVG_TIME(LV_MAT_MULT_DIRECT_MOD);
  double old_elementwise_mult_time = GET_AVG_TIME(ELEM_MULT);
  double elementwise_mult_128_time = GET_AVG_TIME(ELEM_MULT_128);
  double elementwise_mult_direct_mod_time = GET_AVG_TIME(ELEM_MULT_DIRECT_MOD);

  double naive_throughput = db_size / (naive_mat_mult_time * 1000);
  double naive_throughput_128 = db_size / (naive_mat_mult_128_time * 1000);
  double level_mat_mult_throughput = db_size / (level_mat_mult_time * 1000);
  double level_mat_mult_128_throughput = db_size / (level_mat_mult_128_time * 1000);
  double level_mat_mult_direct_mod_throughput = db_size / (level_mat_mult_direct_mod_time * 1000);
  double old_elementwise_mult_throughput = db_size / (old_elementwise_mult_time * 1000); 
  double elementwise_mult_128_throughput = db_size / (elementwise_mult_128_time * 1000);
  double elementwise_mult_direct_mod_throughput = db_size / (elementwise_mult_direct_mod_time * 1000);

  BENCH_PRINT("Matrix size: " << db_size / 1024 / 1024 << " MB");
  BENCH_PRINT("Naive level mat mat throughput: \t" << (size_t)naive_throughput << " MB/s");
  BENCH_PRINT("Naive level mat mat 128 throughput: \t" << (size_t)naive_throughput_128 << " MB/s");
  BENCH_PRINT("Level mat mat throughput: \t\t" << (size_t) level_mat_mult_throughput << " MB/s");
  BENCH_PRINT("Level mat mat 128 throughput: \t\t" << (size_t)level_mat_mult_128_throughput << " MB/s");
  BENCH_PRINT("Level mat mat direct mod throughput: \t" << (size_t)level_mat_mult_direct_mod_throughput << " MB/s");
  BENCH_PRINT("Elementwise mat throughput: \t\t" << (size_t)old_elementwise_mult_throughput << " MB/s");
  BENCH_PRINT("Elementwise mat 128 throughput: \t" << (size_t)elementwise_mult_128_throughput << " MB/s");
  BENCH_PRINT("Elementwise mat direct mod throughput: \t" << (size_t)elementwise_mult_direct_mod_throughput << " MB/s");
  #ifdef HAVE_EIGEN
  double level_mat_mult_eigen_time = GET_AVG_TIME(EIGEN_MULT);
  double level_mat_mult_eigen_throughput = db_size / (level_mat_mult_eigen_time * 1000);
  BENCH_PRINT("Level mat mat Eigen throughput: \t" << (size_t)level_mat_mult_eigen_throughput << " MB/s");
  #endif
}


void PirTest::test_batch_decomp() {
  // I observed that we do external product for each polynomial in the selected
  // database after first dimension multiplication. 
  // Then, I wonder if it is possible to batch the external product. 
  // ! here, batch means vectorized operations. Essentially: 
  // We can first do the vectorized BFV homomorphic subtraction,
  // then we decompose all the ciphertext (which generates $l$ times more ciphertexts)
  // then then the external product matrix multiplication,
  // then the delayed modulus reduction.
  // then the vectorized BFV homomorphic addition.
  // ? If we batch everything, can we save some time?

  // Well, I use thi test to test the batched CRT and NTT operations. 
  // They are the most time-consuming operations in the external product.
  print_func_name(__FUNCTION__);
  CLEAN_TIMER();

  PirParams pir_params;
  const size_t other_dim_sz = pir_params.get_other_dim_sz();
  const auto params = pir_params.get_seal_params();
  auto context_ = seal::SEALContext(params);
  auto context_data = context_.first_context_data();
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);
  auto ntt_tables = context_data->small_ntt_tables();
  seal::util::RNSBase *rns_base = context_data->rns_tool()->base_q();
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  auto pool = seal::MemoryManager::GetPool();


  // create other_dim_sz many ciphertexts with zero plaintext
  std::vector<seal::Ciphertext> ct_vec(other_dim_sz);
  for (size_t i = 0; i < other_dim_sz; i++) {
    seal::Ciphertext ct;
    encryptor_->encrypt_zero_symmetric(ct);
    ct_vec[i] = ct;
  }

  TIME_START("Batch compose");
  for (size_t i = 0; i < other_dim_sz; i++) {
    for (size_t poly_id = 0; poly_id < ct_vec[i].size(); poly_id++) {
      auto ct_ptr = ct_vec[i].data(poly_id);
      rns_base->compose_array(ct_ptr, coeff_count, pool);
    }
  }
  TIME_END("Batch compose");

  // decompose
  TIME_START("Batch decompose");
  for (size_t i = 0; i < other_dim_sz; i++) {
    for (size_t poly_id = 0; poly_id < ct_vec[i].size(); poly_id++) {
      for (size_t p = 0; p < pir_params.get_l(); p++) {
        auto ct_ptr = ct_vec[i].data(poly_id);
        rns_base->decompose_array(ct_ptr, coeff_count, pool);
      }
    }
  }
  TIME_END("Batch decompose");

  TIME_START("Batch NTT");
  for (size_t i = 0; i < other_dim_sz; i++) {
    for (size_t poly_id = 0; poly_id < ct_vec[i].size(); poly_id++) {
      for (size_t p = 0; p < pir_params.get_l(); p++) {
        auto ct_ptr = ct_vec[i].data(poly_id);
        ct_vec[i].is_ntt_form() = true; 
        seal::util::ntt_negacyclic_harvey(ct_ptr, ntt_tables[0]);
        seal::util::ntt_negacyclic_harvey(ct_ptr + coeff_count, ntt_tables[1]);
      }
    }
  }
  TIME_END("Batch NTT");

  // to make sure that the optimization is not too aggressive
  size_t dummy_sum = 0;
  for (size_t i = 0; i < other_dim_sz; i++) {
    for (size_t poly_id = 0; poly_id < ct_vec[i].size(); poly_id++) {
      auto ct_ptr = ct_vec[i].data(poly_id);
      for (size_t j = 0; j < coeff_count; j++) {
        dummy_sum += ct_ptr[j];
      }
    }
  }

  // ============= Profiling the batch compose ==============
  END_EXPERIMENT();
  PRINT_RESULTS(); // uncomment this line to see the actual time elapsed in each function.
  
  // ! And it looks like optimized CRT and NTT saves you dozens of milliseconds. Maybe not worth it. 
  // The problem is: doing many decomposition at a time requires some memory allocations.
}


void PirTest::test_fast_expand_query() {
  print_func_name(__FUNCTION__);
  
  // In this test, I want to make sure if the fast_expand_query is working as expected.
  // There are two ways to order the even and odd parts of a polynomial in the expanding process.
  // One way (the normal way) is to put the even part in it's own location, and the odd part is shifted by expansion tree level size.
  // The other way (the fast way) is to put the even part in 2b and the odd part in 2b + 1.
  // Both of them expand like a binary tree, but the order of the resulting polynomial is different.
  // Here is the access pattern of the normal expansion: https://raw.githubusercontent.com/chenyue42/images-for-notes/master/uPic/expansion.png
  // And the fast expansion will look like a noremal binary tree.


  PirParams pir_params;
  auto params = pir_params.get_seal_params();
  auto context_ = pir_params.get_context();
  auto evaluator_ = seal::Evaluator(context_);
  auto keygen_ = seal::KeyGenerator(context_);
  auto secret_key_ = keygen_.secret_key();
  auto encryptor_ = new seal::Encryptor(context_, secret_key_);
  auto decryptor_ = new seal::Decryptor(context_, secret_key_);  
  const size_t coeff_count = DatabaseConstants::PolyDegree;
  std::stringstream query_stream;
  const size_t fst_dim_sz = DatabaseConstants::MaxFstDimSz;
  const size_t useful_cnt = pir_params.get_fst_dim_sz() + pir_params.get_l() * (pir_params.get_dims().size() - 1);

  PirClient client(pir_params);
  PirServer server(pir_params);
  const size_t client_id = client.get_client_id();

  // ============= setup the server ==============
  std::stringstream galois_key_stream, gsw_stream, data_stream;
  // Client create galois keys and gsw keys and writes to the stream (to the
  // server)
  client.create_galois_keys(galois_key_stream);
  client.write_gsw_to_stream(client.generate_gsw_from_key(), gsw_stream);
  //--------------------------------------------------------------------------------
  // Server receives the gsw keys and galois keys and loads them when needed
  server.set_client_galois_key(client_id, galois_key_stream);
  server.set_client_gsw_key(client_id, gsw_stream);

  // ============= Generate the query ==============
  const size_t query_idx = 253;
  const size_t reversed_idx = utils::bit_reverse(query_idx, pir_params.get_expan_height());
  seal::Ciphertext normal_seeded = client.generate_query(query_idx);
  seal::Ciphertext fast_seeded = client.fast_generate_query(query_idx);

  // ============= Serialize the query ==============
  seal::Ciphertext normal_query, fast_query;
  normal_seeded.save(query_stream);
  normal_query.load(context_, query_stream);
  fast_seeded.save(query_stream);
  fast_query.load(context_, query_stream);
  client.test_budget(normal_query);
  client.test_budget(fast_query);
  // decrypt the query and print it
  auto normal_decrypted = client.decrypt_result({normal_query});
  auto fast_decrypted = client.decrypt_result({fast_query});
  BENCH_PRINT("raw packed query: " << normal_decrypted[0].to_string());
  BENCH_PRINT("fast packed query: " << fast_decrypted[0].to_string());
  PRINT_BAR;

  // ============= Expand the query ==============
  auto normal_exp_q = server.expand_query(client_id, normal_query);
  auto fast_exp_q = server.fast_expand_qry(client_id, fast_query);

  client.test_budget(normal_exp_q[query_idx % fst_dim_sz]);
  client.test_budget(fast_exp_q[query_idx % fst_dim_sz]);

  normal_decrypted = client.decrypt_result(normal_exp_q);
  fast_decrypted = client.decrypt_result(fast_exp_q);
  BENCH_PRINT("normal Expanded query: " << normal_decrypted[query_idx % fst_dim_sz].to_string());
  BENCH_PRINT("fast Expanded query: " << fast_decrypted[query_idx % fst_dim_sz].to_string());

  for (size_t i = 0; i < useful_cnt; i++) {
    // BENCH_PRINT(fast_decrypted[i].to_string());
    if (fast_decrypted[i].to_string() == "1") {
      BENCH_PRINT("equals one at idx: " << i << " ");
    }
  }
}
