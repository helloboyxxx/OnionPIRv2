#pragma once

void run_tests();

// the main test for PIR
void test_pir();

// BFV & GSW tests
void bfv_example();
void test_ct_sub();
void test_ntt_add();
void test_ntt_scalar_mul();
void test_external_product(); 

// SEAL serialization
void serialization_example();

// Other tests
void test_prime_gen();
void test_reading_pt();
void test_reading_ct();
void test_matrix_mult();
void test_pt_size();
void ideal_fst_dim();