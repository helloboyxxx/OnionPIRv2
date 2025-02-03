#pragma once

void run_tests();

// BFV & GSW tests
void bfv_example();
void test_ct_sub();
void test_ntt_add();
void test_ntt_scalar_mul();

void test_external_product(); 

// SEAL serialization
void serialization_example();

// the main test for PIR
void test_pir();


// Other tests
void test_prime_gen();
void test_reading_pt();
void test_reading_ct();
void test_reading_inter();
void test_matrix_mult();
void test_pt_size();