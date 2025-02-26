#pragma once

void run_tests();

// the main test for PIR
void test_pir();

// BFV & GSW tests
void bfv_example();
void test_external_product(); 

// SEAL serialization.
void serialization_example();

// Other tests
void test_prime_gen();

// matrix tests
void test_single_mat_mult();
void test_fst_dim_mult();
void level_mat_mult_demo();
void level_mat_mult128_demo();
void component_wise_mult_demo();