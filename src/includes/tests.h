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
void ideal_fst_dim();

// matrix tests
void test_matrix_mult();
void level_mat_mult_demo();
void component_wise_mult_demo();