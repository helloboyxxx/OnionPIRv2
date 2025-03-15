#pragma once

void run_tests();

// ! the main test for PIR
void test_pir();

// ======================== BFV & GSW tests ========================
void bfv_example();
void test_external_product(); 

// ======================== SEAL Serialization ========================
void serialization_example();

// ======================== Matrix tests ========================
// test the matrix multiplication performance when using only one level/degree
void test_single_mat_mult();
// simulation of the first dimension multiplication
void test_fst_dim_mult();

// ======================== Other tests ========================
void test_prime_gen();
