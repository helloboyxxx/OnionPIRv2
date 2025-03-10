#pragma once
#include "seal/seal.h"
#include "utils.h"
#include <stdint.h>
#include <stddef.h>

#define ALIGN_BYTES 64  // Typical cache line size

// define a structure for a matrix
typedef struct {
    uint64_t *data;
    size_t rows;
    size_t cols;
    size_t levels;
} matrix_t; 

typedef struct {
    uint128_t *data;
    size_t rows;
    size_t cols;
    size_t levels;
} matrix128_t; 



uint64_t* allocate_aligned_uint64(size_t num_elements);
uint128_t* allocate_aligned_uint128(size_t num_elements);

void naive_mat_mult(matrix_t *A, matrix_t *B, matrix_t *out);

void naive_level_mat_mult(matrix_t *A, matrix_t *B, matrix_t *out);

// performing coeff_val_cnt many matrix matrix multiplications. Assumes that the B matrix has only two columns.
void level_mat_mult(matrix_t *A, matrix_t *B, matrix_t *out);

// suitable for the first dimension. The output is in uint128_t.
void level_mat_mult_128(matrix_t *A, matrix_t *B, matrix128_t *out);

void native_mat_mult_128(matrix_t *A, matrix_t *B, matrix128_t *out);

void naive_level_mat_mult_128(matrix_t *A, matrix_t *B, matrix128_t *out);

void level_mat_mult_direct_mod(matrix_t *A, matrix_t *B, matrix_t *out, const seal::Modulus mod);

// Perform the Matrix Multiplication over a direct product over component wise vector multiplication.
void component_wise_mult(matrix_t *A, matrix_t *B, matrix_t *out);

void component_wise_mult_128(matrix_t *A, matrix_t *B, matrix128_t *out);

// This is using intel::hexl::EltwiseMultMod for each component wise multiplication.
void component_wise_mult_direct_mod(matrix_t *A, matrix_t *B, uint64_t *out, const uint64_t mod);

void level_mat_mult_eigen(matrix_t *A, matrix_t *B, matrix_t *out);

void level_mat_mult_arma(matrix_t *A, matrix_t *B, matrix_t *out);


// calculates c = c + a * b mod m using Barrett reduction
inline void mult_add_mod(uint64_t &a, uint64_t &b, uint64_t &c, const seal::Modulus &m) {
    uint128_t tmp = (uint128_t)a * b + c;
    uint64_t raw[2] = {static_cast<uint64_t>(tmp), static_cast<uint64_t>(tmp >> 64)};
    c = util::barrett_reduce_128(raw, m);
}