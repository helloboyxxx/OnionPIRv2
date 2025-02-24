#pragma once
#include <stdint.h>
#include <stddef.h>

// define a structure for a matrix
typedef struct {
    uint64_t *data;
    size_t rows;
    size_t cols;
    size_t levels;
} matrix_t; 

// performing coeff_val_cnt many matrix matrix multiplications. Assumes that the B matrix has only two columns.
void level_mat_mult(matrix_t *A, matrix_t *B, matrix_t *out);

// Perform the Matrix Multiplication over a direct product over component wise vector multiplication.
void component_wise_mult(matrix_t *A, matrix_t *B, matrix_t *out);