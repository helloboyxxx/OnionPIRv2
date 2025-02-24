#include "matrix.h"
#include "seal/util/clang.h"
#include <cstring>

void level_mat_mult(matrix_t *A, matrix_t *B, matrix_t *out) {
  const size_t rows = A->rows; // multiple of 32
  const size_t cols = A->cols; // multiple of 32
  const size_t p = B->cols; // p=2 (assumed)
  const size_t levels = A->levels;
  const uint64_t *A_data = A->data;
  const uint64_t *B_data = B->data;
  uint64_t *out_data = out->data;

  // We always assume p=2. Because the BFV ciphertext has two polynomials. 
  // This assumption keeps the code simple.
  if (p != 2) { return; } 

  // define pointers
  const uint64_t *A_ptr;
  const uint64_t *B_ptr;
  uint64_t *C_ptr;
  uint64_t db0, db1, db2, db3, db4, db5, db6, db7;
  uint64_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  uint64_t tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15;
  size_t i, j, level; 

  // For each "level," we do one standard mat-mat multiplication.
  // A(level) is m-by-n, B(level) is n-by-2, out(level) is m-by-2
  for (level = 0; level < levels; ++level) {
    // Offsets into the flat arrays for this level
    A_ptr = A_data + level * (rows * cols);
    B_ptr = B_data + level * (cols * p);
    C_ptr = out_data + level * (rows * p);

    // Then we can compute a normal matrix multiplication
    // This is a slight variation of the 
    for (i = 0; i < rows; i += 8) {
      tmp0 = 0; tmp1 = 0; tmp2 = 0; tmp3 = 0;
      tmp4 = 0; tmp5 = 0; tmp6 = 0; tmp7 = 0;
      tmp8 = 0; tmp9 = 0; tmp10 = 0; tmp11 = 0;
      tmp12 = 0; tmp13 = 0; tmp14 = 0; tmp15 = 0;
      for (j = 0; j < cols; j++) {
        db0 = A_ptr[i * cols + j];
        db1 = A_ptr[(i + 1) * cols + j];
        db2 = A_ptr[(i + 2) * cols + j];
        db3 = A_ptr[(i + 3) * cols + j];
        db4 = A_ptr[(i + 4) * cols + j];
        db5 = A_ptr[(i + 5) * cols + j];
        db6 = A_ptr[(i + 6) * cols + j];
        db7 = A_ptr[(i + 7) * cols + j];
        tmp0 += db0 * B_ptr[j * p];
        tmp1 += db0 * B_ptr[j * p + 1];
        tmp2 += db1 * B_ptr[j * p];
        tmp3 += db1 * B_ptr[j * p + 1];
        tmp4 += db2 * B_ptr[j * p];
        tmp5 += db2 * B_ptr[j * p + 1];
        tmp6 += db3 * B_ptr[j * p];
        tmp7 += db3 * B_ptr[j * p + 1];
        tmp8 += db4 * B_ptr[j * p];
        tmp9 += db4 * B_ptr[j * p + 1];
        tmp10 += db5 * B_ptr[j * p];
        tmp11 += db5 * B_ptr[j * p + 1];
        tmp12 += db6 * B_ptr[j * p];
        tmp13 += db6 * B_ptr[j * p + 1];
        tmp14 += db7 * B_ptr[j * p];
        tmp15 += db7 * B_ptr[j * p + 1];
      }
      C_ptr[i * p] = tmp0;
      C_ptr[i * p + 1] = tmp1;
      C_ptr[(i + 1) * p] = tmp2;
      C_ptr[(i + 1) * p + 1] = tmp3;
      C_ptr[(i + 2) * p] = tmp4;
      C_ptr[(i + 2) * p + 1] = tmp5;
      C_ptr[(i + 3) * p] = tmp6;
      C_ptr[(i + 3) * p + 1] = tmp7;
      C_ptr[(i + 4) * p] = tmp8;
      C_ptr[(i + 4) * p + 1] = tmp9;
      C_ptr[(i + 5) * p] = tmp10;
      C_ptr[(i + 5) * p + 1] = tmp11;
      C_ptr[(i + 6) * p] = tmp12;
      C_ptr[(i + 6) * p + 1] = tmp13;
      C_ptr[(i + 7) * p] = tmp14;
      C_ptr[(i + 7) * p + 1] = tmp15;
    }
  } // end for(level)
}


void level_mat_mult_128(matrix_t *A, matrix_t *B, matrix128_t *out) {
  const size_t rows = A->rows; // multiple of 32
  const size_t cols = A->cols; // multiple of 32
  const size_t p = B->cols; // p=2 (assumed)
  const size_t levels = A->levels;
  const uint64_t *A_data = A->data;
  const uint64_t *B_data = B->data;
  uint128_t *out_data = out->data;

  // We always assume p=2. Because the BFV ciphertext has two polynomials. 
  // This assumption keeps the code simple.
  if (p != 2) { return; } 

  // define pointers
  const uint64_t *A_ptr;
  const uint64_t *B_ptr;
  uint128_t *C_ptr;
  uint128_t db0, db1, db2, db3, db4, db5, db6, db7;
  uint128_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  uint128_t tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15;
  size_t i, j, level; 

  // For each "level," we do one standard mat-mat multiplication.
  // A(level) is m-by-n, B(level) is n-by-2, out(level) is m-by-2
  for (level = 0; level < levels; ++level) {
    // Offsets into the flat arrays for this level
    A_ptr = A_data + level * (rows * cols);
    B_ptr = B_data + level * (cols * p);
    C_ptr = out_data + level * (rows * p);

    // Then we can compute a normal matrix multiplication
    // This is a slight variation of the 
    for (i = 0; i < rows; i += 8) {
      tmp0 = 0; tmp1 = 0; tmp2 = 0; tmp3 = 0;
      tmp4 = 0; tmp5 = 0; tmp6 = 0; tmp7 = 0;
      tmp8 = 0; tmp9 = 0; tmp10 = 0; tmp11 = 0;
      tmp12 = 0; tmp13 = 0; tmp14 = 0; tmp15 = 0;
      for (j = 0; j < cols; j++) {
        db0 = A_ptr[i * cols + j];
        db1 = A_ptr[(i + 1) * cols + j];
        db2 = A_ptr[(i + 2) * cols + j];
        db3 = A_ptr[(i + 3) * cols + j];
        db4 = A_ptr[(i + 4) * cols + j];
        db5 = A_ptr[(i + 5) * cols + j];
        db6 = A_ptr[(i + 6) * cols + j];
        db7 = A_ptr[(i + 7) * cols + j];
        tmp0 += db0 * B_ptr[j * p];
        tmp1 += db0 * B_ptr[j * p + 1];
        tmp2 += db1 * B_ptr[j * p];
        tmp3 += db1 * B_ptr[j * p + 1];
        tmp4 += db2 * B_ptr[j * p];
        tmp5 += db2 * B_ptr[j * p + 1];
        tmp6 += db3 * B_ptr[j * p];
        tmp7 += db3 * B_ptr[j * p + 1];
        tmp8 += db4 * B_ptr[j * p];
        tmp9 += db4 * B_ptr[j * p + 1];
        tmp10 += db5 * B_ptr[j * p];
        tmp11 += db5 * B_ptr[j * p + 1];
        tmp12 += db6 * B_ptr[j * p];
        tmp13 += db6 * B_ptr[j * p + 1];
        tmp14 += db7 * B_ptr[j * p];
        tmp15 += db7 * B_ptr[j * p + 1];
      }
      C_ptr[i * p] = tmp0;
      C_ptr[i * p + 1] = tmp1;
      C_ptr[(i + 1) * p] = tmp2;
      C_ptr[(i + 1) * p + 1] = tmp3;
      C_ptr[(i + 2) * p] = tmp4;
      C_ptr[(i + 2) * p + 1] = tmp5;
      C_ptr[(i + 3) * p] = tmp6;
      C_ptr[(i + 3) * p + 1] = tmp7;
      C_ptr[(i + 4) * p] = tmp8;
      C_ptr[(i + 4) * p + 1] = tmp9;
      C_ptr[(i + 5) * p] = tmp10;
      C_ptr[(i + 5) * p + 1] = tmp11;
      C_ptr[(i + 6) * p] = tmp12;
      C_ptr[(i + 6) * p + 1] = tmp13;
      C_ptr[(i + 7) * p] = tmp14;
      C_ptr[(i + 7) * p + 1] = tmp15;
    }
  } // end for(level)
}



void component_wise_mult(matrix_t *A, matrix_t *B, matrix_t *out) {
  const size_t m = A->rows; // multiple of 32
  const size_t n = A->cols; // multiple of 32
  const size_t p = B->cols; // p=2 (assumed)
  const size_t levels = A->levels;
  uint64_t *A_data = A->data;
  uint64_t *B_data = B->data;
  uint64_t *out_data = out->data;
  // Safety check (not strictly necessary, but wise):
  if (p != 2) { return; }  
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      uint64_t *db_ptr = A_data + (i * n + j) * levels;
      uint64_t *q0 = B_data + j * 2 * levels;
      uint64_t *q1 = q0 + levels;
      uint64_t *out_0 = out_data + i * 2 * levels;
      uint64_t *out_1 = out_0 + levels;
      #pragma GCC unroll 32
      for (size_t level = 0; level < levels; ++level) {
        out_0[level] += db_ptr[level] * q0[level];
        out_1[level] += db_ptr[level] * q1[level];
      }
    }
  }
}