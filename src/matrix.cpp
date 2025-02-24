#include "matrix.h"
#include <cstring>

void level_mat_mult(matrix_t *A, matrix_t *B, matrix_t *out) {
  const size_t m = A->rows; // multiple of 32
  const size_t n = A->cols; // multiple of 32
  const size_t p = B->cols; // p=2 (assumed)
  const size_t levels = A->levels;
  const uint64_t *A_data = A->data;
  const uint64_t *B_data = B->data;
  uint64_t *out_data = out->data;
  // Safety check (not strictly necessary, but wise):
  if (p != 2) {
    // Handle error or return
    return;
  }

  // For each "level," we do one standard matrix multiplication.
  // A(level) is m-by-n, B(level) is n-by-2, out(level) is m-by-2
  for (size_t level = 0; level < levels; ++level) {
    // Offsets into the flat arrays for this level
    const uint64_t *A_ptr = A_data + level * (m * n);
    const uint64_t *B_ptr = B_data + level * (n * p);
    uint64_t *C_ptr = out_data + level * (m * p);
    uint64_t db0, db1, db2, db3;
    uint64_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    for (size_t i = 0; i < m; i += 4) {
      tmp0 = 0; tmp1 = 0; tmp2 = 0; tmp3 = 0;
      tmp4 = 0; tmp5 = 0; tmp6 = 0; tmp7 = 0;
      for (size_t j = 0; j < n; j++) {
        db0 = A_ptr[i * n + j];
        db1 = A_ptr[(i + 1) * n + j];
        db2 = A_ptr[(i + 2) * n + j];
        db3 = A_ptr[(i + 3) * n + j];
        tmp0 += db0 * B_ptr[j * p];
        tmp1 += db0 * B_ptr[j * p + 1];
        tmp2 += db1 * B_ptr[j * p];
        tmp3 += db1 * B_ptr[j * p + 1];
        tmp4 += db2 * B_ptr[j * p];
        tmp5 += db2 * B_ptr[j * p + 1];
        tmp6 += db3 * B_ptr[j * p];
        tmp7 += db3 * B_ptr[j * p + 1];
      }
      C_ptr[i * p] = tmp0;
      C_ptr[i * p + 1] = tmp1;
      C_ptr[(i + 1) * p] = tmp2;
      C_ptr[(i + 1) * p + 1] = tmp3;
      C_ptr[(i + 2) * p] = tmp4;
      C_ptr[(i + 2) * p + 1] = tmp5;
      C_ptr[(i + 3) * p] = tmp6;
      C_ptr[(i + 3) * p + 1] = tmp7;
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
  if (p != 2) {
    // Handle error or return
    return;
  }  
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








