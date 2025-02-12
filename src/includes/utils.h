#pragma once
#include "pir.h"
#include "seal/seal.h"
#include <iostream>
#include <fstream>


#ifdef _DEBUG
#define PRINT_INT_ARRAY(arr_name, arr, size) \
    do {                                     \
        std::cout << arr_name << ": [";      \
        for (int i = 0; i < size; ++i) {     \
            std::cout << arr[i];             \
            if (i < size - 1)                \
                std::cout << ", ";           \
        }                                    \
        std::cout << "]" << std::endl;       \
    } while (0)
#endif

#ifdef _BENCHMARK
#define PRINT_INT_ARRAY(arr_name, arr, size) ;  // do nothing
#endif


template <typename T> std::string to_string(T x) {
  std::string ret;
  if (x == 0) {
    return "0";
  }
  while (x) {
    ret += (x % 10) + '0';
    x /= 10;
  }
  reverse(ret.begin(), ret.end());
  return ret;
}

namespace utils {

  // global sum for testing
extern uint64_t global_sum;

void print_sum();




/*!
    Helper function for multiply_poly_acum. Multiplies two operands together and
   stores the result in product_acum.
*/
inline void multiply_acum(const uint64_t op1, const uint64_t op2, uint128_t &product_acum) {

  /*
  Uncomment these lines to examine the reading performance.
  This shows that the bottle neck is the memory access.
  The memory reading takes about 85% of the actual computation.
  */

  // asm volatile("" : : "r"(op1) : "memory");  // ct
  // asm volatile("" : : "r"(op2) : "memory"); // pt
  // asm volatile("" : : "r"(product_acum) : "memory"); // intermediate result

  // global_sum += op2;  // pt
  
  // The actual computation.
  product_acum = product_acum + static_cast<uint128_t>(op1) * (op2);
}

/*!
    Multiplies two polynomials in NTT form together and adds the result to a
   third polynomial in NTT form.
    @param ct_ptr - Pointer to the start of the data of the first polynomial
    @param pt_ptr - Pointer to the start of the data of the second polynomial
    @param size - Number of polynomial coefficients
    @param result - Pointer to the start of the data of the result polynomial
*/
inline void multiply_poly_acum(const uint64_t* __restrict ct_ptr, const uint64_t* __restrict pt_ptr, const size_t size,
                               uint128_t* __restrict result) {
    for (size_t cc = 0; cc < size; cc += 16) {
        // Perform the actual computation
        multiply_acum(ct_ptr[cc], pt_ptr[cc], result[cc]);
        multiply_acum(ct_ptr[cc + 1], pt_ptr[cc + 1], result[cc + 1]);
        multiply_acum(ct_ptr[cc + 2], pt_ptr[cc + 2], result[cc + 2]);
        multiply_acum(ct_ptr[cc + 3], pt_ptr[cc + 3], result[cc + 3]);
        multiply_acum(ct_ptr[cc + 4], pt_ptr[cc + 4], result[cc + 4]);
        multiply_acum(ct_ptr[cc + 5], pt_ptr[cc + 5], result[cc + 5]);
        multiply_acum(ct_ptr[cc + 6], pt_ptr[cc + 6], result[cc + 6]);
        multiply_acum(ct_ptr[cc + 7], pt_ptr[cc + 7], result[cc + 7]);
        multiply_acum(ct_ptr[cc + 8], pt_ptr[cc + 8], result[cc + 8]);
        multiply_acum(ct_ptr[cc + 9], pt_ptr[cc + 9], result[cc + 9]);
        multiply_acum(ct_ptr[cc + 10], pt_ptr[cc + 10], result[cc + 10]);
        multiply_acum(ct_ptr[cc + 11], pt_ptr[cc + 11], result[cc + 11]);
        multiply_acum(ct_ptr[cc + 12], pt_ptr[cc + 12], result[cc + 12]);
        multiply_acum(ct_ptr[cc + 13], pt_ptr[cc + 13], result[cc + 13]);
        multiply_acum(ct_ptr[cc + 14], pt_ptr[cc + 14], result[cc + 14]);
        multiply_acum(ct_ptr[cc + 15], pt_ptr[cc + 15], result[cc + 15]);
        // multiply_acum(ct_ptr[cc + 16], pt_ptr[cc + 16], result[cc + 16]);
        // multiply_acum(ct_ptr[cc + 17], pt_ptr[cc + 17], result[cc + 17]);
        // multiply_acum(ct_ptr[cc + 18], pt_ptr[cc + 18], result[cc + 18]);
        // multiply_acum(ct_ptr[cc + 19], pt_ptr[cc + 19], result[cc + 19]);
        // multiply_acum(ct_ptr[cc + 20], pt_ptr[cc + 20], result[cc + 20]);
        // multiply_acum(ct_ptr[cc + 21], pt_ptr[cc + 21], result[cc + 21]);
        // multiply_acum(ct_ptr[cc + 22], pt_ptr[cc + 22], result[cc + 22]);
        // multiply_acum(ct_ptr[cc + 23], pt_ptr[cc + 23], result[cc + 23]);
        // multiply_acum(ct_ptr[cc + 24], pt_ptr[cc + 24], result[cc + 24]);
        // multiply_acum(ct_ptr[cc + 25], pt_ptr[cc + 25], result[cc + 25]);
        // multiply_acum(ct_ptr[cc + 26], pt_ptr[cc + 26], result[cc + 26]);
        // multiply_acum(ct_ptr[cc + 27], pt_ptr[cc + 27], result[cc + 27]);
        // multiply_acum(ct_ptr[cc + 28], pt_ptr[cc + 28], result[cc + 28]);
        // multiply_acum(ct_ptr[cc + 29], pt_ptr[cc + 29], result[cc + 29]);
        // multiply_acum(ct_ptr[cc + 30], pt_ptr[cc + 30], result[cc + 30]);
        // multiply_acum(ct_ptr[cc + 31], pt_ptr[cc + 31], result[cc + 31]);
    }
}




void negacyclic_shift_poly_coeffmod(seal::util::ConstCoeffIter poly, size_t coeff_count,
                                    size_t shift, const seal::Modulus &modulus,
                                    seal::util::CoeffIter result);
void shift_polynomial(seal::EncryptionParameters &params, seal::Ciphertext &encrypted,
                      seal::Ciphertext &destination, size_t index);
} // namespace utils


// Inplace calculate the additive inverse of a seal::Plaintext
void negate_poly_inplace(seal::Plaintext &plain);

// Convert a 128-bit unsigned integer to a string
std::string uint128_to_string(uint128_t value);

/**
 * @brief Construct a RGSW gadget. Notice that the gadget is from large to
 * small, i.e., the first row is B^(log q / log B -1), the final row is 1.
 */
std::vector<std::vector<uint64_t>>
gsw_gadget(size_t l, uint64_t base_log2, size_t coeff_mod_count,
           const std::vector<seal::Modulus> &coeff_modulus);


// Generate a prime that is bit_width long
std::uint64_t generate_prime(size_t bit_width);

void writeIdxToEntry(const uint64_t idx, Entry &entry);

// Extract first 8 bytes of the given Entry and format it as a uint64_t and return.
uint64_t get_entry_idx(const Entry &entry);

void print_entry(const Entry &entry);

bool entry_is_equal(const Entry &entry1, const Entry &entry2);


// logical entry index to the actuall index in the database.
// This is a trick we do for reducing the miss cache rate when evaluating the first dimension.
size_t poly_idx_to_actual(const size_t entry_idx, const size_t fst_dim_sz, const size_t other_dim_sz);

void print_progress(size_t current, size_t total);

/**
* @brief Given an entry id and the length of the entry, generate a random entry using random number generator.
* 
* @param entry_id entry id. Not used in implementation, but can be useful if we want to generate entries with specific ids.
* @param len length(size) of the entry. Each entry is a vector of bytes.
* @param random_file random file for quick entry generation.
* @return Entry 
*/
Entry generate_entry(const uint64_t entry_id, const size_t entry_size, std::ifstream &random_file);


size_t next_pow_of_2(const size_t n);

size_t roundup_div(const size_t numerator, const size_t denominator);
