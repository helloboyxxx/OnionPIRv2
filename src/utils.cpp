#include "utils.h"
#include <fstream>
#include <stdexcept>

void utils::negacyclic_shift_poly_coeffmod(seal::util::ConstCoeffIter poly, size_t coeff_count,
                                           size_t shift, const seal::Modulus &modulus,
                                           seal::util::CoeffIter result) {
  if (shift == 0) {
    set_uint(poly, coeff_count, result);
    return;
  }

  uint64_t index_raw = shift;
  const uint64_t coeff_count_mod_mask = static_cast<uint64_t>(coeff_count) - 1;
  for (size_t i = 0; i < coeff_count; i++, poly++, index_raw++) {
    uint64_t index = index_raw & coeff_count_mod_mask;  // shifted index, possibly wrapping around
    if (!(index_raw & static_cast<uint64_t>(coeff_count)) || !*poly) {
      // for those entries that are not wrapped around
      result[index] = *poly;
    } else {
      // For wrapped around entries, we fill in additive inverse.
      result[index] = modulus.value() - *poly; 
    }
  }
}

void utils::shift_polynomial(seal::EncryptionParameters &params, seal::Ciphertext &encrypted,
                             seal::Ciphertext &destination, size_t index) {
  const auto encrypted_count = encrypted.size();
  const auto coeff_count = params.poly_modulus_degree();
  const auto coeff_mod_count = params.coeff_modulus().size() - 1;
  destination = encrypted;
  for (int i = 0; i < encrypted_count; i++) {
    for (int j = 0; j < coeff_mod_count; j++) {
      negacyclic_shift_poly_coeffmod(encrypted.data(i) + (j * coeff_count), coeff_count, index,
                                     params.coeff_modulus()[j],
                                     destination.data(i) + (j * coeff_count));
    }
  }
}


void negate_poly_inplace(seal::Plaintext &plain) {
  std::cout << "TODO" << std::endl;
}


std::string uint128_to_string(uint128_t value) {
    // Split the 128-bit value into two 64-bit parts
    uint64_t high = value >> 64;
    uint64_t low = static_cast<uint64_t>(value);

    std::ostringstream oss;

    // Print the high part, if it's non-zero, to avoid leading zeros
    if (high != 0) {
        oss << high << " * 2^64 + " << low;
    } else {
        oss << low;
    }
    return oss.str();
}



std::vector<std::vector<uint64_t>> gsw_gadget(size_t l, uint64_t base_log2, size_t coeff_mod_count,
                const std::vector<seal::Modulus> &coeff_modulus) {
  // Create RGSW gadget.
  std::vector<std::vector<uint64_t>> gadget(coeff_mod_count, std::vector<uint64_t>(l));
  for (int i = 0; i < coeff_mod_count; i++) {
    const uint128_t mod = coeff_modulus[i].value();
    uint128_t pow = 1;
    for (int j = l - 1; j >= 0; j--) {
      gadget[i][j] = pow;
      pow = (pow << base_log2) % mod;
    }
  }

// #ifdef _DEBUG
//   for (int i = 0; i < coeff_mod_count; i++) {
//     std::cout << "Gadget for mod " << i << ": " << std::endl;
//     for (int j = 0; j < l; j++) {
//       std::cout << uint128_to_string(gadget[i][j]) << " ";
//     }
//     std::cout << std::endl;
//   }
// #endif

  return gadget;
}


/**
 * @brief Generate the smallest prime that is at least bit_width bits long.
 * Optimized for repeated calls with the same bit_width by caching results in a file.
 * @param bit_width >= 2 and <= 64
 * @return std::uint64_t  
 */
std::uint64_t generate_prime(size_t bit_width) {
  if (bit_width < 2) throw std::invalid_argument("Bit width must be at least 2.");

  // ================= Read from file if it exists
  std::ifstream file("prime_cache.txt");
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      std::stringstream ss(line);
      std::string token;
      std::getline(ss, token, ',');
      size_t cached_bit_width = std::stoul(token);
      if (cached_bit_width == bit_width) {
        std::getline(ss, token, ',');
        return std::stoull(token);
      }
    }
    file.close();
  }

  // Otherwise, generate a new prime
  std::uint64_t candidate = 1ULL << (bit_width - 1);
  do {
      candidate++;
      // Ensure candidate is odd, as even numbers greater than 2 cannot be prime
      candidate |= 1;
  } while (!seal::util::is_prime(seal::Modulus(candidate)));

  // write the bit_width, candidate pair to csv file
  std::string file_name = "prime_cache.txt";
  std::ofstream out_file(file_name);
  if (out_file.is_open()) {
    out_file << bit_width << "," << candidate << std::endl;
    out_file.close();
  } else {
    std::cerr << "Failed to create file " << file_name << std::endl;
  }
  
  return candidate;
}

// converting a uint64_t to a std::vector<uint8_t> of size 8. Assyming the input vector has at least 8 elements.
void writeIdxToEntry(const uint64_t idx, Entry &entry) {
  // Convert id to bytes and write them to the start of the entry.
  for (int i = 0; i < 8; ++i) {
    // Extract the i-th byte from the least significant to the most significant
    entry[7 - i] = static_cast<uint8_t>((idx >> (i * 8)) & 0xFF);
  }
}

uint64_t get_entry_idx(const Entry &entry) {
  uint64_t idx = 0;
  for (int i = 0; i < 8; ++i) {
    idx |= static_cast<uint64_t>(entry[7 - i]) << (i * 8);
  }
  return idx;
}


void print_entry(const Entry &entry) {
  int cnt = 0;
  for (auto &val : entry) {
    if (cnt < 10) { 
      std::cout << (int)val << ", ";
    }
    cnt += 1;
  }
  std::cout << std::endl;
}


bool entry_is_equal(const Entry &entry1, const Entry &entry2) {
  for (size_t i = 0; i < entry1.size(); i++) {
    if (entry1[i] != entry2[i]) {
      std::cerr << "Entries are not equal" << std::endl;
      return false;
    }
  }
  std::cout << "Entries are equal" << std::endl;
  return true;
}

size_t poly_idx_to_actual(const size_t poly_idx, const size_t fst_dim_sz, const size_t other_dim_sz) {
  // we remap (N/N_1) * k + j to k + j * N_1. 
  // Then k + j * N_1 is again remapped to (N/N_1)*k * (k/t) + j*t + k mod t
  size_t k = poly_idx / other_dim_sz;
  size_t j = poly_idx % other_dim_sz;

  auto ts = DatabaseConstants::TileSz;
  auto max_mult_of_ts = (fst_dim_sz / ts) * ts;
  auto remain_fst_dim_sz = fst_dim_sz % ts;
  if (k < max_mult_of_ts) {
    return other_dim_sz * ts * (k / ts) + j * ts + k % ts;
  } else {
    return other_dim_sz * max_mult_of_ts + j * remain_fst_dim_sz + k - max_mult_of_ts;
  }

  // return k + j * fst_dim_sz;
}

void print_progress(size_t current, size_t total) {
  float progress = static_cast<float>(current) / total;
  int bar_width = 70;
  std::cout << "[";
  int pos = static_cast<int>(bar_width * progress);
  for (int i = 0; i < bar_width; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}


Entry generate_entry(const uint64_t entry_id, const size_t entry_size, std::ifstream &random_file) {
  Entry entry(entry_size);

  // write the entry_id to the first 8 bytes of the entry
  writeIdxToEntry(entry_id, entry);

  // fill the rest of the entry with random bytes
  // random_file.read(reinterpret_cast<char *>(entry.data()), entry_size);
  random_file.read(reinterpret_cast<char *>(entry.data() + 8), entry_size - 8);
  return entry;
}


size_t next_pow_of_2(const size_t n) {
  size_t p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

size_t roundup_div(const size_t numerator, const size_t denominator) {
  if (denominator == 0) {
    throw std::invalid_argument("roundup_div division by zero");
  }
  return (numerator + denominator - 1) / denominator;
}





// temp test
uint64_t utils::global_sum = 0; // Definition and initialization
void utils::print_sum() {
  std::cout << "global_sum: " << utils::global_sum << std::endl;
}