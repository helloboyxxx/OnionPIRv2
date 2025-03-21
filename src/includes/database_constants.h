#pragma once
#include <cstddef>
#include <array>

namespace DatabaseConstants {
  // Currently, if the degree is 4096, 256 for the first dimension looks optimal. 
  // If the degree is 2048, 512 for the first dimension looks optimal.
  constexpr size_t MaxFstDimSz = 256; // Maximum size of the first dimension. Actual size can only be smaller.
  // constexpr size_t MaxFstDimSz = 512;

  // ============================================================================================
  // ! THE FOLLOWING FEW CHOICES ARE FOR POLYNOMIAL DEGREE = 4096
  // ============================================================================================

  // ! ========================== 2^19 * 24KB = 12GB (bvest throughput) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 19;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 9;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q

  // ! ========================== 2^16 * 24KB = 1.5GB (normal test) ==========================
  constexpr size_t PolyDegree = 4096;
  constexpr size_t NumEntries = 1 << 16;            // number of entries in the database
  constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  constexpr size_t GSW_L = 4;                       // parameter for GSW scheme
  constexpr size_t GSW_L_KEY = 9;                   // GSW for query expansion
  constexpr size_t PlainMod = 49;
  constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q

  // ! ========================== 2^15 * 24KB = 768MB (for Mac M1 8GB RAM ˚∆˚) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 15;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 4;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 9;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q

  // ! ========================== 2^10 * 24KB = 24MB (quick test) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 10;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 10;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q


// ! ========================== 2^8 * 24KB = 6MB (first dimension size test) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = MaxFstDimSz;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 10;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 30;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q


  // ! ========================== 2^23 * 1KB = 8GB (general large test case when n = 4096) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 23;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 1024;                // 1KB
  // constexpr size_t GSW_L = 4;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 9;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q

  // ! ========================== 2^20 * 1KB = 1GB (general normal test case when n = 4096) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 20;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 1024;                // 1KB
  // constexpr size_t GSW_L = 4;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 9;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q

  // ! ========================== 2^19 * 1KB = 512MB (general small test case when n = 4096) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 19;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 1024;                // 1KB
  // constexpr size_t GSW_L = 4;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 9;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q

  // ! ========================== 2^10 * 1KB = 1MB (general tiny test case when n = 4096) ==========================
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 10;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 1024;                // 1KB
  // constexpr size_t GSW_L = 4;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 9;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q

  // ============================================================================================
  // ! THE FOLLOWING FEW CHOICES ARE FOR POLYNOMIAL DEGREE = 2048
  // ============================================================================================

  // ! ========================== 2^21 * 4KB = 8GB (best throughput) ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 21;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                   // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 2> CoeffMods = {60, 60}; // log q = 60.


  // ! ========================== 2^18 * 4KB = 1GB ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 18;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                   // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 2> CoeffMods = {60, 60}; // log q = 60.


  // ! ========================== 2^16 * 4KB = 256MB (quick test) ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 16;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                   // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 2> CoeffMods = {60, 60}; // log q = 60.

  // ! ========================== 2^23 * 1KB = 8GB (general test case when n = 2048) ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 23;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 1024;                // 1KB
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                   // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 2> CoeffMods = {60, 60}; // log q = 60.

  // ! ========================== 2^20 * 1KB = 1GB (general test case when n = 2048) ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 20;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 1024;                // 1KB
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                   // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 2> CoeffMods = {60, 60}; // log q = 60.

  // ! ========================== 2^22 * 256B = 1GB (general test case when n = 2048) ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 22;             // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 256;                  // 256B
  // constexpr size_t GSW_L = 5;                        // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                   // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 2> CoeffMods = {60, 60}; // log q = 60.


  // ! ========================== 2^21 * 256B = 512MB (general small test case when n = 2048) ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 21;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 256;                 // 256B
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                  // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 2> CoeffMods = {60, 60}; // log q = 60.


  // ======================================== SPECIAL TEST CASES ====================================================
  // ! ========================== Test many small mods ==========================
  // constexpr size_t PolyDegree = 2048;
  // constexpr size_t NumEntries = 1 << 16;            // number of entries in the database
  // constexpr size_t EntrySize = 0;                   // 0 means calculated automatically. Take the largest possible value.
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 20;                   // GSW for query expansion
  // constexpr size_t PlainMod = 17;
  // constexpr std::array<size_t, 4> CoeffMods = {20, 20, 20, 20}; // log q = 60.

  // ! ========================== 2^20 * 1KB = 1GB (special test case when n = 4096) ==========================
  // ! Keep the GSW_L same as the small test case
  // constexpr size_t PolyDegree = 4096;
  // constexpr size_t NumEntries = 1 << 20;            // number of entries in the database. Will be padded to multiples of other dimension size.
  // constexpr size_t EntrySize = 1024;                // 1KB
  // constexpr size_t GSW_L = 5;                       // parameter for GSW scheme
  // constexpr size_t GSW_L_KEY = 15;                   // GSW for query expansion
  // constexpr size_t PlainMod = 49;
  // constexpr std::array<size_t, 3> CoeffMods = { 60, 60, 60 }; // the first two addes up to log q








} // namespace DatabaseConstants