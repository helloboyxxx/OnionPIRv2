# OnionV2

### Huge performance improvements

#### First dimension

- Not doing splitting for the first dimension. 
  - Pro: about 2x speedup in the first dimension. 
  - Con: larger noise growth, must use larger GSW_L and smaller plaintext mods. $60/124 \to 49 / 120$.
- Delayed modulus. 
  - Pro: only one mod required for each result of the first dimension. Hence, $O(N/N_1)$ instead of $O(N)$. About 2x speedup. 
  - Con: addition may brings overflow? Say all plaintext in the database is $2^{\log t -1}$. Multiplied by random int $\in [2^{\log q} - 1]$ can be large. After doing $N_1$ many multiplications and additions, what's the probability of overflowing? Take $q = 60, t = 17$, both packed to `uint64_t`.

#### Flags:

- -O3 flag: 
  - 250MB/s v.s.1000MB/s on Mac M1.
- AVX flag for SEAL. It makes NTT faster.

### Smaller things

- Barrett reduction for calculating modulus (about 10% speedup?)
- AVX 
- GSW $l$ separation.
  - Expansion uses a $\mathsf{RGSW}(s)$. For that, we use large $l = 9$ or $l = 15$. Canceling the noise growth brought by not doing the splitting. 
  - Retrieving stage we use $l = 4$ or $l = 5$. 
- Cache performance for the first dimension. (about 70 MB/s throughput speedup?)
- PRG for galois keys and $\mathsf{RGSW}(s)$. Cut half the storage size.

### What's next 

Modulus switching for a single modulus.

Abusing AVX-512? 

- Can try something very quick: `_mm512_mullo_epi32` for `{30, 30, 60}`

NTT form polynomial addition

memory direct scan 

