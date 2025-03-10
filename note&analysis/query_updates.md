# Updates on query related functionalities and details

The code was different from the pseudocode in OnionPIR paper. Yue made some changes on the code so to align with the pseudocode logic.



### Updates on `generate_query` on client side.

- Previously, if the `query_indices[i] == 0`for dimension $i$, the code packed some "special values" to the coefficients of the query. Correspondingly, in `other_dim_mux` on the server side, if the selection vector is RGSW(0), then it outputs the second half of the given vector. These old code are not consistent with the output of `get_query_indices`, and are againsts the conventional vector order. Hence, the first update is to change the code so that we indeed pack the value 1 when `query_indices[i] == 1`.
- TODO: also encrypt the first $l$ rows for RGSW queries. Worth it? The current `query_to_gsw` takes about 60ms for each dimension. 
  - Looks like this is impossible. Long story in short: it is impossible to put the complete secret key $s$ inside a coefficient..




### Question on exsisting `generate_query` code:

Only $1 / 3$ of the RGSW gadget are used, others are 0. In the following example, $B = 256 = \text{0x100}$. 

<center>
  <figure>
    <img src=" https://raw.githubusercontent.com/helloboyxxx/images-for-notes/master/uPic/image-20240723222155623.png " style="width:50%;" />
    <figcaption> Decrypted RGSW(1) in "problematic" generate_query </figcaption>
  </figure>
</center>

This could be serious because $l$ and $B$ controls the gadget, which affects the ability to decompose and reconstruct values. In the current algorithm,  only with $l \geq 8$ can we guarantee sufficient precision to retrieve the correct data (check the final section for other default values.) 

<center>
  <figure>
    <img src=" https://raw.githubusercontent.com/helloboyxxx/images-for-notes/master/uPic/image-20240723235037356.png" alt="image-20240723235037356 " style="width:40%;" />
    <figcaption> Retrieved result v.s. Actual data </figcaption>
  </figure>
</center>
It is not because the values are not filled in, but is because the first $2/3$ values are all rounded to 0. This is of course weird, but it works now.



### Strange but working code design

There is a reason why we have to "first encrypt the initial plaintext then add RGSW gadgets": RGSW gadget uses the ciphertext modulus for coefficients, which is larger than the plaintext modulus. This is related to how BFV uses these coefficients. One must dive into SEAL to understand how this works. 

One way to learn if the current code is working exactly as we expected is to try to decrypt the GSW ciphertext and check if the gadget looks the same as in the paper.



---

### Default PIR param

```cpp
#define DB_SZ       1 << 15
#define NUM_DIM     8
#define NUM_ENTRIES 1 << 15
#define ENTRY_SZ    12000
#define GSW_L       9
#define GSW_L_KEY   9
```

#### Some expected values:

```
==============================================================
                       PIR PARAMETERS                         
==============================================================
num_entries_											 = 32768
l_																= 9
base_log2_												  = 8
entry_size_												   = 12000
DBSize_ (num plaintexts in database)        = 32768
DBCapacity (max num of entries)      		   = 32768
dimensions_                           						 = [ 256 2 2 2 2 2 2 2 ]
seal_params_.poly_modulus_degree()        = 4096
seal_params_.coeff_modulus().bit_count   = [36 + 36 + 37] bits
seal_params_.coeff_modulus().size()      	  = 3
seal_params_.plain_modulus().bitcount() = 25
==============================================================


plain_modulus = 16777259
inverse of 256 in plain_modulus: 12484640
context_data->parms().coeff_modulus() = {68585185425, 68585013729}
68585185425 * 68585013729 / 16777259
delta = 280374516718114

```



---

### TODO: 

Better profiling method. 

Better unit test.









