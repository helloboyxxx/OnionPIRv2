Here I am listing the basic pipeline of the current OnionPIRv2. I will focus on the places where NTT transforms happen (so that we can try to remove them). 



`test_pir()`in `test.cpp`: 

1. Setup the pir parameters using constants defined in `database_constants.h`

2. server generate data in a specific order. Each entry has its own index.  (`server.gend_data()`)

3. Initiate a new client, create client's Galois keys (key switching keys) and $\mathsf{RGSW}(s)$. These two are sent to the server and stored. 

4. ###### Online phase: client generate the query BFV ciphertext. This BFV is in coefficient form (`client.generate_query`).

5. `server.make_query()` first expand_query and reconstruct query vectors for different dimensions. 

6. Transform the first dimension query vector to NTT form and evaluate the first dimension. 

7. After the matrix-vector multiplication, we transform the resulting database to coefficient form. 

8. For following dimensions, we run `evaluate_gsw_product`, where we do the $\mathsf{RGSW}(b) * (y - x) + x$ trick. The $*$ is `exteral_product`.

   - Inside `external_product`, we first **decompose the $y-x$ BFV ciphertext** in coefficient form and **transform it to NTT form**. Then we perform the matrix-matrix multiplication between the decomposed BFV and the GSW selection vector.

   - After the external product, we transform the resulting BFV ciphertexts to coefficients form again.

9. Finally, we will get a single result in the coefficient form, which will be returned to the client for decryption.



The main issue is at step 8. If we know how to perform the gadget decomposition in NTT form, we should be able to remove the "NTT back and forth" computations. In theory, the decomposition is some divisions of different powers of the "base" $B$. In our current solution, since we are using coefficient form, we compute divisions by right-shifting each coefficient by $p \cdot \log B$, where $1 \leq p \leq l$ is the power of the base. As for the polynomials in NTT form, I have no idea. Especially when they are in RNS.

When using single coefficient modulus, my best guess is that we can perform constant multiplications in NTT form to replace the right shifts. Other than this, I currently don't have a solution for this. 

Without these messy transformations, a quick experiment result shows that the "other dimensions" can be 25%-40% faster. (Not a huge improvement overall).

