Here I am listing the basic pipeline of the current OnionPIRv2. I will focus on the places where NTT transforms happen (so that we can try to remove them). 



`test_pir()`in `test.cpp`: 

1. Setup the pir parameters using constants defined in `database_constants.h`

2. server generate data in a specific order. Each entry has its own index.  (`server.gend_data()`)

3. Initiate a new client, create client's Galois keys (key switching keys) and $\mathsf{RGSW}(s)$. These two are sent to the server and stored. 

4. ###### Online phase: client generate the query BFV ciphertext. This BFV is in coefficient form (`client.generate_query`).

5. `server.make_seeded_query()` first expand_query and reconstruct query vectors for different dimensions. 

6. Transform the first dimension query vector to NTT form and evaluate the first dimension. 

7. After the matrix-vector multiplication, we transform the resulting database to coefficient form. 

8. For following dimensions, we run `evaluate_gsw_product`, where we do the $\mathsf{RGSW}(b) * (y - x) + x$ trick. The $*$ is `exteral_product`.

   - Inside `external_product`, we first **decompose the $y-x$ BFV ciphertext** in coefficient form and **transform it to NTT form**. Then we perform the matrix-matrix multiplication between the decomposed BFV and the GSW selection vector.









