Here I am listing the basic pipeline of the current OnionPIRv2. 

`test_pir()`in `test.cpp`: 

1. Setup the pir parameters using constants defined in `database_constants.h`
2. server generate data in a specific order. Each entry has its own index.  (`server.gend_data()`)
3. Initiate a new client, create client's Galois keys (key switching keys) and $\mathsf{RGSW}(s)$. 