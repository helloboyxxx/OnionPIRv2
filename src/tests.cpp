#include "tests.h"
#include "gsw_eval.h"
#include "pir.h"
#include "server.h"
#include "client.h"
#include "utils.h"
#include <cassert>
#include <iostream>
#include <bitset>



#define EXPERIMENT_ITERATIONS 3

void print_func_name(std::string func_name) {
#ifdef _DEBUG
  std::cout << "                    "<< func_name << "(Debug build)" << std::endl;
#endif
#ifdef _BENCHMARK
  std::cout << "                    "<< func_name << "(Benchmark build)" << std::endl;
#endif
}

void run_tests() {
  DEBUG_PRINT("Running tests");
  PRINT_BAR;
  test_pir();
  PRINT_BAR;
  DEBUG_PRINT("Tests finished");
}


void test_pir() {
  print_func_name(__FUNCTION__);
  auto success_count = 0;
  
  // ============== setting parameters for PIR scheme ==============
  PirParams pir_params;
  pir_params.print_params();
  PirServer server(pir_params); // Initialize the server with the parameters
  
  BENCH_PRINT("Initializing server...");
  // Data to be stored in the database.
  server.gen_data();
  BENCH_PRINT("Server initialized");

  // some global results
  size_t galois_key_size = 0;
  size_t gsw_key_size = 0;
  size_t query_size = 0;
  size_t response_size = 0;

  // Run the query process many times.
  srand(time(0)); // reset the seed for the random number generator
  for (size_t i = 0; i < EXPERIMENT_ITERATIONS; i++) {
    
    // ============= OFFLINE PHASE ==============
    // Initialize the client
    PirClient client(pir_params);
    const size_t client_id = client.get_client_id();
    std::stringstream galois_key_stream, gsw_stream, data_stream;

    // Client create galois keys and gsw keys and writes to the stream (to the server)
    galois_key_size = client.create_galois_keys(galois_key_stream);
    gsw_key_size = client.write_gsw_to_stream(
        client.generate_gsw_from_key(), gsw_stream);
    //--------------------------------------------------------------------------------
    // Server receives the gsw keys and galois keys and loads them when needed
    server.set_client_galois_key(client_id, galois_key_stream);
    server.set_client_gsw_key(client_id, gsw_stream);

    // ===================== ONLINE PHASE =====================
    // Client start generating query
    size_t entry_index = rand() % pir_params.get_num_entries();
  
    // ============= CLIENT ===============
    TIME_START(CLIENT_TOT_TIME);
    PirQuery query = client.generate_query(entry_index);
    query_size = client.write_query_to_stream(query, data_stream);
    TIME_END(CLIENT_TOT_TIME);
    
    // ============= SERVER ===============
    TIME_START(SERVER_TOT_TIME);
    auto result = server.make_query(client_id, data_stream);
    TIME_END(SERVER_TOT_TIME);

    // ============= CLIENT ===============
    TIME_START(CLIENT_TOT_TIME);
    // client gets result from the server and decrypts it
    auto decrypted_result = client.decrypt_result(result);
    Entry result_entry = client.get_entry_from_plaintext(entry_index, decrypted_result[0]);
    TIME_END(CLIENT_TOT_TIME);

    // write the result to the stream to test the size
    std::stringstream result_stream;
    response_size = result[0].save(result_stream);
    result_stream.str(std::string()); // clear the stream

    // Directly get the plaintext from server. Not part of PIR.
    Entry actual_entry = server.direct_get_entry(entry_index);
    // extract and print the actual entry index
    uint64_t actual_entry_idx = get_entry_idx(actual_entry);
    uint64_t result_entry_idx = get_entry_idx(result_entry);
    
    END_EXPERIMENT();
    // ============= PRINTING RESULTS ===============    
    DEBUG_PRINT("\t\tWanted/result/actual idx:\t" << entry_index << " / " << result_entry_idx << " / " << actual_entry_idx);
    #ifdef _DEBUG
    PRINT_RESULTS(i+1);
    #endif

    if (entry_is_equal(result_entry, actual_entry)) {
      // print a green success message
      std::cout << "\033[1;32mSuccess!\033[0m" << std::endl;
      success_count++;
    } else {
      // print a red failure message
      std::cout << "\033[1;31mFailure!\033[0m" << std::endl;
      std::cout << "PIR Result:\t";
      print_entry(result_entry);
      std::cout << "Actual Entry:\t";
      print_entry(actual_entry);
    }
    PRINT_BAR;
  }

  double avg_server_time = GET_AVG_TIME(SERVER_TOT_TIME);
  double throughput = pir_params.get_DBSize_MB() / (avg_server_time / 1000);
  

  // ============= PRINTING FINAL RESULTS ===============]
  BENCH_PRINT("                       FINAL RESULTS                         ")
  PRINT_BAR;
  BENCH_PRINT("Success rate: " << success_count << "/" << EXPERIMENT_ITERATIONS);
  BENCH_PRINT("galois key size: " << galois_key_size << " bytes");
  BENCH_PRINT("gsw key size: " << gsw_key_size << " bytes");
  BENCH_PRINT("total key size: " << static_cast<double>(galois_key_size + gsw_key_size) / 1024 / 1024 << "MB");
  BENCH_PRINT("query size: " << query_size << " bytes");
  BENCH_PRINT("response size: " << response_size << " bytes");
  
  // PRINT_AVERAGE_RESULTS();
  PRETTY_PRINT();
  BENCH_PRINT("Throughput: " << throughput << " MB/s");
}

// TODO: I deleted some small test cases and examples because I am refactoring the code. I will add them back later.
