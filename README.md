# OnionPIR

A rewrite of OnionPIR without the NFLlib dependency. 
You need to have SEAL installed. After installation, set `CMAKE_PREFIX_PATH` to the library's location. Separate versions of the library can be used for debugging and benchmarking. To run as a debug build, set -DCMAKE_BUILD_TYPE=Debug as a cmake option. To run benchmarks, set -DCMAKE_BUILD_TYPE=Benchmark. The benchmark build type is used by default.

Build and run the project like this:
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug .. 
make
./Onion-PIR
```

If using clangd, you can use `cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON .` to generate the `compile_commands.json`.

