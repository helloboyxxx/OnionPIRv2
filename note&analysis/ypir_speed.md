```
Running YPIR (w/ DoublePIR) on a database of 1073741824 bits, and performing cross-client batching over 1 clients. 
The server performance measurement will be averaged over 5 trials.
Measurement completed. See the README for details on what the following fields mean.
Result:
{
  "offline": {
    "uploadBytes": 0,
    "downloadBytes": 0,
    "serverTimeMs": 6178,
    "clientTimeMs": 0,
    "simplepirPrepTimeMs": 3534,
    "simplepirHintBytes": 29360128,
    "doublepirHintBytes": 14680064
  },
  "online": {
    "uploadBytes": 604160,
    "downloadBytes": 12288,
    "simplepirQueryBytes": 65536,
    "doublepirQueryBytes": 57344,
    "simplepirRespBytes": 28672,
    "doublepirRespBytes": 12288,
    "serverTimeMs": 75,
    "clientQueryGenTimeMs": 748,
    "clientDecodeTimeMs": 0,
    "firstPassTimeMs": 15,
    "secondPassTimeMs": 4,
    "ringPackingTimeMs": 54,
    "sqrtNBytes": 8192,
    "allServerTimesMs": [
      76,
      75,
      75,
      75,
      75
    ],
    "stdDevServerTimeMs": 0.4
  }
}
➜  ypir git:(main) cargo run --release -- 10737418240
    Finished release [optimized] target(s) in 0.03s
     Running `target/release/run 10737418240`
Running YPIR (w/ DoublePIR) on a database of 10737418240 bits, and performing cross-client batching over 1 clients. 
The server performance measurement will be averaged over 5 trials.
Measurement completed. See the README for details on what the following fields mean.
Result:
{
  "offline": {
    "uploadBytes": 0,
    "downloadBytes": 0,
    "serverTimeMs": 61650,
    "clientTimeMs": 0,
    "simplepirPrepTimeMs": 54199,
    "simplepirHintBytes": 117440512,
    "doublepirHintBytes": 14680064
  },
  "online": {
    "uploadBytes": 997376,
    "downloadBytes": 12288,
    "simplepirQueryBytes": 262144,
    "doublepirQueryBytes": 229376,
    "simplepirRespBytes": 114688,
    "doublepirRespBytes": 12288,
    "serverTimeMs": 350,
    "clientQueryGenTimeMs": 2995,
    "clientDecodeTimeMs": 0,
    "firstPassTimeMs": 270,
    "secondPassTimeMs": 22,
    "ringPackingTimeMs": 54,
    "sqrtNBytes": 32768,
    "allServerTimesMs": [
      349,
      349,
      354,
      351,
      347
    ],
    "stdDevServerTimeMs": 2.3664319132398464
  }
}
```

