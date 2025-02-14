// Well, thanks ChatGPT for writing this clean logger for me.
#ifndef LOGGING_H
#define LOGGING_H

#include <unordered_map>
#include <vector>
#include <chrono>
#include <string>


// print for debug. Easily turn on/off by defining _DEBUG
#ifdef _DEBUG
#define DEBUG_PRINT(s) std::cout << s << std::endl;
#endif

#ifdef _BENCHMARK
#define DEBUG_PRINT(s) ; // do nothing
#endif

#define BENCH_PRINT(s) std::cout << s << std::endl;

#define PRINT_BAR BENCH_PRINT("==============================================================");

// predefine some name for logging
#define CORE_TIME "Core"
#define FST_DIM_TIME "First dim"
#define OTHER_DIM_TIME "Other dim"
#define EXPAND_TIME "Expand"
#define CONVERT_TIME "Convert"
#define SERVER_TOT_TIME "Server total"
#define CLIENT_TOT_TIME "Client total"
#define EXTERN_PROD_TOT_TIME "External product total"
#define EXTERN_PROD_MAT_MULT_TIME "External product mat mult"
#define FST_NTT_TIME "First dim NTT"
#define OTHER_NTT_TIME "Other dim NTT"
#define DECOMP_RLWE_TIME "Decomp RLWE"
#define FST_DELEY_MOD_TIME "First dim delay mod"


// Hierarchical structure for pretty result
const std::unordered_map<std::string, std::vector<std::string>> LOG_HIERARCHY = {
    {SERVER_TOT_TIME, {EXPAND_TIME, CONVERT_TIME, FST_DIM_TIME, OTHER_DIM_TIME}},
    {FST_DIM_TIME, {CORE_TIME, FST_DELEY_MOD_TIME, FST_NTT_TIME}},
    {OTHER_DIM_TIME, {OTHER_NTT_TIME, EXTERN_PROD_TOT_TIME}},
    {EXTERN_PROD_TOT_TIME, {DECOMP_RLWE_TIME, EXTERN_PROD_MAT_MULT_TIME}}
};



class TimerLogger {
private:
    // Stores start times of active sections
    std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> startTimes;

    // Stores all timing results for multiple experiments
    std::vector<std::unordered_map<std::string, double>> experimentRecords;

    // Stores timing data for the current experiment
    std::unordered_map<std::string, double> currentExperiment;

    // Private constructor for Singleton
    TimerLogger() = default;

    // Recursive helper for pretty printing
    void prettyPrintHelper(const std::string& section, const std::string& prefix, bool isLast) const;

public:
    // Singleton instance
    static TimerLogger& getInstance();

    // Start logging time for a section
    void start(const std::string& sectionName);

    // Stop logging time for a section
    void end(const std::string& sectionName);

    // End the current experiment and start a new one
    void endExperiment();

    // Print results for specific experiment. -1 to print all experiments
    void printResults(size_t expId = -1);

    // Compute and print average time across experiments
    void printAverageResults();

    double getAvgTime(const std::string& sectionName);

    // Pretty print hierarchical results
    void prettyPrint();

    // Prevent copying
    TimerLogger(const TimerLogger&) = delete;
    TimerLogger& operator=(const TimerLogger&) = delete;
};


// Macros for easy time logging
#define TIME_START(sec) TimerLogger::getInstance().start(sec)
#define TIME_END(sec) TimerLogger::getInstance().end(sec)
#define END_EXPERIMENT() TimerLogger::getInstance().endExperiment()
#define PRINT_RESULTS(expId) TimerLogger::getInstance().printResults(expId)
#define PRINT_AVERAGE_RESULTS() TimerLogger::getInstance().printAverageResults()
#define GET_AVG_TIME(sec) TimerLogger::getInstance().getAvgTime(sec)
#define PRETTY_PRINT() TimerLogger::getInstance().prettyPrint()


#endif // LOGGER_H
