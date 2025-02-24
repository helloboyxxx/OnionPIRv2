#include "logging.h"
#include <iostream>


// Singleton instance getter
TimerLogger& TimerLogger::getInstance() {
    static TimerLogger instance;
    return instance;
}

// Start logging time for a section
void TimerLogger::start(const std::string& sectionName) {
    startTimes[sectionName] = std::chrono::high_resolution_clock::now();
}

// Stop logging time for a section
void TimerLogger::end(const std::string& sectionName) {
    auto it = startTimes.find(sectionName);
    if (it != startTimes.end()) {
        double duration = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - it->second).count();
        currentExperiment[sectionName] += duration;
        startTimes.erase(it);
    }
}

// End the current experiment and store the results
void TimerLogger::endExperiment() {
    experimentRecords.push_back(currentExperiment);
    currentExperiment.clear();  // Reset for the next experiment
}

// Print timing results for a specific experiment or all experiments
void TimerLogger::printResults(size_t expId) {
    if (experimentRecords.empty()) {
        std::cout << "No experiments recorded.\n";
        return;
    }

    if (expId == -1) {
        // Print all experiments
        std::cout << "========================== Experiment Timing Results =========================\n";
        for (size_t i = 0; i < experimentRecords.size(); ++i) {
            std::cout << "Experiment " << i + 1 << ":\n";
            for (const auto& [section, time] : experimentRecords[i]) {
                std::cout << "  " << section << ": " << time << " ms\n";
            }
        }
    } else if (expId >= 1 && expId <= experimentRecords.size()) {
        // Print a specific experiment
        std::cout << "========================== Experiment " << expId << " Timing Results =========================\n";
        for (const auto& [section, time] : experimentRecords[expId - 1]) {
            std::cout << section << ": " << time << " ms\n";
        }
    } else {
        std::cout << "Invalid experiment index. Please choose between 1 and " << experimentRecords.size() << ".\n";
    }
}

// Compute and print average timing results across experiments
void TimerLogger::printAverageResults() {
    if (experimentRecords.empty()) {
        std::cout << "No experiments recorded.\n";
        return;
    }

    std::unordered_map<std::string, double> totalTimes;
    std::unordered_map<std::string, int> count;

    // Sum up all times for each section
    for (const auto& record : experimentRecords) {
        for (const auto& [section, time] : record) {
            totalTimes[section] += time;
            count[section]++;
        }
    }

    std::cout << "===== Average Timing Results =====\n";
    for (const auto& [section, totalTime] : totalTimes) {
        double avgTime = totalTime / count[section];
        std::cout << section << ": " << (size_t)avgTime << " ms\n";
    }
}


double TimerLogger::getAvgTime(const std::string& sectionName) {
    double sum = 0;
    for (const auto& record : experimentRecords) {
        sum += record.at(sectionName);
    }
    return sum / experimentRecords.size();
}



// Recursive function to print hierarchical sections with indentation
void TimerLogger::prettyPrintHelper(const std::string& section, const std::string& prefix, bool isLast) const {
    auto it = experimentRecords.back().find(section);
    if (it != experimentRecords.back().end()) {
        std::cout << prefix;
        std::cout << (isLast ? "└── " : "├── ");
        std::cout << section << ": " << (size_t)it->second << " ms\n";
    }

    auto subSections = LOG_HIERARCHY.find(section);
    if (subSections != LOG_HIERARCHY.end()) {
        size_t numSubSections = subSections->second.size();
        for (size_t i = 0; i < numSubSections; i++) {
            std::string newPrefix = prefix + (isLast ? "    " : "│   ");
            prettyPrintHelper(subSections->second[i], newPrefix, i == numSubSections - 1);
        }
    }
}

// Pretty print hierarchical results in tree format
void TimerLogger::prettyPrint() {
    if (experimentRecords.empty()) {
        std::cout << "No experiments recorded.\n";
        return;
    }

    std::cout << "========== Average Results ==========\n";
    prettyPrintHelper(SERVER_TOT_TIME, "", false);
    prettyPrintHelper(CLIENT_TOT_TIME, "", true);
}

void TimerLogger::cleanup() {
    startTimes.clear();
    experimentRecords.clear();
    currentExperiment.clear();
}

