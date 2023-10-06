#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>


struct Timer {
    std::string name;

    std::chrono::time_point<std::chrono::steady_clock> startTime;
    double lastElapsed = 0.;
    size_t numElapsed = 0;
    double sumElapsed = 0;
    double minElapsed, maxElapsed;

    inline explicit Timer(const std::string &name_) : name(name_) {}

    inline void start() {
        startTime = std::chrono::steady_clock::now();
    }

    inline void stop() {
        auto stopTime = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = stopTime - startTime;
        lastElapsed = elapsed.count();

        sumElapsed += lastElapsed;
        if (0 == numElapsed) {
            minElapsed = lastElapsed;
            maxElapsed = lastElapsed;
        } else {
            minElapsed = std::min(minElapsed, lastElapsed);
            maxElapsed = std::max(maxElapsed, lastElapsed);
        }

        ++numElapsed;
    }

    inline double meanElapsed() {
        return 0 == numElapsed ? 0. : sumElapsed / numElapsed;
    }

    inline void print(int mpiRank) {
        std::cout << "Timer " << std::setw(16) << std::left << name << " on rank " << std::right << std::setw(3) << mpiRank << " : "
                  << std::setw(6) << (size_t) (1e3 * sumElapsed) << " ms, "
                  << std::setw(6) << (size_t) (1e3 * meanElapsed()) << " ms, "
                  << std::setw(6) << (size_t) (1e3 * minElapsed) << " ms, "
                  << std::setw(6) << (size_t) (1e3 * maxElapsed) << " ms, "
                  << std::setw(3) << numElapsed
                  << " (total, mean, min, max, timings)"
                  << std::endl;
    }
};