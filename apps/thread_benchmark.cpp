#include "Poisson.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

int main() {

    std::ofstream file("thread_benchmark_result.txt");
    if (!file.is_open())
        return 1;

    size_t runs = 10000;
    int n = 100;
    for (int thd = 0; thd <= 15; thd++) {
        std::cout << "thd = " << thd << std::endl;
        long total_time = 0;
        if (thd > 1) {
            file << "Matrix filling (" << thd << " threads)" << std::endl;
        } else if (thd == 0) {
            file << "Matrix is filling in main thread" << std::endl;
        } else if (thd == 1) {
            file << "Matrix is filling using std::execution::par" << std::endl;
        }
        for (size_t i = 0; i < runs; i++) {
            std::vector<double> vals;
            std::vector<double> rhs;
            std::vector<size_t> addr;
            std::vector<size_t> cols;
            auto start = std::chrono::high_resolution_clock::now();
            Poisson(n, addr, cols, vals, rhs, thd);
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }
        double average_time = (double)total_time / (double)runs;
        file << "Average time (" << thd << " threads): " << average_time << std::endl << std::endl;
    }

    return 0;
}
