#include "poisson.hpp"
#include <vector>
#include <chrono>
#include <iostream>

int main()
{
	size_t runs = 100;
	int n = 250;
	for(int thd = 1; thd <= 15; thd++)
	{
	std::vector<double> vals;
	std::vector<double> rhs;
	std::vector<size_t> addr;
	std::vector<size_t> cols;


	long total_time = 0;
	std::cout << "Matrix filling (" << thd << " threads)" << std::endl;
	for(size_t i = 0; i < runs; i++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		Poisson(n, addr, cols, vals, rhs, thd);
		auto end = std::chrono::high_resolution_clock::now();
		total_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	}
	double average_time = (double) total_time / (double) runs;
	std::cout << "Average time (" << thd << " threads): " << average_time << std::endl << std::endl;
	}

	return 0;
}


