#include <thread>
#include <iostream>
#include <vector>
#include <functional>

void fill_vector(std::vector<double>&  arr, size_t i_begin, size_t i_end)
{
	for(size_t i = i_begin; i < i_end; i++)
	{
		arr[i] = double(i);
	}
}
int main()
{
	size_t n_threads = 4;
	size_t n = 10;
	std::vector<double> arr(n, 0.0);
	size_t pivot = n / n_threads;
	std::vector<std::thread> threads;
	for(size_t i = 0; i < n_threads; i++)
	{
		size_t start = i * pivot;
		size_t finish = (i == n_threads - 1) ? n : (i + 1) * pivot;
		threads.push_back(std::thread(fill_vector, std::ref(arr), start, finish));
	}
	for(size_t i = 0; i < n_threads; i++)
	{
		threads[i].join();
	}
	for(size_t i = 0; i < arr.size(); i++)
	{
		std::cout << arr[i] << ' ';
	}
	return 0;
}
