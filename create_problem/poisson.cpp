#include "poisson.hpp"
#include <cmath>
#include  <thread> 
#include <functional>


double exactFunction(double x, double y) {
    return 0.2 * (std::sin(2 * M_PI * x) - std::cos(12 * M_PI * y) + std::tan(x * y) + std::exp(-x * y));
}

double dirichletBoundaryFunction(double x, double y) {
    return exactFunction(x, y);
}

double sourceFunction(double x, double y) {
    return -0.2 * (-4 * M_PI * M_PI * std::sin(2 * M_PI * x) + 144 * M_PI * M_PI * std::cos(12 * M_PI * y) +
                   2 * (x * x + y * y) * std::pow(1.0 / std::cos(x * y), 2) * std::tan(x * y) +
                   (x * x + y * y) * std::exp(-x * y));
}

void Poisson(int n, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
             std::vector<double>& rhs, size_t thd) {
    size_t n2 = n * n;
    size_t inner_n = (n - 2)*(n - 2);
    size_t bound_n = n2 - inner_n;

    addr.clear();
    addr.reserve(n2 + 1);
    addr.push_back(0);
    size_t temp = 0;
    for(int j = 0; j < n; j++)
    {
	    for(int i = 0; i < n; i++)
	    {
		if (i == 0 || i == n - 1 || j == 0 || j == n - 1)
			temp++;
		else
			temp += 5;
		    
	    addr.push_back(temp);
	    }
    }


    cols.clear();
    cols.resize(bound_n + inner_n*5);

    vals.clear();
    vals.resize(bound_n + inner_n*5);

    rhs.resize(n2);

    double h = 1.0 / (double)(n - 1);
    size_t n_threads = thd;
    size_t pivot = n / n_threads;
    std::vector<std::thread> threads;
    for(size_t i = 0; i < n_threads; i++)
	{
	   size_t start_j = pivot * i;
	   size_t end_j = (i == n_threads - 1) ? n : pivot * (i + 1);
	   threads.push_back(std::thread(full_mat_vectors, n, h, std::ref(addr), std::ref(cols), std::ref(vals), std::ref(rhs), start_j, end_j));
	}
    for(size_t i = 0; i < n_threads; i++)
	   threads[i].join();
}

void full_mat_vectors(int n, double h, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double> &vals, std::vector<double> &rhs,  int start_j, int end_j)
{
    for (int j = start_j; j < end_j; ++j) {
        for (int i = 0; i < n; ++i) {
		int k = i + j * n;
	        size_t ind = addr[k];
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1) {
                cols[ind] = (k);
                vals[ind++] = (1);
                rhs[k] = dirichletBoundaryFunction(i * h, j * h);
            } else {
                cols[ind] = (k - n);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k - 1);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k);
                vals[ind++] = (4.0 / (h * h));

                cols[ind] = (k + 1);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k + n);
                vals[ind++] = (-1.0 / (h * h));

                rhs[k] = sourceFunction(i * h, j * h);
            }
        }
    }
}
