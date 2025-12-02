#ifndef POISSON_HPP
#define POISSON_HPP
#include "shwartz.hpp"
#include <vector>
#include <functional>

double exactFunction(double x, double y);

double dirichletBoundaryFunction(double x, double y);

double sourceFunction(double x, double y);

void Poisson(const Mesh& mesh, std::function<double(double, double)> source_function,
             std::function<double(double, double)> boundary_function, std::vector<size_t>& addr,
             std::vector<size_t>& cols, std::vector<double>& vals, std::vector<double>& rhs); 
// void Poisson(int n, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
//             std::vector<double>& rhs, size_t threads = 0);
             // void full_mat_vectors(int n, double h, std::vector<size_t>& addr, std::vector<size_t>& cols,
             // std::vector<double>& vals,
//                      std::vector<double>& rhs, int j);

#endif
