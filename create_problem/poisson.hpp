#ifndef POISSON_HPP
#define POISSON_HPP
#include <vector>

double exactFunction(double x, double y); 

double dirichletBoundaryFunction(double x, double y);

double sourceFunction(double x, double y);

void Poisson(int n, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
             std::vector<double>& rhs);
void full_mat_vectors(int n, double h, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double> &vals, std::vector<double> &rhs,  int start_j, int end_j);

#endif
