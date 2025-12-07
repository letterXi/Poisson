#ifndef POISSON_HPP
#define POISSON_HPP
#include "poisson_problem/grid/grid.hpp"
#include "poisson_problem/mat/csr_matrix.hpp"
#include <functional>
#include <vector>

class PoissonSlaeBuilder {
private:
    CsrMatrix mat_;
    std::vector<double> rhs_;

public:
    PoissonSlaeBuilder(const Grid& grid);
    const CsrMatrix& get_mat() const;
    std::vector<double>  get_rhs(const Grid &grid_, std::function<double(double, double)> source_function,
                 std::function<double(double, double)> boundary_function);
};
double exactFunction(double x, double y);

double dirichletBoundaryFunction(double x, double y);

double sourceFunction(double x, double y);

void Poisson(const Grid& grid, std::function<double(double, double)> source_function,
             std::function<double(size_t, size_t)> boundary_function, std::vector<size_t>& addr,
             std::vector<size_t>& cols, std::vector<double>& vals, std::vector<double>& rhs);
// void Poisson(int n, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
//             std::vector<double>& rhs, size_t threads = 0);
// void full_mat_vectors(int n, double h, std::vector<size_t>& addr, std::vector<size_t>& cols,
// std::vector<double>& vals,
// std::vector<double>& rhs, int j);

#endif
