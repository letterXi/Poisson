#ifndef AMGCL_SOLVER_HPP
#define AMGCL_SOLVER_HPP

#include "CsrMatrix.hpp"
#include <map>
#include <memory>
#include <string>

class AmgclSolver {
public:
    AmgclSolver(size_t maxit = 1000, double eps = 1e-8);
    ~AmgclSolver();
    AmgclSolver(std::initializer_list<std::pair<std::string, std::string>> amgc_params);
    void solve(const std::vector<double>& rhs, std::vector<double>& ret) const;
    void set_matrix(const CsrMatrix& mat);

private:
    size_t maxit_;
    double tolerance_;
    std::map<std::string, std::string> params_;

    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif