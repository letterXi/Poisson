#ifndef JACOBI_SOLVE_HPP
#define JACOBI_SOLVE_HPP

#include "poisson_problem/ddm/global_solve.hpp"
#include <optional>

class JacobiSolve : public GlobalSolve {
public:
    JacobiSolve(const Grid& gird, size_t overlap, double tolerance, const std::string& glue_strategy,
                size_t maxit = 1000);
    void solve(std::vector<double>& u, size_t& iters) override;
};

#endif