#ifndef JACOBI_SOLVE_HPP
#define JACOBI_SOLVE_HPP

#include "poisson_problem/ddm/global_solve.hpp"

class JacobiSolve : public GlobalSolve {
public:
    JacobiSolve(const Grid& gird, size_t overlap);
    void solve(std::vector<double>& u, size_t& iters) override;
};

#endif