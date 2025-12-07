#ifndef SEIDEL_SOLVE_HPP
#define SEIDEL_SOLVE_HPP
#include "poisson_problem/ddm/global_solve.hpp"

class SeidelSolve : public GlobalSolve {
public:
    SeidelSolve(const Grid& gird, size_t overlap, double tolerance, const std::string& glue_strategy,
                size_t maxit = 1000);
    void solve(std::vector<double>& u, size_t& iters) override;
};

#endif
