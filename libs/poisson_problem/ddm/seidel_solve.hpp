#ifndef SEIDEL_SOLVE_HPP
#define SEIDEL_SOLVE_HPP
#include "poisson_problem/ddm/global_solve.hpp"

class SeidelSolve : public GlobalSolve {
public:
    SeidelSolve(const Grid& gird, size_t overlap);
    void solve(std::vector<double>& u, size_t& iters) override;
};

#endif
