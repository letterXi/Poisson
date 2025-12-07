#include "poisson_problem/ddm/seidel_solve.hpp"
#include "poisson_problem/norms/norms.hpp"

inline std::vector<double> operator-(const std::vector<double*>& u1, const std::vector<double*>& u2) {
    std::vector<double> res(u1.size());
    for (size_t i = 0; i < u1.size(); i++) {
        res[i] = *u1[i] - *u2[i];
    }
    return res;
}

SeidelSolve::SeidelSolve(const Grid& grid, size_t overlap) : GlobalSolve(grid, overlap) {}

void SeidelSolve::solve(std::vector<double>& u, size_t& iters) {
    u.reserve(3);
    for (size_t iter = 1; iter <= 1000; ++iter) {

        right_domain_solve_->give_boundary(*left_domain_solve_);
        left_domain_solve_->solve();
        left_domain_solve_->give_boundary(*right_domain_solve_);
        right_domain_solve_->solve();

        u[0] = 0;
        std::vector<double> residual =
            left_domain_solve_->getIntersectionSolve() - right_domain_solve_->getIntersectionSolve();
        if (inf_norm(residual) < 1e-10) {
            iters = iter;
            break;
        }
    }
}