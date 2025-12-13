#include "poisson_problem/ddm/jacobi_solve.hpp"
#include "poisson_problem/norms/norms.hpp"
#include "poisson_problem/vector_operations/vector_operations.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>

JacobiSolve::JacobiSolve(const Grid& grid, size_t overlap, double tolerance, const std::string& glue_strategy,
                         size_t maxit)
    : GlobalSolve(grid, overlap, tolerance, glue_strategy, maxit) {}

void JacobiSolve::solve(std::vector<double>& u, size_t& iters) {
    std::fill(u.begin(), u.end(), 0.0);
    std::vector<double> u_prev = u;
    for (size_t iter = 1; iter <= maxit_; ++iter) {

        right_domain_solve_->give_boundary(*left_domain_solve_);
        left_domain_solve_->give_boundary(*right_domain_solve_);
        left_domain_solve_->solve();
        right_domain_solve_->solve();

        u.swap(u_prev);
        glueSolve(u, glue_strategy_);

        std::vector<double> residual = u - u_prev;
        double norm_residual = l2_norm(residual);

        convergence_history_.push_back(std::make_pair(iter, norm_residual));

        if (norm_residual < tolerance_) {
            iters = iter;
            return;
        }
    }
    throw std::runtime_error("Jacobi did not converge in " + std::to_string(maxit_) + " iterations");
}