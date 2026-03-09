#include "jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include <algorithm>
#include <execution>

void JacobiSchwarzSolver::iterate(std::vector<double>& u) {
    std::for_each(std::execution::seq, subdomains_.begin(), subdomains_.end(), [&](const auto& sub) {
        sub->get_overlap_boundary(u);
        sub->solve();
    });
    this->connect_solves(u);
}

void JacobiSchwarzSolver::parallel_iterate(std::vector<double>& u) {
    std::for_each(std::execution::par, subdomains_.begin(), subdomains_.end(), [&](const auto& sub) {
        sub->get_overlap_boundary(u);
        sub->solve();
    });
    this->connect_solves(u);
}

void JacobiSchwarzSolver::parallel_solve(std::vector<double>& u, size_t& iters) {
    for (size_t iter = 1; iter <= maxiter_; iter++) {
        std::vector<double> u_old = u;
        this->parallel_iterate(u);
        if (norm_inf(diff_of(u, u_old)) < tolerance_) {
            iters = iter;
            break;
        }
    }
    if (iters >= maxiter_)
      throw std::runtime_error("Convergence failed: residuals too high");
}


std::string JacobiSchwarzSolver::get_name() const {
    return "JacobiSchwarz";
}
