#include "original_schwarz_solver.hpp"
#include <algorithm>

void OriginalSchwarzSolver::iterate(std::vector<double>& u) {
    for (auto& sub : subdomains_) {
        sub->get_overlap_boundary(u);
        sub->solve();

        const auto& indices = sub->get_indices();
        const auto& local_u = sub->get_u();
        for (size_t i = 0; i < indices.size(); i++)
            u[indices[i]] = local_u[i];
    }
}

std::string OriginalSchwarzSolver::get_name() const {
    return "OriginalSchwarz";
}
