#include "poisson_problem_solver/utils/norms.hpp"
#include <stdexcept>
#include "schwarz_solver.hpp"
#include <algorithm>
#include <unordered_set>

SchwarzSolver::SchwarzSolver(size_t N, std::vector<size_t> mask, std::function<double(double, double)> source_function,
                             std::function<double(double, double)> boundary_function, size_t maxiter, double tolerance)
    : N_(N), mask_(mask), source_function_(source_function), boundary_function_(boundary_function), maxiter_(maxiter), tolerance_(tolerance) {
    std::vector<size_t> unique_idx = this->find_unique_indices();
    for (size_t i = 0; i < unique_idx.size(); i++) {
        subdomains_.push_back(std::make_unique<Subdomain>(N));
        id_to_idx_[unique_idx[i]] = i;
    }
    for (size_t point = 0; point < mask_.size(); point++) {
        size_t id = mask_[point];
        subdomains_[id_to_idx_[id]]->add_point(point);
    }
}
const std::vector<std::unique_ptr<Subdomain>>& SchwarzSolver::get_subdomains() const {
    return subdomains_;
}
const std::unordered_map<size_t, size_t>& SchwarzSolver::get_id_to_idx() const {
    return id_to_idx_;
}

std::vector<size_t> SchwarzSolver::find_unique_indices() const {
    std::unordered_set<size_t> all_ids;
    std::vector<size_t> unique_indices;
    for (auto point : mask_)
        if (all_ids.insert(point).second)
            unique_indices.push_back(point);
    return unique_indices;
}

void SchwarzSolver::set_overlap(size_t overlap) {
    for (size_t o = 1; o <= overlap; o++) {
        for (auto& sub : subdomains_) {
            if (o == overlap)
                sub->set_swap_index(sub->get_indices().size());
            std::vector<size_t> temp;
            for (size_t point = 0; point < N_ * N_; point++)
                if (sub->is_neighbor(point))
                    temp.push_back(point);
            for (auto t : temp)
                sub->add_point(t);
        }
    }
    for (size_t point = 0; point < N_ * N_; point++) {
        size_t over = 0;
        for (const auto& sub : subdomains_)
            if (sub->contains(point))
                over++;
        double weight = 1.0 / static_cast<double>(over);
        for (auto& sub : subdomains_)
            if (sub->contains(point))
                sub->set_weight(point, weight);
    }
    this->create_slaes();
}

void SchwarzSolver::create_slaes() {
    double h = 1.0 / static_cast<double>(N_);
    for (auto& sub : subdomains_)
        sub->create_matrix(h, source_function_, boundary_function_);
}

void SchwarzSolver::connect_solves(std::vector<double>& u) {
    std::fill(u.begin(), u.end(), 0.0);

    for (const auto& sub : subdomains_) {
        const auto& indices = sub->get_indices();
        const auto& weights = sub->get_weights();
        const auto& local_u = sub->get_u();

        for (size_t i = 0; i < indices.size(); i++)
            u[indices[i]] += weights[i] * local_u[i];
    }
}

void SchwarzSolver::solve(std::vector<double>& u, size_t& iters) {
    for (size_t iter = 1; iter <= maxiter_; iter++) {
        std::vector<double> u_old = u;
        this->iterate(u);
        if (norm_inf(diff_of(u, u_old)) < tolerance_) {
            iters = iter;
            break;
        }
    }
    if (iters >= maxiter_)
      throw std::runtime_error("Convergence failed: residuals too high");
}

void SchwarzSolver::initialize(const std::vector<double>& global_u) {
    for (auto& sub : subdomains_) {
        const auto& indices = sub->get_indices();
        std::vector<double> local_u(indices.size());

        for (size_t i = 0; i < indices.size(); i++)
            local_u[i] = global_u[indices[i]];
        sub->set_u(local_u);
    }
}
