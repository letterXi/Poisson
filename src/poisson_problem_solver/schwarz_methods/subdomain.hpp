#ifndef SUBDOMAIN_HPP
#define SUBDOMAIN_HPP

#include "poisson_problem_solver/mat/csr_matrix.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include <functional>
#include <vector>

class Subdomain {
public:
    Subdomain(size_t N);
    void add_point(size_t point, double weight = 1);

    bool is_neighbor(size_t point) const;

    bool is_boundary(size_t point) const;

    bool contains(size_t point) const;

    const std::vector<size_t>& get_indices() const;

    const std::vector<double>& get_weights() const;

    const std::vector<int>& get_mask() const;

    const CsrMatrix& get_matrix() const;

    void create_matrix(double h, std::function<double(double, double)> source_f,
                       std::function<double(double, double)> boundary_f);

    void set_weight(size_t point, double weight);

    void set_swap_index(size_t index);

    size_t get_swap_index() const;

    const std::vector<double>& get_rhs() const;

    void get_overlap_boundary(const std::vector<double>& u_global);

    void solve();

    const std::vector<double>& get_u() const;

    void set_u(const std::vector<double>& u);

private:
    size_t N_;
    std::vector<size_t> indices_;
    std::vector<double> u_;
    std::vector<int> mask_;
    std::vector<double> weights_;
    CsrMatrix matrix_;
    std::vector<double> rhs_;
    size_t swap_index_;
    AmgclSolver solver_;
};

#endif
