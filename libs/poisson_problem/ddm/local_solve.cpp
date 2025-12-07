#include "poisson_problem/ddm/local_solve.hpp"
#include "poisson_problem/mat_solver/csr_mat_solver.hpp"
LocalSolve::LocalSolve(size_t id, Grid grid, size_t intersection)
    : id_(id), grid_(grid), intersection_(intersection), u_(grid_.points(), 0.0), boundary_(grid_.get_N_y(), 0.0),
      slae_(grid_) {}
void LocalSolve::give_boundary(LocalSolve& other) const {
    other.set_boundary(slice(intersection_));
}

std::vector<double> LocalSolve::slice(size_t i) const {
    std::vector<double> u_slice;
    for (size_t j = 0; j < grid_.get_N_y(); j++)
        u_slice.push_back(u_[grid_.getK(i, j)]);
    return u_slice;
}
void LocalSolve::set_boundary(const std::vector<double>& boundary) {
    boundary_ = boundary;
}

void LocalSolve::solve() {
    AmgclSolver solver_({{"solver.type", "gmres"}});
    solver_.set_matrix(slae_.get_mat());
    std::vector<double> rhs = slae_.get_rhs(grid_, sourceFunction, [this](double x, double y) {
        if ((grid_.getI(x) == 0 && id_ == 1) || (grid_.getI(x) == grid_.get_N_x() - 1 && id_ == 0)) {
            if ((grid_.getJ(y) > 0 && grid_.getJ(y) < grid_.get_N_y() - 1)) {
                return boundary_[grid_.getJ(y)];
            }
        }
        return dirichletBoundaryFunction(x, y);
    });
    solver_.solve(rhs, u_);
}
const std::vector<double>& LocalSolve::get_solve() const {
    return u_;
}

double LocalSolve::expand(double x, double y) const {
    if (grid_.is_in_domain(x, y))
        return u_[grid_.getK(grid_.getI(x), grid_.getJ(y))];
    return 0.0;
}

double LocalSolve::chiConst(double x, double y) const {
    size_t i = grid_.getI(x);
    if (!grid_.is_in_domain(x, y))
        return 0.0;
    if (id_ == 0) {
        if (i >= intersection_)
            return 0.5;
        return 1.0;
    } else {
        if (i <= intersection_)
            return 0.5;
        return 1.0;
    }
}

double LocalSolve::chiContinuous(double x, double y) const {
    size_t i = grid_.getI(x);
    if (!grid_.is_in_domain(x, y))
        return 0.0;
    if (id_ == 0) {
        if (i >= intersection_)
            return (grid_.getX(grid_.get_N_x() - 1) - x) /
                   (grid_.get_h() * static_cast<double>(grid_.get_N_x() - 1 - intersection_));
        return 1.0;
    } else {
        if (i <= intersection_)
            return (x - grid_.getX(0)) / (grid_.get_h() * static_cast<double>(intersection_));
        return 1.0;
    }
}

const Grid& LocalSolve::get_grid() const {
    return grid_;
}
