#include "poisson_problem/ddm/global_solve.hpp"
#include "poisson_problem/vtk/vtk_save.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

GlobalSolve::GlobalSolve(const Grid& grid, size_t overlap, double tolerance, const std::string& glue_strategy,
                         size_t maxit)
    : global_grid_(std::make_unique<Grid>(grid)), tolerance_(tolerance), glue_strategy_(glue_strategy), maxit_(maxit) {
    splitDomain(overlap);
    overlap_ = overlap;
}

void GlobalSolve::splitDomain(size_t overlap) {
    if (overlap == 0 || overlap > global_grid_->get_N_x() - 1)
        throw std::invalid_argument("invalid overlap");

    size_t N = global_grid_->get_N_x();
    size_t N1, N2;
    size_t temp = 0;
    double left_x0 = 0.0;
    double left_y0 = 0.0;
    double right_x0 = 0.0;
    double right_y0 = 0.0;
    size_t left_intersection = 0;
    size_t right_intersection = N - 1;

    if (overlap != N) {
        if (overlap % 2 == 0) {
            temp = overlap / 2 - 1;
        } else if (overlap % 2 != 0)
            temp = overlap / 2;
        if (N % 2 != 0) {
            N2 = N - N / 2 + 1 + temp;
            N1 = N - N2 + overlap + 1;
        } else if (N % 2 == 0) {
            N1 = N / 2 + 1 + temp;
            N2 = N - N1 + overlap + 1;
        }

        right_x0 = static_cast<double>(N1 - overlap - 1) * global_grid_->get_h();
        right_y0 = 0.0;

        left_intersection = N1 - 1 - overlap;
        right_intersection = overlap;
    }

    left_domain_solve_ = std::make_unique<LocalSolve>(
        0, Grid(left_x0, left_y0, N1, global_grid_->get_N_y(), global_grid_->get_h()), left_intersection);
    right_domain_solve_ = std::make_unique<LocalSolve>(
        1, Grid(right_x0, right_y0, N2, global_grid_->get_N_y(), global_grid_->get_h()), right_intersection);
}

void GlobalSolve::glueSolve(std::vector<double>& u, const std::string& strategy) {
    for (size_t j = 0; j < global_grid_->get_N_y(); j++) {
        for (size_t i = 0; i < global_grid_->get_N_x(); i++) {
            double chi_left;
            double chi_right;
            double x = global_grid_->getX(i);
            double y = global_grid_->getY(j);
            if (strategy == "chi_const") {
                chi_left = left_domain_solve_->chiConst(x, y);
                chi_right = right_domain_solve_->chiConst(x, y);
            } else if (strategy == "chi_continuous") {
                chi_left = left_domain_solve_->chiContinuous(x, y);
                chi_right = right_domain_solve_->chiContinuous(x, y);
            } else {
                throw std::invalid_argument("the strategy can be only chi_const or chi_continuous");
            }
            u[global_grid_->getK(i, j)] =
                left_domain_solve_->expand(x, y) * chi_left + right_domain_solve_->expand(x, y) * chi_right;
        }
    }
};

const std::vector<std::pair<size_t, double>>& GlobalSolve::getConvergenceHistory() const {
    return convergence_history_;
}

void GlobalSolve::clearConvergenceHistory() {
    convergence_history_.clear();
}
void GlobalSolve::changeOverlap(size_t overlap) {
    if (overlap == 0 || overlap > global_grid_->get_N_x())
        throw std::invalid_argument("invalid overlap");
    splitDomain(overlap);
    clearConvergenceHistory();
    overlap_ = overlap;
}
void GlobalSolve::changeGlueStrategy(const std::string& strategy) {
    if (strategy == "chi_const" || strategy == "chi_continuous") {
        glue_strategy_ = strategy;
        left_domain_solve_->resetSolve();
        right_domain_solve_->resetSolve();
        clearConvergenceHistory();
        return;
    }
    throw std::invalid_argument("the strategy can be only chi_const or chi_continuous");
}

LocalSolve GlobalSolve::getLeftLocalSolve() const {
    return *left_domain_solve_;
}

LocalSolve GlobalSolve::getRightLocalSolve() const {
    return *right_domain_solve_;
}