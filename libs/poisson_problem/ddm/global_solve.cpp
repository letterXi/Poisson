#include "poisson_problem/ddm/global_solve.hpp"
#include "poisson_problem/vtk/vtk_save.hpp"
#include <cmath>
#include <iostream>

double norm(const std::vector<double*>& u1, const std::vector<double*>& u2) {
    double result = 0;
    for (size_t i = 0; i < u1.size(); i++) {
        result += (*u1[i] - *u2[i]) * (*u1[i] - *u2[i]);
    }
    return std::sqrt(result);
}
GlobalSolve::GlobalSolve(const Grid& grid, size_t overlap) : global_grid_(std::make_unique<Grid>(grid)) {
    splitDomain(overlap);
}

void GlobalSolve::splitDomain(size_t overlap) {
    size_t temp = 0;
    size_t N = global_grid_->get_N_x();
    size_t N1, N2;
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

    double left_x0 = 0.0;
    double left_y0 = 0.0;

    double right_x0 = static_cast<double>(N1 - overlap - 1) * global_grid_->get_h();
    double right_y0 = 0.0;

    size_t left_intersection = N1 - 1 - overlap;
    size_t right_intersection = overlap;

    left_domain_solve_ = std::make_unique<LocalSolve>(
        0, Grid(left_x0, left_y0, N1, global_grid_->get_N_y(), global_grid_->get_h()), left_intersection);
    right_domain_solve_ = std::make_unique<LocalSolve>(
        1, Grid(right_x0, right_y0, N2, global_grid_->get_N_y(), global_grid_->get_h()), right_intersection);
}

// void GlobalSolve::solve() {

// std::vector<double> u(global_grid_->points());
// vtkWriter sol("sol_ex", "sol", *global_grid_);
// for (int k = 1; k <= 1000; k++) {
//     std::vector<double> ex_u_2;
//     std::vector<double> ex_u_1;
//     for (size_t j = 0; j < global_grid_->get_N_y(); j++) {
//         for (size_t i = 0; i < global_grid_->get_N_x(); i++) {
//             ex_u_2.push_back(right_domain_solve_->chiContinuous(global_grid_->getX(i), global_grid_->getY(j)));
//         }
//     }
//     for (size_t j = 0; j < global_grid_->get_N_y(); j++) {
//         for (size_t i = 0; i < global_grid_->get_N_x(); i++) {
//             ex_u_1.push_back(left_domain_solve_->chiContinuous(global_grid_->getX(i), global_grid_->getY(j)));
//         }
//     }

// for (size_t j = 0; j < global_grid_->points(); j++) {
//     u[j] = (ex_u_1[j] + ex_u_2[j]);
// }

// sol.add_scalars(u, "glue_" + std::to_string(k));
// right_domain_solve_->give_boundary(*left_domain_solve_);
// left_domain_solve_->give_boundary(*right_domain_solve_);
// left_domain_solve_->solve();
// right_domain_solve_->solve();

// std::cout << k << ' '
//           << norm(left_domain_solve_->getIntersectionSolve(), right_domain_solve_->getIntersectionSolve())
//           << std::endl;

// if (norm(left_domain_solve_->getIntersectionSolve(), right_domain_solve_->getIntersectionSolve()) < 1e-10)
//     break;

// sol.write();
// }
// }