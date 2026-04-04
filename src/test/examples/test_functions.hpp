#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include "poisson_problem_solver/grid/grid_regular2d.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

inline double exactFunction(double x, double y) {
    return 0.2 * (std::sin(2 * M_PI * x) - std::cos(12 * M_PI * y) + std::tan(x * y) + std::exp(-x * y));
}

inline void fill_u_exact(std::vector<double>& u, const RegularGrid2D& grid) {
    u.clear();
    const auto& points = grid.get_points();
    u.resize(points.size());
    for (size_t i = 0; i < u.size(); i++) {
        const auto& point = points[i];
        u[i] = exactFunction(point.x, point.y);
    }
}

inline double dirichletBoundaryFunction(double x, double y) {
    return exactFunction(x, y);
}

inline double sourceFunction(double x, double y) {
    return -0.2 * (-4 * M_PI * M_PI * std::sin(2 * M_PI * x) + 144 * M_PI * M_PI * std::cos(12 * M_PI * y) +
                   2 * (x * x + y * y) * (1.0 / std::cos(x * y)) * (1.0 / std::cos(x * y)) * std::tan(x * y) +
                   (x * x + y * y) * std::exp(-x * y));
}

inline void print_centered(const std::string& word, int width) {
    int side = (width - word.length()) / 2;
    if (side <= 0) {
        std::cout << word << std::endl;
        return;
    }
    std::cout << std::setfill('=') << std::setw(side + word.length()) << word << std::setw(side) << "" << std::endl;
    std::cout << std::setfill(' ');
}

inline void solve_with_screenshot(const std::string& stem, std::vector<double>& u, const std::vector<double>& u_exact,
                                  const RegularGrid2D& grid, SchwarzSolver& solver) {
    print_centered(stem, 40);
    double error = std::numeric_limits<double>::max();
    size_t iters;
    VtkWriter::StepManager schwarz_filename_mgr(stem, 1);
    bool force_write = false;
    size_t N = solver.N();

    for (size_t iter = 0; iter < 10000; iter++) {
        std::string path = schwarz_filename_mgr.add(iter, force_write);
        if (!path.empty()) {
            std::ofstream fs(path);
            VtkWriter::append_header(stem, fs);
            VtkWriter::append_points(grid.get_points(), N, fs);
            VtkWriter::append_point_data_header(N * N, fs);
            VtkWriter::add_point_data(u, "U", fs);
            VtkWriter::add_point_data(solver.overlap_point(), "overlap", fs);
            fs.close();
        }
        if (force_write)
            break;

        solver.iterate(u);

        solver.connect_solves(u);
        double current_error = norm_inf_2(u, u_exact);
        if (std::abs(error - current_error) < 1e-11) {
            force_write = true;
        }
        error = current_error;
        iters = iter;
    }
    std::cout << "overlap_ratio = " << solver.overlap_ratio() << "%" << std::endl;
    std::cout << "iters = " << iters << ", error = " << error << std::endl << std::endl;
}

#endif
