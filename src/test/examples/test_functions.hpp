#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include "poisson_problem_solver/grid/grid_regular2d.hpp"
#include "poisson_problem_solver/schwarz_methods/jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/original_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/schwarz_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include <chrono>
#include <cmath>
#include <filesystem>
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

inline void save_overlap_vtk(const std::string& header, const RegularGrid2D& grid, SchwarzSolver& solver) {
    std::ofstream fs(header + ".vtk");
    VtkWriter::append_header(header, fs);
    VtkWriter::append_points(grid.get_points(), solver.N(), fs);
    VtkWriter::append_point_data_header(solver.N() * solver.N(), fs);
    const auto& subs = solver.get_subdomains();
    VtkWriter::add_point_data(solver.overlap_point(), "overlap", fs);
    for (size_t i = 0; i < subs.size(); i++)

    fs.close();
}
inline void solve_with_screenshot(const std::string& stem, std::vector<double>& u, const std::vector<double>& u_exact,
                                  const RegularGrid2D& grid, SchwarzSolver& solver) {
    print_centered(stem, 40);
    double error = norm_inf_2(u, u_exact);
    size_t iters;
    VtkWriter::StepManager schwarz_filename_mgr(stem, 1);
    bool force_write = false;
    size_t N = solver.N();

    for (size_t iter = 1; iter < 10000; iter++) {
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
        if (std::abs(error - current_error) < 1e-15) {
            force_write = true;
        }
        error = current_error;
        iters = iter;
    }
    std::cout << "overlap_ratio = " << solver.overlap_ratio() << "%" << std::endl;
    std::cout << "iters = " << iters << ", error = " << error << std::endl << std::endl;
}

inline void log_convergence_metrics(size_t N, const std::vector<size_t> mask, const std::string& path) {
    std::vector<double> u(N * N);
    std::vector<double> u_exact(N*N);
    RegularGrid2D grid(0,0,1,1, N,N);
    fill_u_exact(u_exact, grid);

    size_t overlap, iters;

    if (std::filesystem::is_directory(path))
        std::filesystem::remove_all(path);
    std::filesystem::create_directory(path);

    overlap = 1;
    std::ofstream f1(path + "/jacobi.csv");
    f1 << "overlap,iters,overhead,v,time" << std::endl;
    while (true) {
        JacobiSchwarzSolver jac_sol(N, mask, sourceFunction, dirichletBoundaryFunction, 100000);
        jac_sol.set_overlap(overlap);
        std::fill(u.begin(), u.end(), 0);
        jac_sol.initialize(u);
        
        auto start = std::chrono::steady_clock::now();
        jac_sol.solve(u, iters);
        auto end = std::chrono::steady_clock::now();
        double elapsed_ms = std::chrono::duration<double,std::milli>(end - start).count();

        f1 << overlap << ',' << iters << ',' << jac_sol.overhead_ratio() << ',' << iters * jac_sol.overhead_ratio()
           << ',' << elapsed_ms << std::endl;

        overlap++;
        if (jac_sol.is_collapse())
            break;
    }
    f1.close();

    overlap = 1;
    std::ofstream f2(path + "/schwarz.csv");
    f2 << "overlap,iters,overhead,v,time" << std::endl;
    while (true) {
        OriginalSchwarzSolver sch_sol(N, mask, sourceFunction, dirichletBoundaryFunction, 10000);
        sch_sol.set_overlap(overlap);
        std::fill(u.begin(), u.end(), 0);
        sch_sol.initialize(u);

        auto start = std::chrono::steady_clock::now();
        sch_sol.solve(u, iters);
        auto end = std::chrono::steady_clock::now();
        double elapsed_ms = std::chrono::duration<double,std::milli>(end - start).count();

        f2 << overlap << ',' << iters << ',' << sch_sol.overhead_ratio() << ',' << iters * sch_sol.overhead_ratio()
           << ',' << elapsed_ms << std::endl;
        overlap++;
        if (sch_sol.is_collapse())
            break;
    }
    f2.close();
}

#endif
