#include "poisson_problem/ddm/jacobi_solve.hpp"
#include "poisson_problem/ddm/seidel_solve.hpp"
#include "poisson_problem/grid/grid.hpp"
#include "poisson_problem/vtk/vtk_save.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

void iter_residual_csv_writer(const std::string& file_name, const std::vector<std::pair<size_t, double>>& iter_res) {
    std::filesystem::create_directories("convergence_methods");
    std::ofstream file("convergence_methods/" + file_name + ".csv");
    file << "iter," + file_name << std::endl;
    for (size_t i = 0; i < iter_res.size(); i++)
        file << iter_res[i].first << ',' << iter_res[i].second << std::endl;
    file.close();
}

void overlap_iter_csv_writer(const std::string& file_name, const std::vector<size_t>& overlaps,
                             const std::vector<size_t>& iters) {
    std::filesystem::create_directories("compare_overlap_methods");
    std::ofstream file("compare_overlap_methods/" + file_name + ".csv");
    file << "overlap," + file_name << std::endl;
    for (size_t i = 0; i < overlaps.size(); i++) {
        file << overlaps[i] << ',' << iters[i] << std::endl;
    }
    file.close();
}

int main() {

    double x_0 = 0.0;
    double y_0 = 0.0;
    size_t N = 100;
    double h = 1.0 / static_cast<double>(N - 1);

    Grid grid(x_0, y_0, N, N, h);

    size_t overlap = 20;
    double tolerance = 1e-6;
    JacobiSolve jacobi_solver(grid, overlap, tolerance, "chi_const");

    vtkWriter jacobi_writer("jacobi_solve", "jacobi_solve", grid);

    SeidelSolve seidel_solver(grid, overlap, tolerance, "chi_const");

    vtkWriter seidel_writer("seidel_solve", "seidel_solve", grid);

    std::vector<double> u(grid.points());
    size_t iters;

    jacobi_solver.solve(u, iters);
    jacobi_writer.add_scalars(u, "U_chi_const");
    seidel_solver.solve(u, iters);
    seidel_writer.add_scalars(u, "U_chi_const");

    iter_residual_csv_writer("jacobi_chi_const_iter", jacobi_solver.getConvergenceHistory());
    iter_residual_csv_writer("seidel_chi_const_iter", seidel_solver.getConvergenceHistory());

    jacobi_solver.changeGlueStrategy("chi_continuous");
    seidel_solver.changeGlueStrategy("chi_continuous");

    seidel_solver.solve(u, iters);
    seidel_writer.add_scalars(u, "U_chi_continuous");

    jacobi_solver.solve(u, iters);
    jacobi_writer.add_scalars(u, "U_chi_continuous");

    jacobi_writer.write();
    seidel_writer.write();

    iter_residual_csv_writer("jacobi_chi_continuous_iter", jacobi_solver.getConvergenceHistory());
    iter_residual_csv_writer("seidel_chi_continuous_iter", seidel_solver.getConvergenceHistory());

    std::vector<size_t> overlaps(N);
    std::vector<size_t> jacobi_iters;
    std::vector<size_t> seidel_iters;
    size_t counter = 1;
    std::generate(overlaps.begin(), overlaps.end(), [&counter]() { return counter++; });
    std::cout << "overlaps[-1] = " << overlaps.back() << std::endl;

    jacobi_solver.changeGlueStrategy("chi_const");
    seidel_solver.changeGlueStrategy("chi_const");
    for (auto overlap : overlaps) {
        jacobi_solver.changeOverlap(overlap);
        seidel_solver.changeOverlap(overlap);
        jacobi_solver.solve(u, iters);
        jacobi_iters.push_back(iters);
        std::cout << iters << std::endl;
        seidel_solver.solve(u, iters);
        seidel_iters.push_back(iters);
        std::cout << iters << std::endl;
        std::cout << "chi_const overlap = " << overlap << std::endl;
    }
    overlap_iter_csv_writer("jacobi_chi_const_overlap", overlaps, jacobi_iters);
    overlap_iter_csv_writer("seidel_chi_const_overlap", overlaps, seidel_iters);
    jacobi_iters.clear();
    seidel_iters.clear();

    jacobi_solver.changeGlueStrategy("chi_continuous");
    seidel_solver.changeGlueStrategy("chi_continuous");
    for (auto overlap : overlaps) {
        jacobi_solver.changeOverlap(overlap);
        seidel_solver.changeOverlap(overlap);
        jacobi_solver.solve(u, iters);
        std::cout << iters << std::endl;
        jacobi_iters.push_back(iters);
        seidel_solver.solve(u, iters);
        std::cout << iters << std::endl;
        seidel_iters.push_back(iters);
        std::cout << "chi_continuous overlap = " << overlap << std::endl;
    }
    overlap_iter_csv_writer("jacobi_chi_continuous_overlap", overlaps, jacobi_iters);
    overlap_iter_csv_writer("seidel_chi_continuous_overlap", overlaps, seidel_iters);
    jacobi_iters.clear();
    seidel_iters.clear();

    return 0;
}