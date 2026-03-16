#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/original_schwarz_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include "test_functions.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

TEST_CASE("Comparison of the performance of parallel and sequential Jacobi methods",
          "[examples][jacobi][heavy_calculation][parallel]") {
    size_t N = 600;
    double h = 1.0 / static_cast<double>(N);
    std::vector<size_t> mask = create_block_mask(N, 2, 2);
    std::vector<double> u_0(N * N);
    std::vector<double> u(N * N);
    std::vector<double> u_exact(N * N);

    std::vector<VtkWriter::Point> points(N * N);
    for (size_t j = 0; j < N; j++) {
        for (size_t i = 0; i < N; i++) {
            double x = static_cast<double>(i) * h;
            double y = static_cast<double>(j) * h;
            size_t k = i + j * N;

            u_0[k] = static_cast<double>(mask[k]);
            u_exact[k] = exactFunction(x, y);
            points[k] = {x, y, 0.0};
        }
    }

    JacobiSchwarzSolver solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    OriginalSchwarzSolver schwarz_solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(100);
    schwarz_solver.set_overlap(100);
    schwarz_solver.initialize(u_0);

    size_t iters;

    solver.initialize(u_0);
    u = u_0;
    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(u, iters);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Jacobi seq_solve time: " << elapsed.count() << " ms" << std::endl;

    solver.initialize(u_0);
    u = u_0;
    start = std::chrono::high_resolution_clock::now();
    solver.parallel_solve(u, iters);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Jacobi par_solve time: " << elapsed.count() << " ms" << std::endl;

    u = u_0;
    start = std::chrono::high_resolution_clock::now();
    schwarz_solver.solve(u, iters);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "OriginalSchwarz solve time: " << elapsed.count() << " ms" << std::endl;
}
