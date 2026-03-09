#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/original_schwarz_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

static double exactFunction(double x, double y) {
    return 0.2 * (std::sin(2 * M_PI * x) - std::cos(12 * M_PI * y) + std::tan(x * y) + std::exp(-x * y));
}

static double dirichletBoundaryFunction(double x, double y) {
    return exactFunction(x, y);
}

static double sourceFunction(double x, double y) {
    return -0.2 * (-4 * M_PI * M_PI * std::sin(2 * M_PI * x) + 144 * M_PI * M_PI * std::cos(12 * M_PI * y) +
                   2 * (x * x + y * y) * (1.0 / std::cos(x * y)) * (1.0 / std::cos(x * y)) * std::tan(x * y) +
                   (x * x + y * y) * std::exp(-x * y));
}
static std::vector<size_t> create_grid_mask(size_t N, size_t rows, size_t cols) {
    std::vector<size_t> mask(N * N);

    double block_h = static_cast<double>(N) / rows;
    double block_w = static_cast<double>(N) / cols;

    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < N; ++i) {
            size_t r = static_cast<size_t>(j / block_h);
            size_t c = static_cast<size_t>(i / block_w);

            if (r >= rows)
                r = rows - 1;
            if (c >= cols)
                c = cols - 1;

            mask[i + j * N] = r * cols + c;
        }
    }
    return mask;
}

TEST_CASE("Comparison of the performance of parallel and sequential Jacobi methods",
          "[examples][jacobi][heavy_calculation][parallel]") {
    size_t N = 1000;
    double h = 1.0 / static_cast<double>(N);
    std::vector<size_t> mask = create_grid_mask(N, 3, 4);
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
