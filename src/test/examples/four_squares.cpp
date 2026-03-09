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
static std::vector<size_t> create_grid_mask(size_t N, size_t k) {
    std::vector<size_t> mask(N * N);

    size_t block_size = N / k;

    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < N; ++i) {
            size_t block_row = j / block_size;
            size_t block_col = i / block_size;

            if (block_row >= k)
                block_row = k - 1;
            if (block_col >= k)
                block_col = k - 1;

            mask[i + j * N] = block_row * k + block_col;
        }
    }
    return mask;
}

TEMPLATE_TEST_CASE("Schwarz methods for four squares on 100x100 grid with overlap 29h", "[examples][schwarz]",
                   JacobiSchwarzSolver, OriginalSchwarzSolver) {
    size_t N = 100;
    double h = 1.0 / static_cast<double>(N);
    std::vector<size_t> mask = create_grid_mask(N, 2);
    std::vector<double> u(N * N);
    std::vector<double> u_exact(N * N);

    std::vector<VtkWriter::Point> points(N * N);
    for (size_t j = 0; j < N; j++) {
        for (size_t i = 0; i < N; i++) {
            double x = static_cast<double>(i) * h;
            double y = static_cast<double>(j) * h;
            size_t k = i + j * N;

            u[k] = static_cast<double>(mask[k]);
            u_exact[k] = exactFunction(x, y);
            points[k] = {x, y, 0.0};
        }
    }

    TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(15);
    solver.initialize(u);

    std::string stem = "four_squares_" + solver.get_name();
    VtkWriter::StepManager schwarz_filename_mgr(stem, 1);

    bool force = false;
    double error = std::numeric_limits<double>::max();
    size_t iters;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t iter = 0; iter < 1000; iter++) {
        std::string path = schwarz_filename_mgr.add(iter, force);
        if (!path.empty()) {
            std::ofstream fs(path);
            VtkWriter::append_header(stem, fs);
            VtkWriter::append_points(points, N, fs);
            VtkWriter::append_point_data_header(N * N, fs);
            VtkWriter::add_point_data(u, "U", fs);
            fs.close();
        }
        if (force)
            break;

        solver.iterate(u);
        double current_error = norm_inf(diff_of(u, u_exact));
        if (std::abs(error - current_error) < 1e-16) {
            force = true;
        }
        error = current_error;
        iters = iter;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << "==========" << stem << "==========" << std::endl;
    std::cout << "iters = " << iters << ", error = " << error << std::endl;
    std::cout << "Iteration time: " << elapsed.count() << " ms" << std::endl << std::endl;
}
