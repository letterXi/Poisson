#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/original_schwarz_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
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

TEMPLATE_TEST_CASE("Overlap convergence study", "[overlap]", OriginalSchwarzSolver, JacobiSchwarzSolver) {
    size_t N = 100;
    double h = 1.0 / static_cast<double>(N);
    std::vector<size_t> mask = create_grid_mask(N, 4, 4);
    std::vector<double> u(N * N);
    std::vector<double> u_0(N * N);
    std::vector<double> u_exact(N * N);

    for (size_t j = 0; j < N; j++) {
        for (size_t i = 0; i < N; i++) {
            double x = static_cast<double>(i) * h;
            double y = static_cast<double>(j) * h;
            size_t k = i + j * N;

            u_0[k] = 100.0 * std::sin(static_cast<double>(i) + std::cos(static_cast<double>(j)));
            u_exact[k] = exactFunction(x, y);
        }
    }

    size_t iters = 0;
    for (size_t o = 1; o <= N; o++) {
        TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction);
        solver.set_overlap(o);
        solver.initialize(u_0);
        u = u_0;

        std::ofstream file;
        if (o == 1) {
            file.open(solver.get_name() + "_convergence.csv");
            file << "overlap,iters" << std::endl;
        } else {
            file.open(solver.get_name() + "_convergence.csv", std::ios::app);
        }
        double error = std::numeric_limits<double>::max();

        for (size_t iter = 0; iter < 10000; iter++) {
            solver.iterate(u);
            double current_error = norm_inf(diff_of(u, u_exact));
            if (std::abs(error - current_error) < 1e-16) {
                iters = iter;
                break;
            }
            error = current_error;
        }

        file << o << ',' << iters << std::endl;
        file.close();
    }
}
