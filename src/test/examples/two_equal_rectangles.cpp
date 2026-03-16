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

TEMPLATE_TEST_CASE("Schwarz methods for two equal rectangles on 100x100 grid with overlap 29h",
                   "[examples][schwarz][two_rect][heavy_calculation]", OriginalSchwarzSolver, JacobiSchwarzSolver) {
    size_t N = 100;
    double h = 1.0 / static_cast<double>(N);
    std::vector<size_t> mask(N * N);
    std::vector<double> u(N * N);
    std::vector<double> u_exact(N * N);

    std::vector<VtkWriter::Point> points(N * N);
    for (size_t j = 0; j < N; j++) {
        for (size_t i = 0; i < N; i++) {
            double x = static_cast<double>(i) * h;
            double y = static_cast<double>(j) * h;
            size_t k = i + j * N;
            if (i >= N / 2)
                mask[k] = 1;
            u[k] = 100.0 * std::sin(static_cast<double>(i) + std::cos(static_cast<double>(j)));
            u_exact[k] = exactFunction(x, y);
            points[k] = {x, y, 0.0};
        }
    }

    TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(15); // 15*2 - 1
    solver.initialize(u);

    std::string stem = "2_eq_rect_" + solver.get_name();
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
            VtkWriter::add_point_data(mask, "mask", fs);
            for (size_t i = 0; i < solver.get_subdomains().size(); i++)
                VtkWriter::add_point_data(solver.get_subdomains()[i]->get_mask(), std::to_string(i), fs);
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
