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


TEMPLATE_TEST_CASE("Schwarz methods for four squares on 100x100 grid with overlap 31h",
                   "[examples][schwarz][heavy_calculation][4_subs]", JacobiSchwarzSolver, OriginalSchwarzSolver) {
    size_t N = 100;
    double h = 1.0 / static_cast<double>(N);
    std::vector<size_t> mask = create_block_mask(N, 2, 2);
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
    solver.set_overlap(16);
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
            VtkWriter::add_point_data(mask, "mask", fs);
            for(size_t i = 0; i < solver.get_subdomains().size(); i++)
              VtkWriter::add_point_data(solver.get_subdomains()[i]->get_mask(), std::to_string(i), fs);
            fs.close();
        }
        if (force)
            break;

        solver.iterate(u);
        double current_error = norm_inf(diff_of(u, u_exact));
        if (std::abs(error - current_error) < 1e-11) {
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
