#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include "test_functions.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>
#include "poisson_problem_solver/utils/make_masks.hpp"

TEMPLATE_TEST_CASE("Schwarz methods for two equal rectangles on 100x100 grid with overlap 29h",
                   "[examples][schwarz][two_rect][heavy_calculation]", OriginalSchwarzSolver, JacobiSchwarzSolver) {
    size_t N = 100;
    RegularGrid2D grid(0, 0, 1, 1, N, N);
    std::vector<size_t> mask = block_mask(grid, 1, 2);
    std::vector<double> u(N * N, 0);
    std::vector<double> u_exact(N * N);

    TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(15); // 15*2 - 1
    solver.initialize(u);
    std::string stem = "2_eq_rect_" + solver.get_name();
    solve_with_screenshot(stem, u, u_exact, grid, solver);
}
