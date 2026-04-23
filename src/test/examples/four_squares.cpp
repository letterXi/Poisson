#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include "test_functions.hpp"
#include "poisson_problem_solver/grid/grid_regular2d.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>
#include "poisson_problem_solver/utils/make_masks.hpp"

TEMPLATE_TEST_CASE("Schwarz methods for four squares on 100x100 grid with overlap 31h",
                   "[examples][schwarz][heavy_calculation][4_subs]", JacobiSchwarzSolver, OriginalSchwarzSolver) {
    size_t N = 100;
    RegularGrid2D grid(0,0,1,1, N,N);
    std::vector<size_t> mask = block_mask(grid, 2, 2);
    std::vector<double> u(N * N, 0.0);
    std::vector<double> u_exact(N * N);
    fill_u_exact(u_exact, grid);

    TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(15);
    solver.initialize(u);

    std::string stem = "four_squares_" + solver.get_name();
    solve_with_screenshot(stem, u,u_exact, grid, solver); 

}
