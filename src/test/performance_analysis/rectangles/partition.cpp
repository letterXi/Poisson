#include "../../examples/test_functions.hpp"
#include "../../main_test.hpp"
#include "poisson_problem_solver/grid/grid_regular2d.hpp"
#include "poisson_problem_solver/schwarz_methods/jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/original_schwarz_solver.hpp"
#include "poisson_problem_solver/utils/make_masks.hpp"

TEMPLATE_TEST_CASE("rectangle partition", "[performance]", JacobiSchwarzSolver, OriginalSchwarzSolver) {
    size_t N = 500;
    size_t iters;
    RegularGrid2D grid(0, 0, 1, 1, N, N);
    std::vector<double> u(grid.npoints(), 0);
    std::vector<double> u_exact(grid.npoints());
    fill_u_exact(u_exact, grid);
    std::vector<size_t> mask = stripes_mask(grid, 2, 30);
    TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(40);
    std::cout << solver.overlap_ratio() << "%" << std::endl;
    solver.initialize(u);
    solver.solve(u, u_exact, iters);
    std::cout << solver.get_name() << ' '<< iters << std::endl;
}
