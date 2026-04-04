#include "../main_test.hpp"
#include "poisson_problem_solver/grid/grid_regular2d.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/original_schwarz_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include "test_functions.hpp"
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

bool is_ring(const Point& point, const Point& center, double r, double R) {
    double x = point.x;
    double y = point.y;
    double xc = center.x;
    double yc = center.y;
    return (((x - xc) * (x - xc) + (y - yc) * (y - yc) <= R * R) &&
            ((x - xc) * (x - xc) + (y - yc) * (y - yc) >= r * r));
}

TEMPLATE_TEST_CASE("Schwarz methods for ring subdomain on 500x500 grid with overlap 49h",
                   "[examples][schwarz][heavy_calculation][ring]", JacobiSchwarzSolver, OriginalSchwarzSolver) {
    size_t N = 100;
    RegularGrid2D grid(0, 0, 1, 1, N, N);
    std::vector<size_t> mask(N * N, 0);
    std::vector<double> u(N * N, 0);
    std::vector<double> u_exact(N * N);

    const auto& points = grid.get_points();
    for (size_t k = 0; k < points.size(); k++) {
        if (is_ring(points[k], {0.4, 0.7, 0}, 0.1, 0.3))
            mask[k] = 1;
    }

    TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(10);
    solver.initialize(u);

    std::string stem = "ring1_" + solver.get_name();
    solve_with_screenshot(stem, u, u_exact, grid, solver);
}
