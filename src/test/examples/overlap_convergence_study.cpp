#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/utils/make_masks.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include "test_functions.hpp"
#include <cmath>
#include <fstream>
#include <limits>
#include <vector>

TEMPLATE_TEST_CASE("Overlap convergence study", "[overlap][heavy_calculation]", OriginalSchwarzSolver,
                   JacobiSchwarzSolver) {
    size_t N = 200;
    RegularGrid2D grid(0, 0, 1, 1, N, N);
    std::vector<size_t> mask = block_mask(grid, 2, 2);
    std::vector<double> u(N * N);
    std::vector<double> u_0(N * N, 0);
    std::vector<double> u_exact(N * N);
    fill_u_exact(u_exact, grid);

    size_t iters = 0;
    for (size_t o = 1; o <= N*N; o++) {
        TestType solver(N, mask, sourceFunction, dirichletBoundaryFunction, 10000, 1e-11);
        solver.set_overlap(o);
        solver.initialize(u_0);
        u = u_0;

        std::ofstream file;
        if (o == 1) {
            file.open(solver.get_name() + "_convergence.csv");
            file << "overlap,iters,k,v" << std::endl;
        } else {
            file.open(solver.get_name() + "_convergence.csv", std::ios::app);
        }
        solver.solve(u, u_exact, iters);
        file << solver.overlap_ratio() << ',' << iters << ',' << solver.overhead_ratio() << ',' << solver.overhead_ratio() * iters << std::endl;
        file.close();
        if (solver.is_collapse())
          break;
    }
}
