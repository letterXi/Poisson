#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include "poisson_problem_solver/utils/tictoc.hpp"
#include "test_functions.hpp"
#include <vector>
#include "poisson_problem_solver/utils/make_masks.hpp"

TEST_CASE("Comparison of the performance of parallel and sequential Jacobi methods",
          "[examples][jacobi][heavy_calculation][parallel]") {
    size_t N = 100; 
    RegularGrid2D grid(0,0,1,1, N,N);
    std::vector<size_t> mask = block_mask(grid, 2, 2);
    std::vector<double> u_0(N * N, 0.0);
    std::vector<double> u_exact(N*N);
    fill_u_exact(u_exact, grid);
    std::vector<double> u(N * N);
    JacobiSchwarzSolver solver(N, mask, sourceFunction, dirichletBoundaryFunction, 10000);
    OriginalSchwarzSolver schwarz_solver(N, mask, sourceFunction, dirichletBoundaryFunction, 10000);
    size_t overlap = 50;
    solver.set_overlap(overlap);
    schwarz_solver.set_overlap(overlap);

    size_t iters;

    solver.initialize(u_0);
    u = u_0;

    std::cout << solver.overlap_ratio() << "%" << std::endl;

    Tic("jacobi_seq");
    solver.solve(u, u_exact, iters);
    Toc("jacobi_seq");
    std::cout << "jacobi_seq iters = " << iters << std::endl;

    solver.initialize(u_0);
    u = u_0;

    Tic("jacobi_par");
    solver.parallel_solve(u, u_exact, iters);
    Toc("jacobi_par");
    std::cout << "jacobi_par iters = " << iters << std::endl;

    schwarz_solver.initialize(u_0);
    u = u_0;

    Tic("schwarz");
    schwarz_solver.solve(u, u_exact, iters);
    Toc("schwarz");
    std::cout << "schwarz iters = " << iters << std::endl;
}
