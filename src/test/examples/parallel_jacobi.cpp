#include "../main_test.hpp"
#include "poisson_problem_solver/mat_solver/csr_mat_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/jacobi_schwarz_solver.hpp"
#include "poisson_problem_solver/schwarz_methods/original_schwarz_solver.hpp"
#include "poisson_problem_solver/utils/norms.hpp"
#include "poisson_problem_solver/utils/vtk.hpp"
#include "poisson_problem_solver/utils/tictoc.hpp"
#include "test_functions.hpp"
#include <vector>

TEST_CASE("Comparison of the performance of parallel and sequential Jacobi methods",
          "[examples][jacobi][heavy_calculation][parallel]") {
    size_t N = 400; 
    std::vector<size_t> mask = create_block_mask(N, 2, 2);
    std::vector<double> u_0(N * N, 0.0);
    std::vector<double> u(N * N);
    JacobiSchwarzSolver solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    OriginalSchwarzSolver schwarz_solver(N, mask, sourceFunction, dirichletBoundaryFunction);
    solver.set_overlap(10);
    schwarz_solver.set_overlap(10);

    size_t iters;

    solver.initialize(u_0);
    u = u_0;

    Tic("jacobi_seq");
    solver.solve(u, iters);
    Toc("jacobi_seq");
    std::cout << "jacobi_seq iters = " << iters << std::endl;

    solver.initialize(u_0);
    u = u_0;

    Tic("jacobi_par");
    solver.parallel_solve(u, iters);
    Toc("jacobi_par");
    std::cout << "jacobi_par iters = " << iters << std::endl;

    schwarz_solver.initialize(u_0);
    u = u_0;

    Tic("schwarz");
    schwarz_solver.solve(u, iters);
    Toc("schwarz");
    std::cout << "schwarz iters = " << iters << std::endl;
}
