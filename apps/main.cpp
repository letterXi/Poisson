#include "poisson_problem/ddm/jacobi_solve.hpp"
#include "poisson_problem/ddm/seidel_solve.hpp"
#include "poisson_problem/grid/grid.hpp"
#include "poisson_problem/vtk/vtk_save.hpp"
#include <fstream>
#include <vector>

#include <iostream>
int main() {

    double x_0 = 0.0;
    double y_0 = 0.0;
    size_t N = 100;
    double h = 1.0 / static_cast<double>(N - 1);

    Grid grid(x_0, y_0, N, N, h);

    size_t overlap = 30;
    double tolerance = 1e-6;
    JacobiSolve jacobi_solver(grid, overlap, tolerance, "chi_const");
    SeidelSolve seidel_solver(grid, overlap, tolerance, "chi_const");

    std::vector<double> u_jacobi(grid.points());
    std::vector<double> u_seidel(grid.points());

    size_t iters_jacobi;
    size_t iters_seidel;

    jacobi_solver.solve(u_jacobi, iters_jacobi);
    seidel_solver.solve(u_seidel, iters_seidel);

    vtkWriter writer("ddm_solve", "Jacobi and Seidel ddm method", grid);

    writer.add_scalars(u_seidel, "U_seidel");
    writer.add_scalars(u_jacobi, "U_jacobi");

    writer.write();

    std::cout << iters_jacobi << ' ' << iters_seidel << std::endl;

    return 0;
}