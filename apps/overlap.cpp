#include "poisson_problem/ddm/jacobi_solve.hpp"
#include "poisson_problem/ddm/local_solve.hpp"
#include "poisson_problem/grid/grid.hpp"
#include "poisson_problem/vtk/vtk_save.hpp"
#include <fstream>

int main() {
    size_t N = 101;
    double h = 1.0 / static_cast<double>(N - 1);
    Grid grid(0.0, 0.0, N, 2, h);

    JacobiSolve solve(grid, 99, 1e-10, "chi_const");

    LocalSolve left = solve.getLeftLocalSolve();
    LocalSolve right = solve.getRightLocalSolve();

    std::ofstream left_file("left_chi.csv");
    std::ofstream right_file("right_chi.csv");

    left_file << "x,chi_left" << std::endl;
    right_file << "x,chi_right" << std::endl;

    for (size_t i = 0; i < N; i++) {
        double x = grid.getX(i);
        double y = 0;

        left_file << x << ',' << left.chiContinuous(x, y) << std::endl;
        right_file << x << ',' << right.chiContinuous(x, y) << std::endl;
    }

    left_file.close();
    right_file.close();

    return 0;
}