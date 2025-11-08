#include "AmgclSolver.hpp"
#include "CsrMatrix.hpp"
#include "poisson.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

void saveVtkFile(std::string filename, std::vector<double> u, std::vector<double> u_ex, int n) {

    double h = 1.0 / (n - 1);
    std::ofstream vtkfile(filename);
    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << "Solve Poisson problem" << std::endl;
    vtkfile << "ASCII" << std::endl << std::endl;
    vtkfile << "DATASET STRUCTURED_GRID" << std::endl;
    vtkfile << "DIMENSIONS" << ' ' << n << ' ' << n << ' ' << 1 << std::endl;
    vtkfile << "POINTS" << ' ' << n * n << ' ' << "float" << std::endl;
    if (!vtkfile.is_open())
        std::cout << "no open" << std::endl;
    else {
        for (int j = 0; j < n; j++) {
            double y = j * h;
            for (int i = 0; i < n; i++) {
                double x = i * h;
                vtkfile << x << ' ' << y << ' ' << 0.0 << std::endl;
            }
        }
        vtkfile << std::endl << "POINT_DATA" << ' ' << n * n << std::endl;

        vtkfile << "SCALARS U float" << std::endl;
        vtkfile << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < u.size(); ++i)
            vtkfile << u[i] << std::endl;

        vtkfile << "SCALARS U_exact float" << std::endl;
        vtkfile << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < u_ex.size(); ++i)
            vtkfile << u_ex[i] << std::endl;

        vtkfile.close();
    }
}

int main() {

    std::vector<size_t> addr;
    std::vector<size_t> cols;
    std::vector<double> vals;
    std::vector<double> rhs;
    int n = 200;

    Poisson(n, addr, cols, vals, rhs, 3);
    std::vector<double> u;

    CsrMatrix mat(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver({{"solver.type", "bicgstab"}});
    solver.set_matrix(mat);
    solver.solve(rhs, u);

    std::vector<double> exact_u;

    double h = 1.0 / n;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            exact_u.push_back(exactFunction(i * h, j * h));
        }
    }

    saveVtkFile("solve.vtk", u, exact_u, n);

    return 0;
}
