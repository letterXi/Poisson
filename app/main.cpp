#include "AmgclSolver.hpp"
#include "CsrMatrix.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

double exactFunction(double x, double y) {
    return 0.2 * (std::sin(2 * M_PI * x) - std::cos(12 * M_PI * y) + std::tan(x * y) + std::exp(-x * y));
}

double dirichletBoundaryFunction(double x, double y) {
    return exactFunction(x, y);
}

double sourceFunction(double x, double y) {
    return -0.2 * (-4 * M_PI * M_PI * std::sin(2 * M_PI * x) + 144 * M_PI * M_PI * std::cos(12 * M_PI * y) +
                   2 * (x * x + y * y) * std::pow(1.0 / std::cos(x * y), 2) * std::tan(x * y) +
                   (x * x + y * y) * std::exp(-x * y));
}

void saveVtkFile(std::string filename, std::vector<double> u, int n) {

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
        vtkfile.close();
    }
}

void Poisson(int n, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
             std::vector<double>& rhs) {
    size_t n2 = n * n;
    double h = 1.0 / (double)(n - 1);

    addr.clear();
    addr.reserve(n2 + 1);
    addr.push_back(0);

    cols.clear();
    cols.reserve(n2 * 5);

    vals.clear();
    vals.reserve(n2 * 5);

    rhs.resize(n2);

    for (int j = 0, k = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i, ++k) {
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1) {
                cols.push_back(k);
                vals.push_back(1);
                rhs[k] = dirichletBoundaryFunction(i * h, j * h);
            } else {
                cols.push_back(k - n);
                vals.push_back(-1.0 / (h * h));

                cols.push_back(k - 1);
                vals.push_back(-1.0 / (h * h));

                cols.push_back(k);
                vals.push_back(4.0 / (h * h));

                cols.push_back(k + 1);
                vals.push_back(-1.0 / (h * h));

                cols.push_back(k + n);
                vals.push_back(-1.0 / (h * h));

                rhs[k] = sourceFunction(i * h, j * h);
            }
            addr.push_back(cols.size());
        }
    }
}

int main() {

    std::vector<size_t> addr;
    std::vector<size_t> cols;
    std::vector<double> vals;
    std::vector<double> rhs;
    int n = 100;
    Poisson(n, addr, cols, vals, rhs);

    std::vector<double> u;

    CsrMatrix mat(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver({{"solver.type", "bicgstab"}});
    solver.set_matrix(mat);
    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(rhs, u);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
    saveVtkFile("numsolve.vtk", u, n);

    n = 300;
    std::vector<double> exact_u;
    double h = 1.0 / (n - 1);

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            exact_u.push_back(exactFunction(i * h, j * h));
        }
    }

    saveVtkFile("exactsolve.vtk", exact_u, n);

    return 0;
}