#include "AmgclSolver.hpp"
#include "CsrMatrix.hpp"
#include "Poisson.hpp"
#include "shwartz.hpp"
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

    int N = 100;
    int d = 0;
    int N1, N2;
    int h = 30;
    if (h % 2 == 0) {
        d = h / 2 - 1;
    } else if (h % 2 != 0)
        d = h / 2;
    if (N % 2 != 0) {
        N2 = N - N / 2 + 1 + d;
        N1 = N - N2 + h + 1;
    } else if (N % 2 == 0) {
        N1 = N / 2 + 1 + d;
        N2 = N - N1 + h + 1;
    }

    double h_ = 1.0 / (N - 1);
    double x01 = 0.0;
    double y01 = 0.0;
    double x02 = (N1 - h - 1.0) * h_;
    double y02 = 0.0;

    size_t in1 = N1 - 1 - h;
    size_t in2 = h; 

    std::cout << in1 << ' ' << in2 << std::endl;

    Mesh mesh1(x01, y01, N1, N, h_);
    Mesh mesh2(x02, y02, N2, N, h_); 
    
    localSolve solve1(0, mesh1, in1);
    localSolve solve2(0, mesh2, in2);  

    for(int i = 1; i < 10; i++)
    {
    solve2.give_boundary(solve1);
    solve1.give_boundary(solve2);

    solve1.solve();
    solve2.solve();
    }



    std::ofstream file1("N1.txt");
    std::ofstream file2("N2.txt");

    if (!file1.is_open() || !file2.is_open())
        std::cout << "dont open" << std::endl;
    file1 << "x,y,u" << std::endl;
    file2 << "x,y,u" << std::endl;

    std::cout << N1 << ' ' << N2 << std::endl;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N1; i++) {
            double x = x01 + h_ * i;
            double y = y01 + h_ * j;
            file1 << x << ',' << y << ',' << solve1.get_solve()[i + j*N1]<< std::endl;
        }
        for (int i = 0; i < N2; i++) {
            double x = x02 + h_ * i;
            double y = y02 + h_ * j;
            file2 << x << ',' << y << ',' << solve2.get_solve()[i + j*N2] << std::endl;
        }
    }

    file1.close();
    file2.close();
    // std::vector<size_t> addr;
    // std::vector<size_t> cols;
    // std::vector<double> vals;
    // std::vector<double> rhs;
    // int n = 100;
    //
    //   Poisson(n, addr, cols, vals, rhs, 3);
    //  std::vector<double> u;
    //
    //   CsrMatrix mat(std::move(addr), std::move(cols), std::move(vals));
    //  AmgclSolver solver({{"solver.type", "bicgstab"}});
    // solver.set_matrix(mat);
    // solver.solve(rhs, u);
    //
    //   std::vector<double> exact_u;
    //
    //   double h = 1.0 / n;
    //
    //   for (int j = 0; j < n; j++) {
    //      for (int i = 0; i < n; i++) {
    //         exact_u.push_back(exactFunction(i * h, j * h));
    //    }
    //}
    //
    //   saveVtkFile("solve.vtk", u, exact_u, n);

    return 0;
}
