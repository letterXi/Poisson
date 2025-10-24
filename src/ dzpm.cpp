#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/idrs.hpp>
#include <amgcl/solver/richardson.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

void createResultsFolder() {
    std::filesystem::create_directory("results");
}

void writeFrame(const std::vector<double>& T, int N, int frame) {
    std::string filename = "results/frame_" + std::to_string(frame) + ".txt";
    std::ofstream file(filename);

    file << "X,Y,T" << std::endl;
    double h = 1.0 / (N - 1);

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            int k = j * N + i;
            file << i * h << "," << j * h << "," << T[k] << std::endl;
        }
    }

    file.close();
}

class SparseMatrix {
public:
    SparseMatrix(std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val)
        : ptr_(ptr), col_(col), val_(val), n((int)ptr.size() - 1) {}
    auto getMat() {
        return std::tie(n, ptr_, col_, val_);
    }
    auto getN() {
        return n;
    }

private:
    std::vector<int> ptr_, col_;
    std::vector<double> val_;
    int n;
};

typedef amgcl::backend::builtin<double> Backend;

template<template<class> class Relaxation, typename IterativeSolver> class AmgSolver {
    typedef amgcl::make_solver<amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, Relaxation>,
                               IterativeSolver>
        Solver;

public:
    AmgSolver(SparseMatrix& mat, std::vector<double>& rhs)
        : sol(mat.getMat()), mat_(mat), rhs_(rhs), x(mat.getN(), 0.0) {}

    std::vector<double> getSolve() {
        return x;
    }

    void solve() {
        std::tie(iterations, res) = sol(rhs_, x);
    }

    void setRhs(std::vector<double>& rhs) {
        rhs_ = rhs;
    }

    double getRes() {
        return res;
    }

    double getIter() {
        return iterations;
    }

private:
    Solver sol;
    SparseMatrix mat_;
    std::vector<double> rhs_;
    std::vector<double> x;
    int iterations;
    double res;
};

void createRhs(int n, int m, std::vector<double>& T, std::vector<double>& rhs) {
    double tau = 1.0 / (m - 1);
    for (int j = 0, k = 0; j < n; j++) {
        for (int i = 0; i < n; i++, k++) {
            if (j == 0 || j == n - 1) {
                rhs[k] = T[k];
            } else if (i == 0 || i == n - 1) {
                rhs[k] = 0.0;
            } else {
                rhs[k] = T[k] / tau;
            }
        }
    }
}

void createEquations(int n, int m, std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val) {
    double h = 1.0 / (n - 1);
    double tau = 1.0 / (m - 1);

    int n2 = n * n;
    ptr.clear();
    ptr.reserve(n2 + 1);
    ptr.push_back(0);
    col.clear();
    col.reserve(n2 * 5);
    val.clear();
    val.reserve(n2 * 5);

    for (int j = 0, k = 0; j < n; j++) {
        for (int i = 0; i < n; i++, k++) {
            if (j == 0 || j == n - 1) {
                col.push_back(k);
                val.push_back(1.0);
            } else if (i == 0) {
                col.push_back(k);
                val.push_back(-1.0);
                col.push_back(k + 1);
                val.push_back(1.0);
            } else if (i == n - 1) {
                col.push_back(k - 1);
                val.push_back(-1.0);
                col.push_back(k);
                val.push_back(1.0);
            } else {
                col.push_back(k - n);
                val.push_back(-1.0 / (h * h));

                col.push_back(k - 1);
                val.push_back(-1.0 / (h * h));

                col.push_back(k);
                val.push_back(1.0 / tau + 4.0 / (h * h));

                col.push_back(k + 1);
                val.push_back(-1.0 / (h * h));

                col.push_back(k + n);
                val.push_back(-1.0 / (h * h));
            }
            ptr.push_back((int)col.size());
        }
    }
}

void fullT(std::vector<double>& T, int N) {
    T.resize(N * N);

    for (int j = 0, k = 0; j < N; j++) {
        for (int i = 0; i < N; i++, k++) {
            if (j == N - 1) {
                T[k] = 1.0;
            } else if (j == 0) {
                T[k] = 0.4;
            } else {
                T[k] = 1e-8;
            }
        }
    }
}

template<typename T> void output(std::vector<T> a) {
    for (size_t i = 0; i < a.size(); i++) {
        std::cout << a[i] << ' ';
    }
    std::cout << std::endl;
}

int main() {
    int N = 200;
    int M = 1000;

    std::vector<int> ptr, col;
    std::vector<double> val, rhs;

    rhs.resize(N * N);

    std::vector<double> T;
    fullT(T, N);

    createEquations(N, M, ptr, col, val);

    SparseMatrix mat(ptr, col, val);

    AmgSolver<amgcl::relaxation::ilu0, amgcl::solver::bicgstab<Backend>> mainsoler(mat, rhs);
    createResultsFolder();

    writeFrame(T, N, 0);

    for (int j = 0; j < M; j++) {
        createRhs(N, M, T, rhs);
        mainsoler.setRhs(rhs);
        mainsoler.solve();
        T = mainsoler.getSolve();

        if (j % 10 == 0) {
            writeFrame(T, N, j / 10 + 1);
            std::cout << "Записан кадр " << j / 10 + 1 << " (шаг " << j << ")"
                      << "" << mainsoler.getIter() << std::endl;
        }
        if ((int)mainsoler.getIter() < 1) {
            break;
        }
    }
    writeFrame(T, N, M / 10 + 1);
    std::cout << "Записан финальный кадр" << std::endl;

    std::cout << "=== ТЕСТИРОВАНИЕ РЕШАТЕЛЕЙ ===" << std::endl;

    // SPAI0 со всеми решателями
    std::cout << "\n--- SPAI0 ---" << std::endl;
    {
        AmgSolver<amgcl::relaxation::spai0, amgcl::solver::bicgstab<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "BiCGSTAB: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::spai0, amgcl::solver::cg<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "CG: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::spai0, amgcl::solver::gmres<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "GMRES: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::spai0, amgcl::solver::idrs<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "IDRS: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::spai0, amgcl::solver::richardson<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "Richardson: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }

    // ILU0 со всеми решателями
    std::cout << "\n--- ILU0 ---" << std::endl;
    {
        AmgSolver<amgcl::relaxation::ilu0, amgcl::solver::bicgstab<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "BiCGSTAB: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::ilu0, amgcl::solver::cg<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "CG: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::ilu0, amgcl::solver::gmres<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "GMRES: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::ilu0, amgcl::solver::idrs<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "IDRS: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::ilu0, amgcl::solver::richardson<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "Richardson: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }

    // Gauss-Seidel со всеми решателями
    std::cout << "\n--- Gauss-Seidel ---" << std::endl;
    {
        AmgSolver<amgcl::relaxation::gauss_seidel, amgcl::solver::bicgstab<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "BiCGSTAB: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::gauss_seidel, amgcl::solver::cg<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "CG: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::gauss_seidel, amgcl::solver::gmres<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "GMRES: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::gauss_seidel, amgcl::solver::idrs<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "IDRS: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::gauss_seidel, amgcl::solver::richardson<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "Richardson: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }

    // Chebyshev со всеми решателями
    std::cout << "\n--- Chebyshev ---" << std::endl;
    {
        AmgSolver<amgcl::relaxation::chebyshev, amgcl::solver::bicgstab<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "BiCGSTAB: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::chebyshev, amgcl::solver::cg<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "CG: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::chebyshev, amgcl::solver::gmres<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "GMRES: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::chebyshev, amgcl::solver::idrs<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "IDRS: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::chebyshev, amgcl::solver::richardson<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "Richardson: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }

    // Damped Jacobi со всеми решателями
    std::cout << "\n--- Damped Jacobi ---" << std::endl;
    {
        AmgSolver<amgcl::relaxation::damped_jacobi, amgcl::solver::bicgstab<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "BiCGSTAB: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::damped_jacobi, amgcl::solver::cg<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "CG: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::damped_jacobi, amgcl::solver::gmres<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "GMRES: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::damped_jacobi, amgcl::solver::idrs<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "IDRS: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }
    {
        AmgSolver<amgcl::relaxation::damped_jacobi, amgcl::solver::richardson<Backend>> solv(mat, rhs);
        solv.solve();
        std::cout << "Richardson: " << solv.getIter() << ' ' << solv.getRes() << std::endl;
    }

    return 0;
}
