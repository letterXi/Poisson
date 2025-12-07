#include "poisson_problem/slae_build/poisson_slae_build.hpp"
#include <algorithm>
#include <cmath>
#include <execution>
#include <iostream>
#include <ranges>
#include <thread>

double exactFunction(double x, double y) {
    return 0.2 * (std::sin(2 * M_PI * x) - std::cos(12 * M_PI * y) + std::tan(x * y) + std::exp(-x * y));
}

double dirichletBoundaryFunction(double x, double y) {
    return exactFunction(x, y);
}

double sourceFunction(double x, double y) {
    return -0.2 * (-4 * M_PI * M_PI * std::sin(2 * M_PI * x) + 144 * M_PI * M_PI * std::cos(12 * M_PI * y) +
                   2 * (x * x + y * y) * (1.0 / std::cos(x * y)) * (1.0 / std::cos(x * y)) * std::tan(x * y) +
                   (x * x + y * y) * std::exp(-x * y));
}

void Poisson(const Grid& grid, std::function<double(double, double)> source_function,
             std::function<double(size_t, size_t)> boundary_function, std::vector<size_t>& addr,
             std::vector<size_t>& cols, std::vector<double>& vals, std::vector<double>& rhs) {
    size_t n2 = grid.points();
    size_t inner_n = (grid.get_N_x() - 2) * (grid.get_N_y() - 2);
    size_t bound_n = n2 - inner_n;

    addr.clear();
    addr.reserve(n2 + 1);
    addr.push_back(0);
    size_t temp = 0;
    for (size_t j = 0; j < grid.get_N_y(); j++) {
        for (size_t i = 0; i < grid.get_N_x(); i++) {
            if (grid.is_boundary(i, j))
                temp++;
            else
                temp += 5;

            addr.push_back(temp);
        }
    }

    cols.clear();
    cols.resize(bound_n + inner_n * 5);

    vals.clear();
    vals.resize(bound_n + inner_n * 5);

    rhs.resize(n2);

    double h = grid.get_h();
    size_t n = grid.get_N_x();
    for (size_t j = 0; j < grid.get_N_y(); j++) {
        for (size_t i = 0; i < grid.get_N_x(); ++i) {
            size_t k = grid.getK(i, j);
            size_t ind = addr[k];
            if (grid.is_boundary(i, j)) {
                cols[ind] = (k);
                vals[ind++] = (1);
                rhs[k] = boundary_function(i, j);
            } else {
                cols[ind] = (k - n);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k - 1);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k);
                vals[ind++] = (4.0 / (h * h));

                cols[ind] = (k + 1);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k + n);
                vals[ind++] = (-1.0 / (h * h));

                rhs[k] = source_function(grid.getX(i), grid.getY(j));
            }
        }
    }
}
/*
void Poisson(int n, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
             std::vector<double>& rhs, size_t thd) {
    size_t n2 = n * n;
    size_t inner_n = (n - 2) * (n - 2);
    size_t bound_n = n2 - inner_n;

    addr.clear();
    addr.reserve(n2 + 1);
    addr.push_back(0);
    size_t temp = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            if (i == 0 || i == n - 1 || j == 0 || j == n - 1)
                temp++;
            else
                temp += 5;

            addr.push_back(temp);
        }
    }

    cols.clear();
    cols.resize(bound_n + inner_n * 5);

    vals.clear();
    vals.resize(bound_n + inner_n * 5);

    rhs.resize(n2);
    double h = 1.0 / (double)(n - 1);

    if (thd == 0) {
        for (int j = 0; j < n; j++)
            full_mat_vectors(n, h, std::ref(addr), std::ref(cols), std::ref(vals), std::ref(rhs), j);
    } else if (thd == 1) {
        auto j_inds = std::views::iota(0, n);
        std::for_each(std::execution::par, j_inds.begin(), j_inds.end(), [n, h, &addr, &cols, &vals, &rhs](int j) {
            full_mat_vectors(n, h, std::ref(addr), std::ref(cols), std::ref(vals), std::ref(rhs), j);
        });
    } else {
        size_t pivot = n / thd;
        std::vector<std::thread> threads;

        for (size_t i = 0; i < thd; i++) {
            size_t start_j = pivot * i;
            size_t end_j = (i == thd - 1) ? n : pivot * (i + 1);

            threads.push_back(std::thread([n, h, &addr, &cols, &vals, &rhs, start_j, end_j]() {
                for (size_t j = start_j; j < end_j; j++)
                    full_mat_vectors(n, h, std::ref(addr), std::ref(cols), std::ref(vals), std::ref(rhs),
                                     static_cast<int>(j));
            }));
        }
        for (size_t i = 0; i < thd; i++)
            threads[i].join();
    }
}

void full_mat_vectors(int n, double h, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
                      std::vector<double>& rhs, int j) {
    for (int i = 0; i < n; ++i) {
        int k = i + j * n;
        size_t ind = addr[k];
        if (i == 0 || j == 0 || i == n - 1 || j == n - 1) {
            cols[ind] = (k);
            vals[ind++] = (1);
            rhs[k] = dirichletBoundaryFunction(i * h, j * h);
        } else {
            cols[ind] = (k - n);
            vals[ind++] = (-1.0 / (h * h));

            cols[ind] = (k - 1);
            vals[ind++] = (-1.0 / (h * h));

            cols[ind] = (k);
            vals[ind++] = (4.0 / (h * h));

            cols[ind] = (k + 1);
            vals[ind++] = (-1.0 / (h * h));

            cols[ind] = (k + n);
            vals[ind++] = (-1.0 / (h * h));

            rhs[k] = sourceFunction(i * h, j * h);
        }
    }
}*/
PoissonSlaeBuilder::PoissonSlaeBuilder(const Grid& grid) {
    std::vector<size_t> addr;
    std::vector<size_t> cols;
    std::vector<double> vals;

    size_t n2 = grid.points();
    size_t inner_n = (grid.get_N_x() - 2) * (grid.get_N_y() - 2);
    size_t bound_n = n2 - inner_n;

    addr.clear();
    addr.reserve(n2 + 1);
    addr.push_back(0);
    size_t temp = 0;
    for (size_t j = 0; j < grid.get_N_y(); j++) {
        for (size_t i = 0; i < grid.get_N_x(); i++) {
            if (grid.is_boundary(i, j))
                temp++;
            else
                temp += 5;

            addr.push_back(temp);
        }
    }

    cols.clear();
    cols.resize(bound_n + inner_n * 5);

    vals.clear();
    vals.resize(bound_n + inner_n * 5);

    double h = grid.get_h();
    size_t n = grid.get_N_x();
    for (size_t j = 0; j < grid.get_N_y(); j++) {
        for (size_t i = 0; i < grid.get_N_x(); ++i) {
            size_t k = grid.getK(i, j);
            size_t ind = addr[k];
            if (grid.is_boundary(i, j)) {
                cols[ind] = (k);
                vals[ind++] = (1);
            } else {
                cols[ind] = (k - n);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k - 1);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k);
                vals[ind++] = (4.0 / (h * h));

                cols[ind] = (k + 1);
                vals[ind++] = (-1.0 / (h * h));

                cols[ind] = (k + n);
                vals[ind++] = (-1.0 / (h * h));
            }
        }
    }

    mat_ = CsrMatrix(std::move(addr), std::move(cols), std::move(vals));
}

const CsrMatrix& PoissonSlaeBuilder::get_mat() const {
    return mat_;
}

std::vector<double> PoissonSlaeBuilder::get_rhs(const Grid& grid_,
                                                std::function<double(double, double)> source_function,
                                                std::function<double(double, double)> boundary_function) {
    rhs_ = std::vector<double>(grid_.points(), 0.0);
    for (size_t j = 0; j < grid_.get_N_y(); j++) {
        for (size_t i = 0; i < grid_.get_N_x(); i++) {
            double x = grid_.getX(i);
            double y = grid_.getY(j);
            if (grid_.is_boundary(i, j)) {
                rhs_[grid_.getK(i, j)] = boundary_function(x, y);
            } else {
                rhs_[grid_.getK(i, j)] = source_function(x, y);
            }
        }
    }
    return rhs_;
}
