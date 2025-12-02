#include "AmgclSolver.hpp"
#include "CsrMatrix.hpp"
#include "Poisson.hpp"
#include "shwartz.hpp"
#include <iostream>

localSolve::localSolve(size_t id, const Mesh& mesh, size_t intersection)
    : id_(id), mesh_(mesh), intersection_(intersection), u_(mesh_.points(), 0.1), boundary_(mesh_.get_N_y(), 0.1) {}
void localSolve::give_boundary(localSolve& other) const {
    other.set_boundary(slice(intersection_));
}

std::vector<double> localSolve::slice(size_t i) const {
    std::vector<double> u_slice;
    for (size_t j = 0; j < mesh_.get_N_y(); j++)
        u_slice.push_back(u_[mesh_.getK(i, j)]);
    return u_slice;
}
void localSolve::set_boundary(std::vector<double> boundary) {
    boundary_ = boundary;
}

void localSolve::solve() {
    std::vector<size_t> addr;
    std::vector<size_t> cols;
    std::vector<double> vals;
    std::vector<double> rhs;

    Poisson(
        mesh_, sourceFunction,
        [this](double x, double y) {
            if ((x - mesh_.getX(intersection_)) < 1e-10) {
                size_t j = static_cast<size_t>((y - mesh_.getY(0)) / static_cast<double>(mesh_.get_h()));
                if (j > 0 && j < mesh_.get_N_y() - 1) {
                    return this->boundary_[j];
                }
            }

            return dirichletBoundaryFunction(x, y);
        },
        addr, cols, vals, rhs);
    CsrMatrix mat(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver({{"solver.type", "bicgstab"}});
    solver.set_matrix(mat);
    solver.solve(rhs, u_);
}
std::vector<double> localSolve::get_solve() {
    return u_;
}

Mesh::Mesh(double x_0, double y_0, size_t N_x, size_t N_y, double h)
    : x_0_(x_0), y_0_(y_0), N_x_(N_x), N_y_(N_y), h_(h) {}
double Mesh::getX(size_t i) const {
    return x_0_ + static_cast<double>(i) * h_;
}
double Mesh::getY(size_t j) const {
    return y_0_ + static_cast<double>(j) * h_;
}

size_t Mesh::getK(size_t i, size_t j) const {
    return i + j * N_x_;
}

std::pair<double, double> Mesh::getPoint(size_t i, size_t j) const {
    return std::pair<double, double>(getX(i), getY(j));
}

size_t Mesh::get_N_x() const {
    return N_x_;
}
size_t Mesh::get_N_y() const {
    return N_y_;
}

double Mesh::get_h() const {
    return h_;
}

size_t Mesh::points() const {
    return N_x_ * N_y_;
}

bool Mesh::is_boundary(size_t i, size_t j) const {
    if (i == 0 || j == 0 || i == N_x_ - 1 || j == N_y_ - 1) {
        return true;
    }

    return false;
}
