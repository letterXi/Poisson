#include "poisson_problem/grid/grid.hpp"
#include <cmath>
#include <utility>

Grid::Grid(double x_0, double y_0, size_t N_x, size_t N_y, double h)
    : x_0_(x_0), y_0_(y_0), N_x_(N_x), N_y_(N_y), h_(h) {}
double Grid::getX(size_t i) const {
    return x_0_ + static_cast<double>(i) * h_;
}
double Grid::getY(size_t j) const {
    return y_0_ + static_cast<double>(j) * h_;
}

size_t Grid::getK(size_t i, size_t j) const {
    return i + j * N_x_;
}


size_t Grid::get_N_x() const {
    return N_x_;
}
size_t Grid::get_N_y() const {
    return N_y_;
}

double Grid::get_h() const {
    return h_;
}

size_t Grid::points() const {
    return N_x_ * N_y_;
}

size_t Grid::getI(double x) const {
    return static_cast<size_t>(std::round((x - x_0_) / h_));
}

size_t Grid::getJ(double y) const {
    return static_cast<size_t>(std::round((y - y_0_) / h_));
}


bool Grid::is_boundary(size_t i, size_t j) const {
    if (i == 0 || j == 0 || i == N_x_ - 1 || j == N_y_ - 1) {
        return true;
    }

    return false;
}

bool Grid::is_in_domain(double x, double y) const
{
    if (x < getX(0) || x > getX(N_x_ - 1) || y < getY(0) || y > getY(N_y_ - 1) )
            return false;
    return true;
}
