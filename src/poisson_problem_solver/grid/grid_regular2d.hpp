#pragma once

#include <vector>

struct Point {
    double x, y, z;
};

class RegularGrid2D {
public:
    RegularGrid2D(double x0, double y0, double x1, double y1, size_t nx, size_t ny) {
        nx_ = nx;
        ny_ = ny;
        points_.reserve(nx_ * ny_);
        double hx = (x1 - x0) / static_cast<double>(nx - 1);
        double hy = (y1 - y0) / static_cast<double>(ny - 1);
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                points_.push_back({x0 + i * hx, y0 + j * hy, 0});
            }
        }
    }

    Point get_point(size_t i, size_t j) const {
        return points_[i + j * nx_];
    }

    Point get_point(size_t k) const {
        return points_[k];
    }

    const std::vector<Point>& get_points() const {
        return points_;
    }

    size_t nx() const {
        return nx_;
    }
    size_t ny() const {
        return ny_;
    }

    size_t npoints() const
    {
      return points_.size();
    }

    double hx() const {
        return (points_[1].x - points_[0].x);
    }

    double hy() const {
        return (points_[1].y - points_[0].y);
    }

private:
    size_t nx_, ny_;
    std::vector<Point> points_;
};
