#ifndef GRID_HPP
#define GRID_HPP
#include <cstddef>
class Grid {
public:
    Grid(double x_0, double y_0, size_t N_x, size_t N_y, double h);
    double getX(size_t i) const;
    double getY(size_t j) const;
    size_t getI(double x) const;
    size_t getJ(double y) const;
    size_t getK(size_t i, size_t j) const;
    size_t get_N_x() const;
    size_t get_N_y() const;
    double get_h() const;
    size_t points() const;
    bool is_boundary(size_t i, size_t j) const;
    bool is_in_domain(double x, double y) const;

private:
    double x_0_;
    double y_0_;
    size_t N_x_;
    size_t N_y_;
    double h_;
};
#endif 
