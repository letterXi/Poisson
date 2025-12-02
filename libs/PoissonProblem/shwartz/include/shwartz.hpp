#ifndef SHWARTZ_HPP
#define SHWARTZ_HPP
#include <utility>
#include <vector>

class Mesh {
public:
    Mesh(double x_0, double y_0, size_t N_x, size_t N_y, double h);
    double getX(size_t i) const;
    double getY(size_t j) const;
    size_t getK(size_t i, size_t j) const;
    std::pair<double, double> getPoint(size_t i, size_t j) const;
    size_t get_N_x() const;
    size_t get_N_y() const;
    double get_h() const;
    size_t points() const;
    bool is_boundary(size_t i, size_t j) const;

private:
    double x_0_;
    double y_0_;
    size_t N_x_;
    size_t N_y_;
    double h_;
};
class localSolve {
public:
    localSolve(size_t id, const Mesh& mesh, size_t intersection);
    void give_boundary(localSolve& other) const;
    void set_boundary(std::vector<double> boundary);
    std::vector<double> slice(size_t i) const;
    void solve();
    std::vector<double> get_solve();

private:
    size_t id_;
    Mesh mesh_;
    size_t intersection_;
    std::vector<double> u_;
    std::vector<double> boundary_;
};

#endif
