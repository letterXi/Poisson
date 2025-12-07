#ifndef LOCAL_SOLVE_HPP
#define LOCAL_SOLVE_HPP
#include "poisson_problem/grid/grid.hpp"
#include "poisson_problem/slae_build/poisson_slae_build.hpp"
#include <vector>

class LocalSolve {
public:
    LocalSolve(size_t id, Grid grid, size_t intersection);
    ~LocalSolve() = default;
    void give_boundary(LocalSolve& other) const;
    void set_boundary(const std::vector<double>& boundary);
    std::vector<double> slice(size_t i) const;
    void solve();
    const std::vector<double>& get_solve() const;
    const Grid& get_grid() const;
    double chiContinuous(double x, double y) const;
    double chiConst(double x, double y) const;
    double expand(double x, double y) const;

private:
    size_t id_;
    Grid grid_;
    size_t intersection_;
    std::vector<double> u_;
    std::vector<double> boundary_;
    PoissonSlaeBuilder slae_;
};
#endif
