#ifndef GLOBAL_SOLVE_HPP
#define GLOBAL_SOLVE_HPP
#include "poisson_problem/ddm/local_solve.hpp"
#include "poisson_problem/grid/grid.hpp"
#include <memory>
#include <vector>

class GlobalSolve {
public:
    virtual ~GlobalSolve() = default;
    void splitDomain(size_t overlap);
    std::vector<double> get_solve() const;
    GlobalSolve(const Grid& gird, size_t overlap);
    virtual void solve(std::vector<double>& u, size_t& iters) = 0;

protected:
    std::unique_ptr<LocalSolve> left_domain_solve_;
    std::unique_ptr<LocalSolve> right_domain_solve_;
    std::unique_ptr<Grid> global_grid_;
};
#endif
