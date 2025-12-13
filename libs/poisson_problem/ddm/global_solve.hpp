#ifndef GLOBAL_SOLVE_HPP
#define GLOBAL_SOLVE_HPP
#include "poisson_problem/ddm/local_solve.hpp"
#include "poisson_problem/grid/grid.hpp"
#include <memory>
#include <string>
#include <utility>
#include <vector>

class GlobalSolve {
public:
    virtual ~GlobalSolve() = default;
    GlobalSolve(const Grid& grid, size_t overlap, double tolerance, const std::string& glue_strategy,
                size_t maxit = 1000);
    virtual void solve(std::vector<double>& u, size_t& iters) = 0;
    const std::vector<std::pair<size_t, double>>& getConvergenceHistory() const;
    void changeOverlap(size_t overlap);
    void changeGlueStrategy(const std::string& strategy);
    void clearConvergenceHistory();
    void splitDomain(size_t overlap);
    LocalSolve getLeftLocalSolve() const;
    LocalSolve getRightLocalSolve() const;

protected:
    std::unique_ptr<LocalSolve> left_domain_solve_;
    std::unique_ptr<LocalSolve> right_domain_solve_;
    std::unique_ptr<Grid> global_grid_;
    std::vector<std::pair<size_t, double>> convergence_history_;
    size_t overlap_;
    double tolerance_;
    std::string glue_strategy_;
    size_t maxit_;

    void glueSolve(std::vector<double>& u, const std::string& strategy);
};

#endif
