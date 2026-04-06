#include "schwarz_solver.hpp"

class JacobiSchwarzSolver : public SchwarzSolver {

public:
    JacobiSchwarzSolver(size_t N, std::vector<size_t> mask, std::function<double(double, double)> source_function,
                        std::function<double(double, double)> boundary_function, size_t maxiter = 1000,
                        double tolerance = 1e-6)
        : SchwarzSolver(N, mask, source_function, boundary_function, maxiter, tolerance) {}

    void iterate(std::vector<double>& u) override;
    void parallel_iterate(std::vector<double>& u);
    void parallel_solve(std::vector<double>& u, size_t& iters);
    void parallel_solve(std::vector<double>& u, const std::vector<double>& u_exact, size_t& iters);
    std::string get_name() const override;
};
