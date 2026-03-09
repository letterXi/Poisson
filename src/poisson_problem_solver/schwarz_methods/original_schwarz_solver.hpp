#include "schwarz_solver.hpp"

class OriginalSchwarzSolver : public SchwarzSolver {

public:
    OriginalSchwarzSolver(size_t N, std::vector<size_t> mask, std::function<double(double, double)> source_function,
                          std::function<double(double, double)> boundary_function, size_t maxiter = 1000, double tolerance = 1e-6)
        : SchwarzSolver(N, mask, source_function, boundary_function, maxiter, tolerance) {}

    void iterate(std::vector<double>& u) override;
    std::string get_name() const override;
};
