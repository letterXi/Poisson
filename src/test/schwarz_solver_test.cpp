#include "main_test.hpp"
#include "poisson_problem_solver/schwarz_methods/schwarz_solver.hpp"
#include <unordered_map>
#include <vector>

inline double source_function(double x, double y) {
    return x - x + y - y + 1;
}

inline double boundary_function(double x, double y) {
    return x - x + y - y - 1;
}

class SchwarzSolverTest : public SchwarzSolver {
public:
    SchwarzSolverTest(size_t N, std::vector<size_t> mask, std::function<double(double, double)> source_function,
                      std::function<double(double, double)> boundary_function)
        : SchwarzSolver(N, mask, source_function, boundary_function, 1000, 1e-6) {};
    void iterate(std::vector<double>& u) override {
        u[0] = u[0];
    }
    std::string get_name() const override {
        return "Schwarz";
    }
};

TEST_CASE("SchwarzSolver") {
    // 6 6 6 0 0 0 0 0
    // 6 0 0 0 0 1 1 0
    // 6 0 0 1 0 1 1 0
    // 0 0 1 1 1 0 0 0
    // 0 0 0 1 0 0 0 0
    // 0 0 0 0 0 0 4 0
    // 0 2 0 0 0 4 4 0
    // 0 0 0 0 0 4 4 0

    size_t N = 8;
    std::vector<size_t> mask = {0, 0, 0, 0, 0, 4, 4, 0, 0, 2, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 0,
                                4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 6, 0, 0, 1,
                                0, 1, 1, 0, 6, 0, 0, 0, 0, 1, 1, 0, 6, 6, 6, 0, 0, 0, 0, 0};

    std::vector<std::vector<size_t>> sub_masks = {{0,  1,  2,  3,  4,  7,  8,  10, 11, 12, 15, 16, 17, 18, 19,
                                                   20, 21, 23, 24, 25, 26, 28, 29, 30, 31, 32, 33, 37, 38, 39,
                                                   41, 42, 44, 47, 49, 50, 51, 52, 55, 59, 60, 61, 62, 63},
                                                  {5, 6, 13, 14, 22},
                                                  {9},
                                                  {27, 34, 35, 36, 43, 45, 46, 53, 54},
                                                  {40, 48, 56, 57, 58}};

    std::unordered_map<size_t, size_t> correct_id_to_idx;
    std::vector<size_t> correct_unique_indices = {0, 4, 2, 1, 6};

    for (size_t i = 0; i < correct_unique_indices.size(); i++)
        correct_id_to_idx[correct_unique_indices[i]] = i;

    SchwarzSolverTest solver(N, mask, source_function, boundary_function);
    std::unordered_map<size_t, size_t> id_to_idx = solver.get_id_to_idx();
    std::vector<size_t> unique_indices = solver.find_unique_indices();

    SECTION("SchwarzSolver::find_unique_indices()") {

        CHECK(unique_indices == correct_unique_indices);
    }

    // check id_to_idx
    SECTION("SchwarzSolver::id_to_idx") {
        CHECK(id_to_idx == correct_id_to_idx);
    }

    SECTION("SchwarzSolver::subdomains") {
        for (size_t i = 0; i < unique_indices.size(); i++)
            CHECK(solver.get_subdomains()[id_to_idx[unique_indices[i]]]->get_indices() == sub_masks[i]);
    }
    SECTION("SchwarzSolver::overlap") {
        std::vector<double> correct_sub_1_weights = {
            1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0, 1.0 / 2.0, 1.0 / 4.0, 1.0 / 2.0,
            1.0 / 3.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0,
            1.0 / 3.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        solver.set_overlap(3);
        CHECK(solver.get_subdomains()[1]->get_swap_index() == 18);
        std::vector<double> sub_1_weights = solver.get_subdomains()[1]->get_weights();
        CHECK(sub_1_weights == correct_sub_1_weights);
    }
}
