#define CATCH_CONFIG_MAIN
#include "AmgclSolver.hpp"
#include "CsrMatrix.hpp"
#include <catch2/catch_all.hpp>
#include <vector>
using Catch::Approx;

TEST_CASE("AmgclSolver creation", "[AmgclSolver][constructor]") {
    SECTION("Default parameters") {
        AmgclSolver solver(100, 1e-6);
        REQUIRE(true);
    }

    SECTION("Custom parameters") {
        AmgclSolver solver({{"solver.type", "cg"}, {"solver.tol", "1e-8"}});
        REQUIRE(true);
    }
}

TEST_CASE("AmgclSolver simple 1x1 system", "[AmgclSolver][solve][1x1]") {
    // 5
    std::vector<size_t> addr = {0, 1};
    std::vector<size_t> cols = {0};
    std::vector<double> vals = {5.0};

    CsrMatrix matrix(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver(10, 1e-10);

    REQUIRE_NOTHROW(solver.set_matrix(matrix));

    std::vector<double> rhs = {10.0}; // 5*x = 10
    std::vector<double> solution;

    REQUIRE_NOTHROW(solver.solve(rhs, solution));
    REQUIRE(solution.size() == 1);
    REQUIRE(solution[0] == Approx(2.0)); // x = 2
}

TEST_CASE("AmgclSolver 2x2 diagonal system", "[AmgclSolver][solve][2x2][diagonal]") {
    // 2 0
    // 0 3
    std::vector<size_t> addr = {0, 1, 2};
    std::vector<size_t> cols = {0, 1};
    std::vector<double> vals = {2.0, 3.0};

    CsrMatrix matrix(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver(10, 1e-10);

    REQUIRE_NOTHROW(solver.set_matrix(matrix));

    std::vector<double> rhs = {4.0, 9.0}; // 2x=4, 3y=9
    std::vector<double> solution;

    REQUIRE_NOTHROW(solver.solve(rhs, solution));
    REQUIRE(solution.size() == 2);
    REQUIRE(solution[0] == Approx(2.0)); // x = 2
    REQUIRE(solution[1] == Approx(3.0)); // y = 3
}

TEST_CASE("AmgclSolver 2x2 full system", "[AmgclSolver][solve][2x2][full]") {
    // 2 1
    // 1 2
    std::vector<size_t> addr = {0, 2, 4};
    std::vector<size_t> cols = {0, 1, 0, 1};
    std::vector<double> vals = {2.0, 1.0, 1.0, 2.0};

    CsrMatrix matrix(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver(50, 1e-8);

    REQUIRE_NOTHROW(solver.set_matrix(matrix));

    std::vector<double> rhs = {5.0, 4.0}; // 2x + y = 5, x + 2y = 4
    std::vector<double> solution;

    REQUIRE_NOTHROW(solver.solve(rhs, solution));
    REQUIRE(solution.size() == 2);
    REQUIRE(solution[0] == Approx(2.0)); // x = 2
    REQUIRE(solution[1] == Approx(1.0)); // y = 1
}

TEST_CASE("AmgclSolver zero rhs", "[AmgclSolver][solve][zero-rhs]") {
    // 2 0
    // 0 2
    std::vector<size_t> addr = {0, 1, 2};
    std::vector<size_t> cols = {0, 1};
    std::vector<double> vals = {2.0, 2.0};

    CsrMatrix matrix(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver(10, 1e-10);

    REQUIRE_NOTHROW(solver.set_matrix(matrix));

    std::vector<double> rhs = {0.0, 0.0};
    std::vector<double> solution;

    REQUIRE_NOTHROW(solver.solve(rhs, solution));
    REQUIRE(solution.size() == 2);
    REQUIRE(solution[0] == Approx(0.0));
    REQUIRE(solution[1] == Approx(0.0));
}
