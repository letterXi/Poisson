#define CATCH_CONFIG_MAIN
#include "AmgclSolver.hpp"
#include "CsrMatrix.hpp"
#include "poisson.hpp"
#include <catch2/catch_all.hpp>
#include <vector>
using Catch::Approx;

TEST_CASE("True fill matrix", "[Poisson]") {
    const int n = 10;

    CsrMatrix matrix_0, matrix_1, matrix_2, matrix_3;

    std::vector<size_t> addr_0, cols_0, addr_1, cols_1, addr_2, cols_2, addr_3, cols_3;
    std::vector<double> vals_0, vals_1, vals_2, vals_3;
    std::vector<double> rhs_0, rhs_1, rhs_2, rhs_3;

    Poisson(n, addr_0, cols_0, vals_0, rhs_0, 0);
    Poisson(n, addr_1, cols_1, vals_1, rhs_1, 1);
    Poisson(n, addr_2, cols_2, vals_2, rhs_2, 2);
    Poisson(n, addr_3, cols_3, vals_3, rhs_3, 3);

    matrix_0 = CsrMatrix(std::move(addr_0), std::move(cols_0), std::move(vals_0));
    matrix_1 = CsrMatrix(std::move(addr_1), std::move(cols_1), std::move(vals_1));
    matrix_2 = CsrMatrix(std::move(addr_2), std::move(cols_2), std::move(vals_2));
    matrix_3 = CsrMatrix(std::move(addr_3), std::move(cols_3), std::move(vals_3));

    REQUIRE(matrix_0.addr() == matrix_1.addr());
    REQUIRE(matrix_0.addr() == matrix_2.addr());
    REQUIRE(matrix_0.addr() == matrix_3.addr());

    REQUIRE(matrix_0.cols() == matrix_1.cols());
    REQUIRE(matrix_0.cols() == matrix_2.cols());
    REQUIRE(matrix_0.cols() == matrix_3.cols());

    REQUIRE(matrix_0.vals() == matrix_1.vals());
    REQUIRE(matrix_0.vals() == matrix_2.vals());
    REQUIRE(matrix_0.vals() == matrix_3.vals());

    REQUIRE(rhs_0 == rhs_1);
    REQUIRE(rhs_1 == rhs_2);
    REQUIRE(rhs_2 == rhs_3);
}