#define CATCH_CONFIG_MAIN
#include "poisson_problem/mat/csr_matrix.hpp"
#include <catch2/catch_all.hpp>
#include <vector>
using Catch::Approx;
TEST_CASE("CsrMatrix with common matrix") {
    // 0 3 0
    // 1 1 6
    // 2 0 0
    std::vector<double> vals = {3, 1, 1, 6, 2};
    std::vector<size_t> cols = {1, 0, 1, 2, 0};
    std::vector<size_t> addr = {0, 1, 4, 5};

    CsrMatrix mat(std::move(addr), std::move(cols), std::move(vals));

    SECTION("find nonzeros values", "[find][nonzeros]") {
        REQUIRE(mat.find_value(0, 1) == Approx(3));
        REQUIRE(mat.find_value(1, 0) == Approx(1));
        REQUIRE(mat.find_value(1, 1) == Approx(1));
        REQUIRE(mat.find_value(1, 2) == Approx(6));
        REQUIRE(mat.find_value(2, 0) == Approx(2));
    }

    SECTION("find zeros", "[find][zeros]") {
        REQUIRE(mat.find_value(0, 0) == Approx(0));
        REQUIRE(mat.find_value(0, 2) == Approx(0));
        REQUIRE(mat.find_value(2, 1) == Approx(0));
        REQUIRE(mat.find_value(2, 2) == Approx(0));
    }

    SECTION("use basic methods", "[nrows][nonzeros]") {
        REQUIRE(mat.n_nonzeros() == 5);
        REQUIRE(mat.n_rows() == 3);
    }

    SECTION("set new values", "[set_values][find]") {
        const std::vector<size_t> old_cols = mat.cols();
        const std::vector<size_t> old_addr = mat.addr();

        std::vector<double> new_vals = {5, 23, -1, 45, 30};
        mat.set_values(std::move(new_vals));
        REQUIRE(mat.find_value(0, 1) == Approx(5));
        REQUIRE(mat.find_value(1, 0) == Approx(23));
        REQUIRE(mat.find_value(1, 1) == Approx(-1));
        REQUIRE(mat.find_value(1, 2) == Approx(45));
        REQUIRE(mat.find_value(2, 0) == Approx(30));

        SECTION("const stencil") {
            REQUIRE(mat.cols() == old_cols);
            REQUIRE(mat.addr() == old_addr);
        }
    }
}

TEST_CASE("CsrMatrix with empty matrix") {
    // 0 0
    // 0 0

    std::vector<double> vals = {};
    std::vector<size_t> cols = {};
    std::vector<size_t> addr = {0, 0, 0};

    CsrMatrix mat(std::move(addr), std::move(cols), std::move(vals));

    REQUIRE(mat.n_rows() == 2);
    REQUIRE(mat.n_nonzeros() == 0);

    REQUIRE(mat.find_value(0, 0) == Approx(0));
    REQUIRE(mat.find_value(0, 1) == Approx(0));
    REQUIRE(mat.find_value(1, 0) == Approx(0));
    REQUIRE(mat.find_value(1, 1) == Approx(0));
}
