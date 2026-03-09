#include "main_test.hpp"
#include "poisson_problem_solver/schwarz_methods/subdomain.hpp"
#include <unordered_map>
#include <vector>

double s_f(double x, double y) {
    return x - x + y - y + 1;
}

double b_f(double x, double y) {
    return x - x + y - y - 1;
}


TEST_CASE("Subdomain::create_matrix()") {
    // 0 0 * * 0
    // 0 * * * *
    // * * * * *
    // 0 * * * 0
    // 0 0 * 0 0
    size_t N = 5;
    double h = 1.0 / N;
    double b = -1.0 / (h * h);
    double t = -1.0 / (h * h);
    double l = -1.0 / (h * h);
    double r = -1.0 / (h * h);
    double c = 4.0 / (h * h);

    Subdomain subdomain(N);
    std::vector<size_t> indices = {7, 11, 12, 13, 17, 18, 2, 6, 8, 10, 14, 16, 19, 22, 23};
    for (auto i : indices)
        subdomain.add_point(i);

    subdomain.create_matrix(h, s_f, b_f);

    std::vector<double> correct_rhs = {1, 1, 1, 1, 1, 1, -1, 0, 0, -1, -1, 0, -1, -1, -1};
    std::vector<size_t> correct_addr = {0, 5, 10, 15, 20, 25, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39};
    std::vector<size_t> correct_cols = {0, 2, 6, 7,  8,  1, 2, 7, 9,  11, 0, 1, 2, 3, 4,  2,  3,  5,  8, 10,
                                        2, 4, 5, 11, 13, 3, 4, 5, 12, 14, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    std::vector<double> correct_vals = {
        c, t, b, l, r, c, r, b, l, t, b,   l,   c,   r,   t,   l,   c,   t,   b,   r,
        b, c, r, l, t, b, l, c, r, t, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    };
    auto matrix = subdomain.get_matrix();
    CHECK(matrix.vals() == correct_vals);
    CHECK(matrix.addr() == correct_addr);
    CHECK(matrix.cols() == correct_cols);
    CHECK(subdomain.get_rhs() == correct_rhs);
}
TEST_CASE("Subdomain") {
    // 0 0 * * * 0 0 0 0 0 0 0
    // 0 0 * * * 0 0 0 0 0 0 0
    // 0 0 0 * 0 0 0 0 * 0 0 0
    // 0 0 0 0 0 0 * 0 0 0 0 0
    // 0 * 0 0 0 * * * 0 0 0 0
    // * * * 0 0 0 * 0 0 0 0 0
    // * * * 0 0 0 0 0 0 0 0 *
    // * * * 0 0 0 0 0 0 0 * *
    // * 0 0 0 0 0 0 0 0 * * *
    // 0 0 0 0 0 0 0 0 * * * *
    // 0 0 0 0 0 0 0 * * * * *
    // 0 0 0 0 0 0 * * * * * *

    size_t N = 12;
    Subdomain subdomain(N);
    std::vector<size_t> correct_indices = {6,  7,  8,  9,  10, 11, 19,  20,  21,  22,  23,  32,  33,  34,  35,
                                           36, 45, 46, 47, 48, 49, 50,  58,  59,  60,  61,  62,  71,  72,  73,
                                           74, 78, 85, 89, 90, 91, 102, 111, 116, 122, 123, 124, 134, 135, 136};
    std::vector<size_t> correct_neighbors = {5,   18,  24,  31,  37,  38,  44,  51,  57,  63,  66,  70,
                                             75,  77,  79,  83,  84,  86,  88,  92,  97,  99,  101, 103,
                                             104, 110, 112, 114, 115, 117, 121, 125, 128, 133, 137};
    std::vector<size_t> correct_boundary = {6,  7,  8,  9,   10,  11,  19,  23,  32,  35,  36, 45,
                                            47, 48, 49, 50,  58,  59,  60,  62,  71,  72,  74, 78,
                                            85, 89, 91, 102, 111, 116, 122, 124, 134, 135, 136};

    for (auto val : correct_indices)
        subdomain.add_point(val);

    SECTION("Subdomain::contains()") {
        std::vector<size_t> all_subdomain_points;
        for (size_t point = 0; point < N * N; point++)
            if (subdomain.contains(point))
                all_subdomain_points.push_back(point);
        CHECK(all_subdomain_points == correct_indices);
    }

    SECTION("Subdomain::get_weights()") {
        for (const auto& weight : subdomain.get_weights())
            CHECK(weight == Approx(1.0));
    }

    SECTION("Subdomain::add_point()") {
        std::vector<size_t> indices = subdomain.get_indices();
        std::vector<int> mask = subdomain.get_mask();
        REQUIRE(indices.size() == correct_indices.size());
        CHECK(indices == correct_indices);

        int temp = 0;
        for (size_t k = 0; k < N * N; k++) {
            if (k == correct_indices[temp])
                CHECK(mask[k] == temp++);
            else
                CHECK(mask[k] == -1);
        }
    }

    SECTION("Subdomain::is_neighbor()") {
        int temp = 0;
        for (size_t i = 0; i < N * N; i++) {
            if (i == correct_neighbors[temp]) {
                CHECK(subdomain.is_neighbor(i) == true);
                temp++;
            } else
                CHECK(subdomain.is_neighbor(i) == false);
        }
    }

    SECTION("Subdomain::is_boundary()") {
        int temp = 0;
        for (size_t i = 0; i < N * N; i++) {
            if (i == correct_boundary[temp]) {
                CHECK(subdomain.is_boundary(i) == true);
                temp++;
            } else
                CHECK(subdomain.is_boundary(i) == false);
        }
    }

    SECTION("Subdomain::set_weight()") {
        std::vector<double> weights(correct_indices.size());
        for (size_t i = 0; i < weights.size(); i++)
            weights[i] = static_cast<double>(i);
        for (size_t i = 0; i < correct_indices.size(); i++)
            subdomain.set_weight(correct_indices[i], weights[i]);
        CHECK(weights == subdomain.get_weights());
    }
}
