#include "poisson_problem/vector_operations/vector_operations.hpp"

std::vector<double> operator-(const std::vector<double>& u1, const std::vector<double>& u2) {
    std::vector<double> res(u1.size());
    for (size_t i = 0; i < u1.size(); i++) {
        res[i] = u1[i] - u2[i];
    }
    return res;
}