#include "poisson_problem/norms/norms.hpp"
#include <cmath>

double inf_norm(const std::vector<double>& u) {
    double temp = 0.0;
    for (size_t i = 0; i < u.size(); i++) {
        if (std::abs(u[i] > temp))
            temp = std::abs(u[i]);
    }
    return temp;
}