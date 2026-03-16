#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include <cmath>


inline double exactFunction(double x, double y) {
    return 0.2 * (std::sin(2 * M_PI * x) - std::cos(12 * M_PI * y) + std::tan(x * y) + std::exp(-x * y));
}

inline double dirichletBoundaryFunction(double x, double y) {
    return exactFunction(x, y);
}

inline double sourceFunction(double x, double y) {
    return -0.2 * (-4 * M_PI * M_PI * std::sin(2 * M_PI * x) + 144 * M_PI * M_PI * std::cos(12 * M_PI * y) +
                   2 * (x * x + y * y) * (1.0 / std::cos(x * y)) * (1.0 / std::cos(x * y)) * std::tan(x * y) +
                   (x * x + y * y) * std::exp(-x * y));
}

#endif
