#ifndef NORMS_HPP
#define NORMS_HPP
#include <vector>
double inf_norm(const std::vector<double>& u);
double l2_norm(const std::vector<double>& u);
double rms_norm(const std::vector<double>& u);
#endif