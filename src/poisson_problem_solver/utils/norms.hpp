#ifndef NORMS_HPP
#define NORMS_HPP

#include <vector>
#include <cmath>     
#include <algorithm> 

template <typename T>
inline std::vector<T> diff_of(const std::vector<T>& v1, const std::vector<T>& v2)
{
    std::vector<T> res;
    res.reserve(v1.size());
    for(size_t i = 0; i < v1.size(); i++)
    {
        res.push_back(v1[i] - v2[i]);
    }
    return res;
}

inline double norm_inf(const std::vector<double>& v)
{
    double max_v = 0.0;
    for(double x : v)
        max_v = std::max(std::abs(x), max_v);
    return max_v;
}

inline double norm_L2(const std::vector<double>& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return std::sqrt(sum / v.size());
}

#endif
