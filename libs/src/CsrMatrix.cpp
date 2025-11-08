#include "CsrMatrix.hpp"
#include <algorithm>

const std::vector<size_t>& CsrMatrix::addr() const {
    return addr_;
}

const std::vector<size_t>& CsrMatrix::cols() const {
    return cols_;
}

std::vector<double>& CsrMatrix::vals() {
    return vals_;
}

const std::vector<double>& CsrMatrix::vals() const {
    return vals_;
}

void CsrMatrix::set_values(std::vector<double>&& vals) {
    vals_ = std::move(vals);
}

size_t CsrMatrix::n_nonzeros() const {
    return cols_.size();
}

size_t CsrMatrix::n_rows() const {
    return addr_.size() - 1;
}

double CsrMatrix::find_value(size_t i, size_t j) const {
    std::vector<size_t>::const_iterator it_begin = cols_.begin() + addr_[i];
    std::vector<size_t>::const_iterator it_end = cols_.begin() + addr_[i + 1];

    std::vector<size_t>::const_iterator fnd = std::lower_bound(it_begin, it_end, j);
    if (fnd != it_end && *fnd == j) {
        size_t a = fnd - cols_.begin();
        return vals_[a];
    } else {
        return 0;
    }
}
