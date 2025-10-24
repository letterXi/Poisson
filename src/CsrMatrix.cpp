#include "CsrMatrix.hpp"

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