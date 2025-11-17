#ifndef CSR_MATRIX_HPP
#define CSR_MATRIX_HPP

#include <vector>

class CsrMatrix {
public:
    CsrMatrix() = default;

    CsrMatrix(std::vector<size_t>&& addr, std::vector<size_t>&& cols, std::vector<double>&& vals)
        : addr_(std::move(addr)), cols_(std::move(cols)), vals_(std::move(vals)) {};

    const std::vector<size_t>& addr() const;

    const std::vector<size_t>& cols() const;

    std::vector<double>& vals();

    const std::vector<double>& vals() const;

    double find_value(size_t i, size_t j) const;

    void set_values(std::vector<double>&& vals);

    size_t n_nonzeros() const;

    size_t n_rows() const;

private:
    std::vector<size_t> addr_;
    std::vector<size_t> cols_;
    std::vector<double> vals_;
};

#endif
