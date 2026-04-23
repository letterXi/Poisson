#include "subdomain.hpp"
#include <algorithm>

Subdomain::Subdomain(size_t N) : N_(N), mask_(N * N, -1) {}

void Subdomain::add_point(size_t point, double weight) {
    mask_[point] = indices_.size();
    indices_.push_back(point);
    weights_.push_back(weight);
}

bool Subdomain::is_neighbor(size_t point) const {
    if (mask_[point] != -1)
        return false;

    size_t j = point / N_;
    size_t i = point % N_;

    // Проверка 4 основных направлений (Крест)
    if (i > 0 && mask_[point - 1] != -1) return true;           // Лево
    if (i < N_ - 1 && mask_[point + 1] != -1) return true;     // Право
    if (j > 0 && mask_[point - N_] != -1) return true;         // Низ
    if (j < N_ - 1 && mask_[point + N_] != -1) return true;    // Верх

    // ДОБАВЛЯЕМ 4 ДИАГОНАЛИ (чтобы не было усечения под 45°)
    if (i > 0 && j > 0 && mask_[point - N_ - 1] != -1) return true;         // Низ-Лево
    if (i < N_ - 1 && j > 0 && mask_[point - N_ + 1] != -1) return true;     // Низ-Право
    if (i > 0 && j < N_ - 1 && mask_[point + N_ - 1] != -1) return true;     // Верх-Лево
    if (i < N_ - 1 && j < N_ - 1 && mask_[point + N_ + 1] != -1) return true; // Верх-Право

    return false;
}

bool Subdomain::is_boundary(size_t point) const {
    if (mask_[point] == -1)
        return false;

    size_t j = point / N_;
    size_t i = point % N_;

    if (i == 0 || i == N_ - 1 || j == 0 || j == N_ - 1)
        return true;

    if (mask_[point + 1] == -1 || mask_[point - 1] == -1 || mask_[point + N_] == -1 || mask_[point - N_] == -1)
        return true;

    return false;
}

const std::vector<size_t>& Subdomain::get_indices() const {
    return indices_;
}
const std::vector<double>& Subdomain::get_weights() const {
    return weights_;
}
const std::vector<int>& Subdomain::get_mask() const {
    return mask_;
}

bool Subdomain::contains(size_t point) const {
    if (mask_[point] == -1)
        return false;
    return true;
}

void Subdomain::set_weight(size_t point, double weight) {
    weights_[mask_[point]] = weight;
}

const CsrMatrix& Subdomain::get_matrix() const {
    return matrix_;
}

void Subdomain::create_matrix(double h, std::function<double(double, double)> source_func,
                              std::function<double(double, double)> boundary_func) {
    std::vector<size_t> addr = {0};
    std::vector<size_t> cols;
    std::vector<double> vals;

    cols.reserve(indices_.size() * 5);
    vals.reserve(indices_.size() * 5);
    addr.reserve(indices_.size() + 1);

    rhs_ = std::vector<double>(indices_.size(), 0.0);

    double coeff_1 = -1.0 / (h * h);
    double coeff_2 = 4.0 / (h * h);
    for (size_t i = 0; i < indices_.size(); i++) {
        size_t I = indices_[i] % N_;
        size_t J = indices_[i] / N_;
        double x = static_cast<double>(I) * h;
        double y = static_cast<double>(J) * h;

        if (this->is_boundary(indices_[i])) {
            cols.push_back(i);
            vals.push_back(1.0);
            if (I == 0 || J == 0 || I == N_ - 1 || J == N_ - 1)
                rhs_[i] = boundary_func(x, y);
        } else {
            rhs_[i] = source_func(x, y);
            std::vector<size_t> row_indices = {
                i, static_cast<size_t>(mask_[indices_[i] - N_]), static_cast<size_t>(mask_[indices_[i] + N_]),
                static_cast<size_t>(mask_[indices_[i] - 1]), static_cast<size_t>(mask_[indices_[i] + 1])};

            std::sort(row_indices.begin(), row_indices.end());

            for (size_t col_idx : row_indices) {
                cols.push_back(col_idx);
                if (col_idx == i) {
                    vals.push_back(coeff_2);
                } else {
                    vals.push_back(coeff_1);
                }
            }
        }
        addr.push_back(cols.size());
    }
    matrix_ = CsrMatrix(std::move(addr), std::move(cols), std::move(vals));
    solver_.set_matrix(matrix_);
    u_ = std::vector<double>(indices_.size());
}

void Subdomain::set_swap_index(size_t index) {
    swap_index_ = index;
}
const std::vector<double>& Subdomain::get_rhs() const {
    return rhs_;
}

size_t Subdomain::get_swap_index() const {
    return swap_index_;
}
void Subdomain::get_overlap_boundary(const std::vector<double>& u_global) {
    for (size_t i = swap_index_; i < indices_.size(); i++)
        rhs_[i] = u_global[indices_[i]];
}

void Subdomain::solve() {
    solver_.solve(rhs_, u_);
}

const std::vector<double>& Subdomain::get_u() const {
    return u_;
}

void Subdomain::set_u(const std::vector<double>& u) {
    u_ = u;
}
size_t Subdomain::npoints() const {
    return indices_.size();
}
bool Subdomain::is_global_collapse() const
{
  return indices_.size() == N_*N_;
}
