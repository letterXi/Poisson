#include "subdomain.hpp"
#include <algorithm>

Subdomain::Subdomain(size_t N) : N_(N) {}

void Subdomain::add_point(size_t point, double weight) {
    mask_[point] = indices_.size();
    indices_.push_back(point);
    weights_.push_back(weight);
}

bool Subdomain::is_neighbor(size_t point) const {
    // Если точка уже в подобласти — она не сосед, она "своя"
    if (this->contains(point)) return false;

    size_t j = point / N_;
    size_t i = point % N_;

    // Проверяем 8 соседей вокруг точки
    for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
            if (di == 0 && dj == 0) continue;
            
            size_t ni = i + di;
            size_t nj = j + dj;
            
            if (ni < N_ && nj < N_) {
                size_t neighbor_idx = nj * N_ + ni;
                if (this->contains(neighbor_idx)) return true;
            }
        }
    }
    return false;
}

bool Subdomain::is_boundary(size_t point) const {
    // 1. Если точки вообще нет в этой подобласти — это не её граница
    auto it = mask_.find(point);
    if (it == mask_.end()) return false;

    size_t j = point / N_;
    size_t i = point % N_;

    // 2. Проверка физических границ всей сетки (Дирихле)
    if (i == 0 || i == N_ - 1 || j == 0 || j == N_ - 1)
        return true;

    // 3. Проверка внутренних границ (интерфейсов)
    // Если хотя бы одного соседа нет в нашей mask_map_ — значит точка граничная
    if (mask_.find(point + 1) == mask_.end() || 
        mask_.find(point - 1) == mask_.end() || 
        mask_.find(point + N_) == mask_.end() || 
        mask_.find(point - N_) == mask_.end()) 
    {
        return true;
    }

    return false;
}

const std::vector<size_t>& Subdomain::get_indices() const {
    return indices_;
}
const std::vector<double>& Subdomain::get_weights() const {
    return weights_;
}

bool Subdomain::contains(size_t point) const {
    return mask_.find(point) != mask_.end();
}

void Subdomain::set_weight(size_t point, double weight) {
    auto it = mask_.find(point);
    if (it != mask_.end()) {
        weights_[it->second] = weight; // it->second — это и есть локальный индекс
    }
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
        size_t global_idx = indices_[i];
        size_t I = global_idx % N_;
        size_t J = global_idx / N_;
        double x = static_cast<double>(I) * h;
        double y = static_cast<double>(J) * h;

        if (this->is_boundary(global_idx)) {
            cols.push_back(i);
            vals.push_back(1.0);
            
            if (I == 0 || J == 0 || I == N_ - 1 || J == N_ - 1) {
                rhs_[i] = boundary_func(x, y);
            }
        } 
        else {
            rhs_[i] = source_func(x, y);

            // Используем .at(), так как для внутренних точек соседи точно есть в маске
            std::vector<size_t> row_indices = {
                i,                          
                mask_.at(global_idx - N_), 
                mask_.at(global_idx + N_), 
                mask_.at(global_idx - 1),  
                mask_.at(global_idx + 1)   
            };

            std::sort(row_indices.begin(), row_indices.end());

            for (size_t col_idx : row_indices) {
                cols.push_back(col_idx);
                vals.push_back((col_idx == i) ? coeff_2 : coeff_1);
            }
        }
        addr.push_back(cols.size());
    }

    matrix_ = CsrMatrix(std::move(addr), std::move(cols), std::move(vals));
    solver_.set_matrix(matrix_);
    u_ = std::vector<double>(indices_.size(), 0.0);
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
