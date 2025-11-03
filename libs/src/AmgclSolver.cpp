#include "AmgclSolver.hpp"
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <boost/property_tree/ptree.hpp>
#include <vector>

class AmgclSolver::Impl {
public:
    using param_t = boost::property_tree::ptree;
    using matrix_t = amgcl::backend::crs<double, size_t>;
    using backend_t = amgcl::backend::builtin<double>;
    using solver_t = amgcl::make_solver<
        amgcl::amg<backend_t, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>,
        amgcl::runtime::solver::wrapper<backend_t>>;

    Impl(matrix_t mat, param_t param)
        : solver_(mat, param), dim_((size_t)mat.nrows), maxit_(param.get("solver_maxiter", -1)) {}

    void solve(const std::vector<double>& rhs, std::vector<double>& ret) const {
        ret.resize(dim_, 0.0);
        int n_iters;
        double err;

        std::tie(n_iters, err) = solver_(rhs, ret);
    }

private:
    solver_t solver_;
    size_t dim_;
    size_t maxit_;
};

void AmgclSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const {
    pimpl_->solve(rhs, ret);
}

AmgclSolver::AmgclSolver(std::initializer_list<std::pair<std::string, std::string>> amgc_params) {
    for (auto it : amgc_params) {
        params_[it.first] = it.second;
    }
}

AmgclSolver::AmgclSolver(size_t maxit, double tolerance) : maxit_(maxit), tolerance_(tolerance) {}

AmgclSolver::~AmgclSolver() = default;

void AmgclSolver::set_matrix(const CsrMatrix& mat) {
    Impl::matrix_t amgcl_matrix;
    amgcl_matrix.own_data = false;
    amgcl_matrix.ncols = amgcl_matrix.nrows = mat.n_rows();
    amgcl_matrix.nnz = mat.n_nonzeros();
    amgcl_matrix.ptr = const_cast<size_t*>(mat.addr().data());
    amgcl_matrix.col = const_cast<size_t*>(mat.cols().data());

    const std::vector<double>& refmat = mat.vals();
    amgcl_matrix.val = const_cast<double*>(refmat.data());

    Impl::param_t prm;
    prm.put("solver.type", "fgmres");
    prm.put("solver.tol", tolerance_);
    prm.put("solver.maxiter", maxit_);
    prm.put("precond.coarsening.type", "smoothed_aggregation");
    prm.put("precond.relax.type", "spai0");

    for (auto it : params_) {
        prm.put(it.first, it.second);
    }

    pimpl_.reset(new Impl(amgcl_matrix, prm));
}