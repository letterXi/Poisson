#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "subdomain.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class SchwarzSolver {
public:
    SchwarzSolver(size_t N, std::vector<size_t> mask, std::function<double(double, double)> source_function,
                  std::function<double(double, double)> boundary_function, size_t maxiter = 1000, double tolerance = 1e-6);

    virtual ~SchwarzSolver() = default;

    const std::vector<std::unique_ptr<Subdomain>>& get_subdomains() const;

    const std::unordered_map<size_t, size_t>& get_id_to_idx() const;

    std::vector<size_t> find_unique_indices() const;

    void set_overlap(size_t overlap);

    virtual void iterate(std::vector<double>& u) = 0;

    void solve(std::vector<double>& u, size_t& iters);

    void create_slaes();

    void connect_solves(std::vector<double>& u);

    void initialize(const std::vector<double>& global_u);

    virtual std::string get_name() const = 0;

protected:
    size_t N_;
    std::vector<size_t> mask_;
    std::vector<std::unique_ptr<Subdomain>> subdomains_;
    std::unordered_map<size_t, size_t> id_to_idx_;
    std::function<double(double, double)> source_function_;
    std::function<double(double, double)> boundary_function_;
    size_t maxiter_;
    double tolerance_;
};

#endif
