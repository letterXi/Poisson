#ifndef VTK_SAVER_HPP
#define VTK_SAVER_HPP
#include "poisson_problem/grid/grid.hpp"
#include <string>
#include <utility>
#include <vector>

class vtkWriter {
public:
    vtkWriter(const std::string& name, const std::string& describe, const Grid& grid);
    void add_scalars(std::vector<double> a, const std::string& scalar_name);
    void write() const;

private:
    std::string file_name_;
    std::string file_describe_;
    Grid grid_;
    std::vector<std::pair<std::string, std::vector<double>>> scalars_;
};
#endif
