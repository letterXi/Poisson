#ifndef VTK_SAVER_HPP
#define VTK_SAVER_HPP
#include "shwartz.hpp"
#include <string>
#include <utility>
#include <vector>

class vtkWriter {
public:
    vtkWriter(const std::string& name, const std::string& describe, const Mesh& mesh);
    void add_scalars(std::vector<double> a, const std::string& scalar_name);
    void write() const;

private:
    std::string file_name_;
    std::string file_describe_;
    Mesh mesh_;
    std::vector<std::pair<std::string, std::vector<double>>> scalars_;
};
#endif
