#include "poisson_problem/vtk/vtk_save.hpp"
#include <fstream>

vtkWriter::vtkWriter(const std::string& name, const std::string& describe, const Grid& grid)
    : file_name_(name), file_describe_(describe), grid_(grid) {}
void vtkWriter::write() const {
    std::ofstream vtkfile(file_name_ + ".vtk");
    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << file_describe_ << std::endl;
    vtkfile << "ASCII" << std::endl << std::endl;
    vtkfile << "DATASET STRUCTURED_GRID" << std::endl;
    vtkfile << "DIMENSIONS" << ' ' << grid_.get_N_x() << ' ' << grid_.get_N_y() << ' ' << 1 << std::endl;
    vtkfile << "POINTS" << ' ' << grid_.points() << ' ' << "float" << std::endl;
    for (size_t j = 0; j < grid_.get_N_y(); j++) {
        for (size_t i = 0; i < grid_.get_N_x(); i++) {
            vtkfile << grid_.getX(i) << ' ' << grid_.getY(j) << ' ' << 0.0 << std::endl;
        }
    }
    vtkfile << std::endl << "POINT_DATA" << ' ' << grid_.points() << std::endl;
    for (auto scalar : scalars_) {
        vtkfile << "SCALARS " << scalar.first << " float" << std::endl;
        vtkfile << "LOOKUP_TABLE default" << std::endl;
        for (auto value : scalar.second)
            vtkfile << value << std::endl;
    }

    vtkfile.close();
}

void vtkWriter::add_scalars(std::vector<double> a, const std::string& scalar_name) {
    scalars_.push_back(std::make_pair(scalar_name, a));
}
