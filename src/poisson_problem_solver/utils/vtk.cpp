#include "vtk.hpp"
#include <filesystem>
#include <iomanip>
#include <sstream>

void VtkWriter::append_header(const std::string& describe, std::ostream& fs) {
    fs << "# vtk DataFile Version 3.0" << std::endl;
    fs << describe << std::endl;
    fs << "ASCII" << std::endl << std::endl;
}
void VtkWriter::append_points(const std::vector<Point> points, size_t N, std::ostream& fs) {
    fs << "DATASET STRUCTURED_GRID" << std::endl;
    fs << "DIMENSIONS" << ' ' << N << ' ' << N << ' ' << 1 << std::endl;
    fs << "POINTS" << ' ' << points.size() << ' ' << "float" << std::endl;
    for (auto point : points)
        fs << point.x << ' ' << point.y << ' ' << point.z << std::endl;
}
void VtkWriter::add_point_data(const std::vector<double>& data, const std::string& data_cap, std::ostream& fs) {
    fs << "SCALARS " << data_cap << " float" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    for (auto value : data)
        fs << value << std::endl;
}

void VtkWriter::append_point_data_header(size_t data_size, std::ostream& fs) {
    fs << "POINT_DATA" << ' ' << data_size << std::endl;
}

VtkWriter::StepManager::StepManager(const std::string& stem, size_t step) : stem_(stem), step_(step) {
    if (std::filesystem::is_directory(stem))
        std::filesystem::remove_all(stem);
    std::filesystem::create_directory(stem);
}

std::string VtkWriter::StepManager::add(size_t iter, bool force) {
    if (!force && iter % this->step_ != 0) {
        return "";
    }

    std::ostringstream fn;
    fn << stem_ << '_' << std::setfill('0') << std::setw(8) << iter << ".vtk";
    std::string ret = stem_ + '/' + fn.str();
    return ret;
}
