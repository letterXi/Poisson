#ifndef VTK_HPP
#define VTK_HPP

#include <fstream>
#include <string>
#include <vector>

namespace VtkWriter {

struct Point {
    double x;
    double y;
    double z;
};

void append_header(const std::string& describe, std::ostream& fs);
void append_points(const std::vector<Point> points, size_t N, std::ostream& fs);
void add_point_data(const std::vector<double>& data, const std::string& data_cap, std::ostream& fs);
void append_point_data_header(size_t data_size, std::ostream& fs);

class StepManager {
public:
    StepManager(const std::string& stem, size_t step);
    std::string add(size_t step, bool force);

private:
    std::string stem_;
    size_t step_;
};
} // namespace VtkWriter

#endif
