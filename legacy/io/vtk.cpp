#include "dtt/io/vtk.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <string_view>
#include <vector>

namespace dtt::io {

bool write_vtp_polydata(const sim::Particles &particles, std::string_view file_path,
                        const sim::ForceField *forces) {
    std::ofstream out{std::string(file_path)};
    if (!out.is_open()) {
        return false;
    }

    const std::size_t n = particles.x.size();
    out << R"(<?xml version="1.0"?>)" << "\n";
    out << R"(<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">)" << "\n";
    out << "<PolyData>\n";
    out << "  <Piece NumberOfPoints=\"" << n << "\" NumberOfVerts=\"" << n
        << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    out << "    <Points>\n";
    out << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << std::setprecision(12);
    for (std::size_t i = 0; i < n; ++i) {
        out << "        " << particles.x[i] << " " << particles.y[i] << " 0.0\n";
    }
    out << "      </DataArray>\n";
    out << "    </Points>\n";

    out << "    <Verts>\n";
    out << "      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) {
        out << "        " << static_cast<int>(i) << "\n";
    }
    out << "      </DataArray>\n";
    out << "      <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) {
        out << "        " << static_cast<int>(i + 1) << "\n";
    }
    out << "      </DataArray>\n";
    out << "    </Verts>\n";

    out << "    <PointData>\n";
    out << "      <DataArray type=\"Float64\" Name=\"mass\" NumberOfComponents=\"1\" "
           "format=\"ascii\">\n";
    for (double m : particles.mass) {
        out << "        " << m << "\n";
    }
    out << "      </DataArray>\n";

    out << "      <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" "
           "format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) {
        out << "        " << particles.vx[i] << " " << particles.vy[i] << " 0.0\n";
    }
    out << "      </DataArray>\n";

    if (forces) {
        out << "      <DataArray type=\"Float64\" Name=\"force\" NumberOfComponents=\"3\" "
               "format=\"ascii\">\n";
        for (const auto &f : *forces) {
            out << "        " << f[0] << " " << f[1] << " 0.0\n";
        }
        out << "      </DataArray>\n";
    }

    out << "    </PointData>\n";
    out << "  </Piece>\n";
    out << "</PolyData>\n";
    out << "</VTKFile>\n";

    return true;
}

bool write_pvd_collection(const std::string &pvd_path,
                          const std::vector<std::pair<double, std::string>> &time_file_pairs) {
    std::ofstream out(pvd_path);
    if (!out.is_open()) {
        return false;
    }

    out << R"(<?xml version="1.0"?>)" << "\n";
    out << R"(<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">)" << "\n";
    out << "<Collection>\n";
    for (const auto &[time, file] : time_file_pairs) {
        out << "  <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << file
            << "\"/>\n";
    }
    out << "</Collection>\n";
    out << "</VTKFile>\n";
    return true;
}

} // namespace dtt::io
