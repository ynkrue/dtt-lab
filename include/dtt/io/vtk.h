#pragma once

#include "dtt/sim/forces.h"
#include "dtt/sim/particles.h"

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace dtt::io {

// Write particles (and optional forces) as legacy ASCII VTK PolyData.
// Points are (x, y, 0). Outputs VERTICES so ParaView treats them as particles.
bool write_legacy_vtk(const sim::Particles& particles,
                      std::string_view file_path,
                      const sim::ForceField* forces = nullptr);

// Write particles as XML .vtp PolyData so they can be used in a .pvd collection.
bool write_vtp_polydata(const sim::Particles& particles,
                        std::string_view file_path,
                        const sim::ForceField* forces = nullptr);

// Write a .pvd collection referencing per-step VTK files.
// time_file_pairs: (time, relative file name).
bool write_pvd_collection(const std::string& pvd_path,
                          const std::vector<std::pair<double, std::string>>& time_file_pairs);

}  // namespace dtt::io
