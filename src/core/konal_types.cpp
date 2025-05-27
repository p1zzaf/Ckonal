#include "konal_types.hpp"
#include <stdexcept>
#include <algorithm> // For std::transform

namespace Ckonal {

std::string to_string(CoordSys cs) {
    switch (cs) {
        case CoordSys::CARTESIAN: return "cartesian";
        case CoordSys::SPHERICAL: return "spherical";
        default: return "unknown";
    }
}

CoordSys coord_sys_from_string(const std::string& s_in) {
    std::string s = s_in;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower); // Convert to lowercase
    if (s == "cartesian") return CoordSys::CARTESIAN;
    if (s == "spherical") return CoordSys::SPHERICAL;
    throw std::invalid_argument("Unknown coordinate system string: " + s_in);
}

} // namespace Ckonal