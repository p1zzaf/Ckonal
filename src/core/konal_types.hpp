#ifndef KONAL_TYPES_HPP
#define KONAL_TYPES_HPP

#include "konal_constants.hpp"
#include <array>
#include <vector>
#include <string>
#include <functional> // For std::hash

namespace Ckonal {

// 3D Index for grids
struct Index3D {
    size_t i1 = 0, i2 = 0, i3 = 0;

    Index3D() = default;
    Index3D(size_t r, size_t c, size_t d) : i1(r), i2(c), i3(d) {}

    bool operator==(const Index3D& other) const {
        return i1 == other.i1 && i2 == other.i2 && i3 == other.i3;
    }
    bool operator!=(const Index3D& other) const {
        return !(*this == other);
    }
    // Add < operator for use in std::map or std::set if needed
    bool operator<(const Index3D& other) const {
        if (i1 != other.i1) return i1 < other.i1;
        if (i2 != other.i2) return i2 < other.i2;
        return i3 < other.i3;
    }
};

// Coordinate system type
enum class CoordSys {
    CARTESIAN,
    SPHERICAL
    // GEOGRAPHICAL could be added if direct support is needed beyond transformations
};

std::string to_string(CoordSys cs);
CoordSys coord_sys_from_string(const std::string& s);


// 3D Point/Vector
using Point3D = std::array<real_t, 3>;

// For multi-dimensional data (e.g., field values)
// Option 1: Nested vectors (simpler, potentially less performant for large data)
// template<typename T>
// using GridData3D = std::vector<std::vector<std::vector<T>>>;
// template<typename T>
// using GridData4D = std::vector<std::vector<std::vector<std::vector<T>>>>;

// Option 2: Flattened vector with manual indexing (more performant)
// This will be managed within the Field classes.

} // namespace Ckonal

// Hash function for Index3D for use in std::unordered_map/set
namespace std {
    template <> struct hash<Ckonal::Index3D> {
        size_t operator()(const Ckonal::Index3D& idx) const {
            size_t h1 = std::hash<size_t>{}(idx.i1);
            size_t h2 = std::hash<size_t>{}(idx.i2);
            size_t h3 = std::hash<size_t>{}(idx.i3);
            // A common way to combine hashes
            return h1 ^ (h2 << 1) ^ (h3 << 2); // Or boost::hash_combine
        }
    };
}


#endif // KONAL_TYPES_HPP