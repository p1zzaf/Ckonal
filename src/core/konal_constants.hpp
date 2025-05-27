#ifndef KONAL_CONSTANTS_HPP
#define KONAL_CONSTANTS_HPP

#include <cstddef> // For size_t
#include <cstdint> // For fixed-width integers
#include <limits>  // For std::numeric_limits

namespace Ckonal {

// Basic data types (mirroring DTYPE_REAL, DTYPE_UINT, etc.)
using real_t = double;
using uint_t = std::uint32_t;
using int_t = std::int32_t;
using bool_t = bool; // Or potentially char/std::byte for memory if packing bools

// Mathematical Constants
constexpr real_t PI = 3.14159265358979323846;
constexpr real_t TWO_PI = 2.0 * PI;
constexpr real_t INF = std::numeric_limits<real_t>::infinity();
constexpr real_t NEG_INF = -std::numeric_limits<real_t>::infinity();
constexpr real_t NaN = std::numeric_limits<real_t>::quiet_NaN();
constexpr real_t EPSILON = 1e-6;

// Physical Constants
constexpr real_t EARTH_RADIUS_KM = 6371.0; // kilometers

// Default fill values or special markers
constexpr int_t HEAP_INVALID_INDEX = -1;

} // namespace pykonal

#endif // KONAL_CONSTANTS_HPP