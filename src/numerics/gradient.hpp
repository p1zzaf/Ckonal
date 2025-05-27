#ifndef Ckonal_GRADIENT_HPP
#define Ckonal_GRADIENT_HPP

#include "../core/konal_types.hpp"
#include <vector>
#include <memory> // For std::unique_ptr

// Forward declarations
namespace Ckonal {
class ScalarField3D; // Gradient is calculated for a scalar field
class VectorField3D; // Gradient result is a vector field
}

namespace Ckonal {
namespace gradient {

// --- Gradient for Cartesian Coordinates ---
// Calculates gradient and stores it in the pre-allocated result_values vector.
// result_values should be sized N0*N1*N2*3.
void calculate_cartesian(
    const ScalarField3D& scalar_field,
    std::vector<real_t>& result_values // Output: flattened vector components
);

// --- Gradient for Spherical Coordinates ---
void calculate_spherical(
    const ScalarField3D& scalar_field,
    std::vector<real_t>& result_values // Output: flattened vector components [gr, gtheta, gphi, ...]
);

} // namespace gradient
} // namespace Ckonal

#endif // Ckonal_GRADIENT_HPP