#ifndef Ckonal_INTERPOLATION_HPP
#define Ckonal_INTERPOLATION_HPP

#include "../core/konal_types.hpp"
#include "../fields/field_base.hpp" // Need FieldBase for grid info
#include <vector>

namespace Ckonal {
// Forward declare field types to avoid circular dependencies if they call interpolation directly.
// However, it's cleaner if interpolation functions are general and take FieldBase + data.
// class ScalarField3D;
// class VectorField3D;

namespace interpolation {

// --- Trilinear Interpolation for Scalar Fields ---
real_t trilinear_scalar(
    const FieldBase& field_grid,        // Provides grid structure (min_coords, intervals, npts, periodicity)
    const std::vector<real_t>& values,  // Flattened scalar values N0*N1*N2
    const Point3D& point,               // Point to interpolate at (in field's coord_sys)
    real_t null_value = NaN
);

// --- Trilinear Interpolation for Vector Fields ---
Point3D trilinear_vector(
    const FieldBase& field_grid,        // Provides grid structure
    const std::vector<real_t>& values,  // Flattened vector values N0*N1*N2*3
    const Point3D& point,               // Point to interpolate at
    const Point3D& null_value = {NaN, NaN, NaN}
);

} // namespace interpolation
} // namespace Ckonal

#endif // Ckonal_INTERPOLATION_HPP