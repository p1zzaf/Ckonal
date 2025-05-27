#include "interpolation.hpp"
#include <cmath> // For std::floor, std::fmod
#include <algorithm> // For std::max, std::min

namespace Ckonal {
namespace interpolation {

real_t trilinear_scalar(
    const FieldBase& field_grid,
    const std::vector<real_t>&_values, // Renamed to avoid conflict with FieldBase::values if it existed
    const Point3D& point,
    real_t null_value)
{
    const auto& min_coords = field_grid.get_min_coords();
    const auto& intervals = field_grid.get_node_intervals();
    const auto& npts = field_grid.get_npts();
    const auto& is_periodic = field_grid.get_axis_is_periodic();
    const auto& is_null_axis = field_grid.get_axis_is_null();

    Point3D normalized_idx_float; // Stores (point_coord - min_coord) / interval for each axis
    std::array<int_t, 3> base_indices;    // Stores floor(normalized_idx_float) for each axis
    Point3D deltas;                // Stores fractional part for each axis

    for (size_t axis = 0; axis < 3; ++axis) {
        if (is_null_axis[axis]) {
            if (npts[axis] != 1 || std::abs(point[axis] - min_coords[axis]) > 1e-6 * intervals[axis]) {
                 // If axis is null, point must be at the single layer, or it's out of simple bounds for null axis.
                 // Ckonal's original value() method used ii[iax][0]=0, ii[iax][1]=0 for null axes.
                 // And delta[iax] = idx[iax] % 1. If idx[iax] is not 0 for null axis, delta != 0.
                 // This implies point[axis] must effectively be min_coords[axis].
                 // Let's be strict: if axis is null, point[axis] should match min_coords[axis].
                 // However, the python code's delta calculation for null axis (idx[axis]%1) implies it can handle
                 // non-zero idx[axis] as long as ii[axis][0/1] are 0. This would make delta[axis] = idx[axis].
                 // This part needs very careful translation of the python logic:
                 // idx[iax]   = (point[iax] - self.cy_min_coords[iax]) / self.cy_node_intervals[iax]
                 // if self.cy_iax_isnull[iax]:
                 //     ii[iax][0] = 0
                 //     ii[iax][1] = 0
                 // delta[iax] = idx[iax] % 1
                 // This means for a null axis, it interpolates between node (0,y,z) and (0,y,z) with weight delta[0].
                 // Which is effectively just using the value at (0,y,z) if delta[0] is an integer.
                 // If delta[0] is not integer, it uses (value + (value - value) * delta[0]), still value.
                 // This implies the result is independent of point[axis] for a null axis IF ii are both 0.
                 // Let's follow:
            }
             // For null axis, normalized_idx_float might not be 0 if point[axis] != min_coords[axis]
            if (intervals[axis] < 1e-9) { // Avoid division by zero for null axis if interval is zero
                normalized_idx_float[axis] = 0.0;
            } else {
                normalized_idx_float[axis] = (point[axis] - min_coords[axis]) / intervals[axis];
            }
            base_indices[axis] = 0; // Index is always 0 for null axis
            deltas[axis] = std::fmod(normalized_idx_float[axis], 1.0);
            if (deltas[axis] < 0.0) deltas[axis] += 1.0; // Ensure delta is in [0,1)
                                                        // Actually, python % 1 gives positive remainder.
                                                        // if idx[iax] can be negative.
                                                        // If idx[iax] is -0.5, idx[iax]%1 = 0.5 in python.
                                                        // C++ fmod(-0.5, 1.0) = -0.5.
                                                        // So, adjustment is needed for negative idx.
        } else {
            // Boundary checks for non-periodic axes
            if (!is_periodic[axis] && (point[axis] < min_coords[axis] - 1e-9 || point[axis] > field_grid.get_max_coords()[axis] + 1e-9)) {
                 // 1e-9 for float tolerance
                return null_value;
            }
            normalized_idx_float[axis] = (point[axis] - min_coords[axis]) / intervals[axis];
            base_indices[axis] = static_cast<int_t>(std::floor(normalized_idx_float[axis]));
            deltas[axis] = normalized_idx_float[axis] - static_cast<real_t>(base_indices[axis]); // Fractional part
        }
    }
    
    // Determine the 8 corner indices (i0,i1), (j0,j1), (k0,k1)
    std::array<std::array<uint_t, 2>, 3> corner_grid_indices; // Stores [axis][0_or_1]
    for (size_t axis = 0; axis < 3; ++axis) {
        if (is_null_axis[axis]) {
            corner_grid_indices[axis][0] = 0;
            corner_grid_indices[axis][1] = 0;
        } else {
            corner_grid_indices[axis][0] = static_cast<uint_t>(base_indices[axis]);
            corner_grid_indices[axis][1] = static_cast<uint_t>(base_indices[axis] + 1);

            // Handle periodicity for indices
            if (is_periodic[axis]) {
                corner_grid_indices[axis][0] = (corner_grid_indices[axis][0] % npts[axis] + npts[axis]) % npts[axis]; // Ensure positive modulo
                corner_grid_indices[axis][1] = (corner_grid_indices[axis][1] % npts[axis] + npts[axis]) % npts[axis];
            } else {
                // Clamp indices for non-periodic boundaries if point is exactly on max_coord
                // or slightly outside due to precision but not enough to trigger null_value earlier.
                if (corner_grid_indices[axis][0] >= npts[axis]) corner_grid_indices[axis][0] = npts[axis] - 1;
                if (corner_grid_indices[axis][1] >= npts[axis]) corner_grid_indices[axis][1] = npts[axis] - 1;
                // This also handles if base_indices[axis] was npts[axis]-1, then base_indices[axis]+1 would be npts[axis].
            }
        }
    }

    // Get values at the 8 corners
    // f000 = values[i0, j0, k0] etc.
    real_t f_ijk[2][2][2]; // Using 0/1 for low/high index along each axis
    for (int i_off = 0; i_off < 2; ++i_off) {
        for (int j_off = 0; j_off < 2; ++j_off) {
            for (int k_off = 0; k_off < 2; ++k_off) {
                Index3D corner_idx(
                    corner_grid_indices[0][i_off],
                    corner_grid_indices[1][j_off],
                    corner_grid_indices[2][k_off]
                );
                // Boundary check for actual access, though clamping should handle it for non-periodic.
                if (corner_idx.i1 >= npts[0] || corner_idx.i2 >= npts[1] || corner_idx.i3 >= npts[2]){
                    // This case should ideally be caught by clamping or periodicity.
                    // If point is outside non-periodic boundary, null_value returned earlier.
                    // If point is on boundary, clamping makes indices valid.
                    // This might occur if logic error in index calculation.
                    return null_value; // Safety return
                }
                f_ijk[i_off][j_off][k_off] = _values[field_grid.get_flat_index(corner_idx)];
            }
        }
    }

    // Interpolate along x (axis 0)
    real_t f_00 = f_ijk[0][0][0] * (1.0 - deltas[0]) + f_ijk[1][0][0] * deltas[0];
    real_t f_01 = f_ijk[0][0][1] * (1.0 - deltas[0]) + f_ijk[1][0][1] * deltas[0];
    real_t f_10 = f_ijk[0][1][0] * (1.0 - deltas[0]) + f_ijk[1][1][0] * deltas[0];
    real_t f_11 = f_ijk[0][1][1] * (1.0 - deltas[0]) + f_ijk[1][1][1] * deltas[0];

    // Interpolate along y (axis 1)
    real_t f_0 = f_00 * (1.0 - deltas[1]) + f_10 * deltas[1];
    real_t f_1 = f_01 * (1.0 - deltas[1]) + f_11 * deltas[1];

    // Interpolate along z (axis 2)
    real_t final_value = f_0 * (1.0 - deltas[2]) + f_1 * deltas[2];

    return final_value;
}


Point3D trilinear_vector(
    const FieldBase& field_grid,
    const std::vector<real_t>& _values, // Flattened N0*N1*N2*3
    const Point3D& point,
    const Point3D& null_value)
{
    // Similar logic to trilinear_scalar, but interpolate each component (x,y,z) of the vector separately.
    const auto& min_coords = field_grid.get_min_coords();
    const auto& intervals = field_grid.get_node_intervals();
    const auto& npts = field_grid.get_npts();
    const auto& is_periodic = field_grid.get_axis_is_periodic();
    const auto& is_null_axis = field_grid.get_axis_is_null();

    Point3D normalized_idx_float;
    std::array<int_t, 3> base_indices;
    Point3D deltas;

    for (size_t axis = 0; axis < 3; ++axis) {
        if (is_null_axis[axis]) {
             if (intervals[axis] < 1e-9) { normalized_idx_float[axis] = 0.0; }
             else { normalized_idx_float[axis] = (point[axis] - min_coords[axis]) / intervals[axis];}
            base_indices[axis] = 0;
            deltas[axis] = std::fmod(normalized_idx_float[axis], 1.0);
            if (normalized_idx_float[axis] < 0.0 && deltas[axis] != 0.0) deltas[axis] += 1.0;

        } else {
            if (!is_periodic[axis] && (point[axis] < min_coords[axis] - 1e-9 || point[axis] > field_grid.get_max_coords()[axis] + 1e-9)) {
                return null_value;
            }
            normalized_idx_float[axis] = (point[axis] - min_coords[axis]) / intervals[axis];
            base_indices[axis] = static_cast<int_t>(std::floor(normalized_idx_float[axis]));
            deltas[axis] = normalized_idx_float[axis] - static_cast<real_t>(base_indices[axis]);
        }
    }

    std::array<std::array<uint_t, 2>, 3> corner_grid_indices;
    for (size_t axis = 0; axis < 3; ++axis) {
        if (is_null_axis[axis]) {
            corner_grid_indices[axis][0] = 0;
            corner_grid_indices[axis][1] = 0;
        } else {
            corner_grid_indices[axis][0] = static_cast<uint_t>(base_indices[axis]);
            corner_grid_indices[axis][1] = static_cast<uint_t>(base_indices[axis] + 1);
            if (is_periodic[axis]) {
                corner_grid_indices[axis][0] = (corner_grid_indices[axis][0] % npts[axis] + npts[axis]) % npts[axis];
                corner_grid_indices[axis][1] = (corner_grid_indices[axis][1] % npts[axis] + npts[axis]) % npts[axis];
            } else {
                if (corner_grid_indices[axis][0] >= npts[axis]) corner_grid_indices[axis][0] = npts[axis] - 1;
                if (corner_grid_indices[axis][1] >= npts[axis]) corner_grid_indices[axis][1] = npts[axis] - 1;
            }
        }
    }
    
    Point3D interpolated_vector_value = {0.0, 0.0, 0.0};

    for (size_t comp = 0; comp < 3; ++comp) { // Interpolate for x, y, z components of the vector
        real_t f_ijk[2][2][2];
        for (int i_off = 0; i_off < 2; ++i_off) {
            for (int j_off = 0; j_off < 2; ++j_off) {
                for (int k_off = 0; k_off < 2; ++k_off) {
                    Index3D corner_idx(
                        corner_grid_indices[0][i_off],
                        corner_grid_indices[1][j_off],
                        corner_grid_indices[2][k_off]
                    );
                    if (corner_idx.i1 >= npts[0] || corner_idx.i2 >= npts[1] || corner_idx.i3 >= npts[2]){
                        return null_value; 
                    }
                    size_t flat_base_idx = field_grid.get_flat_index(corner_idx) * 3; // Each node has 3 components
                    f_ijk[i_off][j_off][k_off] = _values[flat_base_idx + comp];
                }
            }
        }

        real_t f_00 = f_ijk[0][0][0] * (1.0 - deltas[0]) + f_ijk[1][0][0] * deltas[0];
        real_t f_01 = f_ijk[0][0][1] * (1.0 - deltas[0]) + f_ijk[1][0][1] * deltas[0];
        real_t f_10 = f_ijk[0][1][0] * (1.0 - deltas[0]) + f_ijk[1][1][0] * deltas[0];
        real_t f_11 = f_ijk[0][1][1] * (1.0 - deltas[0]) + f_ijk[1][1][1] * deltas[0];

        real_t f_0 = f_00 * (1.0 - deltas[1]) + f_10 * deltas[1];
        real_t f_1 = f_01 * (1.0 - deltas[1]) + f_11 * deltas[1];

        interpolated_vector_value[comp] = f_0 * (1.0 - deltas[2]) + f_1 * deltas[2];
    }

    return interpolated_vector_value;
}


} // namespace interpolation
} // namespace Ckonal