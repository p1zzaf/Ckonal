#include "gradient.hpp"
#include "../fields/scalar_field.hpp" // Need full definition for access
#include "../fields/vector_field.hpp" // For constructing result if needed, though here we fill raw vector
#include <cmath>     // For sin
#include <stdexcept> // For errors

namespace Ckonal {
namespace gradient {

// Helper to get value from scalar_field, handling boundary conditions for gradient stencils
// For gradient, we often need points like (i-1, j, k), (i+1, j, k), etc.
// This helper needs to respect periodicity and grid bounds.
// 'idx_offset' is [-2, -1, 0, 1, 2] typically.
// 'axis' is the dimension along which the offset is applied.
// 'base_idx' is the central point for the stencil.
inline real_t get_value_for_stencil(const ScalarField3D& field, const Index3D& base_idx, int offset, size_t axis_dim) {
    Index3D current_idx = base_idx;
    int_t target_component_idx = static_cast<int_t>(base_idx.i1); // Example for axis 0
    if (axis_dim == 0) target_component_idx = static_cast<int_t>(base_idx.i1) + offset;
    else if (axis_dim == 1) target_component_idx = static_cast<int_t>(base_idx.i2) + offset;
    else if (axis_dim == 2) target_component_idx = static_cast<int_t>(base_idx.i3) + offset;

    const auto& npts = field.get_npts();
    const auto& is_periodic = field.get_axis_is_periodic();

    if (is_periodic[axis_dim]) {
        target_component_idx = (target_component_idx % static_cast<int_t>(npts[axis_dim]) + static_cast<int_t>(npts[axis_dim])) % static_cast<int_t>(npts[axis_dim]);
    } else {
        if (target_component_idx < 0 || target_component_idx >= static_cast<int_t>(npts[axis_dim])) {
            // This should ideally not happen if boundary stencils are chosen correctly.
            // If it does, means we need a value outside the non-periodic domain.
            // Python's np.gradient handles boundaries automatically. We need to emulate.
            // For now, this indicates an issue with stencil choice at boundary.
            // Let's assume the calling logic (boundary vs interior) handles this.
            // If we are at a boundary and need an "outside" point for a central diff,
            // we should be using a one-sided diff instead.
            // So, this path suggests a logic error in the caller for non-periodic.
            // However, np.gradient implicitly uses lower order at boundaries.
            // Ckonal's _gradient_of_spherical uses explicit 2nd order forward/backward at edges.
            throw std::out_of_range("Stencil requires point outside non-periodic boundary.");
        }
    }

    if (axis_dim == 0) current_idx.i1 = static_cast<size_t>(target_component_idx);
    else if (axis_dim == 1) current_idx.i2 = static_cast<size_t>(target_component_idx);
    else if (axis_dim == 2) current_idx.i3 = static_cast<size_t>(target_component_idx);
    
    return field.get_value(current_idx);
}


void calculate_cartesian(const ScalarField3D& scalar_field, std::vector<real_t>& result_values) {
    const auto& npts = scalar_field.get_npts();
    const auto& intervals = scalar_field.get_node_intervals();
    const auto& is_null_axis = scalar_field.get_axis_is_null();

    result_values.assign(scalar_field.get_total_nodes() * 3, 0.0); // Initialize with zeros

    for (uint_t i = 0; i < npts[0]; ++i) {
        for (uint_t j = 0; j < npts[1]; ++j) {
            for (uint_t k = 0; k < npts[2]; ++k) {
                Index3D current_idx(i, j, k);
                size_t flat_grad_idx = scalar_field.get_flat_index(current_idx) * 3;

                for (size_t axis = 0; axis < 3; ++axis) {
                    if (is_null_axis[axis]) {
                        result_values[flat_grad_idx + axis] = 0.0;
                        continue;
                    }

                    real_t val_plus1, val_minus1;
                    real_t h = intervals[axis];

                    // Determine current component index (i, j, or k)
                    uint_t current_comp_val = 0;
                    if (axis == 0) current_comp_val = i;
                    else if (axis == 1) current_comp_val = j;
                    else current_comp_val = k;

                    if (npts[axis] == 1) { // Should be caught by is_null_axis, but double check
                        result_values[flat_grad_idx + axis] = 0.0;
                        continue;
                    }
                    
                    // np.gradient uses central differences in interior, and first/last differences at boundaries.
                    // For 2nd order accuracy:
                    if (current_comp_val == 0 && !scalar_field.get_axis_is_periodic()[axis]) { // Forward difference at lower boundary
                        real_t val_0 = scalar_field.get_value(current_idx);
                        real_t val_1 = get_value_for_stencil(scalar_field, current_idx, 1, axis);
                        if (npts[axis] > 2) { // Need 3 points for 2nd order forward
                            real_t val_2 = get_value_for_stencil(scalar_field, current_idx, 2, axis);
                            result_values[flat_grad_idx + axis] = (-3.0 * val_0 + 4.0 * val_1 - val_2) / (2.0 * h);
                        } else { // Fallback to 1st order forward if only 2 points
                            result_values[flat_grad_idx + axis] = (val_1 - val_0) / h;
                        }
                    } else if (current_comp_val == npts[axis] - 1 && !scalar_field.get_axis_is_periodic()[axis]) { // Backward difference at upper boundary
                        real_t val_0 = scalar_field.get_value(current_idx);
                        real_t val_m1 = get_value_for_stencil(scalar_field, current_idx, -1, axis);
                         if (npts[axis] > 2) { // Need 3 points for 2nd order backward
                            real_t val_m2 = get_value_for_stencil(scalar_field, current_idx, -2, axis);
                            result_values[flat_grad_idx + axis] = (3.0 * val_0 - 4.0 * val_m1 + val_m2) / (2.0 * h);
                        } else { // Fallback to 1st order backward
                             result_values[flat_grad_idx + axis] = (val_0 - val_m1) / h;
                        }
                    } else { // Central difference for interior or periodic boundaries
                        val_plus1 = get_value_for_stencil(scalar_field, current_idx, 1, axis);
                        val_minus1 = get_value_for_stencil(scalar_field, current_idx, -1, axis);
                        result_values[flat_grad_idx + axis] = (val_plus1 - val_minus1) / (2.0 * h);
                    }
                }
            }
        }
    }
}


void calculate_spherical(const ScalarField3D& scalar_field, std::vector<real_t>& result_values) {
    const auto& npts = scalar_field.get_npts();
    const auto& intervals = scalar_field.get_node_intervals(); // dr, dtheta, dphi
    const auto& is_null_axis = scalar_field.get_axis_is_null();

    result_values.assign(scalar_field.get_total_nodes() * 3, 0.0);

    real_t dr = intervals[0];
    real_t dt = intervals[1];
    real_t dp = intervals[2];

    for (uint_t ir = 0; ir < npts[0]; ++ir) { // rho index
        for (uint_t it = 0; it < npts[1]; ++it) { // theta index
            for (uint_t ip = 0; ip < npts[2]; ++ip) { // phi index
                Index3D current_idx(ir, it, ip);
                size_t flat_grad_idx = scalar_field.get_flat_index(current_idx) * 3;
                Point3D node_coords = scalar_field.get_node_coords(current_idx); // [rho, theta, phi]
                real_t rho = node_coords[0];
                real_t theta = node_coords[1];
                // real_t phi = node_coords[2]; // Not directly used in denominators here

                // --- Gradient wrt rho (axis 0) ---
                if (is_null_axis[0]) {
                    result_values[flat_grad_idx + 0] = 0.0;
                } else {
                    if (npts[0] == 1) result_values[flat_grad_idx + 0] = 0.0;
                    else if (ir == 0 && !scalar_field.get_axis_is_periodic()[0]) { // Lower boundary (rho_min)
                        real_t v0 = scalar_field.get_value(current_idx);
                        real_t v1 = get_value_for_stencil(scalar_field, current_idx, 1, 0);
                        if (npts[0] > 2) {
                            real_t v2 = get_value_for_stencil(scalar_field, current_idx, 2, 0);
                            result_values[flat_grad_idx + 0] = (-3.0*v0 + 4.0*v1 - v2) / (2.0*dr);
                        } else {
                             result_values[flat_grad_idx + 0] = (v1 - v0) / dr;
                        }
                    } else if (ir == npts[0] - 1 && !scalar_field.get_axis_is_periodic()[0]) { // Upper boundary (rho_max)
                        real_t v0 = scalar_field.get_value(current_idx);
                        real_t v_1 = get_value_for_stencil(scalar_field, current_idx, -1, 0);
                        if (npts[0] > 2) {
                            real_t v_2 = get_value_for_stencil(scalar_field, current_idx, -2, 0);
                            result_values[flat_grad_idx + 0] = (3.0*v0 - 4.0*v_1 + v_2) / (2.0*dr);
                        } else {
                            result_values[flat_grad_idx + 0] = (v0 - v_1) / dr;
                        }
                    } else { // Interior or periodic
                        real_t v_p1 = get_value_for_stencil(scalar_field, current_idx, 1, 0);
                        real_t v_m1 = get_value_for_stencil(scalar_field, current_idx, -1, 0);
                        result_values[flat_grad_idx + 0] = (v_p1 - v_m1) / (2.0*dr);
                    }
                }

                // --- Gradient wrt theta (axis 1) ---
                // df/dtheta = (1/rho) * dV/d(angle_theta)
                if (is_null_axis[1]) {
                    result_values[flat_grad_idx + 1] = 0.0;
                } else {
                    real_t grad_val_angle_theta;
                    if (npts[1] == 1) grad_val_angle_theta = 0.0;
                    else if (it == 0 && !scalar_field.get_axis_is_periodic()[1]) {
                        real_t v0 = scalar_field.get_value(current_idx);
                        real_t v1 = get_value_for_stencil(scalar_field, current_idx, 1, 1);
                        if (npts[1] > 2) {
                            real_t v2 = get_value_for_stencil(scalar_field, current_idx, 2, 1);
                            grad_val_angle_theta = (-3.0*v0 + 4.0*v1 - v2) / (2.0*dt);
                        } else {
                            grad_val_angle_theta = (v1 - v0) / dt;
                        }
                    } else if (it == npts[1] - 1 && !scalar_field.get_axis_is_periodic()[1]) {
                        real_t v0 = scalar_field.get_value(current_idx);
                        real_t v_1 = get_value_for_stencil(scalar_field, current_idx, -1, 1);
                        if (npts[1] > 2) {
                            real_t v_2 = get_value_for_stencil(scalar_field, current_idx, -2, 1);
                            grad_val_angle_theta = (3.0*v0 - 4.0*v_1 + v_2) / (2.0*dt);
                        } else {
                            grad_val_angle_theta = (v0 - v_1) / dt;
                        }
                    } else {
                        real_t v_p1 = get_value_for_stencil(scalar_field, current_idx, 1, 1);
                        real_t v_m1 = get_value_for_stencil(scalar_field, current_idx, -1, 1);
                        grad_val_angle_theta = (v_p1 - v_m1) / (2.0*dt);
                    }
                    // Add 1/rho scaling, careful if rho is zero (should be caught by FieldBase validation)
                    result_values[flat_grad_idx + 1] = (std::abs(rho) < 1e-9) ? 0.0 : grad_val_angle_theta / rho;
                }

                // --- Gradient wrt phi (axis 2) ---
                // df/dphi = (1 / (rho*sin(theta))) * dV/d(angle_phi)
                if (is_null_axis[2]) {
                    result_values[flat_grad_idx + 2] = 0.0;
                } else {
                    real_t grad_val_angle_phi;
                     if (npts[2] == 1) grad_val_angle_phi = 0.0;
                    else if (ip == 0 && !scalar_field.get_axis_is_periodic()[2]) {
                        real_t v0 = scalar_field.get_value(current_idx);
                        real_t v1 = get_value_for_stencil(scalar_field, current_idx, 1, 2);
                        if (npts[2] > 2) {
                           real_t v2 = get_value_for_stencil(scalar_field, current_idx, 2, 2);
                           grad_val_angle_phi = (-3.0*v0 + 4.0*v1 - v2) / (2.0*dp);
                        } else {
                            grad_val_angle_phi = (v1-v0) / dp;
                        }
                    } else if (ip == npts[2] - 1 && !scalar_field.get_axis_is_periodic()[2]) {
                        real_t v0 = scalar_field.get_value(current_idx);
                        real_t v_1 = get_value_for_stencil(scalar_field, current_idx, -1, 2);
                        if (npts[2] > 2) {
                            real_t v_2 = get_value_for_stencil(scalar_field, current_idx, -2, 2);
                            grad_val_angle_phi = (3.0*v0 - 4.0*v_1 + v_2) / (2.0*dp);
                        } else {
                            grad_val_angle_phi = (v0 - v_1) / dp;
                        }
                    } else { // Interior or periodic
                        real_t v_p1 = get_value_for_stencil(scalar_field, current_idx, 1, 2);
                        real_t v_m1 = get_value_for_stencil(scalar_field, current_idx, -1, 2);
                        grad_val_angle_phi = (v_p1 - v_m1) / (2.0*dp);
                    }

                    real_t sin_theta = std::sin(theta);
                    // Careful with rho*sin(theta) being zero (poles or r=0)
                    real_t denominator = rho * sin_theta;
                    if (std::abs(denominator) < 1e-9) {
                        // At poles (sin_theta=0) or origin (rho=0), d/dphi is ill-defined or component should be 0.
                        // Python code would use the small denominator from `norm` property.
                        // If `norm` component for phi is 0, then aa,bb,cc for phi become 0 in solver.
                        // Here we are calculating the gradient itself.
                        // Let's set to 0 if denominator is too small, consistent with how it might be used.
                        result_values[flat_grad_idx + 2] = 0.0;
                    } else {
                        result_values[flat_grad_idx + 2] = grad_val_angle_phi / denominator;
                    }
                }
            }
        }
    }
}

} // namespace gradient
} // namespace Ckonal