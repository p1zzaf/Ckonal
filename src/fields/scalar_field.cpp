#include "scalar_field.hpp"
#include "vector_field.hpp" // For return type of calculate_gradient
#include "numerics/interpolation.hpp" // For trilinear_scalar
#include "numerics/gradient.hpp"      // For gradient::calculate_cartesian/spherical
#include <stdexcept>
#include <algorithm> // For std::fill

namespace Ckonal {

// --- Constructor ---
ScalarField3D::ScalarField3D(CoordSys coord_sys)
    : FieldBase(coord_sys, FieldType::SCALAR) {
    // Initialize m_values based on default npts from FieldBase
    // FieldBase constructor calls update_derived_properties(), which sets initial m_npts.
    m_values.assign(get_total_nodes(), NaN); // Fill with NaN by default
}

// --- Overriding set_npts to handle m_values resizing ---
void ScalarField3D::set_npts(const std::array<uint_t, 3>& npts_val) {
    FieldBase::set_npts(npts_val); // Call base to update m_npts and derived props
    // Resize m_values and fill new elements with NaN, old data might be lost or partially kept
    // depending on std::vector::resize behavior if shrinking.
    // For simplicity, let's reassign like Python (np.full)
    m_values.assign(get_total_nodes(), NaN);
}

// --- Value Access ---
real_t ScalarField3D::get_value(const Index3D& idx) const {
    // get_flat_index will throw std::out_of_range if idx is invalid
    return m_values[get_flat_index(idx)];
}

real_t ScalarField3D::get_value(size_t i1, size_t i2, size_t i3) const {
    return get_value(Index3D(i1, i2, i3));
}

void ScalarField3D::set_value(const Index3D& idx, real_t value) {
    m_values[get_flat_index(idx)] = value;
}

void ScalarField3D::set_value(size_t i1, size_t i2, size_t i3, real_t value) {
    set_value(Index3D(i1, i2, i3), value);
}

const std::vector<real_t>& ScalarField3D::get_values_raw() const {
    return m_values;
}

std::vector<real_t>& ScalarField3D::get_values_raw_mut() {
    return m_values;
}

void ScalarField3D::set_all_values(const std::vector<real_t>& all_vals) {
    if (all_vals.size() != get_total_nodes()) {
        throw std::invalid_argument("Size of input vector does not match field dimensions for set_all_values.");
    }
    m_values = all_vals;
}

void ScalarField3D::fill(real_t fill_value) {
    std::fill(m_values.begin(), m_values.end(), fill_value);
}

// --- Interpolation ---
real_t ScalarField3D::interpolate_value(const Point3D& point, real_t null_value) const {
    // Delegate to the free function in interpolation namespace
    return interpolation::trilinear_scalar(*this, m_values, point, null_value);
}

std::vector<real_t> ScalarField3D::resample_values(const std::vector<Point3D>& points, real_t null_value) const {
    std::vector<real_t> resampled_data;
    resampled_data.reserve(points.size());
    for (const auto& p : points) {
        resampled_data.push_back(interpolate_value(p, null_value));
    }
    return resampled_data;
}

// --- Gradient ---
std::unique_ptr<VectorField3D> ScalarField3D::calculate_gradient() const {
    // Create a new VectorField3D to store the gradient
    // The new field will have the same grid definition as this scalar field
    auto grad_vector_field = std::make_unique<VectorField3D>(get_coord_sys());
    grad_vector_field->set_min_coords(get_min_coords());
    grad_vector_field->set_node_intervals(get_node_intervals());
    // set_npts will also correctly size the internal m_values of grad_vector_field
    grad_vector_field->set_npts(get_npts());

    // Call the appropriate gradient calculation function based on the coordinate system
    // The gradient functions will fill the m_values of grad_vector_field.
    if (get_coord_sys() == CoordSys::CARTESIAN) {
        gradient::calculate_cartesian(*this, grad_vector_field->get_values_raw_mut());
    } else if (get_coord_sys() == CoordSys::SPHERICAL) {
        gradient::calculate_spherical(*this, grad_vector_field->get_values_raw_mut());
    } else {
        // Should not happen if CoordSys enum is exhaustive and handled
        throw std::runtime_error("Gradient calculation not implemented for the current coordinate system.");
    }

    return grad_vector_field;
}

// --- Ray Tracing (Placeholder) ---
// std::vector<Point3D> ScalarField3D::trace_ray_to(const Point3D& end_point) const {
//     // 1. Calculate gradient field
//     auto grad_field = calculate_gradient();
//     // 2. Get recommended step size
//     real_t step = get_recommended_step_size(); // or grad_field->get_recommended_step_size()
//     // 3. Iteratively step from end_point backwards using gradient and interpolated values
//     // ... implementation ...
//     std::vector<Point3D> ray_path;
//     // ... populate ray_path ...
//     // std::reverse(ray_path.begin(), ray_path.end()); // Python version reverses at the end
//     return ray_path;
// }

} // namespace Ckonal