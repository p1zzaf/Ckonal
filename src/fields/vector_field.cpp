#include "vector_field.hpp"
#include "numerics/interpolation.hpp" // For trilinear_vector
#include <stdexcept>
#include <algorithm> // For std::fill_n (though direct loop is also fine)

namespace Ckonal {

// --- Constructor ---
VectorField3D::VectorField3D(CoordSys coord_sys)
    : FieldBase(coord_sys, FieldType::VECTOR) {
    // Initialize m_values based on default npts from FieldBase
    m_values.assign(get_total_nodes() * 3, NaN); // Each node has 3 components
}

// --- Overriding set_npts to handle m_values resizing ---
void VectorField3D::set_npts(const std::array<uint_t, 3>& npts_val) {
    FieldBase::set_npts(npts_val); // Call base
    m_values.assign(get_total_nodes() * 3, NaN);
}

// --- Value Access ---
Point3D VectorField3D::get_value(const Index3D& idx) const {
    size_t base_flat_idx = get_flat_index(idx) * 3;
    return {m_values[base_flat_idx], m_values[base_flat_idx + 1], m_values[base_flat_idx + 2]};
}

Point3D VectorField3D::get_value(size_t i1, size_t i2, size_t i3) const {
    return get_value(Index3D(i1, i2, i3));
}

void VectorField3D::set_value(const Index3D& idx, const Point3D& value) {
    size_t base_flat_idx = get_flat_index(idx) * 3;
    m_values[base_flat_idx]     = value[0];
    m_values[base_flat_idx + 1] = value[1];
    m_values[base_flat_idx + 2] = value[2];
}

void VectorField3D::set_value(size_t i1, size_t i2, size_t i3, const Point3D& value) {
    set_value(Index3D(i1, i2, i3), value);
}

const std::vector<real_t>& VectorField3D::get_values_raw() const {
    return m_values;
}

std::vector<real_t>& VectorField3D::get_values_raw_mut() {
    return m_values;
}

void VectorField3D::set_all_values(const std::vector<real_t>& all_vals) {
    if (all_vals.size() != get_total_nodes() * 3) {
        throw std::invalid_argument("Size of input vector does not match field dimensions for set_all_values (vector field).");
    }
    m_values = all_vals;
}

void VectorField3D::fill(const Point3D& fill_value) {
    for (size_t i = 0; i < get_total_nodes(); ++i) {
        size_t base_idx = i * 3;
        m_values[base_idx]     = fill_value[0];
        m_values[base_idx + 1] = fill_value[1];
        m_values[base_idx + 2] = fill_value[2];
    }
}

// --- Interpolation ---
Point3D VectorField3D::interpolate_value(const Point3D& point, const Point3D& null_value) const {
    // Delegate to the free function in interpolation namespace
    return interpolation::trilinear_vector(*this, m_values, point, null_value);
}

// std::vector<Point3D> VectorField3D::resample_values(const std::vector<Point3D>& points, const Point3D& null_value) const {
//     std::vector<Point3D> resampled_data;
//     resampled_data.reserve(points.size());
//     for (const auto& p : points) {
//         resampled_data.push_back(interpolate_value(p, null_value));
//     }
//     return resampled_data;
// }

} // namespace Ckonal