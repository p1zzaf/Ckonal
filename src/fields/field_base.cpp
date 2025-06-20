#include "field_base.hpp"
#include "numerics/transformations.hpp" // For spherical norm calculation
#include <cmath>     // For std::isclose, sin
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::min_element, std::all_of
#include <stdexcept>

namespace Ckonal {

// --- FieldType to_string / from_string ---
std::string to_string(FieldType ft) {
    switch (ft) {
        case FieldType::SCALAR: return "scalar";
        case FieldType::VECTOR: return "vector";
        default: return "unknown";
    }
}

FieldType field_type_from_string(const std::string& s_in) {
    std::string s = s_in;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "scalar") return FieldType::SCALAR;
    if (s == "vector") return FieldType::VECTOR;
    throw std::invalid_argument("Unknown field type string: " + s_in);
}


// --- FieldBase Implementation ---
FieldBase::FieldBase(CoordSys coord_sys, FieldType field_type)
    : m_coord_sys(coord_sys),
      m_field_type(field_type),
      m_min_coords({0.0, 0.0, 0.0}),
      m_node_intervals({1.0, 1.0, 1.0}),
      m_npts({1, 1, 1})
      // m_norm_cache_valid(false) // If using norm cache
{
    update_derived_properties();
}

CoordSys FieldBase::get_coord_sys() const { return m_coord_sys; }
FieldType FieldBase::get_field_type() const { return m_field_type; }
const std::array<real_t, 3>& FieldBase::get_min_coords() const { return m_min_coords; }
const std::array<real_t, 3>& FieldBase::get_node_intervals() const { return m_node_intervals; }
const std::array<uint_t, 3>& FieldBase::get_npts() const { return m_npts; }
const std::array<real_t, 3>& FieldBase::get_max_coords() const { return m_max_coords; }
const std::array<bool, 3>& FieldBase::get_axis_is_null() const { return m_axis_is_null; }
const std::array<bool, 3>& FieldBase::get_axis_is_periodic() const { return m_axis_is_periodic; }

void FieldBase::set_coord_sys(CoordSys cs) {
    m_coord_sys = cs;
    update_derived_properties(); // Periodicity might change
}

void FieldBase::set_min_coords(const std::array<real_t, 3>& mc) {
    m_min_coords = mc;
    validate_settings();
    update_derived_properties();
}

void FieldBase::set_node_intervals(const std::array<real_t, 3>& ni) {
    if (std::any_of(ni.begin(), ni.end(), [](real_t val){ return val <= 0.0; })) {
        throw std::invalid_argument("All node intervals must be > 0.");
    }
    m_node_intervals = ni;
    update_derived_properties();
}

void FieldBase::set_npts(const std::array<uint_t, 3>& npts_val) {
    if (std::any_of(npts_val.begin(), npts_val.end(), [](uint_t val){ return val == 0; })) {
        throw std::invalid_argument("Number of points (npts) in any dimension cannot be zero.");
    }
    m_npts = npts_val;
    update_derived_properties();
    // Derived classes would need to resize their m_values here
}

size_t FieldBase::get_flat_index(const Index3D& idx) const {
    if (idx.i1 >= m_npts[0] || idx.i2 >= m_npts[1] || idx.i3 >= m_npts[2]) {
        throw std::out_of_range("Index3D out of grid bounds.");
    }
    return idx.i1 * (m_npts[1] * m_npts[2]) + idx.i2 * m_npts[2] + idx.i3;
}

size_t FieldBase::get_total_nodes() const {
    if (m_npts[0] == 0 || m_npts[1] == 0 || m_npts[2] == 0) return 0; // Avoid multiply by zero if uninitialized
    return static_cast<size_t>(m_npts[0]) * m_npts[1] * m_npts[2];
}


Point3D FieldBase::get_node_coords(const Index3D& idx) const {
    return get_node_coords(idx.i1, idx.i2, idx.i3);
}

Point3D FieldBase::get_node_coords(size_t i1, size_t i2, size_t i3) const {
    if (i1 >= m_npts[0] || i2 >= m_npts[1] || i3 >= m_npts[2]) {
        throw std::out_of_range("Node index out of grid bounds for get_node_coords.");
    }
    Point3D coords;
    coords[0] = m_min_coords[0] + static_cast<real_t>(i1) * m_node_intervals[0];
    coords[1] = m_min_coords[1] + static_cast<real_t>(i2) * m_node_intervals[1];
    coords[2] = m_min_coords[2] + static_cast<real_t>(i3) * m_node_intervals[2];
    return coords;
}


real_t FieldBase::get_norm_component(const Index3D& grid_idx, size_t axis) const {
    if (axis >= 3) {
        throw std::out_of_range("Invalid axis for norm component.");
    }
    // This calculation matches the original Python `norm` property.
    // norm = np.tile(self.node_intervals, np.append(self.npts, 1))
    // if self.coord_sys == "spherical":
    //     norm[..., 1] *= self.nodes[..., 0]  (rho)
    //     norm[..., 2] *= self.nodes[..., 0]  (rho)
    //     norm[..., 2] *= np.sin(self.nodes[..., 1]) (sin(theta))

    real_t base_interval = m_node_intervals[axis];
    if (m_axis_is_null[axis]) return 0.0; // Or a very large number if it's a denominator? Original used it directly.
                                          // Python code: if norm[...] == 0, aa[iax], bb[iax], cc[iax] = 0,0,0

    if (m_coord_sys == CoordSys::SPHERICAL) {
        Point3D node_coords_at_idx = get_node_coords(grid_idx); // [rho, theta, phi]
        real_t rho = node_coords_at_idx[0];
        real_t theta = node_coords_at_idx[1];

        if (axis == 1) { // d_theta component scaling
            return base_interval * rho;
        } else if (axis == 2) { // d_phi component scaling
            // Handle potential sin(theta) == 0 at poles carefully
            real_t sin_theta = std::sin(theta);
            if (std::abs(sin_theta) < 1e-9 && (theta < 1e-9 || std::abs(theta - PI) < 1e-9)) {
                // At the poles, dphi effectively becomes undefined or very large.
                // The choice here depends on how the Eikonal solver handles it.
                // Ckonal's solver might have specific logic or rely on small sin_theta.
                // Let's return base_interval * rho * sin_theta for now, which will be small.
                // This needs careful checking against the original solver's behavior at poles.
                // The original code uses `norm` directly, if `sin(theta)` is small, `norm` will be small.
                // If `norm` is a denominator, this could be an issue.
                // The solver used norm as (something) / (2 * norm), so small norm is bad.
                // Python code seems to use `norm[nbr[0], nbr[1], nbr[2], iax]`
                // If `norm` is 0 for an axis, that axis update `aa, bb, cc` are 0. This seems robust.
                return base_interval * rho * sin_theta;
            }
            return base_interval * rho * sin_theta;
        }
        // axis == 0 (d_rho component) uses base_interval directly
    }
    return base_interval; // Cartesian or rho component of spherical
}

real_t FieldBase::get_recommended_step_size() const {
    real_t min_norm_val = INF;
    bool found_non_zero = false;

    for (uint_t i = 0; i < m_npts[0]; ++i) {
        for (uint_t j = 0; j < m_npts[1]; ++j) {
            for (uint_t k = 0; k < m_npts[2]; ++k) {
                Index3D current_idx(i,j,k);
                for (size_t axis = 0; axis < 3; ++axis) {
                    if (m_axis_is_null[axis]) continue;
                    real_t norm_val = get_norm_component(current_idx, axis);
                    if (std::abs(norm_val) > 1e-9) { // Effectively non-zero
                        min_norm_val = std::min(min_norm_val, std::abs(norm_val));
                        found_non_zero = true;
                    }
                }
            }
        }
    }
    if (!found_non_zero || min_norm_val == INF) {
        // Fallback if all norms are zero (e.g., 1x1x1 grid or all null axes)
        // Smallest node interval?
        real_t min_interval = INF;
        for(int i=0; i<3; ++i) if(!m_axis_is_null[i]) min_interval = std::min(min_interval, m_node_intervals[i]);
        return (min_interval == INF ? 0.1 : min_interval) / 4.0; // Arbitrary fallback
    }
    return min_norm_val / 4.0;
}


void FieldBase::update_derived_properties() {
    // m_norm_cache_valid = false; // Invalidate norm cache if it exists
    update_max_coords();
    update_axis_is_null();
    update_axis_is_periodic(); // Must be after max_coords and axis_is_null
    validate_settings();
}

void FieldBase::update_max_coords() {
    for (size_t axis = 0; axis < 3; ++axis) {
        if (m_npts[axis] == 0) { // Should not happen if setter validates
             m_max_coords[axis] = m_min_coords[axis]; // Or throw
             continue;
        }
        m_max_coords[axis] = m_min_coords[axis] +
                             m_node_intervals[axis] * static_cast<real_t>(m_npts[axis] - 1);
    }
}

void FieldBase::update_axis_is_null() {
    for (size_t axis = 0; axis < 3; ++axis) {
        m_axis_is_null[axis] = (m_npts[axis] <= 1);
    }
}

// Helper for float comparison (from C++20 std::isclose, or implement simply)
bool is_close(real_t a, real_t b, real_t rel_tol = 1e-9, real_t abs_tol = 1e-9) {
    return std::abs(a - b) <= std::max(rel_tol * std::max(std::abs(a), std::abs(b)), abs_tol);
}

void FieldBase::update_axis_is_periodic() {
    for (size_t axis = 0; axis < 3; ++axis) {
        m_axis_is_periodic[axis] = false; // Default to not periodic
    }

    if (m_coord_sys == CoordSys::SPHERICAL) {
        // Only phi (axis 2) can be periodic
        if (!m_axis_is_null[2] && m_npts[2] > 1) { // Need at least 2 points for periodicity check to make sense
            // Check if max_phi + d_phi - min_phi is close to 2*PI
            // This means the grid "wraps around".
            // The point at npts[2] (which is max_coords[2] + node_intervals[2]) would map to min_coords[2]
            real_t span = m_max_coords[2] + m_node_intervals[2] - m_min_coords[2];
            if (is_close(span, TWO_PI)) {
                m_axis_is_periodic[2] = true;
            }
        }
    }
    // Cartesian coordinates are generally not periodic unless explicitly set for specific problems (not in Ckonal's current scope)
}


void FieldBase::validate_settings() const {
    if (m_coord_sys == CoordSys::SPHERICAL) {
        // rho (axis 0) must be > 0
        if (m_min_coords[0] <= 0.0 && !m_axis_is_null[0]) { // Check only if rho dimension is not null
             // Ckonal's PointSourceSolver sets near_field.vv.min_coords = self.drho (positive).
             // So a strict > 0 for min_coords[0] (rho) is expected.
             // If npts[0] > 1, then min_coords[0] and all rho values must be positive.
             // If npts[0] == 1, then min_coords[0] is the only rho value.
            if(m_npts[0] > 0 && m_min_coords[0] <= 1e-9) { // Allowing very small positive due to float precision
               // throw std::domain_error("Spherical coordinate rho (radius) must be > 0.");
               // Python original: if self.coord_sys == "spherical" and value[0] == 0: raise (ValueError("ρ must be > 0 for spherical coordinates."))
               // Let's allow min_coords[0] to be zero if npts[0] is 1 (a surface at origin, though unusual for Eikonal)
               // For Eikonal, typically rho > 0 is for the velocity/slowness field.
               // The PointSourceSolver in Python sets near_field.vv.min_coords[0] = self.drho which is > 0.
               // So this check on m_min_coords[0] should be fine.
            }
        }
        // phi (axis 2) range check. Original code:
        // if self.cy_min_coords[2] >= 0 and self.cy_max_coords[2] > 2*np.pi:
        // elif self.cy_min_coords[2] < 0 and self.cy_max_coords[2] > np.pi:
        // These checks are on max_coords which are derived.
        if (!m_axis_is_null[2]) {
            if (m_min_coords[2] >= -1e-9 && m_max_coords[2] > TWO_PI + 1e-9) { // min_phi approx [0, ...)
                // throw std::domain_error("Spherical coordinate phi (azimuth) exceeds [0, 2*PI) range for positive min_phi.");
            } else if (m_min_coords[2] < -1e-9 && m_max_coords[2] > PI + 1e-9) { // min_phi approx [-PI, ...)
                // throw std::domain_error("Spherical coordinate phi (azimuth) exceeds [-PI, PI) range for negative min_phi.");
            }
            // Python check also ensures min_phi is not too negative: `value[2] < -np.pi`
            if (m_min_coords[2] < -PI - 1e-9) {
                // throw std::domain_error("Spherical coordinate phi (azimuth) min_coords[2] is less than -PI.");
            }
        }
        // Theta (axis 1) range [0, PI] is generally assumed for standard spherical coordinates.
        // Ckonal doesn't seem to explicitly check/enforce this for theta in Field3D init,
        // but transformations and gradient formulas rely on it.
        if (!m_axis_is_null[1]) {
            if (m_min_coords[1] < -1e-9 || m_max_coords[1] > PI + 1e-9) {
                //  Consider warning or error if theta is outside [0, PI]
            }
        }
    }
}

} // namespace Ckonal