#ifndef Ckonal_VECTOR_FIELD_HPP
#define Ckonal_VECTOR_FIELD_HPP

#include "field_base.hpp"
#include <vector>

namespace Ckonal {

class VectorField3D : public FieldBase {
public:
    VectorField3D(CoordSys coord_sys = CoordSys::CARTESIAN);

    // --- Overriding setters to handle m_values resizing ---
    void set_npts(const std::array<uint_t, 3>& npts_val) override;

    // --- Value Access ---
    Point3D get_value(const Index3D& idx) const; // Returns [vx, vy, vz]
    Point3D get_value(size_t i1, size_t i2, size_t i3) const;
    void set_value(const Index3D& idx, const Point3D& value);
    void set_value(size_t i1, size_t i2, size_t i3, const Point3D& value);

    // Get direct access to the underlying data vector
    const std::vector<real_t>& get_values_raw() const; // Data is [v0x,v0y,v0z, v1x,v1y,v1z, ...]
    std::vector<real_t>& get_values_raw_mut();

    // Set all values from an external vector
    void set_all_values(const std::vector<real_t>& all_vals); // Expects flattened N0*N1*N2*3
    void fill(const Point3D& fill_value);

    // --- Interpolation ---
    // `point` is in the field's own coordinate system
    Point3D interpolate_value(const Point3D& point, const Point3D& null_value = {NaN, NaN, NaN}) const;
    // std::vector<Point3D> resample_values(const std::vector<Point3D>& points, const Point3D& null_value) const;

private:
    // Flattened 1D vector for N0xN1xN2x3 vector values.
    // Order: [idx0_x, idx0_y, idx0_z, idx1_x, idx1_y, idx1_z, ...]
    std::vector<real_t> m_values;
};

} // namespace Ckonal

#endif // Ckonal_VECTOR_FIELD_HPP