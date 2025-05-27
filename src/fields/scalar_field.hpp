#ifndef Ckonal_SCALAR_FIELD_HPP
#define Ckonal_SCALAR_FIELD_HPP

#include "field_base.hpp"
#include <vector>
#include <limits> // For NaN
#include <memory>

// Forward declaration for VectorField3D if gradient returns it
namespace Ckonal { class VectorField3D; }

namespace Ckonal {

class ScalarField3D : public FieldBase {
public:
    ScalarField3D(CoordSys coord_sys = CoordSys::CARTESIAN);

    // --- Overriding setters to handle m_values resizing ---
    void set_npts(const std::array<uint_t, 3>& npts_val) override;

    // --- Value Access ---
    real_t get_value(const Index3D& idx) const;
    real_t get_value(size_t i1, size_t i2, size_t i3) const;
    void set_value(const Index3D& idx, real_t value);
    void set_value(size_t i1, size_t i2, size_t i3, real_t value);

    // Get direct access to the underlying data vector (e.g., for Heap or direct manipulation)
    // Be careful with const correctness.
    const std::vector<real_t>& get_values_raw() const;
    std::vector<real_t>& get_values_raw_mut(); // Use with caution

    // Set all values from an external vector
    void set_all_values(const std::vector<real_t>& all_vals);
    void fill(real_t fill_value);


    // --- Interpolation ---
    // `point` is in the field's own coordinate system
    real_t interpolate_value(const Point3D& point, real_t null_value = NaN) const;
    std::vector<real_t> resample_values(const std::vector<Point3D>& points, real_t null_value = NaN) const;

    // --- Gradient ---
    // Returns a new VectorField3D object representing the gradient.
    // The caller owns the returned object.
    std::unique_ptr<VectorField3D> calculate_gradient() const;

    // --- Ray Tracing ---
    // `end_point` is in the field's own coordinate system
    // Returns a vector of Point3D representing the ray path.
    // std::vector<Point3D> trace_ray_to(const Point3D& end_point) const; // trace_ray in python was from end

private:
    std::vector<real_t> m_values; // Flattened 1D vector for N0xN1xN2 scalar values

    // Specific gradient calculation methods
    std::unique_ptr<VectorField3D> gradient_cartesian() const;
    std::unique_ptr<VectorField3D> gradient_spherical() const;
};

} // namespace Ckonal

#endif // Ckonal_SCALAR_FIELD_HPP