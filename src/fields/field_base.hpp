#ifndef Ckonal_FIELD_BASE_HPP
#define Ckonal_FIELD_BASE_HPP

#include "core/konal_types.hpp"
#include "core/konal_constants.hpp"
#include <vector>
#include <string>
#include <stdexcept> // For exceptions
#include <array>

// Forward declaration for serializer if we make it a friend or use it directly
namespace Ckonal { namespace io { class FieldSerializer; } }

namespace Ckonal {

enum class FieldType {
    SCALAR,
    VECTOR,
    UNKNOWN
};

std::string to_string(FieldType ft);
FieldType field_type_from_string(const std::string& s);


class FieldBase {
public:
    // --- Constructors ---
    FieldBase(CoordSys coord_sys = CoordSys::CARTESIAN, FieldType field_type = FieldType::UNKNOWN);
    virtual ~FieldBase() = default; // Important for base class with virtual functions

    // ---禁止拷贝和移动，或显式默认/删除---
    FieldBase(const FieldBase&) = delete;
    FieldBase& operator=(const FieldBase&) = delete;
    FieldBase(FieldBase&&) = default;
    FieldBase& operator=(FieldBase&&) = default;

    // --- Getters ---
    CoordSys get_coord_sys() const;
    FieldType get_field_type() const;
    const std::array<real_t, 3>& get_min_coords() const;
    const std::array<real_t, 3>& get_node_intervals() const;
    const std::array<uint_t, 3>& get_npts() const;
    const std::array<real_t, 3>& get_max_coords() const; // Derived
    const std::array<bool, 3>& get_axis_is_null() const; // Derived
    const std::array<bool, 3>& get_axis_is_periodic() const; // Derived

    // --- Setters ---
    // Setting these will trigger updates to derived properties like max_coords, axis_is_null, etc.
    // And potentially resize internal data arrays if they exist at this level or in derived classes.
    virtual void set_coord_sys(CoordSys cs); // May affect periodicity
    virtual void set_min_coords(const std::array<real_t, 3>& mc);
    virtual void set_node_intervals(const std::array<real_t, 3>& ni);
    virtual void set_npts(const std::array<uint_t, 3>& npts);

    // --- Grid Node Coordinates ---
    // Generates all node coordinates. Could be memory intensive for large grids.
    // Returns a flattened vector of Point3D: [p000, p001, ..., pN0N1N2]
    // Or perhaps it's better to get coordinates for a specific (i,j,k) index?
    // The original `nodes` property returned a (N0,N1,N2,3) numpy array.
    // For C++, returning a huge vector might not be ideal.
    // Let's provide a method to get a single node's coordinates.
    Point3D get_node_coords(const Index3D& idx) const;
    Point3D get_node_coords(size_t i1, size_t i2, size_t i3) const;

    // --- Gradient Scaling Factors (Norm property in Python) ---
    // This is (dx, dy, dz) effectively, but can be more complex for spherical.
    // It's a 4D array in Python (N0,N1,N2,3).
    // We'll need a similar structure or a way to compute these on the fly.
    // For now, let's define a getter that derived classes might use/override.
    // It might be better to compute these within gradient/solver logic where needed.
    // The original `norm` property calculated and cached it.
    // Let's provide access to the components.
    real_t get_norm_component(const Index3D& grid_idx, size_t axis) const;

    // --- Step Size for Ray Tracing ---
    // Originally norm[~np.isclose(norm, 0)].min() / 4
    // This needs access to all norm components.
    virtual real_t get_recommended_step_size() const;

    // --- Utility for flattened index ---
    // Public for convenience if external users need it, or protected if only internal.
    // Making it public is fine if it's a well-defined utility.
    size_t get_flat_index(const Index3D& idx) const;
    size_t get_total_nodes() const;

    // --- Pure virtual function for value access (to be implemented by derived classes) ---
    // This is tricky because scalar returns real_t, vector returns Point3D.
    // We might need separate derived classes for that or use std::variant/any if C++17+
    // For now, let's assume derived classes will have their own `get_value_at(Index3D)`

    // --- Serialization (placeholder) ---
    // Friend class for serializer to access protected/private members
    // friend class io::FieldSerializer;
    // virtual void serialize(io::Serializer& s) const = 0;
    // virtual void deserialize(io::Deserializer& d) = 0;


protected:
    // --- Core Grid Properties ---
    CoordSys m_coord_sys;
    FieldType m_field_type;
    std::array<real_t, 3> m_min_coords;      // [x0, y0, z0] or [r0, theta0, phi0]
    std::array<real_t, 3> m_node_intervals;  // [dx, dy, dz] or [dr, dtheta, dphi]
    std::array<uint_t, 3> m_npts;            // [N0, N1, N2]

    // --- Derived Grid Properties (cached for efficiency) ---
    std::array<real_t, 3> m_max_coords;
    std::array<bool, 3> m_axis_is_null;      // True if npts[axis] == 1
    std::array<bool, 3> m_axis_is_periodic;  // True for phi in spherical if full 2*PI range

    // --- Internal update methods ---
    void update_derived_properties(); // Calls all _update_... methods
    void update_max_coords();
    void update_axis_is_null();
    void update_axis_is_periodic(); // Depends on coord_sys, min_coords, max_coords, node_intervals

    // --- Internal data storage for 'norm' (gradient scaling factors) ---
    // This was N0xN1xN2x3 in Python. Storing it can be memory intensive.
    // Option 1: Store it (std::vector<real_t> flattened).
    // Option 2: Compute on the fly in get_norm_component.
    // Let's try Option 2 first to save memory, compute it when needed.
    // If performance becomes an issue, we can cache it.
    // mutable std::vector<real_t> m_norm_cache; // For Option 1
    // mutable bool m_norm_cache_valid;         // For Option 1

    // --- Validation ---
    void validate_settings() const; // Check for inconsistencies, e.g., spherical rho > 0
};

} // namespace Ckonal

#endif // Ckonal_FIELD_BASE_HPP