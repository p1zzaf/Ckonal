#ifndef Ckonal_TRANSFORMATIONS_HPP
#define Ckonal_TRANSFORMATIONS_HPP

#include "core/konal_types.hpp"
#include "core/konal_constants.hpp"
#include <vector>

namespace Ckonal {
namespace transformations {

// --- Geographic to Spherical and vice-versa ---
// Geographic: lat, lon (degrees), depth (km, positive down from EARTH_RADIUS_KM)
// Spherical: r (km), theta (radians, from Z+), phi (radians, from X+ in XY plane)
Point3D geo_to_sph(const Point3D& geo_coords); // geo_coords: [lat_deg, lon_deg, depth_km]
Point3D sph_to_geo(const Point3D& sph_coords); // sph_coords: [r_km, theta_rad, phi_rad]

std::vector<Point3D> geo_to_sph_batch(const std::vector<Point3D>& geo_coords_vec);
std::vector<Point3D> sph_to_geo_batch(const std::vector<Point3D>& sph_coords_vec);


// --- Cartesian to Spherical and vice-versa ---
// Cartesian: x, y, z
// Spherical: r, theta, phi
// 'origin' is in the 'from' coordinate system and defines the origin of the 'to' system.
Point3D xyz_to_sph(const Point3D& xyz_coords,
                   const Point3D& cartesian_origin = {0.0, 0.0, 0.0},
                   bool force_phi_positive = false);

Point3D sph_to_xyz(const Point3D& sph_coords,
                   const Point3D& spherical_origin_for_new_xyz_frame = {0.0, 0.0, 0.0});
                   // spherical_origin is [r,theta,phi] of the new Cartesian frame's origin,
                   // expressed in the old spherical frame.

std::vector<Point3D> xyz_to_sph_batch(const std::vector<Point3D>& xyz_coords_vec,
                                      const Point3D& cartesian_origin = {0.0, 0.0, 0.0},
                                      bool force_phi_positive = false);

std::vector<Point3D> sph_to_xyz_batch(const std::vector<Point3D>& sph_coords_vec,
                                      const Point3D& spherical_origin_for_new_xyz_frame = {0.0, 0.0, 0.0});


// --- Spherical to Spherical ---
// Transform points from an old spherical system to a new spherical system.
// 'new_origin_in_old_sph_coords' is [r,theta,phi] of the new system's origin,
// expressed in the old spherical coordinate system.
Point3D sph_to_sph(const Point3D& old_sph_coords,
                   const Point3D& new_origin_in_old_sph_coords = {0.0, 0.0, 0.0},
                   bool force_phi_positive = false);

std::vector<Point3D> sph_to_sph_batch(const std::vector<Point3D>& old_sph_coords_vec,
                                      const Point3D& new_origin_in_old_sph_coords = {0.0, 0.0, 0.0},
                                      bool force_phi_positive = false);


// --- Cartesian to Cartesian (simple translation) ---
Point3D xyz_to_xyz(const Point3D& xyz_coords,
                   const Point3D& cartesian_origin_offset = {0.0, 0.0, 0.0});

std::vector<Point3D> xyz_to_xyz_batch(const std::vector<Point3D>& xyz_coords_vec,
                                      const Point3D& cartesian_origin_offset = {0.0, 0.0, 0.0});

// --- Rotation Matrix (if needed later) ---
// std::array<std::array<real_t, 3>, 3> get_rotation_matrix(real_t alpha_rad, real_t beta_rad, real_t gamma_rad);

// Helper: degrees to radians and vice-versa
inline real_t degrees_to_radians(real_t degrees) { return degrees * PI / 180.0; }
inline real_t radians_to_degrees(real_t radians) { return radians * 180.0 / PI; }

} // namespace transformations
} // namespace Ckonal

#endif // Ckonal_TRANSFORMATIONS_HPP