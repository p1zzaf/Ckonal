#include "transformations.hpp"
#include <cmath>    // For sin, cos, acos, atan2, sqrt, fmod
#include <stdexcept> // For potential errors if inputs are out of expected range

namespace Ckonal {
namespace transformations {

// --- Geographic to Spherical ---
Point3D geo_to_sph(const Point3D& geo_coords) {
    // geo_coords: [lat_deg, lon_deg, depth_km]
    // sph_coords: [r_km, theta_rad, phi_rad]
    Point3D sph_coords;
    sph_coords[0] = EARTH_RADIUS_KM - geo_coords[2]; // r = R_earth - depth
    sph_coords[1] = PI / 2.0 - degrees_to_radians(geo_coords[0]); // theta = pi/2 - lat_rad
    sph_coords[2] = degrees_to_radians(geo_coords[1]); // phi = lon_rad
    return sph_coords;
}

std::vector<Point3D> geo_to_sph_batch(const std::vector<Point3D>& geo_coords_vec) {
    std::vector<Point3D> sph_coords_vec;
    sph_coords_vec.reserve(geo_coords_vec.size());
    for (const auto& geo : geo_coords_vec) {
        sph_coords_vec.push_back(geo_to_sph(geo));
    }
    return sph_coords_vec;
}

// --- Spherical to Geographic ---
Point3D sph_to_geo(const Point3D& sph_coords) {
    // sph_coords: [r_km, theta_rad, phi_rad]
    // geo_coords: [lat_deg, lon_deg, depth_km]
    Point3D geo_coords;
    geo_coords[0] = radians_to_degrees(PI / 2.0 - sph_coords[1]); // lat_deg = 90 - theta_deg
    geo_coords[1] = radians_to_degrees(sph_coords[2]); // lon_deg = phi_deg
    geo_coords[2] = EARTH_RADIUS_KM - sph_coords[0]; // depth = R_earth - r
    return geo_coords;
}

std::vector<Point3D> sph_to_geo_batch(const std::vector<Point3D>& sph_coords_vec) {
    std::vector<Point3D> geo_coords_vec;
    geo_coords_vec.reserve(sph_coords_vec.size());
    for (const auto& sph : sph_coords_vec) {
        geo_coords_vec.push_back(sph_to_geo(sph));
    }
    return geo_coords_vec;
}


// --- Cartesian to Spherical ---
Point3D xyz_to_sph(const Point3D& xyz_coords, const Point3D& cartesian_origin, bool force_phi_positive) {
    Point3D rel_xyz = {
        xyz_coords[0] - cartesian_origin[0],
        xyz_coords[1] - cartesian_origin[1],
        xyz_coords[2] - cartesian_origin[2]
    };

    real_t r = std::sqrt(rel_xyz[0] * rel_xyz[0] + rel_xyz[1] * rel_xyz[1] + rel_xyz[2] * rel_xyz[2]);
    real_t theta, phi;

    if (r < 1e-9) { // At origin
        theta = 0.0;
        phi = 0.0;
    } else {
        theta = std::acos(rel_xyz[2] / r); // [0, PI]
        phi = std::atan2(rel_xyz[1], rel_xyz[0]); // [-PI, PI]
    }

    if (force_phi_positive) {
        if (phi < 0.0) {
            phi += TWO_PI; // range [0, 2PI)
        }
    }
    // Default atan2 range is [-PI, PI]. Ckonal's sph2sph has specific phi adjustments.
    // The python code for xyz2sph is simpler for phi.
    // `pp = np.mod(pp, 2*np.pi)` if force_phi_positive.
    // `np.mod(pp, -np.pi, pp, where=pp>np.pi)` otherwise.
    // This means: if not force_phi_positive, map (PI, 2PI) to (-PI, 0).
    // If force_phi_positive is false (default in Python for xyz2sph) it seems phi stays in [-PI, PI].
    // Let's match the Python's sph2sph default for phi:
    // if not force_phi_positive and phi > PI (e.g. from atan2 slightly > PI for negative x, tiny negative y)
    // then map it to the equivalent angle in [-PI,PI).
    // atan2 already gives [-PI, PI]. If force_phi_positive, then map to [0, 2PI).

    return {r, theta, phi};
}

std::vector<Point3D> xyz_to_sph_batch(const std::vector<Point3D>& xyz_coords_vec,
                                      const Point3D& cartesian_origin, bool force_phi_positive) {
    std::vector<Point3D> sph_coords_vec;
    sph_coords_vec.reserve(xyz_coords_vec.size());
    for (const auto& xyz : xyz_coords_vec) {
        sph_coords_vec.push_back(xyz_to_sph(xyz, cartesian_origin, force_phi_positive));
    }
    return sph_coords_vec;
}


// --- Spherical to Cartesian ---
Point3D sph_to_xyz(const Point3D& sph_coords, const Point3D& spherical_origin_for_new_xyz_frame) {
    real_t r_node = sph_coords[0];
    real_t theta_node = sph_coords[1];
    real_t phi_node = sph_coords[2];

    // Convert node to Cartesian in its own (old) frame
    real_t xx_node = r_node * std::sin(theta_node) * std::cos(phi_node);
    real_t yy_node = r_node * std::sin(theta_node) * std::sin(phi_node);
    real_t zz_node = r_node * std::cos(theta_node);

    // Convert the 'spherical_origin_for_new_xyz_frame' (which is given in the old spherical frame)
    // to Cartesian coordinates in the old frame. This Cartesian point is where the new XYZ frame's origin lies.
    real_t r_origin = spherical_origin_for_new_xyz_frame[0];
    real_t theta_origin = spherical_origin_for_new_xyz_frame[1];
    real_t phi_origin = spherical_origin_for_new_xyz_frame[2];

    real_t x0 = r_origin * std::sin(theta_origin) * std::cos(phi_origin);
    real_t y0 = r_origin * std::sin(theta_origin) * std::sin(phi_origin);
    real_t z0 = r_origin * std::cos(theta_origin);
    
    // The result is the node's Cartesian coordinates relative to the new origin
    return {xx_node - x0, yy_node - y0, zz_node - z0};
}

std::vector<Point3D> sph_to_xyz_batch(const std::vector<Point3D>& sph_coords_vec,
                                      const Point3D& spherical_origin_for_new_xyz_frame) {
    std::vector<Point3D> xyz_coords_vec;
    xyz_coords_vec.reserve(sph_coords_vec.size());
    for (const auto& sph : sph_coords_vec) {
        xyz_coords_vec.push_back(sph_to_xyz(sph, spherical_origin_for_new_xyz_frame));
    }
    return xyz_coords_vec;
}

// --- Spherical to Spherical ---
Point3D sph_to_sph(const Point3D& old_sph_coords, const Point3D& new_origin_in_old_sph_coords, bool force_phi_positive) {
    // Strategy: old_sph -> old_xyz -> translate origin -> new_xyz -> new_sph
    // 1. Convert old_sph_coords to Cartesian (in the old frame, origin at 0,0,0 of old frame)
    real_t r_old = old_sph_coords[0];
    real_t theta_old = old_sph_coords[1];
    real_t phi_old = old_sph_coords[2];
    Point3D old_xyz_node = {
        r_old * std::sin(theta_old) * std::cos(phi_old),
        r_old * std::sin(theta_old) * std::sin(phi_old),
        r_old * std::cos(theta_old)
    };

    // 2. Convert new_origin_in_old_sph_coords to Cartesian (this is where the new frame's origin lies in old_xyz)
    real_t r_origin_new = new_origin_in_old_sph_coords[0];
    real_t theta_origin_new = new_origin_in_old_sph_coords[1];
    real_t phi_origin_new = new_origin_in_old_sph_coords[2];
    Point3D new_origin_xyz_in_old_frame = {
        r_origin_new * std::sin(theta_origin_new) * std::cos(phi_origin_new),
        r_origin_new * std::sin(theta_origin_new) * std::sin(phi_origin_new),
        r_origin_new * std::cos(theta_origin_new)
    };

    // 3. Translate: get node's coordinates relative to the new origin (still in old_xyz alignment)
    Point3D node_xyz_in_new_frame_alignment = {
        old_xyz_node[0] - new_origin_xyz_in_old_frame[0],
        old_xyz_node[1] - new_origin_xyz_in_old_frame[1],
        old_xyz_node[2] - new_origin_xyz_in_old_frame[2]
    };

    // 4. Convert these translated XYZ coordinates to spherical (this will be the new_sph_coords)
    //    Here, the 'cartesian_origin' for xyz_to_sph is {0,0,0} because we've already translated.
    Point3D new_sph_coords = xyz_to_sph(node_xyz_in_new_frame_alignment, {0.0, 0.0, 0.0}, force_phi_positive);

    // Python's sph2sph has a specific phi adjustment:
    // pp = np.mod(pp, 2*np.pi)
    // if force_phi_positive is False:
    //     np.mod(pp, -np.pi, pp, where=pp>np.pi) -> maps (PI, 2PI) to (-PI, 0)
    // Our xyz_to_sph with force_phi_positive=false should already give phi in [-PI, PI].
    // If force_phi_positive=true, it gives [0, 2PI).
    // The Python code's logic:
    //   pp = np.arctan2(xyz[...,1], xyz[...,0]) -> pp in [-PI, PI]
    //   pp = np.mod(pp, 2*np.pi) -> pp in [0, 2PI) IF original pp was negative, else no change.
    //                              e.g. -PI/2 becomes 3PI/2. PI stays PI.
    //   if force_phi_positive is False:
    //      np.mod(pp, -np.pi, pp, where=pp>np.pi) -> if pp is (PI, 2PI), then pp = pp - 2PI. (incorrect interpretation)
    // Correct interpretation of np.mod(x, N, out=arr, where=cond):
    //  `np.mod(pp, 2*np.pi)`: This ensures `pp` is in `[0, 2*np.pi)` if `2*np.pi` is positive.
    //                         If pp was -pi/2, it becomes 3pi/2.
    //  `np.mod(pp, -np.pi, pp, where=pp>np.pi)`: This is unusual. `numpy.fmod` is closer to C's `fmod`.
    //                                            `np.mod(a, n)` has sign of divisor `n`.
    //  The intent of `np.mod(pp, -np.pi, pp, where=pp>np.pi)` seems to be a typo or complex way to say:
    //  If `pp > np.pi` (after the first `np.mod(pp, 2*np.pi)` means `pp` is in `(np.pi, 2*np.pi)`),
    //  then set `pp = pp - 2*np.pi` to bring it to `(-np.pi, 0)`.
    //  So, if `force_phi_positive` is false, final range is `[-np.pi, np.pi]`.
    //  Our `xyz_to_sph` already does this if `force_phi_positive` is false for its output `phi`.

    // Let's re-check phi logic from python `sph2sph` directly:
    // pp  = np.arctan2(xyz[...,1], xyz[...,0])  // pp in [-PI, PI]
    // pp = np.mod(pp, 2*np.pi)                 // pp in [0, 2*PI) (e.g. -0.1 becomes ~1.9*PI)
    // if force_phi_positive is False:
    //     np.mod(pp, -np.pi, pp, where=pp>np.pi) // This is trying to map (PI, 2PI) to (-PI, 0)
    //                                           // but np.mod(val, -N) is (val % -N), which is weird.
    //                                           // A more direct way: if pp > PI then pp -= TWO_PI
    // So, the target logic for phi seems to be:
    // 1. Calculate initial phi (e.g. via atan2, result in [-PI, PI]).
    // 2. If force_phi_positive: map to [0, 2PI).
    // 3. If NOT force_phi_positive: map to [-PI, PI]. (atan2 already does this).
    // The python `sph2sph` `pp = np.mod(pp, 2*np.pi)` line is the confusing one if not `force_phi_positive`.
    // Given `xyz_to_sph`'s phi handling, this should be correct.

    return new_sph_coords;
}


std::vector<Point3D> sph_to_sph_batch(const std::vector<Point3D>& old_sph_coords_vec,
                                      const Point3D& new_origin_in_old_sph_coords, bool force_phi_positive) {
    std::vector<Point3D> new_sph_coords_vec;
    new_sph_coords_vec.reserve(old_sph_coords_vec.size());
    for (const auto& old_sph : old_sph_coords_vec) {
        new_sph_coords_vec.push_back(sph_to_sph(old_sph, new_origin_in_old_sph_coords, force_phi_positive));
    }
    return new_sph_coords_vec;
}


// --- Cartesian to Cartesian ---
Point3D xyz_to_xyz(const Point3D& xyz_coords, const Point3D& cartesian_origin_offset) {
    return {
        xyz_coords[0] - cartesian_origin_offset[0],
        xyz_coords[1] - cartesian_origin_offset[1],
        xyz_coords[2] - cartesian_origin_offset[2]
    };
}

std::vector<Point3D> xyz_to_xyz_batch(const std::vector<Point3D>& xyz_coords_vec,
                                      const Point3D& cartesian_origin_offset) {
    std::vector<Point3D> result_vec;
    result_vec.reserve(xyz_coords_vec.size());
    for (const auto& xyz : xyz_coords_vec) {
        result_vec.push_back(xyz_to_xyz(xyz, cartesian_origin_offset));
    }
    return result_vec;
}

/* Rotation matrix (if needed)
std::array<std::array<real_t, 3>, 3> get_rotation_matrix(real_t alpha_rad, real_t beta_rad, real_t gamma_rad) {
    // ... implementation of Z-Y'-Z'' Euler rotation ...
    return rotation_matrix_data;
}
*/

} // namespace transformations
} // namespace Ckonal