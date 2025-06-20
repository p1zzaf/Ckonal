#ifndef Ckonal_EQ_LOCATOR_HPP
#define Ckonal_EQ_LOCATOR_HPP

#include "core/konal_types.hpp"
#include "fields/scalar_field.hpp"
#include "io/field_serializer.hpp" // For loading travel times
#include <string>
#include <vector>
#include <map>
#include <memory>

// --- Optimization library placeholder ---
// #include <some_optimization_library.hpp> // e.g., NLopt, Ceres, dlib

namespace Ckonal {

// For storing arrival data: station_name/phase -> arrival_time
using ArrivalData = std::map<std::string, real_t>;

// For residual statistics (simplified from Ckonal's NormalDistribution for now)
struct ResidualStats {
    real_t mean = 0.0;
    real_t std_dev = 1.0; // Assume normal distribution for log_likelihood
    real_t log_pdf(real_t residual) const;
};
using ResidualRVMap = std::map<std::string, ResidualStats>;

using Point4D = std::array<real_t, 4>;
class EQLocator {
public:
    EQLocator(CoordSys coord_sys = CoordSys::SPHERICAL);

    void set_traveltime_directory(const std::string& dir_path); // Directory where travel time fields are stored

    void add_arrivals(const ArrivalData& new_arrivals);
    void clear_arrivals();
    const ArrivalData& get_arrivals() const;

    void add_residual_stats(const ResidualRVMap& new_stats);
    void clear_residual_stats();

    // --- Location ---
    // hypocenter_guess: [x,y,z,t0] or [r,theta,phi,t0]
    // search_deltas: defines search bounds around the guess
    // Returns located hypocenter [x,y,z,t0] and possibly uncertainty/RMS
    Point4D locate_event(const Point4D& initial_guess, const Point4D& search_deltas);

    // Calculate RMS for a given hypocenter
    real_t calculate_rms(const Point4D& hypocenter) const;

    // Calculate log-likelihood for a given model (hypocenter)
    real_t calculate_log_likelihood(const Point4D& model) const;


private:
    CoordSys m_coord_sys;
    std::string m_traveltime_dir; // Path to load travel time fields from
    ArrivalData m_arrivals;
    ResidualRVMap m_residual_stats;

    // Cache for loaded travel time fields to avoid repeated file I/O
    // Key: station_name/phase (matching arrival data keys)
    mutable std::map<std::string, std::unique_ptr<ScalarField3D>> m_traveltime_cache;

    // Helper to load a travel time field on demand
    const ScalarField3D* get_traveltime_field(const std::string& key) const;

    // --- Optimization (Simplified for now) ---
    // This would internally use an optimization algorithm.
    // For a simple version, it could be a grid search.
    Point4D perform_grid_search(const Point4D& initial_guess, const Point4D& search_deltas,
                                size_t n_steps_spatial, size_t n_steps_time) const;
};

// Using std::array for 4D points (x,y,z,t or r,th,ph,t)



} // namespace Ckonal

#endif // Ckonal_EQ_LOCATOR_HPP