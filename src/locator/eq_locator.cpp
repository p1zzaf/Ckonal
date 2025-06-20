#include "eq_locator.hpp"
#include "utils/math_utils.hpp" // For magnitude etc.
#include <cmath>     // For std::sqrt
#include <stdexcept> // For errors
#include <iostream>  // For debug/logging

namespace Ckonal {

real_t ResidualStats::log_pdf(real_t residual) const {
    if (std_dev <= 1e-9) return (std::abs(residual - mean) < 1e-9) ? 0.0 : NEG_INF; // Delta function
    real_t z = (residual - mean) / std_dev;
    return -0.5 * z * z - std::log(std_dev * std::sqrt(TWO_PI));
}


EQLocator::EQLocator(CoordSys coord_sys) : m_coord_sys(coord_sys) {}

void EQLocator::set_traveltime_directory(const std::string& dir_path) {
    m_traveltime_dir = dir_path;
    // Basic check: does directory exist? (platform dependent or use filesystem library)
    // For now, assume valid path.
    m_traveltime_cache.clear(); // Clear cache if directory changes
}

void EQLocator::add_arrivals(const ArrivalData& new_arrivals) {
    for (const auto& pair : new_arrivals) {
        m_arrivals[pair.first] = pair.second;
    }
}
void EQLocator::clear_arrivals() { m_arrivals.clear(); }
const ArrivalData& EQLocator::get_arrivals() const { return m_arrivals; }

void EQLocator::add_residual_stats(const ResidualRVMap& new_stats) {
     for (const auto& pair : new_stats) {
        m_residual_stats[pair.first] = pair.second;
    }
}
void EQLocator::clear_residual_stats() { m_residual_stats.clear(); }


const ScalarField3D* EQLocator::get_traveltime_field(const std::string& key) const {
    auto it = m_traveltime_cache.find(key);
    if (it != m_traveltime_cache.end()) {
        return it->second.get();
    }

    if (m_traveltime_dir.empty()) {
        // std::cerr << "Warning: Travel time directory not set for EQLocator." << std::endl;
        return nullptr;
    }

    // Construct file path. Assume key is "NET.STA.PHASE" and files are named accordingly, e.g., "NET.STA.PHASE.bin"
    // This naming convention needs to be consistent with how fields are saved.
    std::string filename = key + ".pkb"; // Ckonal Binary
    std::string full_path = m_traveltime_dir + "/" + filename; // Add path separator logic for cross-platform

    try {
        io::BinaryFieldSerializer loader;
        std::unique_ptr<ScalarField3D> tt_field = loader.load_scalar_field(full_path);
        if (tt_field) {
            // Check if tt_field's coordinate system matches m_coord_sys (or handle transform)
            if (tt_field->get_coord_sys() != m_coord_sys) {
                std::cerr << "Warning: Loaded travel time field for " << key
                          << " has different coordinate system than locator. Transformation needed but not yet implemented."
                          << std::endl;
                // For now, don't cache if mismatch, or transform if possible.
                return nullptr; 
            }
            auto* ptr = tt_field.get();
            m_traveltime_cache[key] = std::move(tt_field);
            return ptr;
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to load travel time field for " << key << " from " << full_path << ": " << e.what() << std::endl;
    }
    return nullptr;
}

real_t EQLocator::calculate_rms(const Point4D& hypocenter) const {
    if (m_arrivals.empty()) return 0.0;

    Point3D location = {hypocenter[0], hypocenter[1], hypocenter[2]};
    real_t origin_time_t0 = hypocenter[3];
    real_t sum_sq_residuals = 0.0;
    int valid_arrivals = 0;

    for (const auto& arrival_pair : m_arrivals) {
        const std::string& key = arrival_pair.first;         // e.g., "STA1.P"
        real_t observed_arrival_time = arrival_pair.second;

        const ScalarField3D* tt_field = get_traveltime_field(key);
        if (!tt_field) {
            // std::cerr << "Warning: No travel time field for " << key << ", skipping for RMS." << std::endl;
            continue; // Skip if travel time field is not available
        }

        // Interpolate travel time from the field at the given location
        // Null value for interpolation: if point is outside, tt is effectively infinite
        real_t predicted_travel_time = tt_field->interpolate_value(location, INF);

        if (std::isinf(predicted_travel_time)) {
            // std::cerr << "Warning: Predicted travel time is INF for " << key << " at location, skipping." << std::endl;
            continue; // Cannot calculate residual if predicted TT is invalid
        }
        
        real_t residual = observed_arrival_time - (origin_time_t0 + predicted_travel_time);
        sum_sq_residuals += residual * residual;
        valid_arrivals++;
    }

    if (valid_arrivals == 0) return INF; // Or some other indicator of failure
    return std::sqrt(sum_sq_residuals / static_cast<real_t>(valid_arrivals));
}


real_t EQLocator::calculate_log_likelihood(const Point4D& model) const {
    if (m_arrivals.empty() || m_residual_stats.empty()) return NEG_INF;

    Point3D location = {model[0], model[1], model[2]};
    real_t origin_time_t0 = model[3];
    real_t total_log_likelihood = 0.0;
    bool_t at_least_one_valid = false;

    for (const auto& arrival_pair : m_arrivals) {
        const std::string& key = arrival_pair.first;
        real_t observed_arrival_time = arrival_pair.second;

        auto stats_it = m_residual_stats.find(key);
        if (stats_it == m_residual_stats.end()) {
            // std::cerr << "Warning: No residual stats for " << key << ", skipping for log-likelihood." << std::endl;
            continue;
        }
        const ResidualStats& stats = stats_it->second;

        const ScalarField3D* tt_field = get_traveltime_field(key);
        if (!tt_field) {
            // std::cerr << "Warning: No travel time field for " << key << ", skipping for log-likelihood." << std::endl;
            continue;
        }

        real_t predicted_travel_time = tt_field->interpolate_value(location, INF);
        if (std::isinf(predicted_travel_time)) {
            total_log_likelihood += NEG_INF; // Or a very large negative number
            at_least_one_valid = true; // Technically counts, but with zero probability
            continue;
        }

        real_t residual = observed_arrival_time - (origin_time_t0 + predicted_travel_time);
        total_log_likelihood += stats.log_pdf(residual);
        at_least_one_valid = true;
    }
    return at_least_one_valid ? total_log_likelihood : NEG_INF; // If no arrivals processed
}


Point4D EQLocator::perform_grid_search(const Point4D& initial_guess, const Point4D& search_deltas,
                                     size_t n_steps_spatial, size_t n_steps_time) const {
    Point4D best_hypocenter = initial_guess;
    real_t min_rms = calculate_rms(initial_guess);
    if (std::isinf(min_rms) && !m_arrivals.empty()) { // Initial guess might be bad
        min_rms = INF; // Ensure it can be improved
    }


    std::array<real_t, 4> step_sizes;
    std::array<size_t, 4> n_steps_dim = {n_steps_spatial, n_steps_spatial, n_steps_spatial, n_steps_time};

    for (int i = 0; i < 4; ++i) {
        step_sizes[i] = (n_steps_dim[i] <= 1) ? 0.0 : (2.0 * search_deltas[i]) / static_cast<real_t>(n_steps_dim[i] -1);
    }

    Point4D current_hypo;
    for (size_t i0 = 0; i0 < n_steps_dim[0]; ++i0) {
        current_hypo[0] = (initial_guess[0] - search_deltas[0]) + static_cast<real_t>(i0) * step_sizes[0];
        for (size_t i1 = 0; i1 < n_steps_dim[1]; ++i1) {
            current_hypo[1] = (initial_guess[1] - search_deltas[1]) + static_cast<real_t>(i1) * step_sizes[1];
            for (size_t i2 = 0; i2 < n_steps_dim[2]; ++i2) {
                current_hypo[2] = (initial_guess[2] - search_deltas[2]) + static_cast<real_t>(i2) * step_sizes[2];
                for (size_t i3 = 0; i3 < n_steps_dim[3]; ++i3) {
                    current_hypo[3] = (initial_guess[3] - search_deltas[3]) + static_cast<real_t>(i3) * step_sizes[3];
                    
                    real_t current_rms = calculate_rms(current_hypo);
                    if (current_rms < min_rms) {
                        min_rms = current_rms;
                        best_hypocenter = current_hypo;
                    }
                }
            }
        }
    }
    // std::cout << "Grid search best RMS: " << min_rms << std::endl;
    return best_hypocenter;
}


// In src/locator/eq_locator.cpp
Point4D EQLocator::locate_event(const Point4D& initial_guess, const Point4D& search_deltas) {
    if (m_arrivals.empty()) {
        throw std::runtime_error("No arrival data to perform location.");
    }
    if (m_traveltime_dir.empty() && m_traveltime_cache.empty()){
        throw std::runtime_error("Travel time data source not specified or empty.");
    }

    // --- TODO: Replace with a proper optimization algorithm ---
    // For now, use the simple grid search.
    // Example: 11 steps for spatial, 5 for time.
    size_t n_spatial_steps_for_grid = 11; 
    size_t n_time_steps_for_grid = 5;

    // Ensure steps are at least 1 if they were configurable and could be zero
    // For hardcoded values like above, this isn't strictly necessary unless they might change to 0.
    if (n_spatial_steps_for_grid == 0) n_spatial_steps_for_grid = 1;
    if (n_time_steps_for_grid == 0) n_time_steps_for_grid = 1;


    Point4D located_hypocenter = perform_grid_search(initial_guess, search_deltas, n_spatial_steps_for_grid, n_time_steps_for_grid);
    
    // If using differential_evolution or similar, the call would be here:
    // auto objective_func = [&](const std::vector<real_t>& params) {
    //     Point4D hypo = {params[0], params[1], params[2], params[3]};
    //     return calculate_rms(hypo); // Or -log_likelihood for MLE
    // };
    // std::vector<std::pair<real_t, real_t>> bounds(4);
    // for(int i=0; i<4; ++i) bounds[i] = {initial_guess[i] - search_deltas[i], initial_guess[i] + search_deltas[i]};
    // auto result = some_optimization_library::minimize(objective_func, bounds, ...);
    // located_hypocenter = {result.params[0], ...};

    return located_hypocenter;
}


} // namespace Ckonal