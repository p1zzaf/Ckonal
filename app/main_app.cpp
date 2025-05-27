#include "config_parser.hpp"
#include "vpr_model_loader.hpp"
#include "../src/fields/scalar_field.hpp"    
#include "../src/solver/eikonal_solver.hpp" 
#include "../src/io/field_serializer.hpp"   
#include "../src/core/konal_types.hpp"      // Defines Ckonal::CoordSys, Ckonal::Index3D etc.
#include "../src/core/konal_constants.hpp"  // Defines Ckonal::uint_t, Ckonal::NaN, Ckonal::INF, Ckonal::EPSILON
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept> 
#include <iomanip>   
#include <cmath>     

// Helper to print a small part of the field for verification
void print_field_slice(const Ckonal::ScalarField3D& field, const std::string& name, Ckonal::uint_t k_slice = 0, size_t max_dim = 5) { // <--- MODIFIED: Ckonal::uint_t and size_t for max_dim
    const auto& npts = field.get_npts();
    if (k_slice >= npts[2]) {
        std::cout << name << ": k_slice " << k_slice << " out of bounds." << std::endl;
        return;
    }
    std::cout << "--- " << name << " (Slice at Z-index " << k_slice 
              << ", Z-coord=" << field.get_min_coords()[2] + static_cast<double>(k_slice) * field.get_node_intervals()[2]  // <--- Added static_cast for clarity
              << ", up to " << max_dim << "x" << max_dim << ") ---" << std::endl;

    for (Ckonal::uint_t i = 0; i < std::min(npts[0], static_cast<Ckonal::uint_t>(max_dim)); ++i) { // <--- MODIFIED
        for (Ckonal::uint_t j = 0; j < std::min(npts[1], static_cast<Ckonal::uint_t>(max_dim)); ++j) { // <--- MODIFIED
            double val = field.get_value({i, j, k_slice});
            if (std::isinf(val)) {
                std::cout << std::setw(7) << "inf";
            } else if (std::isnan(val)) {
                std::cout << std::setw(7) << "nan";
            }
            else {
                std::cout << std::fixed << std::setprecision(2) << std::setw(7) << val;
            }
        }
        std::cout << std::endl;
    }
    std::cout << "--- End of " << name << " Slice ---" << std::endl;
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file.ini>" << std::endl;
        return 1;
    }
    std::string config_file_path = argv[1];

    ckonal::app::Config config; // This remains in ckonal::app namespace
    try {
        ckonal::app::ConfigParser config_parser;
        config = config_parser.parse(config_file_path);
        std::cout << "Config file '" << config_file_path << "' parsed successfully." << std::endl;
        std::cout << "  Grid X: min=" << config.grid_x_min << ", step=" << config.grid_x_step << ", count=" << config.grid_x_count << std::endl;
        std::cout << "  Velocity model VPR path: " << config.velocity_model_path << std::endl;
        std::cout << "  Source: x=" << config.source_x << ", y=" << config.source_y << ", z=" << config.source_z
                  << ", time=" << config.source_travel_time << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing config file: " << e.what() << std::endl;
        return 1;
    }

    // --- 1. Setup Velocity Field from VPR and Config ---
    Ckonal::ScalarField3D velocity_model(Ckonal::CoordSys::CARTESIAN); // <--- MODIFIED
    velocity_model.set_min_coords({config.grid_x_min, config.grid_y_min, config.grid_z_min});
    velocity_model.set_node_intervals({config.grid_x_step, config.grid_y_step, config.grid_z_step});
    velocity_model.set_npts({static_cast<Ckonal::uint_t>(config.grid_x_count),  // <--- MODIFIED
                             static_cast<Ckonal::uint_t>(config.grid_y_count), 
                             static_cast<Ckonal::uint_t>(config.grid_z_count)});

    ckonal::app::VprModelLoader vpr_loader; // This is fine, it's in ckonal::app
    double vpr_min_val, vpr_max_val;
    try {
        if (!vpr_loader.load(config.velocity_model_path, velocity_model, vpr_min_val, vpr_max_val)) {
            std::cerr << "Failed to load VPR velocity model." << std::endl;
            return 1;
        }
        std::cout << "VPR velocity model '" << config.velocity_model_path << "' loaded. VPR min/max: "
                  << vpr_min_val << "/" << vpr_max_val << std::endl;
        // print_field_slice(velocity_model, "Loaded Velocity Model", velocity_model.get_npts()[2] / 2);

    } catch (const std::exception& e) {
        std::cerr << "Error loading VPR model: " << e.what() << std::endl;
        return 1;
    }

    // --- 2. Setup EikonalSolver ---
    Ckonal::EikonalSolver solver(Ckonal::CoordSys::CARTESIAN); // <--- MODIFIED

    Ckonal::ScalarField3D& solver_velocity = solver.get_velocity_field(); // <--- MODIFIED
    solver_velocity.set_min_coords(velocity_model.get_min_coords());
    solver_velocity.set_node_intervals(velocity_model.get_node_intervals());
    solver_velocity.set_npts(velocity_model.get_npts());
    solver_velocity.set_all_values(velocity_model.get_values_raw());
    std::cout << "Solver velocity field configured." << std::endl;


    // --- 3. Initialize Solver State and Add Source ---
    try {
        solver.init_solver_state();
        std::cout << "Solver state initialized." << std::endl;

        Ckonal::Point3D source_coords_phys = {config.source_x, config.source_y, config.source_z}; // <--- MODIFIED
        Ckonal::Index3D source_idx; // <--- MODIFIED
        bool source_on_grid = true;
        for (int i=0; i<3; ++i) {
            if (solver_velocity.get_node_intervals()[i] < Ckonal::EPSILON) { // <--- MODIFIED for EPSILON
                // Assign a default index for this dimension if step is zero
                if (i==0) source_idx.i1 = 0;
                else if (i==1) source_idx.i2 = 0;
                else source_idx.i3 = 0;
                
                if (std::abs(source_coords_phys[i] - solver_velocity.get_min_coords()[i]) > Ckonal::EPSILON && solver_velocity.get_npts()[i] > 1) { // <--- MODIFIED
                    std::cerr << "Warning: Source coord " << i << " differs from min_coord for a zero-step/single-point axis." << std::endl;
                }
                continue;
            }
            double normalized_coord = (source_coords_phys[i] - solver_velocity.get_min_coords()[i]) / solver_velocity.get_node_intervals()[i];
            const auto& n_pts_dim = solver_velocity.get_npts()[i]; // This is Ckonal::uint_t

            if (normalized_coord < -Ckonal::EPSILON || normalized_coord > static_cast<double>(n_pts_dim -1) + Ckonal::EPSILON) { // <--- MODIFIED
                std::cerr << "Error: Source coordinate " << i << " (" << source_coords_phys[i]
                          << ") is outside the grid range [" << solver_velocity.get_min_coords()[i]
                          << ", " << solver_velocity.get_max_coords()[i] << "]." << std::endl;
                source_on_grid = false;
                break;
            }
            Ckonal::uint_t snapped_idx_dim = static_cast<Ckonal::uint_t>(std::round(normalized_coord)); // <--- MODIFIED
            if (snapped_idx_dim >= n_pts_dim) snapped_idx_dim = n_pts_dim - 1; 

            if (i==0) source_idx.i1 = snapped_idx_dim;
            else if (i==1) source_idx.i2 = snapped_idx_dim;
            else source_idx.i3 = snapped_idx_dim;
        }


        if (!source_on_grid) {
             std::cerr << "Source location is outside the defined grid. Cannot proceed." << std::endl;
            return 1;
        }
        
        std::cout << "Source physical coords (" << config.source_x << "," << config.source_y << "," << config.source_z
                  << ") snapped to grid index (" << source_idx.i1 << "," << source_idx.i2 << "," << source_idx.i3 << ")." << std::endl;


        if (solver.add_source_point(source_idx, config.source_travel_time)) {
            std::cout << "Source point added successfully with travel time " << config.source_travel_time << std::endl;
        } else {
            std::cerr << "Failed to add source point at grid index ("
                      << source_idx.i1 << "," << source_idx.i2 << "," << source_idx.i3 << ")." << std::endl;
            std::cerr << "Velocity at source index: " << solver_velocity.get_value(source_idx) << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during solver initialization or source addition: " << e.what() << std::endl;
        return 1;
    }

    // --- 4. Solve Eikonal Equation ---
    std::cout << "Solving Eikonal equation..." << std::endl;
    try {
        if (solver.solve()) {
            std::cout << "Eikonal equation solved successfully." << std::endl;

            const Ckonal::ScalarField3D& travel_time_result = solver.get_traveltime_field(); // <--- MODIFIED
            print_field_slice(travel_time_result, "Computed Travel Times", travel_time_result.get_npts()[2] / 2);

            Ckonal::io::BinaryFieldSerializer tt_serializer; // <--- MODIFIED
            std::string tt_output_path = "output_travel_times.pkb"; 
            try {
                tt_serializer.save_scalar_field(travel_time_result, tt_output_path);
                std::cout << "Computed travel time field saved to: " << tt_output_path << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Error saving computed travel time field: " << e.what() << std::endl;
            }

        } else {
            std::cerr << "Eikonal solver failed." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during Eikonal solve: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "--- Ckonal Application Finished ---" << std::endl;
    return 0;
}