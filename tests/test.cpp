#include <iostream>
#include <vector>
#include <string>
#include <cmath> // For std::isinf, std::abs
#include "core/konal_constants.hpp"
#include "core/konal_types.hpp"
#include "fields/scalar_field.hpp"
#include "solver/eikonal_solver.hpp"
#include "io/field_serializer.hpp"

// Helper to print a small field
void print_travel_times(const Ckonal::ScalarField3D& tt_field, size_t max_print_dim = 5) {
    const auto& npts = tt_field.get_npts();
    std::cout << "Travel Time Field (up to " << max_print_dim << "x" << max_print_dim << " slice at k=0):\n";
for (size_t i = 0; i < std::min(static_cast<size_t>(npts[0]), max_print_dim); ++i) {
for (size_t j = 0; j < std::min(static_cast<size_t>(npts[1]), max_print_dim); ++j) {
            // Print for k=0 slice
            if (npts[2] > 0) {
                Ckonal::real_t val = tt_field.get_value(i, j, 0);
                if (std::isinf(val)) {
                    std::cout << " inf   ";
                } else {
                    printf("%6.2f ", val);
                }
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


int main() {
    std::cout << "--- Ckonal C++ Eikonal Solver Demo ---" << std::endl;

    // --- 1. Setup EikonalSolver and Velocity Field ---
    Ckonal::CoordSys coord_system = Ckonal::CoordSys::CARTESIAN;
    Ckonal::EikonalSolver solver(coord_system);

    // Define velocity field properties
    Ckonal::ScalarField3D& velocity_model = solver.get_velocity_field();
    velocity_model.set_min_coords({0.0, 0.0, 0.0});
    velocity_model.set_node_intervals({1.0, 1.0, 1.0}); // Grid spacing of 1 unit
    velocity_model.set_npts({11, 11, 11}); // A small 11x11x11 grid

    // Set velocity values (e.g., constant velocity)
    Ckonal::real_t constant_velocity = 2.0; // units/second
    velocity_model.fill(constant_velocity);
    std::cout << "Velocity model configured with constant velocity: " << constant_velocity << std::endl;

    // --- 2. Initialize Solver State and Add Source ---
    try {
        solver.init_solver_state(); // Initializes travel time field, FMM states, heap
        std::cout << "Solver state initialized." << std::endl;

        Ckonal::Index3D source_location(5, 5, 5); // Source at the center of the grid
        // Ckonal::Index3D source_location(0, 0, 0); // Source at a corner
        if (solver.add_source_point(source_location, 0.0)) {
            std::cout << "Source point added at grid index ("
                      << source_location.i1 << "," << source_location.i2 << "," << source_location.i3 << ")." << std::endl;
        } else {
            std::cerr << "Failed to add source point." << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during solver initialization or source addition: " << e.what() << std::endl;
        return 1;
    }

    // --- 3. Solve Eikonal Equation ---
    std::cout << "Solving Eikonal equation..." << std::endl;
    try {
        if (solver.solve()) {
            std::cout << "Eikonal equation solved successfully." << std::endl;

            // --- 4. Access and Print Results ---
            const Ckonal::ScalarField3D& travel_time_result = solver.get_traveltime_field();
            print_travel_times(travel_time_result);

            // Example: Get travel time at a specific point
            Ckonal::Index3D query_idx(6, 5, 5); // One step away from source (5,5,5) in i1
            Ckonal::real_t tt_at_query = travel_time_result.get_value(query_idx);
            std::cout << "Travel time at (" << query_idx.i1 << "," << query_idx.i2 << "," << query_idx.i3 << "): " << tt_at_query << std::endl;
            // Expected: roughly 1.0 / 2.0 = 0.5 for this simple case if interval is 1.0

            Ckonal::Point3D query_coord = travel_time_result.get_node_coords(query_idx);
            Ckonal::real_t tt_interpolated = travel_time_result.interpolate_value({query_coord[0] + 0.5, query_coord[1], query_coord[2]});
             std::cout << "Interpolated TT at (" << query_coord[0] + 0.5 << "," << query_coord[1] << "," << query_coord[2] << "): " << tt_interpolated << std::endl;


            // --- 5. Save Travel Time Field (Optional) ---
            Ckonal::io::BinaryFieldSerializer serializer;
            std::string tt_output_path = "travel_times.pkb";
            try {
                serializer.save_scalar_field(travel_time_result, tt_output_path);
                std::cout << "Travel time field saved to: " << tt_output_path << std::endl;

                // --- 6. Load Travel Time Field (Test) ---
                std::unique_ptr<Ckonal::ScalarField3D> loaded_tt_field = serializer.load_scalar_field(tt_output_path);
                if (loaded_tt_field) {
                    std::cout << "Travel time field loaded successfully from: " << tt_output_path << std::endl;
                    // print_travel_times(*loaded_tt_field); // Verify
                }
            } catch (const std::exception& e) {
                std::cerr << "Error during field serialization/deserialization: " << e.what() << std::endl;
            }

        } else {
            std::cerr << "Eikonal solver failed or had no sources." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during Eikonal solve: " << e.what() << std::endl;
        return 1;
    }
    
    // In src/main.cpp
// ... (EikonalSolver demo) ...

    // --- PointSourceSolver Demo (Illustrative - requires more setup) ---
    std::cout << "\n--- PointSourceSolver Demo ---" << std::endl;
    try {
        Ckonal::PointSourceSolver pssolver(Ckonal::CoordSys::CARTESIAN); // Far field is Cartesian

        // Far-field velocity setup
        Ckonal::ScalarField3D& far_vv_pss = pssolver.get_velocity_field();
        far_vv_pss.set_min_coords({0.0, 0.0, 0.0});       // km
        far_vv_pss.set_node_intervals({10.0, 10.0, 10.0}); // 10 km spacing
        far_vv_pss.set_npts({21, 21, 11});               // Grid: 0-200km (X), 0-200km (Y), 0-100km (Z)
        far_vv_pss.fill(6.0);                            // Constant velocity 6 km/s
        std::cout << "PSSolver: Far-field velocity configured." << std::endl;

        // Source location (in far-field Cartesian coordinates)
        // Make sure it's well within the far-field grid
        pssolver.set_source_location({100.0, 100.0, 50.0}); // Center of the far-field grid

        // Near-field parameters (can use defaults or set them)
        // Let's use defaults by not calling setters, or explicitly set:
        pssolver.set_near_field_radial_nodes(32);    // Default was 64
        pssolver.set_near_field_polar_nodes(17);     // Default was 33
        pssolver.set_near_field_azimuthal_nodes(32); // Default was 64
        // pssolver.set_near_field_radial_interval(0.5); // Default will be calculated based on far-field interval

        std::cout << "PSSolver: Calling init_solver_state..." << std::endl;
        pssolver.init_solver_state(); 
        std::cout << "PSSolver: init_solver_state finished." << std::endl;
        
        std::cout << "PSSolver: Calling solve()..." << std::endl;
        if (pssolver.solve()) {
            std::cout << "PointSourceSolver solved successfully." << std::endl;
            // You can print parts of pssolver.get_traveltime_field() here
            std::cout << "PointSourceSolver: Travel times in far-field (slice at Z index "
          << pssolver.get_traveltime_field().get_npts()[2] / 2 << "):" << std::endl;
print_travel_times(pssolver.get_traveltime_field(), 7); // Print a 7x7 slice from center Z

// Example query point near source in far-field
Ckonal::Index3D far_query_idx(
    pssolver.get_traveltime_field().get_npts()[0] / 2 + 1, // e.g., 100+10 = 110km if npts[0]/2 maps to 100km
    pssolver.get_traveltime_field().get_npts()[1] / 2,     // 100km
    pssolver.get_traveltime_field().get_npts()[2] / 2      // 50km
);
// Ensure far_query_idx is valid before using
if (far_query_idx.i1 < pssolver.get_traveltime_field().get_npts()[0] &&
    far_query_idx.i2 < pssolver.get_traveltime_field().get_npts()[1] &&
    far_query_idx.i3 < pssolver.get_traveltime_field().get_npts()[2]) {
    Ckonal::Point3D far_query_coords_actual = pssolver.get_traveltime_field().get_node_coords(far_query_idx);
    Ckonal::real_t tt_far_query = pssolver.get_traveltime_field().get_value(far_query_idx);
    std::cout << "PSSolver: Far-field TT at far_idx(" << far_query_idx.i1 << "," << far_query_idx.i2 << "," << far_query_idx.i3 << ")"
              << " corresponding to coords (" << far_query_coords_actual[0] << "," << far_query_coords_actual[1] << "," << far_query_coords_actual[2] << ")"
              << " is " << tt_far_query << std::endl;
}
        } else {
            std::cerr << "PointSourceSolver main: solve() returned false." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error in PointSourceSolver demo: " << e.what() << std::endl;
    }

    std::cout << "\nPress Enter to continue...";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
    std::cin.get(); 
    return 0;
}