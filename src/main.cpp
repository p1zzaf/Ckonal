#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath> // For std::isinf, std::abs
#include "core/konal_constants.hpp"
#include "core/konal_types.hpp"
#include "fields/scalar_field.hpp"
#include "solver/eikonal_solver.hpp"
#include "io/field_serializer.hpp"
#include "io/parmset.h"
#include "io/filecheck.h"
#include "io/vprmodelio.h"
#include "io/geom.h"
#include "io/npy.h"

int main(int argc, char* argv[]) {
    std::cout << "--- Ckonal C++ Eikonal Solver for 3D model under cartesian coordinate ---" << std::endl;

    // run mode
    int mode_run = 0;
    int frunflag = -1;
    // windows
    const string path_seperator("\\");

     // 检查参数数量
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input par file>" << std::endl;
        return 1; // 返回非零值表示错误
    }

    // 获取字符串参数
    std::string parfile = argv[1];

    // check file exist
    if(!fileExists(parfile)){
        std::cerr << "PAR file not exist!!   " << parfile << std::endl;
        return 1;
    }

    // parse PAR FILE
    ParmSet ps;
    std::ifstream parms_ifs(parfile);
    parms_ifs >> ps;
    parms_ifs.close();

    // get running mode
    if(ps.isDef("CKONAL_RUN_MODE")){mode_run = ps.getInt("CKONAL_RUN_MODE");}
    if(mode_run != 0 || mode_run != 1){
        mode_run = 0;
        std::cout << "!! wrong CKONAL_RUN_MODE in PAR file: " << mode_run << "  we set it to 0!" << std::endl;
    }

    std::cout << "check PAR file for Ckonal [mode: " << mode_run << " ] ..." << std::endl;
    // mode 0:
    if(mode_run == 0){
        if(!ps.isDef("IN_3D_VELOCITY_MODEL")){std::cout << "not define IN_3D_VELOCITY_MODEL!!" << std::endl;frunflag = -1;}
        if(!ps.isDef("CKONAL_SOURCE_X")){std::cout << "not define CKONAL_SOURCE_X!!" << std::endl;frunflag = -1;}
        if(!ps.isDef("CKONAL_SOURCE_Y")){std::cout << "not define CKONAL_SOURCE_Y!!" << std::endl;frunflag = -1;}
        if(!ps.isDef("CKONAL_SOURCE_Z")){std::cout << "not define CKONAL_SOURCE_Z!!" << std::endl;frunflag = -1;}
        if(!ps.isDef("OUT_PATH_TRAVEL_TIME")){std::cout << "not define OUT_PATH_TRAVEL_TIME!!" << std::endl;frunflag = -1;}
        if(frunflag == -1){
            std::cerr << "!!!! para missing in PAR file: " << parfile << std::endl;
            return -1;
        }
    }
    // mode 1:
    if(mode_run == 1){
        if(!ps.isDef("IN_3D_VELOCITY_MODEL")){std::cout << "not define IN_3D_VELOCITY_MODEL!!" << std::endl;frunflag = -1;}
        if(!ps.isDef("GEOMETRY_FILE")){std::cout << "not define GEOMETRY_FILE!!" << std::endl;frunflag = -1;}
        if(!ps.isDef("GEOMETRY_CHANGE_XY")){std::cout << "not define GEOMETRY_CHANGE_XY!!" << std::endl;frunflag = -1;}
        if(!ps.isDef("OUT_PATH_TRAVEL_TIME")){std::cout << "not define OUT_PATH_TRAVEL_TIME!!" << std::endl;frunflag = -1;}
        if(frunflag == -1){
            std::cerr << "!!!! para missing in PAR file: " << parfile << std::endl;
            return -1;
        }
    }
    
    // get 3D velocity model
    std::string velocity3D_file = ps.getString("IN_3D_VELOCITY_MODEL");

    // get output path of travel time field
    std::string out_path_traveltime = ps.getString("OUT_PATH_TRAVEL_TIME");

    // read 3D velocity model data to vector
    std::vector<Ckonal::real_t> vmodel_values;
    std::vector<Ckonal::real_t> modelcoords_min;
    std::vector<Ckonal::real_t> modelcoords_step;
    std::vector<Ckonal::uint_t>  modelcoords_num;
    frunflag = readvpr(vmodel_values, modelcoords_min, modelcoords_num, modelcoords_step, velocity3D_file);
    if(frunflag < 0){
        std::cerr << "load 3D model file (vpr format) error!!" << std::endl;
        return -2;
    }

    // --- 1. Setup EikonalSolver and Velocity Field ---
    std::cout << "set coordinate of Ckonal ..." << std::endl;
    Ckonal::CoordSys coord_system = Ckonal::CoordSys::CARTESIAN;
    Ckonal::EikonalSolver solver(coord_system);

    std::cout << "set grids of Ckonal ..." << std::endl;
    // Define velocity field properties
    Ckonal::ScalarField3D& velocity_model = solver.get_velocity_field();
    // set min of coords: x y z
    velocity_model.set_min_coords({modelcoords_min[0], modelcoords_min[1], modelcoords_min[2]});
    // set model step: x y z
    velocity_model.set_node_intervals({modelcoords_step[0], modelcoords_step[1], modelcoords_step[2]});
    // set coords num: x y z
    velocity_model.set_npts({modelcoords_num[0], modelcoords_num[1], modelcoords_num[2]});

    // Set velocity values (e.g., constant velocity)
    // Ckonal::real_t constant_velocity = 2.0; // units/second
    // velocity_model.fill(constant_velocity);
    // std::cout << "Velocity model configured with constant velocity: " << constant_velocity << std::endl;
    
    std::cout << "set velocity of Ckonal ..." << std::endl;
    // Set velocity values
    velocity_model.set_all_values(vmodel_values);

    // ---------- Ckonal mode 1 ----------
    if(mode_run == 1){
        // get geometry file
        int if_change_xy_for_geometry = ps.getInt("GEOMETRY_CHANGE_XY");
        std::string input_geometry_file = ps.getString("GEOMETRY_FILE");
        // read geometry
        Geometry geometry(input_geometry_file);
        std::cout << "calculate travel time field for each station ..." << std::endl;
        // station loop
        for(int i=0; i<geometry.node_number_; i++){
            // generate out travel time field file name
            std::ostringstream oss;
            oss << out_path_traveltime << path_seperator << "tt_";
            oss << std::setw(6) << std::setfill('0') << i+1;
            oss << ".npy";
            std::string ouf_tt_for_station_i = oss.str();
            // --- 2. Initialize Solver State and Add Source ---
            try {
                solver.init_solver_state(); // Initializes travel time field, FMM states, heap
                std::cout << "Solver state initialized." << std::endl;
                
                // set source location
                Ckonal::uint_t source_index_x, source_index_y, source_index_z;
                double source_x, source_y;
                if(if_change_xy_for_geometry > 0){
                    source_x = geometry.node_y_[i];
                    source_y = geometry.node_x_[i];
                }else{
                    source_x = geometry.node_x_[i];
                    source_y = geometry.node_y_[i];
                }
                source_index_x = std::round((source_x - modelcoords_min[0]) / modelcoords_step[0]);
                source_index_y = std::round((source_y - modelcoords_min[1]) / modelcoords_step[1]);
                Ckonal::Index3D source_location(source_index_x, source_index_y, 0); // Source at surface
                // Ckonal::Index3D source_location(0, 0, 0); // Source at a corner
                if (solver.add_source_point(source_location, 0.0)) {
                    std::cout << "Source point added at grid index ("
                            << source_location.i1 << "," << source_location.i2 << "," << source_location.i3 << ")." << std::endl;
                } else {
                    std::cerr << "Failed to add source point." << std::endl;
                    return -1;
                }
            } catch (const std::exception& e) {
                std::cerr << "Error during solver initialization or source addition: " << e.what() << std::endl;
                return -1;
            }

            // solve
            std::cout << "Solving Eikonal equation..." << std::endl;
            try {
                if (solver.solve()) {
                    std::cout << "Eikonal equation solved successfully." << std::endl;
                } else {
                    std::cerr << "Eikonal solver failed or had no sources." << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << "Error during Eikonal solve: " << e.what() << std::endl;
                return 1;
            }
            // solve end

            // --- 5. output travel time field ---
            const Ckonal::ScalarField3D& travel_time_result = solver.get_traveltime_field();
            const std::vector<Ckonal::real_t>& tt1d = travel_time_result.get_values_raw();
            npy::write_npy(ouf_tt_for_station_i, tt1d);
            std::cout << "write out travel time table for station " << i+1 << " as npy format: " << ouf_tt_for_station_i << std::endl;
        }
        // station loop end
    }

    // ---------- Ckonal mode 0 ----------
    if(mode_run == 0){
        // generate out travel time field file name
        std::ostringstream oss;
        oss << out_path_traveltime << path_seperator << "tt.npy";
        std::string ouf_tt = oss.str();
        // --- 2. Initialize Solver State and Add Source ---
        try {
            solver.init_solver_state(); // Initializes travel time field, FMM states, heap
            std::cout << "Solver state initialized." << std::endl;
            
            // set source location
            Ckonal::uint_t source_index_x, source_index_y, source_index_z;
            double source_x, source_y, source_z;
            source_x = ps.getFloat("CKONAL_SOURCE_X");
            source_y = ps.getFloat("CKONAL_SOURCE_Y");
            source_z = ps.getFloat("CKONAL_SOURCE_Z");
            std::cout << "calculate travel time field for source (" << source_x << "," << source_y << "," << source_z << ") ..." << std::endl;

            source_index_x = std::round((source_x - modelcoords_min[0]) / modelcoords_step[0]);
            source_index_y = std::round((source_y - modelcoords_min[1]) / modelcoords_step[1]);
            source_index_z = std::round((source_z - modelcoords_min[2]) / modelcoords_step[2]);

            Ckonal::Index3D source_location(source_index_x, source_index_y, source_index_z); // Source at surface
            // Ckonal::Index3D source_location(0, 0, 0); // Source at a corner
            if (solver.add_source_point(source_location, 0.0)) {
                std::cout << "Source point added at grid index ("
                        << source_location.i1 << "," << source_location.i2 << "," << source_location.i3 << ")." << std::endl;
            } else {
                std::cerr << "Failed to add source point." << std::endl;
                return -1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error during solver initialization or source addition: " << e.what() << std::endl;
            return -1;
        }

        // solve
        std::cout << "Solving Eikonal equation..." << std::endl;
        try {
            if (solver.solve()) {
                std::cout << "Eikonal equation solved successfully." << std::endl;
            } else {
                std::cerr << "Eikonal solver failed or had no sources." << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error during Eikonal solve: " << e.what() << std::endl;
            return 1;
        }
        // solve end

        // --- 5. output travel time field ---
        const Ckonal::ScalarField3D& travel_time_result = solver.get_traveltime_field();
        const std::vector<Ckonal::real_t>& tt1d = travel_time_result.get_values_raw();
        npy::write_npy(ouf_tt, tt1d);
        std::cout << "write out travel time table " << " as npy format: " << ouf_tt << std::endl;
    }
    
    // return
    return 0;
}

