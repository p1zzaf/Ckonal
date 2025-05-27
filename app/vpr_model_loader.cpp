#include "vpr_model_loader.hpp"
#include "../src/core/konal_constants.hpp" // For Ckonal::uint_t if needed directly (though ScalarField3D uses it)
#include <iostream> 
#include <iomanip>  

namespace ckonal {
namespace app {

bool VprModelLoader::read_next_significant_line(std::ifstream& infile, std::string& line, int& current_line_num) {
    while (std::getline(infile, line)) {
        current_line_num++;
        line = trim_string(line);
        if (!line.empty()) {
            return true; 
        }
    }
    return false; 
}


bool VprModelLoader::load(const std::string& vpr_file_path,
                          Ckonal::ScalarField3D& velocity_field, // <--- MODIFIED: Use Ckonal::
                          double& out_value_min,
                          double& out_value_max) {
    std::ifstream infile(vpr_file_path);
    if (!infile.is_open()) {
        std::cerr << "Error: Failed to open VPR model file: " << vpr_file_path << std::endl;
        return false;
    }

    const auto& npts = velocity_field.get_npts();
    Ckonal::uint_t nx = npts[0]; // <--- MODIFIED: Use Ckonal::
    Ckonal::uint_t ny = npts[1]; // <--- MODIFIED: Use Ckonal::
    Ckonal::uint_t nz = npts[2]; // <--- MODIFIED: Use Ckonal::

    if (nx == 0 || ny == 0 || nz == 0) {
        std::cerr << "Error: Velocity field dimensions (nx, ny, nz) are not set or are zero before loading VPR model." << std::endl;
        return false;
    }

    std::string line;
    int line_num = 0;

    if (!read_next_significant_line(infile, line, line_num)) {
        std::cerr << "Error: VPR file " << vpr_file_path << " is empty or contains no data after header." << std::endl;
        return false;
    }
    std::istringstream header_ss(line);
    if (!(header_ss >> out_value_min >> out_value_max)) {
        std::cerr << "Error: Failed to parse header (value_min value_max) from VPR file " << vpr_file_path << " at line " << line_num << std::endl;
        return false;
    }

    for (Ckonal::uint_t i = 0; i < nx; ++i) { // <--- MODIFIED
        for (Ckonal::uint_t j = 0; j < ny; ++j) { // <--- MODIFIED
            if (!read_next_significant_line(infile, line, line_num)) {
                std::cerr << "Error: Premature end of file in VPR model " << vpr_file_path
                          << ". Expected Y-slice " << j << " for X-slice " << i << " at line approx " << line_num << "." << std::endl;
                return false;
            }

            std::istringstream data_ss(line);
            for (Ckonal::uint_t k = 0; k < nz; ++k) { // <--- MODIFIED
                double val;
                if (!(data_ss >> val)) {
                    std::cerr << "Error: Failed to read Z-value " << k << " for (X=" << i << ", Y=" << j
                              << ") from VPR file " << vpr_file_path << " at line " << line_num << "." << std::endl;
                    std::cerr << "Line content: '" << line << "'" << std::endl;
                    return false;
                }
                try {
                    velocity_field.set_value({i, j, k}, val);
                } catch (const std::out_of_range& e) {
                     std::cerr << "Error: Grid index (" << i << "," << j << "," << k 
                               << ") out of bounds for velocity field when setting VPR value. "
                               << "VPR dimensions might not match config. " << e.what() << std::endl;
                    return false;
                }
            }
            double extra_val_check;
            if (data_ss >> extra_val_check) {
                std::cerr << "Warning: Extra data found on line " << line_num << " for (X=" << i << ", Y=" << j
                          << ") in VPR file " << vpr_file_path << ". Expected " << nz << " Z-values." << std::endl;
            }
        }
        if (i < nx - 1) {
            std::string temp_line_for_sep;
            if (std::getline(infile, temp_line_for_sep)) { 
                line_num++;
                temp_line_for_sep = trim_string(temp_line_for_sep);
                if (!temp_line_for_sep.empty()) {
                    std::cerr << "Warning: Expected empty separator line after X-slice " << i
                              << " in VPR file " << vpr_file_path << " at line " << line_num
                              << ", but found: '" << temp_line_for_sep << "'. Will attempt to continue." << std::endl;
                }
            } else if (!infile.eof()){ 
                 std::cerr << "Warning: Error reading potential separator line after X-slice " << i << std::endl;
            }
        }
    }

    if (read_next_significant_line(infile, line, line_num)) {
        std::cerr << "Warning: Extra data found at the end of VPR file " << vpr_file_path << " starting at line " << line_num << ": '" << line << "'" << std::endl;
    }

    infile.close();
    return true;
}

} // namespace app
} // namespace ckonal