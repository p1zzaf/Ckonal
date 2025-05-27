#ifndef CKONAL_VPR_MODEL_LOADER_HPP
#define CKONAL_VPR_MODEL_LOADER_HPP

// Use the correct namespace for ScalarField3D from the core library
#include "../src/fields/scalar_field.hpp" // This defines Ckonal::ScalarField3D
#include "config_parser.hpp" 
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ckonal { // Keep app specific code in ckonal namespace if desired
namespace app {

class VprModelLoader {
public:
    VprModelLoader() = default;

    bool load(const std::string& vpr_file_path,
              Ckonal::ScalarField3D& velocity_field, // <--- MODIFIED: Use Ckonal::
              double& out_value_min,    
              double& out_value_max);   

private:
    bool read_next_significant_line(std::ifstream& infile, std::string& line, int& current_line_num);
};

} // namespace app
} // namespace ckonal

#endif // CKONAL_VPR_MODEL_LOADER_HPP