#include "config_parser.hpp"
#include <iostream> // For potential debug output

namespace ckonal {
namespace app {

Config ConfigParser::parse(const std::string& filename) {
    sections.clear();
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Failed to open config file: " + filename);
    }

    std::string line;
    std::string current_section;
    int line_number = 0;

    while (std::getline(infile, line)) {
        line_number++;
        line = trim_string(line);

        if (line.empty() || line[0] == '#' || line[0] == ';') { // Skip empty lines and comments
            continue;
        }

        if (line[0] == '[' && line.back() == ']') { // Section header
            current_section = trim_string(line.substr(1, line.length() - 2));
            if (current_section.empty()) {
                throw std::runtime_error("Syntax error in config file " + filename + " at line " + std::to_string(line_number) + ": Empty section name.");
            }
            sections[current_section] = {}; // Create section if not exists
        } else if (!current_section.empty()) { // Key-value pair
            size_t delimiter_pos = line.find('=');
            if (delimiter_pos != std::string::npos) {
                std::string key = trim_string(line.substr(0, delimiter_pos));
                std::string value = trim_string(line.substr(delimiter_pos + 1));
                if (key.empty()) {
                    throw std::runtime_error("Syntax error in config file " + filename + " at line " + std::to_string(line_number) + ": Empty key name in section [" + current_section + "].");
                }
                sections[current_section][key] = value;
            } else {
                 throw std::runtime_error("Syntax error in config file " + filename + " at line " + std::to_string(line_number) + ": Missing '=' in key-value pair in section [" + current_section + "].");
            }
        } else {
             throw std::runtime_error("Syntax error in config file " + filename + " at line " + std::to_string(line_number) + ": Key-value pair outside of any section.");
        }
    }
    infile.close();

    // Populate Config struct
    Config config;
    config.grid_x_min   = get_required_value<double>("Grid", "x_min");
    config.grid_x_step  = get_required_value<double>("Grid", "x_step");
    config.grid_x_count = get_required_value<int>("Grid", "x_count");
    config.grid_y_min   = get_required_value<double>("Grid", "y_min");
    config.grid_y_step  = get_required_value<double>("Grid", "y_step");
    config.grid_y_count = get_required_value<int>("Grid", "y_count");
    config.grid_z_min   = get_required_value<double>("Grid", "z_min");
    config.grid_z_step  = get_required_value<double>("Grid", "z_step");
    config.grid_z_count = get_required_value<int>("Grid", "z_count");

    config.velocity_model_path = get_required_value<std::string>("VelocityModel", "file_path");

    config.source_x = get_required_value<double>("Source", "x_coord");
    config.source_y = get_required_value<double>("Source", "y_coord");
    config.source_z = get_required_value<double>("Source", "z_coord");
    config.source_travel_time = get_value<double>("Source", "travel_time", 0.0); // Default to 0.0 if missing

    return config;
}

} // namespace app
} // namespace ckonal