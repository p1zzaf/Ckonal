#ifndef CKONAL_CONFIG_PARSER_HPP
#define CKONAL_CONFIG_PARSER_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept> // For std::runtime_error, std::invalid_argument
#include <algorithm> // For trim

namespace ckonal {
namespace app {

// Helper to trim whitespace from a string
inline std::string trim_string(const std::string& str) {
    const std::string whitespace = " \t\n\r\f\v";
    size_t start = str.find_first_not_of(whitespace);
    if (start == std::string::npos) return ""; // Empty or all whitespace
    size_t end = str.find_last_not_of(whitespace);
    return str.substr(start, end - start + 1);
}

struct Config {
    // Grid parameters
    double grid_x_min, grid_x_step;
    int grid_x_count;
    double grid_y_min, grid_y_step;
    int grid_y_count;
    double grid_z_min, grid_z_step;
    int grid_z_count;

    // Velocity model
    std::string velocity_model_path;

    // Source parameters
    double source_x, source_y, source_z;
    double source_travel_time;

    Config() : grid_x_min(0), grid_x_step(0), grid_x_count(0),
               grid_y_min(0), grid_y_step(0), grid_y_count(0),
               grid_z_min(0), grid_z_step(0), grid_z_count(0),
               source_x(0), source_y(0), source_z(0), source_travel_time(0.0) {}
};

class ConfigParser {
public:
    ConfigParser() = default;

    Config parse(const std::string& filename);

private:
    // Store parsed values temporarily before populating Config struct
    std::map<std::string, std::map<std::string, std::string>> sections;

    template<typename T>
    T get_value(const std::string& section, const std::string& key, const T& default_value) const;

    template<typename T>
    T get_required_value(const std::string& section, const std::string& key) const;
};


// Template specializations for get_value and get_required_value
template<>
inline std::string ConfigParser::get_value<std::string>(const std::string& section, const std::string& key, const std::string& default_value) const {
    auto sec_it = sections.find(section);
    if (sec_it != sections.end()) {
        auto key_it = sec_it->second.find(key);
        if (key_it != sec_it->second.end()) {
            return key_it->second;
        }
    }
    return default_value;
}

template<>
inline double ConfigParser::get_value<double>(const std::string& section, const std::string& key, const double& default_value) const {
    auto sec_it = sections.find(section);
    if (sec_it != sections.end()) {
        auto key_it = sec_it->second.find(key);
        if (key_it != sec_it->second.end()) {
            try {
                return std::stod(key_it->second);
            } catch (const std::invalid_argument& ia) {
                throw std::runtime_error("Invalid double value for " + section + "/" + key + ": " + key_it->second);
            } catch (const std::out_of_range& oor) {
                 throw std::runtime_error("Double value out of range for " + section + "/" + key + ": " + key_it->second);
            }
        }
    }
    return default_value;
}

template<>
inline int ConfigParser::get_value<int>(const std::string& section, const std::string& key, const int& default_value) const {
     auto sec_it = sections.find(section);
    if (sec_it != sections.end()) {
        auto key_it = sec_it->second.find(key);
        if (key_it != sec_it->second.end()) {
            try {
                return std::stoi(key_it->second);
            } catch (const std::invalid_argument& ia) {
                throw std::runtime_error("Invalid integer value for " + section + "/" + key + ": " + key_it->second);
            } catch (const std::out_of_range& oor) {
                 throw std::runtime_error("Integer value out of range for " + section + "/" + key + ": " + key_it->second);
            }
        }
    }
    return default_value;
}


template<typename T>
T ConfigParser::get_required_value(const std::string& section, const std::string& key) const {
    auto sec_it = sections.find(section);
    if (sec_it == sections.end()) {
        throw std::runtime_error("Required section missing in config file: [" + section + "]");
    }
    auto key_it = sec_it->second.find(key);
    if (key_it == sec_it->second.end()) {
        throw std::runtime_error("Required key '" + key + "' missing in section [" + section + "]");
    }
    // Use get_value with a dummy default to leverage its conversion and error handling
    return get_value<T>(section, key, T{}); // T{} is default constructor for T
}


} // namespace app
} // namespace ckonal
#endif // CKONAL_CONFIG_PARSER_HPP