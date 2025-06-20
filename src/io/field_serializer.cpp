#include "field_serializer.hpp"
#include <stdexcept> // For runtime_error

namespace Ckonal {
namespace io {

// Magic number and version for basic file format validation
// In io/field_serializer.cpp
// Example: Choose hex values that are somewhat memorable or distinct
const uint32_t MAGIC_NUMBER_SCALAR = 0x504B3543; // PK5C (ASCII approx)
const uint32_t MAGIC_NUMBER_VECTOR = 0x504B5645; // PKVE (ASCII approx)
// Or simpler ones:
// const uint32_t MAGIC_NUMBER_SCALAR = 0xA1B2C3D0;
// const uint32_t MAGIC_NUMBER_VECTOR = 0xA1B2C3D1;
const uint16_t FILE_FORMAT_VERSION = 1;

void BinaryFieldSerializer::save_scalar_field(const ScalarField3D& field, const std::string& path) const {
    std::ofstream outfile(path, std::ios::binary | std::ios::trunc);
    if (!outfile.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    // 1. Magic number & Version
    write_binary(outfile, MAGIC_NUMBER_SCALAR);
    write_binary(outfile, FILE_FORMAT_VERSION);

    // 2. Coordinate System (as int/enum value)
    CoordSys cs_val = field.get_coord_sys();
    write_binary(outfile, static_cast<int>(cs_val)); // Store enum as underlying type

    // 3. Grid Parameters
    write_binary_array(outfile, field.get_min_coords());
    write_binary_array(outfile, field.get_node_intervals());
    write_binary_array(outfile, field.get_npts());
    // Derived properties (max_coords, is_null, is_periodic) will be recomputed on load.

    // 4. Scalar Values
    write_binary_vector(outfile, field.get_values_raw());

    if (!outfile.good()) {
         throw std::runtime_error("Error occurred during writing to file: " + path);
    }
    outfile.close();
}

std::unique_ptr<ScalarField3D> BinaryFieldSerializer::load_scalar_field(const std::string& path) const {
    std::ifstream infile(path, std::ios::binary);
    if (!infile.is_open()) {
        throw std::runtime_error("Failed to open file for reading: " + path);
    }

    uint32_t magic;
    uint16_t version;
    read_binary(infile, magic);
    read_binary(infile, version);

    if (magic != MAGIC_NUMBER_SCALAR) {
        throw std::runtime_error("Invalid file format (magic number mismatch) for scalar field: " + path);
    }
    if (version != FILE_FORMAT_VERSION) {
        throw std::runtime_error("Unsupported file format version for scalar field: " + path);
    }

    int cs_int;
    read_binary(infile, cs_int);
    CoordSys cs = static_cast<CoordSys>(cs_int);

    auto field = std::make_unique<ScalarField3D>(cs);

    std::array<real_t, 3> min_coords_val;
    std::array<real_t, 3> node_intervals_val;
    std::array<uint_t, 3> npts_val;

    read_binary_array(infile, min_coords_val);
    read_binary_array(infile, node_intervals_val);
    read_binary_array(infile, npts_val);

    field->set_min_coords(min_coords_val);
    field->set_node_intervals(node_intervals_val);
    field->set_npts(npts_val); // This will also size its internal values vector

    std::vector<real_t> values_data;
    read_binary_vector(infile, values_data);
    if (values_data.size() != field->get_total_nodes()) {
        throw std::runtime_error("Data size mismatch when loading scalar field values from: " + path);
    }
    field->set_all_values(values_data);


    if (!infile.good() && !infile.eof()) { // Check for errors other than EOF
         throw std::runtime_error("Error occurred during reading from file: " + path);
    }
    infile.close();
    return field;
}


void BinaryFieldSerializer::save_vector_field(const VectorField3D& field, const std::string& path) const {
    std::ofstream outfile(path, std::ios::binary | std::ios::trunc);
    if (!outfile.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    write_binary(outfile, MAGIC_NUMBER_VECTOR);
    write_binary(outfile, FILE_FORMAT_VERSION);

    CoordSys cs_val = field.get_coord_sys();
    write_binary(outfile, static_cast<int>(cs_val));

    write_binary_array(outfile, field.get_min_coords());
    write_binary_array(outfile, field.get_node_intervals());
    write_binary_array(outfile, field.get_npts());

    // Vector Values (N0*N1*N2*3)
    write_binary_vector(outfile, field.get_values_raw());

    if (!outfile.good()) {
         throw std::runtime_error("Error occurred during writing to file: " + path);
    }
    outfile.close();
}

std::unique_ptr<VectorField3D> BinaryFieldSerializer::load_vector_field(const std::string& path) const {
    std::ifstream infile(path, std::ios::binary);
    if (!infile.is_open()) {
        throw std::runtime_error("Failed to open file for reading: " + path);
    }

    uint32_t magic;
    uint16_t version;
    read_binary(infile, magic);
    read_binary(infile, version);

    if (magic != MAGIC_NUMBER_VECTOR) {
        throw std::runtime_error("Invalid file format (magic number mismatch) for vector field: " + path);
    }
    if (version != FILE_FORMAT_VERSION) {
        throw std::runtime_error("Unsupported file format version for vector field: " + path);
    }

    int cs_int;
    read_binary(infile, cs_int);
    CoordSys cs = static_cast<CoordSys>(cs_int);

    auto field = std::make_unique<VectorField3D>(cs);

    std::array<real_t, 3> min_coords_val;
    std::array<real_t, 3> node_intervals_val;
    std::array<uint_t, 3> npts_val;

    read_binary_array(infile, min_coords_val);
    read_binary_array(infile, node_intervals_val);
    read_binary_array(infile, npts_val);

    field->set_min_coords(min_coords_val);
    field->set_node_intervals(node_intervals_val);
    field->set_npts(npts_val);

    std::vector<real_t> values_data;
    read_binary_vector(infile, values_data);
     if (values_data.size() != field->get_total_nodes() * 3) {
        throw std::runtime_error("Data size mismatch when loading vector field values from: " + path);
    }
    field->set_all_values(values_data);

    if (!infile.good() && !infile.eof()) {
         throw std::runtime_error("Error occurred during reading from file: " + path);
    }
    infile.close();
    return field;
}


} // namespace io
} // namespace Ckonal