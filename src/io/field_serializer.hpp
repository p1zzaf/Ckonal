#ifndef Ckonal_FIELD_SERIALIZER_HPP
#define Ckonal_FIELD_SERIALIZER_HPP

#include "../fields/scalar_field.hpp"
#include "../fields/vector_field.hpp"
#include <string>
#include <fstream> // For file streams
#include <memory>

namespace Ckonal {
namespace io {

// --- Abstract Serializer/Deserializer Interface (Optional but good practice) ---
// class IFieldSerializer {
// public:
//     virtual ~IFieldSerializer() = default;
//     virtual void save_scalar_field(const ScalarField3D& field, const std::string& path) const = 0;
//     virtual std::unique_ptr<ScalarField3D> load_scalar_field(const std::string& path) const = 0;
//     // Add similar for VectorField3D
// };

// --- Simple Binary Field Serializer ---
// For now, a concrete class is fine.
class BinaryFieldSerializer /* : public IFieldSerializer */ {
public:
    void save_scalar_field(const ScalarField3D& field, const std::string& path) const;
    std::unique_ptr<ScalarField3D> load_scalar_field(const std::string& path) const;

    void save_vector_field(const VectorField3D& field, const std::string& path) const;
    std::unique_ptr<VectorField3D> load_vector_field(const std::string& path) const;

private:
    // Helper to write basic types
    template<typename T>
    void write_binary(std::ofstream& out, const T& value) const {
        out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    template<typename T>
    void read_binary(std::ifstream& in, T& value) const {
        in.read(reinterpret_cast<char*>(&value), sizeof(T));
    }

    // Helper to write std::array
    template<typename T, size_t N>
    void write_binary_array(std::ofstream& out, const std::array<T, N>& arr) const {
        for (const auto& val : arr) {
            write_binary(out, val);
        }
    }
    template<typename T, size_t N>
    void read_binary_array(std::ifstream& in, std::array<T, N>& arr) const {
        for (auto& val : arr) {
            read_binary(in, val);
        }
    }

    // Helper to write std::vector<real_t> (for values)
    void write_binary_vector(std::ofstream& out, const std::vector<real_t>& vec) const {
        size_t size = vec.size();
        write_binary(out, size);
        out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(real_t));
    }
    void read_binary_vector(std::ifstream& in, std::vector<real_t>& vec) const {
        size_t size;
        read_binary(in, size);
        vec.resize(size);
        in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(real_t));
    }
};

// --- Placeholder for TravelTimeInventory replacement ---
// If EQLocator needs to save/load multiple fields, this would be the place.
// For now, EQLocator might just use multiple calls to BinaryFieldSerializer.
// class TravelTimeDataStore {
// public:
//     // Methods to add/read fields by a key (e.g., station/phase)
//     // to/from a directory structure or a single container file.
// };


} // namespace io
} // namespace Ckonal

#endif // Ckonal_FIELD_SERIALIZER_HPP