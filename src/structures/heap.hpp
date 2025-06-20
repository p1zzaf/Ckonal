#ifndef Ckonal_HEAP_HPP
#define Ckonal_HEAP_HPP

#include "core/konal_types.hpp"
#include "core/konal_constants.hpp"
#include <vector>
#include <functional> // For std::less, std::greater
#include <stdexcept>  // For std::out_of_range
#include <iostream>   // For potential debugging

namespace Ckonal {

// Forward declaration for Field classes if needed for values_ptr_
// class ScalarField3D; // If Heap directly accesses ScalarField3D internals

class Heap {
public:
    // The 'values' array (e.g., travel times) is external to the heap.
    // The heap stores indices (Index3D) and uses these external values for comparison.
    // We need a way to access these values. A raw pointer or a reference to the
    // data container can work. For simplicity with flattened arrays, a pointer
    // and dimensions might be easier.
    Heap(const real_t* values_ptr, const std::array<size_t, 3>& dims);

    // Disable copy and move semantics for now, or implement properly if needed.
    Heap(const Heap&) = delete;
    Heap& operator=(const Heap&) = delete;
    Heap(Heap&&) = default; // Default move constructor
    Heap& operator=(Heap&&) = default; // Default move assignment

    void push(const Index3D& idx);
    Index3D pop();
    bool empty() const;
    size_t size() const;

    // To be called if the value associated with an index already in the heap changes.
    // 'heap_pos' is the current position of the item in the m_keys vector.
    void update_key_at(size_t heap_pos);

    // Access to the heap_index array (maps grid Index3D to its position in m_keys)
    // This is mutable because pop/push modifies it.
    // The external user of Heap will typically create and pass this array.
    // For better encapsulation, Heap could own this, but original design has it separate.
    // Let's try with Heap owning it for now.
    // If an Index3D is not in heap, its heap_index_data value is HEAP_INVALID_INDEX.
    // This method is to get the heap position of a grid node.
    int_t get_heap_position(const Index3D& grid_idx) const;


private:
    std::vector<Index3D> m_keys; // Stores the actual heap (elements are Index3D)
    
    // Data for sorting:
    // Pointer to the 3D grid data (e.g., travel times) that the heap sorts by.
    // The heap itself doesn't own this data.
    const real_t* m_values_ptr; 
    std::array<size_t, 3> m_value_dims; // Dimensions of the m_values_ptr grid (N0, N1, N2)

    // Mapping from grid Index3D to its position in m_keys.
    // Flattened array: size = N0 * N1 * N2.
    // Stores the index within m_keys, or HEAP_INVALID_INDEX if not in heap.
    std::vector<int_t> m_heap_index_data; 

    // Helper to get flat index for m_values_ptr and m_heap_index_data
    inline size_t get_flat_value_index(const Index3D& idx) const {
        // TODO: Add bounds checking if necessary, or ensure idx is always valid
        return idx.i1 * (m_value_dims[1] * m_value_dims[2]) + idx.i2 * m_value_dims[2] + idx.i3;
    }

    // Comparison function using m_values_ptr
    bool compare_indices(const Index3D& idx1, const Index3D& idx2) const {
        return m_values_ptr[get_flat_value_index(idx1)] < m_values_ptr[get_flat_value_index(idx2)];
    }

    void sift_up(size_t k);   // aka bubble_up, swim. For push.
    void sift_down(size_t k); // aka bubble_down, sink. For pop & update.

    // For sift_up and sift_down, they might need to know the start/end bounds
    // as in the original sift_down(j_start, j) and sift_up(j_start)
    // Let's stick to standard heap sift operations first.
    // The original sift_down(j_start, j) sifts element j towards j_start (root direction).
    // The original sift_up(j_start) sifts element j_start away from root.
    // These names are a bit confusing compared to typical heap algos.
    // Standard sift_up (from child to root):
    void heapify_up(size_t child_idx);
    // Standard sift_down (from parent to child):
    void heapify_down(size_t parent_idx);

    void swap_keys(size_t i, size_t j);
};

} // namespace Ckonal

#endif // Ckonal_HEAP_HPP