#include "heap.hpp"
#include <algorithm> // For std::swap

namespace Ckonal {

Heap::Heap(const real_t* values_ptr, const std::array<size_t, 3>& dims)
    : m_values_ptr(values_ptr), m_value_dims(dims) {
    if (!values_ptr) {
        throw std::invalid_argument("Heap values_ptr cannot be null.");
    }
    if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0) {
        throw std::invalid_argument("Heap dimensions cannot be zero.");
    }
    m_heap_index_data.resize(dims[0] * dims[1] * dims[2], HEAP_INVALID_INDEX);
    m_keys.reserve(1024); // Pre-allocate some space
}

bool Heap::empty() const {
    return m_keys.empty();
}

size_t Heap::size() const {
    return m_keys.size();
}

void Heap::swap_keys(size_t i, size_t j) {
    std::swap(m_keys[i], m_keys[j]);
    // Update heap_index_data for the swapped elements
    m_heap_index_data[get_flat_value_index(m_keys[i])] = static_cast<int_t>(i);
    m_heap_index_data[get_flat_value_index(m_keys[j])] = static_cast<int_t>(j);
}

// Standard sift_up (bubble up from child_idx towards root)
void Heap::heapify_up(size_t child_idx) {
    if (child_idx == 0) return; // Root
    size_t parent_idx = (child_idx - 1) / 2;
    // Min-heap: if child is smaller than parent, swap
    while (child_idx > 0 && compare_indices(m_keys[child_idx], m_keys[parent_idx])) {
        swap_keys(child_idx, parent_idx);
        child_idx = parent_idx;
        if (child_idx == 0) break;
        parent_idx = (child_idx - 1) / 2;
    }
}

// Standard sift_down (bubble down from parent_idx towards leaves)
void Heap::heapify_down(size_t parent_idx) {
    size_t N = m_keys.size();
    if (N == 0) return;

    size_t current_idx = parent_idx;
    while (true) {
        size_t left_child_idx = 2 * current_idx + 1;
        size_t right_child_idx = 2 * current_idx + 2;
        size_t smallest_idx = current_idx;

        if (left_child_idx < N && compare_indices(m_keys[left_child_idx], m_keys[smallest_idx])) {
            smallest_idx = left_child_idx;
        }
        if (right_child_idx < N && compare_indices(m_keys[right_child_idx], m_keys[smallest_idx])) {
            smallest_idx = right_child_idx;
        }

        if (smallest_idx != current_idx) {
            swap_keys(current_idx, smallest_idx);
            current_idx = smallest_idx; // Move down to the swapped child
        } else {
            break; // Heap property satisfied at this subtree
        }
    }
}


void Heap::push(const Index3D& grid_idx) {
    // Check if already in heap? Original code seems to allow re-pushing if it's
    // handled by 'unknown' flag logic in solver, or if update logic is used.
    // For FMM, if a node is re-evaluated and its time improves, it might be
    // "re-pushed" or its key updated.
    // If it's already there and value changes, `update_key_at` should be used.
    // If we are strictly adding new items:
    if (m_heap_index_data[get_flat_value_index(grid_idx)] != HEAP_INVALID_INDEX) {
        // This could be an error, or it means we need to update.
        // For now, let's assume push is for new items or items whose values might have changed
        // and are now being (re)considered.
        // The original `solver.pyx` pushes if `unknown[nbr]` is true,
        // otherwise it calls `trial.sift_down(0, heap_index[nbr])` which implies an update.
        // This heap's `update_key_at` matches the `sift_down` usage.
        // So, push should be for *new* items to the heap.
        // Let's make `push` simpler: always add, then sift up.
        // The caller (Solver) must manage if it's an update or a new push.
    }

    m_keys.push_back(grid_idx);
    size_t new_item_heap_idx = m_keys.size() - 1;
    m_heap_index_data[get_flat_value_index(grid_idx)] = static_cast<int_t>(new_item_heap_idx);
    heapify_up(new_item_heap_idx);
}

Index3D Heap::pop() {
    if (empty()) {
        throw std::out_of_range("Heap is empty, cannot pop.");
    }
    Index3D top_item = m_keys[0];
    
    // Mark popped item as no longer in heap
    m_heap_index_data[get_flat_value_index(top_item)] = HEAP_INVALID_INDEX;

    if (m_keys.size() > 1) {
        m_keys[0] = m_keys.back(); // Move last item to root
        m_heap_index_data[get_flat_value_index(m_keys[0])] = 0; // Update its index
    }
    m_keys.pop_back();

    if (!empty()) {
        heapify_down(0); // Restore heap property from root
    }
    return top_item;
}

// This is the crucial "update" function.
// If the value for m_keys[heap_pos] in m_values_ptr has changed, call this.
void Heap::update_key_at(size_t heap_pos) {
    if (heap_pos >= m_keys.size()) {
        throw std::out_of_range("Invalid heap_pos for update_key_at.");
    }
    // The value has changed. It could be smaller (needs sift_up) or larger (needs sift_down).
    // A common trick is to try sift_up, and if it doesn't move, then try sift_down.
    // Or, compare with parent to decide.
    if (heap_pos > 0) {
        size_t parent_idx = (heap_pos - 1) / 2;
        if (compare_indices(m_keys[heap_pos], m_keys[parent_idx])) {
            heapify_up(heap_pos);
            return; // If it moved up, it's done.
        }
    }
    // If it didn't move up (or it's the root), try sifting down.
    heapify_down(heap_pos);
}


int_t Heap::get_heap_position(const Index3D& grid_idx) const {
    // Add bounds check for grid_idx if necessary
    return m_heap_index_data[get_flat_value_index(grid_idx)];
}


} // namespace Ckonal