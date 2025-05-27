#ifndef Ckonal_EIKONAL_SOLVER_HPP
#define Ckonal_EIKONAL_SOLVER_HPP

#include "../core/konal_types.hpp"
#include "../core/konal_constants.hpp"
#include "../fields/scalar_field.hpp"
#include "../structures/heap.hpp"
#include <vector>
#include <memory> // For std::unique_ptr for fields if solver owns them

namespace Ckonal {

// Enum for FMM node states
enum class FMMState : char {
    UNKNOWN,
    TRIAL,
    KNOWN
};


class EikonalSolver {
public:
    EikonalSolver(CoordSys coord_sys = CoordSys::CARTESIAN);
    virtual ~EikonalSolver() = default;

    EikonalSolver(const EikonalSolver&) = delete;
    EikonalSolver& operator=(const EikonalSolver&) = delete;
    EikonalSolver(EikonalSolver&&) = default;
    EikonalSolver& operator=(EikonalSolver&&) = default;

    ScalarField3D& get_velocity_field();
    const ScalarField3D& get_velocity_field() const;
    ScalarField3D& get_traveltime_field();
    const ScalarField3D& get_traveltime_field() const;

    virtual void init_solver_state();
    bool add_source_point(const Index3D& source_idx, real_t source_time = 0.0);
    
    virtual bool solve();

    // --- Debugging / State Inspection ---
    bool is_trial_heap_empty() const {
        if (m_trial_heap) {
            return m_trial_heap->empty();
        }
        return true; 
    }

    size_t get_trial_heap_size() const {
        if (m_trial_heap) {
            return m_trial_heap->size();
        }
        return 0;
    }

    FMMState get_fmm_state(const Index3D& idx) const;


protected:
    CoordSys m_coord_sys;
    ScalarField3D m_velocity_field;
    ScalarField3D m_traveltime_field;
    std::vector<FMMState> m_fmm_states;
    std::unique_ptr<Heap> m_trial_heap; 
    bool m_is_initialized;

    size_t get_flat_state_index(const Index3D& idx) const;
    void update_neighbor(const Index3D& active_idx, const Index3D& neighbor_idx);
    real_t calculate_updated_time(const Index3D& target_node_idx);
};


class PointSourceSolver : public EikonalSolver {
public:
    PointSourceSolver(CoordSys far_field_coord_sys = CoordSys::CARTESIAN);

    void set_source_location(const Point3D& loc);
    Point3D get_source_location() const;

    void set_near_field_radial_nodes(uint_t n_rho);
    void set_near_field_polar_nodes(uint_t n_theta);
    void set_near_field_azimuthal_nodes(uint_t n_phi);
    void set_near_field_radial_interval(real_t drho);

    void init_solver_state() override; 
    bool solve() override;            

private:
    std::unique_ptr<EikonalSolver> m_near_field_solver;
    Point3D m_source_loc_far_field; 

    uint_t m_n_rho;   
    uint_t m_n_theta; 
    uint_t m_n_phi;   
    real_t m_drho;    

    static constexpr uint_t DEFAULT_N_RHO = 64;
    static constexpr uint_t DEFAULT_N_THETA = 33; 
    static constexpr uint_t DEFAULT_N_PHI = 64;   
    
    bool m_source_loc_set;

    void initialize_near_field_grid_parameters();
    void interpolate_far_field_velocity_to_near_field();
    void initialize_near_field_source_band();
    void interpolate_near_field_traveltime_to_far_field();
    void initialize_far_field_source_band_from_near_field();
};


} // namespace Ckonal

#endif // Ckonal_EIKONAL_SOLVER_HPP