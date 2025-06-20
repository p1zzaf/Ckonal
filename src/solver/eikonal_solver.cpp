#include "eikonal_solver.hpp"
#include "utils/math_utils.hpp" 
#include "numerics/transformations.hpp" 
#include <algorithm> 
#include <cmath>     
#include <stdexcept> 
#include <limits>    
#include <iostream> // For debugging cout/cerr

namespace Ckonal {

// --- EikonalSolver Implementation ---
EikonalSolver::EikonalSolver(CoordSys coord_sys)
    : m_coord_sys(coord_sys),
      m_velocity_field(coord_sys),       
      m_traveltime_field(coord_sys),     
      m_is_initialized(false) {}

ScalarField3D& EikonalSolver::get_velocity_field() { return m_velocity_field; }
const ScalarField3D& EikonalSolver::get_velocity_field() const { return m_velocity_field; }
ScalarField3D& EikonalSolver::get_traveltime_field() { return m_traveltime_field; }
const ScalarField3D& EikonalSolver::get_traveltime_field() const { return m_traveltime_field; }

size_t EikonalSolver::get_flat_state_index(const Index3D& idx) const {
    return m_traveltime_field.get_flat_index(idx);
}

FMMState EikonalSolver::get_fmm_state(const Index3D& idx) const {
    if (!m_is_initialized) return FMMState::UNKNOWN; 
    return m_fmm_states[get_flat_state_index(idx)];
}

void EikonalSolver::init_solver_state() {
    if (m_velocity_field.get_total_nodes() == 0 ||
        m_velocity_field.get_npts()[0] == 0) { 
        throw std::runtime_error("Velocity field grid must be configured (npts, min_coords, intervals) before initializing solver.");
    }

    m_traveltime_field.set_coord_sys(m_velocity_field.get_coord_sys());
    m_traveltime_field.set_min_coords(m_velocity_field.get_min_coords());
    m_traveltime_field.set_node_intervals(m_velocity_field.get_node_intervals());
    m_traveltime_field.set_npts(m_velocity_field.get_npts()); 

    m_traveltime_field.fill(INF);
    m_fmm_states.assign(m_traveltime_field.get_total_nodes(), FMMState::UNKNOWN);

    const auto& npts_uint = m_traveltime_field.get_npts();
    std::array<size_t, 3> npts_size_t = {
        static_cast<size_t>(npts_uint[0]),
        static_cast<size_t>(npts_uint[1]),
        static_cast<size_t>(npts_uint[2])
    };
    m_trial_heap = std::make_unique<Heap>(m_traveltime_field.get_values_raw().data(), npts_size_t);

    m_is_initialized = true;
}

bool EikonalSolver::add_source_point(const Index3D& source_idx, real_t source_time) {
    if (!m_is_initialized) {
        throw std::runtime_error("Solver not initialized. Call init_solver_state() before adding sources.");
    }

    const auto& npts = m_traveltime_field.get_npts();
    if (source_idx.i1 >= npts[0] || source_idx.i2 >= npts[1] || source_idx.i3 >= npts[2]) {
        return false; 
    }

    if (m_velocity_field.get_value(source_idx) <= 1e-9) { 
        return false;
    }

    size_t flat_idx = get_flat_state_index(source_idx);

    m_traveltime_field.set_value(source_idx, source_time);
    m_fmm_states[flat_idx] = FMMState::TRIAL;
    m_trial_heap->push(source_idx); 

    return true;
}


bool EikonalSolver::solve() {
    if (!m_is_initialized) {
        throw std::runtime_error("Solver not initialized. Call init_solver_state() and add sources before solving.");
    }
    if (m_trial_heap->empty()) {
        return false; 
    }

    const auto& npts_field = m_traveltime_field.get_npts();
    const auto& is_periodic_field = m_traveltime_field.get_axis_is_periodic();

    while (!m_trial_heap->empty()) {
        Index3D active_idx = m_trial_heap->pop();
        size_t flat_active_idx = get_flat_state_index(active_idx);

        if (m_fmm_states[flat_active_idx] == FMMState::KNOWN) {
            continue;
        }
        m_fmm_states[flat_active_idx] = FMMState::KNOWN;

        for (int axis = 0; axis < 3; ++axis) {      
            for (int dir = -1; dir <= 1; dir += 2) { 
                Index3D neighbor_idx = active_idx;
                bool is_valid_neighbor = true;

                long long current_comp_val_long = 0; 
                if (axis == 0) current_comp_val_long = static_cast<long long>(active_idx.i1);
                else if (axis == 1) current_comp_val_long = static_cast<long long>(active_idx.i2);
                else current_comp_val_long = static_cast<long long>(active_idx.i3);
                
                long long neighbor_comp_val_long = current_comp_val_long + dir;

                if (is_periodic_field[axis]) {
                    neighbor_comp_val_long = (neighbor_comp_val_long % static_cast<long long>(npts_field[axis]) +
                                             static_cast<long long>(npts_field[axis])) % static_cast<long long>(npts_field[axis]);
                } else {
                    if (neighbor_comp_val_long < 0 || neighbor_comp_val_long >= static_cast<long long>(npts_field[axis])) {
                        is_valid_neighbor = false;
                    }
                }
                
                if (!is_valid_neighbor) continue;

                if (axis == 0) neighbor_idx.i1 = static_cast<size_t>(neighbor_comp_val_long);
                else if (axis == 1) neighbor_idx.i2 = static_cast<size_t>(neighbor_comp_val_long);
                else neighbor_idx.i3 = static_cast<size_t>(neighbor_comp_val_long);

                update_neighbor(active_idx, neighbor_idx);
            }
        }
    }
    return true;
}

void EikonalSolver::update_neighbor(const Index3D& /*active_idx*/, const Index3D& neighbor_idx) {
    size_t flat_neighbor_idx = get_flat_state_index(neighbor_idx);

    if (m_fmm_states[flat_neighbor_idx] == FMMState::KNOWN) {
        return; 
    }

    real_t neighbor_velocity = m_velocity_field.get_value(neighbor_idx);
    if (neighbor_velocity <= 1e-9) { 
        return; 
    }

    real_t old_time = m_traveltime_field.get_value(neighbor_idx); 
    real_t new_time = calculate_updated_time(neighbor_idx);

    if (new_time < old_time) {
        m_traveltime_field.set_value(neighbor_idx, new_time);
        if (m_fmm_states[flat_neighbor_idx] == FMMState::UNKNOWN) {
            m_fmm_states[flat_neighbor_idx] = FMMState::TRIAL;
            m_trial_heap->push(neighbor_idx);
        } else { 
            int_t heap_pos = m_trial_heap->get_heap_position(neighbor_idx);
            if (heap_pos != HEAP_INVALID_INDEX) { 
                m_trial_heap->update_key_at(static_cast<size_t>(heap_pos));
            } else {
                 m_trial_heap->push(neighbor_idx); 
            }
        }
    }
}

real_t EikonalSolver::calculate_updated_time(const Index3D& target_node_idx) {
    std::array<real_t, 3> aa = {0.0, 0.0, 0.0};
    std::array<real_t, 3> bb = {0.0, 0.0, 0.0};
    std::array<real_t, 3> cc = {0.0, 0.0, 0.0};

    const auto& npts = m_traveltime_field.get_npts();
    const auto& is_periodic = m_traveltime_field.get_axis_is_periodic();

    for (int dim_axis = 0; dim_axis < 3; ++dim_axis) {
        real_t h_dim = m_traveltime_field.get_norm_component(target_node_idx, dim_axis);

        if (std::abs(h_dim) < 1e-9) { 
            continue;
        }

        real_t T_candidates[2] = {INF, INF}; 
        bool is_second_order_candidate[2] = {false, false}; // To track if T_candidate was from 2nd order
        real_t T1_for_candidate[2] = {INF, INF};
        real_t T2_for_candidate[2] = {INF, INF};


        for (int dir = -1; dir <= 1; dir += 2) { 
            Index3D neighbor1_idx = target_node_idx;
            Index3D neighbor2_idx = target_node_idx; 
            bool n1_valid = true;
            bool n2_valid = false; 

            long long target_comp_val = 0;
            if (dim_axis == 0) target_comp_val = static_cast<long long>(target_node_idx.i1);
            else if (dim_axis == 1) target_comp_val = static_cast<long long>(target_node_idx.i2);
            else target_comp_val = static_cast<long long>(target_node_idx.i3);

            long long n1_comp_val = target_comp_val + dir;
            if (is_periodic[dim_axis]) {
                n1_comp_val = (n1_comp_val % static_cast<long long>(npts[dim_axis]) + static_cast<long long>(npts[dim_axis])) % static_cast<long long>(npts[dim_axis]);
            } else {
                if (n1_comp_val < 0 || n1_comp_val >= static_cast<long long>(npts[dim_axis])) n1_valid = false;
            }
            if (n1_valid) {
                if (dim_axis == 0) neighbor1_idx.i1 = static_cast<size_t>(n1_comp_val);
                else if (dim_axis == 1) neighbor1_idx.i2 = static_cast<size_t>(n1_comp_val);
                else neighbor1_idx.i3 = static_cast<size_t>(n1_comp_val);
            }

            real_t T1 = INF;
            if (n1_valid && m_fmm_states[get_flat_state_index(neighbor1_idx)] == FMMState::KNOWN) {
                T1 = m_traveltime_field.get_value(neighbor1_idx);
            }

            real_t T2 = INF; 
            if (n1_valid && m_fmm_states[get_flat_state_index(neighbor1_idx)] == FMMState::KNOWN) { 
                long long n2_comp_val = target_comp_val + 2 * dir;
                 if (is_periodic[dim_axis]) {
                    n2_comp_val = (n2_comp_val % static_cast<long long>(npts[dim_axis]) + static_cast<long long>(npts[dim_axis])) % static_cast<long long>(npts[dim_axis]);
                    n2_valid = true; 
                } else {
                    if (n2_comp_val >= 0 && n2_comp_val < static_cast<long long>(npts[dim_axis])) n2_valid = true;
                }

                if (n2_valid) {
                    if (dim_axis == 0) neighbor2_idx.i1 = static_cast<size_t>(n2_comp_val);
                    else if (dim_axis == 1) neighbor2_idx.i2 = static_cast<size_t>(n2_comp_val);
                    else neighbor2_idx.i3 = static_cast<size_t>(n2_comp_val);
                    
                    if (m_fmm_states[get_flat_state_index(neighbor2_idx)] == FMMState::KNOWN) {
                         T2 = m_traveltime_field.get_value(neighbor2_idx);
                    }
                }
            }
            
            int candidate_slot = (dir == -1) ? 0 : 1;
            if (n2_valid && T2 <= T1 && T2 != INF) { 
                T_candidates[candidate_slot] = (4.0 * T1 - T2) / 3.0; 
                is_second_order_candidate[candidate_slot] = true;
                T1_for_candidate[candidate_slot] = T1;
                T2_for_candidate[candidate_slot] = T2;
                                                                     
            } else if (T1 != INF) { 
                T_candidates[candidate_slot] = T1; 
                is_second_order_candidate[candidate_slot] = false;
                T1_for_candidate[candidate_slot] = T1;
                T2_for_candidate[candidate_slot] = INF; // Not used for 1st order
            }
        } 
        
        // Python code's selection logic (fdu[0] > -fdu[1]) is essentially choosing the upwind value
        // that minimizes the resulting T_target. This is equivalent to solving the Eikonal
        // update equation using T_candidates[0] and T_candidates[1] separately (if valid)
        // and picking the one that gives the smallest valid T_target.
        // Or, more simply, for each dimension, we pick the smaller of T_candidates[0] and T_candidates[1]
        // to form the ( (T_target - T_eff_dim) / h_eff_dim )^2 term.
        
        real_t T_eff_dim = INF;
        bool chosen_is_second_order = false;
        real_t T1_chosen = INF;
        real_t T2_chosen = INF;

        if (T_candidates[0] < T_candidates[1]) {
            T_eff_dim = T_candidates[0];
            chosen_is_second_order = is_second_order_candidate[0];
            T1_chosen = T1_for_candidate[0];
            T2_chosen = T2_for_candidate[0];
        } else if (T_candidates[1] < INF) { // T_candidates[1] <= T_candidates[0] or T_cand[0] was INF
            T_eff_dim = T_candidates[1];
            chosen_is_second_order = is_second_order_candidate[1];
            T1_chosen = T1_for_candidate[1];
            T2_chosen = T2_for_candidate[1];
        }


        if (T_eff_dim == INF) { 
            continue;
        }

        if (chosen_is_second_order) { 
            aa[dim_axis] = 9.0 / (4.0 * h_dim * h_dim);
            bb[dim_axis] = (6.0 * T2_chosen - 24.0 * T1_chosen) / (4.0 * h_dim * h_dim);
            cc[dim_axis] = (T2_chosen * T2_chosen - 8.0 * T2_chosen * T1_chosen + 16.0 * T1_chosen * T1_chosen) / (4.0 * h_dim * h_dim);
        } else { 
            aa[dim_axis] = 1.0 / (h_dim * h_dim);
            bb[dim_axis] = -2.0 * T1_chosen / (h_dim * h_dim); // T1_chosen is T_eff_dim here
            cc[dim_axis] = T1_chosen * T1_chosen / (h_dim * h_dim);
        }
    } 

    real_t A_quad = aa[0] + aa[1] + aa[2];
    real_t B_quad = bb[0] + bb[1] + bb[2];
    real_t C_quad_base = cc[0] + cc[1] + cc[2];

    if (std::abs(A_quad) < 1e-12) { // Adjusted tolerance slightly, was 1e-9
        // This case means not enough dimensions contributed with a T^2 term.
        // Could be 1D update (linear equation for T) or no valid update.
        // If, for example, only one dimension (say dim 0) contributes with first order:
        // (T - T1_0)^2 / h0^2 = s^2  => T = T1_0 + s * h0.
        // The quadratic solver might not be appropriate.
        // However, Ckonal's logic with aa,bb,cc seems to handle this implicitly if one or two aa[dim] are zero.
        // If all aa[dim] are zero, A_quad is zero.
        // Then it's B_quad * T + C_quad = 0.  T = -C_quad / B_quad, if B_quad != 0.
        // This scenario needs to be robustly handled by solve_quadratic or here.
        // For now, if A_quad is essentially zero, let's assume no valid quadratic update.
        return INF; 
    }

    real_t vel_at_target = m_velocity_field.get_value(target_node_idx);
    if (vel_at_target <= 1e-9) return INF; // Velocity zero or negative, slowness infinite.
    real_t slowness_sq = 1.0 / (vel_at_target * vel_at_target);
    real_t C_quad = C_quad_base - slowness_sq;

    std::vector<real_t> roots = utils::solve_quadratic(A_quad, B_quad, C_quad);

    real_t updated_T = INF;
    if (!roots.empty()) {
        if (roots.size() == 1) updated_T = roots[0];
        else if (roots.size() == 2) updated_T = std::max(roots[0], roots[1]);
        
        // Basic causality check: new time should not be less than any of the
        // T1_chosen values that contributed to non-zero aa[dim_axis].
        // This is a simplification; a full causality check against all known neighbors is more robust.
        // For now, this matches the general FMM approach of taking the larger root.
    }
    
    return updated_T;
}


// --- PointSourceSolver Implementation ---
PointSourceSolver::PointSourceSolver(CoordSys far_field_coord_sys)
    : EikonalSolver(far_field_coord_sys), 
      m_n_rho(DEFAULT_N_RHO),
      m_n_theta(DEFAULT_N_THETA),
      m_n_phi(DEFAULT_N_PHI),
      m_drho(0.0), 
      m_source_loc_set(false) {
    m_near_field_solver = std::make_unique<EikonalSolver>(CoordSys::SPHERICAL);
    std::cout << "PointSourceSolver: Constructor called. Far-field cs: "
              << to_string(far_field_coord_sys) << std::endl;
}

void PointSourceSolver::set_source_location(const Point3D& loc) {
    m_source_loc_far_field = loc;
    m_source_loc_set = true;
    std::cout << "PointSourceSolver: Source location set to ("
              << loc[0] << ", " << loc[1] << ", " << loc[2] << ") in far-field coordinates." << std::endl;
}
Point3D PointSourceSolver::get_source_location() const {
    if (!m_source_loc_set) throw std::runtime_error("Source location not set for PointSourceSolver.");
    return m_source_loc_far_field;
}

void PointSourceSolver::set_near_field_radial_nodes(uint_t n_rho) { m_n_rho = n_rho; }
void PointSourceSolver::set_near_field_polar_nodes(uint_t n_theta) { m_n_theta = n_theta; }
void PointSourceSolver::set_near_field_azimuthal_nodes(uint_t n_phi) { m_n_phi = n_phi; }
void PointSourceSolver::set_near_field_radial_interval(real_t drho) { m_drho = drho; }


void PointSourceSolver::initialize_near_field_grid_parameters() {
    std::cout << "PointSourceSolver Debug: Initializing near-field grid parameters..." << std::endl;
    if (m_velocity_field.get_total_nodes() == 0) {
        throw std::runtime_error("Far-field velocity grid must be configured for PointSourceSolver.");
    }

    const auto& far_intervals = m_velocity_field.get_node_intervals();
    if (m_drho <= 1e-9) { 
        if (m_coord_sys == CoordSys::CARTESIAN) { 
            real_t min_interval = far_intervals[0];
            if (m_velocity_field.get_npts()[1] > 1 && far_intervals[1] > 1e-9) min_interval = std::min(min_interval, far_intervals[1]);
            else if (m_velocity_field.get_npts()[1] <=1 && min_interval < 1e-9 && far_intervals[1] > 1e-9) min_interval = far_intervals[1]; // if first was zero

            if (m_velocity_field.get_npts()[2] > 1 && far_intervals[2] > 1e-9) min_interval = std::min(min_interval, far_intervals[2]);
            else if (m_velocity_field.get_npts()[2] <=1 && min_interval < 1e-9 && far_intervals[2] > 1e-9) min_interval = far_intervals[2];
            
            m_drho = min_interval / 8.0;
        } else { 
            m_drho = far_intervals[0] / 8.0; 
        }
        std::cout << "PointSourceSolver Debug: Calculated default drho = " << m_drho << std::endl;
    }
    if (m_drho <= 1e-9) {
        m_drho = 0.1; 
        std::cout << "PointSourceSolver Warning: drho was too small or far-field intervals were zero, set to fallback " << m_drho << std::endl;
    }

    m_near_field_solver->get_velocity_field().set_min_coords({m_drho, 0.0, 0.0}); 

    real_t d_theta_near = (m_n_theta > 1) ? (PI / static_cast<real_t>(m_n_theta - 1)) : PI; 
    real_t d_phi_near = (m_n_phi > 0) ? (TWO_PI / static_cast<real_t>(m_n_phi)) : TWO_PI; 
    if (m_n_theta <=1) d_theta_near = PI; // Avoid interval being NaN or Inf if n_theta is 1
    if (m_n_phi <=1) d_phi_near = TWO_PI; // Similar for phi

    m_near_field_solver->get_velocity_field().set_node_intervals({m_drho, d_theta_near, d_phi_near});
    m_near_field_solver->get_velocity_field().set_npts({m_n_rho, m_n_theta, m_n_phi});

    std::cout << "PointSourceSolver Debug: Near-field grid: npts=("
              << m_n_rho << "," << m_n_theta << "," << m_n_phi << "), intervals=("
              << m_drho << "," << d_theta_near << "," << d_phi_near << "), min_coords=("
              << m_drho << ",0,0)." << std::endl;

    m_near_field_solver->init_solver_state(); 
    std::cout << "PointSourceSolver Debug: Near-field solver state initialized." << std::endl;
}


void PointSourceSolver::init_solver_state() {
    std::cout << "PointSourceSolver: Starting init_solver_state..." << std::endl;
    if (!m_source_loc_set) {
        throw std::runtime_error("Source location must be set before initializing PointSourceSolver.");
    }
    EikonalSolver::init_solver_state(); 
    std::cout << "PointSourceSolver: Far-field solver state initialized by base class." << std::endl;

    initialize_near_field_grid_parameters();

    std::cout << "PointSourceSolver: Interpolating far-field velocity to near-field..." << std::endl;
    interpolate_far_field_velocity_to_near_field();
    std::cout << "PointSourceSolver: Far-field velocity interpolation to near-field done." << std::endl;

    std::cout << "PointSourceSolver: Initializing near-field source band..." << std::endl;
    initialize_near_field_source_band();
    std::cout << "PointSourceSolver: Near-field source band initialized." << std::endl;

    m_is_initialized = true;
    std::cout << "PointSourceSolver: init_solver_state finished. Solver is initialized." << std::endl;
}


void PointSourceSolver::interpolate_far_field_velocity_to_near_field() {
    ScalarField3D& near_vv = m_near_field_solver->get_velocity_field();
    const ScalarField3D& far_vv = this->get_velocity_field();

    std::cout << "PointSourceSolver Debug: Interpolating far-field Vv to near-field Vv." << std::endl;
    std::cout << "PointSourceSolver Debug: Far-field Vv npts: ("
              << far_vv.get_npts()[0] << "," << far_vv.get_npts()[1] << "," << far_vv.get_npts()[2] << ")" << std::endl;
    std::cout << "PointSourceSolver Debug: Near-field Vv npts: ("
              << near_vv.get_npts()[0] << "," << near_vv.get_npts()[1] << "," << near_vv.get_npts()[2] << ")" << std::endl;
    std::cout << "PointSourceSolver Debug: Source location (far-field coords): ("
              << m_source_loc_far_field[0] << "," << m_source_loc_far_field[1] << "," << m_source_loc_far_field[2] << ")" << std::endl;


    Point3D source_xyz_abs_in_far_frame; 
    if (this->m_coord_sys == CoordSys::SPHERICAL) { 
        source_xyz_abs_in_far_frame = transformations::sph_to_xyz(m_source_loc_far_field, {0.0, 0.0, 0.0});
    } else { 
        source_xyz_abs_in_far_frame = m_source_loc_far_field;
    }
     std::cout << "PointSourceSolver Debug: Source location (abs Cartesian in far frame): ("
              << source_xyz_abs_in_far_frame[0] << "," << source_xyz_abs_in_far_frame[1] << "," << source_xyz_abs_in_far_frame[2] << ")" << std::endl;


    size_t nan_count = 0;
    for (uint_t i = 0; i < near_vv.get_npts()[0]; ++i) { 
        for (uint_t j = 0; j < near_vv.get_npts()[1]; ++j) { 
            for (uint_t k = 0; k < near_vv.get_npts()[2]; ++k) { 
                Index3D near_idx(i,j,k);
                Point3D near_sph_coords = near_vv.get_node_coords(near_idx); 

                Point3D local_xyz_offset = transformations::sph_to_xyz(near_sph_coords, {0.0, 0.0, 0.0});

                Point3D node_abs_xyz_in_far_frame = {
                    source_xyz_abs_in_far_frame[0] + local_xyz_offset[0],
                    source_xyz_abs_in_far_frame[1] + local_xyz_offset[1],
                    source_xyz_abs_in_far_frame[2] + local_xyz_offset[2]
                };

                Point3D far_field_query_coords;
                if (this->m_coord_sys == CoordSys::SPHERICAL) { 
                    far_field_query_coords = transformations::xyz_to_sph(node_abs_xyz_in_far_frame, {0.0, 0.0, 0.0}, false);
                } else { 
                    far_field_query_coords = node_abs_xyz_in_far_frame;
                }
                
                real_t interpolated_vel = far_vv.interpolate_value(far_field_query_coords, NaN);
                if (std::isnan(interpolated_vel)) {
                    nan_count++;
                    near_vv.set_value(near_idx, 1e-6); 
                } else if (interpolated_vel <= 1e-9) {
                    near_vv.set_value(near_idx, 1e-6); 
                }
                else {
                    near_vv.set_value(near_idx, interpolated_vel);
                }

                if (i==0 && j==0 && k < 3) { 
                     std::cout << "  NearVV interp: near_idx(0,0," << k << ") near_sph(" << near_sph_coords[0] << "," << near_sph_coords[1] << "," << near_sph_coords[2] << ")"
                               << " -> local_xyz_off(" << local_xyz_offset[0] << "," << local_xyz_offset[1] << "," << local_xyz_offset[2] << ")"
                               << " -> node_abs_xyz(" << node_abs_xyz_in_far_frame[0] << "," << node_abs_xyz_in_far_frame[1] << "," << node_abs_xyz_in_far_frame[2] << ")"
                               << " -> far_query(" << far_field_query_coords[0] << "," << far_field_query_coords[1] << "," << far_field_query_coords[2] << ")"
                               << " vel=" << near_vv.get_value(near_idx) << std::endl;
                }
            }
        }
    }
    if (nan_count > 0) {
        std::cout << "PointSourceSolver Warning: " << nan_count << " points in near-field Vv were NaN after interpolation from far-field, filled with tiny positive." << std::endl;
    }
}


void PointSourceSolver::initialize_near_field_source_band() {
    ScalarField3D& near_vv = m_near_field_solver->get_velocity_field();
    const auto& near_npts = near_vv.get_npts();
    std::cout << "PointSourceSolver Debug: Initializing near-field source band with drho = " << m_drho << std::endl;
    int sources_added_to_near_field = 0;

    if (near_npts[0] == 0) { 
        std::cerr << "PointSourceSolver Error: Near field npts[0] is zero in initialize_near_field_source_band." << std::endl;
        return;
    }

    for (uint_t it = 0; it < near_npts[1]; ++it) {
        for (uint_t ip = 0; ip < near_npts[2]; ++ip) {
            Index3D source_shell_idx(0, it, ip); 
            real_t vel_at_node = near_vv.get_value(source_shell_idx);
            if (vel_at_node > 1e-9) { 
                real_t time_to_node = m_drho / vel_at_node;
                m_near_field_solver->add_source_point(source_shell_idx, time_to_node);
                sources_added_to_near_field++;
                 if (it == 0 && ip < 3 && sources_added_to_near_field < 5) { 
                    std::cout << "  Near Source Add: idx(0," << it << "," << ip << ") vel=" << vel_at_node << " time=" << time_to_node << std::endl;
                }
            } else {
                 if (it == 0 && ip < 3) {
                    std::cout << "  Near Source Skip: idx(0," << it << "," << ip << ") vel=" << vel_at_node << std::endl;
                 }
            }
        }
    }
    std::cout << "PointSourceSolver Debug: Added " << sources_added_to_near_field << " sources to near-field trial heap." << std::endl;
    if (sources_added_to_near_field == 0 && near_vv.get_total_nodes() > 0) {
        std::cerr << "PointSourceSolver Warning: No sources added to near-field heap. Check velocities on first rho shell." << std::endl;
    }
}

void PointSourceSolver::interpolate_near_field_traveltime_to_far_field() {
    const ScalarField3D& near_tt = m_near_field_solver->get_traveltime_field();
    ScalarField3D& far_tt = this->get_traveltime_field();
    std::cout << "PointSourceSolver Debug: Interpolating near-field Tt to far-field Tt." << std::endl;

    Point3D source_xyz_abs_in_far_frame;
    if (this->m_coord_sys == CoordSys::SPHERICAL) {
        source_xyz_abs_in_far_frame = transformations::sph_to_xyz(m_source_loc_far_field, {0.0, 0.0, 0.0});
    } else {
        source_xyz_abs_in_far_frame = m_source_loc_far_field;
    }

    size_t updated_far_tt_nodes = 0;
    for (uint_t i = 0; i < far_tt.get_npts()[0]; ++i) {
        for (uint_t j = 0; j < far_tt.get_npts()[1]; ++j) {
            for (uint_t k = 0; k < far_tt.get_npts()[2]; ++k) {
                Index3D far_idx(i,j,k);
                Point3D far_coords_abs_native = far_tt.get_node_coords(far_idx); 

                Point3D far_node_abs_xyz_in_far_frame; 
                if (this->m_coord_sys == CoordSys::SPHERICAL) {
                    far_node_abs_xyz_in_far_frame = transformations::sph_to_xyz(far_coords_abs_native, {0.0, 0.0, 0.0});
                } else {
                    far_node_abs_xyz_in_far_frame = far_coords_abs_native;
                }
                
                Point3D coords_rel_to_source_cartesian = {
                    far_node_abs_xyz_in_far_frame[0] - source_xyz_abs_in_far_frame[0],
                    far_node_abs_xyz_in_far_frame[1] - source_xyz_abs_in_far_frame[1],
                    far_node_abs_xyz_in_far_frame[2] - source_xyz_abs_in_far_frame[2]
                };
                
                Point3D near_field_query_sph_coords = transformations::xyz_to_sph(coords_rel_to_source_cartesian, {0.0, 0.0, 0.0}, false);
                
                real_t interpolated_tt = near_tt.interpolate_value(near_field_query_sph_coords, INF);

                if (interpolated_tt < far_tt.get_value(far_idx)) {
                    far_tt.set_value(far_idx, interpolated_tt);
                    updated_far_tt_nodes++;
                     if (i < 2 && j < 2 && k < 2 && updated_far_tt_nodes < 5) { 
                        std::cout << "  FarTT Update: far_idx(" << i << "," << j << "," << k <<") "
                                  << "far_native(" << far_coords_abs_native[0] << "," << far_coords_abs_native[1] << "," << far_coords_abs_native[2] << ") "
                                  << "-> rel_cart(" << coords_rel_to_source_cartesian[0] << "," << coords_rel_to_source_cartesian[1] << "," << coords_rel_to_source_cartesian[2] << ") "
                                  << "-> near_query_sph(" << near_field_query_sph_coords[0] << "," << near_field_query_sph_coords[1] << "," << near_field_query_sph_coords[2] << ") "
                                  << "tt=" << interpolated_tt << std::endl;
                    }
                }
            }
        }
    }
    std::cout << "PointSourceSolver Debug: Updated " << updated_far_tt_nodes << " nodes in far-field Tt from near-field." << std::endl;
}

void PointSourceSolver::initialize_far_field_source_band_from_near_field() {
    ScalarField3D& far_tt = this->get_traveltime_field(); 
    int sources_added_to_far_field = 0;
    std::cout << "PointSourceSolver Debug: Initializing far-field source band from near-field interpolation." << std::endl;

    for (uint_t i = 0; i < far_tt.get_npts()[0]; ++i) {
        for (uint_t j = 0; j < far_tt.get_npts()[1]; ++j) {
            for (uint_t k = 0; k < far_tt.get_npts()[2]; ++k) {
                Index3D far_idx(i,j,k);
                if (far_tt.get_value(far_idx) < (INF - 1.0)) { 
                    size_t flat_idx = get_flat_state_index(far_idx);
                    if (m_fmm_states[flat_idx] == FMMState::UNKNOWN) {
                        m_fmm_states[flat_idx] = FMMState::TRIAL;
                        m_trial_heap->push(far_idx); 
                        sources_added_to_far_field++;
                    }
                }
            }
        }
    }
    std::cout << "PointSourceSolver Debug: Added " << sources_added_to_far_field << " sources to far-field trial heap from near-field results." << std::endl;
     if (sources_added_to_far_field == 0 && far_tt.get_total_nodes() > 0) {
        std::cerr << "PointSourceSolver Warning: No sources added to far-field heap from near-field interpolation. Far-field solve might fail." << std::endl;
    }
}


bool PointSourceSolver::solve() {
    if (!m_is_initialized) {
        std::cerr << "PointSourceSolver Error: Not initialized!" << std::endl;
        throw std::runtime_error("PointSourceSolver not initialized. Call init_solver_state() first.");
    }
    std::cout << "PointSourceSolver: Starting solve procedure..." << std::endl;

    std::cout << "PointSourceSolver: Solving near-field..." << std::endl;
    bool near_success = m_near_field_solver->solve(); // This calls EikonalSolver::solve() on m_near_field_solver
    if (!near_success) {
        std::cerr << "PointSourceSolver Error: Near-field solve failed or had no sources." << std::endl;
        // Use the new public getter method for m_near_field_solver
        if (m_near_field_solver->is_trial_heap_empty() &&
            m_near_field_solver->get_traveltime_field().get_total_nodes() > 0) {
             std::cerr << "PointSourceSolver Detail: Near-field trial heap was empty or became empty immediately." << std::endl;
        }
        return false;
    }
    std::cout << "PointSourceSolver: Near-field solve completed. Near-field heap size after solve: "
              << m_near_field_solver->get_trial_heap_size() << std::endl; // Use getter

    std::cout << "PointSourceSolver: Interpolating near-field Tt to far-field Tt..." << std::endl;
    interpolate_near_field_traveltime_to_far_field();

    std::cout << "PointSourceSolver: Initializing far-field source band from near-field interpolation..." << std::endl;
    initialize_far_field_source_band_from_near_field();
    // For *this* object's (far-field) heap, we can still access protected m_trial_heap directly
    // or use the public getter for consistency.
    if (this->is_trial_heap_empty()) { // Using public getter
         std::cerr << "PointSourceSolver Warning: Far-field trial heap is EMPTY after band initialization from near-field." << std::endl;
    } else {
         std::cout << "PointSourceSolver: Far-field source band initialized with " << this->get_trial_heap_size() << " trial points." << std::endl;
    }

    std::cout << "PointSourceSolver: Solving far-field (base class EikonalSolver::solve())..." << std::endl;
    bool far_success = EikonalSolver::solve(); // Calls base class solve on *this* object's far-field data
    if (!far_success) {
        std::cerr << "PointSourceSolver Error: Far-field solve failed." << std::endl;
         if (this->is_trial_heap_empty() && this->get_traveltime_field().get_total_nodes() > 0) { // Using public getter
             std::cerr << "PointSourceSolver Detail: Far-field trial heap was empty or became empty immediately during far-solve." << std::endl;
        }
        return false;
    }
    std::cout << "PointSourceSolver: Far-field solve completed. Far-field heap size after solve: "
              << this->get_trial_heap_size() << std::endl; // Using public getter
    std::cout << "PointSourceSolver: Solve procedure finished successfully." << std::endl;
    return true;
}

} // namespace Ckonal