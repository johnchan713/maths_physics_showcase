#ifndef PHYSICS_ADVANCED_HPP
#define PHYSICS_ADVANCED_HPP

/**
 * @file physics_advanced.hpp
 * @brief Central header for all advanced physics modules
 *
 * This is the main entry point for the advanced physics library.
 * Include this file to access all advanced features.
 *
 * Categories:
 * 1. Classical Mechanics (Hamiltonian, Phase Space, Liouville)
 * 2. Cosmology (Friedmann equations, early universe, CMB)
 * 3. Fluid Dynamics (Turbulence, boundary layers, compressible flow)
 * 4. Gauge Theory (Gauge invariance, Higgs mechanism, CP violation)
 * 5. Quantum Field Theory (Particle physics, QGP, supersymmetry)
 * 6. Quantum Mechanics (Advanced QM, quantum chemistry, foundations)
 * 7. Specialized (Loop quantum gravity, operator algebras, nuclear physics)
 *
 * Requirements:
 * - C++17 or later
 * - Standard library only (header-only implementation)
 */

// Category 1: Classical Mechanics
#include "classical_hamiltonian.hpp"
#include "classical_phase_space.hpp"
#include "classical_liouville.hpp"

// Category 2: Cosmology
#include "cosmology_friedmann_equations.hpp"
#include "cosmology_expanding_universe.hpp"
#include "cosmology_energy_density.hpp"
#include "cosmology_early_universe.hpp"

// Category 3: Fluid Dynamics
#include "fluid_dynamics_governing_equations.hpp"
#include "fluid_dynamics_dimensionless_numbers.hpp"
#include "fluid_dynamics_flow_types.hpp"
#include "fluid_dynamics_boundary_layer.hpp"
#include "fluid_dynamics_turbulence.hpp"
#include "fluid_dynamics_vorticity.hpp"
#include "fluid_dynamics_compressible_flow.hpp"

// Category 4: Gauge Theory
#include "gauge_theory_gauge_invariance.hpp"
#include "gauge_theory_higgs_mechanism.hpp"
#include "gauge_theory_running_couplings.hpp"
#include "gauge_theory_symmetries.hpp"
#include "gauge_theory_helicity.hpp"
#include "gauge_theory_cp_violation_kaons.hpp"

// Category 5: Quantum Field Theory
#include "qft_particle_physics.hpp"
#include "qft_spin_statistics.hpp"
#include "qft_antiparticles.hpp"
#include "qft_interactions.hpp"
#include "qft_cross_sections.hpp"
#include "qft_decays.hpp"
#include "qft_quark_gluon_plasma.hpp"
#include "qft_supersymmetry.hpp"

// Category 6: Advanced Quantum Mechanics
#include "advanced_quantum_mechanics.hpp"
#include "quantum_chemistry.hpp"
#include "quantum_foundations.hpp"
#include "relativistic_quantum_mechanics.hpp"

// Category 7: Specialized Topics
#include "loop_quantum_gravity.hpp"
#include "operator_algebras.hpp"
#include "nuclear_physics.hpp"

namespace physics::advanced {

struct Version {
    static constexpr int MAJOR = 2;
    static constexpr int MINOR = 0;
    static constexpr int PATCH = 0;
};

} // namespace physics::advanced

#endif // PHYSICS_ADVANCED_HPP
