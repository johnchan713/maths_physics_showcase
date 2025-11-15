#ifndef PHYSICS_RELATIVISTIC_QUANTUM_MECHANICS_HPP
#define PHYSICS_RELATIVISTIC_QUANTUM_MECHANICS_HPP

#include <vector>
#include <complex>
#include <functional>
#include <cmath>
#include <array>
#include <algorithm>

namespace physics {
namespace relativistic_quantum {

using Complex = std::complex<double>;

/**
 * @brief Physical constants for relativistic quantum mechanics
 */
namespace constants {
    constexpr double h = 6.62607015e-34;        // Planck constant (J·s)
    constexpr double hbar = 1.054571817e-34;    // Reduced Planck constant
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double e = 1.602176634e-19;       // Elementary charge (C)
    constexpr double m_e = 9.1093837015e-31;    // Electron mass (kg)
    constexpr double a_0 = 5.29177210903e-11;   // Bohr radius (m)
    constexpr double alpha = 7.2973525693e-3;   // Fine structure constant ≈ 1/137
    constexpr double mu_B = 9.2740100783e-24;   // Bohr magneton (J/T)
    constexpr double g_e = 2.00231930436256;    // Electron g-factor
}

/**
 * @brief Degenerate Position Eigenstates
 *
 * Handling degeneracy in quantum systems
 */
class DegeneratePositionEigenstates {
public:
    /**
     * @brief Degeneracy of energy level
     *
     * For hydrogen: g_n = n² (including spin: 2n²)
     */
    static int degeneracy_hydrogen(int n, bool include_spin = false) {
        int g = n * n;
        if (include_spin) {
            g *= 2;
        }
        return g;
    }

    /**
     * @brief Degeneracy for 3D harmonic oscillator
     *
     * g_N = (N+1)(N+2)/2
     */
    static int degeneracy_3d_oscillator(int N) {
        return (N + 1) * (N + 2) / 2;
    }

    /**
     * @brief Degenerate perturbation theory
     *
     * Need to diagonalize W within degenerate subspace
     */
    struct DegenerateSubspace {
        int dimension;
        std::vector<std::vector<Complex>> perturbation_matrix;
        std::vector<double> eigenvalues;
        std::vector<std::vector<Complex>> eigenvectors;
    };

    /**
     * @brief Lifted degeneracy
     *
     * Perturbation splits degenerate levels
     */
    static std::vector<double> split_degenerate_level(
        double E_0,
        const std::vector<double>& perturbation_eigenvalues) {

        std::vector<double> split_levels;
        for (double delta_E : perturbation_eigenvalues) {
            split_levels.push_back(E_0 + delta_E);
        }
        return split_levels;
    }

    /**
     * @brief Good quantum numbers
     *
     * Quantum numbers that commute with H
     */
    static bool is_good_quantum_number(
        const std::string& operator_name,
        const std::string& hamiltonian_symmetry) {

        // If [O, H] = 0, then O is a good quantum number
        return operator_name == hamiltonian_symmetry;
    }

    /**
     * @brief Accidental degeneracy
     *
     * Degeneracy not required by symmetry (e.g., hydrogen l-degeneracy)
     */
    static bool is_accidental_degeneracy(
        const std::string& system,
        int n,
        int l1,
        int l2) {

        if (system == "hydrogen") {
            // In hydrogen, different l values are degenerate (accidental)
            return l1 != l2 && l1 < n && l2 < n;
        }
        return false;
    }
};

/**
 * @brief Spin-Half Particles
 *
 * Comprehensive treatment of spin-1/2 systems
 */
class SpinHalfParticles {
public:
    /**
     * @brief Spin angular momentum operators
     *
     * S_x, S_y, S_z with [S_i, S_j] = iℏε_ijk S_k
     */
    struct SpinOperators {
        // Pauli matrices (in units of ℏ/2)
        static std::array<std::array<Complex, 2>, 2> sigma_x() {
            return {{
                {Complex(0, 0), Complex(1, 0)},
                {Complex(1, 0), Complex(0, 0)}
            }};
        }

        static std::array<std::array<Complex, 2>, 2> sigma_y() {
            return {{
                {Complex(0, 0), Complex(0, -1)},
                {Complex(0, 1), Complex(0, 0)}
            }};
        }

        static std::array<std::array<Complex, 2>, 2> sigma_z() {
            return {{
                {Complex(1, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(-1, 0)}
            }};
        }

        // S_i = (ℏ/2) σ_i
        static std::array<std::array<Complex, 2>, 2> S_x() {
            auto sigma = sigma_x();
            double factor = constants::hbar / 2.0;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    sigma[i][j] *= factor;
                }
            }
            return sigma;
        }
    };

    /**
     * @brief Spin eigenstates
     *
     * |↑⟩ = (1, 0)ᵀ, |↓⟩ = (0, 1)ᵀ
     */
    static std::array<Complex, 2> spin_up() {
        return {Complex(1, 0), Complex(0, 0)};
    }

    static std::array<Complex, 2> spin_down() {
        return {Complex(0, 0), Complex(1, 0)};
    }

    /**
     * @brief General spin state
     *
     * |χ⟩ = cos(θ/2)|↑⟩ + e^(iφ)sin(θ/2)|↓⟩
     */
    static std::array<Complex, 2> general_spin_state(double theta, double phi) {
        return {
            Complex(std::cos(theta / 2.0), 0),
            Complex(std::sin(theta / 2.0) * std::cos(phi),
                   std::sin(theta / 2.0) * std::sin(phi))
        };
    }

    /**
     * @brief Spin expectation values
     *
     * ⟨S_x⟩, ⟨S_y⟩, ⟨S_z⟩
     */
    static std::array<double, 3> spin_expectation(
        const std::array<Complex, 2>& state) {

        // ⟨S_z⟩ = (ℏ/2)(|α|² - |β|²)
        double S_z = (constants::hbar / 2.0) *
                    (std::norm(state[0]) - std::norm(state[1]));

        // ⟨S_x⟩ = (ℏ/2)Re(α*β)
        double S_x = (constants::hbar / 2.0) *
                    (std::real(std::conj(state[0]) * state[1]) +
                     std::real(std::conj(state[1]) * state[0]));

        // ⟨S_y⟩ = (ℏ/2)Im(α*β)
        double S_y = (constants::hbar / 2.0) *
                    (std::imag(std::conj(state[0]) * state[1]) -
                     std::imag(std::conj(state[1]) * state[0]));

        return {S_x, S_y, S_z};
    }

    /**
     * @brief Spin precession in magnetic field
     *
     * ω = -γB (Larmor frequency)
     * γ = g_e μ_B/ℏ (gyromagnetic ratio)
     */
    static double larmor_frequency(double magnetic_field) {
        double gamma = constants::g_e * constants::mu_B / constants::hbar;
        return gamma * magnetic_field;
    }

    /**
     * @brief Time evolution of spin state
     *
     * |χ(t)⟩ = e^(-iH t/ℏ)|χ(0)⟩
     */
    static std::array<Complex, 2> evolve_spin_state(
        const std::array<Complex, 2>& initial_state,
        double magnetic_field,
        double time) {

        double omega = larmor_frequency(magnetic_field);
        double phase = omega * time / 2.0;

        Complex phase_up = std::exp(Complex(0, -phase));
        Complex phase_down = std::exp(Complex(0, phase));

        return {
            initial_state[0] * phase_up,
            initial_state[1] * phase_down
        };
    }

    /**
     * @brief Spin-1/2 density matrix
     *
     * ρ = (1/2)(I + r⃗·σ⃗)
     *
     * r⃗: Bloch vector
     */
    static std::array<std::array<Complex, 2>, 2> density_matrix(
        double r_x,
        double r_y,
        double r_z) {

        return {{
            {Complex(0.5 * (1.0 + r_z), 0),
             Complex(0.5 * (r_x - r_y), 0)},
            {Complex(0.5 * (r_x + r_y), 0),
             Complex(0.5 * (1.0 - r_z), 0)}
        }};
    }

    /**
     * @brief Purity of spin state
     *
     * Pure state: Tr(ρ²) = 1
     * Mixed state: Tr(ρ²) < 1
     */
    static double purity(double r_x, double r_y, double r_z) {
        return r_x * r_x + r_y * r_y + r_z * r_z;
    }
};

/**
 * @brief Spin Magnetic Moment and Stern-Gerlach Experiment
 *
 * Interaction of spin with magnetic field
 */
class SpinMagneticMoment {
public:
    /**
     * @brief Magnetic moment operator
     *
     * μ⃗ = -g_e (μ_B/ℏ) S⃗
     */
    static double magnetic_moment_z(int m_s_sign) {
        return -constants::g_e * constants::mu_B * (m_s_sign * 0.5);
    }

    /**
     * @brief Interaction energy with magnetic field
     *
     * H = -μ⃗·B⃗ = g_e μ_B m_s B
     */
    static double interaction_energy(int m_s_sign, double B_field) {
        return constants::g_e * constants::mu_B * (m_s_sign * 0.5) * B_field;
    }

    /**
     * @brief Stern-Gerlach force
     *
     * F_z = μ_z ∂B_z/∂z
     */
    static double stern_gerlach_force(int m_s_sign, double gradient) {
        double mu_z = magnetic_moment_z(m_s_sign);
        return mu_z * gradient;
    }

    /**
     * @brief Beam deflection in Stern-Gerlach apparatus
     *
     * Δz = (1/2)(F/m)t²
     */
    static double beam_deflection(
        int m_s_sign,
        double gradient,
        double apparatus_length,
        double beam_velocity) {

        double force = stern_gerlach_force(m_s_sign, gradient);
        double time = apparatus_length / beam_velocity;
        return 0.5 * (force / constants::m_e) * time * time;
    }

    /**
     * @brief Sequential Stern-Gerlach experiments
     *
     * Demonstrates quantum measurement and collapse
     */
    struct SequentialSG {
        enum class Orientation { X, Y, Z };

        static double probability_after_measurement(
            Orientation first_axis,
            Orientation second_axis,
            bool first_result_up,
            bool second_result_up) {

            // If same axis: probability is 1 if same result, 0 otherwise
            if (first_axis == second_axis) {
                return (first_result_up == second_result_up) ? 1.0 : 0.0;
            }

            // If perpendicular axes: probability is 1/2
            return 0.5;
        }
    };

    /**
     * @brief Landé g-factor
     *
     * g_J = 1 + [J(J+1) - L(L+1) + S(S+1)] / [2J(J+1)]
     */
    static double lande_g_factor(int J, int L, int S) {
        if (J == 0) return 0.0;

        double numerator = J * (J + 1) - L * (L + 1) + S * (S + 1);
        double denominator = 2.0 * J * (J + 1);

        return 1.0 + numerator / denominator;
    }

    /**
     * @brief Magnetic moment for atom
     *
     * μ_J = -g_J μ_B √(J(J+1))
     */
    static double atomic_magnetic_moment(int J, int L, int S) {
        double g_J = lande_g_factor(J, L, S);
        return g_J * constants::mu_B * std::sqrt(J * (J + 1.0));
    }
};

/**
 * @brief Spin-Orbit Coupling
 *
 * Fine structure from relativistic corrections
 */
class SpinOrbitCoupling {
public:
    /**
     * @brief Spin-orbit Hamiltonian
     *
     * H_SO = (1/(2m²c²r)) (dV/dr) L⃗·S⃗
     */
    static double spin_orbit_energy(
        int n,
        int l,
        int j,
        int Z = 1) {

        if (l == 0) return 0.0;

        // For hydrogen-like atoms
        double E_0 = 13.6 * constants::e * Z * Z;  // Rydberg energy

        // Fine structure constant
        double alpha = constants::alpha;

        // Spin-orbit energy
        double energy = (alpha * alpha * E_0) / (n * n * n) *
                       (j * (j + 1) - l * (l + 1) - 0.75) /
                       (l * (l + 0.5) * (l + 1));

        return energy;
    }

    /**
     * @brief Total angular momentum j
     *
     * j = l ± 1/2
     */
    static std::pair<int, int> possible_j_values(int l) {
        // Returns (2j) for integer storage
        if (l == 0) {
            return {1, 1};  // Only j = 1/2
        }
        return {2 * l - 1, 2 * l + 1};  // j = l - 1/2, l + 1/2
    }

    /**
     * @brief Fine structure splitting
     *
     * ΔE = E(j=l+1/2) - E(j=l-1/2)
     */
    static double fine_structure_splitting(int n, int l, int Z = 1) {
        if (l == 0) return 0.0;

        // Use actual j values (stored as 2j)
        auto [j_minus, j_plus] = possible_j_values(l);

        double E_plus = spin_orbit_energy(n, l, j_plus / 2, Z);
        double E_minus = spin_orbit_energy(n, l, j_minus / 2, Z);

        return std::abs(E_plus - E_minus);
    }

    /**
     * @brief L·S expectation value
     *
     * ⟨L⃗·S⃗⟩ = (ℏ²/2)[j(j+1) - l(l+1) - s(s+1)]
     */
    static double LS_expectation(int l, int j_times_2) {
        double j = j_times_2 / 2.0;
        double s = 0.5;

        return (constants::hbar * constants::hbar / 2.0) *
               (j * (j + 1) - l * (l + 1) - s * (s + 1));
    }

    /**
     * @brief Fine structure constant correction
     *
     * Relativistic correction scale: α² ≈ (1/137)² ≈ 5 × 10⁻⁵
     */
    static double fine_structure_scale() {
        return constants::alpha * constants::alpha;
    }

    /**
     * @brief Thomas precession factor
     *
     * Factor of 1/2 from Thomas precession in moving reference frame
     */
    static double thomas_precession_factor() {
        return 0.5;
    }
};

/**
 * @brief Zeeman Effect Revisited
 *
 * Complete treatment including anomalous Zeeman effect
 */
class ZeemanEffectRevisited {
public:
    /**
     * @brief Normal Zeeman effect (no spin)
     *
     * ΔE = μ_B m_l B
     */
    static double normal_zeeman_shift(int m_l, double B_field) {
        return constants::mu_B * m_l * B_field;
    }

    /**
     * @brief Anomalous Zeeman effect (with spin)
     *
     * ΔE = g_J μ_B m_J B
     */
    static double anomalous_zeeman_shift(
        int m_J,
        int J,
        int L,
        int S,
        double B_field) {

        double g_J = SpinMagneticMoment::lande_g_factor(J, L, S);
        return g_J * constants::mu_B * m_J * B_field;
    }

    /**
     * @brief Paschen-Back effect (strong field)
     *
     * When B is strong: decouple L and S
     * ΔE = μ_B(m_l + g_s m_s)B
     */
    static double paschen_back_shift(
        int m_l,
        int m_s,
        double B_field) {

        return constants::mu_B * (m_l + constants::g_e * m_s) * B_field;
    }

    /**
     * @brief Transition between weak and strong field regimes
     *
     * Compare Zeeman energy to fine structure
     */
    static std::string field_regime(double B_field, int n, int l) {
        double zeeman_energy = constants::mu_B * B_field;
        double fine_structure = SpinOrbitCoupling::fine_structure_splitting(n, l);

        if (zeeman_energy < fine_structure) {
            return "Weak field (anomalous Zeeman)";
        } else {
            return "Strong field (Paschen-Back)";
        }
    }

    /**
     * @brief Zeeman pattern for spectral lines
     *
     * Selection rules: ΔJ = 0, ±1 (not 0→0), Δm_J = 0, ±1
     */
    struct ZeemanPattern {
        enum class Polarization { SIGMA_PLUS, SIGMA_MINUS, PI };

        static bool transition_allowed(int J_i, int m_J_i, int J_f, int m_J_f) {
            int delta_J = J_f - J_i;
            int delta_m = m_J_f - m_J_i;

            // ΔJ = 0, ±1 (but not 0→0)
            if (std::abs(delta_J) > 1) return false;
            if (J_i == 0 && J_f == 0) return false;

            // Δm_J = 0, ±1
            if (std::abs(delta_m) > 1) return false;

            return true;
        }

        static Polarization get_polarization(int delta_m) {
            if (delta_m == +1) return Polarization::SIGMA_PLUS;
            if (delta_m == -1) return Polarization::SIGMA_MINUS;
            return Polarization::PI;
        }
    };

    /**
     * @brief Hyperfine structure (nuclear spin interaction)
     *
     * ΔE_hf = A I⃗·J⃗
     */
    static double hyperfine_splitting(
        int I,
        int J,
        int F,
        double A_constant) {

        // ⟨I⃗·J⃗⟩ = (ℏ²/2)[F(F+1) - I(I+1) - J(J+1)]
        return (A_constant / 2.0) *
               (F * (F + 1) - I * (I + 1) - J * (J + 1));
    }

    /**
     * @brief 21 cm line of hydrogen
     *
     * Hyperfine transition: F=1 → F=0
     */
    static double hydrogen_21cm_frequency() {
        return 1420405751.768;  // Hz
    }
};

/**
 * @brief The Klein-Gordon Equation
 *
 * Relativistic equation for spinless particles
 */
class KleinGordonEquation {
public:
    /**
     * @brief Klein-Gordon equation
     *
     * (□ + m²c²/ℏ²)ψ = 0
     * □ = (1/c²)∂²/∂t² - ∇²
     */
    static bool satisfies_klein_gordon(
        double mass,
        double energy,
        double momentum) {

        // Relativistic energy-momentum relation
        // E² = (pc)² + (mc²)²
        double lhs = energy * energy;
        double rhs = momentum * momentum * constants::c * constants::c +
                     mass * mass * constants::c * constants::c * constants::c * constants::c;

        return std::abs(lhs - rhs) < 1e-10;
    }

    /**
     * @brief Plane wave solutions
     *
     * ψ = e^(i(p·x - Et)/ℏ)
     */
    static Complex plane_wave(
        double x,
        double t,
        double momentum,
        double energy) {

        double phase = (momentum * x - energy * t) / constants::hbar;
        return std::exp(Complex(0, phase));
    }

    /**
     * @brief Positive and negative energy solutions
     *
     * E = ±√((pc)² + (mc²)²)
     */
    static std::pair<double, double> energy_solutions(
        double momentum,
        double mass) {

        double E_squared = momentum * momentum * constants::c * constants::c +
                          mass * mass * constants::c * constants::c * constants::c * constants::c;
        double E = std::sqrt(E_squared);

        return {E, -E};  // Positive and negative energy
    }

    /**
     * @brief Klein paradox
     *
     * Transmission through potential barrier can exceed 1
     */
    static double klein_paradox_transmission(
        double energy,
        double barrier_height,
        double mass) {

        if (barrier_height > energy + 2.0 * mass * constants::c * constants::c) {
            // Pair production regime
            return 1.0;  // Perfect transmission
        }

        // Normal tunneling regime
        return 0.0;  // Simplified
    }

    /**
     * @brief Probability density issue
     *
     * Klein-Gordon ρ is not positive definite
     */
    static double probability_density(
        Complex psi,
        Complex psi_t) {

        // ρ = (iℏ/2mc²)(ψ*∂ψ/∂t - ψ∂ψ*/∂t)
        // Can be negative!
        return std::real(Complex(0, constants::hbar / (2.0 * constants::m_e * constants::c * constants::c)) *
                        (std::conj(psi) * psi_t - psi * std::conj(psi_t)));
    }
};

/**
 * @brief The Dirac Equation
 *
 * Relativistic equation for spin-1/2 particles
 */
class DiracEquation {
public:
    /**
     * @brief Dirac equation
     *
     * (iℏ∂/∂t)ψ = (cα⃗·p⃗ + βmc²)ψ
     */

    /**
     * @brief Dirac matrices (4×4)
     *
     * α_i and β matrices
     */
    struct DiracMatrices {
        // In standard representation (Dirac-Pauli)
        // α_i = [[0, σ_i], [σ_i, 0]]
        // β = [[I, 0], [0, -I]]

        static std::vector<std::vector<Complex>> beta() {
            return {
                {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(0, 0), Complex(-1, 0), Complex(0, 0)},
                {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(-1, 0)}
            };
        }

        // Gamma matrices: γ⁰ = β, γⁱ = βα_i
        static std::vector<std::vector<Complex>> gamma0() {
            return beta();
        }
    };

    /**
     * @brief Dirac spinor (4-component)
     *
     * ψ = (ψ_1, ψ_2, ψ_3, ψ_4)ᵀ
     */
    using DiracSpinor = std::array<Complex, 4>;

    /**
     * @brief Free particle solutions (plane waves)
     *
     * u(p) e^(-iEt/ℏ + ip·x/ℏ) for positive energy
     * v(p) e^(+iEt/ℏ + ip·x/ℏ) for negative energy
     */
    static DiracSpinor positive_energy_spinor(
        double momentum,
        int spin_up) {

        double E = std::sqrt(momentum * momentum * constants::c * constants::c +
                            constants::m_e * constants::m_e * constants::c * constants::c * constants::c * constants::c);

        double N = std::sqrt((E + constants::m_e * constants::c * constants::c) / (2.0 * E));

        if (spin_up) {
            return {
                Complex(N, 0),
                Complex(0, 0),
                Complex(N * momentum * constants::c / (E + constants::m_e * constants::c * constants::c), 0),
                Complex(0, 0)
            };
        } else {
            return {
                Complex(0, 0),
                Complex(N, 0),
                Complex(0, 0),
                Complex(N * momentum * constants::c / (E + constants::m_e * constants::c * constants::c), 0)
            };
        }
    }

    /**
     * @brief Positive definite probability density
     *
     * ρ = ψ†ψ (always positive!)
     */
    static double probability_density(const DiracSpinor& psi) {
        double rho = 0.0;
        for (const auto& component : psi) {
            rho += std::norm(component);
        }
        return rho;
    }

    /**
     * @brief Current density
     *
     * j⃗ = cψ†α⃗ψ
     */
    static double current_density_x(const DiracSpinor& psi) {
        // Simplified: j_x = c Re(ψ₁*ψ₃ + ψ₂*ψ₄)
        return constants::c * std::real(std::conj(psi[0]) * psi[2] +
                                        std::conj(psi[1]) * psi[3]);
    }

    /**
     * @brief Continuity equation
     *
     * ∂ρ/∂t + ∇·j⃗ = 0
     */
    static bool satisfies_continuity() {
        // Dirac equation automatically satisfies continuity
        return true;
    }

    /**
     * @brief Non-relativistic limit
     *
     * Recovers Pauli equation with spin
     */
    static std::string nonrelativistic_limit() {
        return "Pauli equation: (p²/2m - eΦ)ψ + (eℏ/2m)σ⃗·B⃗ψ";
    }
};

/**
 * @brief Spin and the Dirac Particle
 *
 * Intrinsic spin from Dirac equation
 */
class DiracParticleSpin {
public:
    /**
     * @brief Spin operator in Dirac theory
     *
     * S⃗ = (ℏ/2)Σ⃗ where Σ_i = [[σ_i, 0], [0, σ_i]]
     */
    static double spin_magnitude() {
        return constants::hbar * std::sqrt(0.5 * 1.5);  // s=1/2
    }

    /**
     * @brief Helicity operator
     *
     * h = Σ⃗·p̂ (spin projection along momentum)
     */
    static int helicity_eigenvalue(bool right_handed) {
        return right_handed ? +1 : -1;
    }

    /**
     * @brief Gyromagnetic ratio
     *
     * Dirac equation predicts g = 2 exactly
     */
    static double dirac_g_factor() {
        return 2.0;
    }

    /**
     * @brief QED correction to g-factor
     *
     * g = 2(1 + α/2π + ...) ≈ 2.00232
     */
    static double qed_g_factor() {
        return 2.0 * (1.0 + constants::alpha / (2.0 * M_PI));
    }

    /**
     * @brief Magnetic moment
     *
     * μ⃗ = -(e/m)S⃗ (from Dirac equation)
     */
    static double magnetic_moment() {
        return constants::e * constants::hbar / (2.0 * constants::m_e);
    }

    /**
     * @brief Zitterbewegung
     *
     * Rapid oscillation with frequency 2mc²/ℏ
     */
    static double zitterbewegung_frequency() {
        return 2.0 * constants::m_e * constants::c * constants::c / constants::hbar;
    }

    /**
     * @brief Zitterbewegung amplitude
     *
     * λ_C = ℏ/(mc) (Compton wavelength)
     */
    static double zitterbewegung_amplitude() {
        return constants::hbar / (constants::m_e * constants::c);
    }

    /**
     * @brief Spin precession in electromagnetic field
     *
     * Thomas-Bargmann-Michel-Telegdi equation
     */
    static double tbmt_precession_frequency(
        double E_field,
        double B_field,
        double velocity) {

        double gamma = 1.0 / std::sqrt(1.0 - velocity * velocity / (constants::c * constants::c));
        double omega_s = -(constants::e / constants::m_e) *
                        (B_field + (gamma - 1) * E_field / constants::c);
        return omega_s;
    }
};

/**
 * @brief Spin-Orbit Coupling in the Dirac Hamiltonian
 *
 * Automatic appearance of spin-orbit coupling
 */
class DiracSpinOrbitCoupling {
public:
    /**
     * @brief Dirac Hamiltonian for central potential
     *
     * Automatic spin-orbit term from Dirac equation
     */
    static double dirac_so_coupling_strength(double r, int Z) {
        // V(r) = -Ze²/(4πε₀r)
        // Leads to L·S coupling with correct Thomas factor
        double V = -Z * constants::e * constants::e * constants::k_e / r;
        double dV_dr = -V / r;

        return dV_dr / (2.0 * constants::m_e * constants::m_e * constants::c * constants::c * r);
    }

    /**
     * @brief Fine structure from Dirac equation
     *
     * E_nj = mc²[1 + (αZ)²/(n - j - 1/2 + √((j+1/2)² - (αZ)²))²]^(-1/2) - mc²
     */
    static double dirac_energy(int n, int j_times_2, int Z) {
        double j = j_times_2 / 2.0;
        double alpha_Z = constants::alpha * Z;

        double denominator_term = n - j - 0.5 +
                                 std::sqrt((j + 0.5) * (j + 0.5) - alpha_Z * alpha_Z);

        double bracket = 1.0 + (alpha_Z * alpha_Z) / (denominator_term * denominator_term);

        return constants::m_e * constants::c * constants::c *
               (1.0 / std::sqrt(bracket) - 1.0);
    }

    /**
     * @brief Darwin term
     *
     * Additional relativistic correction for s-states
     */
    static double darwin_term(int l, int Z) {
        if (l != 0) return 0.0;

        // (πℏ²/2m²c²)|ψ(0)|²(Ze²/4πε₀)
        // Represents Zitterbewegung averaging
        return (M_PI * constants::hbar * constants::hbar * Z * constants::e * constants::e) /
               (2.0 * constants::m_e * constants::m_e * constants::c * constants::c);
    }

    /**
     * @brief Kinetic energy relativistic correction
     *
     * -(p⁴/8m³c²)
     */
    static double kinetic_correction(double momentum) {
        return -std::pow(momentum, 4) / (8.0 * std::pow(constants::m_e, 3) * constants::c * constants::c);
    }

    /**
     * @brief Total fine structure
     *
     * Combination of spin-orbit + Darwin + kinetic
     */
    static double total_fine_structure(int n, int l, int j_times_2, int Z) {
        // Use Dirac energy which includes all corrections automatically
        return dirac_energy(n, j_times_2, Z);
    }
};

/**
 * @brief The Dirac Hydrogen Atom
 *
 * Exact solution including fine structure
 */
class DiracHydrogenAtom {
public:
    /**
     * @brief Quantum number j
     *
     * j = l ± 1/2 (total angular momentum)
     */
    static int quantum_number_j(int l, bool j_plus) {
        // Returns 2j for integer storage
        return j_plus ? (2 * l + 1) : (2 * l - 1);
    }

    /**
     * @brief Dirac energy levels
     *
     * E_nj = mc²[(1 + (αZ)²/(n - j - 1/2 + κ)²)^(-1/2) - 1]
     */
    static double energy_level(int n, int j_times_2, int Z = 1) {
        return DiracSpinOrbitCoupling::dirac_energy(n, j_times_2, Z);
    }

    /**
     * @brief Fine structure constant role
     *
     * α ≈ 1/137 determines size of fine structure
     */
    static double fine_structure_constant() {
        return constants::alpha;
    }

    /**
     * @brief Lamb shift
     *
     * QED correction beyond Dirac equation
     * Mainly affects s-states
     */
    static double lamb_shift(int n, int l) {
        if (l != 0) return 0.0;  // Mainly s-states

        // Approximate formula for hydrogen
        // Δν ≈ 1057.8 MHz for 2S_1/2 - 2P_1/2
        if (n == 2) {
            return 1057.8e6 * constants::h;  // Convert to energy
        }

        return 0.0;
    }

    /**
     * @brief Degeneracy with same j
     *
     * States with same n and j but different l are degenerate in Dirac theory
     * (Lamb shift breaks this)
     */
    static bool is_degenerate(int n1, int l1, int j1, int n2, int l2, int j2) {
        return (n1 == n2) && (j1 == j2);
    }

    /**
     * @brief Spectrum of hydrogen
     *
     * Including fine structure
     */
    struct HydrogenSpectrum {
        int n;
        int l;
        int j_times_2;
        double energy;
        std::string label;

        static HydrogenSpectrum ground_state() {
            return {1, 0, 1, energy_level(1, 1), "1S_1/2"};
        }

        static std::vector<HydrogenSpectrum> n2_states() {
            return {
                {2, 0, 1, energy_level(2, 1), "2S_1/2"},
                {2, 1, 1, energy_level(2, 1), "2P_1/2"},
                {2, 1, 3, energy_level(2, 3), "2P_3/2"}
            };
        }
    };

    /**
     * @brief Radial wave functions
     *
     * Modified from hydrogen to include relativity
     */
    static double radial_wavefunction(double r, int n, int kappa) {
        // Simplified - full solution requires confluent hypergeometric functions
        double a_0 = constants::a_0;
        return std::exp(-r / (n * a_0));
    }
};

/**
 * @brief The Dirac Particle in a Magnetic Field
 *
 * Landau levels and gyromagnetic ratio
 */
class DiracParticleInMagneticField {
public:
    /**
     * @brief Minimal coupling
     *
     * p⃗ → p⃗ - eA⃗
     */
    static double canonical_momentum(double kinetic_momentum, double vector_potential) {
        return kinetic_momentum + constants::e * vector_potential;
    }

    /**
     * @brief Landau levels for Dirac particle
     *
     * E_n = ±√((mc²)² + 2n|e|ℏcB)
     */
    static std::pair<double, double> landau_level_energy(int n, double B_field) {
        double mc2 = constants::m_e * constants::c * constants::c;
        double orbital_term = 2.0 * n * constants::e * constants::hbar * constants::c * B_field;

        double E_squared = mc2 * mc2 + orbital_term;
        double E = std::sqrt(E_squared);

        return {E, -E};
    }

    /**
     * @brief Cyclotron frequency
     *
     * ω_c = eB/m
     */
    static double cyclotron_frequency(double B_field) {
        return constants::e * B_field / constants::m_e;
    }

    /**
     * @brief Magnetic length
     *
     * l_B = √(ℏ/eB)
     */
    static double magnetic_length(double B_field) {
        return std::sqrt(constants::hbar / (constants::e * B_field));
    }

    /**
     * @brief Degeneracy of Landau level
     *
     * g = eBA/(2πℏ) (per unit area A)
     */
    static double landau_degeneracy_per_area(double B_field) {
        return constants::e * B_field / (2.0 * M_PI * constants::hbar);
    }

    /**
     * @brief Pauli term from Dirac equation
     *
     * -(eℏ/2m)σ⃗·B⃗ appears automatically
     */
    static double pauli_interaction_energy(int spin_sign, double B_field) {
        return -(constants::e * constants::hbar / (2.0 * constants::m_e)) *
               spin_sign * B_field;
    }

    /**
     * @brief Anomalous magnetic moment
     *
     * QED correction: g = 2(1 + α/2π) ≈ 2.00232
     */
    static double anomalous_magnetic_moment() {
        return (constants::alpha / (2.0 * M_PI)) * constants::mu_B;
    }

    /**
     * @brief Klein paradox in magnetic field
     *
     * Pair production in strong magnetic fields
     */
    static double critical_magnetic_field() {
        // B_c = m²c³/(eℏ) ≈ 4.4 × 10¹³ G
        return (constants::m_e * constants::m_e * constants::c * constants::c * constants::c) /
               (constants::e * constants::hbar);
    }

    /**
     * @brief Synchrotron radiation
     *
     * Energy loss for relativistic particle in B field
     */
    static double synchrotron_power(double gamma, double B_field) {
        // P = (2e⁴/3m²c³)γ²B²
        double e4 = std::pow(constants::e, 4);
        return (2.0 * e4 / (3.0 * constants::m_e * constants::m_e * constants::c * constants::c * constants::c)) *
               gamma * gamma * B_field * B_field;
    }
};

} // namespace relativistic_quantum
} // namespace physics

#endif // PHYSICS_RELATIVISTIC_QUANTUM_MECHANICS_HPP
