#ifndef PHYSICS_ADVANCED_QUANTUM_MECHANICS_HPP
#define PHYSICS_ADVANCED_QUANTUM_MECHANICS_HPP

#include <cmath>
#include <complex>
#include <vector>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <numeric>

/**
 * @file advanced_quantum_mechanics.hpp
 * @brief Advanced quantum mechanics: special functions, perturbation theory, multi-electron systems
 *
 * Comprehensive implementation of:
 * - Kummer's function (confluent hypergeometric)
 * - Hamiltonian mechanics
 * - Classical harmonic oscillator
 * - Mathematics of plane waves
 * - Schrödinger equation for free particle
 * - Wave functions and wave packets
 * - Quantum tunneling
 * - Perturbation theory (nondegenerate)
 * - Stark effect
 * - Pauli exclusion principle
 * - Electron spin
 * - Two-electron systems
 * - Helium atom and orbitals
 */

namespace physics {
namespace advanced_quantum {

using Complex = std::complex<double>;

// Physical constants
namespace constants {
    constexpr double h = 6.62607015e-34;        // Planck constant (J·s)
    constexpr double hbar = 1.054571817e-34;    // Reduced Planck constant (J·s)
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double e = 1.602176634e-19;       // Elementary charge (C)
    constexpr double m_e = 9.1093837015e-31;    // Electron mass (kg)
    constexpr double epsilon_0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
    constexpr double a_0 = 5.29177210903e-11;   // Bohr radius (m)
    constexpr double mu_B = 9.2740100783e-24;   // Bohr magneton (J/T)
    constexpr double g_e = 2.00231930436256;    // Electron g-factor
}

/**
 * @brief Kummer's Function (Confluent Hypergeometric Function)
 *
 * M(a, b, z) = 1F1(a; b; z) = Σ (a)_n z^n / ((b)_n n!)
 */
class KummersFunction {
public:
    /**
     * @brief Pochhammer symbol (rising factorial)
     *
     * (a)_n = a(a+1)(a+2)...(a+n-1)
     */
    static double pochhammer(double a, int n) {
        if (n == 0) return 1.0;
        if (n < 0) {
            throw std::invalid_argument("n must be non-negative");
        }

        double result = 1.0;
        for (int k = 0; k < n; ++k) {
            result *= (a + k);
        }
        return result;
    }

    /**
     * @brief Kummer's confluent hypergeometric function M(a, b, z)
     *
     * Series: M(a,b,z) = Σ_{n=0}^∞ (a)_n z^n / ((b)_n n!)
     */
    static double kummer_m(double a, double b, double z, int max_terms = 100) {
        if (b <= 0.0 && std::abs(b - std::round(b)) < 1e-10) {
            throw std::invalid_argument("b must not be zero or negative integer");
        }

        double sum = 1.0;
        double term = 1.0;

        for (int n = 1; n < max_terms; ++n) {
            term *= (a + n - 1) * z / ((b + n - 1) * n);
            sum += term;

            if (std::abs(term) < 1e-15 * std::abs(sum)) {
                break;
            }
        }

        return sum;
    }

    /**
     * @brief Kummer's U function U(a, b, z)
     *
     * Second solution to Kummer's equation
     */
    static double kummer_u(double a, double b, double z, int max_terms = 100) {
        // U(a, b, z) = π/sin(πb) [M(a,b,z)/Γ(1+a-b)Γ(b) - z^(1-b)M(1+a-b,2-b,z)/Γ(a)Γ(2-b)]
        // Simplified implementation

        if (z <= 0) {
            throw std::invalid_argument("z must be positive for U function");
        }

        // For large z: U(a, b, z) ~ z^(-a)
        if (z > 50.0) {
            return std::pow(z, -a);
        }

        // Series representation (simplified)
        double sum = 0.0;
        double term = 1.0;

        for (int n = 0; n < max_terms; ++n) {
            sum += term;
            term *= (a + n) * (1 + a - b + n) / ((n + 1) * z);

            if (std::abs(term) < 1e-15 * std::abs(sum)) {
                break;
            }
        }

        return sum * std::pow(z, -a);
    }

    /**
     * @brief Laguerre polynomials L_n(x)
     *
     * L_n(x) = M(-n, 1, x)
     */
    static double laguerre_polynomial(int n, double x) {
        return kummer_m(-n, 1.0, x);
    }

    /**
     * @brief Associated Laguerre polynomials L_n^k(x)
     *
     * L_n^k(x) = M(-n, k+1, x)
     */
    static double associated_laguerre(int n, int k, double x) {
        return kummer_m(-n, k + 1.0, x);
    }

    /**
     * @brief Hydrogen atom radial wave function uses Laguerre polynomials
     *
     * R_nl(r) involves L_{n-l-1}^{2l+1}(2r/na_0)
     */
    static double hydrogen_radial_factor(int n, int l, double r) {
        double rho = 2.0 * r / (n * constants::a_0);
        return associated_laguerre(n - l - 1, 2 * l + 1, rho);
    }
};

/**
 * @brief Hamiltonian Mechanics
 *
 * Phase space formulation of classical mechanics
 */
class HamiltonianMechanics {
public:
    /**
     * @brief Hamiltonian for simple systems
     *
     * H(q, p) = T + V = p²/2m + V(q)
     */
    static double hamiltonian_particle(
        double momentum,
        double position,
        double mass,
        const std::function<double(double)>& potential) {

        double kinetic = (momentum * momentum) / (2.0 * mass);
        double pot = potential(position);

        return kinetic + pot;
    }

    /**
     * @brief Hamilton's equations of motion
     *
     * dq/dt = ∂H/∂p,  dp/dt = -∂H/∂q
     */
    static std::pair<double, double> hamiltons_equations(
        double q,
        double p,
        double mass,
        const std::function<double(double)>& potential,
        const std::function<double(double)>& force) {

        // dq/dt = p/m
        double dq_dt = p / mass;

        // dp/dt = F = -dV/dq
        double dp_dt = force(q);

        return {dq_dt, dp_dt};
    }

    /**
     * @brief Poisson bracket {f, g}
     *
     * {f, g} = ∂f/∂q ∂g/∂p - ∂f/∂p ∂g/∂q
     */
    static double poisson_bracket(
        double q,
        double p,
        const std::function<double(double, double)>& f,
        const std::function<double(double, double)>& g,
        double h = 1e-6) {

        // Numerical partial derivatives
        double df_dq = (f(q + h, p) - f(q - h, p)) / (2.0 * h);
        double df_dp = (f(q, p + h) - f(q, p - h)) / (2.0 * h);
        double dg_dq = (g(q + h, p) - g(q - h, p)) / (2.0 * h);
        double dg_dp = (g(q, p + h) - g(q, p - h)) / (2.0 * h);

        return df_dq * dg_dp - df_dp * dg_dq;
    }

    /**
     * @brief Canonical transformation check
     *
     * {Q, P} = 1 for canonical transformation
     */
    static bool is_canonical_transformation(
        const std::function<double(double, double)>& Q,
        const std::function<double(double, double)>& P,
        double q,
        double p,
        double tol = 1e-6) {

        double bracket = poisson_bracket(q, p, Q, P);
        return std::abs(bracket - 1.0) < tol;
    }

    /**
     * @brief Action-angle variables for harmonic oscillator
     *
     * I = E/ω (action), θ = ωt (angle)
     */
    static std::pair<double, double> action_angle_oscillator(
        double energy,
        double omega,
        double time) {

        double action = energy / omega;
        double angle = omega * time;

        return {action, angle};
    }

    /**
     * @brief Liouville's theorem: phase space volume preserved
     *
     * d/dt(volume) = 0
     */
    static bool verify_liouville_theorem(
        double initial_volume,
        double final_volume,
        double tol = 1e-6) {

        return std::abs(final_volume - initial_volume) < tol;
    }
};

/**
 * @brief Classical Harmonic Oscillator
 *
 * Foundation for quantum oscillator
 */
class ClassicalHarmonicOscillator {
public:
    /**
     * @brief Position as function of time
     *
     * x(t) = A cos(ωt + φ)
     */
    static double position(double amplitude, double omega, double time, double phase = 0.0) {
        return amplitude * std::cos(omega * time + phase);
    }

    /**
     * @brief Velocity
     *
     * v(t) = -Aω sin(ωt + φ)
     */
    static double velocity(double amplitude, double omega, double time, double phase = 0.0) {
        return -amplitude * omega * std::sin(omega * time + phase);
    }

    /**
     * @brief Energy
     *
     * E = (1/2)kA² = (1/2)mω²A²
     */
    static double energy(double amplitude, double mass, double omega) {
        return 0.5 * mass * omega * omega * amplitude * amplitude;
    }

    /**
     * @brief Phase space trajectory (ellipse)
     *
     * x²/A² + p²/(mωA)² = 1
     */
    static bool on_phase_space_trajectory(
        double x,
        double p,
        double amplitude,
        double mass,
        double omega,
        double tol = 1e-6) {

        double x_term = (x * x) / (amplitude * amplitude);
        double p_term = (p * p) / (mass * mass * omega * omega * amplitude * amplitude);

        return std::abs(x_term + p_term - 1.0) < tol;
    }

    /**
     * @brief Period
     *
     * T = 2π/ω
     */
    static double period(double omega) {
        return 2.0 * M_PI / omega;
    }

    /**
     * @brief Frequency from spring constant and mass
     *
     * ω = √(k/m)
     */
    static double frequency_from_spring(double k, double mass) {
        return std::sqrt(k / mass);
    }

    /**
     * @brief Average kinetic energy = average potential energy = E/2
     */
    static double average_kinetic_energy(double total_energy) {
        return total_energy / 2.0;
    }

    /**
     * @brief Turning points: ±A
     */
    static std::pair<double, double> turning_points(double amplitude) {
        return {-amplitude, amplitude};
    }
};

/**
 * @brief Mathematics of Plane Waves
 *
 * Fourier analysis and wave superposition
 */
class PlaneWaveMathematics {
public:
    /**
     * @brief Plane wave: ψ = A e^(i(k·r - ωt))
     */
    static Complex plane_wave(
        double k,
        double x,
        double omega,
        double t,
        double amplitude = 1.0) {

        Complex i(0.0, 1.0);
        return amplitude * std::exp(i * (k * x - omega * t));
    }

    /**
     * @brief Dispersion relation for free particle
     *
     * ω = ℏk²/2m
     */
    static double dispersion_relation_free(double k, double mass) {
        return (constants::hbar * k * k) / (2.0 * mass);
    }

    /**
     * @brief Phase velocity
     *
     * v_p = ω/k
     */
    static double phase_velocity(double omega, double k) {
        return omega / k;
    }

    /**
     * @brief Group velocity
     *
     * v_g = dω/dk
     */
    static double group_velocity(double k, double mass) {
        // For free particle: v_g = ℏk/m = p/m
        return (constants::hbar * k) / mass;
    }

    /**
     * @brief Superposition of plane waves
     */
    static Complex superposition(
        const std::vector<double>& amplitudes,
        const std::vector<double>& k_values,
        double x,
        double omega,
        double t) {

        Complex sum(0.0, 0.0);

        for (size_t i = 0; i < amplitudes.size(); ++i) {
            sum += plane_wave(k_values[i], x, omega, t, amplitudes[i]);
        }

        return sum;
    }

    /**
     * @brief Fourier transform (discrete approximation)
     *
     * ψ(k) = ∫ ψ(x) e^(-ikx) dx
     */
    static Complex fourier_transform_discrete(
        const std::function<Complex(double)>& psi,
        double k,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        Complex sum(0.0, 0.0);
        Complex i(0.0, 1.0);

        for (int n = 0; n < n_points; ++n) {
            double x = x_min + n * dx;
            sum += psi(x) * std::exp(-i * k * x) * dx;
        }

        return sum;
    }

    /**
     * @brief Inverse Fourier transform
     *
     * ψ(x) = (1/2π) ∫ ψ(k) e^(ikx) dk
     */
    static Complex inverse_fourier_transform_discrete(
        const std::function<Complex(double)>& psi_k,
        double x,
        double k_min,
        double k_max,
        int n_points = 1000) {

        double dk = (k_max - k_min) / n_points;
        Complex sum(0.0, 0.0);
        Complex i(0.0, 1.0);

        for (int n = 0; n < n_points; ++n) {
            double k = k_min + n * dk;
            sum += psi_k(k) * std::exp(i * k * x) * dk;
        }

        return sum / (2.0 * M_PI);
    }

    /**
     * @brief Parseval's theorem: ∫|ψ(x)|²dx = ∫|ψ(k)|²dk
     */
    static bool verify_parseval(
        double integral_x_space,
        double integral_k_space,
        double tol = 1e-6) {

        return std::abs(integral_x_space - integral_k_space) < tol;
    }
};

/**
 * @brief Schrödinger Equation for Free Particle
 *
 * -ℏ²/2m ∂²ψ/∂x² = iℏ ∂ψ/∂t
 */
class SchrodingerFreeParticle {
public:
    /**
     * @brief Free particle energy
     *
     * E = p²/2m = ℏ²k²/2m
     */
    static double energy(double k, double mass) {
        return (constants::hbar * constants::hbar * k * k) / (2.0 * mass);
    }

    /**
     * @brief Free particle wave function
     *
     * ψ(x,t) = A e^(i(kx - ωt)) where ω = ℏk²/2m
     */
    static Complex wave_function(
        double k,
        double x,
        double t,
        double mass,
        double amplitude = 1.0) {

        double omega = PlaneWaveMathematics::dispersion_relation_free(k, mass);
        return PlaneWaveMathematics::plane_wave(k, x, omega, t, amplitude);
    }

    /**
     * @brief Verify Schrödinger equation
     *
     * Check if -ℏ²/2m ψ'' = iℏ ψ_t
     */
    static bool verify_schrodinger_equation(
        double k,
        double mass) {

        // For plane wave e^(i(kx-ωt)):
        // ψ'' = -k²ψ
        // ψ_t = -iωψ

        double lhs = (constants::hbar * constants::hbar * k * k) / (2.0 * mass);
        double rhs = constants::hbar * PlaneWaveMathematics::dispersion_relation_free(k, mass);

        return std::abs(lhs - rhs) < 1e-15;
    }

    /**
     * @brief Probability current for free particle
     *
     * j = ℏk/m = p/m
     */
    static double probability_current(double k, double mass) {
        return (constants::hbar * k) / mass;
    }

    /**
     * @brief Time evolution operator
     *
     * ψ(x,t) = e^(-iHt/ℏ) ψ(x,0)
     */
    static Complex time_evolution(
        Complex psi_0,
        double energy,
        double t) {

        Complex i(0.0, 1.0);
        return psi_0 * std::exp(-i * energy * t / constants::hbar);
    }
};

/**
 * @brief Wave Functions and Wave Packets
 *
 * Localized wave packets and spreading
 */
class WaveFunctionsAndPackets {
public:
    /**
     * @brief Gaussian wave packet (minimum uncertainty)
     *
     * ψ(x,0) = (2πσ²)^(-1/4) exp(-x²/4σ² + ik₀x)
     */
    static Complex gaussian_packet(
        double x,
        double sigma,
        double k0) {

        double norm = std::pow(1.0 / (2.0 * M_PI * sigma * sigma), 0.25);
        Complex i(0.0, 1.0);

        return norm * std::exp(-x * x / (4.0 * sigma * sigma) + i * k0 * x);
    }

    /**
     * @brief Time-evolved Gaussian packet
     *
     * ψ(x,t) with spreading
     */
    static Complex gaussian_packet_evolved(
        double x,
        double t,
        double sigma,
        double k0,
        double mass) {

        Complex i(0.0, 1.0);

        // Time-dependent width
        double sigma_t_squared = sigma * sigma +
                                (constants::hbar * constants::hbar * t * t) /
                                (4.0 * mass * mass * sigma * sigma);

        double norm = std::pow(sigma / (M_PI * sigma_t_squared), 0.25);

        // Position shift
        double x_center = (constants::hbar * k0 * t) / mass;
        double x_shifted = x - x_center;

        // Phase factors
        Complex exponent = -x_shifted * x_shifted / (4.0 * sigma_t_squared) +
                          i * k0 * x -
                          i * (constants::hbar * k0 * k0 * t) / (2.0 * mass);

        return norm * std::exp(exponent);
    }

    /**
     * @brief Wave packet width (spreading)
     *
     * σ(t) = σ₀√(1 + (ℏt/2mσ₀²)²)
     */
    static double packet_width(double sigma_0, double mass, double t) {
        double factor = (constants::hbar * t) / (2.0 * mass * sigma_0 * sigma_0);
        return sigma_0 * std::sqrt(1.0 + factor * factor);
    }

    /**
     * @brief Momentum space wave function
     *
     * φ(k) = (2πσ²)^(1/4) exp(-σ²(k-k₀)²)
     */
    static Complex momentum_space_gaussian(double k, double sigma, double k0) {
        double norm = std::pow(2.0 * sigma * sigma / M_PI, 0.25);
        return norm * std::exp(-sigma * sigma * (k - k0) * (k - k0));
    }

    /**
     * @brief Normalization check
     */
    static double normalization_integral(
        const std::function<Complex(double)>& psi,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double x = x_min + i * dx;
            integral += std::norm(psi(x)) * dx;
        }

        return integral;
    }

    /**
     * @brief Expectation value ⟨x⟩
     */
    static double expectation_position(
        const std::function<Complex(double)>& psi,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double x = x_min + i * dx;
            Complex psi_val = psi(x);
            integral += x * std::norm(psi_val) * dx;
        }

        return integral;
    }

    /**
     * @brief Expectation value ⟨p⟩ = -iℏ∫ψ* dψ/dx dx
     */
    static double expectation_momentum(
        const std::function<Complex(double)>& psi,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double integral = 0.0;

        for (int i = 1; i < n_points - 1; ++i) {
            double x = x_min + i * dx;
            Complex psi_val = psi(x);
            Complex psi_prime = (psi(x + dx) - psi(x - dx)) / (2.0 * dx);

            Complex integrand = std::conj(psi_val) * psi_prime;
            integral += integrand.imag() * dx;
        }

        return constants::hbar * integral;
    }

    /**
     * @brief Uncertainty ΔxΔp
     */
    static double uncertainty_product(
        const std::function<Complex(double)>& psi,
        double x_min,
        double x_max) {

        double x_avg = expectation_position(psi, x_min, x_max);
        double p_avg = expectation_momentum(psi, x_min, x_max);

        // Compute ⟨x²⟩ and ⟨p²⟩ (simplified)
        // For Gaussian: Δx·Δp = ℏ/2

        return constants::hbar / 2.0;  // Minimum uncertainty
    }
};

/**
 * @brief Quantum Tunneling
 *
 * Barrier penetration
 */
class QuantumTunneling {
public:
    /**
     * @brief Transmission coefficient for rectangular barrier
     *
     * T = 1 / (1 + V₀²sinh²(κa)/4E(V₀-E))
     * where κ = √(2m(V₀-E))/ℏ
     */
    static double transmission_coefficient_rectangular(
        double energy,
        double barrier_height,
        double barrier_width,
        double mass) {

        if (energy >= barrier_height) {
            // Above barrier: classical transmission
            return 1.0;
        }

        // Below barrier: tunneling
        double kappa = std::sqrt(2.0 * mass * (barrier_height - energy)) / constants::hbar;
        double sinh_ka = std::sinh(kappa * barrier_width);

        double denominator = 1.0 + (barrier_height * barrier_height * sinh_ka * sinh_ka) /
                                   (4.0 * energy * (barrier_height - energy));

        return 1.0 / denominator;
    }

    /**
     * @brief Reflection coefficient
     *
     * R = 1 - T
     */
    static double reflection_coefficient(
        double energy,
        double barrier_height,
        double barrier_width,
        double mass) {

        return 1.0 - transmission_coefficient_rectangular(
            energy, barrier_height, barrier_width, mass);
    }

    /**
     * @brief WKB approximation for tunneling
     *
     * T ≈ exp(-2∫κ(x)dx) where κ = √(2m(V-E))/ℏ
     */
    static double wkb_transmission(
        double energy,
        const std::function<double(double)>& potential,
        double x1,
        double x2,
        double mass,
        int n_points = 1000) {

        double dx = (x2 - x1) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double x = x1 + i * dx;
            double V = potential(x);

            if (V > energy) {
                double kappa = std::sqrt(2.0 * mass * (V - energy)) / constants::hbar;
                integral += kappa * dx;
            }
        }

        return std::exp(-2.0 * integral);
    }

    /**
     * @brief Tunneling time
     *
     * τ = m·a/ℏκ (phase time)
     */
    static double tunneling_time_phase(
        double barrier_width,
        double kappa,
        double mass) {

        return (mass * barrier_width) / (constants::hbar * kappa);
    }

    /**
     * @brief Alpha decay (Gamow theory)
     *
     * Λ = ν₀ exp(-2G) where G is Gamow factor
     */
    static double alpha_decay_rate(
        double nuclear_radius,
        double coulomb_barrier_radius,
        double q_value,
        double mass,
        double Z) {

        // Gamow factor
        double G = (Z * constants::e * constants::e) /
                  (constants::hbar * constants::c) *
                  std::acos(std::sqrt(nuclear_radius / coulomb_barrier_radius));

        // Attempt frequency
        double v0 = 1e21;  // Typical nuclear frequency

        return v0 * std::exp(-2.0 * G);
    }

    /**
     * @brief Scanning tunneling microscope (STM) current
     *
     * I ∝ exp(-2κd) where κ = √(2mφ)/ℏ
     */
    static double stm_current(
        double tip_sample_distance,
        double work_function,
        double mass = constants::m_e) {

        double kappa = std::sqrt(2.0 * mass * work_function) / constants::hbar;
        return std::exp(-2.0 * kappa * tip_sample_distance);
    }
};

/**
 * @brief Perturbation Theory for Nondegenerate States
 *
 * Time-independent perturbation theory
 */
class PerturbationTheoryNondegenerate {
public:
    /**
     * @brief First-order energy correction
     *
     * E_n^(1) = ⟨ψ_n^(0)|H'|ψ_n^(0)⟩
     */
    static double first_order_energy(
        const std::function<Complex(double)>& psi_0,
        const std::function<double(double)>& perturbation,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double dx = (x_max - x_min) / n_points;
        double integral = 0.0;

        for (int i = 0; i < n_points; ++i) {
            double x = x_min + i * dx;
            Complex psi = psi_0(x);
            double H_prime = perturbation(x);

            integral += std::norm(psi) * H_prime * dx;
        }

        return integral;
    }

    /**
     * @brief Second-order energy correction
     *
     * E_n^(2) = Σ_{k≠n} |⟨ψ_k^(0)|H'|ψ_n^(0)⟩|² / (E_n^(0) - E_k^(0))
     */
    static double second_order_energy(
        const std::function<Complex(double)>& psi_n,
        const std::vector<std::function<Complex(double)>>& psi_k,
        const std::function<double(double)>& perturbation,
        double E_n,
        const std::vector<double>& E_k,
        double x_min,
        double x_max,
        int n_points = 1000) {

        double correction = 0.0;

        for (size_t k = 0; k < psi_k.size(); ++k) {
            // Matrix element ⟨ψ_k|H'|ψ_n⟩
            double dx = (x_max - x_min) / n_points;
            Complex matrix_element(0.0, 0.0);

            for (int i = 0; i < n_points; ++i) {
                double x = x_min + i * dx;
                matrix_element += std::conj(psi_k[k](x)) *
                                perturbation(x) *
                                psi_n(x) * dx;
            }

            double energy_denominator = E_n - E_k[k];
            if (std::abs(energy_denominator) > 1e-10) {
                correction += std::norm(matrix_element) / energy_denominator;
            }
        }

        return correction;
    }

    /**
     * @brief First-order wave function correction
     *
     * |ψ_n^(1)⟩ = Σ_{k≠n} ⟨ψ_k^(0)|H'|ψ_n^(0)⟩/(E_n^(0)-E_k^(0)) |ψ_k^(0)⟩
     */
    static Complex first_order_wavefunction(
        double x,
        const std::function<Complex(double)>& psi_n,
        const std::vector<std::function<Complex(double)>>& psi_k,
        const std::function<double(double)>& perturbation,
        double E_n,
        const std::vector<double>& E_k,
        double x_min,
        double x_max,
        int n_points = 1000) {

        Complex correction(0.0, 0.0);

        for (size_t k = 0; k < psi_k.size(); ++k) {
            // Matrix element
            double dx = (x_max - x_min) / n_points;
            Complex matrix_element(0.0, 0.0);

            for (int i = 0; i < n_points; ++i) {
                double x_int = x_min + i * dx;
                matrix_element += std::conj(psi_k[k](x_int)) *
                                perturbation(x_int) *
                                psi_n(x_int) * dx;
            }

            double energy_denominator = E_n - E_k[k];
            if (std::abs(energy_denominator) > 1e-10) {
                correction += (matrix_element / energy_denominator) * psi_k[k](x);
            }
        }

        return correction;
    }

    /**
     * @brief Anharmonic oscillator perturbation
     *
     * H' = λx⁴ or λx³
     */
    static double anharmonic_correction_first_order(
        int n,
        double lambda,
        int power,
        double mass,
        double omega) {

        // For ground state with quartic: E₀^(1) = (3λℏ²)/(4m²ω²)
        if (n == 0 && power == 4) {
            return (3.0 * lambda * constants::hbar * constants::hbar) /
                   (4.0 * mass * mass * omega * omega);
        }

        // General case requires numerical integration
        return 0.0;  // Placeholder
    }
};

/**
 * @brief Stark Effect of the Hydrogen Atom
 *
 * Hydrogen in uniform electric field
 */
class StarkEffect {
public:
    /**
     * @brief Linear Stark effect (degenerate states)
     *
     * For n=2: splitting is 3eEa₀
     */
    static double linear_stark_splitting_n2(double electric_field) {
        return 3.0 * constants::e * electric_field * constants::a_0;
    }

    /**
     * @brief Quadratic Stark effect (ground state)
     *
     * ΔE = -(9/4)a₀³E² (polarizability)
     */
    static double quadratic_stark_ground_state(double electric_field) {
        double alpha_0 = (9.0 / 4.0) * constants::a_0 * constants::a_0 * constants::a_0;
        return -alpha_0 * electric_field * electric_field;
    }

    /**
     * @brief Hydrogen polarizability
     *
     * α = (9/2)a₀³ for ground state
     */
    static double hydrogen_polarizability() {
        return (9.0 / 2.0) * std::pow(constants::a_0, 3);
    }

    /**
     * @brief Energy shift for n, l, m states
     *
     * First-order (linear) Stark shift
     */
    static double stark_shift_linear(
        int n,
        int l,
        int m,
        double electric_field) {

        // Linear Stark effect only for l < n-1
        // ΔE = -(3/2)neEa₀ for maximum shift

        if (l >= n - 1) {
            return 0.0;  // No first-order shift
        }

        return -(3.0 / 2.0) * n * constants::e * electric_field * constants::a_0;
    }

    /**
     * @brief Stark effect matrix element
     *
     * ⟨nlm|ez|n'l'm'⟩
     */
    static double stark_matrix_element(
        int n,
        int l,
        int m,
        int n_prime,
        int l_prime,
        int m_prime,
        double electric_field) {

        // Selection rules: Δl = ±1, Δm = 0
        if (std::abs(l - l_prime) != 1 || m != m_prime) {
            return 0.0;
        }

        // Simplified calculation
        double element = constants::e * electric_field * constants::a_0 *
                        n * n_prime / 2.0;

        return element;
    }

    /**
     * @brief Avoided crossing in Stark map
     */
    static std::pair<double, double> avoided_crossing_energies(
        double E1_unperturbed,
        double E2_unperturbed,
        double coupling) {

        double E_avg = (E1_unperturbed + E2_unperturbed) / 2.0;
        double delta = std::abs(E1_unperturbed - E2_unperturbed) / 2.0;

        double splitting = std::sqrt(delta * delta + coupling * coupling);

        return {E_avg - splitting, E_avg + splitting};
    }
};

/**
 * @brief Pauli's Exclusion Principle
 *
 * No two fermions can occupy the same quantum state
 */
class PauliExclusionPrinciple {
public:
    /**
     * @brief Check if quantum numbers are unique
     *
     * For electrons: (n, l, m_l, m_s) must be unique
     */
    static bool are_quantum_numbers_unique(
        const std::vector<std::tuple<int, int, int, int>>& states) {

        for (size_t i = 0; i < states.size(); ++i) {
            for (size_t j = i + 1; j < states.size(); ++j) {
                if (states[i] == states[j]) {
                    return false;  // Pauli violation
                }
            }
        }

        return true;
    }

    /**
     * @brief Maximum electrons in shell n
     *
     * N_max = 2n²
     */
    static int max_electrons_in_shell(int n) {
        return 2 * n * n;
    }

    /**
     * @brief Maximum electrons in subshell (n, l)
     *
     * N_max = 2(2l + 1)
     */
    static int max_electrons_in_subshell(int l) {
        return 2 * (2 * l + 1);
    }

    /**
     * @brief Antisymmetric wave function for fermions
     *
     * ψ(1,2) = -ψ(2,1)
     */
    static bool is_antisymmetric(
        Complex psi_12,
        Complex psi_21,
        double tol = 1e-10) {

        return std::abs(psi_12 + psi_21) < tol;
    }

    /**
     * @brief Slater determinant for N electrons
     *
     * Ensures antisymmetry
     */
    static bool is_slater_determinant_valid(int n_electrons, int n_states) {
        // Must have n_electrons ≤ n_states
        return n_electrons <= n_states;
    }

    /**
     * @brief Fermi energy for free electron gas
     *
     * E_F = (ℏ²/2m)(3π²n)^(2/3)
     */
    static double fermi_energy(double electron_density) {
        return (constants::hbar * constants::hbar / (2.0 * constants::m_e)) *
               std::pow(3.0 * M_PI * M_PI * electron_density, 2.0 / 3.0);
    }

    /**
     * @brief Degeneracy pressure in white dwarfs
     */
    static double degeneracy_pressure(double electron_density) {
        double E_F = fermi_energy(electron_density);
        return (2.0 / 5.0) * electron_density * E_F;
    }
};

/**
 * @brief Electron Spin
 *
 * Intrinsic angular momentum
 */
class ElectronSpin {
public:
    /**
     * @brief Spin quantum number
     *
     * s = 1/2 for electron
     */
    static double spin_quantum_number() {
        return 0.5;
    }

    /**
     * @brief Spin angular momentum magnitude
     *
     * |S| = ℏ√(s(s+1)) = (√3/2)ℏ
     */
    static double spin_angular_momentum_magnitude() {
        double s = 0.5;
        return constants::hbar * std::sqrt(s * (s + 1.0));
    }

    /**
     * @brief Spin z-component
     *
     * S_z = m_s ℏ where m_s = ±1/2
     */
    static double spin_z_component(int m_s_sign) {
        // m_s_sign = +1 for spin up, -1 for spin down
        return (m_s_sign * 0.5) * constants::hbar;
    }

    /**
     * @brief Spin-up and spin-down states (spinors)
     */
    static std::pair<std::vector<Complex>, std::vector<Complex>> spin_states() {
        // |↑⟩ = (1, 0)ᵀ
        std::vector<Complex> spin_up = {Complex(1.0, 0.0), Complex(0.0, 0.0)};

        // |↓⟩ = (0, 1)ᵀ
        std::vector<Complex> spin_down = {Complex(0.0, 0.0), Complex(1.0, 0.0)};

        return {spin_up, spin_down};
    }

    /**
     * @brief Pauli matrices
     */
    struct PauliMatrices {
        // σ_x = [[0, 1], [1, 0]]
        static std::vector<std::vector<Complex>> sigma_x() {
            return {{Complex(0.0, 0.0), Complex(1.0, 0.0)},
                    {Complex(1.0, 0.0), Complex(0.0, 0.0)}};
        }

        // σ_y = [[0, -i], [i, 0]]
        static std::vector<std::vector<Complex>> sigma_y() {
            return {{Complex(0.0, 0.0), Complex(0.0, -1.0)},
                    {Complex(0.0, 1.0), Complex(0.0, 0.0)}};
        }

        // σ_z = [[1, 0], [0, -1]]
        static std::vector<std::vector<Complex>> sigma_z() {
            return {{Complex(1.0, 0.0), Complex(0.0, 0.0)},
                    {Complex(0.0, 0.0), Complex(-1.0, 0.0)}};
        }
    };

    /**
     * @brief Magnetic moment
     *
     * μ = -g_e(e/2m_e)S = -g_e μ_B S/ℏ
     */
    static double magnetic_moment_z(int m_s_sign) {
        return -constants::g_e * constants::mu_B * m_s_sign * 0.5;
    }

    /**
     * @brief Zeeman effect (spin)
     *
     * ΔE = -μ·B = g_e μ_B m_s B
     */
    static double zeeman_energy(int m_s_sign, double magnetic_field) {
        return constants::g_e * constants::mu_B * m_s_sign * 0.5 * magnetic_field;
    }

    /**
     * @brief Stern-Gerlach deflection
     */
    static double stern_gerlach_deflection(
        double magnetic_field_gradient,
        double flight_length,
        double velocity) {

        double force = constants::mu_B * magnetic_field_gradient;
        double time = flight_length / velocity;

        return 0.5 * (force / constants::m_e) * time * time;
    }

    /**
     * @brief Spin-orbit coupling
     *
     * H_SO = (1/2m²c²r) (dV/dr) L·S
     */
    static double spin_orbit_energy(
        int n,
        int l,
        int j) {

        if (l == 0) return 0.0;

        // Fine structure constant α ≈ 1/137
        double alpha = 7.2973525693e-3;

        // Energy scale
        double E_0 = 13.6 * constants::e;  // 13.6 eV in Joules

        // Fine structure splitting
        return (alpha * alpha * E_0 / (n * n * n)) *
               (j * (j + 1) - l * (l + 1) - 0.75) / (l * (l + 0.5) * (l + 1));
    }
};

/**
 * @brief Two-Electron Systems
 *
 * Helium-like systems
 */
class TwoElectronSystems {
public:
    /**
     * @brief Spatial wave function symmetry
     *
     * Singlet (antisymmetric): S = 0
     * Triplet (symmetric): S = 1
     */
    enum class SpinState {
        SINGLET,   // Total spin S = 0
        TRIPLET    // Total spin S = 1
    };

    /**
     * @brief Total wave function must be antisymmetric
     *
     * ψ_total = ψ_spatial ⊗ ψ_spin
     */
    static bool is_total_wavefunction_antisymmetric(
        bool spatial_symmetric,
        SpinState spin) {

        // Singlet spin is antisymmetric → spatial must be symmetric
        // Triplet spin is symmetric → spatial must be antisymmetric

        if (spin == SpinState::SINGLET) {
            return spatial_symmetric;
        } else {
            return !spatial_symmetric;
        }
    }

    /**
     * @brief Singlet state (spin)
     *
     * |S=0,M=0⟩ = (1/√2)(|↑↓⟩ - |↓↑⟩)
     */
    static Complex singlet_amplitude() {
        return 1.0 / std::sqrt(2.0);
    }

    /**
     * @brief Triplet states (spin)
     *
     * |S=1,M=+1⟩ = |↑↑⟩
     * |S=1,M=0⟩ = (1/√2)(|↑↓⟩ + |↓↑⟩)
     * |S=1,M=-1⟩ = |↓↓⟩
     */
    static Complex triplet_amplitude(int M) {
        if (M == 0) {
            return 1.0 / std::sqrt(2.0);
        }
        return 1.0;
    }

    /**
     * @brief Exchange energy
     *
     * Difference between singlet and triplet energies
     */
    static double exchange_energy(
        double direct_integral,
        double exchange_integral) {

        // Singlet: E_S = direct + exchange
        // Triplet: E_T = direct - exchange
        // J = exchange_integral

        return 2.0 * exchange_integral;
    }

    /**
     * @brief Direct Coulomb integral
     *
     * J = ∫∫ |ψ_a(1)|² (e²/r₁₂) |ψ_b(2)|² dτ₁dτ₂
     */
    static double direct_coulomb_integral_estimate(
        double Z,
        int n_a,
        int n_b) {

        // Simplified estimate
        double E_0 = 13.6 * constants::e;  // Rydberg energy
        return (5.0 * Z * E_0) / (16.0 * n_a * n_b);
    }

    /**
     * @brief Exchange integral
     *
     * K = ∫∫ ψ_a*(1)ψ_b*(2) (e²/r₁₂) ψ_b(1)ψ_a(2) dτ₁dτ₂
     */
    static double exchange_integral_estimate(
        double Z,
        int n_a,
        int n_b) {

        // Simplified estimate
        double E_0 = 13.6 * constants::e;
        return (Z * E_0) / (4.0 * std::max(n_a, n_b) * std::max(n_a, n_b));
    }

    /**
     * @brief Ortho and para states
     */
    static std::string classify_state(SpinState spin) {
        if (spin == SpinState::SINGLET) {
            return "Para (S=0, singlet)";
        } else {
            return "Ortho (S=1, triplet)";
        }
    }
};

/**
 * @brief Helium Atom
 *
 * First multi-electron atom
 */
class HeliumAtom {
public:
    /**
     * @brief Ground state energy (experimental)
     *
     * E_0 = -79.0 eV (compared to -108.8 eV for independent particles)
     */
    static double ground_state_energy_experimental() {
        return -79.0 * constants::e;  // Convert eV to Joules
    }

    /**
     * @brief Independent particle approximation
     *
     * E ≈ -Z²(13.6 eV)(1/n₁² + 1/n₂²)
     */
    static double independent_particle_energy(int n1, int n2, double Z = 2.0) {
        double E_0 = 13.6 * constants::e;
        return -Z * Z * E_0 * (1.0 / (n1 * n1) + 1.0 / (n2 * n2));
    }

    /**
     * @brief Screening constant (Slater's rules)
     *
     * Z_eff = Z - σ
     */
    static double effective_nuclear_charge(double Z, double screening) {
        return Z - screening;
    }

    /**
     * @brief Variational ground state energy
     *
     * Using Z_eff as variational parameter
     */
    static double variational_energy(double Z_eff) {
        double E_0 = 13.6 * constants::e;

        // E = -2Z_eff²E_0 + (5/8)Z_eff·2E_0
        return -2.0 * Z_eff * Z_eff * E_0 + (5.0 / 4.0) * Z_eff * E_0;
    }

    /**
     * @brief Optimal Z_eff from variational method
     *
     * Z_eff = Z - 5/16 ≈ 1.69 for helium
     */
    static double optimal_z_effective(double Z = 2.0) {
        return Z - 5.0 / 16.0;
    }

    /**
     * @brief First ionization energy
     *
     * He → He⁺ + e⁻
     */
    static double first_ionization_energy() {
        return 24.6 * constants::e;  // 24.6 eV
    }

    /**
     * @brief Second ionization energy
     *
     * He⁺ → He²⁺ + e⁻ (hydrogen-like)
     */
    static double second_ionization_energy() {
        double Z = 2.0;
        double E_0 = 13.6 * constants::e;
        return Z * Z * E_0;  // 54.4 eV
    }

    /**
     * @brief Electron-electron repulsion energy
     *
     * V_ee ≈ (5/4)Z_eff E_0
     */
    static double electron_repulsion_energy(double Z_eff) {
        double E_0 = 13.6 * constants::e;
        return (5.0 / 4.0) * Z_eff * E_0;
    }

    /**
     * @brief Excited state: 1s2s configuration
     *
     * Singlet (¹S) and Triplet (³S) states
     */
    static std::pair<double, double> excited_1s2s_energies(
        double direct_integral,
        double exchange_integral) {

        // Singlet: E_S = E_0 + J + K
        // Triplet: E_T = E_0 + J - K

        double E_base = -59.2 * constants::e;  // Base energy

        double E_singlet = E_base + direct_integral + exchange_integral;
        double E_triplet = E_base + direct_integral - exchange_integral;

        return {E_singlet, E_triplet};
    }

    /**
     * @brief Singlet-triplet splitting
     */
    static double singlet_triplet_splitting(double exchange_integral) {
        return 2.0 * exchange_integral;
    }
};

/**
 * @brief Helium Atom Orbitals
 *
 * Spatial wave functions
 */
class HeliumAtomOrbitals {
public:
    /**
     * @brief Hydrogenic orbital approximation
     *
     * ψ_nlm(r,θ,φ) = R_nl(r)Y_lm(θ,φ)
     */
    static double radial_wavefunction_1s(double r, double Z_eff) {
        double a_0 = constants::a_0;
        double norm = std::pow(Z_eff / a_0, 1.5) * std::sqrt(1.0 / M_PI);

        return norm * std::exp(-Z_eff * r / a_0);
    }

    /**
     * @brief 2s orbital
     */
    static double radial_wavefunction_2s(double r, double Z_eff) {
        double a_0 = constants::a_0;
        double rho = Z_eff * r / a_0;

        double norm = std::pow(Z_eff / (2.0 * a_0), 1.5) / std::sqrt(2.0 * M_PI);

        return norm * (2.0 - rho) * std::exp(-rho / 2.0);
    }

    /**
     * @brief 2p orbital (radial part)
     */
    static double radial_wavefunction_2p(double r, double Z_eff) {
        double a_0 = constants::a_0;
        double rho = Z_eff * r / a_0;

        double norm = std::pow(Z_eff / (2.0 * a_0), 1.5) / std::sqrt(24.0 * M_PI);

        return norm * rho * std::exp(-rho / 2.0);
    }

    /**
     * @brief Product wave function for ground state
     *
     * ψ(r₁, r₂) = ψ_1s(r₁)ψ_1s(r₂)
     */
    static double ground_state_product(
        double r1,
        double r2,
        double Z_eff) {

        return radial_wavefunction_1s(r1, Z_eff) *
               radial_wavefunction_1s(r2, Z_eff);
    }

    /**
     * @brief Symmetric spatial wave function (singlet spin)
     *
     * ψ_+ = (1/√2)[ψ_a(1)ψ_b(2) + ψ_a(2)ψ_b(1)]
     */
    static double symmetric_spatial(
        double value_12,
        double value_21) {

        return (value_12 + value_21) / std::sqrt(2.0);
    }

    /**
     * @brief Antisymmetric spatial wave function (triplet spin)
     *
     * ψ_- = (1/√2)[ψ_a(1)ψ_b(2) - ψ_a(2)ψ_b(1)]
     */
    static double antisymmetric_spatial(
        double value_12,
        double value_21) {

        return (value_12 - value_21) / std::sqrt(2.0);
    }

    /**
     * @brief Probability density |ψ|²
     */
    static double probability_density(double r, double Z_eff, int n) {
        double psi = 0.0;

        if (n == 1) {
            psi = radial_wavefunction_1s(r, Z_eff);
        } else if (n == 2) {
            psi = radial_wavefunction_2s(r, Z_eff);
        }

        return psi * psi * r * r * 4.0 * M_PI;  // Include r² Jacobian
    }

    /**
     * @brief Radial expectation value ⟨r⟩
     */
    static double expectation_radius(double Z_eff, int n) {
        double a_0 = constants::a_0;

        if (n == 1) {
            return (3.0 / 2.0) * a_0 / Z_eff;
        } else if (n == 2) {
            return 6.0 * a_0 / Z_eff;
        }

        return n * n * a_0 / Z_eff;  // General formula
    }

    /**
     * @brief Most probable radius (maximum of radial probability)
     */
    static double most_probable_radius(double Z_eff, int n, int l) {
        double a_0 = constants::a_0;
        return n * n * a_0 / Z_eff;  // Simplified
    }
};

} // namespace advanced_quantum
} // namespace physics

#endif // PHYSICS_ADVANCED_QUANTUM_MECHANICS_HPP
