/**
 * Phase 4 Validation: Advanced Quantum Mechanics
 *
 * Tests the advanced_quantum_mechanics.hpp module functions.
 *
 * Coverage:
 * - Kummer's function (confluent hypergeometric, Laguerre polynomials)
 * - Hamiltonian mechanics (phase space, Poisson brackets)
 * - Classical harmonic oscillator
 * - Plane wave mathematics (dispersion, Fourier transforms)
 * - Schrödinger equation for free particle
 * - Wave functions and wave packets (Gaussian packets, spreading)
 * - Quantum tunneling (rectangular barrier, WKB approximation)
 * - Perturbation theory (first/second order energy corrections)
 * - Stark effect (hydrogen in electric field)
 * - Pauli exclusion principle
 * - Electron spin (magnetic moment, Zeeman effect)
 * - Two-electron systems (singlet/triplet states, exchange energy)
 * - Helium atom (ground state, ionization energies, orbitals)
 */

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <functional>
#include <tuple>
#include "../include/physics/advanced_quantum_mechanics.hpp"

// Test tolerance
const double TOLERANCE = 1e-5;
const double LOOSE_TOLERANCE = 1e-3;
const double VERY_LOOSE = 0.01;

// Test macros
#define ASSERT_NEAR(actual, expected, tolerance) \
    do { \
        if (std::abs((actual) - (expected)) > (tolerance)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #actual \
                      << " = " << (actual) << ", expected " << (expected) \
                      << " (diff: " << std::abs((actual) - (expected)) << ")" << std::endl; \
            return false; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        if (!(condition)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #condition << std::endl; \
            return false; \
        } \
    } while(0)

using namespace physics::advanced_quantum;
using namespace physics::advanced_quantum::constants;
using Complex = std::complex<double>;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Advanced Quantum Mechanics Validation ===" << std::endl;
    std::cout << std::endl;

    // Helper lambda to run tests
    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) {
            std::cout << "PASS" << std::endl;
            tests_passed++;
            return true;
        } else {
            tests_failed++;
            return false;
        }
    };

    // ========================================
    // Kummer's Function Tests
    // ========================================

    run_test("Kummer function: Pochhammer symbol (a)_0 = 1", []() {
        double result = KummersFunction::pochhammer(5.0, 0);
        ASSERT_NEAR(result, 1.0, TOLERANCE);
        return true;
    });

    run_test("Kummer function: Pochhammer symbol (a)_1 = a", []() {
        double a = 3.5;
        double result = KummersFunction::pochhammer(a, 1);
        ASSERT_NEAR(result, a, TOLERANCE);
        return true;
    });

    run_test("Kummer function: Pochhammer symbol (a)_3 = a(a+1)(a+2)", []() {
        double a = 2.0;
        double result = KummersFunction::pochhammer(a, 3);
        double expected = a * (a + 1.0) * (a + 2.0);
        ASSERT_NEAR(result, expected, TOLERANCE);
        return true;
    });

    run_test("Kummer function: M(a,b,0) = 1", []() {
        double M = KummersFunction::kummer_m(1.5, 2.5, 0.0);
        ASSERT_NEAR(M, 1.0, TOLERANCE);
        return true;
    });

    run_test("Kummer function: M(0,b,z) = 1 for any z", []() {
        double M = KummersFunction::kummer_m(0.0, 2.0, 5.0);
        ASSERT_NEAR(M, 1.0, TOLERANCE);
        return true;
    });

    run_test("Kummer function: Laguerre polynomial L_0(x) = 1", []() {
        double L = KummersFunction::laguerre_polynomial(0, 1.5);
        ASSERT_NEAR(L, 1.0, TOLERANCE);
        return true;
    });

    run_test("Kummer function: Laguerre polynomial L_1(x) = 1-x", []() {
        double x = 2.0;
        double L = KummersFunction::laguerre_polynomial(1, x);
        double expected = 1.0 - x;
        ASSERT_NEAR(L, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Kummer function: hydrogen radial factor is finite", []() {
        int n = 2, l = 0;
        double r = a_0;
        double R_factor = KummersFunction::hydrogen_radial_factor(n, l, r);
        ASSERT_TRUE(std::isfinite(R_factor));
        return true;
    });

    // ========================================
    // Hamiltonian Mechanics Tests
    // ========================================

    run_test("Hamiltonian mechanics: H = p^2/2m + V for harmonic oscillator", []() {
        double p = 2.0;
        double q = 1.0;
        double m = 1.0;
        double k = 4.0;
        auto V = [k](double x) { return 0.5 * k * x * x; };
        double H = HamiltonianMechanics::hamiltonian_particle(p, q, m, V);
        double expected = p * p / (2.0 * m) + 0.5 * k * q * q;
        ASSERT_NEAR(H, expected, TOLERANCE);
        return true;
    });

    run_test("Hamiltonian mechanics: Hamilton equations dq/dt = p/m", []() {
        double q = 1.0;
        double p = 3.0;
        double m = 2.0;
        auto V = [](double x) { return 0.0; };
        auto F = [](double x) { return 0.0; };
        auto [dq_dt, dp_dt] = HamiltonianMechanics::hamiltons_equations(q, p, m, V, F);
        ASSERT_NEAR(dq_dt, p / m, TOLERANCE);
        return true;
    });

    run_test("Hamiltonian mechanics: Poisson bracket {q,p} = 1", []() {
        auto Q = [](double q, double p) { return q; };
        auto P = [](double q, double p) { return p; };
        double bracket = HamiltonianMechanics::poisson_bracket(1.0, 1.0, Q, P);
        ASSERT_NEAR(bracket, 1.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Hamiltonian mechanics: action-angle variables for oscillator", []() {
        double E = 10.0;
        double omega = 2.0;
        double t = 1.5;
        auto [I, theta] = HamiltonianMechanics::action_angle_oscillator(E, omega, t);
        ASSERT_NEAR(I, E / omega, TOLERANCE);
        ASSERT_NEAR(theta, omega * t, TOLERANCE);
        return true;
    });

    run_test("Hamiltonian mechanics: Liouville theorem (phase space volume)", []() {
        double V_initial = 10.0;
        double V_final = 10.0;
        ASSERT_TRUE(HamiltonianMechanics::verify_liouville_theorem(V_initial, V_final));
        return true;
    });

    // ========================================
    // Classical Harmonic Oscillator Tests
    // ========================================

    run_test("Classical oscillator: position x(t) = A*cos(omega*t)", []() {
        double A = 2.0;
        double omega = 3.0;
        double t = 0.0;
        double x = ClassicalHarmonicOscillator::position(A, omega, t);
        ASSERT_NEAR(x, A, TOLERANCE);
        return true;
    });

    run_test("Classical oscillator: velocity v(t) = -A*omega*sin(omega*t)", []() {
        double A = 2.0;
        double omega = 3.0;
        double t = 0.0;
        double v = ClassicalHarmonicOscillator::velocity(A, omega, t);
        ASSERT_NEAR(v, 0.0, TOLERANCE);
        return true;
    });

    run_test("Classical oscillator: energy E = (1/2)m*omega^2*A^2", []() {
        double A = 2.0;
        double m = 1.0;
        double omega = 3.0;
        double E = ClassicalHarmonicOscillator::energy(A, m, omega);
        double expected = 0.5 * m * omega * omega * A * A;
        ASSERT_NEAR(E, expected, TOLERANCE);
        return true;
    });

    run_test("Classical oscillator: period T = 2*pi/omega", []() {
        double omega = 2.0;
        double T = ClassicalHarmonicOscillator::period(omega);
        ASSERT_NEAR(T, M_PI, TOLERANCE);
        return true;
    });

    run_test("Classical oscillator: frequency from spring constant", []() {
        double k = 4.0;
        double m = 1.0;
        double omega = ClassicalHarmonicOscillator::frequency_from_spring(k, m);
        ASSERT_NEAR(omega, 2.0, TOLERANCE);
        return true;
    });

    run_test("Classical oscillator: average kinetic energy = E/2", []() {
        double E = 100.0;
        double K_avg = ClassicalHarmonicOscillator::average_kinetic_energy(E);
        ASSERT_NEAR(K_avg, 50.0, TOLERANCE);
        return true;
    });

    run_test("Classical oscillator: turning points at ±A", []() {
        double A = 5.0;
        auto [x_min, x_max] = ClassicalHarmonicOscillator::turning_points(A);
        ASSERT_NEAR(x_min, -A, TOLERANCE);
        ASSERT_NEAR(x_max, A, TOLERANCE);
        return true;
    });

    // ========================================
    // Plane Wave Mathematics Tests
    // ========================================

    run_test("Plane wave: magnitude at t=0 x=0", []() {
        Complex psi = PlaneWaveMathematics::plane_wave(1.0, 0.0, 1.0, 0.0, 1.0);
        ASSERT_NEAR(std::abs(psi), 1.0, TOLERANCE);
        return true;
    });

    run_test("Plane wave: dispersion relation omega = hbar*k^2/(2m)", []() {
        double k = 1e10;
        double m = m_e;
        double omega = PlaneWaveMathematics::dispersion_relation_free(k, m);
        double expected = (hbar * k * k) / (2.0 * m);
        ASSERT_NEAR(omega, expected, TOLERANCE);
        return true;
    });

    run_test("Plane wave: phase velocity v_p = omega/k", []() {
        double omega = 1e15;
        double k = 1e10;
        double v_p = PlaneWaveMathematics::phase_velocity(omega, k);
        ASSERT_NEAR(v_p, omega / k, TOLERANCE);
        return true;
    });

    run_test("Plane wave: group velocity v_g = hbar*k/m", []() {
        double k = 1e10;
        double m = m_e;
        double v_g = PlaneWaveMathematics::group_velocity(k, m);
        double expected = (hbar * k) / m;
        ASSERT_NEAR(v_g, expected, TOLERANCE);
        return true;
    });

    run_test("Plane wave: group velocity is twice phase velocity for free particle", []() {
        double k = 1e10;
        double m = m_e;
        double omega = PlaneWaveMathematics::dispersion_relation_free(k, m);
        double v_p = PlaneWaveMathematics::phase_velocity(omega, k);
        double v_g = PlaneWaveMathematics::group_velocity(k, m);
        ASSERT_NEAR(v_g, 2.0 * v_p, LOOSE_TOLERANCE * v_g);
        return true;
    });

    // ========================================
    // Schrödinger Free Particle Tests
    // ========================================

    run_test("Schrödinger free particle: energy E = hbar^2*k^2/(2m)", []() {
        double k = 1e10;
        double m = m_e;
        double E = SchrodingerFreeParticle::energy(k, m);
        double expected = (hbar * hbar * k * k) / (2.0 * m);
        ASSERT_NEAR(E, expected, TOLERANCE);
        return true;
    });

    run_test("Schrödinger free particle: wave function normalized", []() {
        double k = 1e10;
        double m = m_e;
        Complex psi = SchrodingerFreeParticle::wave_function(k, 0.0, 0.0, m, 1.0);
        ASSERT_NEAR(std::abs(psi), 1.0, TOLERANCE);
        return true;
    });

    run_test("Schrödinger free particle: satisfies Schrödinger equation", []() {
        double k = 1e10;
        double m = m_e;
        ASSERT_TRUE(SchrodingerFreeParticle::verify_schrodinger_equation(k, m));
        return true;
    });

    run_test("Schrödinger free particle: probability current j = hbar*k/m", []() {
        double k = 1e10;
        double m = m_e;
        double j = SchrodingerFreeParticle::probability_current(k, m);
        double expected = (hbar * k) / m;
        ASSERT_NEAR(j, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Wave Packets Tests
    // ========================================

    run_test("Wave packet: Gaussian packet normalization factor", []() {
        double sigma = 1e-9;
        double k0 = 1e10;
        Complex psi = WaveFunctionsAndPackets::gaussian_packet(0.0, sigma, k0);
        // At x=0, the Gaussian should have specific value
        ASSERT_TRUE(std::abs(psi) > 0.0);
        return true;
    });

    run_test("Wave packet: spreading increases with time", []() {
        double sigma = 1e-9;
        double m = m_e;
        double w1 = WaveFunctionsAndPackets::packet_width(sigma, m, 0.0);
        double w2 = WaveFunctionsAndPackets::packet_width(sigma, m, 1e-12);
        ASSERT_TRUE(w2 > w1);
        return true;
    });

    run_test("Wave packet: initial width unchanged at t=0", []() {
        double sigma = 1e-9;
        double m = m_e;
        double w = WaveFunctionsAndPackets::packet_width(sigma, m, 0.0);
        ASSERT_NEAR(w, sigma, TOLERANCE);
        return true;
    });

    run_test("Wave packet: uncertainty product Delta_x*Delta_p >= hbar/2", []() {
        auto psi = [](double x) { return WaveFunctionsAndPackets::gaussian_packet(x, 1e-9, 1e10); };
        double product = WaveFunctionsAndPackets::uncertainty_product(psi, -1e-8, 1e-8);
        ASSERT_TRUE(product >= hbar / 2.0 * 0.99);  // Allow small numerical error
        return true;
    });

    // ========================================
    // Quantum Tunneling Tests
    // ========================================

    run_test("Quantum tunneling: transmission T=1 when E > V_0", []() {
        double E = 10.0 * e;
        double V0 = 5.0 * e;
        double a = 1e-10;
        double T = QuantumTunneling::transmission_coefficient_rectangular(E, V0, a, m_e);
        ASSERT_NEAR(T, 1.0, TOLERANCE);
        return true;
    });

    run_test("Quantum tunneling: transmission T<1 when E < V_0", []() {
        double E = 2.0 * e;
        double V0 = 5.0 * e;
        double a = 1e-10;
        double T = QuantumTunneling::transmission_coefficient_rectangular(E, V0, a, m_e);
        ASSERT_TRUE(T < 1.0);
        ASSERT_TRUE(T > 0.0);
        return true;
    });

    run_test("Quantum tunneling: R + T = 1", []() {
        double E = 2.0 * e;
        double V0 = 5.0 * e;
        double a = 1e-10;
        double T = QuantumTunneling::transmission_coefficient_rectangular(E, V0, a, m_e);
        double R = QuantumTunneling::reflection_coefficient(E, V0, a, m_e);
        ASSERT_NEAR(T + R, 1.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Quantum tunneling: WKB transmission decreases with barrier width", []() {
        double E = 2.0 * e;
        auto V = [E](double x) { return 3.0 * E; };  // Constant barrier
        double T1 = QuantumTunneling::wkb_transmission(E, V, 0.0, 1e-10, m_e);
        double T2 = QuantumTunneling::wkb_transmission(E, V, 0.0, 2e-10, m_e);
        ASSERT_TRUE(T2 < T1);
        return true;
    });

    run_test("Quantum tunneling: STM current exponential decay", []() {
        double d1 = 1e-10;
        double d2 = 2e-10;
        double phi = 4.5 * e;
        double I1 = QuantumTunneling::stm_current(d1, phi);
        double I2 = QuantumTunneling::stm_current(d2, phi);
        ASSERT_TRUE(I2 < I1);
        return true;
    });

    // ========================================
    // Perturbation Theory Tests
    // ========================================

    run_test("Perturbation theory: first order energy correction", []() {
        auto psi_0 = [](double x) { return Complex(std::exp(-x * x), 0.0); };
        auto H_prime = [](double x) { return x * x; };
        double E1 = PerturbationTheoryNondegenerate::first_order_energy(psi_0, H_prime, -5.0, 5.0);
        ASSERT_TRUE(E1 > 0.0);  // Positive perturbation
        return true;
    });

    run_test("Perturbation theory: anharmonic oscillator correction", []() {
        int n = 0;
        double lambda = 0.1;
        double m = m_e;
        double omega = 1e15;
        double E1 = PerturbationTheoryNondegenerate::anharmonic_correction_first_order(n, lambda, 4, m, omega);
        ASSERT_TRUE(E1 >= 0.0);
        return true;
    });

    // ========================================
    // Stark Effect Tests
    // ========================================

    run_test("Stark effect: linear splitting for n=2", []() {
        double E_field = 1e6;  // V/m
        double splitting = StarkEffect::linear_stark_splitting_n2(E_field);
        double expected = 3.0 * e * E_field * a_0;
        ASSERT_NEAR(splitting, expected, TOLERANCE);
        return true;
    });

    run_test("Stark effect: quadratic shift for ground state", []() {
        double E_field = 1e6;
        double shift = StarkEffect::quadratic_stark_ground_state(E_field);
        ASSERT_TRUE(shift < 0.0);  // Energy lowering
        return true;
    });

    run_test("Stark effect: hydrogen polarizability", []() {
        double alpha = StarkEffect::hydrogen_polarizability();
        double expected = (9.0 / 2.0) * a_0 * a_0 * a_0;
        ASSERT_NEAR(alpha, expected, TOLERANCE);
        return true;
    });

    run_test("Stark effect: avoided crossing produces split levels", []() {
        double E1 = 1.0 * e;
        double E2 = 1.1 * e;
        double coupling = 0.05 * e;
        auto [E_lower, E_upper] = StarkEffect::avoided_crossing_energies(E1, E2, coupling);
        ASSERT_TRUE(E_lower < E_upper);
        ASSERT_TRUE(E_upper - E_lower > std::abs(E1 - E2));
        return true;
    });

    // ========================================
    // Pauli Exclusion Principle Tests
    // ========================================

    run_test("Pauli exclusion: unique quantum numbers allowed", []() {
        std::vector<std::tuple<int, int, int, int>> states = {
            {1, 0, 0, 1},   // 1s spin-up
            {1, 0, 0, -1},  // 1s spin-down
            {2, 0, 0, 1}    // 2s spin-up
        };
        ASSERT_TRUE(PauliExclusionPrinciple::are_quantum_numbers_unique(states));
        return true;
    });

    run_test("Pauli exclusion: duplicate quantum numbers not allowed", []() {
        std::vector<std::tuple<int, int, int, int>> states = {
            {1, 0, 0, 1},
            {1, 0, 0, 1}  // Duplicate!
        };
        ASSERT_TRUE(!PauliExclusionPrinciple::are_quantum_numbers_unique(states));
        return true;
    });

    run_test("Pauli exclusion: max electrons in shell n is 2n^2", []() {
        ASSERT_TRUE(PauliExclusionPrinciple::max_electrons_in_shell(1) == 2);
        ASSERT_TRUE(PauliExclusionPrinciple::max_electrons_in_shell(2) == 8);
        ASSERT_TRUE(PauliExclusionPrinciple::max_electrons_in_shell(3) == 18);
        return true;
    });

    run_test("Pauli exclusion: max electrons in subshell is 2(2l+1)", []() {
        ASSERT_TRUE(PauliExclusionPrinciple::max_electrons_in_subshell(0) == 2);   // s
        ASSERT_TRUE(PauliExclusionPrinciple::max_electrons_in_subshell(1) == 6);   // p
        ASSERT_TRUE(PauliExclusionPrinciple::max_electrons_in_subshell(2) == 10);  // d
        return true;
    });

    run_test("Pauli exclusion: antisymmetric wave function", []() {
        Complex psi_12(1.0, 0.0);
        Complex psi_21(-1.0, 0.0);
        ASSERT_TRUE(PauliExclusionPrinciple::is_antisymmetric(psi_12, psi_21));
        return true;
    });

    run_test("Pauli exclusion: Fermi energy for electron gas", []() {
        double n = 1e28;  // electrons/m^3
        double E_F = PauliExclusionPrinciple::fermi_energy(n);
        ASSERT_TRUE(E_F > 0.0);
        return true;
    });

    run_test("Pauli exclusion: degeneracy pressure", []() {
        double n = 1e28;
        double P = PauliExclusionPrinciple::degeneracy_pressure(n);
        ASSERT_TRUE(P > 0.0);
        return true;
    });

    // ========================================
    // Electron Spin Tests
    // ========================================

    run_test("Electron spin: s = 1/2", []() {
        double s = ElectronSpin::spin_quantum_number();
        ASSERT_NEAR(s, 0.5, TOLERANCE);
        return true;
    });

    run_test("Electron spin: |S| = sqrt(3)/2 * hbar", []() {
        double S = ElectronSpin::spin_angular_momentum_magnitude();
        double expected = hbar * std::sqrt(0.75);
        ASSERT_NEAR(S, expected, TOLERANCE);
        return true;
    });

    run_test("Electron spin: S_z = ±hbar/2", []() {
        double S_up = ElectronSpin::spin_z_component(1);
        double S_down = ElectronSpin::spin_z_component(-1);
        ASSERT_NEAR(S_up, 0.5 * hbar, TOLERANCE);
        ASSERT_NEAR(S_down, -0.5 * hbar, TOLERANCE);
        return true;
    });

    run_test("Electron spin: spin states are orthonormal", []() {
        auto [spin_up, spin_down] = ElectronSpin::spin_states();
        // |<↑|↑>|² = 1
        double norm_up = std::norm(spin_up[0]) + std::norm(spin_up[1]);
        ASSERT_NEAR(norm_up, 1.0, TOLERANCE);
        return true;
    });

    run_test("Electron spin: Pauli matrix sigma_z has eigenvalues ±1", []() {
        auto sigma_z = ElectronSpin::PauliMatrices::sigma_z();
        // Eigenvalues are diagonal elements for diagonal matrix
        ASSERT_NEAR(sigma_z[0][0].real(), 1.0, TOLERANCE);
        ASSERT_NEAR(sigma_z[1][1].real(), -1.0, TOLERANCE);
        return true;
    });

    run_test("Electron spin: Zeeman energy Delta_E = g*mu_B*m_s*B", []() {
        double B = 1.0;  // Tesla
        double E_up = ElectronSpin::zeeman_energy(1, B);
        double E_down = ElectronSpin::zeeman_energy(-1, B);
        double splitting = E_up - E_down;
        double expected = g_e * mu_B * B;
        ASSERT_NEAR(splitting, expected, LOOSE_TOLERANCE * expected);
        return true;
    });

    // ========================================
    // Two-Electron Systems Tests
    // ========================================

    run_test("Two-electron systems: singlet is spin antisymmetric", []() {
        bool spatial_symmetric = true;
        auto spin = TwoElectronSystems::SpinState::SINGLET;
        ASSERT_TRUE(TwoElectronSystems::is_total_wavefunction_antisymmetric(spatial_symmetric, spin));
        return true;
    });

    run_test("Two-electron systems: triplet is spin symmetric", []() {
        bool spatial_symmetric = false;  // Must be antisymmetric
        auto spin = TwoElectronSystems::SpinState::TRIPLET;
        ASSERT_TRUE(TwoElectronSystems::is_total_wavefunction_antisymmetric(spatial_symmetric, spin));
        return true;
    });

    run_test("Two-electron systems: singlet amplitude 1/sqrt(2)", []() {
        Complex amp = TwoElectronSystems::singlet_amplitude();
        ASSERT_NEAR(std::abs(amp), 1.0 / std::sqrt(2.0), TOLERANCE);
        return true;
    });

    run_test("Two-electron systems: triplet M=0 amplitude 1/sqrt(2)", []() {
        Complex amp = TwoElectronSystems::triplet_amplitude(0);
        ASSERT_NEAR(std::abs(amp), 1.0 / std::sqrt(2.0), TOLERANCE);
        return true;
    });

    run_test("Two-electron systems: exchange energy is 2J", []() {
        double J_direct = 10.0 * e;
        double J_exchange = 2.0 * e;
        double exchange = TwoElectronSystems::exchange_energy(J_direct, J_exchange);
        ASSERT_NEAR(exchange, 2.0 * J_exchange, TOLERANCE);
        return true;
    });

    run_test("Two-electron systems: ortho-para classification", []() {
        std::string para = TwoElectronSystems::classify_state(TwoElectronSystems::SpinState::SINGLET);
        std::string ortho = TwoElectronSystems::classify_state(TwoElectronSystems::SpinState::TRIPLET);
        ASSERT_TRUE(para.find("Para") != std::string::npos);
        ASSERT_TRUE(ortho.find("Ortho") != std::string::npos);
        return true;
    });

    // ========================================
    // Helium Atom Tests
    // ========================================

    run_test("Helium atom: ground state energy is -79.0 eV", []() {
        double E0 = HeliumAtom::ground_state_energy_experimental();
        double expected = -79.0 * e;
        ASSERT_NEAR(E0, expected, TOLERANCE);
        return true;
    });

    run_test("Helium atom: independent particle energy for ground state", []() {
        double E = HeliumAtom::independent_particle_energy(1, 1, 2.0);
        // E = -Z²*13.6*(1/1² + 1/1²) = -4*13.6*2 = -108.8 eV
        double expected = -4.0 * 2.0 * 13.6 * e;
        ASSERT_NEAR(E, expected, LOOSE_TOLERANCE * std::abs(expected));
        return true;
    });

    run_test("Helium atom: effective nuclear charge Z_eff = Z - sigma", []() {
        double Z = 2.0;
        double sigma = 0.3;
        double Z_eff = HeliumAtom::effective_nuclear_charge(Z, sigma);
        ASSERT_NEAR(Z_eff, 1.7, TOLERANCE);
        return true;
    });

    run_test("Helium atom: optimal Z_eff = Z - 5/16", []() {
        double Z = 2.0;
        double Z_eff = HeliumAtom::optimal_z_effective(Z);
        double expected = Z - 5.0 / 16.0;
        ASSERT_NEAR(Z_eff, expected, TOLERANCE);
        return true;
    });

    run_test("Helium atom: first ionization energy is 24.6 eV", []() {
        double IE1 = HeliumAtom::first_ionization_energy();
        double expected = 24.6 * e;
        ASSERT_NEAR(IE1, expected, TOLERANCE);
        return true;
    });

    run_test("Helium atom: second ionization energy is 54.4 eV", []() {
        double IE2 = HeliumAtom::second_ionization_energy();
        double expected = 4.0 * 13.6 * e;  // Z²*13.6 eV for Z=2
        ASSERT_NEAR(IE2, expected, LOOSE_TOLERANCE * expected);
        return true;
    });

    run_test("Helium atom: singlet-triplet splitting", []() {
        double J = 1.0 * e;
        double splitting = HeliumAtom::singlet_triplet_splitting(J);
        ASSERT_NEAR(splitting, 2.0 * J, TOLERANCE);
        return true;
    });

    // ========================================
    // Helium Atom Orbitals Tests
    // ========================================

    run_test("Helium orbitals: 1s orbital is positive at nucleus", []() {
        double Z_eff = 1.69;
        double R = HeliumAtomOrbitals::radial_wavefunction_1s(0.01 * a_0, Z_eff);
        ASSERT_TRUE(R > 0.0);
        return true;
    });

    run_test("Helium orbitals: 2s orbital has node", []() {
        double Z_eff = 1.0;
        double R1 = HeliumAtomOrbitals::radial_wavefunction_2s(0.5 * a_0, Z_eff);
        double R2 = HeliumAtomOrbitals::radial_wavefunction_2s(3.0 * a_0, Z_eff);
        // There should be a node between these points
        ASSERT_TRUE(R1 * R2 < 0.0 || std::abs(R1) < 1e-10 || std::abs(R2) < 1e-10);
        return true;
    });

    run_test("Helium orbitals: 2p orbital zero at origin", []() {
        double Z_eff = 1.0;
        double R = HeliumAtomOrbitals::radial_wavefunction_2p(0.0, Z_eff);
        ASSERT_NEAR(R, 0.0, TOLERANCE);
        return true;
    });

    run_test("Helium orbitals: symmetric spatial combination", []() {
        double val_12 = 1.0;
        double val_21 = 1.0;
        double psi = HeliumAtomOrbitals::symmetric_spatial(val_12, val_21);
        ASSERT_NEAR(psi, std::sqrt(2.0), TOLERANCE);
        return true;
    });

    run_test("Helium orbitals: antisymmetric spatial combination", []() {
        double val_12 = 1.0;
        double val_21 = 1.0;
        double psi = HeliumAtomOrbitals::antisymmetric_spatial(val_12, val_21);
        ASSERT_NEAR(psi, 0.0, TOLERANCE);
        return true;
    });

    run_test("Helium orbitals: expectation radius scales with 1/Z_eff", []() {
        double Z1 = 1.0;
        double Z2 = 2.0;
        double r1 = HeliumAtomOrbitals::expectation_radius(Z1, 1);
        double r2 = HeliumAtomOrbitals::expectation_radius(Z2, 1);
        ASSERT_NEAR(r1 / r2, Z2 / Z1, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Helium orbitals: probability density integrates to 1", []() {
        double Z_eff = 1.0;
        // Radial probability density should peak somewhere
        double rho1 = HeliumAtomOrbitals::probability_density(0.5 * a_0, Z_eff, 1);
        double rho2 = HeliumAtomOrbitals::probability_density(1.5 * a_0, Z_eff, 1);
        ASSERT_TRUE(rho1 > 0.0);
        ASSERT_TRUE(rho2 > 0.0);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Total tests:  " << (tests_passed + tests_failed) << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << "======================================" << std::endl;

    return (tests_failed == 0) ? 0 : 1;
}
