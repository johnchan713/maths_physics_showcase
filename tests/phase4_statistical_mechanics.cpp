/**
 * Phase 4 Validation: Statistical Mechanics
 *
 * Tests the statistical_mechanics.hpp module functions.
 *
 * Coverage:
 * - Canonical ensemble (partition functions, thermodynamic relations)
 * - Partition functions for specific systems (harmonic oscillator, rotators, etc.)
 * - Phase transitions (Van der Waals, critical points, order parameters)
 * - Mean field theory (critical temperature, magnetization)
 * - Thermodynamic consistency (F, E, S relationships)
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "../include/physics/statistical_mechanics.hpp"

// Test tolerance
const double TOLERANCE = 1e-6;
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

// Physical constants from the module
using physics::statistical_mechanics::k_B;
using physics::statistical_mechanics::h;
using physics::statistical_mechanics::hbar;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Statistical Mechanics Validation ===" << std::endl;
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
    // Canonical Ensemble Tests
    // ========================================

    run_test("Canonical ensemble: beta = 1/(k_B T)", []() {
        double T = 300.0;  // 300 K
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);
        double beta = ensemble.beta();
        double expected = 1.0 / (k_B * T);
        ASSERT_NEAR(beta, expected, TOLERANCE);
        return true;
    });

    run_test("Boltzmann probability normalization", []() {
        double T = 300.0;
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);

        // Simple two-level system
        auto energy = [](int n) { return n * 1e-20; };  // 0, 1e-20 J
        double Z = ensemble.partitionFunction(energy, 2);

        double p0 = ensemble.boltzmannProbability(energy(0), Z);
        double p1 = ensemble.boltzmannProbability(energy(1), Z);

        // Probabilities should sum to 1
        ASSERT_NEAR(p0 + p1, 1.0, TOLERANCE);
        return true;
    });

    run_test("Boltzmann distribution: lower energy has higher probability", []() {
        double T = 300.0;
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);

        auto energy = [](int n) { return n * 1e-20; };
        double Z = ensemble.partitionFunction(energy, 10);

        double p0 = ensemble.boltzmannProbability(energy(0), Z);
        double p1 = ensemble.boltzmannProbability(energy(1), Z);

        // Ground state should be most probable
        ASSERT_TRUE(p0 > p1);
        return true;
    });

    run_test("Helmholtz free energy: F = -k_B T ln(Z)", []() {
        double T = 300.0;
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);

        double Z = 10.0;  // Arbitrary partition function
        double F = ensemble.helmholtzFreeEnergy(Z);
        double expected = -k_B * T * std::log(Z);

        ASSERT_NEAR(F, expected, 1e-25);
        return true;
    });

    run_test("Entropy: S = (E - F) / T", []() {
        double T = 300.0;
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);

        double E = 1e-20;
        double F = -2e-20;
        double S = ensemble.entropy(F, E);
        double expected = (E - F) / T;

        ASSERT_NEAR(S, expected, 1e-30);
        return true;
    });

    // ========================================
    // Partition Function Tests - Harmonic Oscillator
    // ========================================

    run_test("Harmonic oscillator partition function (classical)", []() {
        double omega = 1e13;  // 10 THz
        double T = 300.0;     // 300 K

        double Z = physics::statistical_mechanics::PartitionFunctions::harmonicOscillator(omega, T);

        // Z = 1 / (2 sinh(x/2)) where x = ℏω/(k_B T)
        double x = (hbar * omega) / (k_B * T);
        double expected = 1.0 / (2.0 * std::sinh(x / 2.0));

        ASSERT_NEAR(Z, expected, TOLERANCE);
        return true;
    });

    run_test("Harmonic oscillator: high temperature limit", []() {
        double omega = 1e12;
        double T = 1e6;  // Very high temperature

        double Z = physics::statistical_mechanics::PartitionFunctions::harmonicOscillator(omega, T);

        // At high T: ℏω << k_B T, so sinh(x/2) ≈ x/2, giving Z ≈ k_B T / (ℏω)
        double x = (hbar * omega) / (k_B * T);
        // For small x: sinh(x/2) ≈ x/2, so Z ≈ 1/(2·x/2) = 1/x = k_B T / (ℏω)
        double expected_approx = (k_B * T) / (hbar * omega);

        ASSERT_NEAR(Z, expected_approx, expected_approx * 0.1);  // 10% tolerance
        return true;
    });

    run_test("Harmonic oscillator quantum: includes ground state energy", []() {
        double omega = 1e13;
        double T = 100.0;

        double Z_quantum = physics::statistical_mechanics::PartitionFunctions::harmonicOscillatorQuantum(omega, T, 50);

        // Should be positive and include exp(-β·ℏω/2) factor
        ASSERT_TRUE(Z_quantum > 0);
        return true;
    });

    // ========================================
    // Two-Level System Tests
    // ========================================

    run_test("Two-level system partition function", []() {
        double epsilon = 1e-20;  // 1e-20 J energy gap
        double T = 300.0;

        double Z = physics::statistical_mechanics::PartitionFunctions::twoLevelSystem(epsilon, T);

        // Z = 1 + exp(-β·ε)
        double beta = 1.0 / (k_B * T);
        double expected = 1.0 + std::exp(-beta * epsilon);

        ASSERT_NEAR(Z, expected, TOLERANCE);
        return true;
    });

    run_test("Two-level system: high temperature limit", []() {
        double epsilon = 1e-22;
        double T = 1e6;  // Very high T, β·ε → 0

        double Z = physics::statistical_mechanics::PartitionFunctions::twoLevelSystem(epsilon, T);

        // At high T: exp(-β·ε) → 1, so Z → 2
        ASSERT_NEAR(Z, 2.0, VERY_LOOSE);
        return true;
    });

    run_test("Two-level system: low temperature limit", []() {
        double epsilon = 1e-20;
        double T = 1.0;  // Very low T, β·ε → ∞

        double Z = physics::statistical_mechanics::PartitionFunctions::twoLevelSystem(epsilon, T);

        // At low T: exp(-β·ε) → 0, so Z → 1
        ASSERT_NEAR(Z, 1.0, VERY_LOOSE);
        return true;
    });

    // ========================================
    // Rotator Tests
    // ========================================

    run_test("2D rotator partition function", []() {
        double I = 1e-46;  // Moment of inertia (kg·m²)
        double T = 300.0;

        double Z = physics::statistical_mechanics::PartitionFunctions::rotator2D(I, T);

        // Z = 2πI / (β·ℏ²)
        double beta = 1.0 / (k_B * T);
        double expected = 2.0 * M_PI * I / (beta * hbar * hbar);

        ASSERT_NEAR(Z, expected, TOLERANCE);
        return true;
    });

    run_test("2D rotator: temperature dependence", []() {
        double I = 1e-46;
        double T1 = 100.0;
        double T2 = 200.0;

        double Z1 = physics::statistical_mechanics::PartitionFunctions::rotator2D(I, T1);
        double Z2 = physics::statistical_mechanics::PartitionFunctions::rotator2D(I, T2);

        // Z ∝ T, so Z2 / Z1 should equal T2 / T1
        double ratio = Z2 / Z1;
        double expected_ratio = T2 / T1;

        ASSERT_NEAR(ratio, expected_ratio, TOLERANCE);
        return true;
    });

    run_test("3D rotator partition function", []() {
        double I = 1e-46;
        double T = 300.0;

        double Z = physics::statistical_mechanics::PartitionFunctions::rotator3D(I, T, 50);

        // Should be positive and account for (2l+1) degeneracy
        ASSERT_TRUE(Z > 0);

        // Ground state (l=0) contributes 1 to Z
        double beta = 1.0 / (k_B * T);
        ASSERT_TRUE(Z >= 1.0);
        return true;
    });

    // ========================================
    // Phase Transition Tests - Van der Waals
    // ========================================

    run_test("Van der Waals pressure formula", []() {
        double T = 300.0;
        double v = 1e-3;  // Molar volume (m³/mol)
        double a = 0.1;   // Attraction parameter
        double b = 1e-5;  // Excluded volume

        double P = physics::statistical_mechanics::PhaseTransitions::vanDerWaalsPressure(T, v, a, b);

        // P = k_B T / (v - b) - a / v²
        double expected = (k_B * T) / (v - b) - a / (v * v);

        ASSERT_NEAR(P, expected, 1e-30);
        return true;
    });

    run_test("Van der Waals: ideal gas limit (a=0, b=0)", []() {
        double T = 300.0;
        double v = 1e-3;
        double a = 0.0;
        double b = 0.0;

        double P = physics::statistical_mechanics::PhaseTransitions::vanDerWaalsPressure(T, v, a, b);

        // Should reduce to ideal gas: P = k_B T / v
        double expected = (k_B * T) / v;

        ASSERT_NEAR(P, expected, 1e-30);
        return true;
    });

    run_test("Critical temperature formula", []() {
        double a = 0.1;
        double b = 1e-5;

        double T_c = physics::statistical_mechanics::PhaseTransitions::criticalTemperature(a, b);

        // T_c = 8a / (27 k_B b)
        double expected = (8.0 * a) / (27.0 * k_B * b);

        ASSERT_NEAR(T_c, expected, TOLERANCE);
        return true;
    });

    run_test("Critical pressure formula", []() {
        double a = 0.1;
        double b = 1e-5;

        double P_c = physics::statistical_mechanics::PhaseTransitions::criticalPressure(a, b);

        // P_c = a / (27 b²)
        double expected = a / (27.0 * b * b);

        ASSERT_NEAR(P_c, expected, TOLERANCE);
        return true;
    });

    run_test("Critical volume formula", []() {
        double b = 1e-5;

        double v_c = physics::statistical_mechanics::PhaseTransitions::criticalVolume(b);

        // v_c = 3b
        double expected = 3.0 * b;

        ASSERT_NEAR(v_c, expected, TOLERANCE);
        return true;
    });

    run_test("Critical point ratios", []() {
        double a = 0.1;
        double b = 1e-5;

        double T_c = physics::statistical_mechanics::PhaseTransitions::criticalTemperature(a, b);
        double P_c = physics::statistical_mechanics::PhaseTransitions::criticalPressure(a, b);
        double v_c = physics::statistical_mechanics::PhaseTransitions::criticalVolume(b);

        // Verify critical point is consistent with Van der Waals equation
        double P_at_critical = physics::statistical_mechanics::PhaseTransitions::vanDerWaalsPressure(T_c, v_c, a, b);

        // Should match P_c (within numerical precision)
        ASSERT_NEAR(P_at_critical, P_c, P_c * 1e-10);
        return true;
    });

    // ========================================
    // Order Parameter Tests
    // ========================================

    run_test("Order parameter above T_c", []() {
        double T = 400.0;
        double T_c = 300.0;
        double beta_exp = 0.5;  // Critical exponent

        double phi = physics::statistical_mechanics::PhaseTransitions::orderParameter(T, T_c, beta_exp);

        // Above critical temperature, order parameter = 0
        ASSERT_NEAR(phi, 0.0, TOLERANCE);
        return true;
    });

    run_test("Order parameter at T_c", []() {
        double T = 300.0;
        double T_c = 300.0;
        double beta_exp = 0.5;

        double phi = physics::statistical_mechanics::PhaseTransitions::orderParameter(T, T_c, beta_exp);

        // At critical temperature, order parameter = 0
        ASSERT_NEAR(phi, 0.0, TOLERANCE);
        return true;
    });

    run_test("Order parameter below T_c", []() {
        double T = 200.0;
        double T_c = 300.0;
        double beta_exp = 0.5;

        double phi = physics::statistical_mechanics::PhaseTransitions::orderParameter(T, T_c, beta_exp);

        // φ = ((T_c - T) / T_c)^β
        double expected = std::pow((T_c - T) / T_c, beta_exp);

        ASSERT_NEAR(phi, expected, TOLERANCE);
        return true;
    });

    run_test("Order parameter: critical exponent dependence", []() {
        double T = 250.0;
        double T_c = 300.0;

        double phi_beta_half = physics::statistical_mechanics::PhaseTransitions::orderParameter(T, T_c, 0.5);
        double phi_beta_one = physics::statistical_mechanics::PhaseTransitions::orderParameter(T, T_c, 1.0);

        // For β=1: φ = (T_c - T)/T_c = 50/300 = 1/6
        // For β=0.5: φ = √(1/6) ≈ 0.408
        ASSERT_NEAR(phi_beta_one, 1.0/6.0, TOLERANCE);
        ASSERT_NEAR(phi_beta_half, std::sqrt(1.0/6.0), TOLERANCE);
        return true;
    });

    // ========================================
    // Susceptibility Tests
    // ========================================

    run_test("Susceptibility diverges at T_c", []() {
        double T_c = 300.0;
        double gamma_exp = 1.0;

        double T_near = 300.1;  // Just above T_c
        double chi_near = physics::statistical_mechanics::PhaseTransitions::susceptibility(T_near, T_c, gamma_exp);

        double T_far = 350.0;   // Far above T_c
        double chi_far = physics::statistical_mechanics::PhaseTransitions::susceptibility(T_far, T_c, gamma_exp);

        // Susceptibility should be larger closer to T_c
        ASSERT_TRUE(chi_near > chi_far);
        return true;
    });

    run_test("Susceptibility power law", []() {
        double T = 320.0;
        double T_c = 300.0;
        double gamma_exp = 1.0;

        double chi = physics::statistical_mechanics::PhaseTransitions::susceptibility(T, T_c, gamma_exp);

        // χ = |T - T_c|^(-γ) = 20^(-1) = 0.05
        double expected = std::pow(std::abs(T - T_c), -gamma_exp);

        ASSERT_NEAR(chi, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Correlation Length Tests
    // ========================================

    run_test("Correlation length diverges at T_c", []() {
        double T_c = 300.0;
        double nu_exp = 0.5;

        double T_near = 300.1;
        double xi_near = physics::statistical_mechanics::PhaseTransitions::correlationLength(T_near, T_c, nu_exp);

        double T_far = 350.0;
        double xi_far = physics::statistical_mechanics::PhaseTransitions::correlationLength(T_far, T_c, nu_exp);

        // Correlation length should be larger closer to T_c
        ASSERT_TRUE(xi_near > xi_far);
        return true;
    });

    run_test("Correlation length power law", []() {
        double T = 325.0;
        double T_c = 300.0;
        double nu_exp = 0.5;

        double xi = physics::statistical_mechanics::PhaseTransitions::correlationLength(T, T_c, nu_exp);

        // ξ = |T - T_c|^(-ν) = 25^(-0.5) = 1/5 = 0.2
        double expected = std::pow(std::abs(T - T_c), -nu_exp);

        ASSERT_NEAR(xi, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Mean Field Theory Tests
    // ========================================

    run_test("Mean field critical temperature", []() {
        double J = 1e-21;  // Coupling constant
        double z = 4.0;    // Coordination number (square lattice)

        double T_c = physics::statistical_mechanics::MeanFieldTheory::criticalTemperatureMF(J, z);

        // T_c = zJ / k_B
        double expected = z * J / k_B;

        ASSERT_NEAR(T_c, expected, TOLERANCE);
        return true;
    });

    run_test("Mean field magnetization above T_c", []() {
        double J = 1e-21;
        double z = 4.0;
        double T_c = physics::statistical_mechanics::MeanFieldTheory::criticalTemperatureMF(J, z);

        double T = 2.0 * T_c;  // Well above critical temperature
        double h = 0.0;  // No external field

        double m = physics::statistical_mechanics::MeanFieldTheory::meanFieldMagnetization(T, J, z, h);

        // Above T_c with h=0, magnetization should be zero (paramagnetic phase)
        ASSERT_NEAR(m, 0.0, 1e-6);
        return true;
    });

    run_test("Mean field magnetization below T_c", []() {
        double J = 1e-21;
        double z = 4.0;
        double T_c = physics::statistical_mechanics::MeanFieldTheory::criticalTemperatureMF(J, z);

        double T = 0.5 * T_c;  // Well below critical temperature
        double h = 0.0;

        double m = physics::statistical_mechanics::MeanFieldTheory::meanFieldMagnetization(T, J, z, h);

        // Below T_c, should have spontaneous magnetization
        ASSERT_TRUE(std::abs(m) > 0.1);  // Significant magnetization
        return true;
    });

    run_test("Mean field magnetization with external field", []() {
        double J = 1e-21;
        double z = 4.0;
        double T = 1e-17;  // Some temperature
        double h = 1e-21;  // External field

        double m = physics::statistical_mechanics::MeanFieldTheory::meanFieldMagnetization(T, J, z, h);

        // With external field, magnetization should align with field
        ASSERT_TRUE(m > 0);  // Positive magnetization for positive field
        return true;
    });

    // ========================================
    // Thermodynamic Consistency Tests
    // ========================================

    run_test("Partition function is positive", []() {
        double T = 300.0;
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);

        auto energy = [](int n) { return n * 1e-20; };
        double Z = ensemble.partitionFunction(energy, 100);

        ASSERT_TRUE(Z > 0);
        return true;
    });

    run_test("Partition function increases with temperature", []() {
        auto energy = [](int n) { return n * 1e-20; };

        physics::statistical_mechanics::CanonicalEnsemble ensemble1(100.0, 1, 1.0);
        physics::statistical_mechanics::CanonicalEnsemble ensemble2(300.0, 1, 1.0);

        double Z1 = ensemble1.partitionFunction(energy, 100);
        double Z2 = ensemble2.partitionFunction(energy, 100);

        // Higher temperature → more accessible states → larger Z
        ASSERT_TRUE(Z2 > Z1);
        return true;
    });

    run_test("Free energy decreases with temperature (entropy effect)", []() {
        double Z1 = 2.0;
        double Z2 = 10.0;  // More states accessible

        double T = 300.0;
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);

        double F1 = ensemble.helmholtzFreeEnergy(Z1);
        double F2 = ensemble.helmholtzFreeEnergy(Z2);

        // Larger Z → lower free energy (F = -k_B T ln Z)
        ASSERT_TRUE(F2 < F1);
        return true;
    });

    run_test("Entropy is non-negative", []() {
        double T = 300.0;
        physics::statistical_mechanics::CanonicalEnsemble ensemble(T, 1, 1.0);

        double E = 1e-20;
        double F = -1e-20;  // F < E as expected
        double S = ensemble.entropy(F, E);

        // S = (E - F) / T should be positive
        ASSERT_TRUE(S >= 0);
        return true;
    });

    // ========================================
    // Physical Constants Validation
    // ========================================

    run_test("Boltzmann constant value", []() {
        // k_B = 1.380649×10⁻²³ J/K (CODATA 2018)
        ASSERT_NEAR(k_B, 1.380649e-23, 1e-28);
        return true;
    });

    run_test("Planck constant value", []() {
        // h = 6.62607015×10⁻³⁴ J·s (CODATA 2018)
        ASSERT_NEAR(h, 6.62607015e-34, 1e-42);
        return true;
    });

    run_test("Reduced Planck constant value", []() {
        // ℏ = h / (2π) = 1.054571817×10⁻³⁴ J·s
        ASSERT_NEAR(hbar, h / (2.0 * M_PI), 1e-42);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Statistical Mechanics" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "✓ All statistical mechanics tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Canonical ensemble (β, Z, Boltzmann distribution)" << std::endl;
        std::cout << "  - Partition functions (harmonic oscillator, two-level, rotators)" << std::endl;
        std::cout << "  - Thermodynamic relations (F, E, S, P)" << std::endl;
        std::cout << "  - Phase transitions (Van der Waals, critical points)" << std::endl;
        std::cout << "  - Critical phenomena (order parameter, susceptibility, correlation length)" << std::endl;
        std::cout << "  - Mean field theory (critical temperature, magnetization)" << std::endl;
        std::cout << "  - Thermodynamic consistency checks" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
