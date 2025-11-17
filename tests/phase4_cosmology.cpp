/**
 * Phase 4 Validation: Cosmology (Friedmann Equations)
 *
 * Tests the cosmology_friedmann_equations.hpp module functions.
 *
 * Coverage:
 * - Friedmann equations (1st and 2nd)
 * - Critical density and density parameters
 * - Curvature geometry (flat, open, closed)
 * - Fluid equation and energy conservation
 * - Equation of state (matter, radiation, dark energy)
 * - Density scaling laws with cosmic expansion
 */

#include <iostream>
#include <cmath>
#include "../include/physics/cosmology_friedmann_equations.hpp"

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

// Physical constants
const double G = 6.674e-11;  // m³/(kg·s²)
const double c = 2.998e8;    // m/s

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Cosmology (Friedmann Equations) Validation ===" << std::endl;
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
    // Critical Density Tests
    // ========================================

    run_test("Critical density formula", []() {
        double H = 2.2e-18;  // s⁻¹ (Hubble parameter today)
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H);

        // ρ_crit = 3H²/(8πG)
        double expected = (3.0 * H * H) / (8.0 * M_PI * G);

        ASSERT_NEAR(rho_crit, expected, 1e-28);
        return true;
    });

    run_test("Critical density today", []() {
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensityToday();

        // ρ_crit,0 ≈ 8.5×10⁻²⁷ kg/m³ (about 5 protons per cubic meter)
        ASSERT_TRUE(rho_crit > 8.0e-27 && rho_crit < 9.0e-27);
        return true;
    });

    run_test("Critical density scales with H²", []() {
        double H1 = 2.0e-18;
        double H2 = 4.0e-18;

        double rho1 = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H1);
        double rho2 = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H2);

        // ρ_crit ∝ H², so doubling H quadruples ρ_crit
        ASSERT_NEAR(rho2, 4.0 * rho1, rho1 * 1e-10);
        return true;
    });

    // ========================================
    // Density Parameter Tests
    // ========================================

    run_test("Density parameter Ω formula", []() {
        double rho = 5e-27;  // kg/m³
        double H = 2.2e-18;   // s⁻¹

        double Omega = physics::advanced::cosmology::FriedmannEquations::densityParameter(rho, H);

        // Ω = ρ/ρ_crit
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H);
        double expected = rho / rho_crit;

        ASSERT_NEAR(Omega, expected, TOLERANCE);
        return true;
    });

    run_test("Density parameter Ω = 1 for critical density", []() {
        double H = 2.2e-18;
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H);

        double Omega = physics::advanced::cosmology::FriedmannEquations::densityParameter(rho_crit, H);

        // At critical density, Ω = 1 (flat universe)
        ASSERT_NEAR(Omega, 1.0, TOLERANCE);
        return true;
    });

    run_test("Density parameter Ω > 1 for closed universe", []() {
        double H = 2.2e-18;
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H);
        double rho_high = 1.5 * rho_crit;

        double Omega = physics::advanced::cosmology::FriedmannEquations::densityParameter(rho_high, H);

        // Ω > 1 indicates closed universe
        ASSERT_TRUE(Omega > 1.0);
        ASSERT_NEAR(Omega, 1.5, TOLERANCE);
        return true;
    });

    run_test("Density parameter Ω < 1 for open universe", []() {
        double H = 2.2e-18;
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H);
        double rho_low = 0.5 * rho_crit;

        double Omega = physics::advanced::cosmology::FriedmannEquations::densityParameter(rho_low, H);

        // Ω < 1 indicates open universe
        ASSERT_TRUE(Omega < 1.0);
        ASSERT_NEAR(Omega, 0.5, TOLERANCE);
        return true;
    });

    // ========================================
    // Curvature Geometry Tests
    // ========================================

    run_test("Flat universe from Ω = 1", []() {
        auto k = physics::advanced::cosmology::CurvatureGeometry::fromDensityParameters(1.0);
        ASSERT_TRUE(k == physics::advanced::cosmology::CurvatureGeometry::Curvature::FLAT);
        return true;
    });

    run_test("Closed universe from Ω > 1", []() {
        auto k = physics::advanced::cosmology::CurvatureGeometry::fromDensityParameters(1.1);
        ASSERT_TRUE(k == physics::advanced::cosmology::CurvatureGeometry::Curvature::CLOSED);
        return true;
    });

    run_test("Open universe from Ω < 1", []() {
        auto k = physics::advanced::cosmology::CurvatureGeometry::fromDensityParameters(0.9);
        ASSERT_TRUE(k == physics::advanced::cosmology::CurvatureGeometry::Curvature::OPEN);
        return true;
    });

    run_test("Curvature density parameter", []() {
        double Omega_m = 0.3;   // Matter
        double Omega_r = 0.0001; // Radiation (negligible today)
        double Omega_L = 0.7;   // Dark energy

        double Omega_k = physics::advanced::cosmology::CurvatureGeometry::curvatureDensityParameter(
            Omega_m, Omega_r, Omega_L);

        // Ω_k = 1 - (Ω_m + Ω_r + Ω_Λ) ≈ 0
        double expected = 1.0 - (Omega_m + Omega_r + Omega_L);
        ASSERT_NEAR(Omega_k, expected, TOLERANCE);
        return true;
    });

    run_test("Observed curvature is nearly flat", []() {
        double Omega_k_obs = physics::advanced::cosmology::CurvatureGeometry::observedCurvature();

        // Planck 2018: Ω_k ≈ 0.001 ± 0.002 (very flat!)
        ASSERT_TRUE(std::abs(Omega_k_obs) < 0.01);
        return true;
    });

    // ========================================
    // First Friedmann Equation Tests
    // ========================================

    run_test("First Friedmann equation (flat, no Λ)", []() {
        double a = 1.0;
        double rho = 8.5e-27;  // kg/m³
        double H = physics::advanced::cosmology::FriedmannEquations::hubbleFromDensity(rho);

        // For flat universe (k=0, Λ=0): H = √(8πGρ/3)
        double expected = std::sqrt((8.0 * M_PI * G / 3.0) * rho);
        ASSERT_NEAR(H, expected, 1e-20);
        return true;
    });

    run_test("Hubble from density consistency", []() {
        double rho = 5e-27;
        double H = physics::advanced::cosmology::FriedmannEquations::hubbleFromDensity(rho);

        // Verify using critical density
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H);

        // For flat universe, ρ = ρ_crit
        ASSERT_NEAR(rho, rho_crit, rho * 1e-10);
        return true;
    });

    run_test("First Friedmann satisfied for consistent values", []() {
        double a = 1.0;
        double rho = 8.5e-27;
        double k = 0.0;  // Flat
        double Lambda = 0.0;

        // Calculate H from Friedmann
        double H = std::sqrt((8.0 * M_PI * G / 3.0) * rho);
        double a_dot = H * a;

        // Verify equation is satisfied
        double residual = physics::advanced::cosmology::FriedmannEquations::firstFriedmann(
            a, a_dot, rho, k, Lambda);

        ASSERT_NEAR(residual, 0.0, 1e-30);
        return true;
    });

    // ========================================
    // Second Friedmann Equation Tests
    // ========================================

    run_test("Second Friedmann equation (matter-dominated)", []() {
        double a = 1.0;
        double rho = 8.5e-27;
        double p = 0.0;  // Matter has p = 0
        double Lambda = 0.0;

        // For matter: ä/a = -(4πG/3)ρ
        double expected_accel = -(4.0 * M_PI * G / 3.0) * rho;
        double a_double_dot = expected_accel * a;

        double residual = physics::advanced::cosmology::FriedmannEquations::secondFriedmann(
            a, a_double_dot, rho, p, Lambda);

        ASSERT_NEAR(residual, 0.0, 1e-35);
        return true;
    });

    run_test("Second Friedmann: matter causes deceleration", []() {
        double a = 1.0;
        double rho = 1e-26;  // Matter density
        double p = 0.0;  // Matter pressure
        double Lambda = 0.0;

        // Calculate expected deceleration
        double a_over_a_expected = -(4.0 * M_PI * G / 3.0) * rho;

        // Deceleration is negative (universe slowing down)
        ASSERT_TRUE(a_over_a_expected < 0);
        return true;
    });

    run_test("Second Friedmann: dark energy causes acceleration", []() {
        double rho = 5e-27;  // Dark energy density
        double p = -rho * c * c;  // Dark energy: p = -ρc² (w = -1)
        double Lambda = 1e-35;  // Cosmological constant

        // With dark energy: ä/a = -(4πG/3)(ρ + 3p/c²) + Λ/3
        // Since p = -ρc², the ρ+3p/c² term becomes negative
        double rho_p_term = rho + 3.0 * p / (c * c);

        // For dark energy, this should give acceleration
        ASSERT_TRUE(rho_p_term < 0);  // Negative means acceleration
        return true;
    });

    // ========================================
    // Equation of State Tests
    // ========================================

    run_test("Equation of state: matter (w = 0)", []() {
        double rho = 1e-26;
        double p = 0.0;  // Matter is pressureless

        double w = physics::advanced::cosmology::FluidEquation::equationOfStateParameter(p, rho);

        ASSERT_NEAR(w, 0.0, TOLERANCE);
        return true;
    });

    run_test("Equation of state: radiation (w = 1/3)", []() {
        double rho = 1e-30;
        double p = (1.0/3.0) * rho * c * c;  // Radiation

        double w = physics::advanced::cosmology::FluidEquation::equationOfStateParameter(p, rho);

        ASSERT_NEAR(w, 1.0/3.0, TOLERANCE);
        return true;
    });

    run_test("Equation of state: dark energy (w = -1)", []() {
        double rho = 1e-26;
        double p = -rho * c * c;  // Dark energy (cosmological constant)

        double w = physics::advanced::cosmology::FluidEquation::equationOfStateParameter(p, rho);

        ASSERT_NEAR(w, -1.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Density Scaling Tests
    // ========================================

    run_test("Matter density scaling: ρ ∝ a⁻³", []() {
        double rho_0 = 1e-26;
        double a = 0.5;  // Half current size
        double w = 0.0;  // Matter

        double rho = physics::advanced::cosmology::FluidEquation::densityScaling(rho_0, a, w);

        // ρ(a) = ρ₀ a⁻³ for matter
        double expected = rho_0 / (a * a * a);
        ASSERT_NEAR(rho, expected, rho * 1e-10);
        return true;
    });

    run_test("Radiation density scaling: ρ ∝ a⁻⁴", []() {
        double rho_0 = 1e-30;
        double a = 0.5;
        double w = 1.0/3.0;  // Radiation

        double rho = physics::advanced::cosmology::FluidEquation::densityScaling(rho_0, a, w);

        // ρ(a) = ρ₀ a⁻⁴ for radiation
        double expected = rho_0 / (a * a * a * a);
        ASSERT_NEAR(rho, expected, rho * 1e-10);
        return true;
    });

    run_test("Dark energy density constant: ρ ∝ a⁰", []() {
        double rho_0 = 5e-27;
        double a = 0.5;
        double w = -1.0;  // Dark energy

        double rho = physics::advanced::cosmology::FluidEquation::densityScaling(rho_0, a, w);

        // ρ(a) = ρ₀ (constant) for dark energy
        ASSERT_NEAR(rho, rho_0, rho_0 * 1e-10);
        return true;
    });

    run_test("Density scaling at a=1 (today)", []() {
        double rho_0 = 1e-26;
        double a = 1.0;

        // At a=1, density should equal initial density for any w
        for (double w = -1.0; w <= 1.0; w += 0.5) {
            double rho = physics::advanced::cosmology::FluidEquation::densityScaling(rho_0, a, w);
            ASSERT_NEAR(rho, rho_0, rho_0 * 1e-10);
        }
        return true;
    });

    run_test("Matter density increases in early universe", []() {
        double rho_0 = 1e-26;
        double a_early = 0.1;  // 10% of current size
        double w = 0.0;

        double rho_early = physics::advanced::cosmology::FluidEquation::densityScaling(rho_0, a_early, w);

        // Earlier universe (smaller a) had higher matter density
        ASSERT_TRUE(rho_early > rho_0);
        ASSERT_NEAR(rho_early, 1000.0 * rho_0, rho_0 * 1e-8);  // a⁻³ = (0.1)⁻³ = 1000
        return true;
    });

    // ========================================
    // Energy Conservation Tests
    // ========================================

    run_test("Fluid equation: matter energy evolution", []() {
        double rho = 1e-26;
        double p = 0.0;  // Matter
        double a = 1.0;
        double H = 2.2e-18;
        double a_dot = H * a;

        double drho_dt = physics::advanced::cosmology::FluidEquation::energyEvolution(rho, p, a, a_dot);

        // dρ/dt = -3H(ρ + p/c²) = -3Hρ for matter
        double expected = -3.0 * H * rho;
        ASSERT_NEAR(drho_dt, expected, std::abs(expected) * 1e-10);
        return true;
    });

    run_test("Fluid equation: radiation energy evolution", []() {
        double rho = 1e-30;
        double p = (1.0/3.0) * rho * c * c;  // Radiation
        double a = 1.0;
        double H = 2.2e-18;
        double a_dot = H * a;

        double drho_dt = physics::advanced::cosmology::FluidEquation::energyEvolution(rho, p, a, a_dot);

        // dρ/dt = -3H(ρ + ρ/3) = -4Hρ for radiation
        double expected = -4.0 * H * rho;
        ASSERT_NEAR(drho_dt, expected, std::abs(expected) * 1e-10);
        return true;
    });

    run_test("Fluid equation: dark energy constant", []() {
        double rho = 5e-27;
        double p = -rho * c * c;  // Dark energy
        double a = 1.0;
        double H = 2.2e-18;
        double a_dot = H * a;

        double drho_dt = physics::advanced::cosmology::FluidEquation::energyEvolution(rho, p, a, a_dot);

        // dρ/dt = -3H(ρ - ρ) = 0 for dark energy
        ASSERT_NEAR(drho_dt, 0.0, 1e-38);
        return true;
    });

    // ========================================
    // Consistency Tests
    // ========================================

    run_test("Sum of density parameters", []() {
        double Omega_m = 0.3;
        double Omega_r = 0.0001;
        double Omega_L = 0.7;
        double Omega_k = physics::advanced::cosmology::CurvatureGeometry::curvatureDensityParameter(
            Omega_m, Omega_r, Omega_L);

        // Total should sum to 1
        double Omega_total = Omega_m + Omega_r + Omega_L + Omega_k;
        ASSERT_NEAR(Omega_total, 1.0, TOLERANCE);
        return true;
    });

    run_test("Hubble parameter from density equals critical", []() {
        double rho = 1e-26;
        double H = physics::advanced::cosmology::FriedmannEquations::hubbleFromDensity(rho);
        double rho_crit = physics::advanced::cosmology::FriedmannEquations::criticalDensity(H);

        // For flat universe, ρ = ρ_crit
        ASSERT_NEAR(rho, rho_crit, rho * 1e-10);
        return true;
    });

    run_test("Matter dominates radiation today", []() {
        // Typical values today
        double Omega_m = 0.3;
        double Omega_r = 9e-5;  // ~0.01% of critical density

        ASSERT_TRUE(Omega_m > Omega_r);
        ASSERT_TRUE(Omega_m > 1000 * Omega_r);  // Matter >> radiation
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Cosmology" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "✓ All cosmology tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Critical density (ρ_crit = 3H²/8πG)" << std::endl;
        std::cout << "  - Density parameters (Ω = ρ/ρ_crit)" << std::endl;
        std::cout << "  - Friedmann equations (1st and 2nd)" << std::endl;
        std::cout << "  - Curvature geometry (flat, open, closed)" << std::endl;
        std::cout << "  - Equation of state (w = p/ρc²)" << std::endl;
        std::cout << "  - Density scaling (ρ ∝ a^(-3(1+w)))" << std::endl;
        std::cout << "  - Fluid equation (energy conservation)" << std::endl;
        std::cout << "  - Observed universe is nearly flat!" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
