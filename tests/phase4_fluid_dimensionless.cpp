/**
 * Phase 4 Validation: Fluid Dynamics Dimensionless Numbers
 *
 * Tests the fluid_dynamics_dimensionless_numbers.hpp module functions.
 *
 * Coverage:
 * - Reynolds number (inertia/viscous)
 * - Froude number (inertia/gravity)
 * - Mach number (velocity/sound speed)
 * - Prandtl number (momentum/thermal diffusivity)
 * - Grashof and Rayleigh numbers (natural convection)
 * - Nusselt number (heat transfer)
 * - Peclet number (advection/diffusion)
 * - Weber and Capillary numbers (surface tension)
 * - Flow regime classifications
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/fluid_dynamics_dimensionless_numbers.hpp"

// Test tolerances
const double TOLERANCE = 1e-6;
const double NUMERICAL_TOLERANCE = 1e-3;
const double LOOSE_TOLERANCE = 0.01;

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

using namespace physics::advanced::fluid_dynamics;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Fluid Dynamics Dimensionless Numbers Validation ===" << std::endl;
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
    // Reynolds Number Tests
    // ========================================

    run_test("Reynolds: laminar pipe flow", []() {
        double rho = 1000.0;  // Water (kg/m³)
        double U = 0.1;       // m/s
        double L = 0.05;      // 5 cm diameter
        double mu = 1e-3;     // Water viscosity (Pa·s)

        double Re = DimensionlessNumbers::reynoldsNumber(rho, U, L, mu);
        double expected = 1000.0 * 0.1 * 0.05 / 1e-3;  // 5000
        ASSERT_NEAR(Re, expected, TOLERANCE);
        return true;
    });

    run_test("Reynolds: turbulent pipe flow", []() {
        double rho = 1000.0;
        double U = 1.0;
        double L = 0.1;
        double mu = 1e-3;

        double Re = DimensionlessNumbers::reynoldsNumber(rho, U, L, mu);
        ASSERT_TRUE(Re > 4000.0);  // Turbulent regime
        return true;
    });

    run_test("Reynolds: kinematic viscosity form", []() {
        double U = 10.0;       // m/s
        double L = 1.0;        // m
        double nu = 1.5e-5;    // Air kinematic viscosity (m²/s)

        double Re = DimensionlessNumbers::reynoldsNumberKinematic(U, L, nu);
        double expected = 10.0 * 1.0 / 1.5e-5;  // ≈ 666667
        ASSERT_NEAR(Re, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Reynolds: flow around sphere", []() {
        // Golf ball: Re ≈ 100000
        double rho = 1.225;    // Air
        double U = 50.0;       // m/s
        double D = 0.043;      // 43 mm diameter
        double mu = 1.81e-5;

        double Re = DimensionlessNumbers::reynoldsNumber(rho, U, D, mu);
        ASSERT_TRUE(Re > 50000.0 && Re < 150000.0);
        return true;
    });

    run_test("Reynolds: Stokes flow regime", []() {
        double Re = 0.5;
        auto regime = DimensionlessNumbers::classifyFlowRegime(Re);
        ASSERT_TRUE(regime == DimensionlessNumbers::FlowRegime::STOKES);
        return true;
    });

    run_test("Reynolds: laminar regime classification", []() {
        double Re = 1500.0;
        auto regime = DimensionlessNumbers::classifyFlowRegime(Re);
        ASSERT_TRUE(regime == DimensionlessNumbers::FlowRegime::LAMINAR);
        return true;
    });

    run_test("Reynolds: transition regime classification", []() {
        double Re = 3000.0;
        auto regime = DimensionlessNumbers::classifyFlowRegime(Re);
        ASSERT_TRUE(regime == DimensionlessNumbers::FlowRegime::TRANSITION);
        return true;
    });

    run_test("Reynolds: turbulent regime classification", []() {
        double Re = 10000.0;
        auto regime = DimensionlessNumbers::classifyFlowRegime(Re);
        ASSERT_TRUE(regime == DimensionlessNumbers::FlowRegime::TURBULENT);
        return true;
    });

    // ========================================
    // Froude Number Tests
    // ========================================

    run_test("Froude: subcritical open channel flow", []() {
        double U = 1.0;     // m/s
        double g = 9.81;    // m/s²
        double h = 0.5;     // 0.5 m depth

        double Fr = DimensionlessNumbers::froudeNumber(U, g, h);
        double expected = 1.0 / std::sqrt(9.81 * 0.5);  // ≈ 0.451
        ASSERT_NEAR(Fr, expected, TOLERANCE);
        ASSERT_TRUE(Fr < 1.0);  // Subcritical
        return true;
    });

    run_test("Froude: critical flow", []() {
        double U = std::sqrt(9.81 * 1.0);  // √(gh)
        double g = 9.81;
        double h = 1.0;

        double Fr = DimensionlessNumbers::froudeNumber(U, g, h);
        ASSERT_NEAR(Fr, 1.0, TOLERANCE);
        return true;
    });

    run_test("Froude: supercritical flow", []() {
        double U = 5.0;
        double g = 9.81;
        double h = 0.2;

        double Fr = DimensionlessNumbers::froudeNumber(U, g, h);
        double expected = 5.0 / std::sqrt(9.81 * 0.2);
        ASSERT_TRUE(Fr > 1.0);  // Supercritical
        return true;
    });

    run_test("Froude: ship waves", []() {
        // Ship speed 10 knots ≈ 5.14 m/s, length 100 m
        double U = 5.14;
        double g = 9.81;
        double L = 100.0;

        double Fr = DimensionlessNumbers::froudeNumber(U, g, L);
        ASSERT_TRUE(Fr > 0.0 && Fr < 1.0);  // Typical ship Fr
        return true;
    });

    // ========================================
    // Mach Number Tests
    // ========================================

    run_test("Mach: sound speed in air", []() {
        double gamma = 1.4;     // Air
        double R = 287.0;       // J/(kg·K) for air
        double T = 288.15;      // 15°C

        double c = DimensionlessNumbers::soundSpeed(gamma, R, T);
        double expected = std::sqrt(1.4 * 287.0 * 288.15);  // ≈ 340 m/s
        ASSERT_NEAR(c, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Mach: incompressible flow", []() {
        double U = 50.0;       // m/s
        double c = 340.0;      // m/s

        double Ma = DimensionlessNumbers::machNumber(U, c);
        double expected = 50.0 / 340.0;  // ≈ 0.147
        ASSERT_NEAR(Ma, expected, TOLERANCE);

        auto regime = DimensionlessNumbers::classifyCompressibility(Ma);
        ASSERT_TRUE(regime == DimensionlessNumbers::CompressibilityRegime::INCOMPRESSIBLE);
        return true;
    });

    run_test("Mach: subsonic flow", []() {
        double Ma = 0.6;
        auto regime = DimensionlessNumbers::classifyCompressibility(Ma);
        ASSERT_TRUE(regime == DimensionlessNumbers::CompressibilityRegime::SUBSONIC);
        return true;
    });

    run_test("Mach: transonic flow", []() {
        double Ma = 0.95;
        auto regime = DimensionlessNumbers::classifyCompressibility(Ma);
        ASSERT_TRUE(regime == DimensionlessNumbers::CompressibilityRegime::TRANSONIC);
        return true;
    });

    run_test("Mach: supersonic flow", []() {
        double Ma = 2.5;
        auto regime = DimensionlessNumbers::classifyCompressibility(Ma);
        ASSERT_TRUE(regime == DimensionlessNumbers::CompressibilityRegime::SUPERSONIC);
        return true;
    });

    run_test("Mach: hypersonic flow", []() {
        double Ma = 8.0;
        auto regime = DimensionlessNumbers::classifyCompressibility(Ma);
        ASSERT_TRUE(regime == DimensionlessNumbers::CompressibilityRegime::HYPERSONIC);
        return true;
    });

    // ========================================
    // Prandtl Number Tests
    // ========================================

    run_test("Prandtl: air at standard conditions", []() {
        double nu = 1.5e-5;     // m²/s
        double alpha = 2.2e-5;  // m²/s thermal diffusivity

        double Pr = DimensionlessNumbers::prandtlNumber(nu, alpha);
        double expected = 1.5e-5 / 2.2e-5;  // ≈ 0.68
        ASSERT_NEAR(Pr, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Prandtl: water", []() {
        double nu = 1e-6;       // m²/s
        double alpha = 1.4e-7;  // m²/s

        double Pr = DimensionlessNumbers::prandtlNumber(nu, alpha);
        ASSERT_TRUE(Pr > 5.0 && Pr < 10.0);  // Water Pr ≈ 7
        return true;
    });

    run_test("Prandtl: from physical properties", []() {
        double mu = 1.81e-5;     // Pa·s (air)
        double cp = 1005.0;      // J/(kg·K)
        double k = 0.026;        // W/(m·K)

        double Pr = DimensionlessNumbers::prandtlNumberFromProperties(mu, cp, k);
        double expected = (1.81e-5 * 1005.0) / 0.026;  // ≈ 0.70
        ASSERT_NEAR(Pr, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Prandtl: liquid metal low value", []() {
        // Liquid metals have Pr << 1
        double nu = 1e-7;
        double alpha = 1e-5;

        double Pr = DimensionlessNumbers::prandtlNumber(nu, alpha);
        ASSERT_TRUE(Pr < 0.1);  // Very low Pr
        return true;
    });

    run_test("Prandtl: oil high value", []() {
        // Oils have Pr >> 1
        double nu = 1e-4;
        double alpha = 1e-7;

        double Pr = DimensionlessNumbers::prandtlNumber(nu, alpha);
        ASSERT_TRUE(Pr > 100.0);  // Very high Pr
        return true;
    });

    // ========================================
    // Grashof and Rayleigh Number Tests
    // ========================================

    run_test("Grashof: natural convection laminar", []() {
        double g = 9.81;
        double beta = 1.0/300.0;  // Thermal expansion (1/K)
        double dT = 10.0;         // 10 K temperature difference
        double L = 0.1;           // 0.1 m
        double nu = 1.5e-5;       // Air

        double Gr = DimensionlessNumbers::grashofNumber(g, beta, dT, L, nu);
        ASSERT_TRUE(Gr > 0.0);
        return true;
    });

    run_test("Grashof: larger temperature difference", []() {
        double g = 9.81;
        double beta = 1.0/300.0;
        double dT1 = 10.0;
        double dT2 = 20.0;
        double L = 0.1;
        double nu = 1.5e-5;

        double Gr1 = DimensionlessNumbers::grashofNumber(g, beta, dT1, L, nu);
        double Gr2 = DimensionlessNumbers::grashofNumber(g, beta, dT2, L, nu);

        ASSERT_NEAR(Gr2 / Gr1, 2.0, TOLERANCE);  // Gr ∝ ΔT
        return true;
    });

    run_test("Rayleigh: from Grashof and Prandtl", []() {
        double Gr = 1e8;
        double Pr = 0.7;

        double Ra = DimensionlessNumbers::rayleighNumber(Gr, Pr);
        double expected = 1e8 * 0.7;
        ASSERT_NEAR(Ra, expected, TOLERANCE);
        return true;
    });

    run_test("Rayleigh: from physical properties", []() {
        double g = 9.81;
        double beta = 1.0/300.0;
        double dT = 20.0;
        double L = 0.5;
        double nu = 1.5e-5;
        double alpha = 2.2e-5;

        double Ra = DimensionlessNumbers::rayleighNumberFromProperties(
            g, beta, dT, L, nu, alpha);

        double Gr = DimensionlessNumbers::grashofNumber(g, beta, dT, L, nu);
        double Pr = DimensionlessNumbers::prandtlNumber(nu, alpha);
        double expected = Gr * Pr;

        ASSERT_NEAR(Ra, expected, Ra * LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Nusselt Number Tests
    // ========================================

    run_test("Nusselt: pure conduction", []() {
        double h = 25.0;    // W/(m²·K)
        double L = 1.0;     // m
        double k = 25.0;    // W/(m·K)

        double Nu = DimensionlessNumbers::nusseltNumber(h, L, k);
        ASSERT_NEAR(Nu, 1.0, TOLERANCE);  // Nu = 1 for pure conduction
        return true;
    });

    run_test("Nusselt: laminar flat plate correlation", []() {
        double Re = 10000.0;
        double Pr = 0.7;

        double Nu = DimensionlessNumbers::nusseltLaminarFlatPlate(Re, Pr);
        double expected = 0.664 * std::pow(10000.0, 0.5) * std::pow(0.7, 1.0/3.0);
        ASSERT_NEAR(Nu, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Nusselt: turbulent flat plate correlation", []() {
        double Re = 1e6;
        double Pr = 0.7;

        double Nu = DimensionlessNumbers::nusseltTurbulentFlatPlate(Re, Pr);
        double expected = 0.037 * std::pow(1e6, 0.8) * std::pow(0.7, 1.0/3.0);
        ASSERT_NEAR(Nu, expected, expected * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Nusselt: higher Reynolds increases heat transfer", []() {
        double Re1 = 1000.0;
        double Re2 = 10000.0;
        double Pr = 0.7;

        double Nu1 = DimensionlessNumbers::nusseltLaminarFlatPlate(Re1, Pr);
        double Nu2 = DimensionlessNumbers::nusseltLaminarFlatPlate(Re2, Pr);

        ASSERT_TRUE(Nu2 > Nu1);  // Higher Re → higher Nu
        return true;
    });

    // ========================================
    // Peclet Number Tests
    // ========================================

    run_test("Peclet: advection dominated", []() {
        double U = 1.0;
        double L = 1.0;
        double alpha = 1e-5;

        double Pe = DimensionlessNumbers::pecletNumber(U, L, alpha);
        ASSERT_TRUE(Pe > 100.0);  // Advection dominates
        return true;
    });

    run_test("Peclet: from Reynolds and Prandtl", []() {
        double Re = 1000.0;
        double Pr = 0.7;

        double Pe = DimensionlessNumbers::pecletFromReynoldsPrandtl(Re, Pr);
        ASSERT_NEAR(Pe, 700.0, TOLERANCE);
        return true;
    });

    run_test("Peclet: thermal vs mass transfer", []() {
        double U = 0.1;
        double L = 0.5;
        double D_thermal = 2e-5;
        double D_mass = 1e-9;

        double Pe_thermal = DimensionlessNumbers::pecletNumber(U, L, D_thermal);
        double Pe_mass = DimensionlessNumbers::pecletNumber(U, L, D_mass);

        ASSERT_TRUE(Pe_mass > Pe_thermal);  // Mass diffusion slower
        return true;
    });

    // ========================================
    // Schmidt and Sherwood Number Tests
    // ========================================

    run_test("Schmidt: mass transfer analogy to Prandtl", []() {
        double nu = 1.5e-5;     // m²/s
        double D = 2e-5;        // m²/s mass diffusivity

        double Sc = DimensionlessNumbers::schmidtNumber(nu, D);
        double expected = 1.5e-5 / 2e-5;  // 0.75
        ASSERT_NEAR(Sc, expected, TOLERANCE);
        return true;
    });

    run_test("Sherwood: mass transfer analogy to Nusselt", []() {
        double k_m = 0.01;     // m/s mass transfer coefficient
        double L = 1.0;        // m
        double D = 1e-9;       // m²/s

        double Sh = DimensionlessNumbers::sherwoodNumber(k_m, L, D);
        double expected = 0.01 * 1.0 / 1e-9;  // 1e7
        ASSERT_NEAR(Sh, expected, Sh * LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Weber and Capillary Number Tests
    // ========================================

    run_test("Weber: surface tension dominant", []() {
        double rho = 1000.0;   // Water
        double U = 0.1;        // m/s
        double L = 0.001;      // 1 mm droplet
        double sigma = 0.072;  // N/m

        double We = DimensionlessNumbers::weberNumber(rho, U, L, sigma);
        double expected = 1000.0 * 0.1 * 0.1 * 0.001 / 0.072;
        ASSERT_NEAR(We, expected, TOLERANCE);
        ASSERT_TRUE(We < 1.0);  // Surface tension dominant
        return true;
    });

    run_test("Weber: inertia dominant droplet breakup", []() {
        double rho = 1000.0;
        double U = 10.0;       // High velocity
        double L = 0.001;
        double sigma = 0.072;

        double We = DimensionlessNumbers::weberNumber(rho, U, L, sigma);
        ASSERT_TRUE(We > 10.0);  // Inertia dominant
        return true;
    });

    run_test("Capillary: viscous surface tension balance", []() {
        double mu = 1e-3;      // Water
        double U = 0.5;        // m/s
        double sigma = 0.072;  // N/m

        double Ca = DimensionlessNumbers::capillaryNumber(mu, U, sigma);
        double expected = 1e-3 * 0.5 / 0.072;
        ASSERT_NEAR(Ca, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Strouhal Number Tests
    // ========================================

    run_test("Strouhal: vortex shedding frequency", []() {
        double f = 2.0;        // Hz
        double L = 0.1;        // m
        double U = 1.0;        // m/s

        double St = DimensionlessNumbers::strouhalNumber(f, L, U);
        double expected = 2.0 * 0.1 / 1.0;  // 0.2
        ASSERT_NEAR(St, expected, TOLERANCE);
        return true;
    });

    run_test("Strouhal: cylinder vortex shedding", []() {
        // Typical St ≈ 0.2 for cylinder
        double f = 10.0;       // Hz
        double D = 0.05;       // 5 cm diameter
        double U = 2.5;        // m/s

        double St = DimensionlessNumbers::strouhalNumber(f, D, U);
        ASSERT_NEAR(St, 0.2, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Consistency: Reynolds from kinematic equals dynamic form", []() {
        double rho = 1000.0;
        double U = 1.0;
        double L = 0.1;
        double mu = 1e-3;
        double nu = mu / rho;

        double Re1 = DimensionlessNumbers::reynoldsNumber(rho, U, L, mu);
        double Re2 = DimensionlessNumbers::reynoldsNumberKinematic(U, L, nu);

        ASSERT_NEAR(Re1, Re2, TOLERANCE);
        return true;
    });

    run_test("Consistency: Peclet equals Reynolds times Prandtl", []() {
        double U = 1.0;
        double L = 1.0;
        double nu = 1.5e-5;
        double alpha = 2.2e-5;

        double Re = DimensionlessNumbers::reynoldsNumberKinematic(U, L, nu);
        double Pr = DimensionlessNumbers::prandtlNumber(nu, alpha);
        double Pe1 = DimensionlessNumbers::pecletFromReynoldsPrandtl(Re, Pr);
        double Pe2 = DimensionlessNumbers::pecletNumber(U, L, alpha);

        ASSERT_NEAR(Pe1, Pe2, Pe1 * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Consistency: Rayleigh equals Grashof times Prandtl", []() {
        double g = 9.81;
        double beta = 1.0/300.0;
        double dT = 15.0;
        double L = 0.3;
        double nu = 1.5e-5;
        double alpha = 2.2e-5;

        double Gr = DimensionlessNumbers::grashofNumber(g, beta, dT, L, nu);
        double Pr = DimensionlessNumbers::prandtlNumber(nu, alpha);
        double Ra1 = DimensionlessNumbers::rayleighNumber(Gr, Pr);
        double Ra2 = DimensionlessNumbers::rayleighNumberFromProperties(
            g, beta, dT, L, nu, alpha);

        ASSERT_NEAR(Ra1, Ra2, Ra1 * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Physical: Reynolds increases with velocity", []() {
        double rho = 1.225;
        double L = 1.0;
        double mu = 1.81e-5;

        double Re1 = DimensionlessNumbers::reynoldsNumber(rho, 1.0, L, mu);
        double Re2 = DimensionlessNumbers::reynoldsNumber(rho, 10.0, L, mu);

        ASSERT_NEAR(Re2 / Re1, 10.0, TOLERANCE);
        return true;
    });

    run_test("Physical: all dimensionless numbers are positive", []() {
        // Test that physical values give positive dimensionless numbers
        double Re = DimensionlessNumbers::reynoldsNumber(1.2, 10.0, 1.0, 1.8e-5);
        double Fr = DimensionlessNumbers::froudeNumber(5.0, 9.81, 1.0);
        double Ma = DimensionlessNumbers::machNumber(100.0, 340.0);
        double Pr = DimensionlessNumbers::prandtlNumber(1.5e-5, 2.2e-5);

        ASSERT_TRUE(Re > 0.0);
        ASSERT_TRUE(Fr > 0.0);
        ASSERT_TRUE(Ma > 0.0);
        ASSERT_TRUE(Pr > 0.0);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Fluid Dimensionless Numbers" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All fluid dynamics dimensionless numbers tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Reynolds number and flow regimes" << std::endl;
        std::cout << "  - Froude number and open channel flow" << std::endl;
        std::cout << "  - Mach number and compressibility" << std::endl;
        std::cout << "  - Prandtl number and heat transfer" << std::endl;
        std::cout << "  - Grashof and Rayleigh numbers" << std::endl;
        std::cout << "  - Nusselt number correlations" << std::endl;
        std::cout << "  - Peclet and Schmidt numbers" << std::endl;
        std::cout << "  - Weber and Capillary numbers" << std::endl;
        std::cout << "  - Strouhal number and vortex shedding" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
