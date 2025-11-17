/**
 * Phase 4 Validation: Fluid Dynamics Compressible Flow
 *
 * Tests the fluid_dynamics_compressible_flow.hpp module functions.
 *
 * Coverage:
 * - Isentropic flow relations
 * - Stagnation properties
 * - Critical flow conditions
 * - Normal shock relations (Rankine-Hugoniot)
 * - Oblique shock relations
 * - Prandtl-Meyer expansion
 * - Mach number calculations
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/fluid_dynamics_compressible_flow.hpp"

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

    std::cout << "=== Phase 4: Compressible Flow Validation ===" << std::endl;
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

    // Physical constants
    constexpr double GAMMA_AIR = 1.4;

    // ========================================
    // Isentropic Flow Tests
    // ========================================

    run_test("Isentropic: temperature ratio from pressure ratio", []() {
        double p_ratio = 2.0;
        double T_ratio = IsentropicFlow::temperatureRatio(p_ratio, GAMMA_AIR);

        // T2/T1 = (p2/p1)^((gamma-1)/gamma)
        double expected = std::pow(2.0, (GAMMA_AIR - 1.0) / GAMMA_AIR);
        ASSERT_NEAR(T_ratio, expected, TOLERANCE);
        return true;
    });

    run_test("Isentropic: density ratio from pressure ratio", []() {
        double p_ratio = 2.0;
        double rho_ratio = IsentropicFlow::densityRatio(p_ratio, GAMMA_AIR);

        // rho2/rho1 = (p2/p1)^(1/gamma)
        double expected = std::pow(2.0, 1.0 / GAMMA_AIR);
        ASSERT_NEAR(rho_ratio, expected, TOLERANCE);
        return true;
    });

    run_test("Isentropic: pressure from density", []() {
        double p1 = 101325.0;  // Pa
        double rho1 = 1.225;   // kg/m^3
        double rho2 = 2.45;    // Double density

        double p2 = IsentropicFlow::pressure(p1, rho1, rho2, GAMMA_AIR);

        // p2 = p1 * (rho2/rho1)^gamma
        double expected = p1 * std::pow(2.0, GAMMA_AIR);
        ASSERT_NEAR(p2, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Isentropic: speed of sound from pressure and density", []() {
        double p = 101325.0;  // Pa
        double rho = 1.225;   // kg/m^3

        double c = IsentropicFlow::soundSpeed(GAMMA_AIR, p, rho);

        // c = sqrt(gamma * p / rho) ~ 340 m/s for air
        double expected = std::sqrt(GAMMA_AIR * p / rho);
        ASSERT_NEAR(c, expected, TOLERANCE);
        ASSERT_NEAR(c, 340.0, 5.0);  // Approximately 340 m/s
        return true;
    });

    run_test("Isentropic: speed of sound from temperature", []() {
        double T = 288.0;  // K (15°C)
        double R = 287.0;  // J/(kg*K) for air

        double c = IsentropicFlow::soundSpeedFromTemp(GAMMA_AIR, R, T);

        // c = sqrt(gamma * R * T) ~ 340 m/s
        ASSERT_NEAR(c, 340.0, 5.0);
        return true;
    });

    run_test("Isentropic: ideal gas consistency", []() {
        double p = 101325.0;
        double rho = 1.225;
        double T = 288.0;
        double R = 287.0;

        double c1 = IsentropicFlow::soundSpeed(GAMMA_AIR, p, rho);
        double c2 = IsentropicFlow::soundSpeedFromTemp(GAMMA_AIR, R, T);

        // Both methods should give same result
        ASSERT_NEAR(c1, c2, 1.0);
        return true;
    });

    // ========================================
    // Stagnation Properties Tests
    // ========================================

    run_test("Stagnation: temperature at rest", []() {
        double T = 288.0;  // K
        double M = 0.0;    // At rest

        double T0 = StagnationProperties::temperature(T, M, GAMMA_AIR);

        // At rest, T0 = T
        ASSERT_NEAR(T0, T, TOLERANCE);
        return true;
    });

    run_test("Stagnation: temperature at Mach 1", []() {
        double T = 288.0;
        double M = 1.0;

        double T0 = StagnationProperties::temperature(T, M, GAMMA_AIR);

        // T0 = T * (1 + (gamma-1)/2 * M^2) = T * 1.2
        double expected = T * (1.0 + 0.5 * (GAMMA_AIR - 1.0));
        ASSERT_NEAR(T0, expected, TOLERANCE);
        return true;
    });

    run_test("Stagnation: pressure at Mach 2", []() {
        double p = 101325.0;
        double M = 2.0;

        double p0 = StagnationProperties::pressure(p, M, GAMMA_AIR);

        // p0/p = (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1))
        double factor = 1.0 + 0.5 * (GAMMA_AIR - 1.0) * M * M;
        double expected = p * std::pow(factor, GAMMA_AIR / (GAMMA_AIR - 1.0));
        ASSERT_NEAR(p0, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Stagnation: density increases with Mach", []() {
        double rho = 1.225;

        double rho0_M0 = StagnationProperties::density(rho, 0.0, GAMMA_AIR);
        double rho0_M1 = StagnationProperties::density(rho, 1.0, GAMMA_AIR);
        double rho0_M2 = StagnationProperties::density(rho, 2.0, GAMMA_AIR);

        ASSERT_TRUE(rho0_M1 > rho0_M0);
        ASSERT_TRUE(rho0_M2 > rho0_M1);
        return true;
    });

    run_test("Stagnation: static from stagnation temperature", []() {
        double T0 = 345.6;  // K
        double M = 1.0;

        double T = StagnationProperties::staticTemperature(T0, M, GAMMA_AIR);
        double T0_back = StagnationProperties::temperature(T, M, GAMMA_AIR);

        // Should recover T0
        ASSERT_NEAR(T0_back, T0, TOLERANCE);
        return true;
    });

    run_test("Stagnation: static from stagnation pressure", []() {
        double p0 = 200000.0;
        double M = 0.8;

        double p = StagnationProperties::staticPressure(p0, M, GAMMA_AIR);
        double p0_back = StagnationProperties::pressure(p, M, GAMMA_AIR);

        // Should recover p0
        ASSERT_NEAR(p0_back, p0, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Critical Flow Tests
    // ========================================

    run_test("Critical: pressure ratio for air", []() {
        double p_star_ratio = CriticalFlow::criticalPressureRatio(GAMMA_AIR);

        // For gamma = 1.4: (2/(gamma+1))^(gamma/(gamma-1)) ≈ 0.528
        double expected = std::pow(2.0 / (GAMMA_AIR + 1.0), GAMMA_AIR / (GAMMA_AIR - 1.0));
        ASSERT_NEAR(p_star_ratio, expected, TOLERANCE);
        ASSERT_NEAR(p_star_ratio, 0.528, 0.01);
        return true;
    });

    run_test("Critical: temperature ratio", []() {
        double T_star_ratio = CriticalFlow::criticalTemperatureRatio(GAMMA_AIR);

        // T*/T0 = 2/(gamma+1) = 2/2.4 ≈ 0.833
        double expected = 2.0 / (GAMMA_AIR + 1.0);
        ASSERT_NEAR(T_star_ratio, expected, TOLERANCE);
        return true;
    });

    run_test("Critical: density ratio", []() {
        double rho_star_ratio = CriticalFlow::criticalDensityRatio(GAMMA_AIR);

        // rho*/rho0 = (2/(gamma+1))^(1/(gamma-1))
        double expected = std::pow(2.0 / (GAMMA_AIR + 1.0), 1.0 / (GAMMA_AIR - 1.0));
        ASSERT_NEAR(rho_star_ratio, expected, TOLERANCE);
        return true;
    });

    run_test("Critical: sonic velocity", []() {
        double T0 = 300.0;  // K
        double R = 287.0;   // J/(kg*K)

        double c_star = CriticalFlow::criticalVelocity(GAMMA_AIR, R, T0);

        // c* = sqrt(2*gamma*R*T0/(gamma+1))
        double expected = std::sqrt(2.0 * GAMMA_AIR * R * T0 / (GAMMA_AIR + 1.0));
        ASSERT_NEAR(c_star, expected, TOLERANCE);
        return true;
    });

    run_test("Critical: choking condition", []() {
        // At M = 1, critical conditions apply
        double p_ratio = CriticalFlow::criticalPressureRatio(GAMMA_AIR);

        // Critical pressure ratio < 1
        ASSERT_TRUE(p_ratio < 1.0);
        ASSERT_TRUE(p_ratio > 0.5);
        return true;
    });

    // ========================================
    // Normal Shock Tests
    // ========================================

    run_test("Normal shock: downstream Mach at M1 = 2", []() {
        double M1 = 2.0;
        double M2 = NormalShock::downstreamMach(M1, GAMMA_AIR);

        // M2 < 1 (subsonic after shock)
        ASSERT_TRUE(M2 < 1.0);
        ASSERT_NEAR(M2, 0.577, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Normal shock: pressure ratio increases", []() {
        double M1 = 2.0;
        double p_ratio = NormalShock::pressureRatio(M1, GAMMA_AIR);

        // Pressure jumps across shock
        ASSERT_TRUE(p_ratio > 1.0);
        ASSERT_NEAR(p_ratio, 4.5, 0.1);
        return true;
    });

    run_test("Normal shock: density ratio", []() {
        double M1 = 2.0;
        double rho_ratio = NormalShock::densityRatio(M1, GAMMA_AIR);

        // Density increases (but less than pressure)
        ASSERT_TRUE(rho_ratio > 1.0);
        ASSERT_NEAR(rho_ratio, 2.667, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Normal shock: temperature ratio", []() {
        double M1 = 2.0;
        double T_ratio = NormalShock::temperatureRatio(M1, GAMMA_AIR);

        // Temperature increases across shock
        ASSERT_TRUE(T_ratio > 1.0);
        return true;
    });

    run_test("Normal shock: stagnation pressure loss", []() {
        double M1 = 2.0;
        double p0_ratio = NormalShock::stagnationPressureRatio(M1, GAMMA_AIR);

        // Stagnation pressure decreases (entropy increase)
        ASSERT_TRUE(p0_ratio < 1.0);
        ASSERT_NEAR(p0_ratio, 0.7209, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Normal shock: weak shock limit", []() {
        double M1 = 1.1;  // Weak shock
        double M2 = NormalShock::downstreamMach(M1, GAMMA_AIR);

        // M2 should be close to 1
        ASSERT_TRUE(M2 < 1.0 && M2 > 0.9);
        return true;
    });

    run_test("Normal shock: strong shock", []() {
        double M1 = 5.0;
        bool is_strong = NormalShock::isStrongShock(M1);

        ASSERT_TRUE(is_strong);  // M1 > 3
        return true;
    });

    run_test("Normal shock: entropy increase", []() {
        double M1 = 2.0;
        double cp = 1005.0;  // J/(kg*K) for air

        double delta_s = NormalShock::entropyIncrease(M1, GAMMA_AIR, cp);

        // Entropy always increases across shock
        ASSERT_TRUE(delta_s > 0.0);
        return true;
    });

    run_test("Normal shock: shock strength parameter", []() {
        double M1 = 2.0;
        double beta = NormalShock::shockStrength(M1, GAMMA_AIR);

        // beta = (p2 - p1)/p1 > 0
        ASSERT_TRUE(beta > 0.0);
        return true;
    });

    run_test("Normal shock: Rankine-Hugoniot relations consistency", []() {
        double M1 = 3.0;

        double p_ratio = NormalShock::pressureRatio(M1, GAMMA_AIR);
        double rho_ratio = NormalShock::densityRatio(M1, GAMMA_AIR);
        double T_ratio = NormalShock::temperatureRatio(M1, GAMMA_AIR);

        // T_ratio should equal p_ratio / rho_ratio (ideal gas)
        double T_calc = p_ratio / rho_ratio;
        ASSERT_NEAR(T_ratio, T_calc, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Oblique Shock Tests
    // ========================================

    run_test("Oblique shock: normal Mach component", []() {
        double M1 = 2.0;
        double beta = M_PI / 4.0;  // 45 degrees

        double M1n = ObliqueShock::normalMach(M1, beta);

        // M1n = M1 * sin(beta)
        double expected = M1 * std::sin(beta);
        ASSERT_NEAR(M1n, expected, TOLERANCE);
        return true;
    });

    run_test("Oblique shock: shock angle calculation", []() {
        double M1 = 2.0;
        double theta = 10.0 * M_PI / 180.0;  // 10 degrees deflection

        double beta = ObliqueShock::shockAngle(M1, theta, GAMMA_AIR);

        // Shock angle > deflection angle
        ASSERT_TRUE(beta > theta);
        return true;
    });

    run_test("Oblique shock: downstream Mach", []() {
        double M1 = 2.0;
        double beta = 45.0 * M_PI / 180.0;
        double theta = 10.0 * M_PI / 180.0;

        double M2 = ObliqueShock::downstreamMach(M1, beta, theta, GAMMA_AIR);

        // M2 should be subsonic for this case
        ASSERT_TRUE(M2 > 0.0);
        return true;
    });

    run_test("Oblique shock: maximum deflection angle", []() {
        double M1 = 2.0;
        double theta_max = ObliqueShock::maxDeflectionAngle(M1, GAMMA_AIR);

        // Should be positive and less than 90 degrees
        ASSERT_TRUE(theta_max > 0.0);
        ASSERT_TRUE(theta_max < M_PI / 2.0);
        return true;
    });

    run_test("Oblique shock: weak shock solution exists", []() {
        double M1 = 2.0;
        double theta = 5.0 * M_PI / 180.0;  // Small deflection

        double beta = ObliqueShock::shockAngle(M1, theta, GAMMA_AIR);

        // Weak shock has smaller beta
        ASSERT_TRUE(beta < M_PI / 2.0);
        return true;
    });

    // ========================================
    // Prandtl-Meyer Expansion Tests
    // ========================================

    run_test("Prandtl-Meyer: angle at M = 1", []() {
        double M = 1.0;
        double nu = PrandtlMeyerExpansion::prandtlMeyerAngle(M, GAMMA_AIR);

        // At sonic condition, nu = 0
        ASSERT_NEAR(nu, 0.0, TOLERANCE);
        return true;
    });

    run_test("Prandtl-Meyer: angle increases with Mach", []() {
        double nu1 = PrandtlMeyerExpansion::prandtlMeyerAngle(2.0, GAMMA_AIR);
        double nu2 = PrandtlMeyerExpansion::prandtlMeyerAngle(3.0, GAMMA_AIR);

        // nu increases with M
        ASSERT_TRUE(nu2 > nu1);
        return true;
    });

    run_test("Prandtl-Meyer: downstream Mach from expansion", []() {
        double M1 = 2.0;
        double theta = 10.0 * M_PI / 180.0;  // 10 degree expansion

        double M2 = PrandtlMeyerExpansion::downstreamMach(M1, theta, GAMMA_AIR);

        // Expansion increases Mach number
        ASSERT_TRUE(M2 > M1);
        return true;
    });

    run_test("Prandtl-Meyer: maximum expansion angle", []() {
        double nu_max = PrandtlMeyerExpansion::maxPrandtlMeyerAngle(GAMMA_AIR);

        // For gamma = 1.4: nu_max ~ 130 degrees
        ASSERT_TRUE(nu_max > 2.0);  // > 115 degrees
        ASSERT_TRUE(nu_max < 2.5);  // < 143 degrees
        return true;
    });

    run_test("Prandtl-Meyer: isentropic expansion", []() {
        double M1 = 1.5;
        double M2 = 2.0;

        double nu1 = PrandtlMeyerExpansion::prandtlMeyerAngle(M1, GAMMA_AIR);
        double nu2 = PrandtlMeyerExpansion::prandtlMeyerAngle(M2, GAMMA_AIR);

        double theta = nu2 - nu1;

        // Positive expansion angle
        ASSERT_TRUE(theta > 0.0);
        return true;
    });

    run_test("Prandtl-Meyer: expansion to high Mach", []() {
        double M1 = 2.0;
        double theta = 30.0 * M_PI / 180.0;  // Large expansion

        double M2 = PrandtlMeyerExpansion::downstreamMach(M1, theta, GAMMA_AIR);

        // Should reach higher Mach
        ASSERT_TRUE(M2 > 3.0);
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Physical: sound speed increases with temperature", []() {
        double R = 287.0;

        double c1 = IsentropicFlow::soundSpeedFromTemp(GAMMA_AIR, R, 250.0);
        double c2 = IsentropicFlow::soundSpeedFromTemp(GAMMA_AIR, R, 300.0);

        ASSERT_TRUE(c2 > c1);
        return true;
    });

    run_test("Physical: stagnation properties always greater or equal", []() {
        double p = 100000.0;
        double M = 0.8;

        double p0 = StagnationProperties::pressure(p, M, GAMMA_AIR);

        // p0 >= p always
        ASSERT_TRUE(p0 >= p);
        return true;
    });

    run_test("Physical: normal shock always decelerates flow", []() {
        double M1 = 2.5;
        double M2 = NormalShock::downstreamMach(M1, GAMMA_AIR);

        // Supersonic to subsonic
        ASSERT_TRUE(M1 > 1.0);
        ASSERT_TRUE(M2 < 1.0);
        return true;
    });

    run_test("Physical: shock increases pressure and temperature", []() {
        double M1 = 2.0;

        double p_ratio = NormalShock::pressureRatio(M1, GAMMA_AIR);
        double T_ratio = NormalShock::temperatureRatio(M1, GAMMA_AIR);

        ASSERT_TRUE(p_ratio > 1.0);
        ASSERT_TRUE(T_ratio > 1.0);
        return true;
    });

    run_test("Physical: expansion wave decreases pressure and temperature", []() {
        // Expansion is isentropic
        double p1 = 200000.0;
        double rho1 = 2.0;
        double rho2 = 1.0;  // Expansion reduces density

        double p2 = IsentropicFlow::pressure(p1, rho1, rho2, GAMMA_AIR);

        // Pressure decreases
        ASSERT_TRUE(p2 < p1);
        return true;
    });

    run_test("Physical: critical pressure ratio is universal for given gamma", []() {
        double p_star1 = CriticalFlow::criticalPressureRatio(1.4);
        double p_star2 = CriticalFlow::criticalPressureRatio(1.4);

        // Should be identical
        ASSERT_NEAR(p_star1, p_star2, TOLERANCE);
        return true;
    });

    run_test("Physical: oblique shock reduces to normal at 90 degrees", []() {
        double M1 = 2.0;
        double beta = M_PI / 2.0;  // Normal shock

        double M1n = ObliqueShock::normalMach(M1, beta);

        // M1n = M1 when beta = 90 degrees
        ASSERT_NEAR(M1n, M1, TOLERANCE);
        return true;
    });

    run_test("Limit: weak shock has small pressure jump", []() {
        double M1 = 1.05;  // Barely supersonic
        double p_ratio = NormalShock::pressureRatio(M1, GAMMA_AIR);

        // Weak shock: p_ratio ~ 1
        ASSERT_TRUE(p_ratio < 1.2);
        return true;
    });

    run_test("Limit: strong shock density ratio approaches limit", []() {
        double M1 = 10.0;  // Very strong shock
        double rho_ratio = NormalShock::densityRatio(M1, GAMMA_AIR);

        // Limit: rho2/rho1 → (gamma+1)/(gamma-1) = 6 for gamma=1.4
        double limit = (GAMMA_AIR + 1.0) / (GAMMA_AIR - 1.0);
        ASSERT_TRUE(rho_ratio < limit);
        ASSERT_TRUE(rho_ratio > 0.9 * limit);  // Approaching limit
        return true;
    });

    run_test("Limit: Prandtl-Meyer angle approaches maximum", []() {
        double M_high = 20.0;
        double nu = PrandtlMeyerExpansion::prandtlMeyerAngle(M_high, GAMMA_AIR);
        double nu_max = PrandtlMeyerExpansion::maxPrandtlMeyerAngle(GAMMA_AIR);

        // Should be close to maximum
        ASSERT_TRUE(nu < nu_max);
        ASSERT_TRUE(nu > 0.95 * nu_max);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Compressible Flow" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All compressible flow tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Isentropic flow relations" << std::endl;
        std::cout << "  - Stagnation properties calculations" << std::endl;
        std::cout << "  - Critical flow conditions (M = 1)" << std::endl;
        std::cout << "  - Normal shock relations (Rankine-Hugoniot)" << std::endl;
        std::cout << "  - Oblique shock wave theory" << std::endl;
        std::cout << "  - Prandtl-Meyer expansion waves" << std::endl;
        std::cout << "  - Mach number dependencies" << std::endl;
        std::cout << "  - Physical consistency and limiting cases" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
