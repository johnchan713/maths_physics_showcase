/**
 * Phase 4 Validation: Gravitation and Universal Gravitation
 *
 * Tests the gravitation.hpp module functions.
 *
 * Coverage:
 * - Newton's law of universal gravitation F = GMm/r²
 * - Gravitational field strength g = GM/r²
 * - Moon's motion and inverse square law verification
 * - Surface gravity calculation (g ≈ 9.81 m/s²)
 * - Earth's mass determination
 * - Kepler's Third Law T² ∝ r³
 * - Gravitational potential energy U = -GMm/r
 * - Escape velocity v = √(2GM/r)
 * - Binary systems and reduced mass
 * - Center of mass calculations
 */

#include <iostream>
#include <cmath>
#include "../include/physics/gravitation.hpp"

using namespace physics::gravitation;

// Test tolerances
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

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Gravitation Validation ===" << std::endl;
    std::cout << std::endl;

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
    // Universal Gravitation Force
    // ========================================

    run_test("Universal gravitation force formula", []() {
        double m1 = 1000.0;  // 1000 kg
        double m2 = 500.0;   // 500 kg
        double r = 1.0;      // 1 meter

        double F = universalGravitationForce(m1, m2, r);

        double expected = constants::G * m1 * m2 / (r * r);
        ASSERT_NEAR(F, expected, TOLERANCE);
        return true;
    });

    run_test("Gravitation force inverse square law", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double r1 = 1.0;
        double r2 = 2.0;

        double F1 = universalGravitationForce(m1, m2, r1);
        double F2 = universalGravitationForce(m1, m2, r2);

        // F ∝ 1/r², so F1/F2 = (r2/r1)² = 4
        ASSERT_NEAR(F1 / F2, 4.0, TOLERANCE);
        return true;
    });

    run_test("Gravitation force scales with mass product", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double r = 1.0;

        double F1 = universalGravitationForce(m1, m2, r);
        double F2 = universalGravitationForce(2.0 * m1, m2, r);

        ASSERT_NEAR(F2, 2.0 * F1, TOLERANCE);
        return true;
    });

    // ========================================
    // Gravitational Field Strength
    // ========================================

    run_test("Gravitational field strength formula", []() {
        double M = 1e6;  // 1 million kg
        double r = 10.0;  // 10 meters

        double g = gravitationalFieldStrength(M, r);

        double expected = constants::G * M / (r * r);
        ASSERT_NEAR(g, expected, TOLERANCE);
        return true;
    });

    run_test("Field strength inverse square law", []() {
        double M = 1e6;
        double r1 = 10.0;
        double r2 = 20.0;

        double g1 = gravitationalFieldStrength(M, r1);
        double g2 = gravitationalFieldStrength(M, r2);

        ASSERT_NEAR(g1 / g2, 4.0, TOLERANCE);
        return true;
    });

    run_test("Earth surface gravity approximately 9.81", []() {
        double g = calculateSurfaceGravity();

        ASSERT_NEAR(g, 9.81, 0.1);  // Within 0.1 m/s²
        return true;
    });

    run_test("Surface gravity from Earth constants", []() {
        double g = gravitationalFieldStrength(constants::EARTH_MASS,
                                             constants::EARTH_RADIUS);

        ASSERT_NEAR(g, 9.81, 0.1);
        return true;
    });

    // ========================================
    // Moon's Motion and Inverse Square Law
    // ========================================

    run_test("Moon centripetal acceleration", []() {
        double a_c = calculateMoonCentripetalAccel();

        // Should be about 0.00272 m/s²
        ASSERT_TRUE(a_c > 0.002 && a_c < 0.003);
        return true;
    });

    run_test("Moon acceleration from gravitation", []() {
        double g_moon = gravitationalFieldStrength(constants::EARTH_MASS,
                                                   constants::MOON_ORBITAL_RADIUS);

        double a_c = calculateMoonCentripetalAccel();

        // These should match (Newton's insight!)
        ASSERT_NEAR(g_moon, a_c, 1e-4);
        return true;
    });

    run_test("Inverse square law verification", []() {
        double ratio = verifyInverseSquareLaw();

        // Should equal (r_moon/R_earth)² ≈ 3600
        double expected_ratio = std::pow(calculateDistanceRatio(), 2.0);

        ASSERT_NEAR(ratio, expected_ratio, 10.0);
        return true;
    });

    run_test("Distance ratio Moon to Earth", []() {
        double ratio = calculateDistanceRatio();

        // r_moon / R_earth ≈ 60
        ASSERT_TRUE(ratio > 59.0 && ratio < 61.0);
        return true;
    });

    run_test("Surface to Moon gravity ratio", []() {
        double g_surface = calculateSurfaceGravity();
        double g_moon = gravitationalFieldStrength(constants::EARTH_MASS,
                                                   constants::MOON_ORBITAL_RADIUS);

        double ratio = g_surface / g_moon;

        // Should be about 3600
        ASSERT_TRUE(ratio > 3500.0 && ratio < 3700.0);
        return true;
    });

    // ========================================
    // Mass of the Earth
    // ========================================

    run_test("Earth mass from surface gravity", []() {
        double M = calculateEarthMass();

        // Should be close to 5.972e24 kg
        double error = std::abs(M - constants::EARTH_MASS) / constants::EARTH_MASS;
        ASSERT_TRUE(error < 0.01);  // Within 1%
        return true;
    });

    run_test("Earth mass from Moon's orbit", []() {
        double M = calculateEarthMassFromMoon();

        // Should be close to 5.972e24 kg
        double error = std::abs(M - constants::EARTH_MASS) / constants::EARTH_MASS;
        ASSERT_TRUE(error < 0.05);  // Within 5%
        return true;
    });

    run_test("Two methods agree for Earth mass", []() {
        double M1 = calculateEarthMass();
        double M2 = calculateEarthMassFromMoon();

        double diff = std::abs(M1 - M2) / M1;
        ASSERT_TRUE(diff < 0.02);  // Within 2%
        return true;
    });

    run_test("Central mass from arbitrary orbit", []() {
        double r = 1e7;  // 10,000 km
        double T = 5400.0;  // 90 minutes

        double M = calculateCentralMassFromOrbit(r, T);

        ASSERT_TRUE(M > 0);
        return true;
    });

    // ========================================
    // Kepler's Third Law
    // ========================================

    run_test("Kepler Third Law for Moon", []() {
        double T = calculatePeriodKeplerThird(constants::MOON_ORBITAL_RADIUS,
                                             constants::EARTH_MASS);

        // Should match Moon's actual period
        double error = std::abs(T - constants::MOON_ORBITAL_PERIOD) /
                      constants::MOON_ORBITAL_PERIOD;
        ASSERT_TRUE(error < 0.01);
        return true;
    });

    run_test("Kepler constant for Earth", []() {
        double K = calculateKeplerConstant(constants::EARTH_MASS);

        // K = 4π²/(GM)
        double expected = (4.0 * M_PI * M_PI) / (constants::G * constants::EARTH_MASS);
        ASSERT_NEAR(K, expected, TOLERANCE);
        return true;
    });

    run_test("Verify Kepler Third Law Moon orbit", []() {
        bool valid = verifyKeplerThirdLaw(constants::MOON_ORBITAL_RADIUS,
                                         constants::MOON_ORBITAL_PERIOD,
                                         constants::EARTH_MASS,
                                         constants::G,
                                         0.05);  // Looser tolerance

        ASSERT_TRUE(valid);
        return true;
    });

    run_test("Kepler Third Law T squared proportional to r cubed", []() {
        double M = constants::EARTH_MASS;
        double r1 = 1e7;
        double r2 = 2e7;

        double T1 = calculatePeriodKeplerThird(r1, M);
        double T2 = calculatePeriodKeplerThird(r2, M);

        // T1²/T2² = r1³/r2³
        double ratio_T2 = (T1 * T1) / (T2 * T2);
        double ratio_r3 = (r1 * r1 * r1) / (r2 * r2 * r2);

        ASSERT_NEAR(ratio_T2, ratio_r3, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Compare satellite periods", []() {
        double r1 = 1e7;
        double r2 = 8e7;  // 8 times farther

        double ratio = compareSatellitePeriods(r1, r2);

        // T1/T2 = (r1/r2)^(3/2) = (1/8)^(3/2) = 1/(8√8) ≈ 0.044
        double expected = std::pow(1.0 / 8.0, 1.5);

        ASSERT_NEAR(ratio, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Period increases with orbital radius", []() {
        double M = constants::EARTH_MASS;
        double r1 = 1e7;
        double r2 = 2e7;

        double T1 = calculatePeriodKeplerThird(r1, M);
        double T2 = calculatePeriodKeplerThird(r2, M);

        ASSERT_TRUE(T2 > T1);
        return true;
    });

    // ========================================
    // Gravitational Potential Energy
    // ========================================

    run_test("Gravitational potential energy is negative", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double r = 1.0;

        double U = gravitationalPotentialEnergy(m1, m2, r);

        ASSERT_TRUE(U < 0);
        return true;
    });

    run_test("Potential energy formula", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double r = 1.0;

        double U = gravitationalPotentialEnergy(m1, m2, r);

        double expected = -(constants::G * m1 * m2) / r;
        ASSERT_NEAR(U, expected, TOLERANCE);
        return true;
    });

    run_test("Potential energy increases with distance", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double r1 = 1.0;
        double r2 = 2.0;

        double U1 = gravitationalPotentialEnergy(m1, m2, r1);
        double U2 = gravitationalPotentialEnergy(m1, m2, r2);

        // U2 > U1 (both negative, U2 closer to zero)
        ASSERT_TRUE(U2 > U1);
        return true;
    });

    run_test("Potential energy approaches zero at infinity", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double r = 1e10;  // Very large distance

        double U = gravitationalPotentialEnergy(m1, m2, r);

        ASSERT_TRUE(std::abs(U) < 1e-10);  // Much smaller than typical values
        return true;
    });

    // ========================================
    // Escape Velocity
    // ========================================

    run_test("Earth escape velocity", []() {
        double v_esc = calculateEscapeVelocity(constants::EARTH_MASS,
                                              constants::EARTH_RADIUS);

        // Should be about 11.2 km/s = 11200 m/s
        ASSERT_TRUE(v_esc > 11000.0 && v_esc < 11500.0);
        return true;
    });

    run_test("Escape velocity formula", []() {
        double M = 1e24;
        double R = 1e6;

        double v = calculateEscapeVelocity(M, R);

        double expected = std::sqrt(2.0 * constants::G * M / R);
        ASSERT_NEAR(v, expected, TOLERANCE);
        return true;
    });

    run_test("Escape velocity energy conservation", []() {
        double M = constants::EARTH_MASS;
        double R = constants::EARTH_RADIUS;
        double m = 1.0;  // 1 kg test mass

        double v_esc = calculateEscapeVelocity(M, R);

        // Kinetic energy at escape
        double KE = 0.5 * m * v_esc * v_esc;

        // Potential energy at surface
        double PE = gravitationalPotentialEnergy(M, m, R);

        // Total energy should be zero (barely escapes)
        ASSERT_NEAR(KE + PE, 0.0, 1.0);
        return true;
    });

    run_test("Escape velocity scales with mass and radius", []() {
        double M = 1e24;
        double R = 1e6;

        double v1 = calculateEscapeVelocity(M, R);
        double v2 = calculateEscapeVelocity(4.0 * M, 2.0 * R);

        // v ∝ √(M/R), so v2/v1 = √(4M/(2R)) / √(M/R) = √2
        ASSERT_NEAR(v2 / v1, std::sqrt(2.0), LOOSE_TOLERANCE);
        return true;
    });

    run_test("Moon escape velocity smaller than Earth", []() {
        double v_earth = calculateEscapeVelocity(constants::EARTH_MASS,
                                                constants::EARTH_RADIUS);
        double v_moon = calculateEscapeVelocity(constants::MOON_MASS,
                                               constants::MOON_ORBITAL_RADIUS / 60.0);

        ASSERT_TRUE(v_moon < v_earth);
        return true;
    });

    // ========================================
    // Binary Systems
    // ========================================

    run_test("Reduced mass formula", []() {
        double m1 = 1000.0;
        double m2 = 500.0;

        double mu = calculateReducedMass(m1, m2);

        double expected = (m1 * m2) / (m1 + m2);
        ASSERT_NEAR(mu, expected, TOLERANCE);
        return true;
    });

    run_test("Reduced mass less than either mass", []() {
        double m1 = 1000.0;
        double m2 = 500.0;

        double mu = calculateReducedMass(m1, m2);

        ASSERT_TRUE(mu < m1);
        ASSERT_TRUE(mu < m2);
        return true;
    });

    run_test("Reduced mass symmetric", []() {
        double m1 = 1000.0;
        double m2 = 500.0;

        double mu1 = calculateReducedMass(m1, m2);
        double mu2 = calculateReducedMass(m2, m1);

        ASSERT_NEAR(mu1, mu2, TOLERANCE);
        return true;
    });

    run_test("Equal masses reduced mass", []() {
        double m = 1000.0;

        double mu = calculateReducedMass(m, m);

        ASSERT_NEAR(mu, m / 2.0, TOLERANCE);
        return true;
    });

    run_test("Center of mass distance", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double d = 1.5;  // 1.5 meters separation

        double r1 = calculateCMDistance(m1, m2, d);

        // r1 = m2/(m1+m2) × d = 500/1500 × 1.5 = 0.5
        ASSERT_NEAR(r1, 0.5, TOLERANCE);
        return true;
    });

    run_test("Center of mass distances sum to separation", []() {
        double m1 = 1000.0;
        double m2 = 500.0;
        double d = 1.5;

        double r1 = calculateCMDistance(m1, m2, d);
        double r2 = calculateCMDistance(m2, m1, d);

        ASSERT_NEAR(r1 + r2, d, TOLERANCE);
        return true;
    });

    run_test("Equal masses center of mass at midpoint", []() {
        double m = 1000.0;
        double d = 2.0;

        double r = calculateCMDistance(m, m, d);

        ASSERT_NEAR(r, 1.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Test Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Total tests:  " << (tests_passed + tests_failed) << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    return tests_failed == 0 ? 0 : 1;
}
