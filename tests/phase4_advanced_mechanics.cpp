/**
 * Phase 4 Validation: Advanced Mechanics
 *
 * Tests the advanced_mechanics.hpp module functions.
 *
 * Coverage:
 * - Polar coordinates (Cartesian ↔ polar conversion)
 * - Radial and tangential velocities and accelerations
 * - Relative motion and Galilean transformations
 * - Conservative forces (∇ × F = 0)
 * - Potential energy (gravitational, escape velocity)
 * - Orbital mechanics (Kepler's laws, vis-viva equation)
 * - Virial theorem (2⟨T⟩ = -⟨U⟩ for inverse square law)
 * - Variational calculus (Lagrangian, Hamiltonian, action)
 */

#include <iostream>
#include <cmath>
#include <array>
#include <functional>
#include "../include/physics/advanced_mechanics.hpp"

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

using namespace physics::advanced_mechanics;
using namespace physics::advanced_mechanics::constants;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Advanced Mechanics Validation ===" << std::endl;
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
    // Polar Coordinates Tests
    // ========================================

    run_test("Polar coordinates: Cartesian to polar at origin", []() {
        auto [r, theta] = cartesianToPolar(0.0, 0.0);
        ASSERT_NEAR(r, 0.0, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: Cartesian to polar (1,0) = (1,0°)", []() {
        auto [r, theta] = cartesianToPolar(1.0, 0.0);
        ASSERT_NEAR(r, 1.0, TOLERANCE);
        ASSERT_NEAR(theta, 0.0, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: Cartesian to polar (0,1) = (1,90°)", []() {
        auto [r, theta] = cartesianToPolar(0.0, 1.0);
        ASSERT_NEAR(r, 1.0, TOLERANCE);
        ASSERT_NEAR(theta, M_PI / 2.0, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: Cartesian to polar (3,4) = (5,atan(4/3))", []() {
        auto [r, theta] = cartesianToPolar(3.0, 4.0);
        ASSERT_NEAR(r, 5.0, TOLERANCE);
        ASSERT_NEAR(theta, std::atan2(4.0, 3.0), TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: polar to Cartesian (1,0°) = (1,0)", []() {
        auto [x, y] = polarToCartesian(1.0, 0.0);
        ASSERT_NEAR(x, 1.0, TOLERANCE);
        ASSERT_NEAR(y, 0.0, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: polar to Cartesian (1,90°) = (0,1)", []() {
        auto [x, y] = polarToCartesian(1.0, M_PI / 2.0);
        ASSERT_NEAR(x, 0.0, TOLERANCE);
        ASSERT_NEAR(y, 1.0, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: round-trip conversion", []() {
        double x0 = 3.5, y0 = 2.7;
        auto [r, theta] = cartesianToPolar(x0, y0);
        auto [x, y] = polarToCartesian(r, theta);
        ASSERT_NEAR(x, x0, TOLERANCE);
        ASSERT_NEAR(y, y0, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: radial velocity equals dr/dt", []() {
        double dr_dt = 5.0;  // m/s
        double v_r = radialVelocity(dr_dt);
        ASSERT_NEAR(v_r, dr_dt, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: tangential velocity v_theta = r*omega", []() {
        double r = 2.0;      // m
        double omega = 3.0;  // rad/s
        double v_theta = tangentialVelocity(r, omega);
        ASSERT_NEAR(v_theta, r * omega, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: radial acceleration includes centripetal term", []() {
        double d2r_dt2 = 0.0;
        double r = 1.0;
        double omega = 2.0;
        double a_r = radialAcceleration(d2r_dt2, r, omega);
        // For circular motion: a_r = -r*omega^2
        ASSERT_NEAR(a_r, -r * omega * omega, TOLERANCE);
        return true;
    });

    run_test("Polar coordinates: tangential acceleration with Coriolis term", []() {
        double r = 1.0;
        double alpha = 1.0;   // rad/s²
        double dr_dt = 2.0;
        double omega = 3.0;
        double a_theta = tangentialAcceleration(r, alpha, dr_dt, omega);
        double expected = r * alpha + 2.0 * dr_dt * omega;
        ASSERT_NEAR(a_theta, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Relative Motion Tests
    // ========================================

    run_test("Relative motion: velocity difference", []() {
        double v_A = 10.0;  // m/s
        double v_B = 3.0;
        double v_rel = relativeVelocity(v_A, v_B);
        ASSERT_NEAR(v_rel, 7.0, TOLERANCE);
        return true;
    });

    run_test("Relative motion: Galilean transformation", []() {
        double v_lab = 20.0;
        double v_frame = 5.0;
        double v_moving = velocityInMovingFrame(v_lab, v_frame);
        ASSERT_NEAR(v_moving, 15.0, TOLERANCE);
        return true;
    });

    run_test("Relative motion: acceleration difference", []() {
        double a_A = 5.0;
        double a_B = 2.0;
        double a_rel = relativeAcceleration(a_A, a_B);
        ASSERT_NEAR(a_rel, 3.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Conservative Forces Tests
    // ========================================

    run_test("Conservative force: F = -dU/dx", []() {
        double grad_U = 10.0;  // J/m
        double F = forceFromPotential(grad_U);
        ASSERT_NEAR(F, -10.0, TOLERANCE);
        return true;
    });

    run_test("Conservative force: gravitational potential U = -GMm/r", []() {
        double m1 = 1.0e3;   // kg
        double m2 = 2.0e3;   // kg
        double r = 1.0e6;    // m
        double U = gravitationalPotential(m1, m2, r);
        double expected = -GRAVITATIONAL_CONSTANT * m1 * m2 / r;
        ASSERT_NEAR(U, expected, TOLERANCE);
        ASSERT_TRUE(U < 0.0);
        return true;
    });

    run_test("Conservative force: escape velocity v_esc = sqrt(2GM/r)", []() {
        double M = EARTH_MASS;
        double R = 6.371e6;  // Earth radius (m)
        double v_esc = escapeVelocity(M, R);
        double expected = std::sqrt(2.0 * GRAVITATIONAL_CONSTANT * M / R);
        ASSERT_NEAR(v_esc, expected, TOLERANCE);
        // Earth escape velocity ~11.2 km/s
        ASSERT_NEAR(v_esc, 11200.0, 200.0);
        return true;
    });

    run_test("Conservative force: curl = 0 condition", []() {
        // For conservative field: dF_y/dx = dF_x/dy
        double dFy_dx = 5.0;
        double dFx_dy = 5.0;
        ASSERT_TRUE(isConservative(dFy_dx, dFx_dy));
        return true;
    });

    run_test("Conservative force: non-conservative field detection", []() {
        double dFy_dx = 5.0;
        double dFx_dy = 3.0;
        ASSERT_TRUE(!isConservative(dFy_dx, dFx_dy, 1e-10));
        return true;
    });

    run_test("Conservative force: mechanical energy E = K + U", []() {
        double K = 100.0;  // J
        double U = -50.0;  // J
        double E = mechanicalEnergy(K, U);
        ASSERT_NEAR(E, 50.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Orbital Mechanics Tests
    // ========================================

    run_test("Orbital mechanics: angular momentum L = m*r*v", []() {
        double m = 1000.0;   // kg
        double r = 7.0e6;    // m
        double v = 7500.0;   // m/s
        double L = orbitalAngularMomentum(m, r, v);
        ASSERT_NEAR(L, m * r * v, TOLERANCE);
        return true;
    });

    run_test("Orbital mechanics: specific angular momentum h = r*v", []() {
        double r = 7.0e6;
        double v = 7500.0;
        double h = specificAngularMomentum(r, v);
        ASSERT_NEAR(h, r * v, TOLERANCE);
        return true;
    });

    run_test("Orbital mechanics: circular orbit velocity v = sqrt(GM/r)", []() {
        double M = EARTH_MASS;
        double r = 7.0e6;  // Low Earth orbit
        double v = circularOrbitVelocity(M, r);
        double expected = std::sqrt(GRAVITATIONAL_CONSTANT * M / r);
        ASSERT_NEAR(v, expected, TOLERANCE);
        // Should be ~7.5 km/s
        ASSERT_NEAR(v, 7500.0, 100.0);
        return true;
    });

    run_test("Orbital mechanics: circular orbit period T = 2*pi*sqrt(r^3/GM)", []() {
        double M = EARTH_MASS;
        double r = 7.0e6;
        double T = circularOrbitPeriod(r, M);
        double expected = 2.0 * M_PI * std::sqrt(r * r * r / (GRAVITATIONAL_CONSTANT * M));
        ASSERT_NEAR(T, expected, TOLERANCE);
        // Should be ~90 minutes (allow wider tolerance)
        ASSERT_NEAR(T, 5400.0, 500.0);
        return true;
    });

    run_test("Orbital mechanics: specific orbital energy epsilon = -GM/(2a)", []() {
        double M = EARTH_MASS;
        double a = 7.0e6;
        double epsilon = specificOrbitalEnergy(M, a);
        double expected = -GRAVITATIONAL_CONSTANT * M / (2.0 * a);
        ASSERT_NEAR(epsilon, expected, TOLERANCE);
        ASSERT_TRUE(epsilon < 0.0);  // Bound orbit
        return true;
    });

    // ========================================
    // Kepler's Laws Tests
    // ========================================

    run_test("Kepler's third law: T^2 proportional to a^3", []() {
        double M = SOLAR_MASS;
        double a = AU;  // Earth orbit
        double T = keplersThirdLaw(a, M);
        // Should be 1 year = 365.25 days
        double year = 365.25 * 24.0 * 3600.0;
        ASSERT_NEAR(T, year, VERY_LOOSE * year);
        return true;
    });

    run_test("Kepler's third law: inverse calculation of semi-major axis", []() {
        double M = SOLAR_MASS;
        double T = 365.25 * 24.0 * 3600.0;  // 1 year
        double a = semiMajorAxisFromPeriod(T, M);
        ASSERT_NEAR(a, AU, VERY_LOOSE * AU);
        return true;
    });

    run_test("Kepler's laws: eccentricity from perihelion and aphelion", []() {
        double r_p = 1.0e11;  // m
        double r_a = 1.5e11;  // m
        double e = eccentricity(r_a, r_p);
        double expected = (r_a - r_p) / (r_a + r_p);
        ASSERT_NEAR(e, expected, TOLERANCE);
        ASSERT_TRUE(e >= 0.0 && e < 1.0);
        return true;
    });

    run_test("Kepler's laws: circular orbit has e=0", []() {
        double r = 1.0e11;
        double e = eccentricity(r, r);
        ASSERT_NEAR(e, 0.0, TOLERANCE);
        return true;
    });

    run_test("Kepler's laws: semi-major axis a = (r_a + r_p)/2", []() {
        double r_p = 1.0e11;
        double r_a = 1.5e11;
        double a = semiMajorAxis(r_a, r_p);
        ASSERT_NEAR(a, 1.25e11, TOLERANCE);
        return true;
    });

    run_test("Kepler's laws: semi-minor axis b = a*sqrt(1-e^2)", []() {
        double a = 1.0e11;
        double e = 0.6;
        double b = semiMinorAxis(a, e);
        double expected = a * std::sqrt(1.0 - e * e);
        ASSERT_NEAR(b, expected, TOLERANCE);
        return true;
    });

    run_test("Kepler's laws: perihelion distance r_p = a(1-e)", []() {
        double a = 1.0e11;
        double e = 0.2;
        double r_p = perihelionDistance(a, e);
        ASSERT_NEAR(r_p, a * (1.0 - e), TOLERANCE);
        return true;
    });

    run_test("Kepler's laws: aphelion distance r_a = a(1+e)", []() {
        double a = 1.0e11;
        double e = 0.2;
        double r_a = aphelionDistance(a, e);
        ASSERT_NEAR(r_a, a * (1.0 + e), TOLERANCE);
        return true;
    });

    run_test("Kepler's laws: orbital radius at perihelion (nu=0)", []() {
        double a = 1.0e11;
        double e = 0.2;
        double r = orbitalRadius(a, e, 0.0);
        double expected = a * (1.0 - e);
        ASSERT_NEAR(r, expected, TOLERANCE);
        return true;
    });

    run_test("Kepler's laws: orbital radius at aphelion (nu=180°)", []() {
        double a = 1.0e11;
        double e = 0.2;
        double r = orbitalRadius(a, e, M_PI);
        double expected = a * (1.0 + e);
        ASSERT_NEAR(r, expected, TOLERANCE);
        return true;
    });

    run_test("Kepler's laws: vis-viva equation v^2 = GM(2/r - 1/a)", []() {
        double M = EARTH_MASS;
        double r = 7.0e6;
        double a = 8.0e6;
        double v = visVivaEquation(M, r, a);
        double v2_expected = GRAVITATIONAL_CONSTANT * M * (2.0 / r - 1.0 / a);
        ASSERT_NEAR(v * v, v2_expected, LOOSE_TOLERANCE * v2_expected);
        return true;
    });

    run_test("Kepler's laws: vis-viva at perihelion gives maximum velocity", []() {
        double M = SOLAR_MASS;
        double a = AU;
        double e = 0.2;
        double r_p = a * (1.0 - e);
        double r_a = a * (1.0 + e);
        double v_p = visVivaEquation(M, r_p, a);
        double v_a = visVivaEquation(M, r_a, a);
        ASSERT_TRUE(v_p > v_a);
        return true;
    });

    // ========================================
    // Virial Theorem Tests
    // ========================================

    run_test("Virial theorem: 2K + U = 0 for bound orbit", []() {
        double U = -100.0;  // J
        double K = virialTheoremKinetic(U);
        ASSERT_NEAR(K, 50.0, TOLERANCE);
        return true;
    });

    run_test("Virial theorem: total energy E = U/2", []() {
        double U = -100.0;
        double E = virialTheoremEnergy(U);
        ASSERT_NEAR(E, -50.0, TOLERANCE);
        return true;
    });

    run_test("Virial theorem: verification for gravitational system", []() {
        double M = EARTH_MASS;
        double r = 7.0e6;
        double m = 1000.0;
        double v = circularOrbitVelocity(M, r);
        double K = 0.5 * m * v * v;
        double U = gravitationalPotential(M, m, r);
        ASSERT_TRUE(verifyVirialTheorem(K, U, LOOSE_TOLERANCE));
        return true;
    });

    run_test("Virial theorem: K = -E for circular orbit", []() {
        double M = EARTH_MASS;
        double r = 7.0e6;
        double m = 1000.0;
        double v = circularOrbitVelocity(M, r);
        double K = 0.5 * m * v * v;
        double U = gravitationalPotential(M, m, r);
        double E = mechanicalEnergy(K, U);
        ASSERT_NEAR(K, -E, LOOSE_TOLERANCE * K);
        return true;
    });

    // ========================================
    // Variational Calculus Tests
    // ========================================

    run_test("Variational calculus: action for free particle S = (1/2)mv^2*t", []() {
        double m = 1.0;    // kg
        double v = 10.0;   // m/s
        double t = 2.0;    // s
        double S = actionFreeParticle(m, v, t);
        ASSERT_NEAR(S, 0.5 * m * v * v * t, TOLERANCE);
        return true;
    });

    run_test("Variational calculus: Lagrangian L = (1/2)mv^2 for free particle", []() {
        double m = 2.0;
        double v = 5.0;
        double L = lagrangianFreeParticle(m, v);
        ASSERT_NEAR(L, 0.5 * m * v * v, TOLERANCE);
        return true;
    });

    run_test("Variational calculus: Lagrangian L = T - V", []() {
        double T = 100.0;
        double V = 60.0;
        double L = lagrangian(T, V);
        ASSERT_NEAR(L, 40.0, TOLERANCE);
        return true;
    });

    run_test("Variational calculus: Hamiltonian H = T + V", []() {
        double T = 100.0;
        double V = 60.0;
        double H = hamiltonian(T, V);
        ASSERT_NEAR(H, 160.0, TOLERANCE);
        return true;
    });

    run_test("Variational calculus: canonical momentum p = mv", []() {
        double m = 2.0;
        double v = 5.0;
        double p = canonicalMomentum(m, v);
        ASSERT_NEAR(p, m * v, TOLERANCE);
        return true;
    });

    run_test("Variational calculus: Hamiltonian conservation", []() {
        double H1 = 100.0;
        double H2 = 100.0;
        ASSERT_TRUE(isHamiltonianConserved(H1, H2));
        return true;
    });

    run_test("Variational calculus: Hamiltonian not conserved", []() {
        double H1 = 100.0;
        double H2 = 95.0;
        ASSERT_TRUE(!isHamiltonianConserved(H1, H2, 1e-6));
        return true;
    });

    run_test("Variational calculus: reduced mass mu = m1*m2/(m1+m2)", []() {
        double m1 = 2.0;
        double m2 = 3.0;
        double mu = reducedMass(m1, m2);
        double expected = (m1 * m2) / (m1 + m2);
        ASSERT_NEAR(mu, expected, TOLERANCE);
        return true;
    });

    run_test("Variational calculus: reduced mass for equal masses is m/2", []() {
        double m = 4.0;
        double mu = reducedMass(m, m);
        ASSERT_NEAR(mu, m / 2.0, TOLERANCE);
        return true;
    });

    run_test("Variational calculus: effective mass M_eff = m1 + m2", []() {
        double m1 = 2.0;
        double m2 = 3.0;
        double M_eff = effectiveMass(m1, m2);
        ASSERT_NEAR(M_eff, 5.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Tensor Operations Tests
    // ========================================

    run_test("Tensor operations: inertia tensor diagonal element", []() {
        double m = 1.0;
        double x = 3.0;
        double y = 4.0;
        double I_00 = inertiatensorElement2D(m, x, y, 0, 0);
        double r2 = x * x + y * y;
        double expected = m * (r2 - x * x);
        ASSERT_NEAR(I_00, expected, TOLERANCE);
        return true;
    });

    run_test("Tensor operations: inertia tensor off-diagonal element", []() {
        double m = 1.0;
        double x = 3.0;
        double y = 4.0;
        double I_01 = inertiatensorElement2D(m, x, y, 0, 1);
        double expected = -m * x * y;
        ASSERT_NEAR(I_01, expected, TOLERANCE);
        return true;
    });

    run_test("Tensor operations: tensor trace", []() {
        double I_00 = 10.0;
        double I_11 = 15.0;
        double trace = tensorTrace2D(I_00, I_11);
        ASSERT_NEAR(trace, 25.0, TOLERANCE);
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
