/**
 * Phase 4 Validation: Rotational Dynamics
 *
 * Tests the rotational_dynamics.hpp module functions.
 *
 * Coverage:
 * - Angular acceleration α = Δω/Δt
 * - Angular velocity and displacement with constant α
 * - Tangential acceleration a_t = rα
 * - Torque τ = r × F
 * - Angular acceleration from torque α = τ/I
 * - Angular momentum L = Iω
 * - Rotational kinetic energy K_rot = ½Iω²
 * - Moment of inertia formulas (point mass rod cylinder sphere disk hoop)
 * - Parallel axis theorem I = I_cm + Md²
 * - Compound pendulum period T = 2π√(I/mgd)
 * - Center of percussion (sweet spot)
 */

#include <iostream>
#include <cmath>
#include "../include/physics/rotational_dynamics.hpp"

using namespace physics::rotational_dynamics;

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

    std::cout << "=== Phase 4: Rotational Dynamics Validation ===" << std::endl;
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
    // Angular Acceleration
    // ========================================

    run_test("Angular acceleration calculation", []() {
        double omega_i = 10.0;  // 10 rad/s
        double omega_f = 20.0;  // 20 rad/s
        double t = 5.0;         // 5 seconds

        double alpha = calculateAngularAcceleration(omega_i, omega_f, t);

        ASSERT_NEAR(alpha, 2.0, TOLERANCE);  // (20-10)/5 = 2 rad/s²
        return true;
    });

    run_test("Angular acceleration from rest", []() {
        double omega_i = 0.0;
        double omega_f = 30.0;
        double t = 10.0;

        double alpha = calculateAngularAcceleration(omega_i, omega_f, t);

        ASSERT_NEAR(alpha, 3.0, TOLERANCE);
        return true;
    });

    run_test("Final angular velocity", []() {
        double omega_i = 5.0;
        double alpha = 2.0;
        double t = 3.0;

        double omega_f = calculateFinalAngularVelocity(omega_i, alpha, t);

        // ω = ω₀ + αt = 5 + 2×3 = 11
        ASSERT_NEAR(omega_f, 11.0, TOLERANCE);
        return true;
    });

    run_test("Angular displacement with constant acceleration", []() {
        double omega_i = 10.0;
        double alpha = 2.0;
        double t = 5.0;

        double theta = calculateAngularDisplacement(omega_i, alpha, t);

        // θ = ω₀t + ½αt² = 10×5 + ½×2×25 = 75
        ASSERT_NEAR(theta, 75.0, TOLERANCE);
        return true;
    });

    run_test("Tangential acceleration from angular", []() {
        double r = 2.0;     // 2 meters
        double alpha = 3.0;  // 3 rad/s²

        double a_t = calculateTangentialAcceleration(r, alpha);

        ASSERT_NEAR(a_t, 6.0, TOLERANCE);  // a_t = rα
        return true;
    });

    // ========================================
    // Angular Velocity and Direction
    // ========================================

    run_test("Angle turned in time", []() {
        double omega = 2.0;  // 2 rad/s
        double t = 5.0;      // 5 seconds

        double theta = calculateAngleTurned(omega, t);

        ASSERT_NEAR(theta, 10.0, TOLERANCE);  // θ = ωt
        return true;
    });

    run_test("Angular velocity change magnitude", []() {
        double omega = 10.0;
        double angle = M_PI / 3.0;  // 60 degrees

        double delta_omega = calculateAngularVelocityChange(omega, angle);

        // Δω = 2ω sin(θ/2)
        double expected = 2.0 * omega * std::sin(angle / 2.0);
        ASSERT_NEAR(delta_omega, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Torque
    // ========================================

    run_test("Torque calculation perpendicular force", []() {
        double r = 2.0;   // 2 meter lever arm
        double F = 50.0;  // 50 N force

        double tau = calculateTorque(r, F);

        ASSERT_NEAR(tau, 100.0, TOLERANCE);  // τ = rF
        return true;
    });

    run_test("Torque with angle", []() {
        double r = 2.0;
        double F = 50.0;
        double angle = M_PI / 2.0;  // 90 degrees

        double tau = calculateTorqueWithAngle(r, F, angle);

        ASSERT_NEAR(tau, 100.0, TOLERANCE);  // τ = rF sin(90°) = rF
        return true;
    });

    run_test("Torque zero for parallel force", []() {
        double r = 2.0;
        double F = 50.0;
        double angle = 0.0;  // Parallel

        double tau = calculateTorqueWithAngle(r, F, angle);

        ASSERT_NEAR(tau, 0.0, TOLERANCE);
        return true;
    });

    run_test("Torque proportional to lever arm", []() {
        double F = 50.0;
        double tau1 = calculateTorque(1.0, F);
        double tau2 = calculateTorque(2.0, F);

        ASSERT_NEAR(tau2, 2.0 * tau1, TOLERANCE);
        return true;
    });

    // ========================================
    // Angular Acceleration from Torque
    // ========================================

    run_test("Angular acceleration from torque", []() {
        double tau = 100.0;  // 100 N·m
        double I = 50.0;     // 50 kg·m²

        double alpha = calculateAngularAccelFromTorque(tau, I);

        ASSERT_NEAR(alpha, 2.0, TOLERANCE);  // α = τ/I
        return true;
    });

    run_test("Required torque for acceleration", []() {
        double I = 25.0;
        double alpha = 4.0;

        double tau = calculateRequiredTorque(I, alpha);

        ASSERT_NEAR(tau, 100.0, TOLERANCE);  // τ = Iα
        return true;
    });

    run_test("Larger moment needs more torque", []() {
        double alpha = 2.0;
        double tau1 = calculateRequiredTorque(10.0, alpha);
        double tau2 = calculateRequiredTorque(20.0, alpha);

        ASSERT_NEAR(tau2, 2.0 * tau1, TOLERANCE);
        return true;
    });

    // ========================================
    // Angular Momentum
    // ========================================

    run_test("Angular momentum calculation", []() {
        double I = 10.0;     // 10 kg·m²
        double omega = 5.0;  // 5 rad/s

        double L = calculateAngularMomentum(I, omega);

        ASSERT_NEAR(L, 50.0, TOLERANCE);  // L = Iω
        return true;
    });

    run_test("Point mass angular momentum", []() {
        double m = 2.0;    // 2 kg
        double v = 10.0;   // 10 m/s
        double r = 3.0;    // 3 m radius

        double L = calculatePointMassAngularMomentum(m, v, r);

        ASSERT_NEAR(L, 60.0, TOLERANCE);  // L = mvr
        return true;
    });

    run_test("Angular momentum change from torque", []() {
        double tau = 20.0;  // 20 N·m
        double t = 5.0;     // 5 seconds

        double delta_L = calculateAngularMomentumChange(tau, t);

        ASSERT_NEAR(delta_L, 100.0, TOLERANCE);  // ΔL = τΔt
        return true;
    });

    run_test("Angular momentum conservation", []() {
        double I = 10.0;
        double omega = 5.0;

        double L = calculateAngularMomentum(I, omega);

        // If I doubles and ω halves, L should remain constant
        double L2 = calculateAngularMomentum(2.0 * I, omega / 2.0);

        ASSERT_NEAR(L, L2, TOLERANCE);
        return true;
    });

    // ========================================
    // Rotational Kinetic Energy
    // ========================================

    run_test("Rotational kinetic energy", []() {
        double I = 10.0;
        double omega = 4.0;

        double KE = calculateRotationalKE(I, omega);

        ASSERT_NEAR(KE, 80.0, TOLERANCE);  // KE = ½Iω² = ½×10×16
        return true;
    });

    run_test("KE proportional to omega squared", []() {
        double I = 5.0;
        double KE1 = calculateRotationalKE(I, 2.0);
        double KE2 = calculateRotationalKE(I, 4.0);

        ASSERT_NEAR(KE2 / KE1, 4.0, TOLERANCE);  // ω doubles, KE quadruples
        return true;
    });

    run_test("Angular velocity from kinetic energy", []() {
        double KE = 100.0;
        double I = 10.0;

        double omega = calculateAngularVelocityFromKE(KE, I);

        // KE = ½Iω², so ω = √(2KE/I) = √20 ≈ 4.47
        ASSERT_NEAR(omega, std::sqrt(20.0), TOLERANCE);
        return true;
    });

    // ========================================
    // Moment of Inertia Formulas
    // ========================================

    run_test("Point mass moment of inertia", []() {
        double m = 5.0;
        double r = 2.0;

        double I = momentOfInertiaPointMass(m, r);

        ASSERT_NEAR(I, 20.0, TOLERANCE);  // I = mr²
        return true;
    });

    run_test("Rod about center", []() {
        double m = 12.0;
        double L = 2.0;

        double I = momentOfInertiaRodCenter(m, L);

        ASSERT_NEAR(I, 4.0, TOLERANCE);  // I = (1/12)ML² = 12×4/12 = 4
        return true;
    });

    run_test("Rod about end", []() {
        double m = 12.0;
        double L = 2.0;

        double I = momentOfInertiaRodEnd(m, L);

        ASSERT_NEAR(I, 16.0, TOLERANCE);  // I = (1/3)ML² = 12×4/3 = 16
        return true;
    });

    run_test("Rod end is three times center", []() {
        double m = 12.0;
        double L = 2.0;

        double I_center = momentOfInertiaRodCenter(m, L);
        double I_end = momentOfInertiaRodEnd(m, L);

        ASSERT_NEAR(I_end, 4.0 * I_center, TOLERANCE);
        return true;
    });

    run_test("Solid cylinder about axis", []() {
        double m = 10.0;
        double r = 2.0;

        double I = momentOfInertiaCylinder(m, r);

        ASSERT_NEAR(I, 20.0, TOLERANCE);  // I = ½MR² = ½×10×4
        return true;
    });

    run_test("Hollow cylinder about axis", []() {
        double m = 10.0;
        double r = 2.0;

        double I = momentOfInertiaHollowCylinder(m, r);

        ASSERT_NEAR(I, 40.0, TOLERANCE);  // I = MR² = 10×4
        return true;
    });

    run_test("Hollow cylinder twice solid", []() {
        double m = 10.0;
        double r = 2.0;

        double I_solid = momentOfInertiaCylinder(m, r);
        double I_hollow = momentOfInertiaHollowCylinder(m, r);

        ASSERT_NEAR(I_hollow, 2.0 * I_solid, TOLERANCE);
        return true;
    });

    run_test("Solid sphere about diameter", []() {
        double m = 10.0;
        double r = 2.0;

        double I = momentOfInertiaSolidSphere(m, r);

        ASSERT_NEAR(I, 16.0, TOLERANCE);  // I = (2/5)MR² = 0.4×10×4
        return true;
    });

    run_test("Hollow sphere about diameter", []() {
        double m = 10.0;
        double r = 2.0;

        double I = momentOfInertiaHollowSphere(m, r);

        ASSERT_NEAR(I, 26.667, LOOSE_TOLERANCE);  // I = (2/3)MR² = (2/3)×10×4
        return true;
    });

    run_test("Hollow sphere larger than solid", []() {
        double m = 10.0;
        double r = 2.0;

        double I_solid = momentOfInertiaSolidSphere(m, r);
        double I_hollow = momentOfInertiaHollowSphere(m, r);

        ASSERT_TRUE(I_hollow > I_solid);
        return true;
    });

    run_test("Disk about central axis", []() {
        double m = 10.0;
        double r = 2.0;

        double I = momentOfInertiaDisk(m, r);

        ASSERT_NEAR(I, 20.0, TOLERANCE);  // I = ½MR²
        return true;
    });

    run_test("Hoop about central axis", []() {
        double m = 10.0;
        double r = 2.0;

        double I = momentOfInertiaHoop(m, r);

        ASSERT_NEAR(I, 40.0, TOLERANCE);  // I = MR²
        return true;
    });

    run_test("Hoop twice disk moment", []() {
        double m = 10.0;
        double r = 2.0;

        double I_disk = momentOfInertiaDisk(m, r);
        double I_hoop = momentOfInertiaHoop(m, r);

        ASSERT_NEAR(I_hoop, 2.0 * I_disk, TOLERANCE);
        return true;
    });

    // ========================================
    // Parallel Axis Theorem
    // ========================================

    run_test("Parallel axis theorem", []() {
        double I_cm = 10.0;  // Moment about CM
        double m = 5.0;
        double d = 2.0;      // Distance from CM

        double I = parallelAxisTheorem(I_cm, m, d);

        ASSERT_NEAR(I, 30.0, TOLERANCE);  // I = I_cm + md² = 10 + 5×4
        return true;
    });

    run_test("Parallel axis theorem rod center to end", []() {
        double m = 12.0;
        double L = 2.0;

        double I_center = momentOfInertiaRodCenter(m, L);
        double I_end = parallelAxisTheorem(I_center, m, L / 2.0);

        double I_end_direct = momentOfInertiaRodEnd(m, L);

        ASSERT_NEAR(I_end, I_end_direct, TOLERANCE);
        return true;
    });

    run_test("Parallel axis increases moment of inertia", []() {
        double I_cm = 10.0;
        double m = 5.0;
        double d = 1.0;

        double I = parallelAxisTheorem(I_cm, m, d);

        ASSERT_TRUE(I > I_cm);
        return true;
    });

    // ========================================
    // Compound Pendulum
    // ========================================

    run_test("Compound pendulum period", []() {
        double I = 2.0;   // Moment about pivot
        double m = 1.0;
        double d = 0.5;   // Distance to CM

        double T = calculateCompoundPendulumPeriod(I, m, d);

        ASSERT_TRUE(T > 0);
        return true;
    });

    run_test("Compound pendulum equivalent length", []() {
        double I = 2.0;
        double m = 1.0;
        double d = 0.5;

        double L_eq = calculateEquivalentLength(I, m, d);

        // L_eq = I/(md) = 2/(1×0.5) = 4
        ASSERT_NEAR(L_eq, 4.0, TOLERANCE);
        return true;
    });

    run_test("Uniform rod compound pendulum", []() {
        double m = 1.0;
        double L = 1.0;

        // Rod pivoted at end
        double I = momentOfInertiaRodEnd(m, L);
        double d = L / 2.0;  // CM at center

        double T = calculateCompoundPendulumPeriod(I, m, d);

        ASSERT_TRUE(T > 0);
        return true;
    });

    // ========================================
    // Center of Percussion
    // ========================================

    run_test("Center of percussion calculation", []() {
        double I = 2.0;
        double m = 1.0;
        double d = 0.5;

        double d_p = calculateCenterOfPercussion(I, m, d);

        ASSERT_NEAR(d_p, 4.0, TOLERANCE);  // d_p = I/(md)
        return true;
    });

    run_test("Center of percussion equals equivalent length", []() {
        double I = 2.0;
        double m = 1.0;
        double d = 0.5;

        double d_p = calculateCenterOfPercussion(I, m, d);
        double L_eq = calculateEquivalentLength(I, m, d);

        ASSERT_NEAR(d_p, L_eq, TOLERANCE);
        return true;
    });

    run_test("Rod center of percussion", []() {
        double L = 1.0;

        double d_p = calculateRodCenterOfPercussion(L);

        ASSERT_NEAR(d_p, 2.0 / 3.0, TOLERANCE);  // 2L/3
        return true;
    });

    run_test("Rod center of percussion at two thirds length", []() {
        double m = 1.0;
        double L = 1.0;

        // Rod pivoted at end
        double I = momentOfInertiaRodEnd(m, L);
        double d_cm = L / 2.0;

        double d_p = calculateCenterOfPercussion(I, m, d_cm);
        double d_p_formula = calculateRodCenterOfPercussion(L);

        ASSERT_NEAR(d_p, d_p_formula, TOLERANCE);
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
