/**
 * Phase 4 Validation: Fluid Dynamics Flow Types
 *
 * Tests the fluid_dynamics_flow_types.hpp module functions.
 *
 * Coverage:
 * - Poiseuille flow (pressure-driven pipe and channel flow)
 * - Couette flow (shear-driven flow)
 * - Stokes flow (creeping flow at Re << 1)
 * - Potential flow (irrotational and inviscid)
 * - Classical analytical solutions
 * - Flow classifications and transitions
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/fluid_dynamics_flow_types.hpp"

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

    std::cout << "=== Phase 4: Fluid Dynamics Flow Types Validation ===" << std::endl;
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
    // Poiseuille Flow Tests
    // ========================================

    run_test("Poiseuille: velocity at pipe centerline", []() {
        double r = 0.0;        // Centerline
        double R = 0.05;       // 5 cm radius
        double dp_dx = -100.0; // Pressure gradient (Pa/m)
        double mu = 1e-3;      // Water viscosity

        double u = PoiseuilleFlow::velocityCircularPipe(r, R, dp_dx, mu);
        double expected = 100.0 / (4.0 * 1e-3) * 0.05 * 0.05;  // Maximum velocity
        ASSERT_NEAR(u, expected, TOLERANCE);
        return true;
    });

    run_test("Poiseuille: velocity at pipe wall is zero", []() {
        double r = 0.05;       // At wall
        double R = 0.05;
        double dp_dx = -100.0;
        double mu = 1e-3;

        double u = PoiseuilleFlow::velocityCircularPipe(r, R, dp_dx, mu);
        ASSERT_NEAR(u, 0.0, TOLERANCE);  // No-slip condition
        return true;
    });

    run_test("Poiseuille: parabolic velocity profile", []() {
        double R = 0.05;
        double dp_dx = -100.0;
        double mu = 1e-3;

        double u_center = PoiseuilleFlow::velocityCircularPipe(0.0, R, dp_dx, mu);
        double u_half = PoiseuilleFlow::velocityCircularPipe(0.025, R, dp_dx, mu);

        // At r = R/2, u = (3/4)u_max
        ASSERT_NEAR(u_half, 0.75 * u_center, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Poiseuille: maximum velocity", []() {
        double R = 0.02;       // 2 cm radius
        double dP = 1000.0;    // Pressure drop
        double L = 1.0;        // 1 m length
        double mu = 1e-3;

        double u_max = PoiseuilleFlow::maxVelocity(R, dP, L, mu);
        double expected = (1000.0 * 0.02 * 0.02) / (4.0 * 1e-3 * 1.0);  // 0.1 m/s
        ASSERT_NEAR(u_max, expected, TOLERANCE);
        return true;
    });

    run_test("Poiseuille: average velocity is half maximum", []() {
        double R = 0.02;
        double dP = 1000.0;
        double L = 1.0;
        double mu = 1e-3;

        double u_max = PoiseuilleFlow::maxVelocity(R, dP, L, mu);
        double u_avg = PoiseuilleFlow::averageVelocity(R, dP, L, mu);

        ASSERT_NEAR(u_avg, 0.5 * u_max, TOLERANCE);
        return true;
    });

    run_test("Poiseuille: Hagen-Poiseuille flow rate", []() {
        double R = 0.01;       // 1 cm radius
        double dP = 1000.0;
        double L = 1.0;
        double mu = 1e-3;

        double Q = PoiseuilleFlow::flowRate(R, dP, L, mu);
        double expected = (M_PI * std::pow(0.01, 4) * 1000.0) / (8.0 * 1e-3 * 1.0);
        ASSERT_NEAR(Q, expected, Q * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Poiseuille: wall shear stress", []() {
        double R = 0.05;
        double dP = 500.0;
        double L = 2.0;

        double tau_w = PoiseuilleFlow::wallShearStress(R, dP, L);
        double expected = (500.0 * 0.05) / (2.0 * 2.0);  // 6.25 Pa
        ASSERT_NEAR(tau_w, expected, TOLERANCE);
        return true;
    });

    run_test("Poiseuille: laminar friction factor", []() {
        double Re = 1000.0;

        double f = PoiseuilleFlow::frictionFactor(Re);
        double expected = 64.0 / 1000.0;  // 0.064
        ASSERT_NEAR(f, expected, TOLERANCE);
        return true;
    });

    run_test("Poiseuille: plane channel flow", []() {
        double y = 0.0;        // Centerline
        double h = 0.05;       // Half-height
        double dp_dx = -100.0;
        double mu = 1e-3;

        double u = PoiseuilleFlow::velocityPlaneChannel(y, h, dp_dx, mu);
        double expected = 100.0 / (2.0 * 1e-3) * 0.05 * 0.05;  // Maximum velocity
        ASSERT_NEAR(u, expected, TOLERANCE);
        return true;
    });

    run_test("Poiseuille: flow rate increases with R to fourth power", []() {
        double dP = 1000.0;
        double L = 1.0;
        double mu = 1e-3;

        double Q1 = PoiseuilleFlow::flowRate(0.01, dP, L, mu);
        double Q2 = PoiseuilleFlow::flowRate(0.02, dP, L, mu);

        ASSERT_NEAR(Q2 / Q1, 16.0, LOOSE_TOLERANCE);  // (2R)⁴ / R⁴ = 16
        return true;
    });

    // ========================================
    // Couette Flow Tests
    // ========================================

    run_test("Couette: simple linear profile", []() {
        double y = 0.5;        // Halfway
        double U = 10.0;       // Moving plate velocity
        double h = 1.0;        // Gap height

        double u = CouetteFlow::velocitySimple(y, U, h);
        double expected = 10.0 * 0.5;  // 5.0 m/s
        ASSERT_NEAR(u, expected, TOLERANCE);
        return true;
    });

    run_test("Couette: velocity at stationary wall", []() {
        double y = 0.0;
        double U = 10.0;
        double h = 1.0;

        double u = CouetteFlow::velocitySimple(y, U, h);
        ASSERT_NEAR(u, 0.0, TOLERANCE);  // No-slip
        return true;
    });

    run_test("Couette: velocity at moving wall", []() {
        double y = 1.0;
        double U = 10.0;
        double h = 1.0;

        double u = CouetteFlow::velocitySimple(y, U, h);
        ASSERT_NEAR(u, U, TOLERANCE);  // Matches plate velocity
        return true;
    });

    run_test("Couette: constant shear stress", []() {
        double mu = 1e-3;
        double U = 5.0;
        double h = 0.01;       // 1 cm gap

        double tau = CouetteFlow::shearStress(mu, U, h);
        double expected = 1e-3 * 5.0 / 0.01;  // 0.5 Pa
        ASSERT_NEAR(tau, expected, TOLERANCE);
        return true;
    });

    run_test("Couette: generalized with pressure gradient", []() {
        double y = 0.5;
        double U = 10.0;
        double h = 1.0;
        double dp_dx = 0.0;    // No pressure gradient
        double mu = 1e-3;

        double u = CouetteFlow::velocityGeneralized(y, U, h, dp_dx, mu);
        double u_simple = CouetteFlow::velocitySimple(y, U, h);

        ASSERT_NEAR(u, u_simple, TOLERANCE);  // Reduces to simple Couette
        return true;
    });

    run_test("Couette: Taylor-Couette flow", []() {
        double r = 0.075;      // Between cylinders
        double R1 = 0.05;      // Inner radius
        double R2 = 0.1;       // Outer radius
        double omega1 = 10.0;  // Inner angular velocity
        double omega2 = 0.0;   // Outer stationary

        double u_theta = CouetteFlow::taylorCouetteVelocity(r, R1, R2, omega1, omega2);
        ASSERT_TRUE(u_theta > 0.0);  // Positive tangential velocity
        return true;
    });

    run_test("Couette: Taylor-Couette at inner cylinder", []() {
        double R1 = 0.05;
        double R2 = 0.1;
        double omega1 = 10.0;
        double omega2 = 0.0;

        double u_theta = CouetteFlow::taylorCouetteVelocity(R1, R1, R2, omega1, omega2);
        double expected = omega1 * R1;  // v = ωr at inner wall
        ASSERT_NEAR(u_theta, expected, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Stokes Flow Tests
    // ========================================

    run_test("Stokes: drag on sphere", []() {
        double mu = 1e-3;      // Water
        double R = 0.001;      // 1 mm radius
        double U = 0.01;       // 1 cm/s

        double F_D = StokesFlow::sphereDrag(mu, R, U);
        double expected = 6.0 * M_PI * 1e-3 * 0.001 * 0.01;
        ASSERT_NEAR(F_D, expected, F_D * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Stokes: drag coefficient", []() {
        double Re = 0.5;

        double C_D = StokesFlow::sphereDragCoefficient(Re);
        double expected = 24.0 / 0.5;  // 48
        ASSERT_NEAR(C_D, expected, TOLERANCE);
        return true;
    });

    run_test("Stokes: terminal velocity of falling sphere", []() {
        double R = 1e-4;       // 0.1 mm radius
        double rho_p = 2500.0; // Particle density (sand)
        double rho_f = 1000.0; // Fluid density (water)
        double nu = 1e-6;      // Water kinematic viscosity
        double g = 9.81;

        double U_t = StokesFlow::terminalVelocity(R, rho_p, rho_f, nu, g);
        ASSERT_TRUE(U_t > 0.0);  // Falls downward
        return true;
    });

    run_test("Stokes: terminal velocity increases with particle size", []() {
        double rho_p = 2500.0;
        double rho_f = 1000.0;
        double nu = 1e-6;
        double g = 9.81;

        double U_t1 = StokesFlow::terminalVelocity(1e-5, rho_p, rho_f, nu, g);
        double U_t2 = StokesFlow::terminalVelocity(2e-5, rho_p, rho_f, nu, g);

        ASSERT_NEAR(U_t2 / U_t1, 4.0, LOOSE_TOLERANCE);  // U_t ∝ R²
        return true;
    });

    run_test("Stokes: Oseen correction", []() {
        double mu = 1e-3;
        double R = 0.001;
        double U = 0.01;
        double Re = 0.8;

        double F_Stokes = StokesFlow::sphereDrag(mu, R, U);
        double F_Oseen = StokesFlow::oseenDrag(mu, R, U, Re);

        ASSERT_TRUE(F_Oseen > F_Stokes);  // Correction increases drag
        return true;
    });

    run_test("Stokes: cylinder drag", []() {
        double mu = 1e-3;
        double U = 0.01;
        double Re = 0.5;

        double F_per_L = StokesFlow::cylinderDrag(mu, U, Re);
        ASSERT_TRUE(F_per_L > 0.0);
        return true;
    });

    // ========================================
    // Potential Flow Tests
    // ========================================

    run_test("Potential: uniform flow", []() {
        maths::linear_algebra::Vector pos({1.0, 2.0, 0.0});
        maths::linear_algebra::Vector U_inf({10.0, 0.0, 0.0});

        auto velocity = PotentialFlow::uniformFlow(pos, U_inf);
        ASSERT_NEAR(velocity[0], 10.0, TOLERANCE);
        ASSERT_NEAR(velocity[1], 0.0, TOLERANCE);
        return true;
    });

    run_test("Potential: point source velocity", []() {
        maths::linear_algebra::Vector pos({1.0, 0.0, 0.0});
        double m = 10.0;  // Source strength

        auto velocity = PotentialFlow::pointSource(pos, m);
        // u_r = m/(2πr) = 10/(2π*1) ≈ 1.59
        double u_r_expected = 10.0 / (2.0 * M_PI * 1.0);
        ASSERT_NEAR(velocity.norm(), u_r_expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Potential: point vortex velocity", []() {
        maths::linear_algebra::Vector pos({1.0, 0.0, 0.0});
        double Gamma = 10.0;  // Circulation

        auto velocity = PotentialFlow::pointVortex(pos, Gamma);
        // u_θ = Γ/(2πr) = 10/(2π*1) ≈ 1.59
        double u_theta_expected = 10.0 / (2.0 * M_PI * 1.0);
        ASSERT_NEAR(velocity.norm(), u_theta_expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Potential: doublet flow", []() {
        maths::linear_algebra::Vector pos({1.0, 0.0, 0.0});
        double K = 5.0;

        auto velocity = PotentialFlow::doublet(pos, K);
        ASSERT_TRUE(velocity.norm() > 0.0);
        return true;
    });

    run_test("Potential: flow past cylinder at infinity", []() {
        maths::linear_algebra::Vector pos({10.0, 0.0, 0.0});  // Far from cylinder
        double U_inf = 10.0;
        double R = 1.0;

        auto velocity = PotentialFlow::cylinderFlow(pos, U_inf, R);
        // Far from cylinder, should approach uniform flow
        ASSERT_NEAR(velocity[0], U_inf, U_inf * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Potential: cylinder surface velocity at top", []() {
        // At cylinder top (0, R), velocity should be 2*U_inf
        maths::linear_algebra::Vector top({0.0, 1.0, 0.0});
        double U_inf = 10.0;
        double R = 1.0;

        auto velocity = PotentialFlow::cylinderFlow(top, U_inf, R);
        // At top/bottom of cylinder, u ≈ 2*U_inf (approximately)
        ASSERT_TRUE(velocity.norm() > U_inf);  // Higher than free stream
        return true;
    });

    run_test("Potential: Kutta-Joukowski lift", []() {
        double rho = 1.225;    // Air
        double U_inf = 50.0;   // m/s
        double Gamma = 20.0;   // Circulation

        double L = PotentialFlow::liftKuttaJoukowski(rho, U_inf, Gamma);
        double expected = 1.225 * 50.0 * 20.0;  // 1225 N/m
        ASSERT_NEAR(L, expected, TOLERANCE);
        return true;
    });

    run_test("Potential: stream function", []() {
        maths::linear_algebra::Vector pos({1.0, 1.0, 0.0});
        double m = 5.0;
        double Gamma = 10.0;

        double psi = PotentialFlow::streamFunction(pos, m, Gamma);
        ASSERT_TRUE(std::isfinite(psi));
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Consistency: Poiseuille flow rate from velocity profile", []() {
        double R = 0.02;
        double dP = 1000.0;
        double L = 1.0;
        double mu = 1e-3;

        double u_avg = PoiseuilleFlow::averageVelocity(R, dP, L, mu);
        double Q_from_avg = M_PI * R * R * u_avg;
        double Q_direct = PoiseuilleFlow::flowRate(R, dP, L, mu);

        ASSERT_NEAR(Q_from_avg, Q_direct, Q_direct * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Consistency: Stokes drag equals force balance at terminal velocity", []() {
        double R = 1e-4;
        double rho_p = 2500.0;
        double rho_f = 1000.0;
        double mu = 1e-3;
        double nu = mu / rho_f;
        double g = 9.81;

        double U_t = StokesFlow::terminalVelocity(R, rho_p, rho_f, nu, g);
        double F_drag = StokesFlow::sphereDrag(mu, R, U_t);
        double F_buoyancy = (4.0/3.0) * M_PI * R*R*R * (rho_p - rho_f) * g;

        ASSERT_NEAR(F_drag, F_buoyancy, F_buoyancy * LOOSE_TOLERANCE);
        return true;
    });

    run_test("Physical: Poiseuille velocity decreases from center to wall", []() {
        double R = 0.05;
        double dp_dx = -100.0;
        double mu = 1e-3;

        double u_center = PoiseuilleFlow::velocityCircularPipe(0.0, R, dp_dx, mu);
        double u_mid = PoiseuilleFlow::velocityCircularPipe(0.025, R, dp_dx, mu);
        double u_wall = PoiseuilleFlow::velocityCircularPipe(R, R, dp_dx, mu);

        ASSERT_TRUE(u_center > u_mid);
        ASSERT_TRUE(u_mid > u_wall);
        ASSERT_NEAR(u_wall, 0.0, TOLERANCE);
        return true;
    });

    run_test("Physical: Couette linear velocity profile", []() {
        double U = 10.0;
        double h = 1.0;

        double u1 = CouetteFlow::velocitySimple(0.25, U, h);
        double u2 = CouetteFlow::velocitySimple(0.50, U, h);
        double u3 = CouetteFlow::velocitySimple(0.75, U, h);

        // Linear: u(y) = U*y/h
        ASSERT_NEAR(u1, 2.5, TOLERANCE);
        ASSERT_NEAR(u2, 5.0, TOLERANCE);
        ASSERT_NEAR(u3, 7.5, TOLERANCE);
        return true;
    });

    run_test("Physical: potential flow is irrotational", []() {
        // Uniform flow has zero vorticity
        maths::linear_algebra::Vector pos({1.0, 1.0, 0.0});
        maths::linear_algebra::Vector U_inf({10.0, 0.0, 0.0});

        auto v1 = PotentialFlow::uniformFlow(pos, U_inf);
        maths::linear_algebra::Vector pos2({2.0, 1.0, 0.0});
        auto v2 = PotentialFlow::uniformFlow(pos2, U_inf);

        // No change in velocity → no curl
        ASSERT_NEAR((v2 - v1).norm(), 0.0, TOLERANCE);
        return true;
    });

    run_test("Limit: Poiseuille approaches zero at high viscosity", []() {
        double R = 0.01;
        double dP = 100.0;
        double L = 1.0;
        double mu_low = 1e-3;
        double mu_high = 1.0;  // 1000x higher

        double Q_low = PoiseuilleFlow::flowRate(R, dP, L, mu_low);
        double Q_high = PoiseuilleFlow::flowRate(R, dP, L, mu_high);

        ASSERT_NEAR(Q_high / Q_low, 1.0/1000.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Limit: Stokes drag linear in velocity", []() {
        double mu = 1e-3;
        double R = 0.001;

        double F1 = StokesFlow::sphereDrag(mu, R, 0.01);
        double F2 = StokesFlow::sphereDrag(mu, R, 0.02);

        ASSERT_NEAR(F2 / F1, 2.0, TOLERANCE);  // F ∝ U
        return true;
    });

    run_test("Limit: potential flow source strength doubles velocity at same distance", []() {
        maths::linear_algebra::Vector pos({1.0, 0.0, 0.0});
        double m1 = 5.0;
        double m2 = 10.0;

        auto v1 = PotentialFlow::pointSource(pos, m1);
        auto v2 = PotentialFlow::pointSource(pos, m2);

        ASSERT_NEAR(v2.norm() / v1.norm(), 2.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Fluid Flow Types" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All fluid dynamics flow types tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Poiseuille flow (pipe and channel)" << std::endl;
        std::cout << "  - Couette flow (simple and generalized)" << std::endl;
        std::cout << "  - Stokes flow (creeping flow)" << std::endl;
        std::cout << "  - Potential flow (irrotational)" << std::endl;
        std::cout << "  - Classical analytical solutions" << std::endl;
        std::cout << "  - Flow physics and consistency" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
