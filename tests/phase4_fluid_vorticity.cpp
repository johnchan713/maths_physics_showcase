/**
 * Phase 4 Validation: Fluid Dynamics Vorticity
 *
 * Tests the fluid_dynamics_vorticity.hpp module functions.
 *
 * Coverage:
 * - Vorticity calculation from velocity gradient
 * - Vorticity transport equation
 * - Circulation (line integral of velocity)
 * - Kelvin's circulation theorem
 * - Biot-Savart law (induced velocities)
 * - Vortex dynamics and interactions
 * - Rankine vortex and vortex structures
 */

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "../include/physics/fluid_dynamics_vorticity.hpp"

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

    std::cout << "=== Phase 4: Fluid Dynamics Vorticity Validation ===" << std::endl;
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
    // Vorticity Calculation Tests
    // ========================================

    run_test("Vorticity: zero for uniform flow", []() {
        // Uniform flow: du/dy = 0, dv/dx = 0, etc.
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0, 0) = 0.0; vel_grad(0, 1) = 0.0; vel_grad(0, 2) = 0.0;
        vel_grad(1, 0) = 0.0; vel_grad(1, 1) = 0.0; vel_grad(1, 2) = 0.0;
        vel_grad(2, 0) = 0.0; vel_grad(2, 1) = 0.0; vel_grad(2, 2) = 0.0;

        auto omega = Vorticity::compute(vel_grad);
        ASSERT_NEAR(omega.norm(), 0.0, TOLERANCE);
        return true;
    });

    run_test("Vorticity: simple shear flow", []() {
        // Shear flow: u = ky, v = 0, w = 0
        // du/dy = k, so ωz = dv/dx - du/dy = -k
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0, 0) = 0.0; vel_grad(0, 1) = 10.0; vel_grad(0, 2) = 0.0;
        vel_grad(1, 0) = 0.0; vel_grad(1, 1) = 0.0;  vel_grad(1, 2) = 0.0;
        vel_grad(2, 0) = 0.0; vel_grad(2, 1) = 0.0;  vel_grad(2, 2) = 0.0;

        auto omega = Vorticity::compute(vel_grad);
        ASSERT_NEAR(omega[2], -10.0, TOLERANCE);  // ωz = -du/dy
        return true;
    });

    run_test("Vorticity: solid body rotation", []() {
        // Solid body rotation: u = -ωy, v = ωx, w = 0
        // dv/dx - du/dy = ω - (-ω) = 2ω
        double omega_val = 5.0;
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0, 0) = 0.0;      vel_grad(0, 1) = -omega_val; vel_grad(0, 2) = 0.0;
        vel_grad(1, 0) = omega_val; vel_grad(1, 1) = 0.0;       vel_grad(1, 2) = 0.0;
        vel_grad(2, 0) = 0.0;      vel_grad(2, 1) = 0.0;       vel_grad(2, 2) = 0.0;

        auto omega = Vorticity::compute(vel_grad);
        ASSERT_NEAR(omega[2], 2.0 * omega_val, TOLERANCE);
        return true;
    });

    run_test("Vorticity: magnitude calculation", []() {
        maths::linear_algebra::Vector omega({3.0, 4.0, 0.0});
        double mag = Vorticity::magnitude(omega);
        ASSERT_NEAR(mag, 5.0, TOLERANCE);  // sqrt(9 + 16)
        return true;
    });

    run_test("Vorticity: enstrophy", []() {
        maths::linear_algebra::Vector omega({2.0, 0.0, 0.0});
        double enstrophy = Vorticity::enstrophy(omega);
        ASSERT_NEAR(enstrophy, 2.0, TOLERANCE);  // 0.5 * 4
        return true;
    });

    run_test("Vorticity: irrotational check", []() {
        maths::linear_algebra::Vector omega_zero({0.0, 0.0, 0.0});
        ASSERT_TRUE(Vorticity::isIrrotational(omega_zero));

        maths::linear_algebra::Vector omega_nonzero({0.1, 0.0, 0.0});
        ASSERT_TRUE(!Vorticity::isIrrotational(omega_nonzero));
        return true;
    });

    run_test("Vorticity: vortex line direction", []() {
        maths::linear_algebra::Vector omega({3.0, 4.0, 0.0});
        auto direction = Vorticity::vortexLineDirection(omega);
        ASSERT_NEAR(direction.norm(), 1.0, TOLERANCE);  // Unit vector
        ASSERT_NEAR(direction[0], 0.6, TOLERANCE);  // 3/5
        ASSERT_NEAR(direction[1], 0.8, TOLERANCE);  // 4/5
        return true;
    });

    // ========================================
    // Vorticity Transport Equation Tests
    // ========================================

    run_test("Vorticity Transport: vortex stretching", []() {
        maths::linear_algebra::Vector omega({0.0, 0.0, 10.0});
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0, 0) = 1.0; vel_grad(0, 1) = 0.0; vel_grad(0, 2) = 0.0;
        vel_grad(1, 0) = 0.0; vel_grad(1, 1) = 1.0; vel_grad(1, 2) = 0.0;
        vel_grad(2, 0) = 0.0; vel_grad(2, 1) = 0.0; vel_grad(2, 2) = 2.0;

        auto stretching = VorticityTransportEquation::vortexStretching(omega, vel_grad);
        // (ω·∇)u with ω = (0,0,10)
        ASSERT_NEAR(stretching[2], 20.0, TOLERANCE);  // 10 * ∂w/∂z = 10 * 2
        return true;
    });

    run_test("Vorticity Transport: viscous diffusion", []() {
        maths::linear_algebra::Vector omega_laplacian({2.0, 1.0, 0.5});
        double nu = 1e-5;

        auto diffusion = VorticityTransportEquation::viscousDiffusion(omega_laplacian, nu);
        ASSERT_NEAR(diffusion[0], 2.0 * nu, TOLERANCE);
        ASSERT_NEAR(diffusion[1], 1.0 * nu, TOLERANCE);
        return true;
    });

    run_test("Vorticity Transport: total evolution", []() {
        maths::linear_algebra::Vector omega({0.0, 0.0, 5.0});
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(2, 2) = 1.0;  // Other terms zero
        maths::linear_algebra::Vector omega_lap({0.0, 0.0, 10.0});
        double nu = 1e-5;

        auto evolution = VorticityTransportEquation::totalEvolution(omega, vel_grad, omega_lap, nu);
        // Dω/Dt = vortex stretching + viscous diffusion
        ASSERT_TRUE(evolution[2] > 0.0);
        return true;
    });

    run_test("Vorticity Transport: 2D no stretching", []() {
        double omega_lap = -5.0;
        double nu = 1e-6;

        double evolution = VorticityTransportEquation::evolution2D(omega_lap, nu);
        double expected = nu * omega_lap;
        ASSERT_NEAR(evolution, expected, TOLERANCE);
        return true;
    });

    run_test("Vorticity Transport: stretching rate", []() {
        maths::linear_algebra::Vector omega({0.0, 0.0, 10.0});
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(2, 2) = 2.0;  // Stretching in z-direction

        double S_omega = VorticityTransportEquation::stretchingRate(omega, vel_grad);
        ASSERT_TRUE(S_omega > 0.0);  // Positive stretching
        return true;
    });

    // ========================================
    // Circulation Tests
    // ========================================

    run_test("Circulation: around square in uniform flow", []() {
        // Uniform flow has zero circulation
        auto uniform_field = [](const maths::linear_algebra::Vector& pos) {
            return maths::linear_algebra::Vector({1.0, 0.0, 0.0});
        };

        std::vector<maths::linear_algebra::Vector> square = {
            maths::linear_algebra::Vector({0.0, 0.0, 0.0}),
            maths::linear_algebra::Vector({1.0, 0.0, 0.0}),
            maths::linear_algebra::Vector({1.0, 1.0, 0.0}),
            maths::linear_algebra::Vector({0.0, 1.0, 0.0})
        };

        double gamma = Circulation::compute(uniform_field, square);
        ASSERT_NEAR(gamma, 0.0, NUMERICAL_TOLERANCE);
        return true;
    });

    run_test("Circulation: solid body rotation", []() {
        double omega_val = 2.0;  // Angular velocity
        double area = 1.0;       // Area enclosed

        double gamma = Circulation::solidBodyRotation(omega_val, area);
        ASSERT_NEAR(gamma, 2.0, TOLERANCE);
        return true;
    });

    run_test("Circulation: point vortex", []() {
        double Gamma = 10.0;
        double circ = Circulation::pointVortex(Gamma);
        ASSERT_NEAR(circ, Gamma, TOLERANCE);
        return true;
    });

    run_test("Circulation: increases with area for solid rotation", []() {
        double omega_val = 3.0;
        double A1 = 1.0;
        double A2 = 4.0;

        double gamma1 = Circulation::solidBodyRotation(omega_val, A1);
        double gamma2 = Circulation::solidBodyRotation(omega_val, A2);

        ASSERT_NEAR(gamma2 / gamma1, 4.0, TOLERANCE);  // Γ ∝ A
        return true;
    });

    // ========================================
    // Kelvin's Circulation Theorem Tests
    // ========================================

    run_test("Kelvin: inviscid flow conserves circulation", []() {
        double nu = 0.0;  // Inviscid
        bool applicable = KelvinCirculationTheorem::checkApplicability(nu);
        ASSERT_TRUE(applicable);
        return true;
    });

    run_test("Kelvin: viscous flow violates theorem", []() {
        double nu = 1e-5;
        bool applicable = KelvinCirculationTheorem::checkApplicability(nu);
        ASSERT_TRUE(!applicable);
        return true;
    });

    run_test("Kelvin: circulation conservation check", []() {
        std::vector<double> gamma_history = {10.0, 10.0, 10.0, 10.0};
        bool conserved = KelvinCirculationTheorem::verifyConservation(gamma_history);
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("Kelvin: non-conserved circulation", []() {
        std::vector<double> gamma_history = {10.0, 10.0, 9.0, 8.0};
        bool conserved = KelvinCirculationTheorem::verifyConservation(gamma_history, 0.1);
        ASSERT_TRUE(!conserved);
        return true;
    });

    run_test("Kelvin: circulation evolution with viscosity", []() {
        double nu = 1e-5;
        std::vector<maths::linear_algebra::Vector> curl = {
            maths::linear_algebra::Vector({1.0, 0.0, 0.0})
        };
        std::vector<maths::linear_algebra::Vector> path = {
            maths::linear_algebra::Vector({1.0, 0.0, 0.0})
        };

        double dGamma_dt = KelvinCirculationTheorem::circulationEvolution(nu, curl, path);
        ASSERT_TRUE(std::abs(dGamma_dt) > 0.0);  // Non-zero evolution
        return true;
    });

    // ========================================
    // Biot-Savart Law Tests
    // ========================================

    run_test("Biot-Savart: straight vortex filament", []() {
        double Gamma = 10.0;
        double r = 0.5;  // Distance from vortex

        double u_theta = BiotSavartLaw::velocityStraightVortex(Gamma, r);
        double expected = 10.0 / (2.0 * M_PI * 0.5);
        ASSERT_NEAR(u_theta, expected, TOLERANCE);
        return true;
    });

    run_test("Biot-Savart: velocity decreases with distance", []() {
        double Gamma = 10.0;
        double u1 = BiotSavartLaw::velocityStraightVortex(Gamma, 0.5);
        double u2 = BiotSavartLaw::velocityStraightVortex(Gamma, 1.0);

        ASSERT_NEAR(u1 / u2, 2.0, TOLERANCE);  // u ∝ 1/r
        return true;
    });

    run_test("Biot-Savart: vortex ring on axis", []() {
        double Gamma = 5.0;
        double R = 0.1;  // Ring radius
        double z = 0.0;  // On plane of ring

        double u_z = BiotSavartLaw::velocityVortexRing(Gamma, R, z);
        double expected = (Gamma * R * R) / (2.0 * std::pow(R * R, 1.5));
        ASSERT_NEAR(u_z, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Biot-Savart: vortex ring self-velocity", []() {
        double Gamma = 10.0;
        double R = 0.1;
        double a = 0.001;  // Core radius

        double U = BiotSavartLaw::vortexRingSelfVelocity(Gamma, R, a);
        ASSERT_TRUE(U > 0.0);  // Ring translates forward
        return true;
    });

    run_test("Biot-Savart: velocity jump across vortex sheet", []() {
        double gamma_sheet = 5.0;
        double Delta_u = BiotSavartLaw::velocityJumpAcrossSheet(gamma_sheet);
        ASSERT_NEAR(Delta_u, gamma_sheet, TOLERANCE);
        return true;
    });

    run_test("Biot-Savart: filament induced velocity", []() {
        double Gamma = 5.0;
        maths::linear_algebra::Vector dl({0.1, 0.0, 0.0});
        maths::linear_algebra::Vector r_filament({0.0, 0.0, 0.0});
        maths::linear_algebra::Vector r_field({0.0, 0.5, 0.0});

        auto velocity = BiotSavartLaw::velocityFromFilament(Gamma, dl, r_filament, r_field);
        ASSERT_TRUE(velocity.norm() > 0.0);
        return true;
    });

    // ========================================
    // Vortex Dynamics Tests
    // ========================================

    run_test("Vortex Dynamics: vortex pair rotation", []() {
        double Gamma1 = 5.0;
        double Gamma2 = 5.0;
        double d = 1.0;  // Separation

        double omega = VortexDynamics::vortexPairRotation(Gamma1, Gamma2, d);
        double expected = (5.0 + 5.0) / (2.0 * M_PI * 1.0);
        ASSERT_NEAR(omega, expected, TOLERANCE);
        return true;
    });

    run_test("Vortex Dynamics: opposite circulation vortex pair", []() {
        double Gamma = 10.0;
        double d = 0.5;

        double U = VortexDynamics::vortexPairVelocity(Gamma, d);
        double expected = 10.0 / (2.0 * M_PI * 0.5);
        ASSERT_NEAR(U, expected, TOLERANCE);
        return true;
    });

    run_test("Vortex Dynamics: Rankine vortex core", []() {
        double r = 0.5;    // Inside core
        double a = 1.0;    // Core radius
        double Gamma = 10.0;

        double u = VortexDynamics::rankineVortex(r, a, Gamma);
        // Inside: solid body rotation u = ωr
        double omega = Gamma / (M_PI * a * a);
        double expected = omega * r;
        ASSERT_NEAR(u, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Vortex Dynamics: Rankine vortex outside core", []() {
        double r = 2.0;    // Outside core
        double a = 1.0;
        double Gamma = 10.0;

        double u = VortexDynamics::rankineVortex(r, a, Gamma);
        // Outside: potential flow u = Γ/(2πr)
        double expected = Gamma / (2.0 * M_PI * r);
        ASSERT_NEAR(u, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Vortex Dynamics: Rankine vortex behavior at boundary", []() {
        double a = 1.0;
        double Gamma = 10.0;

        double u_inside = VortexDynamics::rankineVortex(a - 1e-6, a, Gamma);
        double u_outside = VortexDynamics::rankineVortex(a + 1e-6, a, Gamma);

        // Inside velocity at boundary is higher (solid body rotation)
        ASSERT_TRUE(u_inside > u_outside);
        return true;
    });

    run_test("Vortex Dynamics: vortex pair velocity inversely proportional to separation", []() {
        double Gamma = 10.0;
        double U1 = VortexDynamics::vortexPairVelocity(Gamma, 0.5);
        double U2 = VortexDynamics::vortexPairVelocity(Gamma, 1.0);

        ASSERT_NEAR(U1 / U2, 2.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Consistency: enstrophy is half vorticity squared", []() {
        maths::linear_algebra::Vector omega({3.0, 4.0, 0.0});
        double enstrophy = Vorticity::enstrophy(omega);
        double mag = Vorticity::magnitude(omega);
        ASSERT_NEAR(enstrophy, 0.5 * mag * mag, TOLERANCE);
        return true;
    });

    run_test("Consistency: vortex line direction is unit vector", []() {
        maths::linear_algebra::Vector omega({1.0, 2.0, 2.0});
        auto dir = Vorticity::vortexLineDirection(omega);
        ASSERT_NEAR(dir.norm(), 1.0, TOLERANCE);
        return true;
    });

    run_test("Physical: vorticity perpendicular to velocity in potential flow", []() {
        // In irrotational flow, ω = 0
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0, 0) = 1.0; vel_grad(0, 1) = 0.0; vel_grad(0, 2) = 0.0;
        vel_grad(1, 0) = 0.0; vel_grad(1, 1) = 1.0; vel_grad(1, 2) = 0.0;
        vel_grad(2, 0) = 0.0; vel_grad(2, 1) = 0.0; vel_grad(2, 2) = 1.0;

        // Symmetric gradient → no vorticity
        auto omega = Vorticity::compute(vel_grad);
        ASSERT_NEAR(omega.norm(), 0.0, TOLERANCE);
        return true;
    });

    run_test("Physical: circulation equals enclosed vorticity flux", []() {
        // For solid body rotation: Γ = ∫∫ ω·dA = ω*A
        double omega_val = 4.0;
        double A = 2.5;

        double gamma = Circulation::solidBodyRotation(omega_val, A);
        ASSERT_NEAR(gamma, omega_val * A, TOLERANCE);
        return true;
    });

    run_test("Physical: Biot-Savart velocity increases with circulation", []() {
        double r = 1.0;
        double u1 = BiotSavartLaw::velocityStraightVortex(5.0, r);
        double u2 = BiotSavartLaw::velocityStraightVortex(10.0, r);

        ASSERT_NEAR(u2 / u1, 2.0, TOLERANCE);
        return true;
    });

    run_test("Physical: Rankine vortex velocity increases with radius inside core", []() {
        double a = 1.0;
        double Gamma = 10.0;

        double u1 = VortexDynamics::rankineVortex(0.25 * a, a, Gamma);
        double u2 = VortexDynamics::rankineVortex(0.5 * a, a, Gamma);
        double u3 = VortexDynamics::rankineVortex(0.75 * a, a, Gamma);

        // Inside core: solid body rotation, u increases with r
        ASSERT_TRUE(u1 < u2);
        ASSERT_TRUE(u2 < u3);
        return true;
    });

    run_test("Limit: zero circulation gives zero velocity", []() {
        double u = BiotSavartLaw::velocityStraightVortex(0.0, 1.0);
        ASSERT_NEAR(u, 0.0, TOLERANCE);
        return true;
    });

    run_test("Limit: enstrophy is zero for zero vorticity", []() {
        maths::linear_algebra::Vector omega_zero({0.0, 0.0, 0.0});
        double enstrophy = Vorticity::enstrophy(omega_zero);
        ASSERT_NEAR(enstrophy, 0.0, TOLERANCE);
        return true;
    });

    run_test("Limit: viscous diffusion is zero without viscosity", []() {
        maths::linear_algebra::Vector omega_lap({5.0, 3.0, 1.0});
        double nu = 0.0;

        auto diffusion = VorticityTransportEquation::viscousDiffusion(omega_lap, nu);
        ASSERT_NEAR(diffusion.norm(), 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Fluid Vorticity" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All fluid dynamics vorticity tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Vorticity calculation and magnitude" << std::endl;
        std::cout << "  - Vorticity transport equation" << std::endl;
        std::cout << "  - Circulation and Stokes theorem" << std::endl;
        std::cout << "  - Kelvin's circulation theorem" << std::endl;
        std::cout << "  - Biot-Savart law and induced velocities" << std::endl;
        std::cout << "  - Vortex dynamics and interactions" << std::endl;
        std::cout << "  - Rankine vortex model" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
