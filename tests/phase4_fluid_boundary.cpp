/**
 * Phase 4 Validation: Fluid Dynamics Boundary Layer
 *
 * Tests the fluid_dynamics_boundary_layer.hpp module functions.
 *
 * Coverage:
 * - Blasius solution (laminar flat plate)
 * - Boundary layer thickness definitions
 * - Von Karman momentum integral
 * - Skin friction coefficients
 * - Turbulent boundary layer correlations
 * - Transition criteria
 * - Velocity profiles (laminar and turbulent)
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/fluid_dynamics_boundary_layer.hpp"

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

    std::cout << "=== Phase 4: Fluid Dynamics Boundary Layer Validation ===" << std::endl;
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
    // Blasius Solution Tests
    // ========================================

    run_test("Blasius: similarity variable", []() {
        double y = 0.005;      // 5 mm from wall
        double x = 1.0;        // 1 m from leading edge
        double U = 10.0;       // m/s
        double nu = 1.5e-5;    // Air kinematic viscosity

        double eta = BlasiusSolution::similarityVariable(y, x, U, nu);
        double expected = 0.005 * std::sqrt(10.0 / (1.5e-5 * 1.0));
        ASSERT_NEAR(eta, expected, TOLERANCE);
        return true;
    });

    run_test("Blasius: velocity profile at wall", []() {
        double eta = 0.0;
        double u_ratio = BlasiusSolution::velocityProfile(eta);
        ASSERT_NEAR(u_ratio, 0.0, TOLERANCE);  // No-slip condition
        return true;
    });

    run_test("Blasius: velocity profile at edge", []() {
        double eta = 6.0;  // Far from wall
        double u_ratio = BlasiusSolution::velocityProfile(eta);
        ASSERT_NEAR(u_ratio, 1.0, TOLERANCE);  // Approaches free stream
        return true;
    });

    run_test("Blasius: velocity profile increases monotonically", []() {
        double u1 = BlasiusSolution::velocityProfile(0.5);
        double u2 = BlasiusSolution::velocityProfile(1.0);
        double u3 = BlasiusSolution::velocityProfile(1.5);

        ASSERT_TRUE(u1 < u2);
        ASSERT_TRUE(u2 < u3);
        return true;
    });

    run_test("Blasius: boundary layer thickness", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double delta = BlasiusSolution::thickness(x, U, nu);
        double expected = 5.0 * std::sqrt(1.5e-5 * 1.0 / 10.0);
        ASSERT_NEAR(delta, expected, TOLERANCE);
        return true;
    });

    run_test("Blasius: displacement thickness", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double delta_star = BlasiusSolution::displacementThickness(x, U, nu);
        double expected = 1.721 * std::sqrt(1.5e-5 * 1.0 / 10.0);
        ASSERT_NEAR(delta_star, expected, TOLERANCE);
        return true;
    });

    run_test("Blasius: momentum thickness", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double theta = BlasiusSolution::momentumThickness(x, U, nu);
        double expected = 0.664 * std::sqrt(1.5e-5 * 1.0 / 10.0);
        ASSERT_NEAR(theta, expected, TOLERANCE);
        return true;
    });

    run_test("Blasius: shape factor", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double delta_star = BlasiusSolution::displacementThickness(x, U, nu);
        double theta = BlasiusSolution::momentumThickness(x, U, nu);

        double H = BlasiusSolution::shapeFactor(delta_star, theta);
        double expected = 1.721 / 0.664;  // ≈ 2.59
        ASSERT_NEAR(H, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Blasius: skin friction coefficient", []() {
        double Re_x = 100000.0;

        double C_f = BlasiusSolution::skinFrictionCoefficient(Re_x);
        double expected = 0.664 / std::sqrt(100000.0);
        ASSERT_NEAR(C_f, expected, TOLERANCE);
        return true;
    });

    run_test("Blasius: wall shear stress", []() {
        double rho = 1.225;    // Air
        double U = 20.0;       // m/s
        double x = 0.5;        // m
        double nu = 1.5e-5;

        double tau_w = BlasiusSolution::wallShearStress(rho, U, x, nu);
        ASSERT_TRUE(tau_w > 0.0);
        return true;
    });

    run_test("Blasius: total drag on flat plate", []() {
        double rho = 1.225;
        double U = 30.0;
        double L = 1.0;
        double W = 0.5;
        double nu = 1.5e-5;

        double F_D = BlasiusSolution::totalDrag(rho, U, L, W, nu);
        ASSERT_TRUE(F_D > 0.0);
        return true;
    });

    // ========================================
    // Boundary Layer Thickness Tests
    // ========================================

    run_test("Boundary Layer: thickness grows with distance", []() {
        double U = 10.0;
        double nu = 1.5e-5;

        double delta1 = BlasiusSolution::thickness(1.0, U, nu);
        double delta2 = BlasiusSolution::thickness(4.0, U, nu);

        ASSERT_NEAR(delta2 / delta1, 2.0, LOOSE_TOLERANCE);  // δ ∝ √x
        return true;
    });

    run_test("Boundary Layer: displacement thickness relation", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double delta = BlasiusSolution::thickness(x, U, nu);
        double delta_star = BlasiusSolution::displacementThickness(x, U, nu);

        // δ*/δ ≈ 1.721/5.0 ≈ 0.344
        double ratio = delta_star / delta;
        ASSERT_NEAR(ratio, 1.721/5.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Boundary Layer: momentum thickness relation", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double delta = BlasiusSolution::thickness(x, U, nu);
        double theta = BlasiusSolution::momentumThickness(x, U, nu);

        // θ/δ ≈ 0.664/5.0 ≈ 0.133
        double ratio = theta / delta;
        ASSERT_NEAR(ratio, 0.664/5.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Boundary Layer: thickness decreases with velocity", []() {
        double x = 1.0;
        double nu = 1.5e-5;

        double delta1 = BlasiusSolution::thickness(x, 10.0, nu);
        double delta2 = BlasiusSolution::thickness(x, 40.0, nu);

        ASSERT_NEAR(delta1 / delta2, 2.0, LOOSE_TOLERANCE);  // δ ∝ 1/√U
        return true;
    });

    // ========================================
    // Von Karman Integral Tests
    // ========================================

    run_test("Von Karman: momentum thickness gradient", []() {
        double theta = 0.001;
        double C_f = 0.005;
        double H = 2.5;
        double U = 10.0;
        double dU_dx = 0.0;  // Zero pressure gradient

        double dtheta_dx = VonKarmanIntegral::momentumThicknessGradient(
            theta, C_f, H, U, dU_dx);

        ASSERT_NEAR(dtheta_dx, C_f/2.0, TOLERANCE);
        return true;
    });

    run_test("Von Karman: skin friction from gradient", []() {
        double dtheta_dx = 0.0025;

        double C_f = VonKarmanIntegral::skinFrictionFromGradient(dtheta_dx);
        ASSERT_NEAR(C_f, 0.005, TOLERANCE);
        return true;
    });

    run_test("Von Karman: separation criterion", []() {
        double dtheta_dx_nosep = 0.01;
        double dtheta_dx_sep = 1e-10;  // Extremely small gradient → separation
        double U = 10.0;
        double nu = 1.5e-5;

        bool sep_nosep = VonKarmanIntegral::isSeparationImminent(dtheta_dx_nosep, U, nu);
        bool sep_sep = VonKarmanIntegral::isSeparationImminent(dtheta_dx_sep, U, nu);

        ASSERT_TRUE(!sep_nosep);  // Normal gradient: no separation
        ASSERT_TRUE(sep_sep);     // Extremely small gradient: separation imminent
        return true;
    });

    // ========================================
    // Turbulent Boundary Layer Tests
    // ========================================

    run_test("Turbulent: boundary layer thickness", []() {
        double x = 1.0;
        double Re_x = 1e6;

        double delta = TurbulentBoundaryLayer::thickness(x, Re_x);
        double expected = 0.37 * 1.0 / std::pow(1e6, 0.2);
        ASSERT_NEAR(delta, expected, TOLERANCE);
        return true;
    });

    run_test("Turbulent: skin friction coefficient", []() {
        double Re_x = 1e6;

        double C_f = TurbulentBoundaryLayer::skinFrictionCoefficient(Re_x);
        double expected = 0.059 / std::pow(1e6, 0.2);
        ASSERT_NEAR(C_f, expected, TOLERANCE);
        return true;
    });

    run_test("Turbulent: average drag coefficient", []() {
        double Re_L = 1e7;

        double C_D = TurbulentBoundaryLayer::avgDragCoefficient(Re_L);
        double expected = 0.074 / std::pow(1e7, 0.2) - 1700.0 / 1e7;
        ASSERT_NEAR(C_D, expected, TOLERANCE);
        return true;
    });

    run_test("Turbulent: power law velocity profile", []() {
        double y = 0.005;
        double delta = 0.01;
        int n = 7;

        double u_ratio = TurbulentBoundaryLayer::powerLawProfile(y, delta, n);
        double expected = std::pow(0.5, 1.0/7.0);  // y/delta = 0.5
        ASSERT_NEAR(u_ratio, expected, TOLERANCE);
        return true;
    });

    run_test("Turbulent: power law at wall", []() {
        double y = 0.0;
        double delta = 0.01;

        double u_ratio = TurbulentBoundaryLayer::powerLawProfile(y, delta);
        ASSERT_NEAR(u_ratio, 0.0, TOLERANCE);
        return true;
    });

    run_test("Turbulent: power law at edge", []() {
        double y = 0.01;
        double delta = 0.01;

        double u_ratio = TurbulentBoundaryLayer::powerLawProfile(y, delta);
        ASSERT_NEAR(u_ratio, 1.0, TOLERANCE);
        return true;
    });

    run_test("Turbulent: log law in viscous sublayer", []() {
        double y_plus = 3.0;  // In viscous sublayer

        double u_plus = TurbulentBoundaryLayer::logLawProfile(y_plus);
        ASSERT_NEAR(u_plus, y_plus, TOLERANCE);  // u⁺ = y⁺
        return true;
    });

    run_test("Turbulent: log law in log layer", []() {
        double y_plus = 100.0;
        double kappa = 0.41;
        double B = 5.0;

        double u_plus = TurbulentBoundaryLayer::logLawProfile(y_plus, kappa, B);
        double expected = (1.0/kappa) * std::log(y_plus) + B;
        ASSERT_NEAR(u_plus, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Turbulent: friction velocity", []() {
        double tau_w = 0.5;    // Pa
        double rho = 1.225;    // Air

        double u_tau = TurbulentBoundaryLayer::frictionVelocity(tau_w, rho);
        double expected = std::sqrt(0.5 / 1.225);
        ASSERT_NEAR(u_tau, expected, TOLERANCE);
        return true;
    });

    run_test("Turbulent: transition Reynolds number", []() {
        double Re_crit = TurbulentBoundaryLayer::transitionReynolds();
        ASSERT_NEAR(Re_crit, 5e5, TOLERANCE);
        return true;
    });

    // ========================================
    // Boundary Layer Parameters Tests
    // ========================================

    run_test("Boundary Layer Parameters: Blasius thicknesses", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        auto t = BoundaryLayerParameters::blasiusThicknesses(x, U, nu);

        ASSERT_TRUE(t.delta > 0.0);
        ASSERT_TRUE(t.delta_star > 0.0);
        ASSERT_TRUE(t.theta > 0.0);
        ASSERT_NEAR(t.H, 2.59, LOOSE_TOLERANCE);  // Blasius shape factor
        return true;
    });

    run_test("Boundary Layer Parameters: turbulent thicknesses", []() {
        double x = 1.0;
        double Re_x = 1e6;

        auto t = BoundaryLayerParameters::turbulentThicknesses(x, Re_x);

        ASSERT_TRUE(t.delta > 0.0);
        ASSERT_TRUE(t.delta_star > 0.0);
        ASSERT_TRUE(t.theta > 0.0);
        ASSERT_TRUE(t.H < 2.0);  // Turbulent H < laminar H
        return true;
    });

    run_test("Boundary Layer Parameters: laminar vs turbulent shape factor", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        auto t_lam = BoundaryLayerParameters::blasiusThicknesses(x, U, nu);
        auto t_turb = BoundaryLayerParameters::turbulentThicknesses(x, 1e6);

        ASSERT_TRUE(t_lam.H > t_turb.H);  // Laminar H ≈ 2.59, Turbulent H ≈ 1.3
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Consistency: skin friction decreases with Reynolds number", []() {
        double C_f1 = BlasiusSolution::skinFrictionCoefficient(1e5);
        double C_f2 = BlasiusSolution::skinFrictionCoefficient(1e6);

        ASSERT_TRUE(C_f1 > C_f2);  // C_f ∝ 1/√Re
        return true;
    });

    run_test("Consistency: displacement thickness less than delta99", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double delta = BlasiusSolution::thickness(x, U, nu);
        double delta_star = BlasiusSolution::displacementThickness(x, U, nu);

        ASSERT_TRUE(delta_star < delta);
        return true;
    });

    run_test("Consistency: momentum thickness less than displacement thickness", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;

        double delta_star = BlasiusSolution::displacementThickness(x, U, nu);
        double theta = BlasiusSolution::momentumThickness(x, U, nu);

        ASSERT_TRUE(theta < delta_star);
        return true;
    });

    run_test("Physical: boundary layer thickness increases downstream", []() {
        double U = 10.0;
        double nu = 1.5e-5;

        double delta1 = BlasiusSolution::thickness(0.5, U, nu);
        double delta2 = BlasiusSolution::thickness(1.0, U, nu);
        double delta3 = BlasiusSolution::thickness(2.0, U, nu);

        ASSERT_TRUE(delta1 < delta2);
        ASSERT_TRUE(delta2 < delta3);
        return true;
    });

    run_test("Physical: turbulent boundary layer thicker than laminar", []() {
        double x = 1.0;
        double U = 10.0;
        double nu = 1.5e-5;
        double Re_x = U * x / nu;

        double delta_lam = BlasiusSolution::thickness(x, U, nu);
        double delta_turb = TurbulentBoundaryLayer::thickness(x, Re_x);

        ASSERT_TRUE(delta_turb > delta_lam);
        return true;
    });

    run_test("Physical: turbulent skin friction higher than laminar", []() {
        double Re_x = 1e6;

        double C_f_lam = BlasiusSolution::skinFrictionCoefficient(Re_x);
        double C_f_turb = TurbulentBoundaryLayer::skinFrictionCoefficient(Re_x);

        ASSERT_TRUE(C_f_turb > C_f_lam);
        return true;
    });

    run_test("Physical: power law profile is monotonic", []() {
        double delta = 0.01;
        double u1 = TurbulentBoundaryLayer::powerLawProfile(0.002, delta);
        double u2 = TurbulentBoundaryLayer::powerLawProfile(0.005, delta);
        double u3 = TurbulentBoundaryLayer::powerLawProfile(0.008, delta);

        ASSERT_TRUE(u1 < u2);
        ASSERT_TRUE(u2 < u3);
        return true;
    });

    run_test("Physical: log law increases with y plus", []() {
        double u1 = TurbulentBoundaryLayer::logLawProfile(50.0);
        double u2 = TurbulentBoundaryLayer::logLawProfile(100.0);
        double u3 = TurbulentBoundaryLayer::logLawProfile(200.0);

        ASSERT_TRUE(u1 < u2);
        ASSERT_TRUE(u2 < u3);
        return true;
    });

    run_test("Physical: wall shear stress decreases downstream", []() {
        double rho = 1.225;
        double U = 20.0;
        double nu = 1.5e-5;

        double tau1 = BlasiusSolution::wallShearStress(rho, U, 0.5, nu);
        double tau2 = BlasiusSolution::wallShearStress(rho, U, 1.0, nu);

        ASSERT_TRUE(tau1 > tau2);  // τ_w ∝ 1/√x
        return true;
    });

    run_test("Limit: boundary layer thickness zero at leading edge", []() {
        double x = 0.0001;  // Very close to leading edge
        double U = 10.0;
        double nu = 1.5e-5;

        double delta = BlasiusSolution::thickness(x, U, nu);
        ASSERT_TRUE(delta < 0.001);  // Very small
        return true;
    });

    run_test("Limit: Blasius velocity profile at large eta", []() {
        double eta_large = 10.0;
        double u_ratio = BlasiusSolution::velocityProfile(eta_large);
        ASSERT_NEAR(u_ratio, 1.0, TOLERANCE);  // Approaches free stream
        return true;
    });

    run_test("Limit: power law decreases as y decreases", []() {
        double delta = 0.01;
        double u1 = TurbulentBoundaryLayer::powerLawProfile(0.005, delta);
        double u2 = TurbulentBoundaryLayer::powerLawProfile(0.001, delta);
        double u3 = TurbulentBoundaryLayer::powerLawProfile(0.0001, delta);

        // Velocity should decrease as y decreases
        ASSERT_TRUE(u2 < u1);
        ASSERT_TRUE(u3 < u2);
        return true;
    });

    run_test("Limit: shape factor consistency", []() {
        double delta_star = 0.01721;
        double theta = 0.00664;

        double H = BlasiusSolution::shapeFactor(delta_star, theta);
        ASSERT_NEAR(H, 2.59, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Fluid Boundary Layer" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All fluid dynamics boundary layer tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Blasius solution for laminar flow" << std::endl;
        std::cout << "  - Boundary layer thickness definitions" << std::endl;
        std::cout << "  - Von Karman momentum integral" << std::endl;
        std::cout << "  - Skin friction coefficients" << std::endl;
        std::cout << "  - Turbulent boundary layer correlations" << std::endl;
        std::cout << "  - Velocity profiles (laminar and turbulent)" << std::endl;
        std::cout << "  - Physical consistency and limiting behavior" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
