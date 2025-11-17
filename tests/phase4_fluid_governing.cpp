/**
 * Phase 4 Validation: Fluid Dynamics Governing Equations
 *
 * Tests the fluid_dynamics_governing_equations.hpp module functions.
 *
 * Coverage:
 * - Continuity equation (mass conservation)
 * - Navier-Stokes equations (momentum conservation)
 * - Euler equations (inviscid flow)
 * - Bernoulli equation (energy conservation)
 * - Energy equation (thermodynamics)
 * - Conservation laws and physical relationships
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/fluid_dynamics_governing_equations.hpp"

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

    std::cout << "=== Phase 4: Fluid Dynamics Governing Equations Validation ===" << std::endl;
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
    // Continuity Equation Tests
    // ========================================

    run_test("Continuity: zero divergence is incompressible", []() {
        bool incomp = ContinuityEquation::isIncompressible(0.0);
        ASSERT_TRUE(incomp);
        return true;
    });

    run_test("Continuity: small divergence is incompressible", []() {
        bool incomp = ContinuityEquation::isIncompressible(1e-8);
        ASSERT_TRUE(incomp);
        return true;
    });

    run_test("Continuity: large divergence is compressible", []() {
        bool comp = !ContinuityEquation::isIncompressible(0.1);
        ASSERT_TRUE(comp);
        return true;
    });

    run_test("Continuity: uniform flow has zero divergence", []() {
        // Uniform flow: u = (1, 0, 0) everywhere
        auto uniform_field = [](const maths::linear_algebra::Vector& pos) {
            return maths::linear_algebra::Vector({1.0, 0.0, 0.0});
        };
        maths::linear_algebra::Vector pos({0.0, 0.0, 0.0});
        double div = ContinuityEquation::velocityDivergence(uniform_field, pos);
        ASSERT_NEAR(div, 0.0, NUMERICAL_TOLERANCE);
        return true;
    });

    run_test("Continuity: source flow has positive divergence", []() {
        // Radial outflow: u = r (expanding)
        auto source_field = [](const maths::linear_algebra::Vector& pos) {
            return pos;  // Velocity proportional to position
        };
        maths::linear_algebra::Vector pos({1.0, 1.0, 1.0});
        double div = ContinuityEquation::velocityDivergence(source_field, pos);
        ASSERT_TRUE(div > 0.0);  // Diverging flow
        return true;
    });

    run_test("Continuity: mass flux through surface", []() {
        FluidState state;
        state.density = 1.225;  // Air kg/m³
        state.velocity = maths::linear_algebra::Vector({10.0, 0.0, 0.0});

        maths::linear_algebra::Vector normal({1.0, 0.0, 0.0});
        double area = 5.0;  // m²

        double flux = ContinuityEquation::massFlux(state, normal, area);
        double expected = 1.225 * 10.0 * 5.0;  // ρ * u * A = 61.25 kg/s
        ASSERT_NEAR(flux, expected, TOLERANCE);
        return true;
    });

    run_test("Continuity: perpendicular flow has zero flux", []() {
        FluidState state;
        state.density = 1000.0;  // Water kg/m³
        state.velocity = maths::linear_algebra::Vector({10.0, 0.0, 0.0});

        maths::linear_algebra::Vector normal({0.0, 1.0, 0.0});  // Perpendicular
        double area = 2.0;

        double flux = ContinuityEquation::massFlux(state, normal, area);
        ASSERT_NEAR(flux, 0.0, TOLERANCE);
        return true;
    });

    run_test("Continuity: density time derivative", []() {
        FluidState state;
        state.density = 1.2;
        state.velocity = maths::linear_algebra::Vector({1.0, 0.0, 0.0});

        double vel_div = 0.5;  // Diverging
        maths::linear_algebra::Vector density_grad({0.1, 0.0, 0.0});

        double drho_dt = ContinuityEquation::densityTimeDerivative(state, vel_div, density_grad);
        double expected = -1.2 * 0.5 - 1.0 * 0.1;  // -0.6 - 0.1 = -0.7
        ASSERT_NEAR(drho_dt, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Navier-Stokes Equation Tests
    // ========================================

    run_test("Navier-Stokes: kinematic viscosity calculation", []() {
        double mu = 1.81e-5;  // Air dynamic viscosity (Pa·s)
        double rho = 1.225;   // Air density (kg/m³)
        double nu = NavierStokesEquations::kinematicViscosity(mu, rho);
        double expected = 1.81e-5 / 1.225;  // ≈ 1.478e-5 m²/s
        ASSERT_NEAR(nu, expected, TOLERANCE);
        return true;
    });

    run_test("Navier-Stokes: water kinematic viscosity", []() {
        double mu = 1e-3;   // Water dynamic viscosity (Pa·s)
        double rho = 1000.0;  // Water density (kg/m³)
        double nu = NavierStokesEquations::kinematicViscosity(mu, rho);
        double expected = 1e-6;  // m²/s
        ASSERT_NEAR(nu, expected, TOLERANCE);
        return true;
    });

    run_test("Navier-Stokes: convective acceleration", []() {
        maths::linear_algebra::Vector velocity({2.0, 1.0, 0.0});
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0, 0) = 0.5; vel_grad(0, 1) = 0.0; vel_grad(0, 2) = 0.0;
        vel_grad(1, 0) = 0.0; vel_grad(1, 1) = 0.3; vel_grad(1, 2) = 0.0;
        vel_grad(2, 0) = 0.0; vel_grad(2, 1) = 0.0; vel_grad(2, 2) = 0.0;

        auto accel = NavierStokesEquations::convectiveAcceleration(velocity, vel_grad);
        // (u·∇)u = [u∂u/∂x + v∂u/∂y, u∂v/∂x + v∂v/∂y, 0]
        // = [2*0.5 + 1*0, 2*0 + 1*0.3, 0] = [1, 0.3, 0]
        ASSERT_NEAR(accel[0], 1.0, TOLERANCE);
        ASSERT_NEAR(accel[1], 0.3, TOLERANCE);
        return true;
    });

    run_test("Navier-Stokes: viscous term", []() {
        maths::linear_algebra::Vector laplacian({2.0, 1.0, 0.0});
        double mu = 1e-3;
        double rho = 1000.0;

        auto viscous = NavierStokesEquations::viscousTerm(laplacian, mu, rho);
        double nu = mu / rho;  // 1e-6
        ASSERT_NEAR(viscous[0], 2.0 * nu, TOLERANCE);
        ASSERT_NEAR(viscous[1], 1.0 * nu, TOLERANCE);
        return true;
    });

    run_test("Navier-Stokes: total acceleration with gravity", []() {
        maths::linear_algebra::Vector pressure_grad({-100.0, 0.0, 0.0});
        maths::linear_algebra::Vector laplacian({0.0, 0.0, 0.0});
        double rho = 1000.0;
        double nu = 1e-6;
        maths::linear_algebra::Vector gravity({0.0, 0.0, -9.81});

        auto accel = NavierStokesEquations::totalAcceleration(
            pressure_grad, laplacian, rho, nu, gravity);

        // a = -∇p/ρ + ν∇²u + g = [100/1000, 0, 0] + [0, 0, 0] + [0, 0, -9.81]
        ASSERT_NEAR(accel[0], 0.1, TOLERANCE);
        ASSERT_NEAR(accel[2], -9.81, TOLERANCE);
        return true;
    });

    run_test("Navier-Stokes: pressure gradient drives flow", []() {
        maths::linear_algebra::Vector pressure_grad({-1000.0, 0.0, 0.0});
        maths::linear_algebra::Vector laplacian({0.0, 0.0, 0.0});
        double rho = 1.225;
        double nu = 1.5e-5;
        maths::linear_algebra::Vector gravity({0.0, 0.0, 0.0});

        auto accel = NavierStokesEquations::totalAcceleration(
            pressure_grad, laplacian, rho, nu, gravity);

        double expected = 1000.0 / 1.225;  // ≈ 816.3 m/s²
        ASSERT_NEAR(accel[0], expected, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Euler Equation Tests (Inviscid)
    // ========================================

    run_test("Euler: inviscid acceleration", []() {
        maths::linear_algebra::Vector pressure_grad({-500.0, 0.0, 0.0});
        double rho = 1.2;
        maths::linear_algebra::Vector gravity({0.0, 0.0, -9.81});

        auto accel = EulerEquations::eulerAcceleration(pressure_grad, rho, gravity);
        double expected_x = 500.0 / 1.2;  // ≈ 416.67 m/s²
        ASSERT_NEAR(accel[0], expected_x, LOOSE_TOLERANCE);
        ASSERT_NEAR(accel[2], -9.81, TOLERANCE);
        return true;
    });

    run_test("Euler: high Reynolds number validates inviscid", []() {
        double Re = 10000.0;
        bool valid = EulerEquations::isInviscidValid(Re);
        ASSERT_TRUE(valid);
        return true;
    });

    run_test("Euler: low Reynolds number invalidates inviscid", []() {
        double Re = 100.0;
        bool valid = EulerEquations::isInviscidValid(Re);
        ASSERT_TRUE(!valid);
        return true;
    });

    run_test("Euler: primitive to conservative variables", []() {
        FluidState state;
        state.density = 1.2;
        state.velocity = maths::linear_algebra::Vector({10.0, 5.0, 0.0});
        state.internal_energy = 1000.0;

        auto U = EulerEquations::primitiveToConservative(state);
        ASSERT_NEAR(U[0], 1.2, TOLERANCE);  // ρ
        ASSERT_NEAR(U[1], 12.0, TOLERANCE);  // ρu
        ASSERT_NEAR(U[2], 6.0, TOLERANCE);   // ρv
        // U[4] = ρ(e + 0.5v²) = 1.2*(1000 + 0.5*125) = 1.2*1062.5
        ASSERT_NEAR(U[4], 1275.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Bernoulli Equation Tests
    // ========================================

    run_test("Bernoulli: dynamic pressure", []() {
        double rho = 1.225;  // Air
        maths::linear_algebra::Vector velocity({30.0, 0.0, 0.0});  // 30 m/s
        double q = BernoulliEquation::dynamicPressure(rho, velocity);
        double expected = 0.5 * 1.225 * 30.0 * 30.0;  // 551.25 Pa
        ASSERT_NEAR(q, expected, TOLERANCE);
        return true;
    });

    run_test("Bernoulli: total pressure", []() {
        double p_static = 101325.0;  // 1 atm
        double rho = 1.225;
        maths::linear_algebra::Vector velocity({20.0, 0.0, 0.0});

        double p_total = BernoulliEquation::totalPressure(p_static, rho, velocity);
        double q = 0.5 * 1.225 * 20.0 * 20.0;  // 245 Pa
        ASSERT_NEAR(p_total, p_static + q, TOLERANCE);
        return true;
    });

    run_test("Bernoulli: velocity from pressure difference", []() {
        double delta_p = 1000.0;  // 1000 Pa
        double rho = 1000.0;  // Water

        double v = BernoulliEquation::velocityFromPressure(delta_p, rho);
        double expected = std::sqrt(2.0 * 1000.0 / 1000.0);  // sqrt(2) ≈ 1.414 m/s
        ASSERT_NEAR(v, expected, TOLERANCE);
        return true;
    });

    run_test("Bernoulli: pressure from velocity change", []() {
        double p1 = 101325.0;
        double rho = 1.225;
        maths::linear_algebra::Vector v1({20.0, 0.0, 0.0});
        maths::linear_algebra::Vector v2({30.0, 0.0, 0.0});
        double h1 = 0.0;
        double h2 = 0.0;

        double p2 = BernoulliEquation::pressureFromVelocity(p1, rho, v1, v2, h1, h2);
        // p2 = p1 + 0.5*rho*(v1² - v2²) = 101325 + 0.5*1.225*(400-900)
        double expected = p1 + 0.5 * 1.225 * (400.0 - 900.0);
        ASSERT_NEAR(p2, expected, TOLERANCE);
        return true;
    });

    run_test("Bernoulli: total head", []() {
        double p = 50000.0;  // 50 kPa
        double rho = 1000.0;  // Water
        maths::linear_algebra::Vector velocity({2.0, 0.0, 0.0});
        double z = 10.0;  // 10 m height
        double g = 9.81;

        double H = BernoulliEquation::totalHead(p, rho, velocity, z, g);
        double expected = 50000.0/(1000.0*9.81) + (4.0)/(2.0*9.81) + 10.0;
        // = 5.097 + 0.204 + 10 ≈ 15.3 m
        ASSERT_NEAR(H, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Bernoulli: Venturi effect", []() {
        // Wide section: low velocity, high pressure
        // Narrow section: high velocity, low pressure
        double p1 = 101325.0;
        double rho = 1000.0;
        maths::linear_algebra::Vector v1({1.0, 0.0, 0.0});
        maths::linear_algebra::Vector v2({4.0, 0.0, 0.0});

        double p2 = BernoulliEquation::pressureFromVelocity(p1, rho, v1, v2, 0.0, 0.0);
        // p2 = p1 + 0.5*ρ*(v1² - v2²) = 101325 + 0.5*1000*(1-16)
        double expected = p1 - 7500.0;  // Pressure drops
        ASSERT_NEAR(p2, expected, TOLERANCE);
        ASSERT_TRUE(p2 < p1);  // Pressure decreases with velocity increase
        return true;
    });

    run_test("Bernoulli: hydrostatic pressure with height", []() {
        double p1 = 101325.0;
        double rho = 1000.0;
        maths::linear_algebra::Vector v1({0.0, 0.0, 0.0});
        maths::linear_algebra::Vector v2({0.0, 0.0, 0.0});
        double h1 = 0.0;
        double h2 = 5.0;  // 5 m higher

        double p2 = BernoulliEquation::pressureFromVelocity(p1, rho, v1, v2, h1, h2);
        // p2 = p1 + ρg(h1 - h2) = 101325 - 1000*9.81*5
        double expected = p1 - 49050.0;
        ASSERT_NEAR(p2, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Energy Equation Tests
    // ========================================

    run_test("Energy: enthalpy calculation", []() {
        double e = 1000.0;  // Internal energy (J/kg)
        double p = 100000.0;  // Pressure (Pa)
        double rho = 1.2;  // Density (kg/m³)

        double h = EnergyEquation::enthalpy(e, p, rho);
        double expected = 1000.0 + 100000.0/1.2;  // ≈ 84333.3 J/kg
        ASSERT_NEAR(h, expected, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Energy: total energy", []() {
        double e = 500.0;  // Internal energy
        maths::linear_algebra::Vector velocity({10.0, 0.0, 0.0});

        double E = EnergyEquation::totalEnergy(e, velocity);
        double expected = 500.0 + 0.5 * 100.0;  // 500 + 50 = 550 J/kg
        ASSERT_NEAR(E, expected, TOLERANCE);
        return true;
    });

    run_test("Energy: total enthalpy stagnation", []() {
        double h = 1000.0;
        maths::linear_algebra::Vector velocity({20.0, 0.0, 0.0});

        double h0 = EnergyEquation::totalEnthalpy(h, velocity);
        double expected = 1000.0 + 0.5 * 400.0;  // 1000 + 200 = 1200 J/kg
        ASSERT_NEAR(h0, expected, TOLERANCE);
        return true;
    });

    run_test("Energy: thermal conduction", []() {
        double temp_laplacian = -10.0;  // ∇²T (K/m²)
        double k = 0.025;  // Air thermal conductivity (W/(m·K))

        double Q_dot = EnergyEquation::thermalConduction(temp_laplacian, k);
        double expected = 0.025 * (-10.0);  // -0.25 W/m³
        ASSERT_NEAR(Q_dot, expected, TOLERANCE);
        return true;
    });

    run_test("Energy: viscous dissipation in simple shear", []() {
        // Simple shear flow: du/dy = constant
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0, 0) = 0.0; vel_grad(0, 1) = 10.0; vel_grad(0, 2) = 0.0;
        vel_grad(1, 0) = 0.0; vel_grad(1, 1) = 0.0;  vel_grad(1, 2) = 0.0;
        vel_grad(2, 0) = 0.0; vel_grad(2, 1) = 0.0;  vel_grad(2, 2) = 0.0;

        double mu = 1e-3;  // Water
        double phi = EnergyEquation::viscousDissipation(vel_grad, mu);
        ASSERT_TRUE(phi > 0.0);  // Dissipation is always positive
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Consistency: Bernoulli energy conservation", []() {
        // Along streamline: p/(ρg) + v²/(2g) + z = constant
        double rho = 1000.0;
        double g = 9.81;

        maths::linear_algebra::Vector v1({5.0, 0.0, 0.0});
        double p1 = 200000.0;
        double z1 = 0.0;

        maths::linear_algebra::Vector v2({10.0, 0.0, 0.0});
        double z2 = 2.0;
        double p2 = BernoulliEquation::pressureFromVelocity(p1, rho, v1, v2, z1, z2, g);

        double H1 = BernoulliEquation::totalHead(p1, rho, v1, z1, g);
        double H2 = BernoulliEquation::totalHead(p2, rho, v2, z2, g);

        ASSERT_NEAR(H1, H2, LOOSE_TOLERANCE);  // Total head conserved
        return true;
    });

    run_test("Consistency: enthalpy and total energy relationship", []() {
        double e = 500.0;
        double p = 100000.0;
        double rho = 1.2;
        maths::linear_algebra::Vector v({15.0, 0.0, 0.0});

        double h = EnergyEquation::enthalpy(e, p, rho);
        double E = EnergyEquation::totalEnergy(e, v);
        double h0 = EnergyEquation::totalEnthalpy(h, v);

        // h0 = h + v²/2 = e + p/ρ + v²/2
        // E = e + v²/2
        // So h0 = E + p/ρ
        double expected_h0 = E + p/rho;
        ASSERT_NEAR(h0, expected_h0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Consistency: stagnation pressure equals total pressure", []() {
        double p_static = 101325.0;
        double rho = 1.225;
        maths::linear_algebra::Vector v({50.0, 0.0, 0.0});

        double p_total = BernoulliEquation::totalPressure(p_static, rho, v);
        double q = BernoulliEquation::dynamicPressure(rho, v);

        ASSERT_NEAR(p_total, p_static + q, TOLERANCE);
        return true;
    });

    run_test("Physical: pressure gradient opposes flow direction", []() {
        // Positive pressure gradient (increasing p) creates negative acceleration
        maths::linear_algebra::Vector pressure_grad({1000.0, 0.0, 0.0});
        maths::linear_algebra::Vector laplacian({0.0, 0.0, 0.0});
        double rho = 1.2;
        double nu = 1.5e-5;
        maths::linear_algebra::Vector gravity({0.0, 0.0, 0.0});

        auto accel = NavierStokesEquations::totalAcceleration(
            pressure_grad, laplacian, rho, nu, gravity);

        ASSERT_TRUE(accel[0] < 0.0);  // Acceleration opposes pressure increase
        return true;
    });

    run_test("Physical: kinematic viscosity is positive", []() {
        double mu = 1.81e-5;
        double rho = 1.225;
        double nu = NavierStokesEquations::kinematicViscosity(mu, rho);
        ASSERT_TRUE(nu > 0.0);
        return true;
    });

    run_test("Physical: dynamic pressure is always positive", []() {
        double rho = 1.225;
        maths::linear_algebra::Vector v({-30.0, 20.0, 10.0});  // Any direction
        double q = BernoulliEquation::dynamicPressure(rho, v);
        ASSERT_TRUE(q > 0.0);
        return true;
    });

    run_test("Physical: total energy exceeds internal energy", []() {
        double e = 1000.0;
        maths::linear_algebra::Vector v({10.0, 5.0, 2.0});
        double E = EnergyEquation::totalEnergy(e, v);
        ASSERT_TRUE(E > e);
        return true;
    });

    run_test("Physical: total enthalpy exceeds static enthalpy", []() {
        double h = 2000.0;
        maths::linear_algebra::Vector v({15.0, 0.0, 0.0});
        double h0 = EnergyEquation::totalEnthalpy(h, v);
        ASSERT_TRUE(h0 > h);
        return true;
    });

    run_test("Limit: Bernoulli at zero velocity gives static pressure", []() {
        double p_static = 101325.0;
        double rho = 1.225;
        maths::linear_algebra::Vector v_zero({0.0, 0.0, 0.0});

        double p_total = BernoulliEquation::totalPressure(p_static, rho, v_zero);
        ASSERT_NEAR(p_total, p_static, TOLERANCE);
        return true;
    });

    run_test("Limit: total energy equals internal energy at rest", []() {
        double e = 800.0;
        maths::linear_algebra::Vector v_zero({0.0, 0.0, 0.0});
        double E = EnergyEquation::totalEnergy(e, v_zero);
        ASSERT_NEAR(E, e, TOLERANCE);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Fluid Governing Equations" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All fluid dynamics governing equations tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Continuity equation (mass conservation)" << std::endl;
        std::cout << "  - Navier-Stokes equations (momentum)" << std::endl;
        std::cout << "  - Euler equations (inviscid flow)" << std::endl;
        std::cout << "  - Bernoulli equation (energy)" << std::endl;
        std::cout << "  - Energy equation (thermodynamics)" << std::endl;
        std::cout << "  - Conservation laws and physical consistency" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
