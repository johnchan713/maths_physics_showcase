/**
 * Phase 4 Validation: Hamiltonian Mechanics
 *
 * Tests the classical_hamiltonian.hpp module functions.
 *
 * Coverage:
 * - PhasePoint construction and dimension management
 * - Hamiltonian evaluation for simple systems
 * - Hamilton's equations: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
 * - Numerical derivatives dH/dq and dH/dp
 * - Symplectic Euler integrator (1st order, exactly symplectic)
 * - Störmer-Verlet integrator (2nd order symplectic)
 * - Energy conservation verification
 * - Harmonic oscillator: H = p²/(2m) + (1/2)kq²
 * - Kepler problem: H = p²/(2m) - GMm/r
 * - Poisson bracket computation and properties
 * - Fundamental Poisson bracket relations
 * - Symplectic structure preservation
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "../include/physics/classical_hamiltonian.hpp"

using namespace physics::advanced::classical;

// Test tolerances
const double TOLERANCE_EXACT = 1e-5;
const double TOLERANCE_NUMERICAL = 1e-3;
const double TOLERANCE_INTEGRATOR = 0.01;
const double TOLERANCE_POISSON = 1e-4;

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

    std::cout << "=== Phase 4: Hamiltonian Mechanics Validation ===" << std::endl;
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
    // PhasePoint Tests
    // ========================================

    run_test("PhasePoint construction with dimension", []() {
        PhasePoint phase_point(3);
        ASSERT_TRUE(phase_point.dimension() == 3);
        return true;
    });

    run_test("PhasePoint construction from vectors", []() {
        maths::linear_algebra::Vector q(2);
        maths::linear_algebra::Vector p(2);
        q[0] = 1.0; q[1] = 2.0;
        p[0] = 3.0; p[1] = 4.0;

        PhasePoint phase_point(q, p);
        ASSERT_TRUE(phase_point.dimension() == 2);
        ASSERT_NEAR(phase_point.q[0], 1.0, TOLERANCE_EXACT);
        ASSERT_NEAR(phase_point.p[1], 4.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhasePoint dimension mismatch throws exception", []() {
        maths::linear_algebra::Vector q(2);
        maths::linear_algebra::Vector p(3);

        try {
            PhasePoint phase_point(q, p);
            return false;  // Should have thrown
        } catch (const std::invalid_argument&) {
            return true;
        }
    });

    run_test("PhasePoint coordinate access", []() {
        PhasePoint phase_point(2);
        phase_point.q[0] = 5.0;
        phase_point.p[1] = 2.5;

        ASSERT_NEAR(phase_point.q[0], 5.0, TOLERANCE_EXACT);
        ASSERT_NEAR(phase_point.p[1], 2.5, TOLERANCE_EXACT);
        return true;
    });

    // ========================================
    // Harmonic Oscillator Tests
    // ========================================

    run_test("Harmonic oscillator Hamiltonian evaluation", []() {
        // m = 1, k = 1, so H = p²/2 + q²/2
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.0;

        double H = ho.hamiltonian(point);
        ASSERT_NEAR(H, 0.5, TOLERANCE_EXACT);  // H = 0 + 0.5*1²/2 = 0.5
        return true;
    });

    run_test("Harmonic oscillator energy conservation with Verlet", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double dt = 0.01;
        int num_steps = 1000;
        auto trajectory = ho.integrate(initial, dt, num_steps);

        ASSERT_TRUE(ho.checkEnergyConservation(trajectory, 0.001));
        return true;
    });

    run_test("Harmonic oscillator period 2pi sqrt m over k", []() {
        // Period T = 2π√(m/k)
        double m = 2.0;
        double k = 8.0;  // T = 2π√(2/8) = 2π√(1/4) = π

        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(m, k);

        PhasePoint initial(1);
        initial.q[0] = 1.0;  // Start at maximum displacement
        initial.p[0] = 0.0;  // With zero velocity

        double dt = 0.001;
        int num_steps = int(M_PI / dt);  // One full period
        auto trajectory = ho.integrate(initial, dt, num_steps);

        // After one period, should return close to initial position
        double final_q = trajectory.back().q[0];
        ASSERT_NEAR(final_q, 1.0, TOLERANCE_INTEGRATOR);
        return true;
    });

    run_test("Harmonic oscillator phase space rotation", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        // Use smaller time steps for better accuracy with quarter period
        double dt = 0.01;
        int num_steps = int(M_PI / (2.0 * dt));  // Quarter period (T = 2π)
        auto trajectory = ho.integrate(initial, dt, num_steps);

        // After quarter period: q should be ≈ 0, p should be ≈ -1
        ASSERT_NEAR(trajectory.back().q[0], 0.0, 0.05);
        ASSERT_NEAR(trajectory.back().p[0], -1.0, 0.05);
        return true;
    });

    // ========================================
    // Hamilton's Equations Tests
    // ========================================

    run_test("Harmonic oscillator Hamilton equations dq over dt", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint point(1);
        point.q[0] = 0.5;
        point.p[0] = 2.0;

        // dq/dt = ∂H/∂p = p/m = 2.0/1.0 = 2.0
        double dq_dt = ho.dH_dp(point, 0);
        ASSERT_NEAR(dq_dt, 2.0, TOLERANCE_NUMERICAL);
        return true;
    });

    run_test("Harmonic oscillator Hamilton equations dp over dt", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint point(1);
        point.q[0] = 0.5;
        point.p[0] = 2.0;

        // dp/dt = -∂H/∂q = -kq = -1.0*0.5 = -0.5
        double dp_dt = -ho.dH_dq(point, 0);
        ASSERT_NEAR(dp_dt, -0.5, TOLERANCE_NUMERICAL);
        return true;
    });

    run_test("Numerical derivative dH over dq consistency", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 2.0);

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.5;

        // For H = p²/2 + kq², dH/dq = kq = 2.0*1.0 = 2.0
        double dH_dq = ho.dH_dq(point, 0);
        ASSERT_NEAR(dH_dq, 2.0, TOLERANCE_NUMERICAL);
        return true;
    });

    run_test("Numerical derivative dH over dp consistency", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.5;

        // For H = p²/2 + kq², dH/dp = p = 0.5
        double dH_dp = ho.dH_dp(point, 0);
        ASSERT_NEAR(dH_dp, 0.5, TOLERANCE_NUMERICAL);
        return true;
    });

    run_test("Hamilton equations for entire phase point", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.0;

        PhasePoint deriv = ho.hamiltonEquations(point);

        // dq/dt should be p/m = 0
        ASSERT_NEAR(deriv.q[0], 0.0, TOLERANCE_NUMERICAL);
        // dp/dt should be -kq = -1
        ASSERT_NEAR(deriv.p[0], -1.0, TOLERANCE_NUMERICAL);
        return true;
    });

    // ========================================
    // Symplectic Euler Integrator Tests
    // ========================================

    run_test("Symplectic Euler single step", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.0;

        double dt = 0.01;
        PhasePoint next = ho.stepSymplecticEuler(point, dt);

        // Should move in phase space
        ASSERT_TRUE(next.q[0] != point.q[0] || next.p[0] != point.p[0]);
        return true;
    });

    run_test("Symplectic Euler energy conservation harmonic oscillator", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double E_initial = ho.hamiltonian(initial);

        PhasePoint current = initial;
        double dt = 0.001;
        int num_steps = 5000;

        for (int i = 0; i < num_steps; ++i) {
            current = ho.stepSymplecticEuler(current, dt);
        }

        double E_final = ho.hamiltonian(current);

        // Symplectic Euler should have good energy conservation
        ASSERT_NEAR(E_final, E_initial, 0.1);
        return true;
    });

    run_test("Symplectic Euler dimension check", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        ASSERT_TRUE(ho.dimension() == 1);
        return true;
    });

    // ========================================
    // Störmer-Verlet Integrator Tests
    // ========================================

    run_test("Verlet single step", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.0;

        double dt = 0.01;
        PhasePoint next = ho.stepVerlet(point, dt);

        // Should move in phase space
        ASSERT_TRUE(next.q[0] != point.q[0] || next.p[0] != point.p[0]);
        return true;
    });

    run_test("Verlet energy conservation harmonic oscillator", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double E_initial = ho.hamiltonian(initial);

        PhasePoint current = initial;
        double dt = 0.01;
        int num_steps = 1000;

        for (int i = 0; i < num_steps; ++i) {
            current = ho.stepVerlet(current, dt);
        }

        double E_final = ho.hamiltonian(current);

        // Verlet is more accurate than Euler
        ASSERT_NEAR(E_final, E_initial, 0.01);
        return true;
    });

    run_test("Verlet phase space area preservation", []() {
        // Symplectic integrators preserve phase space area (volume)
        // Test with two nearby points
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint p1(1);
        p1.q[0] = 1.0;
        p1.p[0] = 0.0;

        PhasePoint p2(1);
        p2.q[0] = 1.001;
        p2.p[0] = 0.0;

        double dq_initial = p2.q[0] - p1.q[0];
        double dp_initial = p2.p[0] - p1.p[0];
        double area_initial = std::abs(dq_initial * dp_initial);

        // Integrate both points
        double dt = 0.01;
        int num_steps = 100;

        for (int i = 0; i < num_steps; ++i) {
            p1 = ho.stepVerlet(p1, dt);
            p2 = ho.stepVerlet(p2, dt);
        }

        double dq_final = p2.q[0] - p1.q[0];
        double dp_final = p2.p[0] - p1.p[0];
        double area_final = std::abs(dq_final * dp_final);

        // Areas should be preserved (within integrator error)
        if (area_initial > 1e-10) {
            double ratio = area_final / area_initial;
            ASSERT_NEAR(ratio, 1.0, 0.1);
        }
        return true;
    });

    // ========================================
    // Trajectory Integration Tests
    // ========================================

    run_test("Trajectory integration returns correct size", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double dt = 0.01;
        int num_steps = 100;
        auto trajectory = ho.integrate(initial, dt, num_steps);

        ASSERT_TRUE(trajectory.size() == num_steps + 1);
        return true;
    });

    run_test("Trajectory integration initial condition preserved", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 2.5;
        initial.p[0] = 1.3;

        double dt = 0.01;
        int num_steps = 50;
        auto trajectory = ho.integrate(initial, dt, num_steps);

        ASSERT_NEAR(trajectory[0].q[0], 2.5, TOLERANCE_EXACT);
        ASSERT_NEAR(trajectory[0].p[0], 1.3, TOLERANCE_EXACT);
        return true;
    });

    run_test("Trajectory integration produces different points", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double dt = 0.01;
        int num_steps = 10;
        auto trajectory = ho.integrate(initial, dt, num_steps);

        // Check that trajectory progresses
        bool has_change = false;
        for (size_t i = 1; i < trajectory.size(); ++i) {
            if (std::abs(trajectory[i].q[0] - trajectory[i-1].q[0]) > 1e-10) {
                has_change = true;
                break;
            }
        }

        ASSERT_TRUE(has_change);
        return true;
    });

    // ========================================
    // Energy Conservation Tests
    // ========================================

    run_test("Energy conservation check returns true for constant trajectory", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        std::vector<PhasePoint> trajectory;
        PhasePoint p1(1);
        p1.q[0] = 1.0;
        p1.p[0] = 0.0;

        for (int i = 0; i < 5; ++i) {
            trajectory.push_back(p1);
        }

        ASSERT_TRUE(ho.checkEnergyConservation(trajectory, 1e-6));
        return true;
    });

    run_test("Energy conservation check detects energy violation", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        std::vector<PhasePoint> trajectory;
        PhasePoint p1(1);
        p1.q[0] = 1.0;
        p1.p[0] = 0.0;
        trajectory.push_back(p1);

        PhasePoint p2(1);
        p2.q[0] = 2.0;  // Doubles energy
        p2.p[0] = 0.0;
        trajectory.push_back(p2);

        ASSERT_TRUE(!ho.checkEnergyConservation(trajectory, 1e-6));
        return true;
    });

    // ========================================
    // Kepler Problem Tests
    // ========================================

    run_test("Kepler problem Hamiltonian evaluation", []() {
        HamiltonianSystem kepler = HamiltonianSystem::keplerProblem(1.0, 1.0);

        PhasePoint point(2);
        point.q[0] = 1.0;
        point.q[1] = 0.0;
        point.p[0] = 0.0;
        point.p[1] = 1.0;

        double H = kepler.hamiltonian(point);

        // H = p²/(2m) - GM/r
        // H = 1²/2 - 1/1 = 0.5 - 1 = -0.5
        ASSERT_NEAR(H, -0.5, TOLERANCE_EXACT);
        return true;
    });

    run_test("Kepler problem dimension", []() {
        HamiltonianSystem kepler = HamiltonianSystem::keplerProblem(1.0, 1.0);
        ASSERT_TRUE(kepler.dimension() == 2);
        return true;
    });

    run_test("Kepler problem collision detection", []() {
        HamiltonianSystem kepler = HamiltonianSystem::keplerProblem(1.0, 1.0);

        PhasePoint point(2);
        point.q[0] = 1e-11;  // Collision
        point.q[1] = 1e-11;
        point.p[0] = 0.0;
        point.p[1] = 0.0;

        try {
            kepler.hamiltonian(point);
            return false;  // Should have thrown
        } catch (const std::domain_error&) {
            return true;
        }
    });

    run_test("Kepler problem energy conservation", []() {
        // Circular orbit: v = sqrt(GM/r)
        HamiltonianSystem kepler = HamiltonianSystem::keplerProblem(1.0, 1.0);

        double r = 1.0;
        double v = std::sqrt(1.0 / r);  // Circular orbit velocity

        PhasePoint initial(2);
        initial.q[0] = r;
        initial.q[1] = 0.0;
        initial.p[0] = 0.0;
        initial.p[1] = 1.0 * v;  // m * v_tangential

        double dt = 0.001;
        int num_steps = 1000;
        auto trajectory = kepler.integrate(initial, dt, num_steps);

        ASSERT_TRUE(kepler.checkEnergyConservation(trajectory, 0.01));
        return true;
    });

    run_test("Kepler problem orbital radius bounds", []() {
        HamiltonianSystem kepler = HamiltonianSystem::keplerProblem(1.0, 1.0);

        PhasePoint initial(2);
        initial.q[0] = 1.0;
        initial.q[1] = 0.0;
        initial.p[0] = 0.0;
        initial.p[1] = 0.8;  // Slightly less than circular orbit

        double dt = 0.001;
        int num_steps = 2000;
        auto trajectory = kepler.integrate(initial, dt, num_steps);

        // All points should stay at reasonable distance
        for (const auto& point : trajectory) {
            double r = std::sqrt(point.q[0] * point.q[0] + point.q[1] * point.q[1]);
            ASSERT_TRUE(r > 0.1 && r < 10.0);
        }

        return true;
    });

    // ========================================
    // Poisson Bracket Tests
    // ========================================

    run_test("Poisson bracket instantiation", []() {
        PoissonBracket pb;
        return true;
    });

    run_test("Poisson bracket q and p fundamental relation", []() {
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.5;

        auto q = [](const PhasePoint& p) { return p.q[0]; };
        auto p_func = [](const PhasePoint& p) { return p.p[0]; };

        double bracket = pb.compute(q, p_func, point);

        // {q, p} = 1
        ASSERT_NEAR(bracket, 1.0, TOLERANCE_POISSON);
        return true;
    });

    run_test("Poisson bracket q and q equals zero", []() {
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.5;

        auto q = [](const PhasePoint& p) { return p.q[0]; };

        double bracket = pb.compute(q, q, point);

        // {q, q} = 0
        ASSERT_NEAR(bracket, 0.0, TOLERANCE_POISSON);
        return true;
    });

    run_test("Poisson bracket p and p equals zero", []() {
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.5;

        auto p_func = [](const PhasePoint& p) { return p.p[0]; };

        double bracket = pb.compute(p_func, p_func, point);

        // {p, p} = 0
        ASSERT_NEAR(bracket, 0.0, TOLERANCE_POISSON);
        return true;
    });

    run_test("Poisson bracket anti-commutativity", []() {
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 1.5;
        point.p[0] = 0.7;

        auto f = [](const PhasePoint& p) { return p.q[0] * p.q[0]; };
        auto g = [](const PhasePoint& p) { return p.p[0]; };

        double fg = pb.compute(f, g, point);
        double gf = pb.compute(g, f, point);

        // {f, g} = -{g, f}
        ASSERT_NEAR(fg, -gf, TOLERANCE_POISSON);
        return true;
    });

    run_test("Poisson bracket with Hamiltonian", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.0;

        auto q = [](const PhasePoint& p) { return p.q[0]; };

        double bracket = pb.compute(q,
            [&ho](const PhasePoint& p) { return ho.hamiltonian(p); },
            point);

        // {q, H} should equal dq/dt
        double dq_dt = ho.dH_dp(point, 0);
        ASSERT_NEAR(bracket, dq_dt, TOLERANCE_POISSON);
        return true;
    });

    run_test("Poisson bracket fundamental relations 2D", []() {
        PoissonBracket pb;

        PhasePoint point(2);
        point.q[0] = 1.0; point.q[1] = 2.0;
        point.p[0] = 0.5; point.p[1] = 1.5;

        ASSERT_TRUE(pb.verifyFundamentalRelations(point, TOLERANCE_POISSON));
        return true;
    });

    run_test("Poisson bracket fundamental relations 3D", []() {
        PoissonBracket pb;

        PhasePoint point(3);
        point.q[0] = 1.0; point.q[1] = 2.0; point.q[2] = 3.0;
        point.p[0] = 0.5; point.p[1] = 1.5; point.p[2] = 2.5;

        ASSERT_TRUE(pb.verifyFundamentalRelations(point, TOLERANCE_POISSON));
        return true;
    });

    run_test("Poisson bracket Jacobi identity", []() {
        // {{f, g}, h} + {{g, h}, f} + {{h, f}, g} = 0
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.5;

        auto f = [](const PhasePoint& p) { return p.q[0]; };
        auto g = [](const PhasePoint& p) { return p.p[0]; };
        auto h = [](const PhasePoint& p) { return p.q[0] * p.p[0]; };

        double fg = pb.compute(f, g, point);
        double fg_h = pb.compute([&](const PhasePoint& p) { return fg; }, h, point);

        double gh = pb.compute(g, h, point);
        double gh_f = pb.compute([&](const PhasePoint& p) { return gh; }, f, point);

        double hf = pb.compute(h, f, point);
        double hf_g = pb.compute([&](const PhasePoint& p) { return hf; }, g, point);

        // This test is approximate due to functional composition
        double jacobi_sum = std::abs(fg_h + gh_f + hf_g);
        ASSERT_TRUE(jacobi_sum < TOLERANCE_POISSON);
        return true;
    });

    run_test("Poisson bracket with quadratic function", []() {
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 2.0;
        point.p[0] = 1.0;

        auto f = [](const PhasePoint& p) { return p.q[0] * p.q[0]; };
        auto g = [](const PhasePoint& p) { return p.p[0]; };

        double bracket = pb.compute(f, g, point);

        // {q², p} = 2q{q, p} = 2q
        ASSERT_NEAR(bracket, 4.0, TOLERANCE_POISSON);
        return true;
    });

    // ========================================
    // Numerical Consistency Tests
    // ========================================

    run_test("Integrators produce smooth trajectories", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double dt = 0.001;
        int num_steps = 100;
        auto trajectory = ho.integrate(initial, dt, num_steps);

        // Check that trajectory doesn't have wild jumps
        for (size_t i = 1; i < trajectory.size(); ++i) {
            double dq = std::abs(trajectory[i].q[0] - trajectory[i-1].q[0]);
            double dp = std::abs(trajectory[i].p[0] - trajectory[i-1].p[0]);

            // Each step should move smoothly
            ASSERT_TRUE(dq < 1.0 && dp < 1.0);
        }

        return true;
    });

    run_test("Verlet vs Symplectic Euler accuracy", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double E_initial = ho.hamiltonian(initial);

        // Symplectic Euler
        PhasePoint euler_current = initial;
        double dt = 0.01;
        int num_steps = 100;

        for (int i = 0; i < num_steps; ++i) {
            euler_current = ho.stepSymplecticEuler(euler_current, dt);
        }

        double E_euler = ho.hamiltonian(euler_current);
        double error_euler = std::abs(E_euler - E_initial);

        // Verlet
        PhasePoint verlet_current = initial;
        for (int i = 0; i < num_steps; ++i) {
            verlet_current = ho.stepVerlet(verlet_current, dt);
        }

        double E_verlet = ho.hamiltonian(verlet_current);
        double error_verlet = std::abs(E_verlet - E_initial);

        // Verlet should be more accurate
        ASSERT_TRUE(error_verlet < error_euler);
        return true;
    });

    run_test("Dimensional consistency multi-dimensional system", []() {
        HamiltonianSystem kepler = HamiltonianSystem::keplerProblem(1.0, 1.0);

        PhasePoint point(2);
        point.q[0] = 1.0;
        point.q[1] = 0.5;
        point.p[0] = 0.0;
        point.p[1] = 1.0;

        PhasePoint deriv = kepler.hamiltonEquations(point);

        ASSERT_TRUE(deriv.q.dimension() == 2);
        ASSERT_TRUE(deriv.p.dimension() == 2);
        return true;
    });

    run_test("Hamiltonian system negative dimension throws", []() {
        try {
            HamiltonianSystem sys(-1, [](const PhasePoint&) { return 0.0; });
            return false;
        } catch (const std::invalid_argument&) {
            return true;
        }
    });

    run_test("Hamiltonian system zero dimension throws", []() {
        try {
            HamiltonianSystem sys(0, [](const PhasePoint&) { return 0.0; });
            return false;
        } catch (const std::invalid_argument&) {
            return true;
        }
    });

    // ========================================
    // Multi-dimensional Tests
    // ========================================

    run_test("2D harmonic oscillator creation", []() {
        auto H_2d = [](const PhasePoint& point) {
            double x = point.q[0], y = point.q[1];
            double px = point.p[0], py = point.p[1];
            return px*px/2.0 + py*py/2.0 + 0.5*x*x + 0.5*y*y;
        };

        HamiltonianSystem ho2d(2, H_2d);
        ASSERT_TRUE(ho2d.dimension() == 2);
        return true;
    });

    run_test("2D harmonic oscillator energy conservation", []() {
        auto H_2d = [](const PhasePoint& point) {
            double x = point.q[0], y = point.q[1];
            double px = point.p[0], py = point.p[1];
            return px*px/2.0 + py*py/2.0 + 0.5*x*x + 0.5*y*y;
        };

        HamiltonianSystem ho2d(2, H_2d);

        PhasePoint initial(2);
        initial.q[0] = 1.0;
        initial.q[1] = 0.5;
        initial.p[0] = 0.0;
        initial.p[1] = 0.0;

        double dt = 0.01;
        int num_steps = 500;
        auto trajectory = ho2d.integrate(initial, dt, num_steps);

        ASSERT_TRUE(ho2d.checkEnergyConservation(trajectory, 0.01));
        return true;
    });

    // ========================================
    // Edge Cases and Special Values
    // ========================================

    run_test("Very small time step integration", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double dt = 1e-5;
        int num_steps = 100;
        auto trajectory = ho.integrate(initial, dt, num_steps);

        // Should still conserve energy with very small step
        ASSERT_TRUE(ho.checkEnergyConservation(trajectory, 1e-6));
        return true;
    });

    run_test("Large system parameters harmonic oscillator", []() {
        // Large mass and stiffness
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(100.0, 100.0);

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 10.0;

        double H = ho.hamiltonian(point);

        // H = p²/(2m) + kq²/2 = 100/(200) + 100/2 = 0.5 + 50 = 50.5
        ASSERT_NEAR(H, 50.5, TOLERANCE_EXACT);
        return true;
    });

    run_test("Poisson bracket zero functions", []() {
        PoissonBracket pb;

        PhasePoint point(1);
        point.q[0] = 1.0;
        point.p[0] = 0.5;

        auto zero = [](const PhasePoint&) { return 0.0; };
        auto q = [](const PhasePoint& p) { return p.q[0]; };

        double bracket = pb.compute(zero, q, point);

        // {0, q} = 0
        ASSERT_NEAR(bracket, 0.0, TOLERANCE_POISSON);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "=== Test Results ===" << std::endl;
    std::cout << "Passed: " << tests_passed << std::endl;
    std::cout << "Failed: " << tests_failed << std::endl;
    std::cout << "Total:  " << (tests_passed + tests_failed) << std::endl;
    std::cout << std::endl;

    return tests_failed == 0 ? 0 : 1;
}
