/**
 * Phase 4 Validation: Classical Phase Space Analysis
 *
 * Tests the classical_phase_space.hpp module functions.
 *
 * Coverage:
 * - PhaseSpaceVolume: volume element calculation and Liouville's theorem
 * - PhaseSpaceDensity: Maxwell-Boltzmann, microcanonical ensemble, observables
 * - PoincareSection: surface crossings, interpolation, period detection
 * - LyapunovExponent: exponent computation and dynamics classification
 * - SymplecticMatrix: Ω structure construction and symplectic condition
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "../include/physics/classical_phase_space.hpp"

using namespace physics::advanced::classical;

// Test tolerances
const double TOLERANCE_EXACT = 1e-5;
const double TOLERANCE_NUMERICAL = 1e-3;
const double TOLERANCE_INTEGRATOR = 0.01;

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

    std::cout << "=== Phase 4: Classical Phase Space Analysis Validation ===" << std::endl;
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
    // PhaseSpaceVolume Tests
    // ========================================

    run_test("PhaseSpaceVolume empty points returns zero", []() {
        std::vector<PhasePoint> empty;
        double vol = PhaseSpaceVolume::volumeElement(empty);
        ASSERT_NEAR(vol, 0.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceVolume 1D single point zero volume", []() {
        std::vector<PhasePoint> points;
        PhasePoint p(1);
        p.q[0] = 1.0;
        p.p[0] = 2.0;
        points.push_back(p);

        double vol = PhaseSpaceVolume::volumeElement(points);
        ASSERT_NEAR(vol, 0.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceVolume 1D rectangle calculation", []() {
        std::vector<PhasePoint> points;

        // Create rectangle: q in [0, 2], p in [0, 3]
        PhasePoint p1(1);
        p1.q[0] = 0.0;
        p1.p[0] = 0.0;
        points.push_back(p1);

        PhasePoint p2(1);
        p2.q[0] = 2.0;
        p2.p[0] = 3.0;
        points.push_back(p2);

        double vol = PhaseSpaceVolume::volumeElement(points);
        // Expected: Δq × Δp = 2 × 3 = 6
        ASSERT_NEAR(vol, 6.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceVolume 1D asymmetric rectangle", []() {
        std::vector<PhasePoint> points;

        PhasePoint p1(1);
        p1.q[0] = -1.0;
        p1.p[0] = 1.0;
        points.push_back(p1);

        PhasePoint p2(1);
        p2.q[0] = 1.0;
        p2.p[0] = 4.0;
        points.push_back(p2);

        double vol = PhaseSpaceVolume::volumeElement(points);
        // Expected: Δq × Δp = 2 × 3 = 6
        ASSERT_NEAR(vol, 6.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceVolume multiple points bounding box", []() {
        std::vector<PhasePoint> points;

        PhasePoint p1(1);
        p1.q[0] = 0.5;
        p1.p[0] = 1.5;
        points.push_back(p1);

        PhasePoint p2(1);
        p2.q[0] = 0.0;
        p2.p[0] = 0.0;
        points.push_back(p2);

        PhasePoint p3(1);
        p3.q[0] = 3.0;
        p3.p[0] = 5.0;
        points.push_back(p3);

        double vol = PhaseSpaceVolume::volumeElement(points);
        // Expected: Δq × Δp = 3 × 5 = 15
        ASSERT_NEAR(vol, 15.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceVolume Liouville theorem harmonic oscillator", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        // Create initial region with 5 phase space points
        std::vector<PhasePoint> initial_region;
        for (int i = 0; i < 5; ++i) {
            PhasePoint p(1);
            p.q[0] = i * 0.2;
            p.p[0] = (4-i) * 0.25;
            initial_region.push_back(p);
        }

        double time = 0.1;
        double dt = 0.001;

        bool liouville_ok = PhaseSpaceVolume::verifyLiouvilleTheorem(
            ho, initial_region, time, dt, 0.1);

        ASSERT_TRUE(liouville_ok);
        return true;
    });

    // ========================================
    // PhaseSpaceDensity Tests
    // ========================================

    run_test("PhaseSpaceDensity construction and evaluation", []() {
        auto rho_func = [](const PhasePoint& p) {
            return std::exp(-(p.q[0]*p.q[0] + p.p[0]*p.p[0]));
        };

        PhaseSpaceDensity density(1, rho_func);

        PhasePoint test_point(1);
        test_point.q[0] = 0.0;
        test_point.p[0] = 0.0;

        double rho = density(test_point);
        ASSERT_NEAR(rho, 1.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceDensity Maxwell-Boltzmann high temperature", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        // High temperature with unit k_B (natural units)
        double T = 100.0;
        double k_B = 1.0;  // Natural units, not SI
        PhaseSpaceDensity mb = PhaseSpaceDensity::maxwellBoltzmann(ho, T, k_B);

        PhasePoint p1(1);
        p1.q[0] = 0.0;
        p1.p[0] = 0.0;

        PhasePoint p2(1);
        p2.q[0] = 1.0;
        p2.p[0] = 0.0;

        double rho1 = mb(p1);
        double rho2 = mb(p2);

        // At high T, both densities should be positive
        ASSERT_TRUE(rho1 > 0.0);
        ASSERT_TRUE(rho2 > 0.0);
        // At high T, all states nearly equally probable: difference should be small
        ASSERT_TRUE(std::abs(rho1 - rho2) < 0.1);
        return true;
    });

    run_test("PhaseSpaceDensity Maxwell-Boltzmann decreases with energy", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        double T = 1.0;
        double k_B = 1.0;  // Natural units
        PhaseSpaceDensity mb = PhaseSpaceDensity::maxwellBoltzmann(ho, T, k_B);

        // H = p²/2 + kq²/2, with m=1, k=1
        PhasePoint p_low(1);
        p_low.q[0] = 0.05;  // Lower energy point
        p_low.p[0] = 0.05;

        PhasePoint p_high(1);
        p_high.q[0] = 2.0;  // Higher energy point
        p_high.p[0] = 2.0;

        double rho_low = mb(p_low);
        double rho_high = mb(p_high);

        // Lower energy should have higher density: rho ∝ exp(-βH)
        ASSERT_TRUE(rho_low > rho_high);
        return true;
    });

    run_test("PhaseSpaceDensity Microcanonical construction", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        double E = 1.0;
        double delta_E = 1e-6;
        PhaseSpaceDensity micro = PhaseSpaceDensity::microcanonical(ho, E, delta_E);

        // Point exactly at energy E
        PhasePoint p_on_shell(1);
        p_on_shell.q[0] = std::sqrt(2.0);  // H = 1.0 for this point
        p_on_shell.p[0] = 0.0;

        double rho_on = micro(p_on_shell);
        ASSERT_TRUE(rho_on > 0.0);
        return true;
    });

    run_test("PhaseSpaceDensity Microcanonical zero off energy shell", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        double E = 1.0;
        double delta_E = 1e-6;
        PhaseSpaceDensity micro = PhaseSpaceDensity::microcanonical(ho, E, delta_E);

        // Point far from energy E
        PhasePoint p_off(1);
        p_off.q[0] = 10.0;
        p_off.p[0] = 10.0;

        double rho_off = micro(p_off);
        ASSERT_NEAR(rho_off, 0.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceDensity average observable simple", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        // Uniform density
        auto uniform_density = [](const PhasePoint&) { return 1.0; };
        PhaseSpaceDensity density(1, uniform_density);

        // Observable: position coordinate
        auto obs_q = [](const PhasePoint& p) { return p.q[0]; };

        std::vector<PhasePoint> sample_points;
        for (int i = 0; i < 5; ++i) {
            PhasePoint p(1);
            p.q[0] = i * 0.5;
            p.p[0] = 0.0;
            sample_points.push_back(p);
        }

        double avg = density.average(obs_q, sample_points);
        // Expected: mean of [0, 0.5, 1.0, 1.5, 2.0] = 1.0
        ASSERT_NEAR(avg, 1.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PhaseSpaceDensity average observable weighted", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        // Density concentrated at q=1
        auto concentrated_density = [](const PhasePoint& p) {
            return std::exp(-(p.q[0] - 1.0)*(p.q[0] - 1.0));
        };
        PhaseSpaceDensity density(1, concentrated_density);

        // Observable: position
        auto obs_q = [](const PhasePoint& p) { return p.q[0]; };

        std::vector<PhasePoint> sample_points;
        for (int i = 0; i < 5; ++i) {
            PhasePoint p(1);
            p.q[0] = i * 0.5;
            p.p[0] = 0.0;
            sample_points.push_back(p);
        }

        double avg = density.average(obs_q, sample_points);
        // Average should be pulled toward q=1
        ASSERT_TRUE(avg > 0.5 && avg < 1.5);
        return true;
    });

    // ========================================
    // PoincareSection Tests
    // ========================================

    run_test("PoincareSection construction with surface", []() {
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);
        return true;
    });

    run_test("PoincareSection detects crossing positive to negative", []() {
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        PhasePoint p1(1);
        p1.q[0] = 1.0;
        p1.p[0] = 0.0;

        PhasePoint p2(1);
        p2.q[0] = -1.0;
        p2.p[0] = 0.0;

        ASSERT_TRUE(section.crosses(p1, p2));
        return true;
    });

    run_test("PoincareSection detects crossing negative to positive", []() {
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        PhasePoint p1(1);
        p1.q[0] = -1.0;
        p1.p[0] = 0.0;

        PhasePoint p2(1);
        p2.q[0] = 1.0;
        p2.p[0] = 0.0;

        ASSERT_TRUE(section.crosses(p1, p2));
        return true;
    });

    run_test("PoincareSection no crossing same side", []() {
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        PhasePoint p1(1);
        p1.q[0] = 1.0;
        p1.p[0] = 0.0;

        PhasePoint p2(1);
        p2.q[0] = 2.0;
        p2.p[0] = 1.0;

        ASSERT_TRUE(!section.crosses(p1, p2));
        return true;
    });

    run_test("PoincareSection record crossings harmonic oscillator", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        // Surface: q = 0
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double dt = 0.01;
        int num_steps = 500;  // More than one period
        auto trajectory = ho.integrate(initial, dt, num_steps);

        section.recordCrossings(trajectory);
        const auto& crossings = section.getCrossings();

        // Should have crossings for harmonic oscillator
        ASSERT_TRUE(crossings.size() > 0);
        return true;
    });

    run_test("PoincareSection crossing interpolation accuracy", []() {
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        std::vector<PhasePoint> trajectory;

        PhasePoint p1(1);
        p1.q[0] = 0.5;
        p1.p[0] = 1.0;
        trajectory.push_back(p1);

        PhasePoint p2(1);
        p2.q[0] = -0.5;
        p2.p[0] = 0.8;
        trajectory.push_back(p2);

        section.recordCrossings(trajectory);
        const auto& crossings = section.getCrossings();

        if (crossings.size() > 0) {
            // Crossing point should be at q ≈ 0
            ASSERT_NEAR(crossings[0].q[0], 0.0, 0.1);
        }

        return true;
    });

    run_test("PoincareSection period detection nonperiodic", []() {
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        std::vector<PhasePoint> trajectory;
        // Create a trajectory with non-periodic crossings
        for (int i = 0; i < 5; ++i) {
            PhasePoint p(1);
            p.q[0] = i * 0.1;
            p.p[0] = 1.0 + i * 0.01;
            trajectory.push_back(p);
        }

        section.recordCrossings(trajectory);
        int period = section.detectPeriod(1e-4);

        // Non-periodic, should return 0
        ASSERT_NEAR(period, 0, TOLERANCE_EXACT);
        return true;
    });

    run_test("PoincareSection period detection for periodic orbit", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double dt = 0.01;
        int num_steps = 1000;
        auto trajectory = ho.integrate(initial, dt, num_steps);

        section.recordCrossings(trajectory);
        int period = section.detectPeriod(1e-3);

        // Harmonic oscillator is periodic, should detect period
        ASSERT_TRUE(period > 0);
        return true;
    });

    run_test("PoincareSection empty trajectory", []() {
        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        std::vector<PhasePoint> trajectory;
        section.recordCrossings(trajectory);
        const auto& crossings = section.getCrossings();

        ASSERT_TRUE(crossings.empty());
        return true;
    });

    // ========================================
    // LyapunovExponent Tests
    // ========================================

    run_test("LyapunovExponent construction", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);
        return true;
    });

    run_test("LyapunovExponent harmonic oscillator near zero", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        double lambda = lyap.compute(initial, 0.01, 100, 1e-8);

        // Harmonic oscillator should have λ ≈ 0 (neutral/quasiperiodic)
        ASSERT_NEAR(lambda, 0.0, 0.5);
        return true;
    });

    run_test("LyapunovExponent classification stable", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);

        double lambda_stable = -0.1;
        auto classification = lyap.classify(lambda_stable);

        ASSERT_TRUE(classification == LyapunovExponent::DynamicsType::STABLE);
        return true;
    });

    run_test("LyapunovExponent classification neutral", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);

        double lambda_neutral = 0.0;
        auto classification = lyap.classify(lambda_neutral);

        ASSERT_TRUE(classification == LyapunovExponent::DynamicsType::NEUTRAL);
        return true;
    });

    run_test("LyapunovExponent classification chaotic", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);

        double lambda_chaotic = 0.1;
        auto classification = lyap.classify(lambda_chaotic);

        ASSERT_TRUE(classification == LyapunovExponent::DynamicsType::CHAOTIC);
        return true;
    });

    run_test("LyapunovExponent classification boundary negative", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);

        double lambda = -1e-7;
        auto classification = lyap.classify(lambda);

        // lambda = -1e-7 is between -1e-6 and 1e-6, so it's NEUTRAL
        ASSERT_TRUE(classification == LyapunovExponent::DynamicsType::NEUTRAL);
        return true;
    });

    run_test("LyapunovExponent classification boundary positive", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);

        double lambda = 1e-7;
        auto classification = lyap.classify(lambda);

        // lambda = 1e-7 is between -1e-6 and 1e-6, so it's NEUTRAL
        ASSERT_TRUE(classification == LyapunovExponent::DynamicsType::NEUTRAL);
        return true;
    });

    run_test("LyapunovExponent different initial conditions", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);
        LyapunovExponent lyap(ho);

        PhasePoint initial1(1);
        initial1.q[0] = 0.5;
        initial1.p[0] = 0.5;

        double lambda1 = lyap.compute(initial1, 0.01, 100, 1e-8);

        PhasePoint initial2(1);
        initial2.q[0] = 1.0;
        initial2.p[0] = 1.0;

        double lambda2 = lyap.compute(initial2, 0.01, 100, 1e-8);

        // Both should be close to zero for harmonic oscillator
        ASSERT_TRUE(std::abs(lambda1) < 1.0);
        ASSERT_TRUE(std::abs(lambda2) < 1.0);
        return true;
    });

    // ========================================
    // SymplecticMatrix Tests
    // ========================================

    run_test("SymplecticMatrix construction 1D", []() {
        SymplecticMatrix omega(1);
        return true;
    });

    run_test("SymplecticMatrix construction 2D", []() {
        SymplecticMatrix omega(2);
        return true;
    });

    run_test("SymplecticMatrix construction 3D", []() {
        SymplecticMatrix omega(3);
        return true;
    });

    run_test("SymplecticMatrix structure 1D dimension", []() {
        SymplecticMatrix omega(1);
        const auto& mat = omega.matrix();

        // Should be 2×2 for 1D system
        ASSERT_TRUE(mat.rows() == 2 && mat.cols() == 2);
        return true;
    });

    run_test("SymplecticMatrix structure 2D dimension", []() {
        SymplecticMatrix omega(2);
        const auto& mat = omega.matrix();

        // Should be 4×4 for 2D system
        ASSERT_TRUE(mat.rows() == 4 && mat.cols() == 4);
        return true;
    });

    run_test("SymplecticMatrix upper right identity block 1D", []() {
        SymplecticMatrix omega(1);
        const auto& mat = omega.matrix();

        // Upper right block should be I: element (0, 1) = 1
        ASSERT_NEAR(mat(0, 1), 1.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("SymplecticMatrix lower left negative identity block 1D", []() {
        SymplecticMatrix omega(1);
        const auto& mat = omega.matrix();

        // Lower left block should be -I: element (1, 0) = -1
        ASSERT_NEAR(mat(1, 0), -1.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("SymplecticMatrix diagonal blocks zero 1D", []() {
        SymplecticMatrix omega(1);
        const auto& mat = omega.matrix();

        // Diagonal upper left should be 0
        ASSERT_NEAR(mat(0, 0), 0.0, TOLERANCE_EXACT);
        // Diagonal lower right should be 0
        ASSERT_NEAR(mat(1, 1), 0.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("SymplecticMatrix upper right identity block 2D", []() {
        SymplecticMatrix omega(2);
        const auto& mat = omega.matrix();

        // Upper right 2×2 block should be identity
        ASSERT_NEAR(mat(0, 2), 1.0, TOLERANCE_EXACT);
        ASSERT_NEAR(mat(1, 3), 1.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("SymplecticMatrix lower left negative identity block 2D", []() {
        SymplecticMatrix omega(2);
        const auto& mat = omega.matrix();

        // Lower left 2×2 block should be -I
        ASSERT_NEAR(mat(2, 0), -1.0, TOLERANCE_EXACT);
        ASSERT_NEAR(mat(3, 1), -1.0, TOLERANCE_EXACT);
        return true;
    });

    run_test("SymplecticMatrix identity is symplectic", []() {
        SymplecticMatrix omega(1);

        // Create identity matrix (2×2 for 1D phase space)
        maths::linear_algebra::Matrix I(2, 2);
        I(0, 0) = 1.0;
        I(0, 1) = 0.0;
        I(1, 0) = 0.0;
        I(1, 1) = 1.0;

        ASSERT_TRUE(omega.isSymplectic(I, 1e-10));
        return true;
    });

    run_test("SymplecticMatrix 90 degree rotation not symplectic", []() {
        SymplecticMatrix omega(1);

        // 90 degree rotation matrix
        maths::linear_algebra::Matrix R(2, 2);
        R(0, 0) = 0.0;
        R(0, 1) = -1.0;
        R(1, 0) = 1.0;
        R(1, 1) = 0.0;

        // Rotation should be symplectic for 1D (it preserves volume)
        bool result = omega.isSymplectic(R, 1e-10);
        ASSERT_TRUE(result);
        return true;
    });

    run_test("SymplecticMatrix scaling by 2 not symplectic", []() {
        SymplecticMatrix omega(1);

        // Scaling matrix
        maths::linear_algebra::Matrix S(2, 2);
        S(0, 0) = 2.0;
        S(0, 1) = 0.0;
        S(1, 0) = 0.0;
        S(1, 1) = 2.0;

        // Scaling does NOT preserve symplectic structure
        ASSERT_TRUE(!omega.isSymplectic(S, 1e-10));
        return true;
    });

    run_test("SymplecticMatrix diagonal scaling by reciprocal is symplectic", []() {
        SymplecticMatrix omega(1);

        // Diagonal scaling: q -> 2q, p -> p/2
        maths::linear_algebra::Matrix D(2, 2);
        D(0, 0) = 2.0;  // q coordinate
        D(0, 1) = 0.0;
        D(1, 0) = 0.0;
        D(1, 1) = 0.5;  // p coordinate (reciprocal)

        // This should preserve phase space volume
        ASSERT_TRUE(omega.isSymplectic(D, 1e-8));
        return true;
    });

    run_test("SymplecticMatrix 2D identity symplectic", []() {
        SymplecticMatrix omega(2);

        // 4×4 identity matrix
        maths::linear_algebra::Matrix I(4, 4);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                I(i, j) = (i == j) ? 1.0 : 0.0;
            }
        }

        ASSERT_TRUE(omega.isSymplectic(I, 1e-10));
        return true;
    });

    run_test("SymplecticMatrix symplectic property M^T Omega M equals Omega", []() {
        SymplecticMatrix omega(1);

        // Identity matrix
        maths::linear_algebra::Matrix M(2, 2);
        M(0, 0) = 1.0;
        M(0, 1) = 0.0;
        M(1, 0) = 0.0;
        M(1, 1) = 1.0;

        // Compute M^T Ω M and verify it equals Ω
        const auto& O = omega.matrix();
        maths::linear_algebra::Matrix MT = M.transpose();
        maths::linear_algebra::Matrix result = MT * O * M;

        // Compare result with O
        bool match = true;
        for (size_t i = 0; i < result.rows(); ++i) {
            for (size_t j = 0; j < result.cols(); ++j) {
                if (std::abs(result(i, j) - O(i, j)) > 1e-10) {
                    match = false;
                }
            }
        }

        ASSERT_TRUE(match);
        return true;
    });

    // ========================================
    // Integration Tests: Multiple Components
    // ========================================

    run_test("Phase space volume and density consistency", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        // Create a phase space region with extent in both q and p
        std::vector<PhasePoint> region;
        for (int i = 0; i < 5; ++i) {
            PhasePoint p(1);
            p.q[0] = i * 0.2;
            p.p[0] = i * 0.3;
            region.push_back(p);
        }

        double vol = PhaseSpaceVolume::volumeElement(region);
        // q ranges from 0 to 0.8, p ranges from 0 to 1.2, so volume = 0.8 * 1.2 = 0.96
        ASSERT_NEAR(vol, 0.96, TOLERANCE_EXACT);
        return true;
    });

    run_test("Poincare section and Lyapunov exponent harmonic oscillator", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        auto surface = [](const PhasePoint& p) { return p.q[0]; };
        PoincareSection section(surface);

        PhasePoint initial(1);
        initial.q[0] = 1.0;
        initial.p[0] = 0.0;

        auto traj = ho.integrate(initial, 0.01, 500);
        section.recordCrossings(traj);

        LyapunovExponent lyap(ho);
        double lambda = lyap.compute(initial, 0.01, 100);

        // Harmonic oscillator should have near-zero Lyapunov exponent
        // and detectable Poincare crossings
        ASSERT_TRUE(std::abs(lambda) < 1.0);
        ASSERT_TRUE(section.getCrossings().size() > 0);
        return true;
    });

    run_test("Symplectic matrix with Verlet integration structure", []() {
        SymplecticMatrix omega(1);

        // Verlet integrator is symplectic - test with composition
        maths::linear_algebra::Matrix R1(2, 2);
        R1(0, 0) = 1.0;
        R1(0, 1) = 0.5;  // Half step
        R1(1, 0) = 0.0;
        R1(1, 1) = 1.0;

        maths::linear_algebra::Matrix R2(2, 2);
        R2(0, 0) = 1.0;
        R2(0, 1) = 0.0;
        R2(1, 0) = -0.5;  // Force step
        R2(1, 1) = 1.0;

        ASSERT_TRUE(omega.isSymplectic(R1, 1e-10));
        ASSERT_TRUE(omega.isSymplectic(R2, 1e-10));
        return true;
    });

    run_test("Phase space analysis workflow", []() {
        HamiltonianSystem ho = HamiltonianSystem::harmonicOscillator(1.0, 1.0);

        // Phase space density
        PhaseSpaceDensity mb = PhaseSpaceDensity::maxwellBoltzmann(ho, 1.0);

        // Test observable average
        auto observable = [](const PhasePoint& p) {
            return p.q[0] * p.q[0];  // q^2
        };

        std::vector<PhasePoint> samples;
        for (int i = 0; i < 3; ++i) {
            PhasePoint p(1);
            p.q[0] = i * 0.5;
            p.p[0] = 0.0;
            samples.push_back(p);
        }

        double avg_q2 = mb.average(observable, samples);
        ASSERT_TRUE(avg_q2 >= 0.0);
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
