/**
 * Phase 4 Validation: Fluid Dynamics Turbulence
 *
 * Tests the fluid_dynamics_turbulence.hpp module functions.
 *
 * Coverage:
 * - Reynolds decomposition and averaging
 * - Reynolds stress tensor
 * - Turbulent kinetic energy and intensity
 * - k-epsilon turbulence model
 * - Mixing length theory
 * - Boussinesq hypothesis
 * - Kolmogorov scales
 * - Wall functions
 */

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "../include/physics/fluid_dynamics_turbulence.hpp"

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

    std::cout << "=== Phase 4: Fluid Dynamics Turbulence Validation ===" << std::endl;
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
    // Reynolds Decomposition Tests
    // ========================================

    run_test("Reynolds decomposition: create decomposed field", []() {
        maths::linear_algebra::Vector instantaneous(3);
        instantaneous[0] = 10.5;
        instantaneous[1] = 5.2;
        instantaneous[2] = 3.1;

        maths::linear_algebra::Vector mean(3);
        mean[0] = 10.0;
        mean[1] = 5.0;
        mean[2] = 3.0;

        auto decomposed = ReynoldsDecomposition::decompose(instantaneous, mean);

        ASSERT_NEAR(decomposed.mean[0], 10.0, TOLERANCE);
        ASSERT_NEAR(decomposed.fluctuation[0], 0.5, TOLERANCE);
        ASSERT_NEAR(decomposed.fluctuation[1], 0.2, TOLERANCE);
        return true;
    });

    run_test("Reynolds decomposition: verify averaging of fluctuations", []() {
        std::vector<maths::linear_algebra::Vector> fluctuations;
        for (int i = 0; i < 100; ++i) {
            maths::linear_algebra::Vector u_prime(3);
            u_prime[0] = std::sin(i * 0.1);
            u_prime[1] = std::cos(i * 0.1);
            u_prime[2] = std::sin(i * 0.2);
            fluctuations.push_back(u_prime);
        }

        // Balanced fluctuations should average to ~zero
        bool verified = ReynoldsDecomposition::verifyAveraging(fluctuations, 0.5);
        ASSERT_TRUE(verified);
        return true;
    });

    run_test("Reynolds decomposition: fluctuation orthogonality", []() {
        maths::linear_algebra::Vector mean(3);
        mean[0] = 10.0;
        mean[1] = 0.0;
        mean[2] = 0.0;

        maths::linear_algebra::Vector inst(3);
        inst[0] = 10.5;
        inst[1] = 0.0;
        inst[2] = 0.0;

        auto decomp = ReynoldsDecomposition::decompose(inst, mean);

        // Mean and fluctuation reconstruction
        double recon = decomp.mean[0] + decomp.fluctuation[0];
        ASSERT_NEAR(recon, inst[0], TOLERANCE);
        return true;
    });

    // ========================================
    // Reynolds Stress Tensor Tests
    // ========================================

    run_test("Reynolds stress: compute from fluctuations", []() {
        std::vector<maths::linear_algebra::Vector> fluctuations;

        // Create simple fluctuation pattern
        for (int i = 0; i < 10; ++i) {
            maths::linear_algebra::Vector u(3);
            u[0] = (i % 2 == 0) ? 1.0 : -1.0;
            u[1] = 0.0;
            u[2] = 0.0;
            fluctuations.push_back(u);
        }

        double rho = 1.225;
        auto stress = ReynoldsStressTensor::compute(rho, fluctuations);

        // Diagonal term should be -rho * <u'^2>
        double expected = -rho * 1.0;  // <u'^2> = 1
        ASSERT_NEAR(stress(0,0), expected, NUMERICAL_TOLERANCE);
        return true;
    });

    run_test("Reynolds stress: turbulent kinetic energy", []() {
        std::vector<maths::linear_algebra::Vector> fluctuations;

        // Unit fluctuations
        maths::linear_algebra::Vector u(3);
        u[0] = 1.0; u[1] = 1.0; u[2] = 1.0;
        fluctuations.push_back(u);

        double k = ReynoldsStressTensor::turbulentKineticEnergy(fluctuations);

        // k = 0.5 * (1^2 + 1^2 + 1^2) = 1.5
        double expected = 0.5 * 3.0;
        ASSERT_NEAR(k, expected, TOLERANCE);
        return true;
    });

    run_test("Reynolds stress: turbulence intensity", []() {
        double k = 1.0;  // m^2/s^2
        double U = 10.0;  // m/s

        double I = ReynoldsStressTensor::turbulenceIntensity(k, U);

        // I = sqrt(2k/3) / U = sqrt(2/3) / 10
        double expected = std::sqrt(2.0/3.0) / 10.0;
        ASSERT_NEAR(I, expected, TOLERANCE);
        return true;
    });

    run_test("Reynolds stress: low turbulence intensity", []() {
        double k = 0.01;
        double U = 20.0;
        double I = ReynoldsStressTensor::turbulenceIntensity(k, U);

        // I < 1% is low turbulence
        ASSERT_TRUE(I < 0.01);
        return true;
    });

    run_test("Reynolds stress: high turbulence intensity", []() {
        double k = 5.0;
        double U = 10.0;
        double I = ReynoldsStressTensor::turbulenceIntensity(k, U);

        // I > 5% is high turbulence
        ASSERT_TRUE(I > 0.05);
        return true;
    });

    run_test("Reynolds stress: anisotropy tensor for isotropic turbulence", []() {
        // For isotropic turbulence, Reynolds stress is diagonal with equal components
        maths::linear_algebra::Matrix R(3, 3);
        double k = 1.0;

        // Normalized Reynolds stress for isotropic case
        R(0,0) = 2.0*k/3.0; R(1,1) = 2.0*k/3.0; R(2,2) = 2.0*k/3.0;

        auto b = ReynoldsStressTensor::anisotropyTensor(R, k);

        // For isotropic turbulence, b_ij should be ~0
        ASSERT_NEAR(b(0,0), 0.0, TOLERANCE);
        ASSERT_NEAR(b(1,1), 0.0, TOLERANCE);
        ASSERT_NEAR(b(2,2), 0.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Boussinesq Hypothesis Tests
    // ========================================

    run_test("Boussinesq: strain rate tensor symmetry", []() {
        maths::linear_algebra::Matrix vel_grad(3, 3);
        vel_grad(0,1) = 2.0;
        vel_grad(1,0) = 1.0;

        auto S = BoussinesqHypothesis::strainRateTensor(vel_grad);

        // Should be symmetric: S_ij = S_ji
        ASSERT_NEAR(S(0,1), S(1,0), TOLERANCE);
        return true;
    });

    run_test("Boussinesq: strain rate magnitude", []() {
        maths::linear_algebra::Matrix S(3, 3);
        S(0,0) = 1.0;
        S(1,1) = 1.0;
        S(2,2) = 1.0;

        double magnitude = BoussinesqHypothesis::strainRateMagnitude(S);

        // |S| = sqrt(2 * sum(S_ij^2))
        double expected = std::sqrt(2.0 * 3.0);
        ASSERT_NEAR(magnitude, expected, TOLERANCE);
        return true;
    });

    run_test("Boussinesq: Reynolds stress from eddy viscosity", []() {
        double mu_t = 0.1;  // Pa*s
        maths::linear_algebra::Matrix S(3, 3);
        S(0,0) = 2.0;

        double rho = 1.0;
        double k = 1.0;

        auto tau_R = BoussinesqHypothesis::reynoldsStress(mu_t, S, rho, k);

        // tau_R = mu_t * S - (2/3)*rho*k*I
        double expected_diag = mu_t * 2.0 - (2.0/3.0) * rho * k;
        ASSERT_NEAR(tau_R(0,0), expected_diag, TOLERANCE);
        return true;
    });

    // ========================================
    // Mixing Length Model Tests
    // ========================================

    run_test("Mixing length: eddy viscosity from mixing length", []() {
        double rho = 1.225;
        double l = 0.01;  // 1 cm
        double dU_dy = 100.0;  // 1/s

        double mu_t = MixingLengthModel::eddyViscosity(rho, l, dU_dy);

        // mu_t = rho * l^2 * |dU/dy|
        double expected = rho * l * l * std::abs(dU_dy);
        ASSERT_NEAR(mu_t, expected, TOLERANCE);
        return true;
    });

    run_test("Mixing length: wall mixing length at wall", []() {
        double y = 0.0;
        double delta = 0.1;

        double l = MixingLengthModel::mixingLengthWall(y, delta);

        // At wall, l = 0
        ASSERT_NEAR(l, 0.0, TOLERANCE);
        return true;
    });

    run_test("Mixing length: wall mixing length at center", []() {
        double y = 0.05;
        double delta = 0.1;
        double kappa = 0.41;

        double l = MixingLengthModel::mixingLengthWall(y, delta, kappa);

        // l = kappa * y * (1 - y/delta) = 0.41 * 0.05 * 0.5
        double expected = kappa * y * (1.0 - y/delta);
        ASSERT_NEAR(l, expected, TOLERANCE);
        return true;
    });

    run_test("Mixing length: free shear layer", []() {
        double delta = 0.05;
        double C = 0.07;

        double l = MixingLengthModel::mixingLengthShear(delta, C);

        // l = C * delta
        ASSERT_NEAR(l, C * delta, TOLERANCE);
        return true;
    });

    run_test("Mixing length: eddy viscosity scales with velocity gradient", []() {
        double rho = 1.0;
        double l = 0.01;

        double mu_t1 = MixingLengthModel::eddyViscosity(rho, l, 10.0);
        double mu_t2 = MixingLengthModel::eddyViscosity(rho, l, 20.0);

        // mu_t should double when gradient doubles
        ASSERT_NEAR(mu_t2 / mu_t1, 2.0, TOLERANCE);
        return true;
    });

    // ========================================
    // k-epsilon Model Tests
    // ========================================

    run_test("k-epsilon: model constants default values", []() {
        KEpsilonModel::Constants c;

        // Check standard k-epsilon constants
        ASSERT_NEAR(c.C_mu, 0.09, TOLERANCE);
        ASSERT_NEAR(c.C_eps1, 1.44, TOLERANCE);
        ASSERT_NEAR(c.C_eps2, 1.92, TOLERANCE);
        ASSERT_NEAR(c.sigma_k, 1.0, TOLERANCE);
        ASSERT_NEAR(c.sigma_eps, 1.3, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: eddy viscosity from k and epsilon", []() {
        double rho = 1.0;
        double k = 1.0;  // m^2/s^2
        double epsilon = 0.1;  // m^2/s^3

        double mu_t = KEpsilonModel::eddyViscosity(rho, k, epsilon);

        // mu_t = rho * C_mu * k^2 / epsilon
        double expected = rho * 0.09 * k * k / epsilon;
        ASSERT_NEAR(mu_t, expected, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: turbulence production", []() {
        double mu_t = 0.1;
        double S = 10.0;  // 1/s

        double P_k = KEpsilonModel::turbulenceProduction(mu_t, S);

        // P_k = mu_t * S^2
        ASSERT_NEAR(P_k, mu_t * S * S, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: turbulent length scale", []() {
        double k = 1.0;
        double epsilon = 0.1;

        double l_t = KEpsilonModel::turbulentLengthScale(k, epsilon);

        // l_t = k^(3/2) / epsilon
        double expected = std::pow(k, 1.5) / epsilon;
        ASSERT_NEAR(l_t, expected, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: turbulent time scale", []() {
        double k = 1.0;
        double epsilon = 0.5;

        double tau_t = KEpsilonModel::turbulentTimeScale(k, epsilon);

        // tau_t = k / epsilon
        ASSERT_NEAR(tau_t, k / epsilon, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: Kolmogorov length scale", []() {
        double nu = 1.5e-5;  // m^2/s (air)
        double epsilon = 1.0;  // m^2/s^3

        double eta = KEpsilonModel::kolmogorovLengthScale(nu, epsilon);

        // eta = (nu^3 / epsilon)^(1/4)
        double expected = std::pow(nu * nu * nu / epsilon, 0.25);
        ASSERT_NEAR(eta, expected, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: Kolmogorov time scale", []() {
        double nu = 1.5e-5;
        double epsilon = 1.0;

        double tau_eta = KEpsilonModel::kolmogorovTimeScale(nu, epsilon);

        // tau_eta = sqrt(nu / epsilon)
        double expected = std::sqrt(nu / epsilon);
        ASSERT_NEAR(tau_eta, expected, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: turbulent Reynolds number", []() {
        double k = 1.0;
        double nu = 1e-5;
        double epsilon = 0.1;

        double Re_t = KEpsilonModel::turbulentReynolds(k, nu, epsilon);

        // Re_t = k^2 / (nu * epsilon)
        double expected = k * k / (nu * epsilon);
        ASSERT_NEAR(Re_t, expected, TOLERANCE);
        return true;
    });

    run_test("k-epsilon: k equation step preserves positivity", []() {
        double k = 1.0;
        double production = 0.5;
        double epsilon = 0.3;
        double diffusion = 0.1;
        double dt = 0.01;

        double k_new = KEpsilonModel::stepKEquation(k, production, epsilon, diffusion, dt);

        // k should remain positive
        ASSERT_TRUE(k_new > 0.0);
        return true;
    });

    run_test("k-epsilon: epsilon equation step preserves positivity", []() {
        double eps = 0.5;
        double k = 1.0;
        double production = 0.3;
        double diffusion = 0.05;
        double dt = 0.01;

        double eps_new = KEpsilonModel::stepEpsilonEquation(eps, k, production, diffusion, dt);

        // epsilon should remain positive
        ASSERT_TRUE(eps_new > 0.0);
        return true;
    });

    run_test("k-epsilon: dissipation rate calculation", []() {
        double mu_t = 0.09;
        double k = 1.0;
        double rho = 1.0;

        double epsilon = KEpsilonModel::dissipationRate(mu_t, k, rho);

        // From mu_t = rho * C_mu * k^2 / epsilon
        double expected = rho * 0.09 * k * k / mu_t;
        ASSERT_NEAR(epsilon, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Wall Functions Tests
    // ========================================

    run_test("Wall functions: y plus calculation", []() {
        double y = 0.001;  // 1 mm
        double u_tau = 0.1;  // m/s
        double nu = 1.5e-5;  // m^2/s

        double y_plus = WallFunctions::yPlus(y, u_tau, nu);

        // y+ = y * u_tau / nu
        double expected = y * u_tau / nu;
        ASSERT_NEAR(y_plus, expected, TOLERANCE);
        return true;
    });

    run_test("Wall functions: viscous sublayer detection", []() {
        double y_plus = 3.0;
        bool is_viscous = WallFunctions::isViscousSublayer(y_plus);

        ASSERT_TRUE(is_viscous);  // y+ < 5
        return true;
    });

    run_test("Wall functions: log layer detection", []() {
        double y_plus = 100.0;
        bool is_log = WallFunctions::isLogLayer(y_plus);

        ASSERT_TRUE(is_log);  // 30 < y+ < 500
        return true;
    });

    run_test("Wall functions: buffer layer not in viscous or log", []() {
        double y_plus = 15.0;  // Buffer layer

        bool not_viscous = !WallFunctions::isViscousSublayer(y_plus);
        bool not_log = !WallFunctions::isLogLayer(y_plus);

        ASSERT_TRUE(not_viscous && not_log);
        return true;
    });

    run_test("Wall functions: wall shear stress from k-epsilon", []() {
        double rho = 1.225;
        double k = 0.5;
        double u_P = 10.0;

        double tau_w = WallFunctions::wallShearStress(rho, k, u_P);

        // tau_w = rho * C_mu^(1/4) * sqrt(k) * u_P
        double C_mu_quarter = std::pow(0.09, 0.25);
        double expected = rho * C_mu_quarter * std::sqrt(k) * u_P;
        ASSERT_NEAR(tau_w, expected, TOLERANCE);
        return true;
    });

    run_test("Wall functions: y plus transitions", []() {
        double nu = 1.5e-5;
        double u_tau = 0.1;

        double y1 = 5.0 * nu / u_tau;  // y+ = 5
        double y2 = 30.0 * nu / u_tau;  // y+ = 30

        double yp1 = WallFunctions::yPlus(y1, u_tau, nu);
        double yp2 = WallFunctions::yPlus(y2, u_tau, nu);

        ASSERT_TRUE(WallFunctions::isViscousSublayer(yp1));
        ASSERT_TRUE(WallFunctions::isLogLayer(yp2));
        return true;
    });

    // ========================================
    // Physical Consistency Tests
    // ========================================

    run_test("Physical: turbulent kinetic energy is positive", []() {
        std::vector<maths::linear_algebra::Vector> flucts;
        maths::linear_algebra::Vector u(3);
        u[0] = 1.0; u[1] = 0.5; u[2] = 0.3;
        flucts.push_back(u);

        double k = ReynoldsStressTensor::turbulentKineticEnergy(flucts);
        ASSERT_TRUE(k > 0.0);
        return true;
    });

    run_test("Physical: eddy viscosity is positive", []() {
        double rho = 1.0;
        double k = 1.0;
        double epsilon = 0.1;

        double mu_t = KEpsilonModel::eddyViscosity(rho, k, epsilon);
        ASSERT_TRUE(mu_t > 0.0);
        return true;
    });

    run_test("Physical: Kolmogorov scale smaller than integral scale", []() {
        double nu = 1.5e-5;
        double k = 1.0;
        double epsilon = 0.1;

        double eta = KEpsilonModel::kolmogorovLengthScale(nu, epsilon);
        double l_t = KEpsilonModel::turbulentLengthScale(k, epsilon);

        ASSERT_TRUE(eta < l_t);  // Kolmogorov << integral
        return true;
    });

    run_test("Physical: mixing length increases from wall", []() {
        double delta = 0.1;

        double l1 = MixingLengthModel::mixingLengthWall(0.01, delta);
        double l2 = MixingLengthModel::mixingLengthWall(0.02, delta);

        ASSERT_TRUE(l2 > l1);
        return true;
    });

    run_test("Physical: mixing length decreases near edge", []() {
        double delta = 0.1;

        double l_mid = MixingLengthModel::mixingLengthWall(0.05, delta);  // Center
        double l_edge = MixingLengthModel::mixingLengthWall(0.09, delta);  // Near edge

        ASSERT_TRUE(l_edge < l_mid);  // l(1-y/delta) decreases
        return true;
    });

    run_test("Physical: turbulence intensity bounds", []() {
        double k = 0.5;
        double U = 10.0;

        double I = ReynoldsStressTensor::turbulenceIntensity(k, U);

        // Physical turbulence intensity should be < 50% typically
        ASSERT_TRUE(I > 0.0 && I < 0.5);
        return true;
    });

    run_test("Physical: turbulent Reynolds number scaling", []() {
        double nu = 1e-5;
        double k = 1.0;

        double Re_t1 = KEpsilonModel::turbulentReynolds(k, nu, 0.1);
        double Re_t2 = KEpsilonModel::turbulentReynolds(k, nu, 0.05);

        // Re_t ~ 1/epsilon, so smaller epsilon → larger Re_t
        ASSERT_TRUE(Re_t2 > Re_t1);
        return true;
    });

    run_test("Physical: k-epsilon constants are positive", []() {
        KEpsilonModel::Constants c;

        ASSERT_TRUE(c.C_mu > 0.0);
        ASSERT_TRUE(c.C_eps1 > 0.0);
        ASSERT_TRUE(c.C_eps2 > 0.0);
        ASSERT_TRUE(c.sigma_k > 0.0);
        ASSERT_TRUE(c.sigma_eps > 0.0);
        return true;
    });

    run_test("Limit: high epsilon reduces turbulent scales", []() {
        double k = 1.0;
        double epsilon_low = 0.1;
        double epsilon_high = 10.0;

        double l_low = KEpsilonModel::turbulentLengthScale(k, epsilon_low);
        double l_high = KEpsilonModel::turbulentLengthScale(k, epsilon_high);

        ASSERT_TRUE(l_high < l_low);
        return true;
    });

    run_test("Limit: Kolmogorov scale depends on viscosity", []() {
        double epsilon = 1.0;

        double eta1 = KEpsilonModel::kolmogorovLengthScale(1e-5, epsilon);
        double eta2 = KEpsilonModel::kolmogorovLengthScale(1e-4, epsilon);

        // Higher viscosity → larger Kolmogorov scale
        ASSERT_TRUE(eta2 > eta1);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Fluid Turbulence" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All fluid dynamics turbulence tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Reynolds decomposition and averaging" << std::endl;
        std::cout << "  - Reynolds stress tensor and turbulent kinetic energy" << std::endl;
        std::cout << "  - Turbulence intensity calculations" << std::endl;
        std::cout << "  - Boussinesq eddy viscosity hypothesis" << std::endl;
        std::cout << "  - Mixing length models" << std::endl;
        std::cout << "  - k-epsilon turbulence model" << std::endl;
        std::cout << "  - Kolmogorov scales" << std::endl;
        std::cout << "  - Wall functions for near-wall modeling" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
