/**
 * Phase 4 Validation: Gauge Theory CP Violation in Kaons
 */

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "../include/physics/gauge_theory_cp_violation_kaons.hpp"

const double TOLERANCE = 1e-6;
const double NUMERICAL_TOLERANCE = 1e-3;
const double LOOSE_TOLERANCE = 0.01;

#define ASSERT_NEAR(actual, expected, tolerance) \
    do { if (std::abs((actual) - (expected)) > (tolerance)) { std::cerr << "FAIL: " << __LINE__ << std::endl; return false; } } while(0)

#define ASSERT_TRUE(condition) \
    do { if (!(condition)) { std::cerr << "FAIL: " << __LINE__ << std::endl; return false; } } while(0)

using namespace physics::advanced::gauge_theory;

int main() {
    int tests_passed = 0, tests_failed = 0;
    std::cout << "=== Phase 4: CP Violation Kaons Validation ===" << std::endl << std::endl;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) { std::cout << "PASS" << std::endl; tests_passed++; return true; }
        else { tests_failed++; return false; }
    };

    run_test("Neutral kaon: flavor eigenstates", []() {
        std::string flavor = NeutralKaonSystem::flavorEigenstates();
        ASSERT_TRUE(flavor.find("K⁰") != std::string::npos);
        return true;
    });

    run_test("Neutral kaon: K_S mass", []() {
        double m = NeutralKaonSystem::kShortMass();
        ASSERT_NEAR(m, 0.4976, 0.001);
        return true;
    });

    run_test("Neutral kaon: K_L mass", []() {
        double m = NeutralKaonSystem::kLongMass();
        ASSERT_NEAR(m, 0.4976, 0.001);
        return true;
    });

    run_test("Neutral kaon: K_S lifetime", []() {
        double tau = NeutralKaonSystem::kShortLifetime();
        ASSERT_TRUE(tau > 0.0 && tau < 1e-9);
        return true;
    });

    run_test("Neutral kaon: K_L lifetime much longer", []() {
        double tau_s = NeutralKaonSystem::kShortLifetime();
        double tau_l = NeutralKaonSystem::kLongLifetime();
        ASSERT_TRUE(tau_l > 100.0 * tau_s);
        return true;
    });

    run_test("Neutral kaon: mass difference", []() {
        double dm = NeutralKaonSystem::massDifference();
        ASSERT_TRUE(dm > 0.0 && dm < 1e-14);
        return true;
    });

    run_test("CP violation: Cronin-Fitch experiment", []() {
        std::string cf = CPViolationDiscovery::croninFitchExperiment();
        ASSERT_TRUE(cf.find("1964") != std::string::npos);
        ASSERT_TRUE(cf.find("Nobel") != std::string::npos);
        return true;
    });

    run_test("CP violation: epsilon parameter", []() {
        double eps = CPViolationDiscovery::epsilonParameter();
        ASSERT_NEAR(eps, 2.2e-3, 0.5e-3);
        return true;
    });

    run_test("Direct vs indirect: epsilon prime", []() {
        double eps_prime = DirectVsIndirectCPViolation::epsilonPrimeOverEpsilon();
        ASSERT_TRUE(eps_prime > 0.0 && eps_prime < 0.01);
        return true;
    });

    run_test("CKM: matrix elements size", []() {
        auto matrix = CKMMatrix::matrixElements();
        ASSERT_TRUE(matrix.size() == 3);
        ASSERT_TRUE(matrix[0].size() == 3);
        return true;
    });

    run_test("CKM: V_ud magnitude", []() {
        auto matrix = CKMMatrix::matrixElements();
        ASSERT_NEAR(matrix[0][0], 0.974, 0.01);
        return true;
    });

    run_test("CKM: V_us magnitude", []() {
        auto matrix = CKMMatrix::matrixElements();
        ASSERT_NEAR(matrix[0][1], 0.224, 0.01);
        return true;
    });

    run_test("CKM: V_ub magnitude small", []() {
        auto matrix = CKMMatrix::matrixElements();
        ASSERT_TRUE(matrix[0][2] < 0.01);
        return true;
    });

    run_test("Unitarity triangle: rhobar", []() {
        double rho = UnitarityTriangle::rhoBar();
        ASSERT_TRUE(rho > 0.0 && rho < 0.5);
        return true;
    });

    run_test("Unitarity triangle: etabar", []() {
        double eta = UnitarityTriangle::etaBar();
        ASSERT_TRUE(eta > 0.0 && eta < 0.5);
        return true;
    });

    run_test("Unitarity triangle: alpha angle", []() {
        double alpha = UnitarityTriangle::alpha();
        ASSERT_TRUE(alpha > 50.0 && alpha < 130.0);
        return true;
    });

    run_test("Unitarity triangle: beta angle", []() {
        double beta = UnitarityTriangle::beta();
        ASSERT_TRUE(beta > 10.0 && beta < 40.0);
        return true;
    });

    run_test("Unitarity triangle: gamma angle", []() {
        double gamma = UnitarityTriangle::gamma();
        ASSERT_TRUE(gamma > 40.0 && gamma < 100.0);
        return true;
    });

    run_test("Unitarity triangle: angles sum to 180", []() {
        double alpha = UnitarityTriangle::alpha();
        double beta = UnitarityTriangle::beta();
        double gamma = UnitarityTriangle::gamma();
        ASSERT_NEAR(alpha + beta + gamma, 180.0, 5.0);
        return true;
    });

    run_test("Unitarity triangle: Jarlskog invariant", []() {
        double J = UnitarityTriangle::jarlskogInvariant();
        ASSERT_TRUE(J > 1e-6 && J < 1e-4);
        return true;
    });

    run_test("B meson: sin 2beta", []() {
        double sin2b = BMesonCPViolation::sin2Beta();
        ASSERT_NEAR(sin2b, 0.7, 0.1);
        return true;
    });

    run_test("CP in SM: sources description", []() {
        std::string sources = CPViolationInSM::sourcesInSM();
        ASSERT_TRUE(sources.find("CKM") != std::string::npos);
        return true;
    });

    run_test("CP in SM: baryogenesis insufficiency", []() {
        std::string bary = CPViolationInSM::baryogenesisInsufficiency();
        ASSERT_TRUE(bary.find("insufficient") != std::string::npos);
        return true;
    });

    std::cout << std::endl << "======================================" << std::endl;
    std::cout << "Phase 4 Results: CP Violation Kaons" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    if (tests_failed == 0) {
        std::cout << std::endl << "All CP violation kaons tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << std::endl << "Some tests FAILED." << std::endl;
        return 1;
    }
}
