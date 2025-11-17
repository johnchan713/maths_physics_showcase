/**
 * Phase 4 Validation: Gauge Theory Running Couplings
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/gauge_theory_running_couplings.hpp"

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
    std::cout << "=== Phase 4: Running Couplings Validation ===" << std::endl << std::endl;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) { std::cout << "PASS" << std::endl; tests_passed++; return true; }
        else { tests_failed++; return false; }
    };

    run_test("RG: beta function definition", []() {
        std::string beta = RenormalizationGroup::betaFunction();
        ASSERT_TRUE(beta.find("β(g)") != std::string::npos);
        return true;
    });

    run_test("RG: one-loop beta", []() {
        std::string one = RenormalizationGroup::oneLoopBeta();
        ASSERT_TRUE(one.find("b₀") != std::string::npos);
        return true;
    });

    run_test("RG: running coupling calculation", []() {
        double alpha_ref = 0.1;
        double Q2 = 10000.0;
        double mu2 = 100.0;
        double b0 = 7.0;
        double alpha = RenormalizationGroup::runningCoupling(alpha_ref, Q2, mu2, b0);
        ASSERT_TRUE(alpha > 0.0 && alpha < alpha_ref);
        return true;
    });

    run_test("QED: beta coefficient", []() {
        int nf = 3;
        double b0 = QEDRunning::betaCoefficient(nf);
        ASSERT_NEAR(b0, -4.0, TOLERANCE);
        return true;
    });

    run_test("QED: alpha at electron mass", []() {
        double alpha = QEDRunning::alphaAtElectronMass();
        ASSERT_NEAR(alpha, 1.0/137.0, 0.0001);
        return true;
    });

    run_test("QED: alpha at Z mass", []() {
        double alpha = QEDRunning::alphaAtZMass();
        ASSERT_NEAR(alpha, 1.0/128.0, 0.001);
        return true;
    });

    run_test("QED: alpha increases with energy", []() {
        double alpha_me = QEDRunning::alphaAtElectronMass();
        double alpha_mz = QEDRunning::alphaAtZMass();
        ASSERT_TRUE(alpha_mz > alpha_me);
        return true;
    });

    run_test("QCD: beta coefficient for 5 flavors", []() {
        double b0 = QCDRunning::betaCoefficient(5);
        ASSERT_NEAR(b0, 7.667, 0.01);
        return true;
    });

    run_test("QCD: alpha_s at Z mass", []() {
        double alpha_s = QCDRunning::alphaSAtZMass();
        ASSERT_NEAR(alpha_s, 0.1179, 0.01);
        return true;
    });

    run_test("QCD: alpha_s decreases with energy", []() {
        double alpha_s1 = QCDRunning::runningAlphaS(10.0);
        double alpha_s2 = QCDRunning::runningAlphaS(100.0);
        ASSERT_TRUE(alpha_s2 < alpha_s1);
        return true;
    });

    run_test("QCD: Lambda_QCD value", []() {
        double lambda = QCDRunning::lambdaQCD();
        ASSERT_NEAR(lambda, 0.213, 0.05);
        return true;
    });

    run_test("EW: sin2theta_W at Z mass", []() {
        double sin2 = ElectroweakRunning::sin2ThetaWAtZMass();
        ASSERT_NEAR(sin2, 0.231, 0.01);
        return true;
    });

    run_test("EW: rho parameter", []() {
        double rho = ElectroweakRunning::rhoParameter();
        ASSERT_NEAR(rho, 1.0, 0.01);
        return true;
    });

    run_test("GUT: GUT scale", []() {
        double mgut = GrandUnification::gutScale();
        ASSERT_TRUE(mgut > 1e15);
        return true;
    });

    run_test("GUT: GUT coupling", []() {
        double alpha_gut = GrandUnification::gutCoupling();
        ASSERT_NEAR(alpha_gut, 0.04, 0.01);
        return true;
    });

    run_test("Physical: QCD asymptotic freedom", []() {
        double alpha_high = QCDRunning::runningAlphaS(1000.0);
        ASSERT_TRUE(alpha_high < 0.1);
        return true;
    });

    std::cout << std::endl << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Running Couplings" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    if (tests_failed == 0) {
        std::cout << std::endl << "All running couplings tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << std::endl << "Some tests FAILED." << std::endl;
        return 1;
    }
}
