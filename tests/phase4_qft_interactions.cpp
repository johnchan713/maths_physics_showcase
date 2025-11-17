/**
 * Phase 4 Validation: QFT Interactions
 */

#include <iostream>
#include <cmath>
#include <string>
#include <complex>
#include "../include/physics/qft_interactions.hpp"

const double TOLERANCE = 1e-6;
const double NUMERICAL_TOLERANCE = 1e-3;
const double LOOSE_TOLERANCE = 0.01;

#define ASSERT_NEAR(actual, expected, tolerance) \
    do { if (std::abs((actual) - (expected)) > (tolerance)) { std::cerr << "FAIL: " << __LINE__ << std::endl; return false; } } while(0)

#define ASSERT_TRUE(condition) \
    do { if (!(condition)) { std::cerr << "FAIL: " << __LINE__ << std::endl; return false; } } while(0)

using namespace physics::advanced::qft;

int main() {
    int tests_passed = 0, tests_failed = 0;
    std::cout << "=== Phase 4: QFT Interactions Validation ===" << std::endl << std::endl;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) { std::cout << "PASS" << std::endl; tests_passed++; return true; }
        else { tests_failed++; return false; }
    };

    run_test("Boson exchange: force range for massless boson", []() {
        double range = BosonExchange::forceRange(0.0);
        ASSERT_TRUE(std::isinf(range));
        return true;
    });

    run_test("Boson exchange: force range for W boson", []() {
        double mW = 80.4;  // GeV
        double range = BosonExchange::forceRange(mW);
        ASSERT_TRUE(range > 0.0 && range < 1e-17);
        return true;
    });

    run_test("Coupling: fine structure constant", []() {
        double alpha = CouplingConstants::fineStructure();
        ASSERT_NEAR(alpha, 1.0/137.0, 0.0001);
        return true;
    });

    run_test("Coupling: strong coupling at MZ", []() {
        double alpha_s = CouplingConstants::strong(91.2);
        ASSERT_NEAR(alpha_s, 0.1179, 0.01);
        return true;
    });

    run_test("Coupling: weak coupling", []() {
        double alpha_w = CouplingConstants::weak();
        ASSERT_NEAR(alpha_w, 1.0/30.0, 0.01);
        return true;
    });

    run_test("Coupling: gravitational is tiny", []() {
        double alpha_g = CouplingConstants::gravitational();
        ASSERT_TRUE(alpha_g < 1e-38);
        return true;
    });

    run_test("Running: QCD asymptotic freedom", []() {
        bool free = RunningCouplings::isAsymptoticallyFree(1000.0);
        ASSERT_TRUE(free);
        return true;
    });

    run_test("Running: QCD Lambda", []() {
        double lambda = RunningCouplings::lambdaQCD();
        ASSERT_NEAR(lambda, 0.2, 0.1);
        return true;
    });

    run_test("Running: confinement at low scale", []() {
        bool conf = RunningCouplings::isConfined(0.1);
        ASSERT_TRUE(conf);
        return true;
    });

    run_test("Fermion-boson: QED vertex charge dependence", []() {
        auto v1 = FermionBosonCoupling::qedVertex(1.0);
        auto v2 = FermionBosonCoupling::qedVertex(2.0);
        ASSERT_NEAR(std::abs(v2) / std::abs(v1), 2.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Fermion-boson: QCD vertex", []() {
        auto v = FermionBosonCoupling::qcdVertex(91.2);
        ASSERT_TRUE(std::abs(v) > 0.0);
        return true;
    });

    run_test("Fermion-boson: Yukawa top quark", []() {
        double yt = FermionBosonCoupling::topYukawa();
        ASSERT_NEAR(yt, 1.0, 0.1);
        return true;
    });

    run_test("Fermion-boson: Yukawa scales with mass", []() {
        double y1 = FermionBosonCoupling::yukawaCoupling(1.0);
        double y2 = FermionBosonCoupling::yukawaCoupling(2.0);
        ASSERT_NEAR(y2 / y1, 2.0, TOLERANCE);
        return true;
    });

    run_test("Interaction range: EM is infinite", []() {
        double range = InteractionRanges::electromagnetic();
        ASSERT_TRUE(std::isinf(range));
        return true;
    });

    run_test("Interaction range: weak is short", []() {
        double range = InteractionRanges::weak();
        ASSERT_TRUE(range < 1e-17);
        return true;
    });

    run_test("Interaction range: strong is nuclear scale", []() {
        double range = InteractionRanges::strong();
        ASSERT_NEAR(range, 1e-15, 0.5e-15);
        return true;
    });

    run_test("Physical: strong stronger than EM", []() {
        double alpha = CouplingConstants::fineStructure();
        double alpha_s = CouplingConstants::strong(91.2);
        ASSERT_TRUE(alpha_s > alpha);
        return true;
    });

    run_test("Physical: weak is suppressed", []() {
        double alpha = CouplingConstants::fineStructure();
        double alpha_w = CouplingConstants::weak();
        ASSERT_TRUE(alpha_w > alpha && alpha_w < 0.1);
        return true;
    });

    run_test("Physical: gravity is weakest", []() {
        double alpha_g = CouplingConstants::gravitational();
        double alpha = CouplingConstants::fineStructure();
        ASSERT_TRUE(alpha_g < alpha * 1e-30);
        return true;
    });

    std::cout << std::endl << "======================================" << std::endl;
    std::cout << "Phase 4 Results: QFT Interactions" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    if (tests_failed == 0) {
        std::cout << std::endl << "All QFT interactions tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << std::endl << "Some tests FAILED." << std::endl;
        return 1;
    }
}
