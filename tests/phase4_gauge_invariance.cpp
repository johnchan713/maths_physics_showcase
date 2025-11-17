/**
 * Phase 4 Validation: Gauge Theory Gauge Invariance
 *
 * Tests the gauge_theory_gauge_invariance.hpp module functions.
 *
 * Coverage:
 * - U(1) gauge theory (QED)
 * - SU(2) gauge theory (Weak isospin)
 * - Electroweak gauge invariance
 * - SU(3) gauge theory (QCD)
 * - Standard Model gauge group
 * - Gauge transformations and covariant derivatives
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/gauge_theory_gauge_invariance.hpp"

const double TOLERANCE = 1e-6;
const double NUMERICAL_TOLERANCE = 1e-3;
const double LOOSE_TOLERANCE = 0.01;

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

using namespace physics::advanced::gauge_theory;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Gauge Invariance Validation ===" << std::endl;
    std::cout << std::endl;

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

    // Gauge Transformation Tests
    run_test("Gauge: global vs local symmetry description", []() {
        std::string desc = GaugeTransformation::globalVsLocal();
        ASSERT_TRUE(desc.find("Global") != std::string::npos);
        ASSERT_TRUE(desc.find("Local") != std::string::npos);
        return true;
    });

    run_test("Gauge: motivation for gauge symmetry", []() {
        std::string motiv = GaugeTransformation::motivation();
        ASSERT_TRUE(motiv.find("boson") != std::string::npos);
        return true;
    });

    run_test("Gauge: redundancy description", []() {
        std::string redun = GaugeTransformation::gaugeRedundancy();
        ASSERT_TRUE(redun.find("A_μ") != std::string::npos);
        return true;
    });

    // U(1) Gauge Theory Tests
    run_test("U1: gauge transformation formula", []() {
        std::string trans = U1GaugeTheory::gaugeTransformation();
        ASSERT_TRUE(trans.find("ψ(x) → ψ'(x)") != std::string::npos);
        return true;
    });

    run_test("U1: covariant derivative", []() {
        std::string cov = U1GaugeTheory::covariantDerivative();
        ASSERT_TRUE(cov.find("D_μ") != std::string::npos);
        return true;
    });

    run_test("U1: minimal coupling", []() {
        std::string min = U1GaugeTheory::minimalCoupling();
        ASSERT_TRUE(min.find("interaction") != std::string::npos);
        return true;
    });

    run_test("U1: field strength tensor", []() {
        std::string field = U1GaugeTheory::fieldStrength();
        ASSERT_TRUE(field.find("F_μν") != std::string::npos);
        return true;
    });

    run_test("U1: QED Lagrangian", []() {
        std::string lag = U1GaugeTheory::lagrangian();
        ASSERT_TRUE(lag.find("invariant") != std::string::npos);
        return true;
    });

    run_test("U1: charge quantization", []() {
        std::string quant = U1GaugeTheory::chargeQuantization();
        ASSERT_TRUE(quant.find("quantized") != std::string::npos);
        return true;
    });

    // SU(2) Gauge Theory Tests
    run_test("SU2: gauge transformation", []() {
        std::string trans = SU2GaugeTheory::gaugeTransformation();
        ASSERT_TRUE(trans.find("SU(2)") != std::string::npos);
        return true;
    });

    run_test("SU2: Lie algebra commutation", []() {
        std::string lie = SU2GaugeTheory::lieAlgebra();
        ASSERT_TRUE(lie.find("τᵃ") != std::string::npos);
        return true;
    });

    run_test("SU2: covariant derivative", []() {
        std::string cov = SU2GaugeTheory::covariantDerivative();
        ASSERT_TRUE(cov.find("D_μ") != std::string::npos);
        return true;
    });

    run_test("SU2: field strength tensor", []() {
        std::string field = SU2GaugeTheory::fieldStrength();
        ASSERT_TRUE(field.find("self-interact") != std::string::npos);
        return true;
    });

    run_test("SU2: Yang-Mills Lagrangian", []() {
        std::string ym = SU2GaugeTheory::yangMillsLagrangian();
        ASSERT_TRUE(ym.find("vertices") != std::string::npos);
        return true;
    });

    run_test("SU2: number of gauge bosons", []() {
        int N = 2;
        int num = SU2GaugeTheory::numberOfGaugeBosons(N);
        ASSERT_NEAR(num, 3, TOLERANCE);  // SU(2): 2^2-1=3
        return true;
    });

    run_test("SU3: number of gauge bosons", []() {
        int N = 3;
        int num = SU2GaugeTheory::numberOfGaugeBosons(N);
        ASSERT_NEAR(num, 8, TOLERANCE);  // SU(3): 3^2-1=8
        return true;
    });

    // Electroweak Gauge Theory Tests
    run_test("EW: gauge group description", []() {
        std::string group = ElectroweakGaugeInvariance::gaugeGroup();
        ASSERT_TRUE(group.find("SU(2)_L × U(1)_Y") != std::string::npos);
        return true;
    });

    run_test("EW: gauge fields", []() {
        std::string fields = ElectroweakGaugeInvariance::gaugeFields();
        ASSERT_TRUE(fields.find("W_μ") != std::string::npos);
        ASSERT_TRUE(fields.find("B_μ") != std::string::npos);
        return true;
    });

    run_test("EW: weak hypercharge formula", []() {
        std::string hyper = ElectroweakGaugeInvariance::hypercharge();
        ASSERT_TRUE(hyper.find("Y = 2(Q - T₃)") != std::string::npos);
        return true;
    });

    run_test("EW: physical gauge bosons", []() {
        std::string phys = ElectroweakGaugeInvariance::physicalBosons();
        ASSERT_TRUE(phys.find("W±") != std::string::npos);
        ASSERT_TRUE(phys.find("Z⁰") != std::string::npos);
        ASSERT_TRUE(phys.find("γ") != std::string::npos);
        return true;
    });

    run_test("EW: Weinberg angle value", []() {
        double theta_W = ElectroweakGaugeInvariance::weinbergAngle();
        ASSERT_TRUE(theta_W > 0.0);
        ASSERT_TRUE(theta_W < M_PI / 2.0);
        return true;
    });

    run_test("EW: sin squared theta W", []() {
        double sin2 = ElectroweakGaugeInvariance::sin2ThetaW();
        ASSERT_NEAR(sin2, 0.231, 0.01);
        return true;
    });

    run_test("EW: electromagnetic coupling", []() {
        std::string em = ElectroweakGaugeInvariance::electromagneticCoupling();
        ASSERT_TRUE(em.find("e = g sin θ_W") != std::string::npos);
        return true;
    });

    run_test("EW: custodial symmetry", []() {
        std::string cust = ElectroweakGaugeInvariance::custodialSymmetry();
        ASSERT_TRUE(cust.find("ρ") != std::string::npos);
        return true;
    });

    // SU(3) QCD Tests
    run_test("QCD: gauge group description", []() {
        std::string group = SU3GaugeTheory::gaugeGroup();
        ASSERT_TRUE(group.find("SU(3)_C") != std::string::npos);
        ASSERT_TRUE(group.find("8 gluons") != std::string::npos);
        return true;
    });

    run_test("QCD: Gell-Mann matrices", []() {
        std::string gm = SU3GaugeTheory::gellMannMatrices();
        ASSERT_TRUE(gm.find("λ") != std::string::npos);
        return true;
    });

    run_test("QCD: covariant derivative", []() {
        std::string cov = SU3GaugeTheory::covariantDerivative();
        ASSERT_TRUE(cov.find("g_s") != std::string::npos);
        return true;
    });

    run_test("QCD: field strength", []() {
        std::string field = SU3GaugeTheory::fieldStrength();
        ASSERT_TRUE(field.find("gluon self-coupling") != std::string::npos);
        return true;
    });

    run_test("QCD: Lagrangian", []() {
        std::string lag = SU3GaugeTheory::lagrangian();
        ASSERT_TRUE(lag.find("confinement") != std::string::npos);
        return true;
    });

    run_test("QCD: confinement description", []() {
        std::string conf = SU3GaugeTheory::confinement();
        ASSERT_TRUE(conf.find("color-singlet") != std::string::npos);
        return true;
    });

    run_test("QCD: asymptotic freedom", []() {
        std::string asymp = SU3GaugeTheory::asymptoticFreedom();
        ASSERT_TRUE(asymp.find("Nobel 2004") != std::string::npos);
        return true;
    });

    // Standard Model Gauge Group Tests
    run_test("SM: complete gauge group", []() {
        std::string complete = StandardModelGaugeGroup::completeGaugeGroup();
        ASSERT_TRUE(complete.find("SU(3)_C × SU(2)_L × U(1)_Y") != std::string::npos);
        return true;
    });

    run_test("SM: coupling constants", []() {
        std::string coup = StandardModelGaugeGroup::couplingConstants();
        ASSERT_TRUE(coup.find("α_s") != std::string::npos);
        return true;
    });

    run_test("SM: quantum numbers", []() {
        std::string qn = StandardModelGaugeGroup::quantumNumbers();
        ASSERT_TRUE(qn.find("SU(3)") != std::string::npos);
        return true;
    });

    run_test("SM: anomaly cancellation", []() {
        std::string anom = StandardModelGaugeGroup::anomalyCancellation();
        ASSERT_TRUE(anom.find("cancel") != std::string::npos);
        return true;
    });

    run_test("SM: why three colors", []() {
        std::string why = StandardModelGaugeGroup::whyThreeColors();
        ASSERT_TRUE(why.find("N_c = 3") != std::string::npos);
        return true;
    });

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Gauge Invariance" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All gauge invariance tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - U(1) gauge theory (QED)" << std::endl;
        std::cout << "  - SU(2) gauge theory (Weak)" << std::endl;
        std::cout << "  - Electroweak unification SU(2)×U(1)" << std::endl;
        std::cout << "  - SU(3) gauge theory (QCD)" << std::endl;
        std::cout << "  - Standard Model gauge group" << std::endl;
        std::cout << "  - Gauge transformations and covariant derivatives" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
