/**
 * Phase 4 Validation: Gauge Theory Symmetries
 *
 * Tests the gauge_theory_symmetries.hpp module functions.
 *
 * Coverage:
 * - Parity (P) transformation and violation
 * - Charge conjugation (C)
 * - Time reversal (T)
 * - CPT theorem
 * - CP violation discovery
 * - Discrete symmetry conservation
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/gauge_theory_symmetries.hpp"

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

using namespace physics::advanced::gauge_theory;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Gauge Theory Symmetries Validation ===" << std::endl;
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

    // ========================================
    // Parity Transformation Tests
    // ========================================

    run_test("Parity: position transformation string", []() {
        std::string transform = ParityTransformation::positionTransformation();
        ASSERT_TRUE(transform.find("r → -r") != std::string::npos);
        return true;
    });

    run_test("Parity: vector transformation", []() {
        std::string transform = ParityTransformation::vectorTransformation();
        ASSERT_TRUE(transform.find("P(v) = -v") != std::string::npos);
        ASSERT_TRUE(transform.find("P(a) = +a") != std::string::npos);
        return true;
    });

    run_test("Parity: fermion parity value", []() {
        int parity = ParityTransformation::fermionParity();
        ASSERT_NEAR(parity, +1, TOLERANCE);
        return true;
    });

    run_test("Parity: antifermion parity value", []() {
        int parity = ParityTransformation::antifermionParity();
        ASSERT_NEAR(parity, -1, TOLERANCE);
        return true;
    });

    run_test("Parity: photon has odd parity", []() {
        int parity = ParityTransformation::photonParity();
        ASSERT_NEAR(parity, -1, TOLERANCE);
        return true;
    });

    run_test("Parity: conserved in strong interaction", []() {
        bool conserved = ParityTransformation::conservedInStrong();
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("Parity: conserved in electromagnetic", []() {
        bool conserved = ParityTransformation::conservedInEM();
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("Parity: violated in weak interaction", []() {
        bool conserved = ParityTransformation::conservedInWeak();
        ASSERT_TRUE(!conserved);  // Violated!
        return true;
    });

    run_test("Parity: Wu experiment description", []() {
        std::string wu = ParityTransformation::wuExperiment();
        ASSERT_TRUE(wu.find("1956") != std::string::npos);
        ASSERT_TRUE(wu.find("Co") != std::string::npos);
        return true;
    });

    run_test("Parity: Lee-Yang prediction", []() {
        std::string ly = ParityTransformation::leeYangPrediction();
        ASSERT_TRUE(ly.find("1957") != std::string::npos);
        ASSERT_TRUE(ly.find("Nobel") != std::string::npos);
        return true;
    });

    // ========================================
    // Charge Conjugation Tests
    // ========================================

    run_test("Charge conjugation: transformation description", []() {
        std::string transform = ChargeConjugation::transformation();
        ASSERT_TRUE(transform.find("particle ↔ antiparticle") != std::string::npos);
        return true;
    });

    run_test("Charge conjugation: photon C eigenvalue", []() {
        int C_gamma = ChargeConjugation::photonEigenvalue();
        ASSERT_NEAR(C_gamma, -1, TOLERANCE);
        return true;
    });

    run_test("Charge conjugation: neutral pion C eigenvalue", []() {
        int C_pi0 = ChargeConjugation::neutralPionEigenvalue();
        ASSERT_NEAR(C_pi0, +1, TOLERANCE);
        return true;
    });

    run_test("Charge conjugation: pi0 C parity", []() {
        int C = ChargeConjugation::cParity("pi0");
        ASSERT_NEAR(C, +1, TOLERANCE);
        return true;
    });

    run_test("Charge conjugation: eta C parity", []() {
        int C = ChargeConjugation::cParity("eta");
        ASSERT_NEAR(C, +1, TOLERANCE);
        return true;
    });

    run_test("Charge conjugation: omega C parity", []() {
        int C = ChargeConjugation::cParity("omega");
        ASSERT_NEAR(C, -1, TOLERANCE);
        return true;
    });

    run_test("Charge conjugation: violated in weak", []() {
        bool conserved = ChargeConjugation::conservedInWeak();
        ASSERT_TRUE(!conserved);
        return true;
    });

    run_test("Charge conjugation: conserved in strong", []() {
        bool conserved = ChargeConjugation::conservedInStrong();
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("Charge conjugation: conserved in EM", []() {
        bool conserved = ChargeConjugation::conservedInEM();
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("Charge conjugation: weak violation description", []() {
        std::string violation = ChargeConjugation::weakCViolation();
        ASSERT_TRUE(violation.find("left-handed") != std::string::npos);
        return true;
    });

    // ========================================
    // Time Reversal Tests
    // ========================================

    run_test("Time reversal: transformation description", []() {
        std::string transform = TimeReversal::transformation();
        ASSERT_TRUE(transform.find("t → -t") != std::string::npos);
        ASSERT_TRUE(transform.find("p → -p") != std::string::npos);
        return true;
    });

    run_test("Time reversal: antiunitarity", []() {
        std::string anti = TimeReversal::antiunitarity();
        ASSERT_TRUE(anti.find("antiunitary") != std::string::npos);
        return true;
    });

    run_test("Time reversal: electric and magnetic fields", []() {
        std::string fields = TimeReversal::electricMagneticFields();
        ASSERT_TRUE(fields.find("E → E") != std::string::npos);
        ASSERT_TRUE(fields.find("B → -B") != std::string::npos);
        return true;
    });

    run_test("Time reversal: violated in weak", []() {
        bool conserved = TimeReversal::conservedInWeak();
        ASSERT_TRUE(!conserved);
        return true;
    });

    run_test("Time reversal: conserved in strong", []() {
        bool conserved = TimeReversal::conservedInStrong();
        ASSERT_TRUE(conserved);
        return true;
    });

    run_test("Time reversal: T violation description", []() {
        std::string violation = TimeReversal::tViolation();
        ASSERT_TRUE(violation.find("kaon") != std::string::npos);
        return true;
    });

    // ========================================
    // CPT Theorem Tests
    // ========================================

    run_test("CPT: theorem statement", []() {
        std::string statement = CPTTheorem::statement();
        ASSERT_TRUE(statement.find("CPT") != std::string::npos);
        ASSERT_TRUE(statement.find("Lorentz-invariant") != std::string::npos);
        return true;
    });

    run_test("CPT: transformation description", []() {
        std::string transform = CPTTheorem::transformation();
        ASSERT_TRUE(transform.find("antiparticle") != std::string::npos);
        return true;
    });

    run_test("CPT: consequences description", []() {
        std::string conseq = CPTTheorem::consequences();
        ASSERT_TRUE(conseq.find("masses") != std::string::npos);
        return true;
    });

    run_test("CPT: electron-positron mass equality", []() {
        double test_limit = CPTTheorem::massEqualityTest();
        ASSERT_NEAR(test_limit, 1e-8, TOLERANCE);
        return true;
    });

    run_test("CPT: proton-antiproton charge equality", []() {
        double test_limit = CPTTheorem::chargeEqualityTest();
        ASSERT_NEAR(test_limit, 1e-21, 1e-22);
        return true;
    });

    run_test("CPT: is not violated", []() {
        bool violated = CPTTheorem::isViolated();
        ASSERT_TRUE(!violated);
        return true;
    });

    run_test("CPT: experimental status", []() {
        std::string status = CPTTheorem::experimentalStatus();
        ASSERT_TRUE(status.find("precision") != std::string::npos);
        return true;
    });

    run_test("CPT: CP violation implies T violation", []() {
        std::string implication = CPTTheorem::cpViolationImpliesTViolation();
        ASSERT_TRUE(implication.find("T violation") != std::string::npos);
        return true;
    });

    // ========================================
    // CP Violation Tests
    // ========================================

    run_test("CP violation: transformation", []() {
        std::string transform = CPViolation::transformation();
        ASSERT_TRUE(transform.find("CP") != std::string::npos);
        return true;
    });

    run_test("CP violation: discovery description", []() {
        std::string discovery = CPViolation::discovery();
        ASSERT_TRUE(discovery.find("Cronin") != std::string::npos);
        ASSERT_TRUE(discovery.find("Fitch") != std::string::npos);
        ASSERT_TRUE(discovery.find("1964") != std::string::npos);
        return true;
    });

    run_test("CP violation: kaon epsilon parameter", []() {
        double epsilon = CPViolation::kaonCPViolation();
        ASSERT_NEAR(epsilon, 2.2e-3, 0.5e-3);
        return true;
    });

    run_test("CP violation: B meson sin 2beta", []() {
        double sin2beta = CPViolation::bMesonCPViolation();
        ASSERT_NEAR(sin2beta, 0.7, 0.1);
        return true;
    });

    run_test("CP violation: sources in SM", []() {
        std::string sources = CPViolation::sourcesInSM();
        ASSERT_TRUE(sources.find("CKM") != std::string::npos);
        return true;
    });

    run_test("CP violation: cosmological importance", []() {
        std::string importance = CPViolation::cosmologicalImportance();
        ASSERT_TRUE(importance.find("baryogenesis") != std::string::npos);
        return true;
    });

    run_test("CP violation: electron EDM limit", []() {
        double edm = CPViolation::electronEDMLimit();
        ASSERT_TRUE(edm < 1e-28);
        return true;
    });

    run_test("CP violation: EDM significance", []() {
        std::string sig = CPViolation::edmSignificance();
        ASSERT_TRUE(sig.find("EDM") != std::string::npos);
        return true;
    });

    // ========================================
    // T Violation Tests
    // ========================================

    run_test("T violation: direct observation", []() {
        std::string obs = TViolation::directObservation();
        ASSERT_TRUE(obs.find("CPLEAR") != std::string::npos);
        return true;
    });

    run_test("T violation: magnitude", []() {
        double mag = TViolation::magnitude();
        ASSERT_NEAR(mag, 2.2e-3, 0.5e-3);
        return true;
    });

    run_test("T violation: thermodynamic arrow", []() {
        std::string arrow = TViolation::thermodynamicArrow();
        ASSERT_TRUE(arrow.find("arrow") != std::string::npos);
        return true;
    });

    // ========================================
    // Combined Symmetries Tests
    // ========================================

    run_test("Combined: conservation table", []() {
        std::string table = CombinedSymmetries::conservationTable();
        ASSERT_TRUE(table.find("Strong") != std::string::npos);
        ASSERT_TRUE(table.find("Weak") != std::string::npos);
        return true;
    });

    run_test("Combined: historical timeline", []() {
        std::string timeline = CombinedSymmetries::timeline();
        ASSERT_TRUE(timeline.find("1957") != std::string::npos);
        ASSERT_TRUE(timeline.find("1964") != std::string::npos);
        return true;
    });

    run_test("Combined: physical significance", []() {
        std::string sig = CombinedSymmetries::significance();
        ASSERT_TRUE(sig.find("CPT") != std::string::npos);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Gauge Symmetries" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "All gauge theory symmetry tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Parity (P) transformation and violation" << std::endl;
        std::cout << "  - Charge conjugation (C) symmetry" << std::endl;
        std::cout << "  - Time reversal (T) symmetry" << std::endl;
        std::cout << "  - CPT theorem and its consequences" << std::endl;
        std::cout << "  - CP violation discovery and magnitude" << std::endl;
        std::cout << "  - T violation observations" << std::endl;
        std::cout << "  - Discrete symmetry conservation laws" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
