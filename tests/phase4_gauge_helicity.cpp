/**
 * Phase 4 Validation: Gauge Theory Helicity
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/gauge_theory_helicity.hpp"

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
    std::cout << "=== Phase 4: Gauge Helicity Validation ===" << std::endl << std::endl;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) { std::cout << "PASS" << std::endl; tests_passed++; return true; }
        else { tests_failed++; return false; }
    };

    run_test("Helicity: definition string", []() {
        std::string def = Helicity::definition();
        ASSERT_TRUE(def.find("S · p̂") != std::string::npos);
        return true;
    });

    run_test("Helicity: right-handed value", []() {
        double h = Helicity::rightHanded();
        ASSERT_NEAR(h, 0.5, TOLERANCE);
        return true;
    });

    run_test("Helicity: left-handed value", []() {
        double h = Helicity::leftHanded();
        ASSERT_NEAR(h, -0.5, TOLERANCE);
        return true;
    });

    run_test("Helicity: operator form", []() {
        std::string op = Helicity::operator_form();
        ASSERT_TRUE(op.find("Σ") != std::string::npos);
        return true;
    });

    run_test("Helicity: Lorentz invariance for massless", []() {
        bool inv = Helicity::isLorentzInvariant(0.0);
        ASSERT_TRUE(inv);
        return true;
    });

    run_test("Helicity: not Lorentz invariant for massive", []() {
        bool inv = Helicity::isLorentzInvariant(1.0);
        ASSERT_TRUE(!inv);
        return true;
    });

    run_test("Helicity: photon helicity states", []() {
        std::string photon = Helicity::photonHelicity();
        ASSERT_TRUE(photon.find("h = ±1") != std::string::npos);
        return true;
    });

    run_test("Helicity: gluon helicity", []() {
        std::string gluon = Helicity::gluonHelicity();
        ASSERT_TRUE(gluon.find("±1") != std::string::npos);
        return true;
    });

    run_test("Chirality: gamma5 operator", []() {
        std::string g5 = Chirality::gammaFiveOperator();
        ASSERT_TRUE(g5.find("γ⁵") != std::string::npos);
        return true;
    });

    run_test("Chirality: projection operators", []() {
        std::string proj = Chirality::projectionOperators();
        ASSERT_TRUE(proj.find("P_L") != std::string::npos);
        ASSERT_TRUE(proj.find("P_R") != std::string::npos);
        return true;
    });

    run_test("Chirality: Weyl spinors", []() {
        std::string weyl = Chirality::weylSpinors();
        ASSERT_TRUE(weyl.find("ψ_L") != std::string::npos);
        return true;
    });

    run_test("Chirality: is Lorentz invariant", []() {
        bool inv = Chirality::isLorentzInvariant();
        ASSERT_TRUE(inv);
        return true;
    });

    run_test("Chirality: massless limit equals helicity", []() {
        std::string limit = Chirality::masslessLimit();
        ASSERT_TRUE(limit.find("chirality = helicity") != std::string::npos);
        return true;
    });

    run_test("Chirality: massive case differs", []() {
        std::string massive = Chirality::massiveCase();
        ASSERT_TRUE(massive.find("≠") != std::string::npos);
        return true;
    });

    run_test("Helicity conservation: conserved for massless", []() {
        bool cons = HelicityConservation::conservedForMassless();
        ASSERT_TRUE(cons);
        return true;
    });

    run_test("Helicity conservation: not conserved for massive", []() {
        bool cons = HelicityConservation::conservedForMassive();
        ASSERT_TRUE(!cons);
        return true;
    });

    run_test("Helicity conservation: flip rate formula", []() {
        double rate = HelicityConservation::helicityFlipRate(0.001, 1.0);
        ASSERT_NEAR(rate, 1e-6, 1e-7);
        return true;
    });

    run_test("Weak interaction: V-A structure", []() {
        std::string va = WeakInteractionChirality::vMinusAStructure();
        ASSERT_TRUE(va.find("V-A") != std::string::npos);
        return true;
    });

    run_test("Weak interaction: charged current", []() {
        std::string cc = WeakInteractionChirality::chargedCurrent();
        ASSERT_TRUE(cc.find("W±") != std::string::npos);
        return true;
    });

    run_test("Weak interaction: neutral current", []() {
        std::string nc = WeakInteractionChirality::neutralCurrent();
        ASSERT_TRUE(nc.find("Z⁰") != std::string::npos);
        return true;
    });

    run_test("Weak interaction: neutrino coupling", []() {
        std::string nu = WeakInteractionChirality::neutrinoCoupling();
        ASSERT_TRUE(nu.find("left-handed") != std::string::npos);
        return true;
    });

    run_test("Weak interaction: beta decay", []() {
        std::string beta = WeakInteractionChirality::betaDecay();
        ASSERT_TRUE(beta.find("n → p") != std::string::npos);
        return true;
    });

    run_test("Weak interaction: parity violation", []() {
        std::string pv = WeakInteractionChirality::parityViolation();
        ASSERT_TRUE(pv.find("parity violation") != std::string::npos);
        return true;
    });

    run_test("Helicity amplitudes: spinor formalism", []() {
        std::string spin = HelicityAmplitudes::spinorHelicityFormalism();
        ASSERT_TRUE(spin.find("spinors") != std::string::npos);
        return true;
    });

    run_test("Chiral anomalies: axial anomaly", []() {
        std::string ax = ChiralAnomalies::axialAnomaly();
        ASSERT_TRUE(ax.find("anomaly") != std::string::npos);
        return true;
    });

    run_test("Chiral anomalies: pion decay", []() {
        std::string pi = ChiralAnomalies::pionDecay();
        ASSERT_TRUE(pi.find("π⁰ → γγ") != std::string::npos);
        return true;
    });

    run_test("Chiral anomalies: anomaly cancellation", []() {
        std::string cancel = ChiralAnomalies::anomalyCancellation();
        ASSERT_TRUE(cancel.find("cancel") != std::string::npos);
        return true;
    });

    std::cout << std::endl << "======================================" << std::endl;
    std::cout << "Phase 4 Results: Gauge Helicity" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    if (tests_failed == 0) {
        std::cout << std::endl << "All gauge helicity tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << std::endl << "Some tests FAILED." << std::endl;
        return 1;
    }
}
