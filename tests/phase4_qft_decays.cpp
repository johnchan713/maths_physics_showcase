/**
 * Phase 4 Validation: QFT Decays
 */

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "../include/physics/qft_decays.hpp"

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
    std::cout << "=== Phase 4: QFT Decays Validation ===" << std::endl << std::endl;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) { std::cout << "PASS" << std::endl; tests_passed++; return true; }
        else { tests_failed++; return false; }
    };

    run_test("Decay rate: lifetime from width", []() {
        double tau = DecayRate::lifetimeFromWidth(1e-15);
        ASSERT_TRUE(tau > 0.0);
        return true;
    });

    run_test("Decay rate: width from lifetime", []() {
        double width = DecayRate::widthFromLifetime(1e-10);
        ASSERT_TRUE(width > 0.0);
        return true;
    });

    run_test("Decay rate: exponential decay", []() {
        double N = DecayRate::exponentialDecay(1000.0, 1.0, 1.0);
        ASSERT_NEAR(N, 1000.0 / M_E, 1.0);
        return true;
    });

    run_test("Decay rate: half-life calculation", []() {
        double t_half = DecayRate::halfLife(1.0);
        ASSERT_NEAR(t_half, std::log(2.0), TOLERANCE);
        return true;
    });

    run_test("Decay rate: decay length calculation", []() {
        double L = DecayRate::decayLength(0.99, 10.0, 1e-10);
        ASSERT_TRUE(L > 0.0);
        return true;
    });

    run_test("Branching ratio: calculation", []() {
        double br = BranchingRatio::branchingFraction(0.5, 1.0);
        ASSERT_NEAR(br, 0.5, TOLERANCE);
        return true;
    });

    run_test("Branching ratio: total width from partials", []() {
        std::vector<double> widths = {0.3, 0.5, 0.2};
        double total = BranchingRatio::totalWidth(widths);
        ASSERT_NEAR(total, 1.0, TOLERANCE);
        return true;
    });

    run_test("Branching ratio: normalization check", []() {
        std::vector<double> brs = {0.5, 0.3, 0.2};
        bool normalized = BranchingRatio::checkNormalization(brs);
        ASSERT_TRUE(normalized);
        return true;
    });

    run_test("Fermis Golden Rule: two-body rate", []() {
        double rate = FermisGoldenRule::twoBodyRate(1.0, 100.0, 50.0);
        ASSERT_TRUE(rate > 0.0);
        return true;
    });

    run_test("Fermis Golden Rule: two-body momentum", []() {
        double p = FermisGoldenRule::twoBodyMomentum(100.0, 10.0, 20.0);
        ASSERT_TRUE(p > 0.0);
        return true;
    });

    run_test("Fermis Golden Rule: kinematically forbidden", []() {
        double p = FermisGoldenRule::twoBodyMomentum(10.0, 20.0, 30.0);
        ASSERT_NEAR(p, 0.0, TOLERANCE);
        return true;
    });

    run_test("Resonances: Breit-Wigner at peak", []() {
        double bw = Resonances::breitWigner(91.2, 91.2, 2.5);
        ASSERT_TRUE(bw > 0.0);
        return true;
    });

    run_test("Resonances: Breit-Wigner off peak", []() {
        double bw_peak = Resonances::breitWigner(91.2, 91.2, 2.5);
        double bw_off = Resonances::breitWigner(100.0, 91.2, 2.5);
        ASSERT_TRUE(bw_off < bw_peak);
        return true;
    });

    run_test("Resonances: FWHM equals width", []() {
        double fwhm = Resonances::fwhm(2.5);
        ASSERT_NEAR(fwhm, 2.5, TOLERANCE);
        return true;
    });

    run_test("Weak decays: muon decay width", []() {
        double width = WeakDecays::muonDecayWidth();
        ASSERT_TRUE(width > 0.0);
        return true;
    });

    run_test("Weak decays: muon lifetime", []() {
        double tau = WeakDecays::muonLifetime();
        ASSERT_NEAR(tau, 2.2e-6, 0.5e-6);
        return true;
    });

    run_test("Weak decays: neutron lifetime", []() {
        double tau = WeakDecays::neutronLifetime();
        ASSERT_NEAR(tau, 880.0, 20.0);
        return true;
    });

    run_test("Weak decays: pion lifetime", []() {
        double tau = WeakDecays::pionLifetime();
        ASSERT_TRUE(tau > 0.0 && tau < 1e-7);
        return true;
    });

    run_test("Weak decays: kaon lifetime", []() {
        double tau = WeakDecays::kaonLifetime();
        ASSERT_TRUE(tau > 0.0 && tau < 1e-7);
        return true;
    });

    run_test("Weak decays: CKM Vud", []() {
        double Vud = WeakDecays::ckmElement(1, 1);
        ASSERT_NEAR(Vud, 0.974, 0.01);
        return true;
    });

    run_test("Weak decays: CKM Vus", []() {
        double Vus = WeakDecays::ckmElement(1, 2);
        ASSERT_NEAR(Vus, 0.225, 0.05);
        return true;
    });

    run_test("Strong decays: rho meson width", []() {
        double width = StrongDecays::rhoMesonWidth();
        ASSERT_NEAR(width, 0.150, 0.05);
        return true;
    });

    run_test("Strong decays: Delta resonance width", []() {
        double width = StrongDecays::deltaResonanceWidth();
        ASSERT_NEAR(width, 0.120, 0.05);
        return true;
    });

    run_test("Strong decays: typical lifetime is short", []() {
        double tau = StrongDecays::typicalLifetime();
        ASSERT_TRUE(tau < 1e-22);
        return true;
    });

    run_test("EM decays: pi0 lifetime", []() {
        double tau = ElectromagneticDecays::neutralPionLifetime();
        ASSERT_TRUE(tau > 1e-18 && tau < 1e-15);
        return true;
    });

    run_test("EM decays: eta meson lifetime", []() {
        double tau = ElectromagneticDecays::etaMesonLifetime();
        ASSERT_TRUE(tau > 0.0);
        return true;
    });

    run_test("EM decays: suppression factor", []() {
        double supp = ElectromagneticDecays::emSuppressionFactor();
        ASSERT_TRUE(supp < 1e-4);
        return true;
    });

    run_test("Physical: strong decays fastest", []() {
        double tau_strong = StrongDecays::typicalLifetime();
        double tau_em = ElectromagneticDecays::neutralPionLifetime();
        double tau_weak = WeakDecays::pionLifetime();
        ASSERT_TRUE(tau_strong < tau_em);
        ASSERT_TRUE(tau_em < tau_weak);
        return true;
    });

    run_test("Physical: lifetime inversely proportional to width", []() {
        double width1 = 1e-15;
        double width2 = 2e-15;
        double tau1 = DecayRate::lifetimeFromWidth(width1);
        double tau2 = DecayRate::lifetimeFromWidth(width2);
        ASSERT_NEAR(tau1 / tau2, 2.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Physical: branching ratios sum to 1", []() {
        std::vector<double> brs = {0.3, 0.4, 0.2, 0.1};
        bool normalized = BranchingRatio::checkNormalization(brs);
        ASSERT_TRUE(normalized);
        return true;
    });

    std::cout << std::endl << "======================================" << std::endl;
    std::cout << "Phase 4 Results: QFT Decays" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    if (tests_failed == 0) {
        std::cout << std::endl << "All QFT decays tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << std::endl << "Some tests FAILED." << std::endl;
        return 1;
    }
}
