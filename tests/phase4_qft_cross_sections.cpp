/**
 * Phase 4 Validation: QFT Cross Sections
 */

#include <iostream>
#include <cmath>
#include <string>
#include "../include/physics/qft_cross_sections.hpp"

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
    std::cout << "=== Phase 4: QFT Cross Sections Validation ===" << std::endl << std::endl;

    auto run_test = [&](const char* name, bool (*test_func)()) {
        std::cout << "Running: " << name << "... ";
        if (test_func()) { std::cout << "PASS" << std::endl; tests_passed++; return true; }
        else { tests_failed++; return false; }
    };

    run_test("CrossSection: unit conversions", []() {
        double barn = CrossSection::barn_to_m2();
        ASSERT_NEAR(barn, 1e-28, TOLERANCE);
        return true;
    });

    run_test("CrossSection: LHC luminosity", []() {
        double lum = CrossSection::lhc_luminosity();
        ASSERT_NEAR(lum, 1e34, TOLERANCE);
        return true;
    });

    run_test("Rutherford: differential cross section", []() {
        double dsigma = RutherfordScattering::differential(2, 79, 10.0, M_PI/4.0);
        ASSERT_TRUE(dsigma > 0.0);
        return true;
    });

    run_test("Rutherford: Mott correction", []() {
        double mott = RutherfordScattering::mottCorrection(M_PI/4.0);
        ASSERT_TRUE(mott > 0.0 && mott <= 1.0);
        return true;
    });

    run_test("QED: e+e- to mu+mu- cross section", []() {
        double sigma = QEDProcesses::electronMuonScattering(91.2);
        ASSERT_TRUE(sigma > 0.0);
        return true;
    });

    run_test("QED: Bhabha scattering", []() {
        double sigma = QEDProcesses::bhabhaScattering(91.2);
        ASSERT_TRUE(sigma > 0.0);
        return true;
    });

    run_test("QED: Compton scattering at low energy", []() {
        double sigma = QEDProcesses::comptonScattering(0.0001);
        ASSERT_TRUE(sigma > 0.0);
        return true;
    });

    run_test("Hadronic: R ratio below charm", []() {
        double R = HadronicCrossSections::Rratio(2.0);
        ASSERT_NEAR(R, 2.0, 0.5);
        return true;
    });

    run_test("Hadronic: R ratio above charm", []() {
        double R = HadronicCrossSections::Rratio(5.0);
        ASSERT_NEAR(R, 10.0/3.0, 0.5);
        return true;
    });

    run_test("Hadronic: R ratio above bottom", []() {
        double R = HadronicCrossSections::Rratio(15.0);
        ASSERT_NEAR(R, 11.0/3.0, 0.5);
        return true;
    });

    run_test("Hadronic: geometric cross section", []() {
        double sigma = HadronicCrossSections::geometricCrossSection(1e-15);
        ASSERT_TRUE(sigma > 0.0);
        return true;
    });

    run_test("Weak: neutrino-nucleon cross section", []() {
        double sigma = WeakProcesses::neutrinoNucleon(1.0);
        ASSERT_TRUE(sigma > 0.0);
        return true;
    });

    run_test("Weak: neutrino cross section grows with energy", []() {
        double sigma1 = WeakProcesses::neutrinoNucleon(1.0);
        double sigma2 = WeakProcesses::neutrinoNucleon(10.0);
        ASSERT_TRUE(sigma2 > sigma1);
        return true;
    });

    run_test("Weak: W boson production", []() {
        double sigma_W = WeakProcesses::wBosonProduction();
        ASSERT_TRUE(sigma_W > 0.0);
        return true;
    });

    run_test("Weak: Z boson production", []() {
        double sigma_Z = WeakProcesses::zBosonProduction();
        ASSERT_TRUE(sigma_Z > 0.0);
        return true;
    });

    run_test("PDF: Bjorken x calculation", []() {
        double x = PartonDistributionFunctions::bjorkenX(10.0, 100.0);
        ASSERT_TRUE(x > 0.0 && x < 1.0);
        return true;
    });

    run_test("PDF: valence quark distribution", []() {
        double pdf = PartonDistributionFunctions::valenceQuarkPDF(0.3);
        ASSERT_TRUE(pdf > 0.0);
        return true;
    });

    run_test("PDF: valence vanishes at boundaries", []() {
        double pdf_0 = PartonDistributionFunctions::valenceQuarkPDF(0.0);
        double pdf_1 = PartonDistributionFunctions::valenceQuarkPDF(1.0);
        ASSERT_NEAR(pdf_0, 0.0, TOLERANCE);
        ASSERT_NEAR(pdf_1, 0.0, TOLERANCE);
        return true;
    });

    run_test("PDF: gluon distribution at small x", []() {
        double pdf = PartonDistributionFunctions::gluonPDF(0.01);
        ASSERT_TRUE(pdf > 0.0);
        return true;
    });

    run_test("Physical: cross section decreases with energy", []() {
        double s1 = QEDProcesses::electronMuonScattering(10.0);
        double s2 = QEDProcesses::electronMuonScattering(100.0);
        ASSERT_TRUE(s2 < s1);
        return true;
    });

    std::cout << std::endl << "======================================" << std::endl;
    std::cout << "Phase 4 Results: QFT Cross Sections" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    if (tests_failed == 0) {
        std::cout << std::endl << "All QFT cross sections tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << std::endl << "Some tests FAILED." << std::endl;
        return 1;
    }
}
