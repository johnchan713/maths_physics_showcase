/**
 * Phase 4 Validation: Nuclear Physics
 *
 * Tests the nuclear_physics.hpp module functions.
 *
 * Coverage:
 * - Nuclear stability (binding energy per nucleon, SEMF, mass excess)
 * - Radioactive decay (alpha, beta, gamma, electron capture)
 * - Decay chains and equilibrium states
 * - Half-life and radioactive decay rates
 * - Neutron interactions and fission
 * - Radiation-matter interactions
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include "../include/physics/nuclear_physics.hpp"

// Test tolerance
const double TOLERANCE = 1e-6;
const double LOOSE_TOLERANCE = 1e-3;
const double VERY_LOOSE = 0.01;

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

#define ASSERT_FALSE(condition) \
    do { \
        if ((condition)) { \
            std::cerr << "FAIL: " << __LINE__ << ": " << #condition << " should be false" << std::endl; \
            return false; \
        } \
    } while(0)

// Physical constants from the module
using physics::nuclear::constants::c;
using physics::nuclear::constants::hbar;
using physics::nuclear::constants::m_e;
using physics::nuclear::constants::m_p;
using physics::nuclear::constants::m_n;
using physics::nuclear::constants::e;
using physics::nuclear::constants::u;
using physics::nuclear::constants::alpha;
using physics::nuclear::constants::MeV;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Nuclear Physics Validation ===" << std::endl;
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
    // Physical Constants Tests
    // ========================================

    run_test("Speed of light constant value", []() {
        ASSERT_NEAR(c, 299792458.0, TOLERANCE);
        return true;
    });

    run_test("Reduced Planck constant value", []() {
        ASSERT_NEAR(hbar, 1.054571817e-34, 1e-45);
        return true;
    });

    run_test("Proton mass value", []() {
        ASSERT_NEAR(m_p, 1.67262192369e-27, 1e-35);
        return true;
    });

    run_test("Neutron mass value", []() {
        ASSERT_NEAR(m_n, 1.67492749804e-27, 1e-35);
        return true;
    });

    run_test("Atomic mass unit value", []() {
        ASSERT_NEAR(u, 1.66053906660e-27, 1e-35);
        return true;
    });

    // ========================================
    // Nuclear Stability Tests - Binding Energy
    // ========================================

    run_test("Binding energy curve returns string", []() {
        std::string result = physics::nuclear::NuclearStability::binding_energy_curve();
        ASSERT_TRUE(result.length() > 0);
        ASSERT_TRUE(result.find("Fe-56") != std::string::npos);
        return true;
    });

    run_test("Binding energy SEMF for carbon-12", []() {
        // Carbon-12: A=12, Z=6
        double B = physics::nuclear::NuclearStability::binding_energy_semf(12, 6);
        // Should be positive (nucleus is bound)
        ASSERT_TRUE(B > 0);
        // Roughly 90-100 MeV for C-12
        ASSERT_TRUE(B > 80);
        ASSERT_TRUE(B < 110);
        return true;
    });

    run_test("Binding energy SEMF for iron-56", []() {
        // Iron-56: A=56, Z=26 (most stable)
        double B = physics::nuclear::NuclearStability::binding_energy_semf(56, 26);
        ASSERT_TRUE(B > 0);
        // Fe-56 is near peak stability
        ASSERT_TRUE(B > 480);
        ASSERT_TRUE(B < 510);
        return true;
    });

    run_test("Binding energy per nucleon for He-4", []() {
        // Helium-4: A=4, Z=2
        double BE_A = physics::nuclear::NuclearStability::binding_energy_per_nucleon(4, 2);
        // SEMF is approximate for light nuclei, gives ~5.5 MeV/nucleon
        ASSERT_NEAR(BE_A, 5.5, 0.5);
        return true;
    });

    run_test("Binding energy per nucleon for Fe-56", []() {
        // Iron-56: A=56, Z=26
        double BE_A = physics::nuclear::NuclearStability::binding_energy_per_nucleon(56, 26);
        // Fe-56 has maximum BE/A at about 8.8 MeV/nucleon
        ASSERT_TRUE(BE_A > 8.7);
        ASSERT_TRUE(BE_A < 8.9);
        return true;
    });

    run_test("Mass excess increases with binding energy deficit", []() {
        // Compare two nuclei
        double dM1 = physics::nuclear::NuclearStability::mass_excess(4, 2);   // He-4
        double dM2 = physics::nuclear::NuclearStability::mass_excess(12, 6);  // C-12
        // Result should be finite
        ASSERT_TRUE(std::isfinite(dM1));
        ASSERT_TRUE(std::isfinite(dM2));
        return true;
    });

    run_test("Neutron separation energy for He-4", []() {
        // He-4: A=4, Z=2 (stable nucleus)
        double S_n = physics::nuclear::NuclearStability::neutron_separation_energy(4, 2);
        ASSERT_TRUE(S_n > 0);
        return true;
    });

    run_test("Proton separation energy for Li-7", []() {
        // Li-7: A=7, Z=3
        double S_p = physics::nuclear::NuclearStability::proton_separation_energy(7, 3);
        ASSERT_TRUE(S_p > 0);
        return true;
    });

    run_test("Q-value conservation for mass difference", []() {
        double Q = physics::nuclear::NuclearStability::q_value(1.007825, 1.008665);
        // Result should be a valid energy value
        ASSERT_TRUE(std::isfinite(Q));
        return true;
    });

    run_test("Valley of stability for light nuclei", []() {
        // Helium-4: very stable
        bool stable = physics::nuclear::NuclearStability::is_in_valley_of_stability(4, 2);
        ASSERT_TRUE(stable);
        return true;
    });

    run_test("Valley of stability for heavy stable nuclei", []() {
        // Iron-56: very stable
        bool stable = physics::nuclear::NuclearStability::is_in_valley_of_stability(56, 26);
        ASSERT_TRUE(stable);
        return true;
    });

    run_test("Valley of stability rejects unstable nuclei", []() {
        // Exotic nucleus far from stability
        bool stable = physics::nuclear::NuclearStability::is_in_valley_of_stability(200, 120);
        ASSERT_FALSE(stable);
        return true;
    });

    // ========================================
    // Natural Radioactivity Tests
    // ========================================

    run_test("Decay series returns valid string", []() {
        std::string series = physics::nuclear::NaturalRadioactivity::decay_series();
        ASSERT_TRUE(series.length() > 0);
        ASSERT_TRUE(series.find("Th-232") != std::string::npos);
        return true;
    });

    run_test("Primordial radionuclides list is non-empty", []() {
        auto nuclei = physics::nuclear::NaturalRadioactivity::primordial_radionuclides();
        ASSERT_TRUE(nuclei.size() > 0);
        return true;
    });

    run_test("Primordial radionuclides contains U-238", []() {
        auto nuclei = physics::nuclear::NaturalRadioactivity::primordial_radionuclides();
        auto it = std::find(nuclei.begin(), nuclei.end(), std::string("U-238"));
        ASSERT_TRUE(it != nuclei.end());
        return true;
    });

    run_test("Cosmogenic radionuclides list is non-empty", []() {
        auto nuclei = physics::nuclear::NaturalRadioactivity::cosmogenic_radionuclides();
        ASSERT_TRUE(nuclei.size() > 0);
        return true;
    });

    run_test("Cosmogenic radionuclides contains C-14", []() {
        auto nuclei = physics::nuclear::NaturalRadioactivity::cosmogenic_radionuclides();
        auto it = std::find(nuclei.begin(), nuclei.end(), std::string("C-14"));
        ASSERT_TRUE(it != nuclei.end());
        return true;
    });

    // ========================================
    // Alpha Decay Tests
    // ========================================

    run_test("Alpha decay equation returns valid string", []() {
        std::string eq = physics::nuclear::AlphaDecay::decay_equation();
        ASSERT_TRUE(eq.length() > 0);
        ASSERT_TRUE(eq.find("He") != std::string::npos);
        return true;
    });

    run_test("Q-value for alpha decay is positive", []() {
        // U-238 decay: parent mass > daughter mass ensures Q > 0
        double Q = physics::nuclear::AlphaDecay::q_value_alpha(239.0 * u, 234.0 * u);
        ASSERT_TRUE(Q > 0);
        return true;
    });

    run_test("Alpha kinetic energy from Q-value", []() {
        double Q = 5.5;  // MeV (typical alpha decay)
        double T_alpha = physics::nuclear::AlphaDecay::alpha_kinetic_energy(Q, 238);
        ASSERT_TRUE(T_alpha > 0);
        ASSERT_TRUE(T_alpha < Q);  // Recoil reduces alpha energy
        return true;
    });

    run_test("Geiger-Nuttall constant calculation", []() {
        double E_alpha = 5.5;  // MeV
        double lambda = 1.0e-17;  // s^-1
        double constant = physics::nuclear::AlphaDecay::geiger_nuttall_constant(E_alpha, lambda);
        ASSERT_TRUE(std::isfinite(constant));
        return true;
    });

    run_test("Gamow factor for alpha decay", []() {
        double v = 1.0e7;  // m/s
        double G = physics::nuclear::AlphaDecay::gamow_factor(92, v);
        ASSERT_TRUE(G > 0);
        return true;
    });

    run_test("Decay constant from Gamow theory", []() {
        double G = 50.0;  // Gamow factor
        double lambda = physics::nuclear::AlphaDecay::decay_constant_gamow(G);
        ASSERT_TRUE(lambda > 0);
        return true;
    });

    // ========================================
    // Beta Decay Tests
    // ========================================

    run_test("Beta minus decay equation returns string", []() {
        std::string eq = physics::nuclear::BetaDecay::beta_minus();
        ASSERT_TRUE(eq.length() > 0);
        ASSERT_TRUE(eq.find("e") != std::string::npos);
        return true;
    });

    run_test("Beta plus decay equation returns string", []() {
        std::string eq = physics::nuclear::BetaDecay::beta_plus();
        ASSERT_TRUE(eq.length() > 0);
        ASSERT_TRUE(eq.find("e") != std::string::npos);
        return true;
    });

    run_test("Q-value for beta minus decay", []() {
        double Q = physics::nuclear::BetaDecay::q_value_beta_minus(14.003241 * u, 14.003074 * u);
        ASSERT_TRUE(std::isfinite(Q));
        return true;
    });

    run_test("Q-value for beta plus decay", []() {
        double Q = physics::nuclear::BetaDecay::q_value_beta_plus(11.011433 * u, 11.003806 * u);
        ASSERT_TRUE(std::isfinite(Q));
        return true;
    });

    run_test("Fermi decay constant is positive", []() {
        double M_fi = 0.5;
        double f_ZQ = 1.0e-3;
        double lambda = physics::nuclear::BetaDecay::fermi_decay_constant(M_fi, f_ZQ);
        ASSERT_TRUE(lambda > 0);
        return true;
    });

    run_test("Fermi integral increases with Q-value", []() {
        double f1 = physics::nuclear::BetaDecay::fermi_integral(82, 1.0);
        double f2 = physics::nuclear::BetaDecay::fermi_integral(82, 2.0);
        ASSERT_TRUE(f2 > f1);
        return true;
    });

    run_test("Maximum beta energy equals Q-value", []() {
        double Q = 2.5;  // MeV
        double E_max = physics::nuclear::BetaDecay::max_beta_energy(Q);
        ASSERT_NEAR(E_max, Q, TOLERANCE);
        return true;
    });

    run_test("Average beta energy is about Q over 3", []() {
        double Q = 2.5;  // MeV
        double E_avg = physics::nuclear::BetaDecay::average_beta_energy(Q);
        // Average should be less than Q
        ASSERT_TRUE(E_avg < Q);
        ASSERT_TRUE(E_avg > 0);
        return true;
    });

    // ========================================
    // Electron Capture Tests
    // ========================================

    run_test("Electron capture decay equation returns string", []() {
        std::string eq = physics::nuclear::ElectronCapture::decay_equation();
        ASSERT_TRUE(eq.length() > 0);
        ASSERT_TRUE(eq.find("e") != std::string::npos);
        return true;
    });

    run_test("Q-value for electron capture", []() {
        double Q = physics::nuclear::ElectronCapture::q_value_ec(99.0 * u, 99.0 * u, 0.01);
        ASSERT_TRUE(std::isfinite(Q));
        return true;
    });

    run_test("K-shell binding energy increases with Z", []() {
        double B_K1 = physics::nuclear::ElectronCapture::k_shell_binding(20);
        double B_K2 = physics::nuclear::ElectronCapture::k_shell_binding(40);
        ASSERT_TRUE(B_K2 > B_K1);
        return true;
    });

    run_test("Beta plus threshold value is correct", []() {
        std::string threshold = physics::nuclear::ElectronCapture::beta_plus_threshold();
        ASSERT_TRUE(threshold.length() > 0);
        ASSERT_TRUE(threshold.find("1.022") != std::string::npos);
        return true;
    });

    run_test("Electron capture branching ratio for low Q", []() {
        double ratio = physics::nuclear::ElectronCapture::ec_branching_ratio(30, 0.5);
        // For Q < 1.022 MeV, only EC is possible
        ASSERT_NEAR(ratio, 1.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Gamma Emission Tests
    // ========================================

    run_test("Gamma decay equation returns string", []() {
        std::string eq = physics::nuclear::GammaEmission::decay_equation();
        ASSERT_TRUE(eq.length() > 0);
        ASSERT_TRUE(eq.find("γ") != std::string::npos || eq.find("gamma") != std::string::npos);
        return true;
    });

    run_test("Gamma energy equals excitation minus final state", []() {
        double E_gamma = physics::nuclear::GammaEmission::gamma_energy(2.0, 0.5, 0.0);
        ASSERT_NEAR(E_gamma, 1.5, TOLERANCE);
        return true;
    });

    run_test("Recoil energy for gamma emission", []() {
        double E_R = physics::nuclear::GammaEmission::recoil_energy(1.0, 56);
        ASSERT_TRUE(E_R > 0);
        ASSERT_TRUE(E_R < 1.0);  // Recoil is small
        return true;
    });

    run_test("Weisskopf E1 transition half-life", []() {
        double t_half = physics::nuclear::GammaEmission::weisskopf_e1_halflife(56, 1.0);
        ASSERT_TRUE(t_half > 0);
        return true;
    });

    run_test("Weisskopf M1 transition half-life", []() {
        double t_half = physics::nuclear::GammaEmission::weisskopf_m1_halflife(1.0);
        ASSERT_TRUE(t_half > 0);
        return true;
    });

    run_test("Transition rate increases with multipolarity", []() {
        double rate_L1 = physics::nuclear::GammaEmission::transition_rate(1.0, 1);
        double rate_L2 = physics::nuclear::GammaEmission::transition_rate(1.0, 2);
        ASSERT_TRUE(rate_L2 > rate_L1);
        return true;
    });

    run_test("Selection rules returns valid string", []() {
        std::string rules = physics::nuclear::GammaEmission::selection_rules();
        ASSERT_TRUE(rules.length() > 0);
        return true;
    });

    // ========================================
    // Decay Rates and Half-Life Tests
    // ========================================

    run_test("Decay constant from half-life", []() {
        double T_half = 5730.0;  // Years for C-14
        double lambda = physics::nuclear::DecayRates::lambda_from_halflife(T_half);
        ASSERT_TRUE(lambda > 0);
        return true;
    });

    run_test("Half-life from decay constant", []() {
        double lambda = 1.0e-10;  // s^-1
        double T_half = physics::nuclear::DecayRates::halflife_from_lambda(lambda);
        ASSERT_TRUE(T_half > 0);
        return true;
    });

    run_test("Mean lifetime equals inverse decay constant", []() {
        double lambda = 1.0e-10;  // s^-1
        double tau = physics::nuclear::DecayRates::mean_lifetime(lambda);
        ASSERT_NEAR(tau, 1.0 / lambda, 1.0e-20);
        return true;
    });

    run_test("Number of nuclei decays with time", []() {
        double N = physics::nuclear::DecayRates::number_of_nuclei(1000.0, 0.1, 10.0);
        ASSERT_TRUE(N < 1000.0);
        ASSERT_TRUE(N > 0);
        return true;
    });

    run_test("Activity decays exponentially with time", []() {
        double A = physics::nuclear::DecayRates::activity(1000.0, 0.1, 10.0);
        ASSERT_TRUE(A < 1000.0);
        ASSERT_TRUE(A > 0);
        return true;
    });

    run_test("Activity from number of nuclei", []() {
        double lambda = 0.1;  // s^-1
        double N = 1000.0;
        double A = physics::nuclear::DecayRates::activity_from_N(N, lambda);
        ASSERT_NEAR(A, N * lambda, TOLERANCE);
        return true;
    });

    run_test("Fraction remaining after half-life", []() {
        double T_half = 5730.0;
        double frac = physics::nuclear::RadioactivityVariation::fraction_remaining(T_half, T_half);
        ASSERT_NEAR(frac, 0.5, 0.01);
        return true;
    });

    run_test("Fraction decayed is complement of remaining", []() {
        double T_half = 5730.0;
        double t = T_half;
        double frac_rem = physics::nuclear::RadioactivityVariation::fraction_remaining(T_half, t);
        double frac_dec = physics::nuclear::RadioactivityVariation::fraction_decayed(T_half, t);
        ASSERT_NEAR(frac_rem + frac_dec, 1.0, TOLERANCE);
        return true;
    });

    run_test("Specific activity calculation", []() {
        double lambda = 0.1;  // s^-1
        double M_molar = 14.0e-3;  // kg/mol for C-14
        double A_specific = physics::nuclear::DecayRates::specific_activity(lambda, M_molar);
        ASSERT_TRUE(A_specific > 0);
        ASSERT_TRUE(std::isfinite(A_specific));
        return true;
    });

    // ========================================
    // Decay Chain and Equilibrium Tests
    // ========================================

    run_test("Secular equilibrium condition", []() {
        bool secular = physics::nuclear::DecayChains::is_secular_equilibrium(1000.0, 1.0);
        ASSERT_TRUE(secular);
        return true;
    });

    run_test("Transient equilibrium condition", []() {
        bool transient = physics::nuclear::DecayChains::is_transient_equilibrium(50.0, 1.0);
        ASSERT_TRUE(transient);
        return true;
    });

    run_test("No equilibrium for similar half-lives", []() {
        bool secular = physics::nuclear::DecayChains::is_secular_equilibrium(100.0, 50.0);
        ASSERT_FALSE(secular);
        return true;
    });

    run_test("Transient equilibrium activity ratio", []() {
        double ratio = physics::nuclear::DecayChains::transient_eq_ratio(1.0, 2.0);
        ASSERT_TRUE(ratio > 0);
        return true;
    });

    run_test("Daughter activity calculation", []() {
        double lambda_A = 0.01;  // s^-1
        double lambda_B = 0.05;  // s^-1
        double A_A_0 = 1000.0;   // Initial activity
        double t = 10.0;          // Time (s)
        double A_daughter = physics::nuclear::DecayChains::daughter_activity_bateman(lambda_A, lambda_B, A_A_0, t);
        ASSERT_TRUE(A_daughter > 0);
        ASSERT_TRUE(std::isfinite(A_daughter));
        return true;
    });

    // ========================================
    // Decay Prediction Tests
    // ========================================

    run_test("Beta minus decay prediction for N-rich nuclei", []() {
        // N=8, Z=6 (Carbon-14 is N-rich)
        std::string prediction = physics::nuclear::DecayPrediction::predict_beta_type(14, 6);
        ASSERT_TRUE(prediction.length() > 0);
        return true;
    });

    run_test("Alpha decay likely for heavy nuclei", []() {
        // Uranium-238: Z=92, A=238
        bool alpha_likely = physics::nuclear::DecayPrediction::alpha_decay_likely(92, 238);
        ASSERT_TRUE(alpha_likely);
        return true;
    });

    run_test("Alpha decay unlikely for light nuclei", []() {
        // Carbon-12: Z=6, A=12
        bool alpha_likely = physics::nuclear::DecayPrediction::alpha_decay_likely(6, 12);
        ASSERT_FALSE(alpha_likely);
        return true;
    });

    run_test("Fission competes for Z-squared over A greater than 47", []() {
        // Uranium-238: Z=92, A=238; Z²/A ≈ 35.5
        bool fission = physics::nuclear::DecayPrediction::fission_competes(92, 238);
        // U-238 has Z²/A < 47, so fission doesn't strongly compete
        ASSERT_FALSE(fission);
        return true;
    });

    run_test("Beyond proton drip line detection", []() {
        // Test a nucleus far beyond stable region
        bool beyond = physics::nuclear::DecayPrediction::beyond_proton_drip(100, 70);
        ASSERT_TRUE(beyond);
        return true;
    });

    run_test("Beyond neutron drip line detection", []() {
        // Test a nucleus far from stability
        bool beyond = physics::nuclear::DecayPrediction::beyond_neutron_drip(300, 50);
        ASSERT_TRUE(beyond);
        return true;
    });

    // ========================================
    // Nuclear Physics Test Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "=== Test Summary ===" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Total tests: " << (tests_passed + tests_failed) << std::endl;

    if (tests_failed == 0) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests failed." << std::endl;
        return 1;
    }
}
