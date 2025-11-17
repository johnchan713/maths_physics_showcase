#include <stdexcept>  // Required for std::invalid_argument used in cosmology_early_universe.hpp
#include "physics/cosmology_early_universe.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <iomanip>

// Test counter
int tests_passed = 0;
int tests_failed = 0;

// Tolerance for floating-point comparisons
const double TOLERANCE = 1e-6;
const double LOOSE_TOLERANCE = 1e-3;
const double VERY_LOOSE_TOLERANCE = 0.01;

// Test macros
#define ASSERT_NEAR(actual, expected, tol) \
    do { \
        double diff = std::abs((actual) - (expected)); \
        if (diff <= (tol)) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cerr << "FAIL: " << __LINE__ << ": " \
                      << #actual << " = " << (actual) \
                      << ", expected " << (expected) \
                      << " (diff: " << diff << ")" << std::endl; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        if (condition) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cerr << "FAIL: " << __LINE__ << ": " \
                      << #condition << " is false" << std::endl; \
        } \
    } while(0)

#define RUN_TEST(name, code) \
    do { \
        std::cout << "Running: " << name << "... "; \
        code; \
        std::cout << "PASS" << std::endl; \
    } while(0)

using namespace physics::advanced::cosmology;

int main() {
    std::cout << "=== Phase 4: Early Universe Cosmology ===" << std::endl;
    std::cout << std::endl;

    // ========================================
    // Cosmic Microwave Background Tests
    // ========================================

    RUN_TEST("CMB temperature today (2.7255 K)", {
        double T_CMB = CosmicMicrowaveBackground::temperatureToday();
        ASSERT_NEAR(T_CMB, 2.7255, TOLERANCE);
    });

    RUN_TEST("CMB temperature at redshift", {
        double z = 1.0;
        double T = CosmicMicrowaveBackground::temperatureAtRedshift(z);
        double expected = 2.7255 * 2.0;  // T(z) = T₀(1+z)
        ASSERT_NEAR(T, expected, TOLERANCE);
    });

    RUN_TEST("CMB temperature scales with (1+z)", {
        double z = 10.0;
        double T = CosmicMicrowaveBackground::temperatureAtRedshift(z);
        double T0 = CosmicMicrowaveBackground::temperatureToday();
        ASSERT_NEAR(T / T0, 1.0 + z, TOLERANCE);
    });

    RUN_TEST("CMB decoupling redshift z ~ 1100", {
        double z_dec = CosmicMicrowaveBackground::decouplingRedshift();
        ASSERT_NEAR(z_dec, 1100.0, TOLERANCE);
    });

    RUN_TEST("CMB temperature at decoupling ~3000 K", {
        double T_dec = CosmicMicrowaveBackground::temperatureAtDecoupling();
        double expected = 2.7255 * 1101.0;  // T₀(1 + z_dec)
        ASSERT_NEAR(T_dec, expected, 1.0);  // ~3000 K
    });

    RUN_TEST("CMB decoupling time 380,000 years", {
        double t_dec = CosmicMicrowaveBackground::decouplingTime();
        ASSERT_NEAR(t_dec, 380000.0, TOLERANCE);
    });

    RUN_TEST("CMB energy density today (Stefan-Boltzmann)", {
        double rho_gamma = CosmicMicrowaveBackground::energyDensityToday();
        ASSERT_TRUE(rho_gamma > 0.0);
        ASSERT_TRUE(rho_gamma < 1e-12);  // Very small today (J/m³)
    });

    RUN_TEST("Planck spectrum is positive", {
        double T = 2.7255;  // K
        double freq = 160e9;  // Hz (peak of CMB)
        double B = CosmicMicrowaveBackground::planckSpectrum(freq, T);
        ASSERT_TRUE(B > 0.0);
    });

    RUN_TEST("Planck spectrum goes to zero at high frequency", {
        double T = 2.7255;
        double freq = 1e15;  // Very high frequency
        double B = CosmicMicrowaveBackground::planckSpectrum(freq, T);
        ASSERT_TRUE(B < 1e-50);  // Essentially zero
    });

    RUN_TEST("CMB anisotropies ΔT/T ~ 10⁻⁵", {
        double aniso = CosmicMicrowaveBackground::anisotropyLevel();
        ASSERT_NEAR(aniso, 1e-5, TOLERANCE);
    });

    RUN_TEST("CMB acoustic peaks exist", {
        auto peaks = CosmicMicrowaveBackground::acousticPeakAngularScales();
        ASSERT_TRUE(peaks.size() >= 3);
        ASSERT_TRUE(peaks[0] > peaks[1]);  // First peak at larger scale
        ASSERT_TRUE(peaks[1] > peaks[2]);
    });

    RUN_TEST("Silk damping scale ~1 Mpc", {
        double lambda_D = CosmicMicrowaveBackground::silkDampingScale();
        ASSERT_NEAR(lambda_D, 1.0, TOLERANCE);
    });

    // ========================================
    // Radiation Era Tests
    // ========================================

    RUN_TEST("Matter-radiation equality redshift z_eq ~ 3400", {
        double z_eq = RadiationEra::equalityRedshift();
        ASSERT_NEAR(z_eq, 3400.0, TOLERANCE);
    });

    RUN_TEST("Matter-radiation equality time t_eq ~ 47,000 years", {
        double t_eq = RadiationEra::equalityTime();
        ASSERT_NEAR(t_eq, 47000.0, TOLERANCE);
    });

    RUN_TEST("Radiation era: a(t) ∝ t^(1/2)", {
        double t_ratio = 4.0;
        double a_ratio = RadiationEra::scaleFactorEvolution(t_ratio);
        double expected = std::sqrt(4.0);  // 2.0
        ASSERT_NEAR(a_ratio, expected, TOLERANCE);
    });

    RUN_TEST("Radiation era: T ∝ t^(-1/2)", {
        double t = 1.0;  // seconds
        double T_ref = 10.0e9;  // 10 GK
        double t_ref = 1.0;
        double T = RadiationEra::temperatureEvolution(t, T_ref, t_ref);
        ASSERT_NEAR(T, T_ref, TOLERANCE);

        // At t = 4 s, T should be half
        double t2 = 4.0;
        double T2 = RadiationEra::temperatureEvolution(t2, T_ref, t_ref);
        ASSERT_NEAR(T2, T_ref / 2.0, LOOSE_TOLERANCE);
    });

    RUN_TEST("Radiation era: H = 1/(2t)", {
        double t = 1.0;  // seconds
        double H = RadiationEra::hubbleParameter(t);
        double expected = 0.5;  // s⁻¹
        ASSERT_NEAR(H, expected, TOLERANCE);
    });

    RUN_TEST("Radiation era: age = 1/(2H)", {
        double t = 100.0;
        double H = RadiationEra::hubbleParameter(t);
        double age = 1.0 / (2.0 * H);
        ASSERT_NEAR(age, t, TOLERANCE);
    });

    RUN_TEST("Radiation energy density positive", {
        double T = 1e10;  // 10 GK
        double rho_r = RadiationEra::energyDensity(T);
        ASSERT_TRUE(rho_r > 0.0);
    });

    RUN_TEST("Radiation energy density ∝ T⁴", {
        double T1 = 1e10;
        double T2 = 2e10;
        double rho1 = RadiationEra::energyDensity(T1);
        double rho2 = RadiationEra::energyDensity(T2);
        double ratio = rho2 / rho1;
        ASSERT_NEAR(ratio, 16.0, LOOSE_TOLERANCE);  // (T2/T1)⁴ = 2⁴ = 16
    });

    RUN_TEST("Effective degrees of freedom at high T", {
        double g_high = RadiationEra::effectiveDegreesOfFreedom(1000.0);  // 1 TeV
        ASSERT_NEAR(g_high, 106.75, TOLERANCE);  // All SM particles
    });

    RUN_TEST("Effective degrees of freedom today", {
        double g_today = RadiationEra::effectiveDegreesOfFreedom(0.0001);  // < 0.5 MeV
        ASSERT_NEAR(g_today, 3.36, TOLERANCE);  // Photons + neutrinos
    });

    RUN_TEST("g_* decreases as universe cools", {
        double g_high = RadiationEra::effectiveDegreesOfFreedom(301.0);  // > 300 MeV
        double g_mid = RadiationEra::effectiveDegreesOfFreedom(10.0);    // 1-300 MeV
        double g_low = RadiationEra::effectiveDegreesOfFreedom(0.1);     // < 0.5 MeV
        ASSERT_TRUE(g_high > g_mid);
        ASSERT_TRUE(g_mid > g_low);
    });

    // ========================================
    // Matter Era Tests
    // ========================================

    RUN_TEST("Matter era starts at z_eq", {
        double z_start = MatterEra::startRedshift();
        double z_eq = RadiationEra::equalityRedshift();
        ASSERT_NEAR(z_start, z_eq, TOLERANCE);
    });

    RUN_TEST("Matter era ends at z ~ 0.4", {
        double z_end = MatterEra::endRedshift();
        ASSERT_NEAR(z_end, 0.4, TOLERANCE);
    });

    RUN_TEST("Matter era: a(t) ∝ t^(2/3)", {
        double t_ratio = 8.0;
        double a_ratio = MatterEra::scaleFactorEvolution(t_ratio);
        double expected = std::pow(8.0, 2.0/3.0);  // 4.0
        ASSERT_NEAR(a_ratio, expected, TOLERANCE);
    });

    RUN_TEST("Matter era: H = 2/(3t)", {
        double t = 1.0;  // seconds
        double H = MatterEra::hubbleParameter(t);
        double expected = 2.0 / 3.0;  // s⁻¹
        ASSERT_NEAR(H, expected, TOLERANCE);
    });

    RUN_TEST("Matter era: age = 2/(3H)", {
        double t = 100.0;
        double H = MatterEra::hubbleParameter(t);
        double age = 2.0 / (3.0 * H);
        ASSERT_NEAR(age, t, TOLERANCE);
    });

    RUN_TEST("Matter era: density perturbations grow as δ ∝ a", {
        double a = 0.5;
        double delta = MatterEra::densityPerturbationGrowth(a);
        ASSERT_NEAR(delta, a, TOLERANCE);
    });

    // ========================================
    // Big Bang Nucleosynthesis Tests
    // ========================================

    RUN_TEST("BBN temperature range 1-0.1 MeV", {
        auto T_range = BigBangNucleosynthesis::temperatureRange();
        ASSERT_NEAR(T_range.first, 1.0, TOLERANCE);
        ASSERT_NEAR(T_range.second, 0.1, TOLERANCE);
    });

    RUN_TEST("BBN time range 1-180 seconds", {
        auto t_range = BigBangNucleosynthesis::timeRange();
        ASSERT_NEAR(t_range.first, 1.0, TOLERANCE);
        ASSERT_NEAR(t_range.second, 180.0, TOLERANCE);
    });

    RUN_TEST("BBN redshift z ~ 4×10⁸", {
        double z_BBN = BigBangNucleosynthesis::redshift();
        ASSERT_NEAR(z_BBN, 4e8, 1e6);  // Approximate
    });

    RUN_TEST("Primordial hydrogen abundance ~75%", {
        auto abundances = BigBangNucleosynthesis::primordialAbundances();
        ASSERT_NEAR(abundances.at("H"), 0.75, 0.01);
    });

    RUN_TEST("Primordial helium-4 abundance ~25%", {
        auto abundances = BigBangNucleosynthesis::primordialAbundances();
        ASSERT_NEAR(abundances.at("He-4"), 0.25, 0.01);
    });

    RUN_TEST("Primordial deuterium abundance ~2.5×10⁻⁵", {
        auto abundances = BigBangNucleosynthesis::primordialAbundances();
        ASSERT_NEAR(abundances.at("D"), 2.5e-5, 1e-6);
    });

    RUN_TEST("H + He-4 = 100% (approximately)", {
        auto abundances = BigBangNucleosynthesis::primordialAbundances();
        double sum = abundances.at("H") + abundances.at("He-4");
        ASSERT_NEAR(sum, 1.0, 0.01);  // ~100% by mass
    });

    RUN_TEST("Deuterium binding energy 2.224 MeV", {
        double B_D = BigBangNucleosynthesis::deuteriumBindingEnergy();
        ASSERT_NEAR(B_D, 2.224, LOOSE_TOLERANCE);
    });

    RUN_TEST("Neutron-proton ratio n/p ~ 1/7", {
        double n_over_p = BigBangNucleosynthesis::neutronProtonRatio();
        ASSERT_NEAR(n_over_p, 1.0/7.0, LOOSE_TOLERANCE);
    });

    RUN_TEST("Helium-4 mass fraction Y_p ~ 0.25", {
        double Y_p = BigBangNucleosynthesis::helium4MassFraction();
        ASSERT_NEAR(Y_p, 0.25, 0.01);
    });

    RUN_TEST("Y_p from n/p ratio: Y_p = 2(n/p)/(1 + n/p)", {
        double n_over_p = BigBangNucleosynthesis::neutronProtonRatio();
        double Y_p_calc = 2.0 * n_over_p / (1.0 + n_over_p);
        double Y_p = BigBangNucleosynthesis::helium4MassFraction();
        ASSERT_NEAR(Y_p, Y_p_calc, TOLERANCE);
    });

    RUN_TEST("Deuterium abundance depends on η", {
        double D_H_6 = BigBangNucleosynthesis::deuteriumAbundance(6.0);
        double D_H_10 = BigBangNucleosynthesis::deuteriumAbundance(10.0);
        ASSERT_TRUE(D_H_10 < D_H_6);  // Higher η → less D (more burned)
    });

    RUN_TEST("Neutrino species from BBN N_ν ~ 3", {
        double N_nu = BigBangNucleosynthesis::neutrinoSpeciesConstraint();
        ASSERT_NEAR(N_nu, 2.99, 0.5);  // Within uncertainty
    });

    // ========================================
    // Baryogenesis Tests
    // ========================================

    RUN_TEST("Baryon-to-photon ratio η ~ 6×10⁻¹⁰", {
        double eta = Baryogenesis::baryonToPhotonRatio();
        ASSERT_NEAR(eta, 6.1e-10, 1e-11);
    });

    RUN_TEST("Baryon asymmetry η_B ~ 8.7×10⁻¹¹", {
        double eta_B = Baryogenesis::baryonAsymmetry();
        ASSERT_NEAR(eta_B, 8.7e-11, 1e-12);
    });

    RUN_TEST("Sakharov conditions exist", {
        auto conditions = Baryogenesis::sakharovConditions();
        ASSERT_TRUE(conditions.size() == 3);
    });

    RUN_TEST("Electroweak baryogenesis description exists", {
        std::string desc = Baryogenesis::electroweakBaryogenesis();
        ASSERT_TRUE(desc.length() > 0);
        ASSERT_TRUE(desc.find("100 GeV") != std::string::npos);
    });

    RUN_TEST("GUT baryogenesis description exists", {
        std::string desc = Baryogenesis::gutBaryogenesis();
        ASSERT_TRUE(desc.length() > 0);
        ASSERT_TRUE(desc.find("10¹⁵ GeV") != std::string::npos);
    });

    RUN_TEST("Leptogenesis description exists", {
        std::string desc = Baryogenesis::leptogenesis();
        ASSERT_TRUE(desc.length() > 0);
        ASSERT_TRUE(desc.find("Lepton") != std::string::npos);
    });

    RUN_TEST("Sphaleron processes description exists", {
        std::string desc = Baryogenesis::sphalerons();
        ASSERT_TRUE(desc.length() > 0);
        ASSERT_TRUE(desc.find("B+L") != std::string::npos);
    });

    RUN_TEST("Annihilation era description exists", {
        std::string desc = Baryogenesis::annihilationEra();
        ASSERT_TRUE(desc.length() > 0);
        ASSERT_TRUE(desc.find("10⁹") != std::string::npos);
    });

    // ========================================
    // Thermal History Tests
    // ========================================

    RUN_TEST("Planck epoch T ~ 10¹⁹ GeV, t ~ 10⁻⁴³ s", {
        auto epoch = ThermalHistory::planckEpoch();
        ASSERT_NEAR(epoch.first, 1e19, 1e18);  // GeV
        ASSERT_NEAR(epoch.second, 1e-43, 1e-44);  // seconds
    });

    RUN_TEST("GUT epoch T ~ 10¹⁵ GeV, t ~ 10⁻³⁷ s", {
        auto epoch = ThermalHistory::gutEpoch();
        ASSERT_NEAR(epoch.first, 1e15, 1e14);
        ASSERT_NEAR(epoch.second, 1e-37, 1e-38);
    });

    RUN_TEST("Electroweak epoch T ~ 100 GeV, t ~ 10⁻¹¹ s", {
        auto epoch = ThermalHistory::electroweakEpoch();
        ASSERT_NEAR(epoch.first, 100.0, 10.0);
        ASSERT_NEAR(epoch.second, 1e-11, 1e-12);
    });

    RUN_TEST("QCD epoch T ~ 200 MeV, t ~ 10⁻⁵ s", {
        auto epoch = ThermalHistory::qcdEpoch();
        ASSERT_NEAR(epoch.first, 0.2, 0.05);  // 200 MeV = 0.2 GeV
        ASSERT_NEAR(epoch.second, 1e-5, 1e-6);
    });

    RUN_TEST("Neutrino decoupling T ~ 1 MeV, t ~ 1 s", {
        auto epoch = ThermalHistory::neutrinoDecoupling();
        ASSERT_NEAR(epoch.first, 0.001, 0.0005);  // 1 MeV = 0.001 GeV
        ASSERT_NEAR(epoch.second, 1.0, 0.5);
    });

    RUN_TEST("BBN epoch T ~ 1 MeV, t ~ 180 s", {
        auto epoch = ThermalHistory::bbnEpoch();
        ASSERT_NEAR(epoch.first, 0.001, 0.0005);
        ASSERT_NEAR(epoch.second, 180.0, 10.0);
    });

    RUN_TEST("Matter-radiation equality in thermal history", {
        auto epoch = ThermalHistory::matterRadiationEquality();
        ASSERT_TRUE(epoch.first > 0.0);  // Temperature (GeV)
        ASSERT_TRUE(epoch.second > 1e12);  // ~47,000 years in seconds
    });

    RUN_TEST("Recombination epoch t ~ 380,000 years", {
        auto epoch = ThermalHistory::recombination();
        double t_years = epoch.second / (365.25 * 24 * 3600);
        ASSERT_NEAR(t_years, 380000.0, 10000.0);
    });

    RUN_TEST("Dark ages duration", {
        auto epoch = ThermalHistory::darkAges();
        double duration = epoch.second - epoch.first;  // years
        ASSERT_TRUE(duration > 1e7);  // > 10 million years
    });

    RUN_TEST("Reionization era t ~ 150 Myr - 1 Gyr", {
        auto epoch = ThermalHistory::reionization();
        ASSERT_NEAR(epoch.first, 150e6, 50e6);  // years
        ASSERT_NEAR(epoch.second, 1e9, 0.5e9);
    });

    RUN_TEST("Chronological order: Planck < GUT < EW < QCD", {
        auto planck = ThermalHistory::planckEpoch();
        auto gut = ThermalHistory::gutEpoch();
        auto ew = ThermalHistory::electroweakEpoch();
        auto qcd = ThermalHistory::qcdEpoch();
        ASSERT_TRUE(planck.second < gut.second);
        ASSERT_TRUE(gut.second < ew.second);
        ASSERT_TRUE(ew.second < qcd.second);
    });

    RUN_TEST("Chronological order: BBN < equality < recombination", {
        auto bbn = ThermalHistory::bbnEpoch();
        auto eq = ThermalHistory::matterRadiationEquality();
        auto rec = ThermalHistory::recombination();
        ASSERT_TRUE(bbn.second < eq.second);
        ASSERT_TRUE(eq.second < rec.second);
    });

    // ========================================
    // Cross-module consistency tests
    // ========================================

    RUN_TEST("CMB decoupling consistent with recombination", {
        double t_dec_CMB = CosmicMicrowaveBackground::decouplingTime();
        auto rec_epoch = ThermalHistory::recombination();
        double t_rec = rec_epoch.second / (365.25 * 24 * 3600);  // Convert to years
        ASSERT_NEAR(t_dec_CMB, t_rec, 10000.0);  // Should match
    });

    RUN_TEST("Radiation-matter equality consistent across modules", {
        double z_eq_rad = RadiationEra::equalityRedshift();
        double z_eq_mat = MatterEra::startRedshift();
        ASSERT_NEAR(z_eq_rad, z_eq_mat, TOLERANCE);
    });

    RUN_TEST("BBN time consistent across modules", {
        auto bbn_range = BigBangNucleosynthesis::timeRange();
        auto thermal_bbn = ThermalHistory::bbnEpoch();
        ASSERT_NEAR(bbn_range.second, thermal_bbn.second, 10.0);
    });

    // ========================================
    // Physical law tests
    // ========================================

    RUN_TEST("Matter era growth faster than radiation era", {
        double t_ratio = 4.0;
        double a_rad = RadiationEra::scaleFactorEvolution(t_ratio);  // t^(1/2)
        double a_mat = MatterEra::scaleFactorEvolution(t_ratio);     // t^(2/3)
        ASSERT_TRUE(a_mat > a_rad);  // Matter era expands faster
    });

    RUN_TEST("Hubble parameter decreases with time (radiation era)", {
        double H_1 = RadiationEra::hubbleParameter(1.0);
        double H_100 = RadiationEra::hubbleParameter(100.0);
        ASSERT_TRUE(H_1 > H_100);
    });

    RUN_TEST("Hubble parameter decreases with time (matter era)", {
        double H_1 = MatterEra::hubbleParameter(1.0);
        double H_100 = MatterEra::hubbleParameter(100.0);
        ASSERT_TRUE(H_1 > H_100);
    });

    RUN_TEST("Temperature decreases as universe expands", {
        double T_z0 = CosmicMicrowaveBackground::temperatureAtRedshift(0.0);
        double T_z10 = CosmicMicrowaveBackground::temperatureAtRedshift(10.0);
        ASSERT_TRUE(T_z10 > T_z0);  // Higher z = earlier = hotter
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4: Early Universe Cosmology" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "✓ All early universe tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - CMB temperature T₀ = 2.7255 K" << std::endl;
        std::cout << "  - CMB decoupling z ~ 1100, t ~ 380,000 years" << std::endl;
        std::cout << "  - Radiation era: a ∝ t^(1/2), H = 1/(2t)" << std::endl;
        std::cout << "  - Matter era: a ∝ t^(2/3), H = 2/(3t)" << std::endl;
        std::cout << "  - Matter-radiation equality z ~ 3400" << std::endl;
        std::cout << "  - BBN abundances: H (75%), He-4 (25%)" << std::endl;
        std::cout << "  - Baryon asymmetry η ~ 6×10⁻¹⁰" << std::endl;
        std::cout << "  - Thermal history from Planck to reionization" << std::endl;
        return 0;
    } else {
        return 1;
    }
}
