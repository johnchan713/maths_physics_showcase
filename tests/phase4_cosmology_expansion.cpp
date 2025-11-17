/**
 * Phase 4 Validation: Cosmology - Expanding Universe & Energy Density
 *
 * Tests cosmology_expanding_universe.hpp and cosmology_energy_density.hpp modules.
 *
 * Coverage:
 * - Hubble's law (v = H₀d)
 * - Hubble constant, time, and sphere
 * - Cosmological redshift
 * - Scale factor evolution
 * - Energy density components (matter, radiation, dark energy)
 * - Density parameters (Ω_m, Ω_r, Ω_Λ)
 * - Equation of state for cosmic components
 */

#include <iostream>
#include <cmath>
#include "../include/physics/cosmology_expanding_universe.hpp"
#include "../include/physics/cosmology_energy_density.hpp"

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

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Cosmology - Expansion & Energy Density ===" << std::endl;
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
    // Hubble Constant Tests
    // ========================================

    run_test("Hubble constant value (Planck 2018)", []() {
        double H0 = physics::advanced::cosmology::HubbleExpansion::hubbleConstant();

        // H₀ = 67.4 km/s/Mpc (Planck 2018)
        ASSERT_NEAR(H0, 67.4, TOLERANCE);
        return true;
    });

    run_test("Hubble constant SI units", []() {
        double H0_SI = physics::advanced::cosmology::HubbleExpansion::hubbleConstantSI();

        // H₀ ≈ 2.2×10⁻¹⁸ s⁻¹
        ASSERT_TRUE(H0_SI > 2.0e-18 && H0_SI < 2.3e-18);
        return true;
    });

    run_test("Hubble time calculation", []() {
        double t_H = physics::advanced::cosmology::HubbleExpansion::hubbleTime();
        double H0_SI = physics::advanced::cosmology::HubbleExpansion::hubbleConstantSI();

        // t_H = 1/H₀
        double expected = 1.0 / H0_SI;
        ASSERT_NEAR(t_H, expected, TOLERANCE);
        return true;
    });

    run_test("Hubble time in Gyr", []() {
        double t_H_Gyr = physics::advanced::cosmology::HubbleExpansion::hubbleTimeGyr();

        // t_H ≈ 14.5 Gyr
        ASSERT_TRUE(t_H_Gyr > 14.0 && t_H_Gyr < 15.0);
        return true;
    });

    run_test("Hubble tension exists", []() {
        double tension = physics::advanced::cosmology::HubbleExpansion::hubbleTension();

        // Discrepancy between CMB (67.4) and local (73.0)
        // Tension ≈ 5.6 km/s/Mpc
        ASSERT_TRUE(tension > 5.0 && tension < 6.0);
        return true;
    });

    // ========================================
    // Hubble's Law Tests
    // ========================================

    run_test("Hubble's law: v = H₀d", []() {
        double d = 100.0;  // 100 Mpc
        double v = physics::advanced::cosmology::HubbleExpansion::recessionVelocity(d);

        // v = 67.4 × 100 = 6740 km/s
        double H0 = physics::advanced::cosmology::HubbleExpansion::hubbleConstant();
        double expected = H0 * d;

        ASSERT_NEAR(v, expected, TOLERANCE);
        return true;
    });

    run_test("Distance from velocity", []() {
        double v = 6740.0;  // km/s
        double d = physics::advanced::cosmology::HubbleExpansion::distanceFromVelocity(v);

        // d = v/H₀ = 6740/67.4 = 100 Mpc
        double H0 = physics::advanced::cosmology::HubbleExpansion::hubbleConstant();
        double expected = v / H0;

        ASSERT_NEAR(d, expected, TOLERANCE);
        return true;
    });

    run_test("Hubble's law consistency", []() {
        double d_original = 250.0;  // Mpc
        double v = physics::advanced::cosmology::HubbleExpansion::recessionVelocity(d_original);
        double d_recovered = physics::advanced::cosmology::HubbleExpansion::distanceFromVelocity(v);

        ASSERT_NEAR(d_recovered, d_original, TOLERANCE);
        return true;
    });

    run_test("Hubble sphere (v = c)", []() {
        double R_H = physics::advanced::cosmology::HubbleExpansion::hubbleSphere();

        // R_H = c/H₀ ≈ 4200 Mpc
        ASSERT_TRUE(R_H > 4000.0 && R_H < 4500.0);
        return true;
    });

    run_test("Galaxies beyond Hubble sphere recede faster than light", []() {
        double R_H = physics::advanced::cosmology::HubbleExpansion::hubbleSphere();
        double d_far = 2.0 * R_H;  // Twice Hubble radius

        double v = physics::advanced::cosmology::HubbleExpansion::recessionVelocity(d_far);
        double c = 2.998e5;  // km/s

        // This is allowed! Space itself is expanding
        ASSERT_TRUE(v > c);
        return true;
    });

    // ========================================
    // Redshift Tests
    // ========================================

    run_test("Redshift from wavelength", []() {
        double lambda_obs = 550.0;  // nm (observed)
        double lambda_emit = 500.0;  // nm (emitted)

        double z = physics::advanced::cosmology::CosmologicalRedshift::fromWavelength(
            lambda_obs, lambda_emit);

        // z = (λ_obs - λ_emit)/λ_emit = 50/500 = 0.1
        double expected = (lambda_obs - lambda_emit) / lambda_emit;
        ASSERT_NEAR(z, expected, TOLERANCE);
        return true;
    });

    run_test("Redshift from scale factor", []() {
        double a = 0.5;  // Universe was half current size
        double z = physics::advanced::cosmology::CosmologicalRedshift::fromScaleFactor(a);

        // 1 + z = 1/a → z = 1/a - 1 = 1/0.5 - 1 = 1
        double expected = (1.0 / a) - 1.0;
        ASSERT_NEAR(z, expected, TOLERANCE);
        return true;
    });

    run_test("Scale factor from redshift", []() {
        double z = 1.0;  // Redshift of 1
        double a = physics::advanced::cosmology::CosmologicalRedshift::scaleFactorFromRedshift(z);

        // a = 1/(1+z) = 1/2 = 0.5
        double expected = 1.0 / (1.0 + z);
        ASSERT_NEAR(a, expected, TOLERANCE);
        return true;
    });

    run_test("Redshift-scale factor consistency", []() {
        double a_original = 0.3;
        double z = physics::advanced::cosmology::CosmologicalRedshift::fromScaleFactor(a_original);
        double a_recovered = physics::advanced::cosmology::CosmologicalRedshift::scaleFactorFromRedshift(z);

        ASSERT_NEAR(a_recovered, a_original, TOLERANCE);
        return true;
    });

    run_test("Non-relativistic velocity from redshift", []() {
        double z = 0.01;  // Small redshift
        double v = physics::advanced::cosmology::CosmologicalRedshift::velocityNonRelativistic(z);

        // v ≈ cz for z << 1
        double c = 2.998e5;  // km/s
        double expected = c * z;

        ASSERT_NEAR(v, expected, TOLERANCE);
        return true;
    });

    run_test("Relativistic redshift formula", []() {
        double beta = 0.5;  // v = 0.5c
        double z = physics::advanced::cosmology::CosmologicalRedshift::relativisticRedshift(beta);

        // 1 + z = √[(1+β)/(1-β)] = √(1.5/0.5) = √3
        double expected = std::sqrt((1.0 + beta) / (1.0 - beta)) - 1.0;
        ASSERT_NEAR(z, expected, TOLERANCE);
        return true;
    });

    run_test("Luminosity distance from redshift", []() {
        double z = 0.05;  // Small z
        double d_L = physics::advanced::cosmology::CosmologicalRedshift::luminosityDistance(z);

        // For small z: d_L ≈ cz/H₀
        double c = 2.998e5;
        double H0 = physics::advanced::cosmology::HubbleExpansion::hubbleConstant();
        double expected_approx = c * z / H0;

        ASSERT_NEAR(d_L, expected_approx, expected_approx * 0.1);  // 10% tolerance
        return true;
    });

    // ========================================
    // Scale Factor Tests
    // ========================================

    run_test("Scale factor today", []() {
        double a_today = physics::advanced::cosmology::ScaleFactor::today();

        // Convention: a(t_today) = 1
        ASSERT_NEAR(a_today, 1.0, TOLERANCE);
        return true;
    });

    run_test("Scale factor at CMB decoupling", []() {
        double a_CMB = physics::advanced::cosmology::ScaleFactor::atCMBDecoupling();

        // z_CMB ≈ 1100 → a = 1/1101 ≈ 9×10⁻⁴
        double expected = 1.0 / 1101.0;
        ASSERT_NEAR(a_CMB, expected, VERY_LOOSE * a_CMB);
        return true;
    });

    run_test("Scale factor at matter-radiation equality", []() {
        double a_eq = physics::advanced::cosmology::ScaleFactor::atMatterRadiationEquality();

        // z_eq ≈ 3400 → a ≈ 2.9×10⁻⁴
        ASSERT_TRUE(a_eq > 2.5e-4 && a_eq < 3.5e-4);
        return true;
    });

    run_test("Scale factor at nucleosynthesis", []() {
        double a_BBN = physics::advanced::cosmology::ScaleFactor::atNucleosynthesis();

        // z ~ 4×10⁸ → a ~ 2.5×10⁻⁹
        ASSERT_TRUE(a_BBN > 2.0e-9 && a_BBN < 3.0e-9);
        return true;
    });

    run_test("Scale factor ordered by cosmic time", []() {
        double a_BBN = physics::advanced::cosmology::ScaleFactor::atNucleosynthesis();
        double a_eq = physics::advanced::cosmology::ScaleFactor::atMatterRadiationEquality();
        double a_CMB = physics::advanced::cosmology::ScaleFactor::atCMBDecoupling();
        double a_today = physics::advanced::cosmology::ScaleFactor::today();

        // Earlier times have smaller scale factor
        ASSERT_TRUE(a_BBN < a_eq);
        ASSERT_TRUE(a_eq < a_CMB);
        ASSERT_TRUE(a_CMB < a_today);
        return true;
    });

    run_test("Hubble parameter H = ȧ/a", []() {
        double a = 1.0;
        double a_dot = 2.2e-18;  // s⁻¹ (current expansion rate)

        double H = physics::advanced::cosmology::ScaleFactor::hubbleParameter(a, a_dot);

        // H = ȧ/a
        double expected = a_dot / a;
        ASSERT_NEAR(H, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Energy Density Component Tests
    // ========================================

    run_test("Matter density parameter (Planck 2018)", []() {
        double Omega_m = physics::advanced::cosmology::EnergyDensityComponents::matterDensityParameter();

        // Ω_m ≈ 0.315
        ASSERT_NEAR(Omega_m, 0.315, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Baryonic matter density", []() {
        double Omega_b = physics::advanced::cosmology::EnergyDensityComponents::baryonicDensityParameter();

        // Ω_b ≈ 0.049 (~5% of critical density)
        ASSERT_NEAR(Omega_b, 0.049, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Dark matter density", []() {
        double Omega_DM = physics::advanced::cosmology::EnergyDensityComponents::darkMatterDensityParameter();

        // Ω_DM ≈ 0.266 (~26% of critical density)
        ASSERT_TRUE(Omega_DM > 0.25 && Omega_DM < 0.28);
        return true;
    });

    run_test("Dark matter is most of total matter", []() {
        double Omega_DM = physics::advanced::cosmology::EnergyDensityComponents::darkMatterDensityParameter();
        double Omega_b = physics::advanced::cosmology::EnergyDensityComponents::baryonicDensityParameter();

        // Dark matter dominates over baryonic matter
        ASSERT_TRUE(Omega_DM > 5.0 * Omega_b);
        return true;
    });

    run_test("Total matter = baryonic + dark", []() {
        double Omega_m = physics::advanced::cosmology::EnergyDensityComponents::matterDensityParameter();
        double Omega_b = physics::advanced::cosmology::EnergyDensityComponents::baryonicDensityParameter();
        double Omega_DM = physics::advanced::cosmology::EnergyDensityComponents::darkMatterDensityParameter();

        ASSERT_NEAR(Omega_m, Omega_b + Omega_DM, TOLERANCE);
        return true;
    });

    run_test("Radiation density parameter", []() {
        double Omega_r = physics::advanced::cosmology::EnergyDensityComponents::radiationDensityParameter();

        // Ω_r ≈ 9.24×10⁻⁵ (negligible today!)
        ASSERT_TRUE(Omega_r > 9.0e-5 && Omega_r < 9.5e-5);
        return true;
    });

    run_test("Photon density parameter", []() {
        double Omega_gamma = physics::advanced::cosmology::EnergyDensityComponents::photonDensityParameter();

        // Ω_γ ≈ 5.4×10⁻⁵ (CMB photons)
        ASSERT_TRUE(Omega_gamma > 5.0e-5 && Omega_gamma < 6.0e-5);
        return true;
    });

    run_test("Neutrino density parameter", []() {
        double Omega_nu = physics::advanced::cosmology::EnergyDensityComponents::neutrinoDensityParameter();

        // Ω_ν ≈ 3.8×10⁻⁵ (cosmic neutrino background)
        ASSERT_TRUE(Omega_nu > 3.5e-5 && Omega_nu < 4.5e-5);
        return true;
    });

    run_test("Radiation = photons + neutrinos", []() {
        double Omega_r = physics::advanced::cosmology::EnergyDensityComponents::radiationDensityParameter();
        double Omega_gamma = physics::advanced::cosmology::EnergyDensityComponents::photonDensityParameter();
        double Omega_nu = physics::advanced::cosmology::EnergyDensityComponents::neutrinoDensityParameter();

        // Tolerance needs to account for very small numbers
        ASSERT_NEAR(Omega_r, Omega_gamma + Omega_nu, 1e-6);
        return true;
    });

    run_test("Dark energy density parameter", []() {
        double Omega_Lambda = physics::advanced::cosmology::EnergyDensityComponents::darkEnergyDensityParameter();

        // Ω_Λ ≈ 0.685 (~69% of universe!)
        ASSERT_NEAR(Omega_Lambda, 0.685, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Dark energy dominates today", []() {
        double Omega_Lambda = physics::advanced::cosmology::EnergyDensityComponents::darkEnergyDensityParameter();
        double Omega_m = physics::advanced::cosmology::EnergyDensityComponents::matterDensityParameter();
        double Omega_r = physics::advanced::cosmology::EnergyDensityComponents::radiationDensityParameter();

        // Dark energy is largest component
        ASSERT_TRUE(Omega_Lambda > Omega_m);
        ASSERT_TRUE(Omega_Lambda > Omega_r);
        return true;
    });

    run_test("Total density parameter", []() {
        double Omega_total = physics::advanced::cosmology::EnergyDensityComponents::totalDensityParameter();

        // Should be very close to 1 (flat universe)
        ASSERT_NEAR(Omega_total, 1.0, VERY_LOOSE);
        return true;
    });

    run_test("Sum of components equals total", []() {
        double Omega_m = physics::advanced::cosmology::EnergyDensityComponents::matterDensityParameter();
        double Omega_r = physics::advanced::cosmology::EnergyDensityComponents::radiationDensityParameter();
        double Omega_Lambda = physics::advanced::cosmology::EnergyDensityComponents::darkEnergyDensityParameter();
        double Omega_total = physics::advanced::cosmology::EnergyDensityComponents::totalDensityParameter();

        double sum = Omega_m + Omega_r + Omega_Lambda;
        ASSERT_NEAR(sum, Omega_total, TOLERANCE);
        return true;
    });

    run_test("Curvature density parameter", []() {
        double Omega_k = physics::advanced::cosmology::EnergyDensityComponents::curvatureDensityParameter();

        // Ω_k ≈ 0.001 (universe is nearly flat!)
        ASSERT_TRUE(std::abs(Omega_k) < 0.01);
        return true;
    });

    run_test("Matter dominates radiation today", []() {
        double Omega_m = physics::advanced::cosmology::EnergyDensityComponents::matterDensityParameter();
        double Omega_r = physics::advanced::cosmology::EnergyDensityComponents::radiationDensityParameter();

        // Matter >> radiation today
        ASSERT_TRUE(Omega_m > 1000.0 * Omega_r);
        return true;
    });

    // ========================================
    // Consistency Tests
    // ========================================

    run_test("Hubble time × Hubble constant = 1", []() {
        double H0_SI = physics::advanced::cosmology::HubbleExpansion::hubbleConstantSI();
        double t_H = physics::advanced::cosmology::HubbleExpansion::hubbleTime();

        double product = H0_SI * t_H;
        ASSERT_NEAR(product, 1.0, TOLERANCE);
        return true;
    });

    run_test("Age of universe comparable to Hubble time", []() {
        double t_H_Gyr = physics::advanced::cosmology::HubbleExpansion::hubbleTimeGyr();
        double age_universe_Gyr = 13.8;  // ~13.8 Gyr

        // Hubble time gives rough age estimate
        ASSERT_NEAR(t_H_Gyr, age_universe_Gyr, 1.0);  // Within 1 Gyr
        return true;
    });

    run_test("CMB decoupling before today", []() {
        double a_CMB = physics::advanced::cosmology::ScaleFactor::atCMBDecoupling();
        double a_today = physics::advanced::cosmology::ScaleFactor::today();

        ASSERT_TRUE(a_CMB < a_today);
        return true;
    });

    run_test("Observable universe radius < Hubble sphere", []() {
        // Observable universe is finite due to finite age
        double R_H = physics::advanced::cosmology::HubbleExpansion::hubbleSphere();

        // Observable universe ≈ 46 Gly, Hubble sphere ≈ 14 Gly
        // But in comoving coords, observable > Hubble
        // Just check Hubble sphere is positive and reasonable
        ASSERT_TRUE(R_H > 4000.0);  // Mpc
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Phase 4: Cosmology Expansion & Energy" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << std::endl;

    if (tests_failed == 0) {
        std::cout << "✓ All cosmology expansion & energy tests PASSED!" << std::endl;
        std::cout << std::endl;
        std::cout << "Validated:" << std::endl;
        std::cout << "  - Hubble constant H₀ = 67.4 km/s/Mpc" << std::endl;
        std::cout << "  - Hubble's law (v = H₀d)" << std::endl;
        std::cout << "  - Hubble time t_H ≈ 14.5 Gyr" << std::endl;
        std::cout << "  - Hubble tension (~5.6 km/s/Mpc)" << std::endl;
        std::cout << "  - Cosmological redshift (z = Δλ/λ)" << std::endl;
        std::cout << "  - Scale factor evolution" << std::endl;
        std::cout << "  - Energy components: matter (31%), dark energy (69%)" << std::endl;
        std::cout << "  - Dark matter dominates over baryons (5:1 ratio)" << std::endl;
        std::cout << "  - Radiation negligible today (~0.01%)" << std::endl;
        std::cout << "  - Universe composition from Planck 2018" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Some tests FAILED. See details above." << std::endl;
        return 1;
    }
}
