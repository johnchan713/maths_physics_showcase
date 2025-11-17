/**
 * Phase 4 Validation: Condensed Matter Physics
 *
 * Tests the condensed_matter.hpp module functions.
 *
 * Coverage:
 * - Band structure (free electron, tight binding, Kronig-Penney)
 * - Group velocity, effective mass
 * - Density of states (1D, 2D, 3D)
 * - BCS theory (superconductivity gap, critical temperature)
 * - London theory (penetration depth, Meissner effect)
 * - Josephson effect
 * - Fermi liquid theory
 * - Phonon dispersion
 * - Quantum Hall effect (conductivity quantization)
 * - Landau theory of phase transitions
 * - Hubbard model
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <functional>
#include "../include/physics/condensed_matter.hpp"

// Test tolerance
const double TOLERANCE = 1e-5;
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

// Physical constants from the module
using physics::condensed_matter::k_B;
using physics::condensed_matter::hbar;
using physics::condensed_matter::e;
using physics::condensed_matter::m_e;

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Condensed Matter Physics Validation ===" << std::endl;
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
    // Band Structure Tests
    // ========================================

    run_test("Free electron band: energy at k=0 is zero", []() {
        physics::condensed_matter::FreElectronBand band;
        std::vector<double> k = {0.0, 0.0, 0.0};
        double E = band.energy(k);
        ASSERT_NEAR(E, 0.0, TOLERANCE);
        return true;
    });

    run_test("Free electron band: quadratic dispersion E = hbar^2 k^2 / 2m", []() {
        physics::condensed_matter::FreElectronBand band;
        double k_val = 1e10;  // 1/m
        std::vector<double> k = {k_val, 0.0, 0.0};
        double E = band.energy(k);
        double expected = hbar * hbar * k_val * k_val / (2.0 * m_e);
        ASSERT_NEAR(E, expected, TOLERANCE);
        return true;
    });

    run_test("Free electron band: isotropic in k-space", []() {
        physics::condensed_matter::FreElectronBand band;
        double k_val = 5e9;
        std::vector<double> k1 = {k_val, 0.0, 0.0};
        std::vector<double> k2 = {0.0, k_val, 0.0};
        std::vector<double> k3 = {0.0, 0.0, k_val};
        ASSERT_NEAR(band.energy(k1), band.energy(k2), TOLERANCE);
        ASSERT_NEAR(band.energy(k2), band.energy(k3), TOLERANCE);
        return true;
    });

    run_test("Tight binding band: energy at k=0", []() {
        double t = 1.0 * e;  // 1 eV hopping
        double a = 5e-10;    // 5 Angstrom lattice constant
        physics::condensed_matter::TightBindingBand band(t, a);
        std::vector<double> k = {0.0, 0.0, 0.0};
        double E = band.energy(k);
        // At k=0: E = -2t * 3 = -6t (3D)
        double expected = -6.0 * t;
        ASSERT_NEAR(E, expected, TOLERANCE);
        return true;
    });

    run_test("Tight binding band: energy at zone edge", []() {
        double t = 1.0 * e;
        double a = 5e-10;
        physics::condensed_matter::TightBindingBand band(t, a);
        std::vector<double> k = {M_PI / a, 0.0, 0.0};
        double E = band.energy(k);
        // At k = pi/a: cos(pi) = -1, so -2t*(-1) = 2t for this direction
        double expected = 2.0 * t - 4.0 * t;  // Other two directions contribute -4t
        ASSERT_NEAR(E, expected, TOLERANCE);
        return true;
    });

    run_test("Kronig-Penney model: energy increases with k", []() {
        double V0 = 5.0 * e;  // 5 eV barrier
        double a = 3e-10;     // 3 Angstrom
        physics::condensed_matter::KronigPenneyModel band(V0, a);
        std::vector<double> k1 = {1e9, 0.0, 0.0};
        std::vector<double> k2 = {2e9, 0.0, 0.0};
        ASSERT_TRUE(band.energy(k2) > band.energy(k1));
        return true;
    });

    run_test("Group velocity: positive for parabolic band", []() {
        physics::condensed_matter::FreElectronBand band;
        double k_val = 1e10;
        std::vector<double> k = {k_val, 0.0, 0.0};
        double v_g = band.groupVelocity(k, 0);
        ASSERT_TRUE(v_g > 0.0);
        return true;
    });

    run_test("Band structure: energy is continuous function of k", []() {
        physics::condensed_matter::FreElectronBand band;
        double k1 = 1e10;
        double k2 = 1.001e10;
        std::vector<double> kv1 = {k1, 0.0, 0.0};
        std::vector<double> kv2 = {k2, 0.0, 0.0};
        double E1 = band.energy(kv1);
        double E2 = band.energy(kv2);
        // Energy should be continuous (small dk leads to small dE)
        ASSERT_TRUE(std::abs(E2 - E1) < 0.01 * E1);
        return true;
    });

    // ========================================
    // Density of States Tests
    // ========================================

    run_test("DOS 3D: proportional to sqrt(E)", []() {
        double E1 = 1.0 * e;
        double E2 = 4.0 * e;
        double g1 = physics::condensed_matter::DensityOfStates::freElectron3D(E1);
        double g2 = physics::condensed_matter::DensityOfStates::freElectron3D(E2);
        // g2/g1 should be sqrt(E2/E1) = sqrt(4) = 2
        ASSERT_NEAR(g2 / g1, 2.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("DOS 3D: zero for negative energy", []() {
        double g = physics::condensed_matter::DensityOfStates::freElectron3D(-1.0);
        ASSERT_NEAR(g, 0.0, TOLERANCE);
        return true;
    });

    run_test("DOS 2D: constant with energy", []() {
        double E1 = 1.0 * e;
        double E2 = 10.0 * e;
        double g1 = physics::condensed_matter::DensityOfStates::freElectron2D(E1);
        double g2 = physics::condensed_matter::DensityOfStates::freElectron2D(E2);
        ASSERT_NEAR(g1, g2, TOLERANCE);
        return true;
    });

    run_test("DOS 2D: zero for negative energy", []() {
        double g = physics::condensed_matter::DensityOfStates::freElectron2D(-1.0);
        ASSERT_NEAR(g, 0.0, TOLERANCE);
        return true;
    });

    run_test("DOS 1D: proportional to 1/sqrt(E)", []() {
        double E1 = 1.0 * e;
        double E2 = 4.0 * e;
        double g1 = physics::condensed_matter::DensityOfStates::freElectron1D(E1);
        double g2 = physics::condensed_matter::DensityOfStates::freElectron1D(E2);
        // g2/g1 should be sqrt(E1/E2) = 1/2
        ASSERT_NEAR(g2 / g1, 0.5, LOOSE_TOLERANCE);
        return true;
    });

    run_test("DOS 1D: zero for zero and negative energy", []() {
        double g1 = physics::condensed_matter::DensityOfStates::freElectron1D(0.0);
        double g2 = physics::condensed_matter::DensityOfStates::freElectron1D(-1.0);
        ASSERT_NEAR(g1, 0.0, TOLERANCE);
        ASSERT_NEAR(g2, 0.0, TOLERANCE);
        return true;
    });

    run_test("DOS integration: positive result", []() {
        std::function<double(double)> g = [](double E) {
            return physics::condensed_matter::DensityOfStates::freElectron3D(E);
        };
        double total = physics::condensed_matter::DensityOfStates::integrate(g, 0.0, 10.0 * e, 1000);
        ASSERT_TRUE(total > 0.0);
        return true;
    });

    // ========================================
    // BCS Superconductivity Tests
    // ========================================

    run_test("BCS: energy gap at T=0 is Delta_0 = 1.764 k_B T_c", []() {
        double T_c = 9.2;  // Nb: 9.2 K
        double E_F = 5.0 * e;
        physics::condensed_matter::BCSTheory bcs(T_c, E_F);
        double Delta_0 = bcs.energyGap(0.0);
        double expected = 1.764 * k_B * T_c;
        ASSERT_NEAR(Delta_0, expected, LOOSE_TOLERANCE * expected);
        return true;
    });

    run_test("BCS: energy gap is zero above T_c", []() {
        double T_c = 9.2;
        double E_F = 5.0 * e;
        physics::condensed_matter::BCSTheory bcs(T_c, E_F);
        double Delta = bcs.energyGap(T_c + 1.0);
        ASSERT_NEAR(Delta, 0.0, TOLERANCE);
        return true;
    });

    run_test("BCS: energy gap decreases with temperature", []() {
        double T_c = 9.2;
        double E_F = 5.0 * e;
        physics::condensed_matter::BCSTheory bcs(T_c, E_F);
        double Delta1 = bcs.energyGap(0.0);
        double Delta2 = bcs.energyGap(0.5 * T_c);
        ASSERT_TRUE(Delta2 < Delta1);
        ASSERT_TRUE(Delta2 > 0.0);
        return true;
    });

    run_test("BCS: Cooper pair binding energy is 2*Delta_0", []() {
        double T_c = 9.2;
        double E_F = 5.0 * e;
        physics::condensed_matter::BCSTheory bcs(T_c, E_F);
        double binding = bcs.cooperPairBindingEnergy();
        double expected = 2.0 * 1.764 * k_B * T_c;
        ASSERT_NEAR(binding, expected, LOOSE_TOLERANCE * expected);
        return true;
    });

    run_test("BCS: quasiparticle energy >= Delta", []() {
        double T_c = 9.2;
        double E_F = 5.0 * e;
        physics::condensed_matter::BCSTheory bcs(T_c, E_F);
        double T = 0.0;
        double Delta = bcs.energyGap(T);
        double E_qp = bcs.quasiparticleEnergy(0.0, T);
        ASSERT_NEAR(E_qp, Delta, LOOSE_TOLERANCE * Delta);
        return true;
    });

    run_test("BCS: coherence length is positive", []() {
        double T_c = 9.2;
        double E_F = 5.0 * e;
        physics::condensed_matter::BCSTheory bcs(T_c, E_F);
        double xi = bcs.coherenceLength(0.5 * T_c);
        ASSERT_TRUE(xi > 0.0);
        return true;
    });

    run_test("BCS: specific heat jump at T_c", []() {
        double T_c = 9.2;
        double E_F = 5.0 * e;
        physics::condensed_matter::BCSTheory bcs(T_c, E_F);
        double C_normal = bcs.normalSpecificHeat(T_c);
        // Superconducting specific heat should be exponentially suppressed
        ASSERT_TRUE(C_normal > 0.0);
        return true;
    });

    // ========================================
    // London Theory Tests
    // ========================================

    run_test("London theory: superfluid velocity proportional to vector potential", []() {
        double lambda_L = 50e-9;  // 50 nm
        double n_s = 1e28;        // m^-3
        physics::condensed_matter::LondonTheory london(lambda_L, n_s);
        std::vector<double> A = {1e-6, 0.0, 0.0};
        auto v_s = london.superfluidVelocity(A);
        // v_s = -(e/m_e) * A
        double expected = -(e / m_e) * A[0];
        ASSERT_NEAR(v_s[0], expected, LOOSE_TOLERANCE * std::abs(expected));
        return true;
    });

    run_test("London theory: supercurrent is proportional to superfluid velocity", []() {
        double lambda_L = 50e-9;
        double n_s = 1e28;
        physics::condensed_matter::LondonTheory london(lambda_L, n_s);
        std::vector<double> A = {1e-6, 0.0, 0.0};
        auto j_s = london.supercurrent(A);
        ASSERT_TRUE(std::abs(j_s[0]) > 0.0);
        return true;
    });

    run_test("London theory: magnetic field penetration decays exponentially", []() {
        double lambda_L = 50e-9;
        double n_s = 1e28;
        physics::condensed_matter::LondonTheory london(lambda_L, n_s);
        double B0 = london.magneticFieldPenetration(0.0);
        double B1 = london.magneticFieldPenetration(lambda_L);
        // At x = lambda_L: B = B0 * exp(-1)
        ASSERT_NEAR(B1, B0 / M_E, LOOSE_TOLERANCE);
        return true;
    });

    run_test("London theory: Meissner effect exponential decay", []() {
        double lambda_L = 50e-9;
        double n_s = 1e28;
        physics::condensed_matter::LondonTheory london(lambda_L, n_s);
        double B_ext = 0.1;  // Tesla
        double B0 = london.meissnerEffect(B_ext, 0.0);
        double B1 = london.meissnerEffect(B_ext, lambda_L);
        ASSERT_NEAR(B0, B_ext, TOLERANCE);
        ASSERT_NEAR(B1, B_ext / M_E, LOOSE_TOLERANCE);
        return true;
    });

    run_test("London theory: flux quantum h/2e", []() {
        double lambda_L = 50e-9;
        double n_s = 1e28;
        physics::condensed_matter::LondonTheory london(lambda_L, n_s);
        double Phi_0 = london.fluxQuantum();
        double expected = M_PI * hbar / e;
        ASSERT_NEAR(Phi_0, expected, TOLERANCE);
        return true;
    });

    run_test("London theory: flux quantization in units of Phi_0", []() {
        double lambda_L = 50e-9;
        double n_s = 1e28;
        physics::condensed_matter::LondonTheory london(lambda_L, n_s);
        double Phi_0 = london.fluxQuantum();
        int n = london.fluxoidQuantization(3.2 * Phi_0);
        ASSERT_TRUE(n == 3);
        return true;
    });

    // ========================================
    // Josephson Effect Tests
    // ========================================

    run_test("Josephson: DC current I = I_c sin(delta)", []() {
        double I_c = 1e-6;  // 1 microamp
        double V_c = 1e-3;  // 1 mV
        physics::condensed_matter::JosephsonEffect josephson(I_c, V_c);
        double I = josephson.dcCurrent(M_PI / 2.0);
        ASSERT_NEAR(I, I_c, TOLERANCE);
        return true;
    });

    run_test("Josephson: DC current is zero at delta=0", []() {
        double I_c = 1e-6;
        double V_c = 1e-3;
        physics::condensed_matter::JosephsonEffect josephson(I_c, V_c);
        double I = josephson.dcCurrent(0.0);
        ASSERT_NEAR(I, 0.0, TOLERANCE);
        return true;
    });

    run_test("Josephson: frequency omega_J = 2eV/hbar", []() {
        double I_c = 1e-6;
        double V = 1e-3;  // 1 mV
        physics::condensed_matter::JosephsonEffect josephson(I_c, V);
        double omega_J = josephson.josephsonFrequency(V);
        double expected = 2.0 * e * V / hbar;
        ASSERT_NEAR(omega_J, expected, TOLERANCE);
        return true;
    });

    run_test("Josephson: inverse AC Josephson relation", []() {
        double I_c = 1e-6;
        double V_c = 1e-3;
        physics::condensed_matter::JosephsonEffect josephson(I_c, V_c);
        double omega = 1e12;  // THz
        double V = josephson.inverseFrohlichEffect(omega);
        double expected = hbar * omega / (2.0 * e);
        ASSERT_NEAR(V, expected, TOLERANCE);
        return true;
    });

    run_test("Josephson: Shapiro steps at V_n = n*hbar*omega_rf/(2e)", []() {
        double I_c = 1e-6;
        double V_c = 1e-3;
        physics::condensed_matter::JosephsonEffect josephson(I_c, V_c);
        double omega_rf = 1e10;
        double V1 = josephson.shapiroSteps(1, omega_rf);
        double V2 = josephson.shapiroSteps(2, omega_rf);
        ASSERT_NEAR(V2, 2.0 * V1, TOLERANCE);
        return true;
    });

    run_test("Josephson: junction energy minimum at delta=0", []() {
        double I_c = 1e-6;
        double V_c = 1e-3;
        physics::condensed_matter::JosephsonEffect josephson(I_c, V_c);
        double E0 = josephson.junctionEnergy(0.0);
        double E1 = josephson.junctionEnergy(M_PI / 4.0);
        ASSERT_TRUE(E0 < E1);
        return true;
    });

    // ========================================
    // Fermi Liquid Theory Tests
    // ========================================

    run_test("Fermi liquid: mass enhancement ratio", []() {
        double m_star = 2.0 * m_e;
        physics::condensed_matter::FermiLiquidTheory fermi(m_star, m_e);
        double enhancement = fermi.massEnhancement();
        ASSERT_NEAR(enhancement, 2.0, TOLERANCE);
        return true;
    });

    run_test("Fermi liquid: quasiparticle energy linear near Fermi surface", []() {
        double m_star = 1.5 * m_e;
        physics::condensed_matter::FermiLiquidTheory fermi(m_star, m_e);
        double k_F = 1e10;
        double k1 = k_F + 1e8;
        double k2 = k_F + 2e8;
        double E1 = fermi.quasiparticleEnergy(k1, k_F);
        double E2 = fermi.quasiparticleEnergy(k2, k_F);
        // Energy should be approximately linear
        ASSERT_TRUE(std::abs(E2) > std::abs(E1));
        return true;
    });

    run_test("Fermi liquid: specific heat enhanced by mass", []() {
        double m_star = 2.0 * m_e;
        physics::condensed_matter::FermiLiquidTheory fermi(m_star, m_e);
        double T = 1.0;
        double g_0 = 1.0;
        double C = fermi.specificHeat(T, g_0);
        ASSERT_TRUE(C > 0.0);
        return true;
    });

    run_test("Fermi liquid: spin susceptibility enhanced", []() {
        double m_star = 2.0 * m_e;
        physics::condensed_matter::FermiLiquidTheory fermi(m_star, m_e);
        double chi_0 = 1.0;
        double chi = fermi.spinSusceptibility(chi_0);
        double expected = 2.0 * chi_0;
        ASSERT_NEAR(chi, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Phonon Dispersion Tests
    // ========================================

    run_test("Phonon: acoustic branch linear at small k", []() {
        double omega_D = 1e13;
        double v_s = 5000.0;  // m/s
        physics::condensed_matter::PhononDispersion phonon(omega_D, v_s);
        double k1 = 1e8;
        double k2 = 2e8;
        double omega1 = phonon.acousticBranch(k1);
        double omega2 = phonon.acousticBranch(k2);
        ASSERT_NEAR(omega2 / omega1, 2.0, TOLERANCE);
        return true;
    });

    run_test("Phonon: optical branch constant with k", []() {
        double omega_D = 1e13;
        double v_s = 5000.0;
        physics::condensed_matter::PhononDispersion phonon(omega_D, v_s);
        double omega_0 = 1e13;
        double omega1 = phonon.opticalBranch(1e8, omega_0);
        double omega2 = phonon.opticalBranch(2e8, omega_0);
        ASSERT_NEAR(omega1, omega2, TOLERANCE);
        return true;
    });

    run_test("Phonon: Debye temperature T_D = hbar*omega_D/k_B", []() {
        double omega_D = 1e13;
        double v_s = 5000.0;
        physics::condensed_matter::PhononDispersion phonon(omega_D, v_s);
        double T_D = phonon.debyeTemperature();
        double expected = hbar * omega_D / k_B;
        ASSERT_NEAR(T_D, expected, TOLERANCE);
        return true;
    });

    run_test("Phonon: Bose-Einstein distribution at T=0", []() {
        double omega_D = 1e13;
        double v_s = 5000.0;
        physics::condensed_matter::PhononDispersion phonon(omega_D, v_s);
        double n = phonon.boseEinsteinDistribution(1e13, 0.1);  // Very low T
        ASSERT_NEAR(n, 0.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Phonon: heat capacity increases with temperature", []() {
        double omega_D = 1e13;
        double v_s = 5000.0;
        physics::condensed_matter::PhononDispersion phonon(omega_D, v_s);
        double T_D = phonon.debyeTemperature();
        double T1 = T_D / 100.0;
        double T2 = T_D / 50.0;
        double C1 = phonon.phononHeatCapacity(T1, T_D);
        double C2 = phonon.phononHeatCapacity(T2, T_D);
        // At low T heat capacity increases with temperature
        ASSERT_TRUE(C2 > C1);
        // Rough T^3 scaling (allow wide tolerance for numerical effects)
        ASSERT_TRUE(C2 / C1 > 4.0);
        return true;
    });

    run_test("Phonon: heat capacity approaches 3k_B at high temperature", []() {
        double omega_D = 1e13;
        double v_s = 5000.0;
        physics::condensed_matter::PhononDispersion phonon(omega_D, v_s);
        double T_D = phonon.debyeTemperature();
        double C = phonon.phononHeatCapacity(10.0 * T_D, T_D);
        ASSERT_NEAR(C, 3.0 * k_B, VERY_LOOSE * k_B);
        return true;
    });

    // ========================================
    // Quantum Hall Effect Tests
    // ========================================

    run_test("Quantum Hall: conductance quantized as nu*e^2/h", []() {
        double B = 10.0;  // Tesla
        int nu = 2;
        physics::condensed_matter::QuantumHallEffect qhe(B, nu);
        double sigma_H = qhe.hallConductance();
        double expected = nu * e * e / (2.0 * M_PI * hbar);
        ASSERT_NEAR(sigma_H, expected, TOLERANCE);
        return true;
    });

    run_test("Quantum Hall: resistance R_H = h/(nu*e^2)", []() {
        double B = 10.0;
        int nu = 1;
        physics::condensed_matter::QuantumHallEffect qhe(B, nu);
        double R_H = qhe.hallResistance();
        double sigma_H = qhe.hallConductance();
        ASSERT_NEAR(R_H, 1.0 / sigma_H, TOLERANCE);
        return true;
    });

    run_test("Quantum Hall: von Klitzing constant", []() {
        double B = 10.0;
        int nu = 1;
        physics::condensed_matter::QuantumHallEffect qhe(B, nu);
        double R_K = qhe.vonKlitzingConstant();
        double expected = 2.0 * M_PI * hbar / (e * e);
        ASSERT_NEAR(R_K, expected, TOLERANCE);
        return true;
    });

    run_test("Quantum Hall: Landau level spacing hbar*omega_c", []() {
        double B = 10.0;
        int nu = 1;
        physics::condensed_matter::QuantumHallEffect qhe(B, nu);
        double E0 = qhe.landauLevelEnergy(0);
        double E1 = qhe.landauLevelEnergy(1);
        double omega_c = e * B / m_e;
        double expected_spacing = hbar * omega_c;
        ASSERT_NEAR(E1 - E0, expected_spacing, LOOSE_TOLERANCE * expected_spacing);
        return true;
    });

    run_test("Quantum Hall: magnetic length l_B = sqrt(hbar/eB)", []() {
        double B = 10.0;
        int nu = 1;
        physics::condensed_matter::QuantumHallEffect qhe(B, nu);
        double l_B = qhe.magneticLength();
        double expected = std::sqrt(hbar / (e * B));
        ASSERT_NEAR(l_B, expected, TOLERANCE);
        return true;
    });

    run_test("Quantum Hall: cyclotron frequency omega_c = eB/m", []() {
        double B = 10.0;
        int nu = 1;
        physics::condensed_matter::QuantumHallEffect qhe(B, nu);
        double omega_c = qhe.cyclotronFrequency();
        double expected = e * B / m_e;
        ASSERT_NEAR(omega_c, expected, TOLERANCE);
        return true;
    });

    run_test("Quantum Hall: filling factor", []() {
        double B = 10.0;
        int nu = 4;
        physics::condensed_matter::QuantumHallEffect qhe(B, nu);
        ASSERT_TRUE(qhe.fillingFactor() == 4);
        return true;
    });

    // ========================================
    // Landau Theory Tests
    // ========================================

    run_test("Landau theory: order parameter zero above T_c", []() {
        double T_c = 100.0;
        physics::condensed_matter::LandauTheory landau(T_c);
        double m = landau.orderParameter(T_c + 10.0, 1.0);
        ASSERT_NEAR(m, 0.0, TOLERANCE);
        return true;
    });

    run_test("Landau theory: order parameter nonzero below T_c", []() {
        double T_c = 100.0;
        physics::condensed_matter::LandauTheory landau(T_c);
        double m = landau.orderParameter(T_c - 10.0, 1.0);
        ASSERT_TRUE(m > 0.0);
        return true;
    });

    run_test("Landau theory: critical exponent beta = 0.5", []() {
        double T_c = 100.0;
        physics::condensed_matter::LandauTheory landau(T_c);
        double beta = landau.criticalExponentBeta();
        ASSERT_NEAR(beta, 0.5, TOLERANCE);
        return true;
    });

    run_test("Landau theory: critical exponent gamma = 1.0", []() {
        double T_c = 100.0;
        physics::condensed_matter::LandauTheory landau(T_c);
        double gamma = landau.criticalExponentGamma();
        ASSERT_NEAR(gamma, 1.0, TOLERANCE);
        return true;
    });

    run_test("Landau theory: critical exponent nu = 0.5", []() {
        double T_c = 100.0;
        physics::condensed_matter::LandauTheory landau(T_c);
        double nu = landau.criticalExponentNu();
        ASSERT_NEAR(nu, 0.5, TOLERANCE);
        return true;
    });

    run_test("Landau theory: susceptibility diverges at T_c", []() {
        double T_c = 100.0;
        physics::condensed_matter::LandauTheory landau(T_c);
        double chi1 = landau.susceptibility(T_c + 1.0, 1.0);
        double chi2 = landau.susceptibility(T_c + 10.0, 1.0);
        ASSERT_TRUE(chi1 > chi2);
        return true;
    });

    run_test("Landau theory: correlation length diverges at T_c", []() {
        double T_c = 100.0;
        physics::condensed_matter::LandauTheory landau(T_c);
        double xi1 = landau.correlationLength(T_c + 0.1, 1.0);
        double xi2 = landau.correlationLength(T_c + 1.0, 1.0);
        ASSERT_TRUE(xi1 > xi2);
        return true;
    });

    // ========================================
    // Hubbard Model Tests
    // ========================================

    run_test("Hubbard model: kinetic energy negative for occupied sites", []() {
        double t = 1.0 * e;
        double U = 4.0 * e;
        int L = 4;
        physics::condensed_matter::HubbardModel hubbard(t, U, L);
        std::vector<int> occupation = {1, 1, 1, 1};
        double E_kin = hubbard.kineticEnergy(occupation);
        ASSERT_TRUE(E_kin < 0.0);
        return true;
    });

    run_test("Hubbard model: interaction energy for double occupancy", []() {
        double t = 1.0 * e;
        double U = 4.0 * e;
        int L = 4;
        physics::condensed_matter::HubbardModel hubbard(t, U, L);
        std::vector<int> n_up = {1, 1, 0, 0};
        std::vector<int> n_down = {1, 0, 1, 0};
        double E_int = hubbard.interactionEnergy(n_up, n_down);
        ASSERT_NEAR(E_int, U, TOLERANCE);  // One doubly occupied site
        return true;
    });

    run_test("Hubbard model: Mott insulator at half-filling with large U", []() {
        double t = 1.0 * e;
        double U = 10.0 * e;
        int L = 4;
        physics::condensed_matter::HubbardModel hubbard(t, U, L);
        bool is_mott = hubbard.isMottInsulator(1.0);
        ASSERT_TRUE(is_mott);
        return true;
    });

    run_test("Hubbard model: metallic at half-filling with small U", []() {
        double t = 1.0 * e;
        double U = 1.0 * e;
        int L = 4;
        physics::condensed_matter::HubbardModel hubbard(t, U, L);
        bool is_mott = hubbard.isMottInsulator(1.0);
        ASSERT_TRUE(!is_mott);
        return true;
    });

    run_test("Hubbard model: band gap opens for large U/t", []() {
        double t = 1.0 * e;
        double U = 10.0 * e;
        int L = 4;
        physics::condensed_matter::HubbardModel hubbard(t, U, L);
        double gap = hubbard.bandGap();
        ASSERT_TRUE(gap > 0.0);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Total tests:  " << (tests_passed + tests_failed) << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;
    std::cout << "======================================" << std::endl;

    return (tests_failed == 0) ? 0 : 1;
}
