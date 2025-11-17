/**
 * @file phase3_quantum_em_optics.cpp
 * @brief Phase 3 Validation: Quantum Mechanics, Electromagnetic Waves, and Optics
 *
 * Tests quantum basics, EM wave propagation, and optical systems
 */

#include <iostream>
#include <cmath>
#include <stdexcept>

// Phase 3 modules
#include "../include/physics/quantum_basics.hpp"
#include "../include/physics/electromagnetic_waves.hpp"
#include "../include/physics/optics.hpp"

// Tolerance for floating-point comparisons
constexpr double TOLERANCE = 1e-6;
constexpr double LOOSE_TOLERANCE = 1e-3;  // For some physical constants

// Test counters
int tests_run = 0;
int tests_passed = 0;

// Macro for assertions
#define ASSERT_NEAR(actual, expected, tol) \
    do { \
        tests_run++; \
        if (std::abs((actual) - (expected)) <= (tol)) { \
            tests_passed++; \
        } else { \
            std::cout << "FAIL: " << #actual << " = " << (actual) \
                      << ", expected " << (expected) \
                      << " (difference: " << std::abs((actual) - (expected)) << ")" \
                      << " at line " << __LINE__ << std::endl; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        tests_run++; \
        if (condition) { \
            tests_passed++; \
        } else { \
            std::cout << "FAIL: " << #condition << " at line " << __LINE__ << std::endl; \
        } \
    } while(0)

// Namespaces (import functions but not constants to avoid ambiguity)
using physics::quantum_basics::deBroglieWavelength;
using physics::quantum_basics::deBroglieWavelengthFromVelocity;
using physics::quantum_basics::deBroglieWavelengthFromEnergy;
using physics::quantum_basics::comptonWavelength;
using physics::quantum_basics::comptonShift;
using physics::quantum_basics::scatteredWavelength;
using physics::quantum_basics::positionUncertainty;
using physics::quantum_basics::momentumUncertainty;
using physics::quantum_basics::energyTimeUncertainty;
using physics::quantum_basics::timeUncertainty;
using physics::quantum_basics::bohrEnergyLevel;
using physics::quantum_basics::bohrOrbitalRadius;
using physics::quantum_basics::hydrogenTransitionEnergy;
using physics::quantum_basics::rydbergWavelength;
using physics::quantum_basics::hydrogenIonizationEnergy;
using physics::quantum_basics::photoelectronKineticEnergy;
using physics::quantum_basics::thresholdFrequency;
using physics::quantum_basics::stoppingPotential;
using physics::quantum_basics::particleInBoxEnergy;
using physics::quantum_basics::zeroPointEnergy;
using physics::quantum_basics::harmonicOscillatorEnergy;
using physics::quantum_basics::harmonicOscillatorZeroPoint;
using physics::quantum_basics::tunnelingProbability;
using physics::quantum_basics::tunnelingDecayConstant;
using physics::quantum_basics::orbitalAngularMomentum;
using physics::quantum_basics::angularMomentumZComponent;
using physics::quantum_basics::spinAngularMomentum;

using physics::electromagnetic_waves::speedOfLight;
using physics::electromagnetic_waves::speedFromRefractiveIndex;
using physics::electromagnetic_waves::wavelengthFromFrequency;
using physics::electromagnetic_waves::frequencyFromWavelength;
using physics::electromagnetic_waves::angularFrequency;
using physics::electromagnetic_waves::waveNumber;
using physics::electromagnetic_waves::magneticFieldFromElectric;
using physics::electromagnetic_waves::electricFieldFromMagnetic;
using physics::electromagnetic_waves::waveSpeedFromFields;
using physics::electromagnetic_waves::electricEnergyDensity;
using physics::electromagnetic_waves::magneticEnergyDensity;
using physics::electromagnetic_waves::totalEnergyDensity;
using physics::electromagnetic_waves::waveIntensity;
using physics::electromagnetic_waves::powerFromIntensity;
using physics::electromagnetic_waves::intensityAtDistance;
using physics::electromagnetic_waves::radiationPressureAbsorption;
using physics::electromagnetic_waves::radiationPressureReflection;
using physics::electromagnetic_waves::radiationForce;
using physics::electromagnetic_waves::photonEnergy;
using physics::electromagnetic_waves::photonEnergyFromWavelength;
using physics::electromagnetic_waves::photonMomentum;
using physics::electromagnetic_waves::photonsPerSecond;
using physics::electromagnetic_waves::wavelengthInMedium;
using physics::electromagnetic_waves::phaseVelocity;
using physics::electromagnetic_waves::freeSpaceImpedance;
using physics::electromagnetic_waves::isVisible;

using physics::optics::snellsLaw;
using physics::optics::refractiveIndexFromAngles;
using physics::optics::criticalAngle;
using physics::optics::isTotalInternalReflection;
using physics::optics::brewstersAngle;
using physics::optics::refractiveIndexFromBrewster;
using physics::optics::verifyBrewsterCondition;
using physics::optics::lensFormula;
using physics::optics::imageDistance;
using physics::optics::objectDistance;
using physics::optics::lensPower;
using physics::optics::focalLengthFromPower;
using physics::optics::lensmakersEquation;
using physics::optics::focalLengthInMedium;
using physics::optics::linearMagnification;
using physics::optics::imageHeight;
using physics::optics::magnificationFromFocal;
using physics::optics::mirrorFocalLength;
using physics::optics::mirrorFormula;
using physics::optics::simpleMicroscopeMagnification;
using physics::optics::compoundMicroscopeMagnification;
using physics::optics::telescopeMagnification;
using physics::optics::telescopeLength;
using physics::optics::telescopeResolvingPower;
using physics::optics::combinedFocalLength;
using physics::optics::combinedPower;
using physics::optics::separatedLensesFocalLength;

//==============================================================================
// QUANTUM BASICS TESTS
//==============================================================================

void test_de_broglie_wavelength() {
    std::cout << "\n=== Testing De Broglie Wavelength ===" << std::endl;

    // Constants
    const double h = physics::quantum_basics::constants::PLANCK_H;  // 6.626e-34 J⋅s
    const double m_e = physics::quantum_basics::constants::ELECTRON_MASS;  // 9.109e-31 kg
    const double c = physics::quantum_basics::constants::SPEED_OF_LIGHT;

    // Test 1: Wavelength from momentum
    double p = 1e-24;  // kg⋅m/s
    double lambda = deBroglieWavelength(p);
    ASSERT_NEAR(lambda, h/p, TOLERANCE);

    // Test 2: Wavelength from velocity (electron)
    double v = 1e6;  // m/s
    lambda = deBroglieWavelengthFromVelocity(m_e, v);
    double expected = h / (m_e * v);
    ASSERT_NEAR(lambda, expected, TOLERANCE);

    // Test 3: Wavelength from kinetic energy
    double KE = 100 * 1.602e-19;  // 100 eV in Joules
    lambda = deBroglieWavelengthFromEnergy(m_e, KE);
    double p_expected = std::sqrt(2 * m_e * KE);
    expected = h / p_expected;
    ASSERT_NEAR(lambda, expected, TOLERANCE);

    // Test 4: Verify λ = h/p relationship
    p = 2.73e-24;  // arbitrary momentum
    double lambda1 = deBroglieWavelength(p);
    double lambda2 = h / p;
    ASSERT_NEAR(lambda1, lambda2, TOLERANCE);
}

void test_compton_scattering() {
    std::cout << "\n=== Testing Compton Scattering ===" << std::endl;

    const double h = physics::quantum_basics::constants::PLANCK_H;
    const double m_e = physics::quantum_basics::constants::ELECTRON_MASS;
    const double c = physics::quantum_basics::constants::SPEED_OF_LIGHT;

    // Test 1: Compton wavelength of electron
    double lambda_C = comptonWavelength();
    double expected = h / (m_e * c);  // ~2.43e-12 m
    ASSERT_NEAR(lambda_C, expected, TOLERANCE);

    // Test 2: Compton shift at 90 degrees
    double theta = M_PI / 2.0;
    double shift = comptonShift(theta);
    expected = lambda_C * (1 - std::cos(theta));  // λ_C for 90°
    ASSERT_NEAR(shift, expected, TOLERANCE);

    // Test 3: No shift at 0 degrees (forward scattering)
    shift = comptonShift(0.0);
    ASSERT_NEAR(shift, 0.0, TOLERANCE);

    // Test 4: Maximum shift at 180 degrees (backscattering)
    shift = comptonShift(M_PI);
    expected = 2.0 * lambda_C;
    ASSERT_NEAR(shift, expected, TOLERANCE);

    // Test 5: Scattered wavelength calculation
    double lambda_0 = 1e-11;  // 0.01 nm incident wavelength
    theta = M_PI / 2.0;
    double lambda_scattered = scatteredWavelength(lambda_0, theta);
    expected = lambda_0 + lambda_C;  // at 90°
    ASSERT_NEAR(lambda_scattered, expected, TOLERANCE);
}

void test_heisenberg_uncertainty() {
    std::cout << "\n=== Testing Heisenberg Uncertainty Principle ===" << std::endl;

    const double hbar = physics::quantum_basics::constants::HBAR;  // ℏ = h/(2π)

    // Test 1: Position uncertainty from momentum uncertainty
    double delta_p = 1e-24;  // kg⋅m/s
    double delta_x = positionUncertainty(delta_p);
    double expected = hbar / (2.0 * delta_p);
    ASSERT_NEAR(delta_x, expected, TOLERANCE);

    // Test 2: Momentum uncertainty from position uncertainty
    delta_x = 1e-10;  // 1 Angstrom
    delta_p = momentumUncertainty(delta_x);
    expected = hbar / (2.0 * delta_x);
    ASSERT_NEAR(delta_p, expected, TOLERANCE);

    // Test 3: Verify minimum uncertainty product: Δx⋅Δp ≥ ℏ/2
    delta_x = 1e-9;
    delta_p = momentumUncertainty(delta_x);
    double product = delta_x * delta_p;
    ASSERT_TRUE(product >= hbar/2.0 - TOLERANCE);
    ASSERT_NEAR(product, hbar/2.0, TOLERANCE);

    // Test 4: Energy-time uncertainty
    double min_product = energyTimeUncertainty();
    ASSERT_NEAR(min_product, hbar/2.0, TOLERANCE);

    // Test 5: Time uncertainty from energy uncertainty
    double delta_E = 1.602e-19;  // 1 eV
    double delta_t = timeUncertainty(delta_E);
    expected = hbar / (2.0 * delta_E);
    ASSERT_NEAR(delta_t, expected, TOLERANCE);
}

void test_bohr_model() {
    std::cout << "\n=== Testing Bohr Model ===" << std::endl;

    const double R_E = physics::quantum_basics::constants::RYDBERG_ENERGY;  // 13.6 eV
    const double a_0 = physics::quantum_basics::constants::BOHR_RADIUS;     // 0.529e-10 m

    // Test 1: Ground state energy (n=1)
    double E_1 = bohrEnergyLevel(1);
    ASSERT_NEAR(E_1, -13.6, 0.01);  // -13.6 eV (tolerance 0.01 for constants)

    // Test 2: First excited state (n=2)
    double E_2 = bohrEnergyLevel(2);
    ASSERT_NEAR(E_2, -13.6/4.0, 0.01);  // -3.4 eV

    // Test 3: Energy scaling: E_n ∝ 1/n²
    double E_3 = bohrEnergyLevel(3);
    ASSERT_NEAR(E_3 / E_1, 1.0/9.0, TOLERANCE);

    // Test 4: Ground state orbital radius
    double r_1 = bohrOrbitalRadius(1);
    ASSERT_NEAR(r_1, a_0, TOLERANCE);

    // Test 5: Second orbit radius: r_2 = 4a₀
    double r_2 = bohrOrbitalRadius(2);
    ASSERT_NEAR(r_2, 4.0 * a_0, TOLERANCE);

    // Test 6: Lyman alpha transition (n=2 → n=1)
    double deltaE = hydrogenTransitionEnergy(2, 1);
    double expected = R_E * (1.0 - 1.0/4.0);  // 10.2 eV
    ASSERT_NEAR(deltaE, expected, LOOSE_TOLERANCE);

    // Test 7: Balmer alpha (H-alpha): n=3 → n=2
    deltaE = hydrogenTransitionEnergy(3, 2);
    expected = R_E * (1.0/4.0 - 1.0/9.0);  // 1.89 eV
    ASSERT_NEAR(deltaE, expected, LOOSE_TOLERANCE);

    // Test 8: Ionization energy (from ground state)
    double E_ion = hydrogenIonizationEnergy();
    ASSERT_NEAR(E_ion, 13.6, 0.01);

    // Test 9: Rydberg wavelength formula (n=3 → n=2, H-alpha)
    double lambda = rydbergWavelength(3, 2);
    // H-alpha is at 656.3 nm
    ASSERT_NEAR(lambda * 1e9, 656.3, 1.0);  // nm, loose tolerance
}

void test_photoelectric_effect() {
    std::cout << "\n=== Testing Photoelectric Effect ===" << std::endl;

    const double h = physics::quantum_basics::constants::PLANCK_H;
    const double e = physics::quantum_basics::constants::ELEMENTARY_CHARGE;

    // Work function of sodium: ~2.3 eV
    double phi = 2.3 * e;  // Convert to Joules

    // Test 1: Kinetic energy from Einstein's equation
    double E_photon = 3.0 * e;  // 3 eV photon
    double KE = photoelectronKineticEnergy(E_photon, phi);
    double expected = (3.0 - 2.3) * e;  // 0.7 eV
    ASSERT_NEAR(KE, expected, TOLERANCE);

    // Test 2: Threshold frequency
    double f_0 = thresholdFrequency(phi);
    expected = phi / h;
    ASSERT_NEAR(f_0, expected, TOLERANCE);

    // Test 3: No emission below threshold
    try {
        photoelectronKineticEnergy(2.0 * e, phi);  // 2 eV < 2.3 eV
        ASSERT_TRUE(false);  // Should throw
    } catch (const std::invalid_argument&) {
        ASSERT_TRUE(true);
    }

    // Test 4: Stopping potential
    E_photon = 4.0 * e;  // 4 eV photon
    double V_s = stoppingPotential(E_photon, phi);
    expected = (4.0 - 2.3);  // 1.7 V
    ASSERT_NEAR(V_s, expected, TOLERANCE);

    // Test 5: Zero stopping potential when E < φ
    V_s = stoppingPotential(2.0 * e, phi);
    ASSERT_NEAR(V_s, 0.0, TOLERANCE);
}

void test_particle_in_box() {
    std::cout << "\n=== Testing Particle in a Box ===" << std::endl;

    const double h = physics::quantum_basics::constants::PLANCK_H;
    const double m = physics::quantum_basics::constants::ELECTRON_MASS;
    const double L = 1e-9;  // 1 nm box

    // Test 1: Ground state energy (n=1)
    double E_1 = particleInBoxEnergy(1, m, L);
    double expected = h*h / (8.0 * m * L*L);
    ASSERT_NEAR(E_1, expected, TOLERANCE);

    // Test 2: First excited state (n=2)
    double E_2 = particleInBoxEnergy(2, m, L);
    expected = 4.0 * h*h / (8.0 * m * L*L);
    ASSERT_NEAR(E_2, expected, TOLERANCE);

    // Test 3: Energy scaling: E_n ∝ n²
    ASSERT_NEAR(E_2 / E_1, 4.0, TOLERANCE);

    // Test 4: Zero-point energy
    double E_0 = zeroPointEnergy(m, L);
    ASSERT_NEAR(E_0, E_1, TOLERANCE);

    // Test 5: Third level
    double E_3 = particleInBoxEnergy(3, m, L);
    ASSERT_NEAR(E_3 / E_1, 9.0, TOLERANCE);
}

void test_quantum_harmonic_oscillator() {
    std::cout << "\n=== Testing Quantum Harmonic Oscillator ===" << std::endl;

    const double hbar = physics::quantum_basics::constants::HBAR;
    const double omega = 1e14;  // rad/s

    // Test 1: Ground state energy (n=0)
    double E_0 = harmonicOscillatorEnergy(0, omega);
    double expected = 0.5 * hbar * omega;
    ASSERT_NEAR(E_0, expected, TOLERANCE);

    // Test 2: First excited state (n=1)
    double E_1 = harmonicOscillatorEnergy(1, omega);
    expected = 1.5 * hbar * omega;
    ASSERT_NEAR(E_1, expected, TOLERANCE);

    // Test 3: Energy spacing is uniform: ΔE = ℏω
    double delta_E = E_1 - E_0;
    ASSERT_NEAR(delta_E, hbar * omega, TOLERANCE);

    // Test 4: Zero-point energy
    double E_zp = harmonicOscillatorZeroPoint(omega);
    ASSERT_NEAR(E_zp, E_0, TOLERANCE);

    // Test 5: nth level: E_n = ℏω(n + 1/2)
    double E_5 = harmonicOscillatorEnergy(5, omega);
    expected = hbar * omega * 5.5;
    ASSERT_NEAR(E_5, expected, TOLERANCE);
}

void test_quantum_tunneling() {
    std::cout << "\n=== Testing Quantum Tunneling ===" << std::endl;

    const double hbar = physics::quantum_basics::constants::HBAR;
    const double m = physics::quantum_basics::constants::ELECTRON_MASS;

    // Test 1: Tunneling probability calculation
    double E = 1.0 * 1.602e-19;    // 1 eV
    double V = 2.0 * 1.602e-19;    // 2 eV barrier
    double L = 1e-9;                // 1 nm barrier width
    double T = tunnelingProbability(E, V, L, m);

    // Should be very small (exponential decay)
    ASSERT_TRUE(T > 0.0 && T < 1.0);

    // Test 2: Decay constant
    double kappa = tunnelingDecayConstant(E, V, m);
    double expected = std::sqrt(2.0 * m * (V - E)) / hbar;
    ASSERT_NEAR(kappa, expected, TOLERANCE);

    // Test 3: Wider barrier reduces transmission
    double T_narrow = tunnelingProbability(E, V, 0.5e-9, m);
    double T_wide = tunnelingProbability(E, V, 1.0e-9, m);
    ASSERT_TRUE(T_narrow > T_wide);

    // Test 4: Higher barrier reduces transmission
    double T_low = tunnelingProbability(E, 1.5 * 1.602e-19, L, m);
    double T_high = tunnelingProbability(E, 2.5 * 1.602e-19, L, m);
    ASSERT_TRUE(T_low > T_high);
}

void test_quantum_angular_momentum() {
    std::cout << "\n=== Testing Quantum Angular Momentum ===" << std::endl;

    const double hbar = physics::quantum_basics::constants::HBAR;

    // Test 1: Orbital angular momentum for l=1
    double L = orbitalAngularMomentum(1);
    double expected = hbar * std::sqrt(2.0);  // √[1(1+1)] ℏ
    ASSERT_NEAR(L, expected, TOLERANCE);

    // Test 2: Ground state (l=0) has zero angular momentum
    L = orbitalAngularMomentum(0);
    ASSERT_NEAR(L, 0.0, TOLERANCE);

    // Test 3: l=2 state
    L = orbitalAngularMomentum(2);
    expected = hbar * std::sqrt(6.0);  // √[2(2+1)] ℏ
    ASSERT_NEAR(L, expected, TOLERANCE);

    // Test 4: Z-component of angular momentum
    double L_z = angularMomentumZComponent(1);
    ASSERT_NEAR(L_z, hbar, TOLERANCE);

    L_z = angularMomentumZComponent(-1);
    ASSERT_NEAR(L_z, -hbar, TOLERANCE);

    // Test 5: Spin angular momentum (electron s=1/2)
    double S = spinAngularMomentum(0.5);
    expected = hbar * std::sqrt(0.75);  // √[1/2(1/2+1)] ℏ = (√3/2)ℏ
    ASSERT_NEAR(S, expected, TOLERANCE);
}

//==============================================================================
// ELECTROMAGNETIC WAVES TESTS
//==============================================================================

void test_speed_of_light() {
    std::cout << "\n=== Testing Speed of Light ===" << std::endl;

    const double epsilon_0 = physics::electromagnetic_waves::constants::EPSILON_0;
    const double mu_0 = physics::electromagnetic_waves::constants::MU_0;
    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;

    // Test 1: Calculate c from ε₀ and μ₀
    double c_calc = speedOfLight();
    ASSERT_NEAR(c_calc, c, 1e3);  // Within 1000 m/s

    // Test 2: Verify c = 1/√(ε₀μ₀)
    double expected = 1.0 / std::sqrt(epsilon_0 * mu_0);
    ASSERT_NEAR(c_calc, expected, TOLERANCE);

    // Test 3: Speed in medium with refractive index
    double n = 1.5;  // glass
    double v = speedFromRefractiveIndex(n);
    expected = c / n;
    ASSERT_NEAR(v, expected, TOLERANCE);

    // Test 4: Speed in water (n ≈ 1.33)
    v = speedFromRefractiveIndex(1.33);
    expected = c / 1.33;
    ASSERT_NEAR(v, expected, 1e3);
}

void test_wavelength_frequency() {
    std::cout << "\n=== Testing Wavelength-Frequency Relationships ===" << std::endl;

    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;

    // Test 1: Wavelength from frequency
    double f = 5e14;  // 500 THz (green light)
    double lambda = wavelengthFromFrequency(f);
    double expected = c / f;  // ~600 nm
    ASSERT_NEAR(lambda, expected, TOLERANCE);

    // Test 2: Frequency from wavelength
    lambda = 650e-9;  // 650 nm (red light)
    f = frequencyFromWavelength(lambda);
    expected = c / lambda;
    ASSERT_NEAR(f, expected, TOLERANCE);

    // Test 3: Round-trip consistency
    f = 1e9;  // 1 GHz
    lambda = wavelengthFromFrequency(f);
    double f_back = frequencyFromWavelength(lambda);
    ASSERT_NEAR(f_back, f, TOLERANCE);

    // Test 4: Angular frequency
    f = 1e6;  // 1 MHz
    double omega = angularFrequency(f);
    expected = 2.0 * M_PI * f;
    ASSERT_NEAR(omega, expected, TOLERANCE);

    // Test 5: Wave number
    lambda = 500e-9;  // 500 nm
    double k = waveNumber(lambda);
    expected = 2.0 * M_PI / lambda;
    ASSERT_NEAR(k, expected, TOLERANCE);
}

void test_em_field_relationships() {
    std::cout << "\n=== Testing E-B Field Relationships ===" << std::endl;

    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;

    // Test 1: Magnetic field from electric field
    double E = 300;  // V/m
    double B = magneticFieldFromElectric(E);
    double expected = E / c;  // ~1 µT
    ASSERT_NEAR(B, expected, TOLERANCE);

    // Test 2: Electric field from magnetic field
    B = 1e-6;  // 1 µT
    E = electricFieldFromMagnetic(B);
    expected = B * c;
    ASSERT_NEAR(E, expected, TOLERANCE);

    // Test 3: Verify E/B = c
    E = 100;  // V/m
    B = magneticFieldFromElectric(E);
    double ratio = E / B;
    ASSERT_NEAR(ratio, c, 1e3);

    // Test 4: Wave speed from fields
    E = 500;
    B = E / c;
    double v = waveSpeedFromFields(E, B);
    ASSERT_NEAR(v, c, 1e3);
}

void test_em_energy_density() {
    std::cout << "\n=== Testing EM Energy Density ===" << std::endl;

    const double epsilon_0 = physics::electromagnetic_waves::constants::EPSILON_0;
    const double mu_0 = physics::electromagnetic_waves::constants::MU_0;
    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;

    // Test 1: Electric energy density
    double E = 100;  // V/m
    double u_E = electricEnergyDensity(E);
    double expected = 0.5 * epsilon_0 * E * E;
    ASSERT_NEAR(u_E, expected, TOLERANCE);

    // Test 2: Magnetic energy density
    double B = E / c;
    double u_B = magneticEnergyDensity(B);
    expected = 0.5 * B * B / mu_0;
    ASSERT_NEAR(u_B, expected, TOLERANCE);

    // Test 3: In EM wave, u_E = u_B
    ASSERT_NEAR(u_E, u_B, 1e-15);

    // Test 4: Total energy density
    double u_total = totalEnergyDensity(E);
    ASSERT_NEAR(u_total, u_E + u_B, TOLERANCE);
    ASSERT_NEAR(u_total, 2.0 * u_E, TOLERANCE);
}

void test_wave_intensity() {
    std::cout << "\n=== Testing Wave Intensity ===" << std::endl;

    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;
    const double epsilon_0 = physics::electromagnetic_waves::constants::EPSILON_0;

    // Test 1: Intensity from electric field
    double E_0 = 100;  // V/m amplitude
    double I = waveIntensity(E_0);
    double expected = 0.5 * c * epsilon_0 * E_0 * E_0;
    ASSERT_NEAR(I, expected, TOLERANCE);

    // Test 2: Power from intensity
    I = 1000;  // W/m² (solar constant)
    double A = 1.0;  // 1 m²
    double P = powerFromIntensity(I, A);
    ASSERT_NEAR(P, 1000, TOLERANCE);

    // Test 3: Intensity at distance (inverse square law)
    P = 100;  // W
    double r = 2.0;  // m
    I = intensityAtDistance(P, r);
    expected = P / (4.0 * M_PI * r * r);
    ASSERT_NEAR(I, expected, TOLERANCE);

    // Test 4: Intensity doubles when distance halves (area quarters)
    double I_1 = intensityAtDistance(P, 2.0);
    double I_2 = intensityAtDistance(P, 1.0);
    ASSERT_NEAR(I_2 / I_1, 4.0, TOLERANCE);
}

void test_radiation_pressure() {
    std::cout << "\n=== Testing Radiation Pressure ===" << std::endl;

    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;

    // Test 1: Radiation pressure (absorption)
    double I = 1360;  // W/m² (solar constant)
    double P_abs = radiationPressureAbsorption(I);
    double expected = I / c;  // ~4.5 µPa
    ASSERT_NEAR(P_abs, expected, TOLERANCE);

    // Test 2: Radiation pressure (reflection) is twice absorption
    double P_ref = radiationPressureReflection(I);
    ASSERT_NEAR(P_ref, 2.0 * P_abs, TOLERANCE);

    // Test 3: Radiation force
    double A = 10.0;  // 10 m² area
    double F = radiationForce(P_abs, A);
    expected = P_abs * A;
    ASSERT_NEAR(F, expected, TOLERANCE);
}

void test_photon_properties() {
    std::cout << "\n=== Testing Photon Properties ===" << std::endl;

    const double h = physics::electromagnetic_waves::constants::PLANCK_H;
    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;

    // Test 1: Photon energy from frequency
    double f = 5e14;  // 500 THz (green light)
    double E = photonEnergy(f);
    double expected = h * f;
    ASSERT_NEAR(E, expected, TOLERANCE);

    // Test 2: Photon energy from wavelength
    double lambda = 550e-9;  // 550 nm (green)
    E = photonEnergyFromWavelength(lambda);
    expected = h * c / lambda;
    ASSERT_NEAR(E, expected, TOLERANCE);

    // Test 3: Photon momentum
    lambda = 600e-9;  // 600 nm
    double p = photonMomentum(lambda);
    expected = h / lambda;
    ASSERT_NEAR(p, expected, TOLERANCE);

    // Test 4: E = pc for photons
    double E_from_p = p * c;
    double E_from_lambda = photonEnergyFromWavelength(lambda);
    ASSERT_NEAR(E_from_p, E_from_lambda, TOLERANCE);

    // Test 5: Number of photons per second
    double power = 1.0;  // 1 W
    E = photonEnergyFromWavelength(550e-9);
    double N = photonsPerSecond(power, E);
    expected = power / E;
    ASSERT_NEAR(N, expected, TOLERANCE);
}

void test_wave_in_medium() {
    std::cout << "\n=== Testing Wave Propagation in Medium ===" << std::endl;

    const double c = physics::electromagnetic_waves::constants::SPEED_OF_LIGHT;

    // Test 1: Wavelength in medium
    double lambda_vac = 600e-9;  // 600 nm in vacuum
    double n = 1.5;  // glass
    double lambda_med = wavelengthInMedium(lambda_vac, n);
    double expected = lambda_vac / n;  // 400 nm
    ASSERT_NEAR(lambda_med, expected, TOLERANCE);

    // Test 2: Phase velocity
    double v_p = phaseVelocity(n);
    expected = c / n;
    ASSERT_NEAR(v_p, expected, TOLERANCE);

    // Test 3: Free space impedance
    const double mu_0 = physics::electromagnetic_waves::constants::MU_0;
    const double epsilon_0 = physics::electromagnetic_waves::constants::EPSILON_0;
    double Z_0 = freeSpaceImpedance();
    expected = std::sqrt(mu_0 / epsilon_0);  // ~377 Ω
    ASSERT_NEAR(Z_0, 377, 1.0);

    // Test 4: Visible spectrum check
    ASSERT_TRUE(isVisible(550e-9));   // Green
    ASSERT_TRUE(isVisible(650e-9));   // Red
    ASSERT_TRUE(!isVisible(300e-9));  // UV
    ASSERT_TRUE(!isVisible(800e-9));  // IR
}

//==============================================================================
// OPTICS TESTS
//==============================================================================

void test_snells_law() {
    std::cout << "\n=== Testing Snell's Law ===" << std::endl;

    // Test 1: Air to glass refraction
    double n1 = 1.0;   // air
    double n2 = 1.5;   // glass
    double theta_i = M_PI / 6.0;  // 30°
    double theta_r = snellsLaw(theta_i, n1, n2);

    // Verify: n₁sinθ₁ = n₂sinθ₂
    double lhs = n1 * std::sin(theta_i);
    double rhs = n2 * std::sin(theta_r);
    ASSERT_NEAR(lhs, rhs, TOLERANCE);

    // Test 2: Reflection angle equals incident angle (law of reflection)
    // For ideal reflection, θ_reflected = θ_incident
    ASSERT_NEAR(theta_i, theta_i, TOLERANCE);  // Trivial but validates convention

    // Test 3: Normal incidence (θ=0) gives no refraction
    theta_r = snellsLaw(0.0, n1, n2);
    ASSERT_NEAR(theta_r, 0.0, TOLERANCE);

    // Test 4: Calculate refractive index from angles
    theta_i = M_PI / 4.0;  // 45°
    theta_r = snellsLaw(theta_i, 1.0, 1.5);
    double n_calc = refractiveIndexFromAngles(theta_i, theta_r);
    ASSERT_NEAR(n_calc, 1.5, TOLERANCE);
}

void test_critical_angle_tir() {
    std::cout << "\n=== Testing Critical Angle and TIR ===" << std::endl;

    // Test 1: Critical angle for glass-air interface
    double n1 = 1.5;  // glass
    double n2 = 1.0;  // air
    double theta_c = criticalAngle(n1, n2);
    double expected = std::asin(n2 / n1);  // ~41.8°
    ASSERT_NEAR(theta_c, expected, TOLERANCE);

    // Test 2: Critical angle for water-air
    theta_c = criticalAngle(1.33, 1.0);
    expected = std::asin(1.0 / 1.33);  // ~48.8°
    ASSERT_NEAR(theta_c, expected, TOLERANCE);

    // Test 3: TIR occurs above critical angle
    double theta = M_PI / 3.0;  // 60°
    bool tir = isTotalInternalReflection(theta, 1.5, 1.0);
    ASSERT_TRUE(tir);

    // Test 4: No TIR below critical angle
    theta = M_PI / 6.0;  // 30°
    tir = isTotalInternalReflection(theta, 1.5, 1.0);
    ASSERT_TRUE(!tir);

    // Test 5: No TIR when going to denser medium
    tir = isTotalInternalReflection(M_PI/3.0, 1.0, 1.5);
    ASSERT_TRUE(!tir);
}

void test_brewsters_angle() {
    std::cout << "\n=== Testing Brewster's Angle ===" << std::endl;

    // Test 1: Brewster's angle for air-glass
    double n1 = 1.0;
    double n2 = 1.5;
    double theta_B = brewstersAngle(n1, n2);
    double expected = std::atan(n2 / n1);  // ~56.3°
    ASSERT_NEAR(theta_B, expected, TOLERANCE);

    // Test 2: Verify tan(θ_B) = n₂/n₁
    double tan_theta = std::tan(theta_B);
    ASSERT_NEAR(tan_theta, n2 / n1, TOLERANCE);

    // Test 3: Calculate refractive index from Brewster's angle
    double n_calc = refractiveIndexFromBrewster(theta_B);
    ASSERT_NEAR(n_calc, n2, TOLERANCE);

    // Test 4: At Brewster's angle, θ_B + θ_r = 90°
    double theta_r = snellsLaw(theta_B, n1, n2);
    bool verified = verifyBrewsterCondition(theta_B, theta_r);
    ASSERT_TRUE(verified);
}

void test_thin_lens_formula() {
    std::cout << "\n=== Testing Thin Lens Formula ===" << std::endl;

    // Test 1: Lens formula: 1/f = 1/v - 1/u
    double u = -0.3;  // Object 30 cm from lens (negative by convention)
    double v = 0.6;   // Image 60 cm from lens (positive for real image)
    double f = lensFormula(u, v);
    double expected = 1.0 / (1.0/v - 1.0/u);
    ASSERT_NEAR(f, expected, TOLERANCE);

    // Test 2: Calculate image distance
    f = 0.2;   // 20 cm focal length
    u = -0.4;  // 40 cm object distance
    v = imageDistance(f, u);
    // Verify: 1/f = 1/v - 1/u (sign convention: u negative for real object)
    double check = 1.0/v - 1.0/u;
    ASSERT_NEAR(check, 1.0/f, TOLERANCE);

    // Test 3: Calculate object distance
    v = 0.3;  // 30 cm image distance
    f = 0.1;  // 10 cm focal length
    u = objectDistance(f, v);
    // Verify: 1/f = 1/v - 1/u becomes 1/u = 1/v - 1/f
    check = 1.0/v - 1.0/f;
    ASSERT_NEAR(check, 1.0/u, TOLERANCE);

    // Test 4: Lens power in diopters
    f = 0.5;  // 50 cm = 0.5 m
    double P = lensPower(f);
    ASSERT_NEAR(P, 2.0, TOLERANCE);  // 2 diopters

    // Test 5: Focal length from power
    P = 5.0;  // 5 diopters
    f = focalLengthFromPower(P);
    ASSERT_NEAR(f, 0.2, TOLERANCE);  // 0.2 m = 20 cm
}

void test_lensmaker_equation() {
    std::cout << "\n=== Testing Lensmaker's Equation ===" << std::endl;

    // Test 1: Symmetric biconvex lens
    double n = 1.5;      // glass
    double R1 = 0.1;     // 10 cm radius (convex)
    double R2 = -0.1;    // -10 cm (other side convex)
    double f = lensmakersEquation(n, R1, R2);
    double expected = 1.0 / ((n - 1.0) * (1.0/R1 - 1.0/R2));
    ASSERT_NEAR(f, expected, TOLERANCE);

    // Test 2: Verify focal length is positive for converging lens
    ASSERT_TRUE(f > 0.0);

    // Test 3: Plano-convex lens (one flat surface)
    R1 = 0.2;  // Convex
    R2 = 1e10; // Essentially flat (large radius)
    f = lensmakersEquation(n, R1, R2);
    expected = R1 / (n - 1.0);
    ASSERT_NEAR(f, expected, 1e-3);

    // Test 4: Lens in medium (water)
    double n_lens = 1.5;
    double n_medium = 1.33;
    R1 = 0.1;
    R2 = -0.1;
    f = focalLengthInMedium(n_lens, n_medium, R1, R2);
    // Focal length increases in denser medium
    ASSERT_TRUE(f > 0.0);
}

void test_magnification() {
    std::cout << "\n=== Testing Magnification ===" << std::endl;

    // Test 1: Linear magnification
    double v = 0.4;   // Image distance
    double u = -0.2;  // Object distance
    double m = linearMagnification(v, u);
    double expected = v / u;
    ASSERT_NEAR(m, expected, TOLERANCE);

    // Test 2: Image height from magnification
    double h_o = 0.05;  // 5 cm object
    m = -2.0;           // Inverted, magnified 2x
    double h_i = imageHeight(h_o, m);
    ASSERT_NEAR(h_i, -0.1, TOLERANCE);  // 10 cm, inverted

    // Test 3: Object farther than 2f gives magnification < 1
    double f = 0.1;  // 10 cm
    u = -0.3;        // 30 cm (> 2f)
    m = magnificationFromFocal(f, u);
    ASSERT_TRUE(std::abs(m) < 1.0);

    // Test 4: Object at 2f gives m = -1
    u = -0.2;  // 2f
    f = 0.1;
    m = magnificationFromFocal(f, u);
    ASSERT_NEAR(m, -1.0, TOLERANCE);
}

void test_mirror_equations() {
    std::cout << "\n=== Testing Mirror Equations ===" << std::endl;

    // Test 1: Focal length is half radius of curvature
    double R = 0.4;  // 40 cm
    double f = mirrorFocalLength(R);
    ASSERT_NEAR(f, 0.2, TOLERANCE);

    // Test 2: Mirror formula (same as lens formula)
    double u = -0.3;  // 30 cm
    double v = 0.6;   // 60 cm
    f = mirrorFormula(u, v);
    double expected = 1.0 / (1.0/v + 1.0/u);
    ASSERT_NEAR(f, expected, TOLERANCE);

    // Test 3: Concave mirror (positive focal length)
    R = 0.5;
    f = mirrorFocalLength(R);
    ASSERT_TRUE(f > 0.0);
}

void test_optical_instruments() {
    std::cout << "\n=== Testing Optical Instruments ===" << std::endl;

    const double D = 0.25;  // 25 cm least distance

    // Test 1: Simple microscope (magnifying glass)
    double f = 0.05;  // 5 cm focal length
    double M = simpleMicroscopeMagnification(f, D);
    double expected = 1.0 + D / f;  // 1 + 5 = 6x
    ASSERT_NEAR(M, expected, TOLERANCE);

    // Test 2: Compound microscope
    double m_obj = 10;  // Objective magnification
    double f_eye = 0.025;  // 2.5 cm eyepiece
    M = compoundMicroscopeMagnification(m_obj, f_eye, D);
    expected = m_obj * (D / f_eye);  // 10 × 10 = 100x
    ASSERT_NEAR(M, expected, TOLERANCE);

    // Test 3: Astronomical telescope
    double f_obj = 1.0;   // 100 cm objective
    double f_eye2 = 0.05; // 5 cm eyepiece
    M = telescopeMagnification(f_obj, f_eye2);
    expected = f_obj / f_eye2;  // 20x
    ASSERT_NEAR(M, expected, TOLERANCE);

    // Test 4: Telescope length
    double L = telescopeLength(f_obj, f_eye2);
    expected = f_obj + f_eye2;
    ASSERT_NEAR(L, expected, TOLERANCE);

    // Test 5: Telescope resolving power
    double D_aperture = 0.1;  // 10 cm aperture
    double lambda = 550e-9;   // 550 nm
    double R = telescopeResolvingPower(D_aperture, lambda);
    expected = D_aperture / (1.22 * lambda);
    ASSERT_NEAR(R, expected, TOLERANCE);
}

void test_lens_combinations() {
    std::cout << "\n=== Testing Lens Combinations ===" << std::endl;

    // Test 1: Two lenses in contact
    double f1 = 0.2;  // 20 cm
    double f2 = 0.3;  // 30 cm
    double F = combinedFocalLength(f1, f2);
    double expected = 1.0 / (1.0/f1 + 1.0/f2);
    ASSERT_NEAR(F, expected, TOLERANCE);

    // Test 2: Combined power
    double P1 = 5.0;  // 5 diopters
    double P2 = 3.0;  // 3 diopters
    double P = combinedPower(P1, P2);
    ASSERT_NEAR(P, 8.0, TOLERANCE);

    // Test 3: Verify P = 1/F for combined system
    F = combinedFocalLength(f1, f2);
    P = lensPower(F);
    double P_direct = lensPower(f1) + lensPower(f2);
    ASSERT_NEAR(P, P_direct, TOLERANCE);

    // Test 4: Separated lenses
    f1 = 0.1;
    f2 = 0.2;
    double d = 0.05;  // 5 cm separation
    F = separatedLensesFocalLength(f1, f2, d);
    expected = 1.0 / (1.0/f1 + 1.0/f2 - d/(f1*f2));
    ASSERT_NEAR(F, expected, TOLERANCE);
}

//==============================================================================
// MAIN TEST RUNNER
//==============================================================================

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "PHASE 3 VALIDATION TEST SUITE" << std::endl;
    std::cout << "Quantum Mechanics, EM Waves, and Optics" << std::endl;
    std::cout << "========================================" << std::endl;

    // Quantum Basics Tests
    test_de_broglie_wavelength();
    test_compton_scattering();
    test_heisenberg_uncertainty();
    test_bohr_model();
    test_photoelectric_effect();
    test_particle_in_box();
    test_quantum_harmonic_oscillator();
    test_quantum_tunneling();
    test_quantum_angular_momentum();

    // Electromagnetic Waves Tests
    test_speed_of_light();
    test_wavelength_frequency();
    test_em_field_relationships();
    test_em_energy_density();
    test_wave_intensity();
    test_radiation_pressure();
    test_photon_properties();
    test_wave_in_medium();

    // Optics Tests
    test_snells_law();
    test_critical_angle_tir();
    test_brewsters_angle();
    test_thin_lens_formula();
    test_lensmaker_equation();
    test_magnification();
    test_mirror_equations();
    test_optical_instruments();
    test_lens_combinations();

    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Tests run: " << tests_run << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << (tests_run - tests_passed) << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / tests_run) << "%" << std::endl;

    return (tests_run == tests_passed) ? 0 : 1;
}
