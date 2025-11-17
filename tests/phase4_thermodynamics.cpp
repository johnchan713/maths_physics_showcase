/**
 * Phase 4 Validation: Thermodynamics and Gas Laws
 *
 * Tests the thermodynamics.hpp module functions.
 *
 * Coverage:
 * - Boyle's Law P₁V₁ = P₂V₂ (constant T)
 * - Charles's Law V₁/T₁ = V₂/T₂ (constant P)
 * - Gay-Lussac's Law P₁/T₁ = P₂/T₂ (constant V)
 * - Combined Gas Law (P₁V₁)/T₁ = (P₂V₂)/T₂
 * - Ideal Gas Law PV = nRT
 * - Kinetic theory and RMS velocity
 * - Pressure from molecular motion
 * - Bulk modulus isothermal and adiabatic
 * - Temperature conversions (Celsius Kelvin Fahrenheit)
 * - Standard conditions and molar volume
 * - Gas density calculations
 */

#include <iostream>
#include <cmath>
#include "../include/physics/thermodynamics.hpp"

using namespace physics::thermodynamics;

// Test tolerances
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

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    std::cout << "=== Phase 4: Thermodynamics Validation ===" << std::endl;
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
    // Boyle's Law
    // ========================================

    run_test("Boyle's Law formula", []() {
        double P1 = 100000.0;  // 100 kPa
        double V1 = 2.0;       // 2 m³
        double V2 = 1.0;       // 1 m³

        double P2 = boylesLaw(P1, V1, V2);

        // P₁V₁ = P₂V₂, so P₂ = P₁V₁/V₂ = 200 kPa
        ASSERT_NEAR(P2, 200000.0, TOLERANCE);
        return true;
    });

    run_test("Boyle's Law inverse relationship", []() {
        double P1 = 100000.0;
        double V1 = 1.0;
        double V2 = 0.5;  // Halve volume

        double P2 = boylesLaw(P1, V1, V2);

        // Pressure should double
        ASSERT_NEAR(P2, 2.0 * P1, TOLERANCE);
        return true;
    });

    run_test("Boyle's Law compression", []() {
        double P1 = 101325.0;  // 1 atm
        double V1 = 10.0;
        double V2 = 2.0;  // Compress to 1/5 volume

        double P2 = boylesLaw(P1, V1, V2);

        ASSERT_NEAR(P2, 5.0 * P1, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Charles's Law
    // ========================================

    run_test("Charles's Law formula", []() {
        double V1 = 1.0;    // 1 m³
        double T1 = 300.0;  // 300 K
        double T2 = 600.0;  // 600 K

        double V2 = charlesLaw(V1, T1, T2);

        // V₁/T₁ = V₂/T₂, so V₂ = V₁T₂/T₁ = 2.0
        ASSERT_NEAR(V2, 2.0, TOLERANCE);
        return true;
    });

    run_test("Charles's Law heating", []() {
        double V1 = 5.0;
        double T1 = 273.15;  // 0°C
        double T2 = 373.15;  // 100°C

        double V2 = charlesLaw(V1, T1, T2);

        ASSERT_NEAR(V2 / V1, T2 / T1, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Charles's Law cooling", []() {
        double V1 = 10.0;
        double T1 = 400.0;
        double T2 = 200.0;  // Halve temperature

        double V2 = charlesLaw(V1, T1, T2);

        ASSERT_NEAR(V2, 5.0, TOLERANCE);  // Volume halves
        return true;
    });

    // ========================================
    // Gay-Lussac's Law
    // ========================================

    run_test("Gay-Lussac's Law formula", []() {
        double P1 = 100000.0;  // 100 kPa
        double T1 = 300.0;     // 300 K
        double T2 = 600.0;     // 600 K

        double P2 = gayLussacsLaw(P1, T1, T2);

        // P₁/T₁ = P₂/T₂, so P₂ = P₁T₂/T₁ = 200 kPa
        ASSERT_NEAR(P2, 200000.0, TOLERANCE);
        return true;
    });

    run_test("Gay-Lussac's Law heating", []() {
        double P1 = 101325.0;
        double T1 = 273.15;
        double T2 = 373.15;

        double P2 = gayLussacsLaw(P1, T1, T2);

        ASSERT_NEAR(P2 / P1, T2 / T1, LOOSE_TOLERANCE);
        return true;
    });

    run_test("Gay-Lussac's Law temperature doubling", []() {
        double P1 = 50000.0;
        double T1 = 200.0;
        double T2 = 400.0;

        double P2 = gayLussacsLaw(P1, T1, T2);

        ASSERT_NEAR(P2, 2.0 * P1, TOLERANCE);
        return true;
    });

    // ========================================
    // Combined Gas Law
    // ========================================

    run_test("Combined Gas Law formula", []() {
        double P1 = 100000.0;
        double V1 = 2.0;
        double T1 = 300.0;
        double V2 = 1.0;
        double T2 = 600.0;

        double P2 = combinedGasLaw(P1, V1, T1, V2, T2);

        // (P₁V₁)/T₁ = (P₂V₂)/T₂
        // P₂ = P₁(V₁/V₂)(T₂/T₁) = 100000 × 2 × 2 = 400000
        ASSERT_NEAR(P2, 400000.0, TOLERANCE);
        return true;
    });

    run_test("Combined Gas Law reduces to Boyle", []() {
        double P1 = 100000.0;
        double V1 = 2.0;
        double T = 300.0;
        double V2 = 1.0;

        double P2_combined = combinedGasLaw(P1, V1, T, V2, T);
        double P2_boyle = boylesLaw(P1, V1, V2);

        ASSERT_NEAR(P2_combined, P2_boyle, TOLERANCE);
        return true;
    });

    run_test("Combined Gas Law reduces to Charles", []() {
        double P = 100000.0;
        double V1 = 1.0;
        double T1 = 300.0;
        double T2 = 600.0;

        // At constant P, if we start with V1 at T1, we can calculate P2 at V2 and T2
        // Then check P2 = P1 (constant pressure condition)
        double V2_charles = charlesLaw(V1, T1, T2);
        double P2 = combinedGasLaw(P, V1, T1, V2_charles, T2);

        ASSERT_NEAR(P2, P, LOOSE_TOLERANCE);
        return true;
    });

    // ========================================
    // Ideal Gas Law - Pressure
    // ========================================

    run_test("Ideal Gas Law pressure calculation", []() {
        double n = 1.0;      // 1 mole
        double V = 0.0224;   // 22.4 liters = 0.0224 m³
        double T = 273.15;   // 0°C

        double P = idealGasLawPressure(n, V, T);

        // At STP: P ≈ 101325 Pa
        ASSERT_NEAR(P, constants::STP_PRESSURE, 1000.0);
        return true;
    });

    run_test("Ideal Gas Law PV = nRT", []() {
        double n = 2.0;
        double V = 1.0;
        double T = 300.0;

        double P = idealGasLawPressure(n, V, T);

        double expected = (n * constants::R * T) / V;
        ASSERT_NEAR(P, expected, TOLERANCE);
        return true;
    });

    run_test("Pressure doubles with moles", []() {
        double n1 = 1.0;
        double n2 = 2.0;
        double V = 1.0;
        double T = 300.0;

        double P1 = idealGasLawPressure(n1, V, T);
        double P2 = idealGasLawPressure(n2, V, T);

        ASSERT_NEAR(P2, 2.0 * P1, TOLERANCE);
        return true;
    });

    // ========================================
    // Ideal Gas Law - Volume
    // ========================================

    run_test("Ideal Gas Law volume calculation", []() {
        double n = 1.0;
        double T = 273.15;
        double P = 101325.0;

        double V = idealGasLawVolume(n, T, P);

        // At STP: V ≈ 0.0224 m³ = 22.4 liters
        ASSERT_NEAR(V, 0.0224, 0.001);
        return true;
    });

    run_test("Volume from ideal gas law", []() {
        double n = 2.0;
        double T = 300.0;
        double P = 100000.0;

        double V = idealGasLawVolume(n, T, P);

        double expected = (n * constants::R * T) / P;
        ASSERT_NEAR(V, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Ideal Gas Law - Temperature
    // ========================================

    run_test("Ideal Gas Law temperature calculation", []() {
        double P = 101325.0;
        double V = 0.0224;
        double n = 1.0;

        double T = idealGasLawTemperature(P, V, n);

        // Should be STP temperature
        ASSERT_NEAR(T, constants::STP_TEMPERATURE, 1.0);
        return true;
    });

    run_test("Temperature from ideal gas law", []() {
        double P = 200000.0;
        double V = 1.0;
        double n = 2.0;

        double T = idealGasLawTemperature(P, V, n);

        double expected = (P * V) / (n * constants::R);
        ASSERT_NEAR(T, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Ideal Gas Law - Moles
    // ========================================

    run_test("Ideal Gas Law moles calculation", []() {
        double P = 101325.0;
        double V = 0.0224;
        double T = 273.15;

        double n = idealGasLawMoles(P, V, T);

        // Should be 1 mole
        ASSERT_NEAR(n, 1.0, 0.01);
        return true;
    });

    run_test("Moles from ideal gas law", []() {
        double P = 100000.0;
        double V = 2.0;
        double T = 300.0;

        double n = idealGasLawMoles(P, V, T);

        double expected = (P * V) / (constants::R * T);
        ASSERT_NEAR(n, expected, TOLERANCE);
        return true;
    });

    // ========================================
    // Kinetic Theory
    // ========================================

    run_test("Pressure from kinetic theory", []() {
        double n_over_V = 2.5e25;  // Number density (molecules/m³)
        double m = 5e-26;          // Molecular mass (kg)
        double v2 = 1e6;           // Mean square velocity (m²/s²)

        double P = pressureFromKineticTheory(n_over_V, m, v2);

        ASSERT_TRUE(P > 0);
        return true;
    });

    run_test("Pressure proportional to temperature", []() {
        double n_over_V = 2.5e25;
        double m = 5e-26;
        double v2_1 = 1e6;
        double v2_2 = 2e6;

        double P1 = pressureFromKineticTheory(n_over_V, m, v2_1);
        double P2 = pressureFromKineticTheory(n_over_V, m, v2_2);

        // P ∝ <v²>, so P₂/P₁ = 2
        ASSERT_NEAR(P2 / P1, 2.0, TOLERANCE);
        return true;
    });

    run_test("RMS velocity calculation", []() {
        double T = 300.0;   // 300 K
        double m = 5e-26;   // ~N₂ molecule mass

        double v_rms = rmsVelocity(T, m);

        // Should be ~500 m/s for nitrogen
        ASSERT_TRUE(v_rms > 400.0 && v_rms < 600.0);
        return true;
    });

    run_test("RMS velocity proportional to sqrt T", []() {
        double m = 5e-26;
        double T1 = 300.0;
        double T2 = 1200.0;  // 4 times hotter

        double v1 = rmsVelocity(T1, m);
        double v2 = rmsVelocity(T2, m);

        // v ∝ √T, so v₂/v₁ = 2
        ASSERT_NEAR(v2 / v1, 2.0, LOOSE_TOLERANCE);
        return true;
    });

    run_test("RMS velocity using molar mass", []() {
        double T = 300.0;
        double M = 0.028;  // N₂ molar mass (kg/mol)

        double v_rms = rmsVelocityMolar(T, M);

        ASSERT_TRUE(v_rms > 400.0 && v_rms < 600.0);
        return true;
    });

    // ========================================
    // Elasticity of Gases
    // ========================================

    run_test("Isothermal bulk modulus equals pressure", []() {
        double P = 100000.0;

        double K_T = isothermalBulkModulus(P);

        ASSERT_NEAR(K_T, P, TOLERANCE);
        return true;
    });

    run_test("Adiabatic bulk modulus gamma times pressure", []() {
        double P = 100000.0;
        double gamma = 1.4;  // Air

        double K_S = adiabaticBulkModulus(P, gamma);

        ASSERT_NEAR(K_S, gamma * P, TOLERANCE);
        return true;
    });

    run_test("Adiabatic modulus greater than isothermal", []() {
        double P = 100000.0;
        double gamma = 1.4;

        double K_T = isothermalBulkModulus(P);
        double K_S = adiabaticBulkModulus(P, gamma);

        ASSERT_TRUE(K_S > K_T);
        return true;
    });

    run_test("Isothermal compressibility inverse of pressure", []() {
        double P = 100000.0;

        double kappa = isothermalCompressibility(P);

        ASSERT_NEAR(kappa, 1.0 / P, TOLERANCE);
        return true;
    });

    // ========================================
    // Temperature Conversions
    // ========================================

    run_test("Celsius to Kelvin", []() {
        double T_K = celsiusToKelvin(0.0);

        ASSERT_NEAR(T_K, 273.15, TOLERANCE);
        return true;
    });

    run_test("Kelvin to Celsius", []() {
        double T_C = kelvinToCelsius(273.15);

        ASSERT_NEAR(T_C, 0.0, TOLERANCE);
        return true;
    });

    run_test("Temperature conversion round trip", []() {
        double T_C = 25.0;

        double T_K = celsiusToKelvin(T_C);
        double T_C2 = kelvinToCelsius(T_K);

        ASSERT_NEAR(T_C2, T_C, TOLERANCE);
        return true;
    });

    run_test("Fahrenheit to Kelvin", []() {
        double T_K = fahrenheitToKelvin(32.0);  // Freezing point

        ASSERT_NEAR(T_K, 273.15, TOLERANCE);
        return true;
    });

    run_test("Kelvin to Fahrenheit", []() {
        double T_F = kelvinToFahrenheit(273.15);

        ASSERT_NEAR(T_F, 32.0, TOLERANCE);
        return true;
    });

    run_test("Fahrenheit conversion round trip", []() {
        double T_F = 68.0;

        double T_K = fahrenheitToKelvin(T_F);
        double T_F2 = kelvinToFahrenheit(T_K);

        ASSERT_NEAR(T_F2, T_F, TOLERANCE);
        return true;
    });

    // ========================================
    // Standard Conditions
    // ========================================

    run_test("Molar volume at STP", []() {
        double V_m = molarVolumeAtSTP();

        // Should be 0.0224 m³/mol = 22.4 L/mol
        ASSERT_NEAR(V_m, 0.0224, 0.001);
        return true;
    });

    run_test("Gas density calculation", []() {
        double P = 101325.0;
        double M = 0.029;  // Air molar mass (kg/mol)
        double T = 273.15;

        double rho = gasDensity(P, M, T);

        // Air density at STP ≈ 1.29 kg/m³
        ASSERT_NEAR(rho, 1.29, 0.1);
        return true;
    });

    run_test("Density proportional to pressure", []() {
        double P1 = 100000.0;
        double P2 = 200000.0;
        double M = 0.029;
        double T = 300.0;

        double rho1 = gasDensity(P1, M, T);
        double rho2 = gasDensity(P2, M, T);

        ASSERT_NEAR(rho2 / rho1, 2.0, TOLERANCE);
        return true;
    });

    run_test("Density inversely proportional to temperature", []() {
        double P = 100000.0;
        double M = 0.029;
        double T1 = 300.0;
        double T2 = 600.0;

        double rho1 = gasDensity(P, M, T1);
        double rho2 = gasDensity(P, M, T2);

        ASSERT_NEAR(rho1 / rho2, 2.0, TOLERANCE);
        return true;
    });

    // ========================================
    // Summary
    // ========================================

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Test Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    std::cout << "Total tests:  " << (tests_passed + tests_failed) << std::endl;
    std::cout << "Success rate: " << (100.0 * tests_passed / (tests_passed + tests_failed)) << "%" << std::endl;

    return tests_failed == 0 ? 0 : 1;
}
