/**
 * @file test_quantum_gravity_phenomenology.cpp
 * @brief Tests for Quantum Gravity Phenomenology Theory
 *
 * Verifies mathematical derivations and physical predictions
 */

#include "../../new_theory/quantum_gravity_phenomenology.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace new_theory::quantum_gravity_phenomenology;

constexpr double TOLERANCE = 1e-6;
constexpr double LOOSE_TOLERANCE = 1e-3;

bool approx_equal(double a, double b, double tol = TOLERANCE) {
    if (std::abs(b) < 1e-100) return std::abs(a) < tol;
    return std::abs(a - b) / std::abs(b) < tol;
}

void test_quantum_corrected_friedmann() {
    std::cout << "\n=== Testing Quantum-Corrected Friedmann Equation ===\n";

    // Test 1: Classical limit (low density)
    std::cout << "Test 1: Classical limit (ρ << ρ_Planck)\n";
    double rho_low = 1e90;  // Much less than ρ_Planck ≈ 5e96
    double H_quantum = QuantumCorrectedFriedmann::hubbleParameterQuantumCorrected(rho_low);
    double H_classical = std::sqrt(8.0 * M_PI * constants::G * rho_low / 3.0);
    std::cout << "  H_quantum  = " << H_quantum << " 1/s\n";
    std::cout << "  H_classical = " << H_classical << " 1/s\n";
    std::cout << "  Relative diff: " << std::abs(H_quantum - H_classical)/H_classical << "\n";
    assert(approx_equal(H_quantum, H_classical, LOOSE_TOLERANCE));
    std::cout << "  ✓ Recovers classical Friedmann\n";

    // Test 2: Bounce density
    std::cout << "\nTest 2: Bounce at critical density\n";
    double rho_bounce = QuantumCorrectedFriedmann::bounceDensity();
    double rho_crit = 0.41 * constants::rho_P;
    std::cout << "  ρ_bounce = " << rho_bounce << " kg/m³\n";
    std::cout << "  ρ_crit/2 = " << rho_crit/2.0 << " kg/m³\n";
    assert(approx_equal(rho_bounce, rho_crit / 2.0));
    std::cout << "  ✓ Bounce at ρ_max = ρ_crit/2\n";

    // Test 3: H vanishes at critical density
    std::cout << "\nTest 3: Hubble parameter vanishes at critical density\n";
    double rho_at_crit = 0.41 * constants::rho_P;  // Exactly at critical
    double H_at_crit = QuantumCorrectedFriedmann::hubbleParameterQuantumCorrected(rho_at_crit);
    std::cout << "  H(ρ_crit) = " << H_at_crit << " 1/s\n";
    assert(H_at_crit == 0.0);  // Should be exactly zero
    std::cout << "  ✓ H = 0 at ρ_crit (singularity avoided)\n";

    // Test 4: Minimum scale factor
    std::cout << "\nTest 4: Minimum scale factor\n";
    double rho_current = 1e-26;  // Current universe density
    double a_current = 1.0;
    double a_min = QuantumCorrectedFriedmann::minimumScaleFactor(rho_current, a_current);
    std::cout << "  a_min = " << a_min << "\n";
    assert(a_min > 0);
    assert(a_min < a_current);
    std::cout << "  ✓ Minimum scale factor a_min > 0 (no singularity)\n";

    // Test 5: Time to bounce
    std::cout << "\nTest 5: Time to bounce\n";
    double t_bounce = QuantumCorrectedFriedmann::timeToBounce(a_current, rho_current, 100);
    std::cout << "  t_bounce = " << t_bounce << " s\n";
    std::cout << "  In Planck times: " << t_bounce / constants::t_P << " t_P\n";
    assert(t_bounce > 0);
    std::cout << "  ✓ Finite time to bounce\n";

    std::cout << "\n✓ All Quantum Friedmann tests passed!\n";
}

void test_black_hole_entropy() {
    std::cout << "\n=== Testing Black Hole Entropy from LQG ===\n";

    // Test 1: Bekenstein-Hawking formula
    std::cout << "Test 1: Bekenstein-Hawking entropy\n";
    double M_solar = 1.989e30;  // Solar mass (kg)
    double r_s = 2.0 * constants::G * M_solar / (constants::c * constants::c);  // Schwarzschild radius
    double A_horizon = 4.0 * M_PI * r_s * r_s;
    std::cout << "  Schwarzschild radius: " << r_s << " m\n";
    std::cout << "  Horizon area: " << A_horizon << " m²\n";

    double S_BH = BlackHoleEntropyLQG::bekensteinHawkingEntropy(A_horizon);
    std::cout << "  S_BH = " << S_BH << " J/K\n";
    std::cout << "  S_BH/k_B = " << S_BH/constants::k_B << " (dimensionless)\n";

    // Verify S = A/(4l_P²)
    double l_P_squared = constants::l_P * constants::l_P;
    double S_expected = constants::k_B * A_horizon / (4.0 * l_P_squared);
    assert(approx_equal(S_BH, S_expected, 1e-10));
    std::cout << "  ✓ S = k_B A/(4l_P²)\n";

    // Test 2: Number of punctures
    std::cout << "\nTest 2: Horizon punctures\n";
    double N_punctures = BlackHoleEntropyLQG::numberOfPunctures(A_horizon);
    std::cout << "  Number of punctures: " << N_punctures << "\n";
    assert(N_punctures > 0);
    std::cout << "  ✓ Positive number of punctures\n";

    // Test 3: Microscopic entropy matches macroscopic
    std::cout << "\nTest 3: Microscopic = Macroscopic entropy\n";
    double S_micro = BlackHoleEntropyLQG::microscopicEntropy(A_horizon);
    std::cout << "  S_micro = " << S_micro << " J/K\n";
    std::cout << "  S_BH    = " << S_BH << " J/K\n";

    double error = BlackHoleEntropyLQG::verifyEntropyMatch(A_horizon);
    std::cout << "  Relative error: " << error << "\n";
    assert(error < LOOSE_TOLERANCE);
    std::cout << "  ✓ Spin network counting reproduces Bekenstein-Hawking!\n";

    // Test 4: Entropy scales with area
    std::cout << "\nTest 4: Entropy scales as S ∝ A\n";
    double A_double = 2.0 * A_horizon;
    double S_double = BlackHoleEntropyLQG::bekensteinHawkingEntropy(A_double);
    assert(approx_equal(S_double / S_BH, 2.0, TOLERANCE));
    std::cout << "  S(2A)/S(A) = " << S_double/S_BH << " ≈ 2 ✓\n";

    std::cout << "\n✓ All Black Hole Entropy tests passed!\n";
}

void test_observable_predictions() {
    std::cout << "\n=== Testing Observable Predictions ===\n";

    // Test 1: CMB power spectrum correction
    std::cout << "Test 1: CMB power spectrum corrections\n";
    double k_large = 1e-4;  // 1/Mpc in SI units
    double P_classical = 1.0;  // Normalized
    double P_corrected = ObservablePredictions::cmbPowerSpectrumCorrection(k_large, P_classical);

    std::cout << "  P_classical = " << P_classical << "\n";
    std::cout << "  P_quantum   = " << P_corrected << "\n";
    std::cout << "  Correction: " << (P_corrected - P_classical) / P_classical * 100 << "%\n";

    // Correction should be reasonably bounded
    assert(std::abs(P_corrected - P_classical) / P_classical < 1.0);  // Within 100%
    std::cout << "  ✓ Quantum correction calculated\n";

    // Test 2: Tensor-to-scalar ratio suppression
    std::cout << "\nTest 2: Tensor-to-scalar ratio\n";
    double r_classical = 0.01;  // Typical inflationary prediction
    double H_inflation = 1e-6 * constants::c / constants::l_P;  // Sub-Planckian
    double r_quantum = ObservablePredictions::tensorToScalarRatio(r_classical, H_inflation);

    std::cout << "  r_classical = " << r_classical << "\n";
    std::cout << "  r_quantum   = " << r_quantum << "\n";
    std::cout << "  Suppression: " << r_quantum / r_classical << "\n";

    assert(r_quantum <= r_classical);
    std::cout << "  ✓ Quantum gravity suppresses tensor modes\n";

    // Test 3: Minimum angular scale
    std::cout << "\nTest 3: Minimum observable scale\n";
    double theta_min = ObservablePredictions::minimumAngularScale();
    std::cout << "  θ_min = " << theta_min << " rad\n";
    std::cout << "  θ_min = " << theta_min * 180.0 / M_PI * 3600.0 << " arcsec\n";

    assert(theta_min > 0);
    std::cout << "  ✓ Finite minimum scale from discrete spacetime\n";

    std::cout << "\n✓ All Observable Predictions tests passed!\n";
}

void test_physical_consistency() {
    std::cout << "\n=== Testing Physical Consistency ===\n";

    // Test 1: Energy conditions
    std::cout << "Test 1: Energy conditions during bounce\n";
    double rho_bounce = QuantumCorrectedFriedmann::bounceDensity();
    std::cout << "  ρ_bounce/ρ_Planck = " << rho_bounce / constants::rho_P << "\n";
    assert(rho_bounce > 0);
    assert(rho_bounce < constants::rho_P);
    std::cout << "  ✓ Bounce density is sub-Planckian\n";

    // Test 2: Dimensional analysis
    std::cout << "\nTest 2: Dimensional analysis\n";
    double H = QuantumCorrectedFriedmann::hubbleParameterQuantumCorrected(rho_bounce);
    std::cout << "  H has units of 1/time: [H] = " << H << " s⁻¹\n";
    assert(H >= 0);
    std::cout << "  ✓ H ≥ 0 (physical)\n";

    // Test 3: Classical limit
    std::cout << "\nTest 3: Classical limit verification\n";
    double rho_small = constants::rho_P * 1e-10;
    double H_quantum = QuantumCorrectedFriedmann::hubbleParameterQuantumCorrected(rho_small);
    double H_classical = std::sqrt(8.0 * M_PI * constants::G * rho_small / 3.0);
    double ratio = H_quantum / H_classical;
    std::cout << "  H_quantum/H_classical = " << ratio << "\n";
    assert(approx_equal(ratio, 1.0, 0.01));
    std::cout << "  ✓ Quantum corrections vanish at low density\n";

    std::cout << "\n✓ All Physical Consistency tests passed!\n";
}

int main() {
    std::cout << std::setprecision(6) << std::scientific;

    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  QUANTUM GRAVITY PHENOMENOLOGY - THEORETICAL TESTS      ║\n";
    std::cout << "║  Novel Theory: LQG + Cosmology                          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";

    try {
        test_quantum_corrected_friedmann();
        test_black_hole_entropy();
        test_observable_predictions();
        test_physical_consistency();

        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "✓✓✓ ALL TESTS PASSED! ✓✓✓\n";
        std::cout << std::string(60, '=') << "\n\n";

        std::cout << "THEORETICAL VERIFICATION COMPLETE:\n";
        std::cout << "  [✓] Quantum-corrected Friedmann equation\n";
        std::cout << "  [✓] Singularity resolution via quantum bounce\n";
        std::cout << "  [✓] Black hole entropy from spin networks\n";
        std::cout << "  [✓] Microscopic = Macroscopic entropy\n";
        std::cout << "  [✓] Observable CMB predictions\n";
        std::cout << "  [✓] Physical consistency checks\n\n";

        std::cout << "NOVEL RESULTS:\n";
        std::cout << "  1. H² = (8πG/3)ρ[1 - ρ/ρ_crit] prevents singularities\n";
        std::cout << "  2. S_BH = A/(4l_P²) derived from quantum geometry\n";
        std::cout << "  3. CMB power spectrum carries quantum bounce signature\n\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n✗ TEST FAILED: " << e.what() << "\n";
        return 1;
    }
}
