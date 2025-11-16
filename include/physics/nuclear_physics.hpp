#ifndef PHYSICS_NUCLEAR_PHYSICS_HPP
#define PHYSICS_NUCLEAR_PHYSICS_HPP

#include <vector>
#include <complex>
#include <functional>
#include <cmath>
#include <array>
#include <string>
#include <map>

namespace physics {
namespace nuclear {

using Complex = std::complex<double>;

/**
 * @brief Physical constants for nuclear physics
 */
namespace constants {
    constexpr double c = 299792458.0;              // Speed of light (m/s)
    constexpr double hbar = 1.054571817e-34;       // Reduced Planck constant (J·s)
    constexpr double m_e = 9.1093837015e-31;       // Electron mass (kg)
    constexpr double m_p = 1.67262192369e-27;      // Proton mass (kg)
    constexpr double m_n = 1.67492749804e-27;      // Neutron mass (kg)
    constexpr double e = 1.602176634e-19;          // Elementary charge (C)
    constexpr double u = 1.66053906660e-27;        // Atomic mass unit (kg)
    constexpr double N_A = 6.02214076e23;          // Avogadro's number (1/mol)
    constexpr double alpha = 7.2973525693e-3;      // Fine structure constant
    constexpr double MeV = 1.602176634e-13;        // MeV in Joules
}

/**
 * @brief Stability of Nuclei
 *
 * Nuclear binding energy and stability criteria
 */
class NuclearStability {
public:
    /**
     * @brief Binding energy per nucleon
     *
     * BE/A peaks at Fe-56 (~8.8 MeV/nucleon)
     */
    static std::string binding_energy_curve() {
        return "BE/A: increases to Fe-56 (8.8 MeV/A), then decreases";
    }

    /**
     * @brief Semi-empirical mass formula (SEMF)
     *
     * B(A,Z) = a_v A - a_s A^(2/3) - a_c Z²/A^(1/3) - a_a (A-2Z)²/A ± δ(A,Z)
     */
    static double binding_energy_semf(int A, int Z) {
        // Coefficients (MeV)
        double a_v = 15.75;  // Volume term
        double a_s = 17.8;   // Surface term
        double a_c = 0.711;  // Coulomb term
        double a_a = 23.7;   // Asymmetry term
        double a_p = 11.18;  // Pairing term

        double N = A - Z;

        // Pairing term δ
        double delta = 0.0;
        if (Z % 2 == 0 && N % 2 == 0) {
            delta = a_p / std::sqrt(A);  // Even-even: +
        } else if (Z % 2 == 1 && N % 2 == 1) {
            delta = -a_p / std::sqrt(A);  // Odd-odd: -
        }
        // Odd-A: delta = 0

        double B = a_v * A
                 - a_s * std::pow(A, 2.0/3.0)
                 - a_c * Z * Z / std::pow(A, 1.0/3.0)
                 - a_a * (A - 2*Z) * (A - 2*Z) / A
                 + delta;

        return B;  // MeV
    }

    /**
     * @brief Binding energy per nucleon
     *
     * BE/A = B/A
     */
    static double binding_energy_per_nucleon(int A, int Z) {
        return binding_energy_semf(A, Z) / A;
    }

    /**
     * @brief Mass excess
     *
     * Δ = M - A (atomic mass - mass number)
     */
    static double mass_excess(int A, int Z) {
        double B = binding_energy_semf(A, Z) * constants::MeV;  // Convert to Joules
        // M = Z*m_p + N*m_n - B/c²
        int N = A - Z;
        double M = Z * constants::m_p + N * constants::m_n - B / (constants::c * constants::c);
        return (M - A * constants::u) / constants::u;  // In atomic mass units
    }

    /**
     * @brief Separation energy (neutron)
     *
     * S_n = B(A,Z) - B(A-1,Z)
     */
    static double neutron_separation_energy(int A, int Z) {
        return binding_energy_semf(A, Z) - binding_energy_semf(A - 1, Z);
    }

    /**
     * @brief Separation energy (proton)
     *
     * S_p = B(A,Z) - B(A-1,Z-1)
     */
    static double proton_separation_energy(int A, int Z) {
        return binding_energy_semf(A, Z) - binding_energy_semf(A - 1, Z - 1);
    }

    /**
     * @brief Q-value for nuclear reaction
     *
     * Q = (m_initial - m_final)c²
     */
    static double q_value(double mass_initial, double mass_final) {
        return (mass_initial - mass_final) * constants::c * constants::c / constants::MeV;  // MeV
    }

    /**
     * @brief Valley of stability
     *
     * Stable nuclei follow N ≈ Z for light, N > Z for heavy
     */
    static bool is_in_valley_of_stability(int A, int Z) {
        int N = A - Z;
        if (A < 40) {
            return std::abs(N - Z) <= 2;  // N ≈ Z for light nuclei
        } else {
            double expected_N = Z * (1.0 + 0.015 * std::pow(A, 2.0/3.0));
            return std::abs(N - expected_N) < 5;
        }
    }
};

/**
 * @brief Natural Radioactivity
 *
 * Types and properties of natural radioactive decay
 */
class NaturalRadioactivity {
public:
    /**
     * @brief Decay series
     *
     * Four natural decay series
     */
    static std::string decay_series() {
        return "4n (Th-232), 4n+1 (Np-237), 4n+2 (U-238), 4n+3 (U-235)";
    }

    /**
     * @brief Primordial radionuclides
     *
     * Long-lived isotopes from nucleosynthesis
     */
    static std::vector<std::string> primordial_radionuclides() {
        return {"U-238", "U-235", "Th-232", "K-40", "Rb-87"};
    }

    /**
     * @brief Cosmogenic radionuclides
     *
     * Produced by cosmic rays
     */
    static std::vector<std::string> cosmogenic_radionuclides() {
        return {"C-14", "Be-10", "Be-7", "H-3"};
    }
};

/**
 * @brief Alpha Decay (α)
 *
 * Emission of He-4 nucleus
 */
class AlphaDecay {
public:
    /**
     * @brief Alpha decay equation
     *
     * ᴬ_Z X → ᴬ⁻⁴_{Z-2} Y + ⁴_2 He
     */
    static std::string decay_equation() {
        return "ᴬ_Z X → ᴬ⁻⁴_{Z-2} Y + α (⁴He)";
    }

    /**
     * @brief Q-value for alpha decay
     *
     * Q_α = [M(A,Z) - M(A-4,Z-2) - M(He-4)]c²
     */
    static double q_value_alpha(double M_parent, double M_daughter) {
        double M_alpha = 4.002603 * constants::u;  // He-4 mass
        return (M_parent - M_daughter - M_alpha) * constants::c * constants::c / constants::MeV;
    }

    /**
     * @brief Kinetic energy of alpha particle
     *
     * T_α ≈ Q(A-4)/A (recoil approximation)
     */
    static double alpha_kinetic_energy(double Q_alpha, int A) {
        return Q_alpha * (A - 4.0) / A;  // MeV
    }

    /**
     * @brief Geiger-Nuttall law
     *
     * log₁₀(λ) = a + b/√E_α (empirical relation)
     */
    static double geiger_nuttall_constant(double E_alpha, double lambda) {
        return std::log10(lambda) * std::sqrt(E_alpha);
    }

    /**
     * @brief Gamow factor
     *
     * G = 2π(Z-2)e²/(ℏv) (barrier penetration)
     */
    static double gamow_factor(int Z, double v) {
        return 2.0 * M_PI * (Z - 2) * constants::e * constants::e / (constants::hbar * v);
    }

    /**
     * @brief Decay constant from Gamow theory
     *
     * λ ≈ ν₀ exp(-2G) where ν₀ ~ 10²¹ s⁻¹
     */
    static double decay_constant_gamow(double G) {
        double nu_0 = 1.0e21;  // Attempt frequency (s⁻¹)
        return nu_0 * std::exp(-2.0 * G);
    }
};

/**
 * @brief Beta Decay (β)
 *
 * β⁻, β⁺, and electron capture
 */
class BetaDecay {
public:
    /**
     * @brief Beta-minus decay
     *
     * n → p + e⁻ + ν̄_e
     */
    static std::string beta_minus() {
        return "ᴬ_Z X → ᴬ_{Z+1} Y + e⁻ + ν̄_e (n → p)";
    }

    /**
     * @brief Beta-plus decay
     *
     * p → n + e⁺ + ν_e
     */
    static std::string beta_plus() {
        return "ᴬ_Z X → ᴬ_{Z-1} Y + e⁺ + ν_e (p → n)";
    }

    /**
     * @brief Q-value for β⁻ decay
     *
     * Q_β⁻ = [M(A,Z) - M(A,Z+1)]c²
     */
    static double q_value_beta_minus(double M_parent, double M_daughter) {
        return (M_parent - M_daughter) * constants::c * constants::c / constants::MeV;
    }

    /**
     * @brief Q-value for β⁺ decay
     *
     * Q_β⁺ = [M(A,Z) - M(A,Z-1) - 2m_e]c²
     */
    static double q_value_beta_plus(double M_parent, double M_daughter) {
        double two_me = 2.0 * constants::m_e;
        return (M_parent - M_daughter - two_me) * constants::c * constants::c / constants::MeV;
    }

    /**
     * @brief Fermi theory decay constant
     *
     * λ = (G_F²/2π³ℏ⁷c⁶) |M_fi|² f(Z,Q)
     */
    static double fermi_decay_constant(double M_fi, double f_ZQ) {
        double G_F = 1.1663787e-5;  // GeV⁻² (Fermi constant)
        double factor = (G_F * G_F) / (2.0 * M_PI * M_PI * M_PI);
        return factor * M_fi * M_fi * f_ZQ;  // Approximate
    }

    /**
     * @brief Fermi integral
     *
     * f(Z,Q) for allowed transitions (approximate)
     */
    static double fermi_integral(int Z, double Q) {
        // Simplified approximation
        return std::pow(Q, 5.0) / 60.0;  // Non-relativistic approximation
    }

    /**
     * @brief Maximum beta energy
     *
     * E_max = Q_β (endpoint energy)
     */
    static double max_beta_energy(double Q_beta) {
        return Q_beta;  // MeV
    }

    /**
     * @brief Average beta energy
     *
     * ⟨E_β⟩ ≈ Q/3 (rule of thumb)
     */
    static double average_beta_energy(double Q_beta) {
        return Q_beta / 3.0;  // MeV
    }
};

/**
 * @brief Electron Capture (EC)
 *
 * p + e⁻ → n + ν_e
 */
class ElectronCapture {
public:
    /**
     * @brief EC decay equation
     *
     * ᴬ_Z X + e⁻ → ᴬ_{Z-1} Y + ν_e
     */
    static std::string decay_equation() {
        return "ᴬ_Z X + e⁻ → ᴬ_{Z-1} Y + ν_e (K-capture)";
    }

    /**
     * @brief Q-value for EC
     *
     * Q_EC = [M(A,Z) - M(A,Z-1)]c² - B_K
     */
    static double q_value_ec(double M_parent, double M_daughter, double B_K = 0.0) {
        return (M_parent - M_daughter) * constants::c * constants::c / constants::MeV - B_K;
    }

    /**
     * @brief K-shell binding energy
     *
     * B_K ≈ 13.6 Z² eV (hydrogen-like approximation)
     */
    static double k_shell_binding(int Z) {
        return 13.6 * Z * Z / 1000.0;  // keV
    }

    /**
     * @brief Competition β⁺ vs EC
     *
     * β⁺ requires Q > 2m_e c² = 1.022 MeV
     */
    static std::string beta_plus_threshold() {
        return "β⁺ requires Q > 1.022 MeV; EC possible for lower Q";
    }

    /**
     * @brief EC/β⁺ branching ratio
     *
     * Ratio depends on Q-value and Z
     */
    static double ec_branching_ratio(int Z, double Q) {
        // Simplified: EC favored for low Q and high Z
        if (Q < 1.022) return 1.0;  // Only EC possible
        double ratio = std::pow(Z / 30.0, 2.0) / (Q - 1.022);
        return ratio / (1.0 + ratio);
    }
};

/**
 * @brief Gamma Emission (γ)
 *
 * Electromagnetic transition between nuclear states
 */
class GammaEmission {
public:
    /**
     * @brief Gamma decay
     *
     * ᴬ_Z X* → ᴬ_Z X + γ
     */
    static std::string decay_equation() {
        return "ᴬ_Z X* → ᴬ_Z X + γ (E_γ = ΔE)";
    }

    /**
     * @brief Gamma energy
     *
     * E_γ = E_initial - E_final - E_recoil
     */
    static double gamma_energy(double E_initial, double E_final, double E_recoil = 0.0) {
        return E_initial - E_final - E_recoil;  // MeV
    }

    /**
     * @brief Nuclear recoil energy
     *
     * E_R = E_γ²/(2Mc²)
     */
    static double recoil_energy(double E_gamma, int A) {
        double M = A * constants::u;
        return (E_gamma * constants::MeV) * (E_gamma * constants::MeV) /
               (2.0 * M * constants::c * constants::c) / constants::MeV;  // MeV
    }

    /**
     * @brief Weisskopf estimate for E1 transition
     *
     * T_1/2(E1) ≈ 6.8×10⁻¹⁵ A^(-2/3) E_γ^(-3) seconds
     */
    static double weisskopf_e1_halflife(int A, double E_gamma) {
        return 6.8e-15 * std::pow(A, -2.0/3.0) * std::pow(E_gamma, -3.0);  // seconds
    }

    /**
     * @brief Weisskopf estimate for M1 transition
     *
     * T_1/2(M1) ≈ 2.2×10⁻¹⁴ E_γ^(-3) seconds
     */
    static double weisskopf_m1_halflife(double E_gamma) {
        return 2.2e-14 * std::pow(E_gamma, -3.0);  // seconds
    }

    /**
     * @brief Multipole transition rate
     *
     * λ(σL) ~ (E_γ/ℏc)^(2L+1)
     */
    static double transition_rate(double E_gamma, int L) {
        double k = E_gamma * constants::MeV / (constants::hbar * constants::c);
        return std::pow(k, 2.0 * L + 1.0);  // Relative units
    }

    /**
     * @brief Selection rules
     *
     * ΔJ, parity change
     */
    static std::string selection_rules() {
        return "E(L): ΔJ ≤ L, π_i π_f = (-1)^L; M(L): ΔJ ≤ L, π_i π_f = (-1)^(L+1)";
    }
};

/**
 * @brief Internal Conversion
 *
 * Nuclear excitation energy transferred to atomic electron
 */
class InternalConversion {
public:
    /**
     * @brief IC process
     *
     * ᴬ_Z X* + e⁻ → ᴬ_Z X + e⁻ (kinetic)
     */
    static std::string process() {
        return "ᴬ_Z X* + e⁻(bound) → ᴬ_Z X + e⁻(free) (IC competes with γ)";
    }

    /**
     * @brief Conversion coefficient
     *
     * α = λ_IC / λ_γ
     */
    static double conversion_coefficient(double lambda_ic, double lambda_gamma) {
        return lambda_ic / lambda_gamma;
    }

    /**
     * @brief IC electron energy
     *
     * E_e = E* - B_n (B_n = binding energy)
     */
    static double ic_electron_energy(double E_excitation, double B_binding) {
        return E_excitation - B_binding;  // MeV
    }

    /**
     * @brief K-shell conversion coefficient (estimate)
     *
     * α_K increases with Z and decreases with E_γ
     */
    static double alpha_k_estimate(int Z, double E_gamma, int L) {
        // Simplified empirical formula
        double factor = std::pow(Z / 100.0, 3.0) / std::pow(E_gamma, 3.0);
        return factor * std::pow(10.0, L);  // Rough estimate
    }

    /**
     * @brief Branching ratio
     *
     * BR_γ = 1/(1 + α_total)
     */
    static double gamma_branching_ratio(double alpha_total) {
        return 1.0 / (1.0 + alpha_total);
    }
};

/**
 * @brief Isomers and Isomeric Transition
 *
 * Metastable excited nuclear states
 */
class NuclearIsomers {
public:
    /**
     * @brief Isomer
     *
     * Metastable state with T_1/2 > 10⁻⁹ s
     */
    static std::string isomer_definition() {
        return "Isomer: excited state with T_1/2 > 1 ns (high spin trap)";
    }

    /**
     * @brief Isomeric transition (IT)
     *
     * ᴬ_Z Xᵐ → ᴬ_Z X + γ (or IC)
     */
    static std::string isomeric_transition() {
        return "IT: ᴬ_Z Xᵐ → ᴬ_Z X + γ (m = metastable state)";
    }

    /**
     * @brief Spin trap
     *
     * Large ΔJ → highly forbidden transition → long T_1/2
     */
    static std::string spin_trap() {
        return "High-spin isomers: large ΔJ inhibits decay (forbidden transition)";
    }

    /**
     * @brief Example isomers
     */
    static std::vector<std::string> examples() {
        return {"Tc-99m (6 hr)", "Co-60m (10.5 min)", "Ta-180m (>10¹⁵ yr)"};
    }
};

/**
 * @brief Decay Chains
 *
 * Sequential radioactive decays
 */
class DecayChains {
public:
    /**
     * @brief Bateman equations
     *
     * Series decay: A → B → C
     */
    static std::string bateman_equations() {
        return "dN_B/dt = λ_A N_A - λ_B N_B (sequential decay)";
    }

    /**
     * @brief Activity of daughter (no initial daughter)
     *
     * A_B(t) = (λ_B/(λ_B - λ_A)) A_A(0) (e^(-λ_A t) - e^(-λ_B t))
     */
    static double daughter_activity_bateman(double lambda_A, double lambda_B,
                                            double A_A_0, double t) {
        if (std::abs(lambda_B - lambda_A) < 1e-10) {
            // Special case: λ_A ≈ λ_B
            return lambda_A * A_A_0 * t * std::exp(-lambda_A * t);
        }
        return (lambda_B / (lambda_B - lambda_A)) * A_A_0 *
               (std::exp(-lambda_A * t) - std::exp(-lambda_B * t));
    }

    /**
     * @brief Secular equilibrium
     *
     * T_1/2(parent) >> T_1/2(daughter): A_B = A_A
     */
    static bool is_secular_equilibrium(double T_half_parent, double T_half_daughter) {
        return T_half_parent > 100.0 * T_half_daughter;
    }

    /**
     * @brief Transient equilibrium
     *
     * T_1/2(parent) > T_1/2(daughter): A_B/A_A = λ_B/(λ_B - λ_A)
     */
    static bool is_transient_equilibrium(double T_half_parent, double T_half_daughter) {
        return (T_half_parent > T_half_daughter) && (T_half_parent < 100.0 * T_half_daughter);
    }

    /**
     * @brief Activity ratio at transient equilibrium
     *
     * A_B/A_A = T_1/2(B) / (T_1/2(B) - T_1/2(A))
     */
    static double transient_eq_ratio(double T_half_A, double T_half_B) {
        return T_half_B / (T_half_B - T_half_A);
    }
};

/**
 * @brief Predicting Type of Decay
 *
 * Systematics for decay mode prediction
 */
class DecayPrediction {
public:
    /**
     * @brief N/Z ratio criterion
     *
     * N/Z > (N/Z)_stable → β⁻ decay
     * N/Z < (N/Z)_stable → β⁺/EC decay
     */
    static std::string predict_beta_type(int A, int Z) {
        int N = A - Z;
        double ratio = static_cast<double>(N) / Z;
        double stable_ratio = (A < 40) ? 1.0 : 1.0 + 0.015 * std::pow(A, 2.0/3.0);

        if (ratio > stable_ratio + 0.1) {
            return "β⁻ decay (neutron-rich)";
        } else if (ratio < stable_ratio - 0.1) {
            return "β⁺/EC decay (proton-rich)";
        } else {
            return "Stable or near-stable";
        }
    }

    /**
     * @brief Alpha decay criterion
     *
     * Z > 82: α decay likely
     */
    static bool alpha_decay_likely(int Z, int A) {
        return (Z > 82) || (A > 209);
    }

    /**
     * @brief Spontaneous fission
     *
     * Z² / A > 47: fission competes with α
     */
    static bool fission_competes(int Z, int A) {
        return (Z * Z / static_cast<double>(A)) > 47.0;
    }

    /**
     * @brief Proton drip line
     *
     * S_p < 0: proton unbound
     */
    static bool beyond_proton_drip(int A, int Z) {
        double S_p = NuclearStability::proton_separation_energy(A, Z);
        return S_p < 0;
    }

    /**
     * @brief Neutron drip line
     *
     * S_n < 0: neutron unbound
     */
    static bool beyond_neutron_drip(int A, int Z) {
        double S_n = NuclearStability::neutron_separation_energy(A, Z);
        return S_n < 0;
    }
};

/**
 * @brief Radioactive Decay Rates
 *
 * Decay constant, activity, and decay law
 */
class DecayRates {
public:
    /**
     * @brief Radioactive decay law
     *
     * N(t) = N₀ exp(-λt)
     */
    static double number_of_nuclei(double N_0, double lambda, double t) {
        return N_0 * std::exp(-lambda * t);
    }

    /**
     * @brief Activity
     *
     * A(t) = λN(t) = A₀ exp(-λt)
     */
    static double activity(double A_0, double lambda, double t) {
        return A_0 * std::exp(-lambda * t);
    }

    /**
     * @brief Activity from number of nuclei
     *
     * A = λN
     */
    static double activity_from_N(double N, double lambda) {
        return lambda * N;
    }

    /**
     * @brief Decay constant from half-life
     *
     * λ = ln(2) / T_1/2
     */
    static double lambda_from_halflife(double T_half) {
        return std::log(2.0) / T_half;
    }

    /**
     * @brief Half-life from decay constant
     *
     * T_1/2 = ln(2) / λ
     */
    static double halflife_from_lambda(double lambda) {
        return std::log(2.0) / lambda;
    }

    /**
     * @brief Mean lifetime
     *
     * τ = 1/λ = T_1/2 / ln(2)
     */
    static double mean_lifetime(double lambda) {
        return 1.0 / lambda;
    }

    /**
     * @brief Number of decays in time interval
     *
     * ΔN = N₀(1 - e^(-λΔt))
     */
    static double decays_in_interval(double N_0, double lambda, double delta_t) {
        return N_0 * (1.0 - std::exp(-lambda * delta_t));
    }

    /**
     * @brief Specific activity
     *
     * a = A/m = λN_A/M (activity per unit mass)
     */
    static double specific_activity(double lambda, double M_molar) {
        return lambda * constants::N_A / M_molar;  // Bq/g
    }
};

/**
 * @brief Units of Measurement for Radioactivity
 *
 * Becquerel, Curie, and conversions
 */
class RadioactivityUnits {
public:
    /**
     * @brief Becquerel (SI unit)
     *
     * 1 Bq = 1 decay/s
     */
    static std::string becquerel() {
        return "1 Bq = 1 disintegration/second (SI unit)";
    }

    /**
     * @brief Curie (traditional unit)
     *
     * 1 Ci = 3.7×10¹⁰ Bq
     */
    static constexpr double curie_to_becquerel(double Ci) {
        return Ci * 3.7e10;
    }

    /**
     * @brief Becquerel to Curie
     */
    static constexpr double becquerel_to_curie(double Bq) {
        return Bq / 3.7e10;
    }

    /**
     * @brief Rutherford (obsolete)
     *
     * 1 Rd = 10⁶ Bq
     */
    static constexpr double rutherford_to_becquerel(double Rd) {
        return Rd * 1.0e6;
    }

    /**
     * @brief Activity in mCi (common unit)
     */
    static double millicurie_to_becquerel(double mCi) {
        return mCi * 3.7e7;  // 1 mCi = 37 MBq
    }

    /**
     * @brief Activity in μCi
     */
    static double microcurie_to_becquerel(double uCi) {
        return uCi * 3.7e4;  // 1 μCi = 37 kBq
    }
};

/**
 * @brief Variation of Radioactivity Over Time
 *
 * Time-dependent decay behavior
 */
class RadioactivityVariation {
public:
    /**
     * @brief Activity as function of time
     *
     * A(t) = A₀ e^(-λt)
     */
    static double activity_at_time(double A_0, double T_half, double t) {
        double lambda = std::log(2.0) / T_half;
        return A_0 * std::exp(-lambda * t);
    }

    /**
     * @brief Time to reach target activity
     *
     * t = (1/λ) ln(A₀/A)
     */
    static double time_to_activity(double A_0, double A_target, double T_half) {
        double lambda = std::log(2.0) / T_half;
        return std::log(A_0 / A_target) / lambda;
    }

    /**
     * @brief Activity after n half-lives
     *
     * A(n T_1/2) = A₀ / 2ⁿ
     */
    static double activity_after_n_halflives(double A_0, int n) {
        return A_0 / std::pow(2.0, n);
    }

    /**
     * @brief Fraction remaining after time t
     *
     * f(t) = e^(-λt)
     */
    static double fraction_remaining(double T_half, double t) {
        double lambda = std::log(2.0) / T_half;
        return std::exp(-lambda * t);
    }

    /**
     * @brief Fraction decayed after time t
     *
     * 1 - f(t) = 1 - e^(-λt)
     */
    static double fraction_decayed(double T_half, double t) {
        return 1.0 - fraction_remaining(T_half, t);
    }
};

/**
 * @brief Radioactive Half-Life
 *
 * Characteristic decay time
 */
class HalfLife {
public:
    /**
     * @brief Definition
     *
     * T_1/2: time for N → N/2
     */
    static std::string definition() {
        return "T_1/2: time for half of nuclei to decay (N → N/2)";
    }

    /**
     * @brief Relation to decay constant
     *
     * T_1/2 = ln(2)/λ ≈ 0.693/λ
     */
    static constexpr double ln2() {
        return 0.693147180559945309417;
    }

    /**
     * @brief Effective half-life
     *
     * T_eff = T_physical × T_biological / (T_physical + T_biological)
     */
    static double effective_halflife(double T_phys, double T_bio) {
        return (T_phys * T_bio) / (T_phys + T_bio);
    }

    /**
     * @brief Half-life range in nature
     */
    static std::string natural_range() {
        return "From 10⁻²³ s (Be-8) to >10¹⁸ yr (Te-128, Xe-136)";
    }

    /**
     * @brief Number of half-lives for 99% decay
     *
     * 0.01 = (1/2)ⁿ → n ≈ 6.64
     */
    static double halflives_for_99_percent_decay() {
        return std::log(100.0) / std::log(2.0);  // ≈ 6.64
    }

    /**
     * @brief Number of half-lives for decay to target fraction
     */
    static double halflives_to_fraction(double fraction) {
        return std::log(1.0 / fraction) / std::log(2.0);
    }
};

/**
 * @brief Plotting Radioactive Decay
 *
 * Visualization and data representation
 */
class DecayPlotting {
public:
    /**
     * @brief Semi-log plot
     *
     * log(A) vs t gives straight line with slope -λ
     */
    static std::string semilog_plot() {
        return "ln(A) vs t: straight line with slope -λ, intercept ln(A₀)";
    }

    /**
     * @brief Slope from semi-log plot
     *
     * λ = -Δ(ln A) / Δt
     */
    static double lambda_from_semilog(double ln_A1, double ln_A2, double t1, double t2) {
        return -(ln_A2 - ln_A1) / (t2 - t1);
    }

    /**
     * @brief Generate decay curve data points
     *
     * Returns A(t) at n equally-spaced time points
     */
    static std::vector<double> generate_decay_curve(double A_0, double lambda,
                                                     double t_max, int n_points) {
        std::vector<double> activities;
        for (int i = 0; i < n_points; ++i) {
            double t = (i * t_max) / (n_points - 1);
            activities.push_back(A_0 * std::exp(-lambda * t));
        }
        return activities;
    }

    /**
     * @brief Decay chain visualization data
     *
     * Returns {A_parent, A_daughter} at time t
     */
    static std::array<double, 2> decay_chain_point(double A_A_0, double lambda_A,
                                                    double lambda_B, double t) {
        double A_A = A_A_0 * std::exp(-lambda_A * t);
        double A_B = DecayChains::daughter_activity_bateman(lambda_A, lambda_B, A_A_0, t);
        return {A_A, A_B};
    }
};

/**
 * @brief Radioactive Equilibrium
 *
 * Secular and transient equilibrium
 */
class RadioactiveEquilibrium {
public:
    /**
     * @brief Secular equilibrium condition
     *
     * λ_parent << λ_daughter: A_daughter = A_parent
     */
    static std::string secular_condition() {
        return "λ_A << λ_B (T_A >> T_B): A_B → A_A (equilibrium)";
    }

    /**
     * @brief Secular equilibrium activity
     *
     * A_B = A_A for t >> 1/λ_B
     */
    static double secular_daughter_activity(double A_parent) {
        return A_parent;  // In equilibrium
    }

    /**
     * @brief Transient equilibrium condition
     *
     * λ_parent < λ_daughter: A_B/A_A = λ_B/(λ_B - λ_A)
     */
    static std::string transient_condition() {
        return "λ_A < λ_B (T_A > T_B): A_B/A_A = λ_B/(λ_B - λ_A) > 1";
    }

    /**
     * @brief Transient equilibrium ratio
     *
     * A_B/A_A at equilibrium
     */
    static double transient_activity_ratio(double lambda_A, double lambda_B) {
        return lambda_B / (lambda_B - lambda_A);
    }

    /**
     * @brief Time to reach equilibrium
     *
     * t_eq ≈ 5/λ_daughter (rule of thumb)
     */
    static double time_to_equilibrium(double lambda_daughter) {
        return 5.0 / lambda_daughter;
    }

    /**
     * @brief Maximum daughter activity in transient equilibrium
     *
     * Occurs at t_max = ln(λ_B/λ_A) / (λ_B - λ_A)
     */
    static double time_of_max_daughter(double lambda_A, double lambda_B) {
        return std::log(lambda_B / lambda_A) / (lambda_B - lambda_A);
    }

    /**
     * @brief Maximum daughter activity value
     */
    static double max_daughter_activity(double A_A_0, double lambda_A, double lambda_B) {
        double t_max = time_of_max_daughter(lambda_A, lambda_B);
        return DecayChains::daughter_activity_bateman(lambda_A, lambda_B, A_A_0, t_max);
    }

    /**
     * @brief No equilibrium
     *
     * λ_parent > λ_daughter: daughter outlives parent
     */
    static std::string no_equilibrium() {
        return "λ_A > λ_B (T_A < T_B): no equilibrium, daughter accumulates";
    }
};

/**
 * @brief Neutron Scattering
 *
 * Elastic and inelastic scattering of neutrons
 */
class NeutronScattering {
public:
    /**
     * @brief Elastic scattering overview
     *
     * n + A → n + A (KE conserved in CM frame)
     */
    static std::string elastic_description() {
        return "Elastic: (n,n) - neutron scatters, nucleus recoils, no excitation";
    }

    /**
     * @brief Inelastic scattering overview
     *
     * n + A → n + A* (nucleus left in excited state)
     */
    static std::string inelastic_description() {
        return "Inelastic: (n,n') - neutron loses energy, nucleus excited";
    }

    // =========================================================================
    // COMPUTATIONAL FUNCTIONS
    // =========================================================================

    /**
     * @brief Average scattering cosine (elastic)
     *
     * μ_avg = 2/(3A) for A >> 1
     *
     * @param A Mass number
     * @return Average cosine of scattering angle
     */
    static double average_scattering_cosine(int A) {
        if (A <= 0) return 0.0;
        return 2.0 / (3.0 * A);
    }

    /**
     * @brief Average logarithmic energy decrement (elastic)
     *
     * ξ = 1 - [(A-1)²/(2A)] ln[(A+1)/(A-1)]
     *
     * @param A Mass number
     * @return Average log energy loss per collision
     */
    static double average_log_energy_decrement(int A) {
        if (A <= 1) return 1.0;  // Special case for hydrogen

        double A_plus = A + 1.0;
        double A_minus = A - 1.0;
        double term = (A_minus * A_minus) / (2.0 * A) * std::log(A_plus / A_minus);
        return 1.0 - term;
    }

    /**
     * @brief Neutron energy after elastic collision
     *
     * E' = E [(A² + 1 + 2A cos(θ_CM)) / (A+1)²]
     *
     * @param E_initial Initial neutron energy (eV)
     * @param A Mass number
     * @param theta_cm Scattering angle in CM frame (radians)
     * @return Final neutron energy (eV)
     */
    static double energy_after_elastic(double E_initial, int A, double theta_cm) {
        if (A <= 0) return E_initial;

        double A_plus_1 = A + 1.0;
        double numerator = A * A + 1.0 + 2.0 * A * std::cos(theta_cm);
        double denominator = A_plus_1 * A_plus_1;

        return E_initial * numerator / denominator;
    }

    /**
     * @brief Minimum energy after elastic scatter (backscatter)
     *
     * E_min = E [(A-1)/(A+1)]²
     *
     * @param E_initial Initial energy (eV)
     * @param A Mass number
     * @return Minimum possible energy (eV)
     */
    static double minimum_energy_elastic(double E_initial, int A) {
        if (A <= 0) return E_initial;

        double ratio = (A - 1.0) / (A + 1.0);
        return E_initial * ratio * ratio;
    }

    /**
     * @brief Maximum energy after elastic scatter (forward)
     *
     * E_max = E (always E for elastic)
     *
     * @param E_initial Initial energy (eV)
     * @return Maximum energy (eV)
     */
    static double maximum_energy_elastic(double E_initial) {
        return E_initial;  // Forward scatter, no energy loss
    }

    /**
     * @brief Number of collisions to slow down
     *
     * n ≈ ln(E_initial/E_final) / ξ
     *
     * @param E_initial Initial energy (eV)
     * @param E_final Final energy (eV)
     * @param A Mass number
     * @return Average number of collisions
     */
    static double collisions_to_slow(double E_initial, double E_final, int A) {
        if (E_initial <= E_final || E_final <= 0.0) return 0.0;

        double xi = average_log_energy_decrement(A);
        if (xi <= 0.0) return 0.0;

        return std::log(E_initial / E_final) / xi;
    }

    /**
     * @brief Inelastic scattering threshold energy
     *
     * E_threshold = Q [(A+1)/A] where Q = excitation energy
     *
     * @param Q_value Excitation energy of level (MeV)
     * @param A Mass number
     * @return Threshold neutron energy (MeV)
     */
    static double inelastic_threshold(double Q_value, int A) {
        if (A <= 0) return 0.0;
        return std::abs(Q_value) * (A + 1.0) / A;
    }

    /**
     * @brief Neutron energy after inelastic scatter
     *
     * E_n' ≈ E_n - Q - E_recoil
     *
     * @param E_initial Initial neutron energy (MeV)
     * @param Q_value Excitation energy (MeV, positive)
     * @param A Mass number
     * @return Final neutron energy (MeV)
     */
    static double energy_after_inelastic(double E_initial, double Q_value, int A) {
        // Simple approximation: E' ≈ E - Q(A+1)/A
        double threshold = inelastic_threshold(Q_value, A);
        if (E_initial < threshold) return 0.0;  // Below threshold

        double energy_loss = Q_value * (A + 1.0) / A;
        return std::max(0.0, E_initial - energy_loss);
    }

    /**
     * @brief Scattering cross-section energy dependence (simple model)
     *
     * σ_s ≈ σ_0 for E > resonance region
     * More complex near resonances
     *
     * @param sigma_0 Constant cross-section (barns)
     * @param energy Neutron energy (eV)
     * @return Scattering cross-section (barns)
     */
    static double scattering_cross_section(double sigma_0, double energy) {
        // Simplified: constant cross-section at high energy
        // Real cross-sections have resonances and energy dependence
        return sigma_0;
    }

    /**
     * @brief Slowing-down power
     *
     * ξΣ_s = moderating power
     *
     * @param xi Average log energy decrement
     * @param sigma_s Macroscopic scattering cross-section (cm⁻¹)
     * @return Slowing-down power (cm⁻¹)
     */
    static double slowing_down_power(double xi, double sigma_s) {
        return xi * sigma_s;
    }

    /**
     * @brief Moderating ratio
     *
     * ξΣ_s / Σ_a = figure of merit for moderators
     *
     * @param xi Average log energy decrement
     * @param sigma_s Scattering cross-section (cm⁻¹)
     * @param sigma_a Absorption cross-section (cm⁻¹)
     * @return Moderating ratio
     */
    static double moderating_ratio(double xi, double sigma_s, double sigma_a) {
        if (sigma_a <= 0.0) return 1e10;  // Perfect moderator
        return (xi * sigma_s) / sigma_a;
    }
};

/**
 * @brief Neutron Absorption
 *
 * Radiative capture reactions
 */
class NeutronAbsorption {
public:
    /**
     * @brief Radiative capture overview
     *
     * n + A → (A+1) + γ
     */
    static std::string radiative_capture_description() {
        return "Radiative capture: (n,γ) - neutron absorbed, gamma ray emitted";
    }

    // =========================================================================
    // COMPUTATIONAL FUNCTIONS
    // =========================================================================

    /**
     * @brief Q-value for radiative capture
     *
     * Q = (M_A + m_n - M_{A+1})c² = binding energy of neutron
     *
     * @param A Mass number of target
     * @param Z Atomic number
     * @return Q-value (MeV, typically 6-8 MeV)
     */
    static double q_value_radiative_capture(int A, int Z) {
        // Approximate neutron binding energy using SEMF
        // B(A+1,Z) - B(A,Z) ≈ neutron separation energy

        // Simplified: typical values are 6-8 MeV
        // More accurate would use actual binding energies
        return 7.0;  // Typical value in MeV
    }

    /**
     * @brief Capture cross-section (1/v law)
     *
     * σ_capture ∝ 1/v for thermal neutrons
     * σ(E) = σ_0 √(E_0/E)
     *
     * @param sigma_0 Cross-section at reference energy (barns)
     * @param E_0 Reference energy (eV, typically 0.0253 eV)
     * @param E Neutron energy (eV)
     * @return Capture cross-section (barns)
     */
    static double capture_cross_section_thermal(double sigma_0, double E_0, double E) {
        if (E <= 0.0) return 0.0;
        return sigma_0 * std::sqrt(E_0 / E);
    }

    /**
     * @brief Breit-Wigner resonance cross-section
     *
     * σ(E) = σ_max [Γ²/4] / [(E - E_R)² + Γ²/4]
     *
     * @param sigma_max Peak cross-section (barns)
     * @param E_resonance Resonance energy (eV)
     * @param gamma_total Total width (eV)
     * @param E Neutron energy (eV)
     * @return Cross-section at energy E (barns)
     */
    static double breit_wigner_resonance(double sigma_max, double E_resonance,
                                          double gamma_total, double E) {
        double delta_E = E - E_resonance;
        double gamma_half = gamma_total / 2.0;
        double denominator = delta_E * delta_E + gamma_half * gamma_half;

        return sigma_max * (gamma_half * gamma_half) / denominator;
    }

    /**
     * @brief Resonance integral
     *
     * I = ∫ σ(E) dE/E over epithermal range
     * Approximation for single resonance
     *
     * @param sigma_0 Thermal cross-section (barns)
     * @param E_resonance Resonance energy (eV)
     * @param gamma Width (eV)
     * @return Resonance integral (barns)
     */
    static double resonance_integral_single(double sigma_0, double E_resonance, double gamma) {
        // Simplified formula: I ≈ (π/2) σ_peak Γ / E_R
        // For detailed calculations, numerical integration needed
        return (M_PI / 2.0) * sigma_0 * gamma / E_resonance;
    }

    /**
     * @brief Activation rate
     *
     * R = φ N σ_a (reactions/cm³/s)
     *
     * @param flux Neutron flux (n/cm²/s)
     * @param number_density Atom density (atoms/cm³)
     * @param sigma_absorption Cross-section (barns = 10⁻²⁴ cm²)
     * @return Reaction rate (reactions/cm³/s)
     */
    static double activation_rate(double flux, double number_density, double sigma_absorption) {
        double sigma_cm2 = sigma_absorption * 1.0e-24;  // Convert barns to cm²
        return flux * number_density * sigma_cm2;
    }

    /**
     * @brief Saturated activity from irradiation
     *
     * A_sat = R (1 - e^(-λt))
     *
     * @param reaction_rate Production rate (atoms/s)
     * @param decay_constant Decay constant (s⁻¹)
     * @param irradiation_time Time (s)
     * @return Activity (Bq)
     */
    static double activity_from_irradiation(double reaction_rate, double decay_constant,
                                             double irradiation_time) {
        return reaction_rate * (1.0 - std::exp(-decay_constant * irradiation_time));
    }

    /**
     * @brief Saturation activity (infinite irradiation)
     *
     * A_sat = R
     *
     * @param reaction_rate Production rate (atoms/s)
     * @return Maximum activity (Bq)
     */
    static double saturation_activity(double reaction_rate) {
        return reaction_rate;
    }

    /**
     * @brief Effective cross-section in reactor spectrum
     *
     * σ_eff = ∫ σ(E) φ(E) dE / ∫ φ(E) dE
     *
     * @param sigma_thermal Thermal cross-section (barns)
     * @param resonance_integral Epithermal contribution (barns)
     * @param thermal_fraction Fraction of thermal flux
     * @return Effective cross-section (barns)
     */
    static double effective_cross_section(double sigma_thermal, double resonance_integral,
                                           double thermal_fraction) {
        double epithermal_fraction = 1.0 - thermal_fraction;
        return sigma_thermal * thermal_fraction + resonance_integral * epithermal_fraction;
    }

    /**
     * @brief Self-shielding factor
     *
     * f = σ_eff / σ_infinite (< 1 for large samples)
     *
     * @param optical_thickness τ = Σt × dimension
     * @return Self-shielding factor
     */
    static double self_shielding_factor(double optical_thickness) {
        // Simplified: f ≈ 1/(1 + τ) for slab geometry
        return 1.0 / (1.0 + optical_thickness);
    }
};

/**
 * @brief Particle Ejection Reactions
 *
 * (n,p), (n,α), (n,2n) and other particle emission
 */
class ParticleEjection {
public:
    /**
     * @brief (n,p) reaction overview
     *
     * n + A(Z) → p + A(Z-1)
     */
    static std::string np_description() {
        return "(n,p): neutron in, proton out - changes Z by -1";
    }

    /**
     * @brief (n,α) reaction overview
     *
     * n + A(Z) → α + (A-3)(Z-2)
     */
    static std::string nalpha_description() {
        return "(n,α): neutron in, alpha out - changes Z by -2, A by -3";
    }

    /**
     * @brief (n,2n) reaction overview
     *
     * n + A → 2n + (A-1)
     */
    static std::string n2n_description() {
        return "(n,2n): neutron multiplication - high energy threshold";
    }

    // =========================================================================
    // COMPUTATIONAL FUNCTIONS
    // =========================================================================

    /**
     * @brief Q-value for (n,p) reaction
     *
     * Q = [M(A,Z) + m_n - M(A,Z-1) - m_p]c²
     *
     * @param A Mass number
     * @param Z Atomic number
     * @return Q-value (MeV, typically negative, endothermic)
     */
    static double q_value_np(int A, int Z) {
        // Approximation: Q ≈ (m_n - m_p)c² + binding energy difference
        // m_n - m_p ≈ 1.293 MeV
        // Usually endothermic for most nuclei
        return -0.8;  // Typical value, usually negative
    }

    /**
     * @brief Q-value for (n,α) reaction
     *
     * Q = [M(A,Z) + m_n - M(A-3,Z-2) - m_α]c²
     *
     * @param A Mass number
     * @param Z Atomic number
     * @return Q-value (MeV)
     */
    static double q_value_nalpha(int A, int Z) {
        // Highly dependent on nucleus
        // Often exothermic for light nuclei, endothermic for heavy
        return 0.0;  // Placeholder, varies widely
    }

    /**
     * @brief Q-value for (n,2n) reaction
     *
     * Q = [M(A,Z) + m_n - M(A-1,Z) - 2m_n]c²
     * Q = -S_n (negative of neutron separation energy)
     *
     * @param A Mass number
     * @param Z Atomic number
     * @return Q-value (MeV, always negative)
     */
    static double q_value_n2n(int A, int Z) {
        // Q ≈ -8 MeV typically (neutron separation energy)
        return -8.0;
    }

    /**
     * @brief Threshold energy for particle ejection
     *
     * E_threshold = -Q [(A_product + m_product)/(A_target)] for Q < 0
     *
     * @param Q_value Q-value (MeV)
     * @param A_target Mass number of target
     * @param A_product Total mass of products
     * @return Threshold energy (MeV)
     */
    static double threshold_energy(double Q_value, int A_target, double A_product) {
        if (Q_value >= 0.0) return 0.0;  // Exothermic, no threshold

        return -Q_value * (A_product / A_target);
    }

    /**
     * @brief (n,p) threshold energy
     *
     * Typically 1-5 MeV
     *
     * @param A Mass number
     * @return Threshold (MeV)
     */
    static double np_threshold(int A) {
        double Q = q_value_np(A, 0);  // Z not critical for estimate
        return threshold_energy(Q, A, A + 1.0);
    }

    /**
     * @brief (n,2n) threshold energy
     *
     * Typically 8-10 MeV
     *
     * @param A Mass number
     * @return Threshold (MeV)
     */
    static double n2n_threshold(int A) {
        double Q = q_value_n2n(A, 0);
        // For (n,2n): threshold ≈ -Q(A+2)/(A-1)
        return -Q * (A + 2.0) / (A - 1.0);
    }

    /**
     * @brief (n,p) cross-section estimate
     *
     * σ(n,p) peaks around 2-5 MeV, then decreases
     *
     * @param energy Neutron energy (MeV)
     * @param sigma_max Peak cross-section (barns)
     * @param E_peak Peak energy (MeV)
     * @return Cross-section (barns)
     */
    static double np_cross_section(double energy, double sigma_max, double E_peak) {
        if (energy < np_threshold(50)) return 0.0;  // Below threshold

        // Simplified Gaussian-like peak
        double exponent = -(energy - E_peak) * (energy - E_peak) / (E_peak * E_peak);
        return sigma_max * std::exp(exponent);
    }

    /**
     * @brief (n,2n) cross-section estimate
     *
     * Rises from threshold, peaks around 14 MeV
     *
     * @param energy Neutron energy (MeV)
     * @param sigma_max Peak cross-section (barns, typically 0.5-2 b)
     * @return Cross-section (barns)
     */
    static double n2n_cross_section(double energy, double sigma_max) {
        double threshold = 8.0;  // Typical threshold

        if (energy < threshold) return 0.0;

        // Simplified: rises from threshold, peaks near 14 MeV
        double E_peak = 14.0;
        if (energy < E_peak) {
            // Rising from threshold
            return sigma_max * (energy - threshold) / (E_peak - threshold);
        } else {
            // Decreasing after peak
            return sigma_max * std::exp(-(energy - E_peak) / 5.0);
        }
    }

    /**
     * @brief Neutron multiplication from (n,2n)
     *
     * ν = average neutrons out per neutron in
     * ν = 1 + σ(n,2n)/σ_total
     *
     * @param sigma_n2n (n,2n) cross-section (barns)
     * @param sigma_total Total cross-section (barns)
     * @return Neutron multiplication factor
     */
    static double neutron_multiplication(double sigma_n2n, double sigma_total) {
        if (sigma_total <= 0.0) return 1.0;
        return 1.0 + sigma_n2n / sigma_total;
    }
};

/**
 * @brief Neutron-Induced Fission
 *
 * Fission reactions and characteristics
 */
class NeutronInducedFission {
public:
    /**
     * @brief Fission overview
     *
     * n + A → fission fragments + neutrons + energy
     */
    static std::string fission_description() {
        return "Fission: nucleus splits into 2 fragments + 2-3 neutrons + ~200 MeV";
    }

    // =========================================================================
    // COMPUTATIONAL FUNCTIONS
    // =========================================================================

    /**
     * @brief Fission Q-value
     *
     * Q ≈ 200 MeV for U-235, Pu-239
     *
     * @param A Mass number (typically 235, 239)
     * @return Energy release (MeV)
     */
    static double fission_q_value(int A) {
        // Approximately 0.85 MeV per nucleon
        return 0.85 * A;
    }

    /**
     * @brief Fission cross-section for U-235 (thermal)
     *
     * σ_f = 585 barns at 0.0253 eV
     *
     * @return Thermal fission cross-section (barns)
     */
    static double u235_fission_cross_section_thermal() {
        return 585.0;  // barns
    }

    /**
     * @brief Fission cross-section for Pu-239 (thermal)
     *
     * σ_f = 747 barns at 0.0253 eV
     *
     * @return Thermal fission cross-section (barns)
     */
    static double pu239_fission_cross_section_thermal() {
        return 747.0;  // barns
    }

    /**
     * @brief Fission cross-section for U-238 (fast)
     *
     * Threshold at ~1 MeV, rises to ~0.5 barns at 14 MeV
     *
     * @param energy Neutron energy (MeV)
     * @return Fission cross-section (barns)
     */
    static double u238_fission_cross_section(double energy) {
        double threshold = 1.0;  // MeV
        if (energy < threshold) return 0.0;

        // Simplified: rises roughly linearly
        return 0.5 * (energy - threshold) / 13.0;
    }

    /**
     * @brief Average neutrons per fission
     *
     * ν̄ for different fissile isotopes
     *
     * @param isotope "U235", "U238", "Pu239"
     * @param energy Neutron energy (MeV)
     * @return Average neutron yield
     */
    static double average_neutrons_per_fission(const std::string& isotope, double energy) {
        // Approximate values, slightly energy dependent
        // ν̄ = ν̄_0 + slope × E

        if (isotope == "U235") {
            return 2.42 + 0.066 * energy;  // Thermal: 2.42, increases with E
        } else if (isotope == "U238") {
            return 2.47 + 0.13 * energy;   // Higher slope for fast
        } else if (isotope == "Pu239") {
            return 2.87 + 0.09 * energy;   // Highest yield
        }
        return 2.5;  // Default
    }

    /**
     * @brief Prompt neutron fraction
     *
     * Fraction of neutrons emitted promptly (< 10⁻¹⁴ s)
     *
     * @return Prompt fraction (typically 0.993-0.997)
     */
    static double prompt_neutron_fraction() {
        return 0.993;  // ~99.3% prompt for U-235
    }

    /**
     * @brief Delayed neutron fraction (β)
     *
     * Critical for reactor control
     *
     * @param isotope "U235", "U238", "Pu239"
     * @return Delayed fraction β
     */
    static double delayed_neutron_fraction(const std::string& isotope) {
        if (isotope == "U235") {
            return 0.0065;  // 0.65%
        } else if (isotope == "U238") {
            return 0.0148;  // 1.48%
        } else if (isotope == "Pu239") {
            return 0.0021;  // 0.21%
        }
        return 0.0065;  // Default to U-235
    }

    /**
     * @brief Fission energy distribution
     *
     * ~200 MeV total, distributed as:
     * Fission fragments: ~168 MeV
     * Neutrons: ~5 MeV
     * Prompt gammas: ~7 MeV
     * Beta decay: ~7 MeV
     * Neutrinos: ~10 MeV (lost)
     * Gamma decay: ~6 MeV
     *
     * @param component "fragments", "neutrons", "prompt_gamma", "beta", "neutrinos", "delayed_gamma"
     * @return Energy (MeV)
     */
    static double fission_energy_component(const std::string& component) {
        if (component == "fragments") return 168.0;
        if (component == "neutrons") return 5.0;
        if (component == "prompt_gamma") return 7.0;
        if (component == "beta") return 7.0;
        if (component == "neutrinos") return 10.0;  // Lost, not deposited
        if (component == "delayed_gamma") return 6.0;
        return 0.0;
    }

    /**
     * @brief Recoverable energy per fission
     *
     * Total minus neutrinos
     *
     * @return Recoverable energy (MeV)
     */
    static double recoverable_energy_per_fission() {
        return 200.0 - 10.0;  // 190 MeV (minus neutrinos)
    }

    /**
     * @brief Fission product yield (mass distribution)
     *
     * Bimodal distribution with peaks at A~95 and A~135
     *
     * @param A_fragment Mass number of fragment
     * @return Yield (percent per fission)
     */
    static double fission_yield_mass(int A_fragment) {
        // Simplified bimodal Gaussian distribution
        // Light peak around A=95, heavy peak around A=135

        double light_peak = 95.0;
        double heavy_peak = 135.0;
        double width = 5.0;

        double light_contrib = 0.03 * std::exp(-(A_fragment - light_peak) * (A_fragment - light_peak) / (2.0 * width * width));
        double heavy_contrib = 0.04 * std::exp(-(A_fragment - heavy_peak) * (A_fragment - heavy_peak) / (2.0 * width * width));

        return (light_contrib + heavy_contrib) * 100.0;  // Percent
    }

    /**
     * @brief Multiplication factor (k-infinity)
     *
     * k∞ = (ν σ_f) / (σ_a) = (ν σ_f) / (σ_f + σ_c)
     *
     * @param nu Average neutrons per fission
     * @param sigma_fission Fission cross-section (barns)
     * @param sigma_capture Capture cross-section (barns)
     * @return Infinite multiplication factor
     */
    static double k_infinity(double nu, double sigma_fission, double sigma_capture) {
        double sigma_absorption = sigma_fission + sigma_capture;
        if (sigma_absorption <= 0.0) return 0.0;

        return (nu * sigma_fission) / sigma_absorption;
    }

    /**
     * @brief Reproduction factor (η - eta)
     *
     * η = ν σ_f / σ_a (neutrons produced per absorption in fuel)
     *
     * @param nu Average neutrons per fission
     * @param sigma_fission Fission cross-section
     * @param sigma_absorption Total absorption cross-section
     * @return Eta factor
     */
    static double eta_factor(double nu, double sigma_fission, double sigma_absorption) {
        if (sigma_absorption <= 0.0) return 0.0;
        return (nu * sigma_fission) / sigma_absorption;
    }

    /**
     * @brief Alpha ratio
     *
     * α = σ_c / σ_f (capture to fission ratio)
     *
     * @param sigma_capture Capture cross-section
     * @param sigma_fission Fission cross-section
     * @return Alpha ratio
     */
    static double alpha_ratio(double sigma_capture, double sigma_fission) {
        if (sigma_fission <= 0.0) return 0.0;
        return sigma_capture / sigma_fission;
    }
};

} // namespace nuclear
} // namespace physics

#endif // PHYSICS_NUCLEAR_PHYSICS_HPP
