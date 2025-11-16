#ifndef PHYSICS_STATISTICAL_MECHANICS_HPP
#define PHYSICS_STATISTICAL_MECHANICS_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>
#include <random>
#include <map>

namespace physics {
namespace statistical_mechanics {

constexpr double k_B = 1.380649e-23;
constexpr double h = 6.62607015e-34;
constexpr double hbar = 1.054571817e-34;

// ============================================================================
// MICROCANONICAL ENSEMBLE
// ============================================================================

class MicrocanonicalEnsemble {
private:
    double E_;
    int N_;
    double V_;

public:
    MicrocanonicalEnsemble(double energy, int particles, double volume)
        : E_(energy), N_(particles), V_(volume) {}

    double multiplicityIdealGas(double deltaE) const {
        double f = 3.0 * N_ / 2.0;
        double C = std::pow(V_, N_) / std::tgamma(f + 1);
        double m = 1.0;
        return C * std::pow(2.0 * M_PI * m * E_ / (h * h), f) * deltaE;
    }

    double entropyBoltzmann(double Omega) const {
        return k_B * std::log(Omega);
    }

    double temperature(double dS_dE) const {
        return 1.0 / dS_dE;
    }

    double pressure(double dS_dV) const {
        return temperature(1.0) * dS_dV;
    }

    double entropyIdealGas() const {
        double lambda_thermal = h / std::sqrt(2.0 * M_PI * 1.0 * k_B * temperature(1.0));
        double S = N_ * k_B * (std::log(V_ / (N_ * std::pow(lambda_thermal, 3))) + 2.5);
        return S;
    }
};

// ============================================================================
// CANONICAL ENSEMBLE
// ============================================================================

class CanonicalEnsemble {
private:
    double T_;
    int N_;
    double V_;

public:
    CanonicalEnsemble(double temperature, int particles, double volume)
        : T_(temperature), N_(particles), V_(volume) {}

    double beta() const { return 1.0 / (k_B * T_); }

    double partitionFunction(std::function<double(int)> energy_levels,
                            int max_level) const {
        double Z = 0.0;
        for (int n = 0; n < max_level; ++n) {
            Z += std::exp(-beta() * energy_levels(n));
        }
        return Z;
    }

    double partitionFunctionIdealGas() const {
        double lambda = h / std::sqrt(2.0 * M_PI * 1.0 * k_B * T_);
        return std::pow(V_ / std::pow(lambda, 3), N_) / std::tgamma(N_ + 1);
    }

    double helmholtzFreeEnergy(double Z) const {
        return -k_B * T_ * std::log(Z);
    }

    double internalEnergy(double Z, double dZ_dT) const {
        return -std::pow(1.0 / (k_B * T_), 2) * (1.0 / Z) * dZ_dT;
    }

    double entropy(double F, double E) const {
        return (E - F) / T_;
    }

    double pressure(double F, double dF_dV) const {
        return -dF_dV;
    }

    double heatCapacity(double dE_dT) const {
        return dE_dT;
    }

    double boltzmannProbability(double E_i, double Z) const {
        return std::exp(-beta() * E_i) / Z;
    }

    template<typename T>
    T ensembleAverage(std::function<T(int)> observable,
                     std::function<double(int)> energy_levels,
                     int max_level) const {
        double Z = partitionFunction(energy_levels, max_level);
        T avg = T();
        for (int n = 0; n < max_level; ++n) {
            double prob = boltzmannProbability(energy_levels(n), Z);
            avg = avg + observable(n) * prob;
        }
        return avg;
    }
};

// ============================================================================
// GRAND CANONICAL ENSEMBLE
// ============================================================================

class GrandCanonicalEnsemble {
private:
    double T_;
    double mu_;
    double V_;

public:
    GrandCanonicalEnsemble(double temperature, double chemical_potential,
                          double volume)
        : T_(temperature), mu_(mu_), V_(volume) {}

    double beta() const { return 1.0 / (k_B * T_); }

    double grandPartitionFunction(
        std::function<double(int, int)> energy_levels,
        int max_N, int max_level) const {

        double Xi = 0.0;
        for (int N = 0; N <= max_N; ++N) {
            for (int n = 0; n < max_level; ++n) {
                double E = energy_levels(N, n);
                Xi += std::exp(-beta() * (E - mu_ * N));
            }
        }
        return Xi;
    }

    double grandPartitionFunctionIdealGas() const {
        double lambda = h / std::sqrt(2.0 * M_PI * 1.0 * k_B * T_);
        double z = std::exp(beta() * mu_);
        return std::exp(z * V_ / std::pow(lambda, 3));
    }

    double grandPotential(double Xi) const {
        return -k_B * T_ * std::log(Xi);
    }

    double averageParticleNumber(double Xi, double dXi_dmu) const {
        return k_B * T_ * (1.0 / Xi) * dXi_dmu;
    }

    double pressure(double Omega) const {
        return -Omega / V_;
    }

    double fluctuationN(double avg_N, double avg_N2) const {
        return std::sqrt(avg_N2 - avg_N * avg_N);
    }

    double compressibility(double dN_dmu) const {
        return V_ / (k_B * T_) * dN_dmu;
    }
};

// ============================================================================
// PARTITION FUNCTIONS FOR SPECIFIC SYSTEMS
// ============================================================================

class PartitionFunctions {
public:
    static double harmonicOscillator(double omega, double T) {
        double beta = 1.0 / (k_B * T);
        double x = beta * hbar * omega;
        return 1.0 / (2.0 * std::sinh(x / 2.0));
    }

    static double harmonicOscillatorQuantum(double omega, double T, int max_n = 100) {
        double beta = 1.0 / (k_B * T);
        double Z = 0.0;
        for (int n = 0; n < max_n; ++n) {
            double E_n = hbar * omega * (n + 0.5);
            Z += std::exp(-beta * E_n);
        }
        return Z;
    }

    static double rotator2D(double I, double T) {
        double beta = 1.0 / (k_B * T);
        return 2.0 * M_PI * I / (beta * hbar * hbar);
    }

    static double rotator3D(double I, double T, int max_l = 50) {
        double beta = 1.0 / (k_B * T);
        double Z = 0.0;
        for (int l = 0; l < max_l; ++l) {
            double E_l = hbar * hbar * l * (l + 1) / (2.0 * I);
            Z += (2 * l + 1) * std::exp(-beta * E_l);
        }
        return Z;
    }

    static double twoLevelSystem(double epsilon, double T) {
        double beta = 1.0 / (k_B * T);
        return 1.0 + std::exp(-beta * epsilon);
    }

    static double paramagneticSpin(double B, double mu_B, double S, double T) {
        double beta = 1.0 / (k_B * T);
        double x = beta * mu_B * B;
        int two_S_plus_1 = static_cast<int>(2 * S + 1);

        double Z = 0.0;
        for (int m = -static_cast<int>(S); m <= static_cast<int>(S); ++m) {
            Z += std::exp(beta * mu_B * B * m);
        }
        return Z;
    }
};

// ============================================================================
// PHASE TRANSITIONS
// ============================================================================

class PhaseTransitions {
public:
    static double vanDerWaalsPressure(double T, double v, double a, double b) {
        return (k_B * T) / (v - b) - a / (v * v);
    }

    static double criticalTemperature(double a, double b) {
        return (8.0 * a) / (27.0 * k_B * b);
    }

    static double criticalPressure(double a, double b) {
        return a / (27.0 * b * b);
    }

    static double criticalVolume(double b) {
        return 3.0 * b;
    }

    static double maxwellConstruction(
        std::function<double(double)> pressure_vs_volume,
        double v_liquid, double v_gas, int n_points = 1000) {

        double dv = (v_gas - v_liquid) / n_points;
        double area_above = 0.0, area_below = 0.0;
        double p_coexist = pressure_vs_volume((v_liquid + v_gas) / 2.0);

        for (int i = 0; i < n_points; ++i) {
            double v = v_liquid + i * dv;
            double p = pressure_vs_volume(v);
            if (p > p_coexist) {
                area_above += (p - p_coexist) * dv;
            } else {
                area_below += (p_coexist - p) * dv;
            }
        }

        if (std::abs(area_above - area_below) < 1e-6) {
            return p_coexist;
        }
        return -1.0;
    }

    static double orderParameter(double T, double T_c, double beta_exponent) {
        if (T >= T_c) return 0.0;
        return std::pow((T_c - T) / T_c, beta_exponent);
    }

    static double susceptibility(double T, double T_c, double gamma_exponent) {
        return std::pow(std::abs(T - T_c), -gamma_exponent);
    }

    static double correlationLength(double T, double T_c, double nu_exponent) {
        return std::pow(std::abs(T - T_c), -nu_exponent);
    }
};

// ============================================================================
// ISING MODEL
// ============================================================================

class IsingModel {
private:
    int L_;
    std::vector<std::vector<int>> spins_;
    double J_;
    double h_;

public:
    IsingModel(int lattice_size, double coupling, double field)
        : L_(lattice_size), J_(coupling), h_(field) {
        spins_.resize(L_, std::vector<int>(L_, 1));
    }

    void randomize() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 1);

        for (int i = 0; i < L_; ++i) {
            for (int j = 0; j < L_; ++j) {
                spins_[i][j] = dis(gen) == 0 ? -1 : 1;
            }
        }
    }

    double energy() const {
        double E = 0.0;
        for (int i = 0; i < L_; ++i) {
            for (int j = 0; j < L_; ++j) {
                int s = spins_[i][j];
                int s_right = spins_[i][(j + 1) % L_];
                int s_down = spins_[(i + 1) % L_][j];
                E += -J_ * s * (s_right + s_down);
                E += -h_ * s;
            }
        }
        return E;
    }

    double magnetization() const {
        int sum = 0;
        for (int i = 0; i < L_; ++i) {
            for (int j = 0; j < L_; ++j) {
                sum += spins_[i][j];
            }
        }
        return static_cast<double>(sum) / (L_ * L_);
    }

    void metropolisStep(double T) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis_site(0, L_ - 1);
        std::uniform_real_distribution<> dis_prob(0.0, 1.0);

        int i = dis_site(gen);
        int j = dis_site(gen);

        int s = spins_[i][j];
        int s_left = spins_[i][(j - 1 + L_) % L_];
        int s_right = spins_[i][(j + 1) % L_];
        int s_up = spins_[(i - 1 + L_) % L_][j];
        int s_down = spins_[(i + 1) % L_][j];

        double dE = 2.0 * J_ * s * (s_left + s_right + s_up + s_down) + 2.0 * h_ * s;

        if (dE < 0 || dis_prob(gen) < std::exp(-dE / (k_B * T))) {
            spins_[i][j] = -s;
        }
    }

    void simulate(double T, int n_steps, int n_equilibrate = 1000) {
        for (int step = 0; step < n_equilibrate; ++step) {
            metropolisStep(T);
        }
    }

    double averageMagnetization(double T, int n_steps, int n_measure = 100) {
        simulate(T, n_steps / 2);

        double M_sum = 0.0;
        for (int m = 0; m < n_measure; ++m) {
            for (int s = 0; s < n_steps / n_measure; ++s) {
                metropolisStep(T);
            }
            M_sum += std::abs(magnetization());
        }
        return M_sum / n_measure;
    }

    double heatCapacity(double T, int n_steps, int n_measure = 100) {
        simulate(T, n_steps / 2);

        double E_sum = 0.0, E2_sum = 0.0;
        for (int m = 0; m < n_measure; ++m) {
            for (int s = 0; s < n_steps / n_measure; ++s) {
                metropolisStep(T);
            }
            double E = energy();
            E_sum += E;
            E2_sum += E * E;
        }

        double avg_E = E_sum / n_measure;
        double avg_E2 = E2_sum / n_measure;
        double var_E = avg_E2 - avg_E * avg_E;

        return var_E / (k_B * T * T);
    }

    double exactSolution1D(double T) const {
        double beta = 1.0 / (k_B * T);
        double Z_1site = 2.0 * std::cosh(beta * h_);
        return std::pow(Z_1site, L_);
    }
};

// ============================================================================
// CORRELATION FUNCTIONS
// ============================================================================

class CorrelationFunctions {
public:
    template<typename Observable>
    static double twoPointCorrelation(
        const std::vector<Observable>& measurements,
        int lag) {

        int N = measurements.size();
        if (lag >= N) return 0.0;

        double avg = std::accumulate(measurements.begin(), measurements.end(), 0.0) / N;

        double correlation = 0.0;
        for (int t = 0; t < N - lag; ++t) {
            correlation += (measurements[t] - avg) * (measurements[t + lag] - avg);
        }
        return correlation / (N - lag);
    }

    static double correlationLength(const std::vector<double>& correlations) {
        if (correlations.size() < 2) return 0.0;

        for (size_t i = 1; i < correlations.size(); ++i) {
            if (correlations[i] < correlations[0] / std::exp(1.0)) {
                return static_cast<double>(i);
            }
        }
        return correlations.size();
    }

    static std::vector<double> autocorrelationFunction(
        const std::vector<double>& time_series, int max_lag) {

        std::vector<double> acf;
        for (int lag = 0; lag < max_lag; ++lag) {
            acf.push_back(twoPointCorrelation(time_series, lag));
        }
        return acf;
    }
};

// ============================================================================
// FLUCTUATION-DISSIPATION
// ============================================================================

class FluctuationDissipation {
public:
    template<typename T>
    static double fluctuation(const std::vector<T>& measurements) {
        if (measurements.empty()) return 0.0;

        double avg = std::accumulate(measurements.begin(), measurements.end(), 0.0)
                    / measurements.size();

        double var = 0.0;
        for (const auto& m : measurements) {
            var += (m - avg) * (m - avg);
        }
        return std::sqrt(var / measurements.size());
    }

    static double responseFunction(double fluctuation_squared, double T) {
        return fluctuation_squared / (k_B * T);
    }

    static double einsteinRelation(double diffusion, double mobility, double T) {
        return std::abs(diffusion - k_B * T * mobility);
    }

    static double greenKubo(const std::vector<double>& current_autocorrelation,
                           double dt) {
        double integral = 0.0;
        for (size_t i = 0; i < current_autocorrelation.size(); ++i) {
            integral += current_autocorrelation[i] * dt;
        }
        return integral;
    }
};

// ============================================================================
// MEAN FIELD THEORY
// ============================================================================

class MeanFieldTheory {
public:
    static double meanFieldMagnetization(double T, double J, double z,
                                        double h = 0.0, int max_iter = 1000,
                                        double tol = 1e-8) {
        double m = 0.5;
        double beta = 1.0 / (k_B * T);

        for (int iter = 0; iter < max_iter; ++iter) {
            double m_new = std::tanh(beta * (z * J * m + h));
            if (std::abs(m_new - m) < tol) {
                return m_new;
            }
            m = m_new;
        }
        return m;
    }

    static double criticalTemperatureMF(double J, double z) {
        return z * J / k_B;
    }

    static double freeEnergyMF(double T, double m, double J, double z, double h) {
        double beta = 1.0 / (k_B * T);
        double arg = std::cosh(beta * (z * J * m + h));
        return -k_B * T * std::log(arg) + 0.5 * z * J * m * m;
    }
};

} // namespace statistical_mechanics
} // namespace physics

#endif // PHYSICS_STATISTICAL_MECHANICS_HPP
