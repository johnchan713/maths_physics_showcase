#ifndef NS_CASCADE_GALERKIN_HPP
#define NS_CASCADE_GALERKIN_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ns_cascade {

using Complex = std::complex<double>;

struct WaveVector {
    int x;
    int y;
    int z;

    WaveVector(int x_value = 0, int y_value = 0, int z_value = 0)
        : x(x_value), y(y_value), z(z_value) {}

    bool operator<(const WaveVector& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }

    WaveVector operator+(const WaveVector& other) const {
        return WaveVector(x + other.x, y + other.y, z + other.z);
    }

    WaveVector operator-() const { return WaveVector(-x, -y, -z); }
    int normSquared() const { return x * x + y * y + z * z; }
};

struct ComplexVector {
    Complex x;
    Complex y;
    Complex z;

    ComplexVector(Complex x_value = Complex(),
                  Complex y_value = Complex(),
                  Complex z_value = Complex())
        : x(x_value), y(y_value), z(z_value) {}

    ComplexVector& operator+=(const ComplexVector& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
};

inline ComplexVector operator+(ComplexVector left, const ComplexVector& right) {
    left += right;
    return left;
}

inline ComplexVector operator-(const ComplexVector& left, const ComplexVector& right) {
    return ComplexVector(left.x - right.x, left.y - right.y, left.z - right.z);
}

inline ComplexVector operator*(const ComplexVector& value, Complex scalar) {
    return ComplexVector(value.x * scalar, value.y * scalar, value.z * scalar);
}

inline ComplexVector operator*(Complex scalar, const ComplexVector& value) {
    return value * scalar;
}

inline ComplexVector operator*(const ComplexVector& value, double scalar) {
    return value * Complex(scalar, 0.0);
}

inline ComplexVector conjugate(const ComplexVector& value) {
    return ComplexVector(std::conj(value.x), std::conj(value.y), std::conj(value.z));
}

inline Complex dot(const WaveVector& wave, const ComplexVector& value) {
    return static_cast<double>(wave.x) * value.x +
           static_cast<double>(wave.y) * value.y +
           static_cast<double>(wave.z) * value.z;
}

inline Complex dot(const ComplexVector& value, const WaveVector& wave) {
    return dot(wave, value);
}

inline Complex innerProduct(const ComplexVector& left, const ComplexVector& right) {
    return std::conj(left.x) * right.x +
           std::conj(left.y) * right.y +
           std::conj(left.z) * right.z;
}

inline double normSquared(const ComplexVector& value) {
    return std::norm(value.x) + std::norm(value.y) + std::norm(value.z);
}

inline double norm(const ComplexVector& value) {
    return std::sqrt(normSquared(value));
}

inline ComplexVector cross(const WaveVector& wave, const ComplexVector& value) {
    return ComplexVector(
        static_cast<double>(wave.y) * value.z - static_cast<double>(wave.z) * value.y,
        static_cast<double>(wave.z) * value.x - static_cast<double>(wave.x) * value.z,
        static_cast<double>(wave.x) * value.y - static_cast<double>(wave.y) * value.x);
}

inline ComplexVector lerayProject(const WaveVector& wave, const ComplexVector& value) {
    const int wave_norm_squared = wave.normSquared();
    if (wave_norm_squared == 0) {
        throw std::invalid_argument("Leray projection is undefined for the zero mode");
    }

    const Complex parallel = dot(wave, value) /
                             static_cast<double>(wave_norm_squared);
    return ComplexVector(
        value.x - static_cast<double>(wave.x) * parallel,
        value.y - static_cast<double>(wave.y) * parallel,
        value.z - static_cast<double>(wave.z) * parallel);
}

class GalerkinSystem {
public:
    using State = std::vector<ComplexVector>;

    struct Diagnostics {
        double energy;
        double enstrophy;
        double palinstrophy;
        double critical_l3_sample;
        double sampled_vorticity_max;
        double vorticity_sup_upper_bound;
        double spectral_centroid;
        double high_shell_energy_fraction;
        double divergence_defect;
        double reality_defect;
        double energy_balance_residual;
    };

    GalerkinSystem(int cutoff, double viscosity)
        : cutoff_(cutoff), viscosity_(viscosity) {
        if (cutoff_ < 1) {
            throw std::invalid_argument("Fourier cutoff must be at least one");
        }
        if (viscosity_ < 0.0) {
            throw std::invalid_argument("Viscosity cannot be negative");
        }

        for (int kx = -cutoff_; kx <= cutoff_; ++kx) {
            for (int ky = -cutoff_; ky <= cutoff_; ++ky) {
                for (int kz = -cutoff_; kz <= cutoff_; ++kz) {
                    const WaveVector wave(kx, ky, kz);
                    if (wave.normSquared() == 0) continue;
                    index_[wave] = modes_.size();
                    modes_.push_back(wave);
                }
            }
        }
    }

    int cutoff() const { return cutoff_; }
    double viscosity() const { return viscosity_; }
    std::size_t modeCount() const { return modes_.size(); }
    const std::vector<WaveVector>& modes() const { return modes_; }

    State zeroState() const { return State(modes_.size()); }

    State deterministicLowModeState(double target_energy = 1.0) const {
        if (target_energy <= 0.0) {
            throw std::invalid_argument("Target energy must be positive");
        }

        State state = zeroState();
        for (std::size_t i = 0; i < modes_.size(); ++i) {
            const WaveVector& wave = modes_[i];
            if (!isCanonicalHalf(wave) || maximumComponent(wave) > 1) continue;

            const double phase = 0.73 * wave.x - 0.41 * wave.y + 0.29 * wave.z;
            const ComplexVector raw(
                Complex(std::sin(phase + 0.2), std::cos(phase - 0.1)),
                Complex(std::cos(phase + 0.7), std::sin(phase + 0.4)),
                Complex(std::sin(phase - 0.5), std::cos(phase + 0.9)));
            const ComplexVector projected = lerayProject(wave, raw);
            state[i] = projected;

            const std::map<WaveVector, std::size_t>::const_iterator negative =
                index_.find(-wave);
            if (negative != index_.end()) {
                state[negative->second] = conjugate(projected);
            }
        }

        const double current_energy = energy(state);
        if (current_energy == 0.0) {
            throw std::runtime_error("Failed to construct non-zero initial data");
        }
        const double scale = std::sqrt(target_energy / current_energy);
        for (std::size_t i = 0; i < state.size(); ++i) state[i] = state[i] * scale;
        return state;
    }

    State rightHandSide(const State& state, bool include_viscosity = true) const {
        requireCompatible(state);
        State result = zeroState();

        if (include_viscosity && viscosity_ != 0.0) {
            for (std::size_t i = 0; i < modes_.size(); ++i) {
                const double damping = -viscosity_ * modes_[i].normSquared();
                result[i] += state[i] * damping;
            }
        }

        const Complex minus_i(0.0, -1.0);
        for (std::size_t p_index = 0; p_index < modes_.size(); ++p_index) {
            const ComplexVector& u_p = state[p_index];
            if (normSquared(u_p) == 0.0) continue;

            for (std::size_t q_index = 0; q_index < modes_.size(); ++q_index) {
                const ComplexVector& u_q = state[q_index];
                if (normSquared(u_q) == 0.0) continue;

                const WaveVector output_wave = modes_[p_index] + modes_[q_index];
                const std::map<WaveVector, std::size_t>::const_iterator output =
                    index_.find(output_wave);
                if (output == index_.end()) continue;

                // Fourier coefficient of -P[(u . grad)u]:
                // -i P_k sum_{p+q=k} (u_p . q) u_q.
                const Complex advection = dot(u_p, modes_[q_index]);
                const ComplexVector interaction = u_q * (minus_i * advection);
                result[output->second] += lerayProject(output_wave, interaction);
            }
        }

        return result;
    }

    void stepRungeKutta4(State& state, double time_step) const {
        requireCompatible(state);
        if (time_step <= 0.0) {
            throw std::invalid_argument("Time step must be positive");
        }

        const State k1 = rightHandSide(state);
        const State k2 = rightHandSide(addScaled(state, k1, 0.5 * time_step));
        const State k3 = rightHandSide(addScaled(state, k2, 0.5 * time_step));
        const State k4 = rightHandSide(addScaled(state, k3, time_step));

        for (std::size_t i = 0; i < state.size(); ++i) {
            state[i] += (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) *
                        (time_step / 6.0);
        }
    }

    double energy(const State& state) const {
        requireCompatible(state);
        double value = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            value += 0.5 * normSquared(state[i]);
        }
        return value;
    }

    double enstrophy(const State& state) const {
        requireCompatible(state);
        double value = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            value += 0.5 * modes_[i].normSquared() * normSquared(state[i]);
        }
        return value;
    }

    double palinstrophy(const State& state) const {
        requireCompatible(state);
        double value = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            const double wave_squared = modes_[i].normSquared();
            value += 0.5 * wave_squared * wave_squared * normSquared(state[i]);
        }
        return value;
    }

    double energyDerivative(const State& state, const State& derivative) const {
        requireCompatible(state);
        requireCompatible(derivative);
        double value = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            value += std::real(innerProduct(state[i], derivative[i]));
        }
        return value;
    }

    Diagnostics diagnostics(const State& state, int sample_points_per_axis = 0) const {
        requireCompatible(state);
        if (sample_points_per_axis == 0) {
            sample_points_per_axis = std::max(2 * cutoff_ + 1, 7);
        }
        if (sample_points_per_axis < 2 * cutoff_ + 1) {
            throw std::invalid_argument(
                "Sampling grid must contain at least 2*cutoff+1 points per axis");
        }

        Diagnostics values;
        values.energy = energy(state);
        values.enstrophy = enstrophy(state);
        values.palinstrophy = palinstrophy(state);
        values.divergence_defect = divergenceDefect(state);
        values.reality_defect = realityDefect(state);
        values.vorticity_sup_upper_bound = vorticitySupremumUpperBound(state);
        values.spectral_centroid = spectralCentroid(state);
        values.high_shell_energy_fraction = highShellEnergyFraction(state);

        const std::pair<double, double> samples = sampledNorms(state,
                                                               sample_points_per_axis);
        values.critical_l3_sample = samples.first;
        values.sampled_vorticity_max = samples.second;

        const State derivative = rightHandSide(state);
        values.energy_balance_residual = std::abs(
            energyDerivative(state, derivative) + 2.0 * viscosity_ * values.enstrophy);
        return values;
    }

    double divergenceDefect(const State& state) const {
        requireCompatible(state);
        double defect = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            defect = std::max(defect, std::abs(dot(modes_[i], state[i])));
        }
        return defect;
    }

    double realityDefect(const State& state) const {
        requireCompatible(state);
        double defect = 0.0;
        for (std::size_t i = 0; i < modes_.size(); ++i) {
            const std::map<WaveVector, std::size_t>::const_iterator negative =
                index_.find(-modes_[i]);
            defect = std::max(defect,
                              norm(state[negative->second] - conjugate(state[i])));
        }
        return defect;
    }

private:
    int cutoff_;
    double viscosity_;
    std::vector<WaveVector> modes_;
    std::map<WaveVector, std::size_t> index_;

    static bool isCanonicalHalf(const WaveVector& wave) {
        if (wave.x != 0) return wave.x > 0;
        if (wave.y != 0) return wave.y > 0;
        return wave.z > 0;
    }

    static int maximumComponent(const WaveVector& wave) {
        return std::max(std::abs(wave.x),
                        std::max(std::abs(wave.y), std::abs(wave.z)));
    }

    void requireCompatible(const State& state) const {
        if (state.size() != modes_.size()) {
            throw std::invalid_argument("State size does not match Galerkin system");
        }
    }

    State addScaled(const State& state, const State& increment, double scale) const {
        State result = state;
        for (std::size_t i = 0; i < result.size(); ++i) {
            result[i] += increment[i] * scale;
        }
        return result;
    }

    double vorticitySupremumUpperBound(const State& state) const {
        double bound = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            bound += norm(cross(modes_[i], state[i]));
        }
        return bound;
    }

    double spectralCentroid(const State& state) const {
        double weighted_energy = 0.0;
        double total_energy = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            const double mode_energy = 0.5 * normSquared(state[i]);
            total_energy += mode_energy;
            weighted_energy += std::sqrt(static_cast<double>(modes_[i].normSquared())) *
                               mode_energy;
        }
        return total_energy == 0.0 ? 0.0 : weighted_energy / total_energy;
    }

    double highShellEnergyFraction(const State& state) const {
        double high_energy = 0.0;
        double total_energy = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            const double mode_energy = 0.5 * normSquared(state[i]);
            total_energy += mode_energy;
            if (maximumComponent(modes_[i]) == cutoff_) high_energy += mode_energy;
        }
        return total_energy == 0.0 ? 0.0 : high_energy / total_energy;
    }

    std::pair<double, double> sampledNorms(const State& state, int points) const {
        const double two_pi = 6.283185307179586476925286766559;
        double velocity_cube_sum = 0.0;
        double vorticity_max = 0.0;

        for (int ix = 0; ix < points; ++ix) {
            const double x = two_pi * ix / points;
            for (int iy = 0; iy < points; ++iy) {
                const double y = two_pi * iy / points;
                for (int iz = 0; iz < points; ++iz) {
                    const double z = two_pi * iz / points;
                    ComplexVector velocity;
                    ComplexVector vorticity;

                    for (std::size_t mode = 0; mode < modes_.size(); ++mode) {
                        const WaveVector& wave = modes_[mode];
                        const double phase = wave.x * x + wave.y * y + wave.z * z;
                        const Complex exponential(std::cos(phase), std::sin(phase));
                        velocity += state[mode] * exponential;
                        vorticity += cross(wave, state[mode]) *
                                     (Complex(0.0, 1.0) * exponential);
                    }

                    const double velocity_magnitude = norm(velocity);
                    velocity_cube_sum += velocity_magnitude * velocity_magnitude *
                                         velocity_magnitude;
                    vorticity_max = std::max(vorticity_max, norm(vorticity));
                }
            }
        }

        const double sample_count = static_cast<double>(points) * points * points;
        return std::make_pair(std::pow(velocity_cube_sum / sample_count, 1.0 / 3.0),
                              vorticity_max);
    }
};

}  // namespace ns_cascade

#endif
