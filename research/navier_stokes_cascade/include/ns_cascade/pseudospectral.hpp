#ifndef NS_CASCADE_PSEUDOSPECTRAL_HPP
#define NS_CASCADE_PSEUDOSPECTRAL_HPP

#include "ns_cascade/galerkin.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ns_cascade {

class PseudospectralSystem {
public:
    using State = std::vector<ComplexVector>;
    using Diagnostics = GalerkinSystem::Diagnostics;
    using ShellDiagnostics = GalerkinSystem::ShellDiagnostics;

    PseudospectralSystem(int grid_size,
                         double viscosity,
                         int requested_cutoff = 0)
        : grid_size_(grid_size),
          cutoff_(0),
          safe_dealias_cutoff_(0),
          viscosity_(viscosity),
          grid_point_count_(0) {
        if (grid_size_ < 8 || !isPowerOfTwo(grid_size_)) {
            throw std::invalid_argument(
                "FFT grid size must be a power of two and at least eight");
        }
        if (viscosity_ < 0.0) {
            throw std::invalid_argument("Viscosity cannot be negative");
        }

        // If |p_i|,|q_i| <= K and K < N/3, no wrapped quadratic product can
        // alias back into a retained component. (N-1)/3 enforces K < N/3.
        safe_dealias_cutoff_ = (grid_size_ - 1) / 3;
        cutoff_ = requested_cutoff == 0 ? safe_dealias_cutoff_ : requested_cutoff;
        if (cutoff_ < 1 || cutoff_ > safe_dealias_cutoff_) {
            throw std::invalid_argument(
                "Requested cutoff exceeds the strict 2/3 de-aliasing limit");
        }

        grid_point_count_ = static_cast<std::size_t>(grid_size_) *
                            static_cast<std::size_t>(grid_size_) *
                            static_cast<std::size_t>(grid_size_);
        modes_.resize(grid_point_count_);

        for (int ix = 0; ix < grid_size_; ++ix) {
            for (int iy = 0; iy < grid_size_; ++iy) {
                for (int iz = 0; iz < grid_size_; ++iz) {
                    const std::size_t index = flatIndex(ix, iy, iz);
                    const WaveVector wave(
                        waveNumber(ix), waveNumber(iy), waveNumber(iz));
                    modes_[index] = wave;
                    if (wave.normSquared() != 0 &&
                        maximumComponent(wave) <= cutoff_) {
                        retained_indices_.push_back(index);
                    }
                }
            }
        }
    }

    int gridSize() const { return grid_size_; }
    int cutoff() const { return cutoff_; }
    int safeDealiasCutoff() const { return safe_dealias_cutoff_; }
    double viscosity() const { return viscosity_; }
    std::size_t gridPointCount() const { return grid_point_count_; }
    std::size_t modeCount() const { return retained_indices_.size(); }

    State zeroState() const { return State(grid_point_count_); }

    State initialState(InitialCondition condition,
                       double target_energy = 1.0) const {
        const GalerkinSystem compact_system(cutoff_, viscosity_);
        const GalerkinSystem::State compact_state =
            compact_system.initialState(condition, target_energy);
        State state = zeroState();
        for (std::size_t i = 0; i < compact_system.modes().size(); ++i) {
            state[indexOf(compact_system.modes()[i])] = compact_state[i];
        }
        return state;
    }

    std::size_t indexOf(const WaveVector& wave) const {
        if (wave.normSquared() == 0 || maximumComponent(wave) > cutoff_) {
            throw std::out_of_range("Wave vector is not a retained FFT mode");
        }
        return flatIndex(gridCoordinate(wave.x),
                         gridCoordinate(wave.y),
                         gridCoordinate(wave.z));
    }

    std::vector<Complex> forwardTransform(
        const std::vector<Complex>& physical_values) const {
        requireScalarGrid(physical_values);
        std::vector<Complex> result = physical_values;
        transform3D(result, false);
        const double normalization = 1.0 / static_cast<double>(grid_point_count_);
        for (std::size_t i = 0; i < result.size(); ++i) {
            result[i] *= normalization;
        }
        return result;
    }

    std::vector<Complex> inverseTransform(
        const std::vector<Complex>& fourier_coefficients) const {
        requireScalarGrid(fourier_coefficients);
        std::vector<Complex> result = fourier_coefficients;
        transform3D(result, true);
        return result;
    }

    State rightHandSide(const State& state,
                        bool include_viscosity = true) const {
        requireCompatible(state);
        std::vector<Complex> velocity_x(grid_point_count_);
        std::vector<Complex> velocity_y(grid_point_count_);
        std::vector<Complex> velocity_z(grid_point_count_);
        std::vector<Complex> vorticity_x(grid_point_count_);
        std::vector<Complex> vorticity_y(grid_point_count_);
        std::vector<Complex> vorticity_z(grid_point_count_);
        fillSpectralFields(state,
                           velocity_x,
                           velocity_y,
                           velocity_z,
                           vorticity_x,
                           vorticity_y,
                           vorticity_z);

        velocity_x = inverseTransform(velocity_x);
        velocity_y = inverseTransform(velocity_y);
        velocity_z = inverseTransform(velocity_z);
        vorticity_x = inverseTransform(vorticity_x);
        vorticity_y = inverseTransform(vorticity_y);
        vorticity_z = inverseTransform(vorticity_z);

        std::vector<Complex> nonlinear_x(grid_point_count_);
        std::vector<Complex> nonlinear_y(grid_point_count_);
        std::vector<Complex> nonlinear_z(grid_point_count_);
        for (std::size_t i = 0; i < grid_point_count_; ++i) {
            // P[u x curl(u)] is exactly -P[(u . grad)u].
            nonlinear_x[i] =
                velocity_y[i] * vorticity_z[i] -
                velocity_z[i] * vorticity_y[i];
            nonlinear_y[i] =
                velocity_z[i] * vorticity_x[i] -
                velocity_x[i] * vorticity_z[i];
            nonlinear_z[i] =
                velocity_x[i] * vorticity_y[i] -
                velocity_y[i] * vorticity_x[i];
        }

        nonlinear_x = forwardTransform(nonlinear_x);
        nonlinear_y = forwardTransform(nonlinear_y);
        nonlinear_z = forwardTransform(nonlinear_z);

        State result = zeroState();
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const WaveVector& wave = modes_[index];
            const ComplexVector nonlinear(
                nonlinear_x[index], nonlinear_y[index], nonlinear_z[index]);
            result[index] = lerayProject(wave, nonlinear);
            if (include_viscosity && viscosity_ != 0.0) {
                result[index] +=
                    state[index] *
                    (-viscosity_ * static_cast<double>(wave.normSquared()));
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
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t i = retained_indices_[n];
            state[i] +=
                (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) *
                (time_step / 6.0);
        }
    }

    double energy(const State& state) const {
        requireCompatible(state);
        double value = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            value += 0.5 * normSquared(state[retained_indices_[n]]);
        }
        return value;
    }

    double enstrophy(const State& state) const {
        requireCompatible(state);
        double value = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            value += 0.5 * static_cast<double>(modes_[index].normSquared()) *
                     normSquared(state[index]);
        }
        return value;
    }

    double palinstrophy(const State& state) const {
        requireCompatible(state);
        double value = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const double wave_squared =
                static_cast<double>(modes_[index].normSquared());
            value += 0.5 * wave_squared * wave_squared * normSquared(state[index]);
        }
        return value;
    }

    double energyDerivative(const State& state,
                            const State& derivative) const {
        requireCompatible(state);
        requireCompatible(derivative);
        double value = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            value += std::real(innerProduct(state[index], derivative[index]));
        }
        return value;
    }

    Diagnostics diagnostics(const State& state,
                            int sample_points_per_axis = 0) const {
        requireCompatible(state);
        if (sample_points_per_axis != 0 &&
            sample_points_per_axis != grid_size_) {
            throw std::invalid_argument(
                "FFT diagnostics use the native grid size for spatial samples");
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

        const std::pair<double, double> samples = sampledNorms(state);
        values.critical_l3_sample = samples.first;
        values.sampled_vorticity_max = samples.second;
        const State derivative = rightHandSide(state);
        values.energy_balance_residual = std::abs(
            energyDerivative(state, derivative) +
            2.0 * viscosity_ * values.enstrophy);
        return values;
    }

    std::vector<ShellDiagnostics> shellDiagnostics(const State& state) const {
        requireCompatible(state);
        const State nonlinear_derivative = rightHandSide(state, false);
        const int shell_count = static_cast<int>(std::ceil(
            std::sqrt(3.0) * static_cast<double>(cutoff_) - 1e-12));
        std::vector<ShellDiagnostics> shells(
            static_cast<std::size_t>(shell_count));
        for (int shell = 1; shell <= shell_count; ++shell) {
            ShellDiagnostics& values = shells[static_cast<std::size_t>(shell - 1)];
            values.shell = shell;
            values.lower_radius = static_cast<double>(shell - 1);
            values.upper_radius = static_cast<double>(shell);
            values.energy = 0.0;
            values.nonlinear_transfer = 0.0;
            values.viscous_dissipation = 0.0;
            values.forward_flux = 0.0;
        }

        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const double wave_squared =
                static_cast<double>(modes_[index].normSquared());
            const int shell = std::max(
                1,
                static_cast<int>(std::ceil(std::sqrt(wave_squared) - 1e-12)));
            ShellDiagnostics& values =
                shells[static_cast<std::size_t>(shell - 1)];
            const double coefficient_energy = normSquared(state[index]);
            values.energy += 0.5 * coefficient_energy;
            values.nonlinear_transfer += std::real(
                innerProduct(state[index], nonlinear_derivative[index]));
            values.viscous_dissipation +=
                viscosity_ * wave_squared * coefficient_energy;
        }

        double cumulative_transfer = 0.0;
        for (std::size_t i = 0; i < shells.size(); ++i) {
            cumulative_transfer += shells[i].nonlinear_transfer;
            shells[i].forward_flux = -cumulative_transfer;
        }
        return shells;
    }

    double divergenceDefect(const State& state) const {
        requireCompatible(state);
        double defect = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            defect = std::max(
                defect, std::abs(dot(modes_[index], state[index])));
        }
        return defect;
    }

    double realityDefect(const State& state) const {
        requireCompatible(state);
        double defect = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const std::size_t negative = indexOf(-modes_[index]);
            defect = std::max(
                defect, norm(state[negative] - conjugate(state[index])));
        }
        return defect;
    }

private:
    int grid_size_;
    int cutoff_;
    int safe_dealias_cutoff_;
    double viscosity_;
    std::size_t grid_point_count_;
    std::vector<WaveVector> modes_;
    std::vector<std::size_t> retained_indices_;

    static bool isPowerOfTwo(int value) {
        return value > 0 && (value & (value - 1)) == 0;
    }

    static int maximumComponent(const WaveVector& wave) {
        return std::max(std::abs(wave.x),
                        std::max(std::abs(wave.y), std::abs(wave.z)));
    }

    int waveNumber(int coordinate) const {
        return coordinate <= grid_size_ / 2 ? coordinate : coordinate - grid_size_;
    }

    int gridCoordinate(int wave_number) const {
        return wave_number < 0 ? grid_size_ + wave_number : wave_number;
    }

    std::size_t flatIndex(int x, int y, int z) const {
        return (static_cast<std::size_t>(x) * grid_size_ +
                static_cast<std::size_t>(y)) *
                   grid_size_ +
               static_cast<std::size_t>(z);
    }

    void requireCompatible(const State& state) const {
        if (state.size() != grid_point_count_) {
            throw std::invalid_argument("State size does not match FFT grid");
        }
    }

    void requireScalarGrid(const std::vector<Complex>& values) const {
        if (values.size() != grid_point_count_) {
            throw std::invalid_argument("Scalar field size does not match FFT grid");
        }
    }

    static void transformLine(std::vector<Complex>& values, bool inverse) {
        const std::size_t size = values.size();
        for (std::size_t i = 1, j = 0; i < size; ++i) {
            std::size_t bit = size >> 1;
            while (j & bit) {
                j ^= bit;
                bit >>= 1;
            }
            j ^= bit;
            if (i < j) std::swap(values[i], values[j]);
        }

        const double two_pi = 6.283185307179586476925286766559;
        for (std::size_t length = 2; length <= size; length <<= 1) {
            const double angle =
                (inverse ? 1.0 : -1.0) * two_pi / static_cast<double>(length);
            const Complex root(std::cos(angle), std::sin(angle));
            for (std::size_t block = 0; block < size; block += length) {
                Complex weight(1.0, 0.0);
                for (std::size_t offset = 0; offset < length / 2; ++offset) {
                    const Complex even = values[block + offset];
                    const Complex odd =
                        values[block + offset + length / 2] * weight;
                    values[block + offset] = even + odd;
                    values[block + offset + length / 2] = even - odd;
                    weight *= root;
                }
            }
            if (length == size) break;
        }
    }

    void transform3D(std::vector<Complex>& values, bool inverse) const {
        std::vector<Complex> line(static_cast<std::size_t>(grid_size_));

        for (int x = 0; x < grid_size_; ++x) {
            for (int y = 0; y < grid_size_; ++y) {
                for (int z = 0; z < grid_size_; ++z) {
                    line[static_cast<std::size_t>(z)] =
                        values[flatIndex(x, y, z)];
                }
                transformLine(line, inverse);
                for (int z = 0; z < grid_size_; ++z) {
                    values[flatIndex(x, y, z)] = line[static_cast<std::size_t>(z)];
                }
            }
        }
        for (int x = 0; x < grid_size_; ++x) {
            for (int z = 0; z < grid_size_; ++z) {
                for (int y = 0; y < grid_size_; ++y) {
                    line[static_cast<std::size_t>(y)] =
                        values[flatIndex(x, y, z)];
                }
                transformLine(line, inverse);
                for (int y = 0; y < grid_size_; ++y) {
                    values[flatIndex(x, y, z)] = line[static_cast<std::size_t>(y)];
                }
            }
        }
        for (int y = 0; y < grid_size_; ++y) {
            for (int z = 0; z < grid_size_; ++z) {
                for (int x = 0; x < grid_size_; ++x) {
                    line[static_cast<std::size_t>(x)] =
                        values[flatIndex(x, y, z)];
                }
                transformLine(line, inverse);
                for (int x = 0; x < grid_size_; ++x) {
                    values[flatIndex(x, y, z)] = line[static_cast<std::size_t>(x)];
                }
            }
        }
    }

    void fillSpectralFields(const State& state,
                            std::vector<Complex>& velocity_x,
                            std::vector<Complex>& velocity_y,
                            std::vector<Complex>& velocity_z,
                            std::vector<Complex>& vorticity_x,
                            std::vector<Complex>& vorticity_y,
                            std::vector<Complex>& vorticity_z) const {
        const Complex imaginary_unit(0.0, 1.0);
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const ComplexVector& velocity = state[index];
            const ComplexVector vorticity =
                cross(modes_[index], velocity) * imaginary_unit;
            velocity_x[index] = velocity.x;
            velocity_y[index] = velocity.y;
            velocity_z[index] = velocity.z;
            vorticity_x[index] = vorticity.x;
            vorticity_y[index] = vorticity.y;
            vorticity_z[index] = vorticity.z;
        }
    }

    State addScaled(const State& state,
                    const State& increment,
                    double scale) const {
        State result = state;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            result[index] += increment[index] * scale;
        }
        return result;
    }

    double vorticitySupremumUpperBound(const State& state) const {
        double bound = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            bound += norm(cross(modes_[index], state[index]));
        }
        return bound;
    }

    double spectralCentroid(const State& state) const {
        double weighted_energy = 0.0;
        double total_energy = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const double mode_energy = 0.5 * normSquared(state[index]);
            total_energy += mode_energy;
            weighted_energy +=
                std::sqrt(static_cast<double>(modes_[index].normSquared())) *
                mode_energy;
        }
        return total_energy == 0.0 ? 0.0 : weighted_energy / total_energy;
    }

    double highShellEnergyFraction(const State& state) const {
        double high_energy = 0.0;
        double total_energy = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const double mode_energy = 0.5 * normSquared(state[index]);
            total_energy += mode_energy;
            if (maximumComponent(modes_[index]) == cutoff_) {
                high_energy += mode_energy;
            }
        }
        return total_energy == 0.0 ? 0.0 : high_energy / total_energy;
    }

    std::pair<double, double> sampledNorms(const State& state) const {
        std::vector<Complex> velocity_x(grid_point_count_);
        std::vector<Complex> velocity_y(grid_point_count_);
        std::vector<Complex> velocity_z(grid_point_count_);
        std::vector<Complex> vorticity_x(grid_point_count_);
        std::vector<Complex> vorticity_y(grid_point_count_);
        std::vector<Complex> vorticity_z(grid_point_count_);
        fillSpectralFields(state,
                           velocity_x,
                           velocity_y,
                           velocity_z,
                           vorticity_x,
                           vorticity_y,
                           vorticity_z);
        velocity_x = inverseTransform(velocity_x);
        velocity_y = inverseTransform(velocity_y);
        velocity_z = inverseTransform(velocity_z);
        vorticity_x = inverseTransform(vorticity_x);
        vorticity_y = inverseTransform(vorticity_y);
        vorticity_z = inverseTransform(vorticity_z);

        double velocity_cube_sum = 0.0;
        double vorticity_max = 0.0;
        for (std::size_t i = 0; i < grid_point_count_; ++i) {
            const double velocity_magnitude = std::sqrt(
                std::norm(velocity_x[i]) + std::norm(velocity_y[i]) +
                std::norm(velocity_z[i]));
            const double vorticity_magnitude = std::sqrt(
                std::norm(vorticity_x[i]) + std::norm(vorticity_y[i]) +
                std::norm(vorticity_z[i]));
            velocity_cube_sum += velocity_magnitude * velocity_magnitude *
                                 velocity_magnitude;
            vorticity_max = std::max(vorticity_max, vorticity_magnitude);
        }
        return std::make_pair(
            std::pow(velocity_cube_sum /
                         static_cast<double>(grid_point_count_),
                     1.0 / 3.0),
            vorticity_max);
    }
};

}  // namespace ns_cascade

#endif
