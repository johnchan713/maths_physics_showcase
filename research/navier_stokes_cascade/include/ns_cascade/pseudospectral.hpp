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

struct VortexTubeParameters {
    double core_radius;
    double separation;
    double bend_amplitude;
    int axial_wavenumber;

    VortexTubeParameters(double core_radius_value = 0.55,
                         double separation_value = 1.6,
                         double bend_amplitude_value = 0.25,
                         int axial_wavenumber_value = 1)
        : core_radius(core_radius_value),
          separation(separation_value),
          bend_amplitude(bend_amplitude_value),
          axial_wavenumber(axial_wavenumber_value) {}
};

struct VortexBundleParameters {
    VortexTubeParameters tubes;
    double orthogonal_pair_weight;
    double phase_offset;

    VortexBundleParameters(
        const VortexTubeParameters& tube_values = VortexTubeParameters(),
        double orthogonal_pair_weight_value = 0.75,
        double phase_offset_value = 1.0471975511965977461542144610932)
        : tubes(tube_values),
          orthogonal_pair_weight(orthogonal_pair_weight_value),
          phase_offset(phase_offset_value) {}
};

struct AdaptiveStepInfo {
    double time_step;
    double velocity_supremum_bound;
    double advective_cfl_upper_bound;
    double viscous_stability_number;
};

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
        if (!std::isfinite(viscosity_) || viscosity_ < 0.0) {
            throw std::invalid_argument("Viscosity must be finite and non-negative");
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
    const std::vector<WaveVector>& gridModes() const { return modes_; }

    State zeroState() const { return State(grid_point_count_); }

    State initialState(InitialCondition condition,
                       double target_energy = 1.0) const {
        if (condition == InitialCondition::VortexTubes) {
            return vortexTubePairState(VortexTubeParameters(), target_energy);
        }
        const GalerkinSystem compact_system(cutoff_, viscosity_);
        const GalerkinSystem::State compact_state =
            compact_system.initialState(condition, target_energy);
        State state = zeroState();
        for (std::size_t i = 0; i < compact_system.modes().size(); ++i) {
            state[indexOf(compact_system.modes()[i])] = compact_state[i];
        }
        return state;
    }

    State vortexTubePairState(const VortexTubeParameters& parameters,
                              double target_energy = 1.0) const {
        validateVortexTubeParameters(parameters);
        if (!std::isfinite(target_energy) || target_energy <= 0.0) {
            throw std::invalid_argument("Target energy must be finite and positive");
        }

        std::vector<Complex> physical_x(grid_point_count_);
        std::vector<Complex> physical_y(grid_point_count_);
        std::vector<Complex> physical_z(grid_point_count_);
        const double two_pi = 6.283185307179586476925286766559;
        const double inverse_core_squared =
            1.0 / (parameters.core_radius * parameters.core_radius);

        for (int ix = 0; ix < grid_size_; ++ix) {
            const double x = two_pi * ix / grid_size_ - 0.5 * two_pi;
            for (int iy = 0; iy < grid_size_; ++iy) {
                const double y = two_pi * iy / grid_size_ - 0.5 * two_pi;
                for (int iz = 0; iz < grid_size_; ++iz) {
                    const double z = two_pi * iz / grid_size_;
                    const double bend_phase = parameters.axial_wavenumber * z;
                    const double bend_x =
                        parameters.bend_amplitude * std::cos(bend_phase);
                    const double bend_y =
                        parameters.bend_amplitude * std::sin(bend_phase);
                    const double center_1_x =
                        -0.5 * parameters.separation + bend_x;
                    const double center_1_y = bend_y;
                    const double center_2_x =
                        0.5 * parameters.separation - bend_x;
                    const double center_2_y = -bend_y;

                    const double delta_1_x = x - center_1_x;
                    const double delta_1_y = y - center_1_y;
                    const double delta_2_x = x - center_2_x;
                    const double delta_2_y = y - center_2_y;
                    const double distance_1_squared =
                        2.0 * (1.0 - std::cos(delta_1_x)) +
                        2.0 * (1.0 - std::cos(delta_1_y));
                    const double distance_2_squared =
                        2.0 * (1.0 - std::cos(delta_2_x)) +
                        2.0 * (1.0 - std::cos(delta_2_y));
                    const double gaussian_1 = std::exp(
                        -0.5 * inverse_core_squared * distance_1_squared);
                    const double gaussian_2 = std::exp(
                        -0.5 * inverse_core_squared * distance_2_squared);
                    const double derivative_1_x =
                        -inverse_core_squared * std::sin(delta_1_x) * gaussian_1;
                    const double derivative_1_y =
                        -inverse_core_squared * std::sin(delta_1_y) * gaussian_1;
                    const double derivative_2_x =
                        -inverse_core_squared * std::sin(delta_2_x) * gaussian_2;
                    const double derivative_2_y =
                        -inverse_core_squared * std::sin(delta_2_y) * gaussian_2;
                    const std::size_t index = flatIndex(ix, iy, iz);

                    // A = (0,0,G1-G2), u = curl(A). The opposite signs form
                    // a counter-rotating pair; helical centre lines make it 3D.
                    physical_x[index] = derivative_1_y - derivative_2_y;
                    physical_y[index] = -(derivative_1_x - derivative_2_x);
                    physical_z[index] = Complex();
                }
            }
        }

        const std::vector<Complex> fourier_x = forwardTransform(physical_x);
        const std::vector<Complex> fourier_y = forwardTransform(physical_y);
        const std::vector<Complex> fourier_z = forwardTransform(physical_z);
        State state = zeroState();
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            state[index] = lerayProject(
                modes_[index],
                ComplexVector(fourier_x[index],
                              fourier_y[index],
                              fourier_z[index]));
        }

        const double current_energy = energy(state);
        if (current_energy == 0.0) {
            throw std::runtime_error("Vortex-tube construction produced zero energy");
        }
        const double scale = std::sqrt(target_energy / current_energy);
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            state[retained_indices_[n]] = state[retained_indices_[n]] * scale;
        }
        return state;
    }

    State vortexBundleState(const VortexBundleParameters& parameters,
                            double target_energy = 1.0) const {
        validateVortexTubeParameters(parameters.tubes);
        if (!std::isfinite(parameters.orthogonal_pair_weight) ||
            parameters.orthogonal_pair_weight < 0.0) {
            throw std::invalid_argument(
                "Orthogonal-pair weight must be finite and non-negative");
        }
        if (!std::isfinite(parameters.phase_offset)) {
            throw std::invalid_argument("Bundle phase offset must be finite");
        }
        if (!std::isfinite(target_energy) || target_energy <= 0.0) {
            throw std::invalid_argument("Target energy must be finite and positive");
        }

        std::vector<Complex> physical_x(grid_point_count_);
        std::vector<Complex> physical_y(grid_point_count_);
        std::vector<Complex> physical_z(grid_point_count_);
        const double two_pi = 6.283185307179586476925286766559;

        for (int ix = 0; ix < grid_size_; ++ix) {
            const double x = two_pi * ix / grid_size_ - 0.5 * two_pi;
            for (int iy = 0; iy < grid_size_; ++iy) {
                const double y = two_pi * iy / grid_size_ - 0.5 * two_pi;
                for (int iz = 0; iz < grid_size_; ++iz) {
                    const double z = two_pi * iz / grid_size_;
                    const double coordinates[3] = {x, y, z};
                    double velocity[3] = {0.0, 0.0, 0.0};

                    // Each contribution is curl(A) for a vector potential
                    // aligned with its tube axis. The z pair is the original
                    // family; the x and y pairs introduce orthogonal contacts.
                    addVortexPairCurl(
                        2, coordinates, parameters.tubes, 0.0, 1.0, velocity);
                    addVortexPairCurl(
                        0,
                        coordinates,
                        parameters.tubes,
                        parameters.phase_offset,
                        parameters.orthogonal_pair_weight,
                        velocity);
                    addVortexPairCurl(
                        1,
                        coordinates,
                        parameters.tubes,
                        -parameters.phase_offset,
                        parameters.orthogonal_pair_weight,
                        velocity);

                    const std::size_t index = flatIndex(ix, iy, iz);
                    physical_x[index] = Complex(velocity[0], 0.0);
                    physical_y[index] = Complex(velocity[1], 0.0);
                    physical_z[index] = Complex(velocity[2], 0.0);
                }
            }
        }

        const std::vector<Complex> fourier_x = forwardTransform(physical_x);
        const std::vector<Complex> fourier_y = forwardTransform(physical_y);
        const std::vector<Complex> fourier_z = forwardTransform(physical_z);
        State state = zeroState();
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            state[index] = lerayProject(
                modes_[index],
                ComplexVector(fourier_x[index],
                              fourier_y[index],
                              fourier_z[index]));
        }

        const double current_energy = energy(state);
        if (current_energy == 0.0) {
            throw std::runtime_error("Vortex-bundle construction produced zero energy");
        }
        const double scale = std::sqrt(target_energy / current_energy);
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            state[retained_indices_[n]] = state[retained_indices_[n]] * scale;
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
        if (!std::isfinite(time_step) || time_step <= 0.0) {
            throw std::invalid_argument("Time step must be finite and positive");
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

    double velocitySupremumUpperBound(const State& state) const {
        requireCompatible(state);
        double bound = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            bound += norm(state[retained_indices_[n]]);
        }
        return bound;
    }

    AdaptiveStepInfo chooseAdaptiveTimeStep(
        const State& state,
        double maximum_time_step,
        double target_cfl = 0.4,
        double diffusion_safety = 2.0) const {
        requireCompatible(state);
        if (!std::isfinite(maximum_time_step) || maximum_time_step <= 0.0) {
            throw std::invalid_argument(
                "Maximum time step must be finite and positive");
        }
        if (!std::isfinite(target_cfl) || !std::isfinite(diffusion_safety) ||
            target_cfl <= 0.0 || diffusion_safety <= 0.0) {
            throw std::invalid_argument(
                "CFL and diffusion safety values must be positive");
        }

        AdaptiveStepInfo information;
        information.velocity_supremum_bound =
            velocitySupremumUpperBound(state);
        if (!std::isfinite(information.velocity_supremum_bound)) {
            throw std::runtime_error(
                "Cannot choose a time step for a non-finite velocity state");
        }
        const double largest_wave_number =
            std::sqrt(3.0) * static_cast<double>(cutoff_);
        const double advective_limit =
            information.velocity_supremum_bound == 0.0
                ? maximum_time_step
                : target_cfl /
                      (information.velocity_supremum_bound * largest_wave_number);
        const double maximum_wave_squared =
            3.0 * static_cast<double>(cutoff_) * cutoff_;
        const double viscous_limit =
            viscosity_ == 0.0
                ? maximum_time_step
                : diffusion_safety / (viscosity_ * maximum_wave_squared);
        information.time_step =
            std::min(maximum_time_step,
                     std::min(advective_limit, viscous_limit));
        if (!std::isfinite(information.time_step) ||
            information.time_step <= 0.0) {
            throw std::runtime_error("Adaptive time-step selection failed");
        }
        information.advective_cfl_upper_bound =
            information.time_step * information.velocity_supremum_bound *
            largest_wave_number;
        information.viscous_stability_number =
            information.time_step * viscosity_ * maximum_wave_squared;
        return information;
    }

    AdaptiveStepInfo stepAdaptiveRungeKutta4(
        State& state,
        double maximum_time_step,
        double target_cfl = 0.4,
        double diffusion_safety = 2.0) const {
        const AdaptiveStepInfo information = chooseAdaptiveTimeStep(
            state, maximum_time_step, target_cfl, diffusion_safety);
        stepRungeKutta4(state, information.time_step);
        return information;
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

    double cutoffShellEnergyFraction(const State& state) const {
        requireCompatible(state);
        return highShellEnergyFraction(state);
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

    double criticalHOneHalf(const State& state) const {
        requireCompatible(state);
        double norm_squared = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const double wave_magnitude = std::sqrt(
                static_cast<double>(modes_[index].normSquared()));
            norm_squared += wave_magnitude * normSquared(state[index]);
        }
        return std::sqrt(norm_squared);
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

    double enstrophyDerivative(const State& state,
                               const State& derivative) const {
        requireCompatible(state);
        requireCompatible(derivative);
        double value = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            value += static_cast<double>(modes_[index].normSquared()) *
                     std::real(innerProduct(state[index], derivative[index]));
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
        values.critical_h_half = criticalHOneHalf(state);
        values.divergence_defect = divergenceDefect(state);
        values.reality_defect = realityDefect(state);
        values.vorticity_sup_upper_bound = vorticitySupremumUpperBound(state);
        values.spectral_centroid = spectralCentroid(state);
        values.high_shell_energy_fraction = highShellEnergyFraction(state);

        const std::pair<double, double> samples = sampledNorms(state);
        values.critical_l3_sample = samples.first;
        values.sampled_vorticity_max = samples.second;
        const State nonlinear_derivative = rightHandSide(state, false);
        values.nonlinear_enstrophy_production =
            enstrophyDerivative(state, nonlinear_derivative);
        values.viscous_enstrophy_destruction =
            2.0 * viscosity_ * values.palinstrophy;
        values.net_enstrophy_rate =
            values.nonlinear_enstrophy_production -
            values.viscous_enstrophy_destruction;
        values.enstrophy_production_to_dissipation =
            values.viscous_enstrophy_destruction == 0.0
                ? 0.0
                : values.nonlinear_enstrophy_production /
                      values.viscous_enstrophy_destruction;
        values.energy_balance_residual =
            std::abs(energyDerivative(state, nonlinear_derivative));
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

    static void addVortexPairCurl(
        int axis,
        const double coordinates[3],
        const VortexTubeParameters& parameters,
        double phase_offset,
        double weight,
        double velocity[3]) {
        if (weight == 0.0) return;

        // For cyclic coordinates (axis,b,c), A_axis = G1-G2 gives
        // curl(A)_b = d_c A_axis and curl(A)_c = -d_b A_axis.
        const int transverse_b = (axis + 1) % 3;
        const int transverse_c = (axis + 2) % 3;
        const double inverse_core_squared =
            1.0 / (parameters.core_radius * parameters.core_radius);
        const double bend_phase =
            parameters.axial_wavenumber * coordinates[axis] + phase_offset;
        const double bend_b =
            parameters.bend_amplitude * std::cos(bend_phase);
        const double bend_c =
            parameters.bend_amplitude * std::sin(bend_phase);
        const double center_1_b = -0.5 * parameters.separation + bend_b;
        const double center_1_c = bend_c;
        const double center_2_b = 0.5 * parameters.separation - bend_b;
        const double center_2_c = -bend_c;

        const double delta_1_b = coordinates[transverse_b] - center_1_b;
        const double delta_1_c = coordinates[transverse_c] - center_1_c;
        const double delta_2_b = coordinates[transverse_b] - center_2_b;
        const double delta_2_c = coordinates[transverse_c] - center_2_c;
        const double distance_1_squared =
            2.0 * (1.0 - std::cos(delta_1_b)) +
            2.0 * (1.0 - std::cos(delta_1_c));
        const double distance_2_squared =
            2.0 * (1.0 - std::cos(delta_2_b)) +
            2.0 * (1.0 - std::cos(delta_2_c));
        const double gaussian_1 = std::exp(
            -0.5 * inverse_core_squared * distance_1_squared);
        const double gaussian_2 = std::exp(
            -0.5 * inverse_core_squared * distance_2_squared);
        const double derivative_1_b =
            -inverse_core_squared * std::sin(delta_1_b) * gaussian_1;
        const double derivative_1_c =
            -inverse_core_squared * std::sin(delta_1_c) * gaussian_1;
        const double derivative_2_b =
            -inverse_core_squared * std::sin(delta_2_b) * gaussian_2;
        const double derivative_2_c =
            -inverse_core_squared * std::sin(delta_2_c) * gaussian_2;

        velocity[transverse_b] +=
            weight * (derivative_1_c - derivative_2_c);
        velocity[transverse_c] -=
            weight * (derivative_1_b - derivative_2_b);
    }

    void validateVortexTubeParameters(
        const VortexTubeParameters& parameters) const {
        const double pi = 3.1415926535897932384626433832795;
        const double two_pi = 2.0 * pi;
        if (!std::isfinite(parameters.core_radius) ||
            parameters.core_radius <= 0.0 || parameters.core_radius > pi) {
            throw std::invalid_argument(
                "Vortex-tube core radius must be finite and in (0, pi]");
        }
        if (!std::isfinite(parameters.separation) ||
            parameters.separation <= 0.0 || parameters.separation >= two_pi) {
            throw std::invalid_argument(
                "Vortex-tube separation must be finite and in (0, 2*pi)");
        }
        if (!std::isfinite(parameters.bend_amplitude) ||
            parameters.bend_amplitude < 0.0 ||
            parameters.bend_amplitude >= pi) {
            throw std::invalid_argument(
                "Vortex-tube bend amplitude must be finite and in [0, pi)");
        }
        if (parameters.axial_wavenumber < 1 ||
            parameters.axial_wavenumber > cutoff_) {
            throw std::invalid_argument(
                "Vortex-tube axial wavenumber must be between 1 and the cutoff");
        }
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
