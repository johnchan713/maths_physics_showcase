#ifndef NS_CASCADE_FFTW_REFERENCE_HPP
#define NS_CASCADE_FFTW_REFERENCE_HPP

#include "ns_cascade/galerkin.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

// Keep this small ABI declaration local to the optional reference backend. It
// lets the research harness use a system FFTW runtime even when development
// headers are unavailable, while CMake still requires the FFTW library before
// building any target that includes this file.
extern "C" {
typedef double ns_cascade_fftw_complex[2];
typedef struct fftw_plan_s* ns_cascade_fftw_plan;
void* fftw_malloc(std::size_t size);
void fftw_free(void* pointer);
ns_cascade_fftw_plan fftw_plan_dft_3d(int first_size,
                                      int second_size,
                                      int third_size,
                                      ns_cascade_fftw_complex* input,
                                      ns_cascade_fftw_complex* output,
                                      int sign,
                                      unsigned flags);
void fftw_execute(const ns_cascade_fftw_plan plan);
void fftw_destroy_plan(ns_cascade_fftw_plan plan);
}

namespace ns_cascade {

struct FftwReferenceDiagnostics {
    double energy;
    double enstrophy;
    double critical_h_half;
    double critical_l3_sample;
    double sampled_vorticity_max;
    double cutoff_shell_energy_fraction;
    double divergence_defect;
    double reality_defect;
};

struct FftwReferenceStepInfo {
    double time_step;
    double velocity_supremum_bound;
    double advective_cfl_upper_bound;
    double viscous_stability_number;
};

namespace fftw_reference_detail {

const int kForward = -1;
const int kBackward = 1;
const unsigned kEstimate = 1U << 6;

class Buffer {
public:
    explicit Buffer(std::size_t size)
        : data_(static_cast<ns_cascade_fftw_complex*>(
              fftw_malloc(sizeof(ns_cascade_fftw_complex) * size))) {
        if (data_ == NULL) throw std::bad_alloc();
    }

    ~Buffer() { fftw_free(data_); }

    ns_cascade_fftw_complex* data() { return data_; }

private:
    Buffer(const Buffer&);
    Buffer& operator=(const Buffer&);
    ns_cascade_fftw_complex* data_;
};

class Plan {
public:
    Plan(int grid_size,
         ns_cascade_fftw_complex* input,
         ns_cascade_fftw_complex* output,
         int sign)
        : plan_(fftw_plan_dft_3d(grid_size,
                                 grid_size,
                                 grid_size,
                                 input,
                                 output,
                                 sign,
                                 kEstimate)) {
        if (plan_ == NULL) {
            throw std::runtime_error("FFTW failed to create a 3D plan");
        }
    }

    ~Plan() { fftw_destroy_plan(plan_); }

    void execute() const { fftw_execute(plan_); }

private:
    Plan(const Plan&);
    Plan& operator=(const Plan&);
    ns_cascade_fftw_plan plan_;
};

}  // namespace fftw_reference_detail

// An intentionally separate time-evolution implementation. It shares the plain
// Fourier coefficient data types and elementary vector arithmetic with the
// main solver. FFTs are executed by FFTW, grid modes are rebuilt here, the
// Leray projection is written out independently, and RK4 stages are assembled
// independently.
//
// One instance owns mutable FFTW work buffers, so concurrent callers must use
// separate instances.
class FftwReferenceSystem {
public:
    using State = std::vector<ComplexVector>;

    FftwReferenceSystem(int grid_size,
                        double viscosity,
                        int requested_cutoff = 0)
        : grid_size_(grid_size),
          cutoff_(0),
          viscosity_(viscosity),
          grid_point_count_(0),
          input_(validatedGridPointCount(grid_size)),
          output_(validatedGridPointCount(grid_size)),
          forward_plan_(grid_size,
                        input_.data(),
                        output_.data(),
                        fftw_reference_detail::kForward),
          backward_plan_(grid_size,
                         input_.data(),
                         output_.data(),
                         fftw_reference_detail::kBackward) {
        if (!std::isfinite(viscosity_) || viscosity_ < 0.0) {
            throw std::invalid_argument(
                "Reference viscosity must be finite and non-negative");
        }
        const int safe_cutoff = (grid_size_ - 1) / 3;
        cutoff_ = requested_cutoff == 0 ? safe_cutoff : requested_cutoff;
        if (cutoff_ < 1 || cutoff_ > safe_cutoff) {
            throw std::invalid_argument(
                "Reference cutoff exceeds the strict 2/3 de-aliasing limit");
        }

        grid_point_count_ = validatedGridPointCount(grid_size_);
        modes_.resize(grid_point_count_);
        for (int x = 0; x < grid_size_; ++x) {
            for (int y = 0; y < grid_size_; ++y) {
                for (int z = 0; z < grid_size_; ++z) {
                    const std::size_t index = flatIndex(x, y, z);
                    const WaveVector wave(
                        waveNumber(x), waveNumber(y), waveNumber(z));
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
    double viscosity() const { return viscosity_; }
    std::size_t gridPointCount() const { return grid_point_count_; }
    const std::vector<WaveVector>& gridModes() const { return modes_; }

    State zeroState() const { return State(grid_point_count_); }

    std::vector<Complex> forwardTransform(
        const std::vector<Complex>& physical_values) const {
        return transform(physical_values, true);
    }

    std::vector<Complex> inverseTransform(
        const std::vector<Complex>& fourier_coefficients) const {
        return transform(fourier_coefficients, false);
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
            nonlinear_x[i] = velocity_y[i] * vorticity_z[i] -
                             velocity_z[i] * vorticity_y[i];
            nonlinear_y[i] = velocity_z[i] * vorticity_x[i] -
                             velocity_x[i] * vorticity_z[i];
            nonlinear_z[i] = velocity_x[i] * vorticity_y[i] -
                             velocity_y[i] * vorticity_x[i];
        }

        // These physical fields have no remaining readers. Releasing them
        // before the forward transforms/derivative reduces peak storage;
        // transform inputs and the independent projection are unchanged.
        std::vector<Complex>().swap(velocity_x);
        std::vector<Complex>().swap(velocity_y);
        std::vector<Complex>().swap(velocity_z);
        std::vector<Complex>().swap(vorticity_x);
        std::vector<Complex>().swap(vorticity_y);
        std::vector<Complex>().swap(vorticity_z);

        nonlinear_x = forwardTransform(nonlinear_x);
        nonlinear_y = forwardTransform(nonlinear_y);
        nonlinear_z = forwardTransform(nonlinear_z);

        State derivative = zeroState();
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const WaveVector& wave = modes_[index];
            const ComplexVector nonlinear(
                nonlinear_x[index], nonlinear_y[index], nonlinear_z[index]);

            // Written independently instead of calling lerayProject().
            const double wave_squared =
                static_cast<double>(wave.normSquared());
            const Complex parallel =
                (static_cast<double>(wave.x) * nonlinear.x +
                 static_cast<double>(wave.y) * nonlinear.y +
                 static_cast<double>(wave.z) * nonlinear.z) /
                wave_squared;
            derivative[index] = ComplexVector(
                nonlinear.x - static_cast<double>(wave.x) * parallel,
                nonlinear.y - static_cast<double>(wave.y) * parallel,
                nonlinear.z - static_cast<double>(wave.z) * parallel);
            if (include_viscosity && viscosity_ != 0.0) {
                derivative[index] += state[index] * (-viscosity_ * wave_squared);
            }
        }
        return derivative;
    }

    void stepRungeKutta4(State& state, double time_step) const {
        requireCompatible(state);
        if (!std::isfinite(time_step) || time_step <= 0.0) {
            throw std::invalid_argument(
                "Reference time step must be finite and positive");
        }

        const State first = rightHandSide(state);
        const State second = rightHandSide(
            addScaled(state, first, 0.5 * time_step));
        const State third = rightHandSide(
            addScaled(state, second, 0.5 * time_step));
        const State fourth = rightHandSide(
            addScaled(state, third, time_step));
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            state[index] +=
                (first[index] + second[index] * 2.0 +
                 third[index] * 2.0 + fourth[index]) *
                (time_step / 6.0);
        }
    }

    // Preserve the original RK4 arithmetic order while retaining only one
    // weighted derivative per active mode. The original four-derivative path
    // remains available to the oracle. This avoids a large peak allocation at
    // N=128 without changing the equation, timestep or Fourier truncation.
    void stepRungeKutta4Compact(State& state, double time_step) const {
        requireCompatible(state);
        if (!std::isfinite(time_step) || time_step <= 0.0) {
            throw std::invalid_argument("Reference time step must be finite and positive");
        }
        State stage = rightHandSide(state);
        std::vector<ComplexVector> accumulated(retained_indices_.size());
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            accumulated[n] = stage[index];
            stage[index] = state[index] + stage[index] * (0.5 * time_step);
        }
        for (int stage_number = 0; stage_number < 3; ++stage_number) {
            const State derivative = rightHandSide(stage);
            for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
                const std::size_t index = retained_indices_[n];
                if (stage_number == 2) {
                    state[index] += (accumulated[n] + derivative[index]) * (time_step / 6.0);
                } else {
                    accumulated[n] += derivative[index] * 2.0;
                    stage[index] = state[index] + derivative[index] *
                        (stage_number == 0 ? 0.5 * time_step : time_step);
                }
            }
        }
    }

    FftwReferenceStepInfo chooseAdaptiveTimeStep(
        const State& state,
        double maximum_time_step,
        double target_cfl = 0.4,
        double diffusion_safety = 2.0) const {
        requireCompatible(state);
        if (!std::isfinite(maximum_time_step) || maximum_time_step <= 0.0) {
            throw std::invalid_argument(
                "Reference maximum time step must be finite and positive");
        }
        if (!std::isfinite(target_cfl) || !std::isfinite(diffusion_safety) ||
            target_cfl <= 0.0 || diffusion_safety <= 0.0) {
            throw std::invalid_argument(
                "Reference CFL controls must be finite and positive");
        }

        FftwReferenceStepInfo information;
        information.velocity_supremum_bound = velocitySupremumUpperBound(state);
        if (!std::isfinite(information.velocity_supremum_bound)) {
            throw std::runtime_error(
                "Reference state has a non-finite velocity bound");
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
        information.time_step = std::min(
            maximum_time_step, std::min(advective_limit, viscous_limit));
        if (!std::isfinite(information.time_step) ||
            information.time_step <= 0.0) {
            throw std::runtime_error(
                "Reference adaptive time-step selection failed");
        }
        information.advective_cfl_upper_bound =
            information.time_step * information.velocity_supremum_bound *
            largest_wave_number;
        information.viscous_stability_number =
            information.time_step * viscosity_ * maximum_wave_squared;
        return information;
    }

    FftwReferenceDiagnostics diagnostics(const State& state) const {
        requireCompatible(state);
        FftwReferenceDiagnostics result;
        result.energy = 0.0;
        result.enstrophy = 0.0;
        double critical_norm_squared = 0.0;
        double cutoff_energy = 0.0;
        result.divergence_defect = 0.0;
        result.reality_defect = 0.0;

        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const WaveVector& wave = modes_[index];
            const double coefficient_energy = 0.5 * normSquared(state[index]);
            const double wave_squared =
                static_cast<double>(wave.normSquared());
            result.energy += coefficient_energy;
            result.enstrophy += wave_squared * coefficient_energy;
            critical_norm_squared +=
                std::sqrt(wave_squared) * normSquared(state[index]);
            if (maximumComponent(wave) == cutoff_) {
                cutoff_energy += coefficient_energy;
            }

            const Complex divergence =
                static_cast<double>(wave.x) * state[index].x +
                static_cast<double>(wave.y) * state[index].y +
                static_cast<double>(wave.z) * state[index].z;
            result.divergence_defect = std::max(
                result.divergence_defect, std::abs(divergence));
            const std::size_t negative = flatIndex(
                gridCoordinate(-wave.x),
                gridCoordinate(-wave.y),
                gridCoordinate(-wave.z));
            result.reality_defect = std::max(
                result.reality_defect,
                norm(state[negative] - conjugate(state[index])));
        }
        result.critical_h_half = std::sqrt(critical_norm_squared);
        result.cutoff_shell_energy_fraction =
            result.energy == 0.0 ? 0.0 : cutoff_energy / result.energy;

        const std::pair<double, double> sampled = sampledNorms(state);
        result.critical_l3_sample = sampled.first;
        result.sampled_vorticity_max = sampled.second;
        return result;
    }

private:
    int grid_size_;
    int cutoff_;
    double viscosity_;
    std::size_t grid_point_count_;
    std::vector<WaveVector> modes_;
    std::vector<std::size_t> retained_indices_;
    mutable fftw_reference_detail::Buffer input_;
    mutable fftw_reference_detail::Buffer output_;
    fftw_reference_detail::Plan forward_plan_;
    fftw_reference_detail::Plan backward_plan_;

    FftwReferenceSystem(const FftwReferenceSystem&);
    FftwReferenceSystem& operator=(const FftwReferenceSystem&);

    static std::size_t validatedGridPointCount(int grid_size) {
        if (grid_size < 8 || (grid_size & (grid_size - 1)) != 0) {
            throw std::invalid_argument(
                "Reference FFT grid must be a power of two and at least eight");
        }
        const std::size_t size = static_cast<std::size_t>(grid_size);
        if (size > std::numeric_limits<std::size_t>::max() / size ||
            size * size > std::numeric_limits<std::size_t>::max() / size) {
            throw std::overflow_error("Reference FFT grid is too large");
        }
        return size * size * size;
    }

    static int maximumComponent(const WaveVector& wave) {
        return std::max(std::abs(wave.x),
                        std::max(std::abs(wave.y), std::abs(wave.z)));
    }

    int waveNumber(int coordinate) const {
        return coordinate <= grid_size_ / 2
                   ? coordinate
                   : coordinate - grid_size_;
    }

    int gridCoordinate(int wave_number) const {
        return wave_number >= 0 ? wave_number : wave_number + grid_size_;
    }

    std::size_t flatIndex(int x, int y, int z) const {
        return (static_cast<std::size_t>(x) * grid_size_ +
                static_cast<std::size_t>(y)) *
                   grid_size_ +
               static_cast<std::size_t>(z);
    }

    void requireCompatible(const State& state) const {
        if (state.size() != grid_point_count_) {
            throw std::invalid_argument(
                "Reference state has the wrong grid size");
        }
    }

    void requireScalarGrid(const std::vector<Complex>& values) const {
        if (values.size() != grid_point_count_) {
            throw std::invalid_argument(
                "Reference scalar field has the wrong grid size");
        }
    }

    std::vector<Complex> transform(const std::vector<Complex>& values,
                                   bool forward) const {
        requireScalarGrid(values);
        for (std::size_t i = 0; i < grid_point_count_; ++i) {
            input_.data()[i][0] = values[i].real();
            input_.data()[i][1] = values[i].imag();
        }
        if (forward) {
            forward_plan_.execute();
        } else {
            backward_plan_.execute();
        }

        const double scale = forward
                                 ? 1.0 / static_cast<double>(grid_point_count_)
                                 : 1.0;
        std::vector<Complex> transformed(grid_point_count_);
        for (std::size_t i = 0; i < grid_point_count_; ++i) {
            transformed[i] =
                Complex(output_.data()[i][0], output_.data()[i][1]) * scale;
        }
        return transformed;
    }

    void fillSpectralFields(const State& state,
                            std::vector<Complex>& velocity_x,
                            std::vector<Complex>& velocity_y,
                            std::vector<Complex>& velocity_z,
                            std::vector<Complex>& vorticity_x,
                            std::vector<Complex>& vorticity_y,
                            std::vector<Complex>& vorticity_z) const {
        const Complex imaginary(0.0, 1.0);
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            const std::size_t index = retained_indices_[n];
            const WaveVector& wave = modes_[index];
            const ComplexVector& velocity = state[index];
            velocity_x[index] = velocity.x;
            velocity_y[index] = velocity.y;
            velocity_z[index] = velocity.z;

            // omega_k = i k cross u_k, written without the shared cross().
            vorticity_x[index] = imaginary *
                (static_cast<double>(wave.y) * velocity.z -
                 static_cast<double>(wave.z) * velocity.y);
            vorticity_y[index] = imaginary *
                (static_cast<double>(wave.z) * velocity.x -
                 static_cast<double>(wave.x) * velocity.z);
            vorticity_z[index] = imaginary *
                (static_cast<double>(wave.x) * velocity.y -
                 static_cast<double>(wave.y) * velocity.x);
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

    double velocitySupremumUpperBound(const State& state) const {
        double bound = 0.0;
        for (std::size_t n = 0; n < retained_indices_.size(); ++n) {
            bound += norm(state[retained_indices_[n]]);
        }
        return bound;
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
