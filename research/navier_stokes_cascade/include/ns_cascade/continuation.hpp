#ifndef NS_CASCADE_CONTINUATION_HPP
#define NS_CASCADE_CONTINUATION_HPP

#include "ns_cascade/checkpoint.hpp"
#include "ns_cascade/fftw_reference.hpp"
#include "ns_cascade/optimization_state_csv.hpp"

#include <functional>
#include <memory>

namespace ns_cascade {

// This format is deliberately separate from the benchmark-initializer restart
// format. It evolves the saved Fourier coefficients without reconstructing or
// normalizing them. The backend and the physical observation clock are frozen.
struct ContinuationConfiguration {
    int grid = 32;
    int cutoff = 10;
    int sampling_grid = 64;
    int dense_sampling_grid = 128;
    bool fftw = true;
    double viscosity = 0.02;
    double maximum_time_step = 0.000125;
    double target_cfl = 0.4;
    double diffusion_safety = 2.0;
    double observation_interval = 0.01;
    double cutoff_limit = 0.008;
};

struct ContinuationObservation {
    double time = 0.0;
    double energy = 0.0;
    double enstrophy = 0.0;
    double palinstrophy = 0.0;
    double h_half = 0.0;
    double l3_sample = 0.0;
    double vorticity_sample = 0.0;
    double l3_dense = 0.0;
    double vorticity_dense = 0.0;
    double vorticity_fourier_bound = 0.0;
    double k_rms = 0.0;
    double stretching = 0.0;
    double viscous_destruction = 0.0;
    double net_enstrophy_rate = 0.0;
    double production_to_dissipation = 0.0;
    double cutoff_fraction = 0.0;
    double divergence_defect = 0.0;
    double reality_defect = 0.0;
    double nonlinear_energy_residual = 0.0;
};

struct ContinuationCheckpoint {
    ContinuationConfiguration configuration;
    std::uint64_t step = 0;
    std::uint64_t completed_observations = 0;
    // A cutoff failure is retained as evidence and cannot be restarted.
    bool cutoff_stopped = false;
    double peak_cutoff_fraction = 0.0;
    double minimum_time_step = 0.0;
    double maximum_time_step = 0.0;
    double maximum_cfl = 0.0;
    double maximum_viscous_number = 0.0;
    double bkm_sampled_integral = 0.0;
    ContinuationObservation initial;
    ContinuationObservation latest;
    OptimizationState state;
};

namespace continuation_detail {

inline void validateConfiguration(const ContinuationConfiguration& c) {
    checkpoint_detail::expectedStateSize(c.grid);
    checkpoint_detail::expectedStateSize(c.sampling_grid);
    checkpoint_detail::expectedStateSize(c.dense_sampling_grid);
    if (c.cutoff < 1 || c.cutoff > (c.grid - 1) / 3 ||
        c.sampling_grid < c.grid || c.dense_sampling_grid < c.sampling_grid) {
        throw std::invalid_argument("Continuation grids/cutoff would drop modes");
    }
    const double values[] = {c.viscosity, c.maximum_time_step, c.target_cfl,
        c.diffusion_safety, c.observation_interval, c.cutoff_limit};
    for (double value : values) {
        if (!std::isfinite(value) || value <= 0.0) {
            throw std::invalid_argument("Continuation settings must be positive and finite");
        }
    }
    if (c.target_cfl > 0.4 || c.diffusion_safety > 2.0 || c.cutoff_limit >= 1.0) {
        throw std::invalid_argument("Continuation settings exceed supported safety limits");
    }
}

inline bool sameConfiguration(const ContinuationConfiguration& a,
                               const ContinuationConfiguration& b) {
    return a.grid == b.grid && a.cutoff == b.cutoff &&
        a.sampling_grid == b.sampling_grid && a.dense_sampling_grid == b.dense_sampling_grid &&
        a.fftw == b.fftw && a.viscosity == b.viscosity &&
        a.maximum_time_step == b.maximum_time_step && a.target_cfl == b.target_cfl &&
        a.diffusion_safety == b.diffusion_safety && a.observation_interval == b.observation_interval &&
        a.cutoff_limit == b.cutoff_limit;
}

inline std::uint64_t observationIndex(double time, double interval) {
    const double index = time / interval;
    if (!std::isfinite(index) || index < 0.0 || index > 1e12 ||
        std::abs(index - std::round(index)) >
            64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, index)) {
        throw std::invalid_argument("Final time must lie on the frozen observation clock");
    }
    return static_cast<std::uint64_t>(std::round(index));
}

inline std::vector<double> observationValues(const ContinuationObservation& d) {
    return {d.time, d.energy, d.enstrophy, d.palinstrophy, d.h_half,
        d.l3_sample, d.vorticity_sample, d.l3_dense, d.vorticity_dense,
        d.vorticity_fourier_bound, d.k_rms, d.stretching,
        d.viscous_destruction, d.net_enstrophy_rate, d.production_to_dissipation,
        d.cutoff_fraction, d.divergence_defect, d.reality_defect,
        d.nonlinear_energy_residual};
}

inline void appendObservation(std::vector<unsigned char>& bytes,
                              const ContinuationObservation& d) {
    for (double value : observationValues(d)) checkpoint_detail::appendDouble(bytes, value);
}

inline ContinuationObservation readObservation(const std::vector<unsigned char>& bytes,
                                                std::size_t& offset) {
    ContinuationObservation d;
    double* values[] = {&d.time, &d.energy, &d.enstrophy, &d.palinstrophy, &d.h_half,
        &d.l3_sample, &d.vorticity_sample, &d.l3_dense, &d.vorticity_dense,
        &d.vorticity_fourier_bound, &d.k_rms, &d.stretching,
        &d.viscous_destruction, &d.net_enstrophy_rate, &d.production_to_dissipation,
        &d.cutoff_fraction, &d.divergence_defect, &d.reality_defect,
        &d.nonlinear_energy_residual};
    for (double* value : values) *value = checkpoint_detail::readDouble(bytes, offset);
    return d;
}

inline void validateCheckpoint(const ContinuationCheckpoint& p) {
    validateConfiguration(p.configuration);
    if (p.state.size() != checkpoint_detail::expectedStateSize(p.configuration.grid)) {
        throw std::runtime_error("Continuation checkpoint has wrong state size");
    }
    for (const auto& d : {p.initial, p.latest}) {
        for (double value : observationValues(d)) {
            if (!std::isfinite(value)) throw std::runtime_error("Nonfinite continuation observation");
        }
        if (d.time < 0.0 || d.energy <= 0.0 || d.h_half <= 0.0 ||
            d.l3_sample <= 0.0 || d.vorticity_sample <= 0.0 ||
            d.l3_dense <= 0.0 || d.vorticity_dense <= 0.0 ||
            d.cutoff_fraction < 0.0 || d.cutoff_fraction > 1.0) {
            throw std::runtime_error("Invalid continuation observation");
        }
    }
    if (p.initial.time != 0.0 || p.completed_observations > UINT64_C(1000000000000)) {
        throw std::runtime_error("Invalid continuation clock origin");
    }
    const double clock_time = p.completed_observations * p.configuration.observation_interval;
    if ((!p.cutoff_stopped && p.latest.time != clock_time) ||
        (p.cutoff_stopped && (p.latest.time < clock_time ||
          p.latest.time > (p.completed_observations + 1) * p.configuration.observation_interval))) {
        throw std::runtime_error("Continuation checkpoint clock is inconsistent");
    }
    const double progress[] = {p.peak_cutoff_fraction, p.minimum_time_step,
        p.maximum_time_step, p.maximum_cfl, p.maximum_viscous_number, p.bkm_sampled_integral};
    for (double value : progress) {
        if (!std::isfinite(value) || value < 0.0) {
            throw std::runtime_error("Invalid continuation progress");
        }
    }
    if (p.peak_cutoff_fraction < p.latest.cutoff_fraction ||
        p.cutoff_stopped != (p.peak_cutoff_fraction > p.configuration.cutoff_limit)) {
        throw std::runtime_error("Continuation cutoff gate is inconsistent");
    }
    const PseudospectralSystem system(p.configuration.grid, p.configuration.viscosity,
                                      p.configuration.cutoff);
    requireOptimizationState(system, p.state);
    const auto& modes = system.gridModes();
    for (std::size_t i = 0; i < p.state.size(); ++i) {
        const auto& value = p.state[i];
        for (double part : {value.x.real(), value.x.imag(), value.y.real(), value.y.imag(), value.z.real(), value.z.imag()}) {
            if (!std::isfinite(part)) throw std::runtime_error("Nonfinite continuation coefficient");
        }
        if ((modes[i].normSquared() == 0 || stateMaximumComponent(modes[i]) > p.configuration.cutoff) &&
            normSquared(value) != 0.0) throw std::runtime_error("Continuation contains an unretained mode");
    }
    if (system.divergenceDefect(p.state) > 1e-9 || system.realityDefect(p.state) > 1e-9) {
        throw std::runtime_error("Invalid continuation Fourier constraints");
    }
    if (std::abs(system.energy(p.state) - p.latest.energy) > 1e-12 * p.latest.energy) {
        throw std::runtime_error("Checkpoint state does not match its energy observation");
    }
}

}  // namespace continuation_detail

class ContinuationSystem {
public:
    explicit ContinuationSystem(const ContinuationConfiguration& configuration)
        : c_(configuration), system_(c_.grid, c_.viscosity, c_.cutoff),
          fftw_(c_.grid, c_.viscosity, c_.cutoff) {
        continuation_detail::validateConfiguration(c_);
        if (c_.sampling_grid != c_.grid) {
            sampler_grid_.reset(new PseudospectralSystem(c_.sampling_grid, c_.viscosity));
            sampler_.reset(new FftwReferenceSystem(c_.sampling_grid, c_.viscosity));
        }
        if (c_.dense_sampling_grid != c_.sampling_grid) {
            dense_grid_.reset(new PseudospectralSystem(c_.dense_sampling_grid, c_.viscosity));
            dense_.reset(new FftwReferenceSystem(c_.dense_sampling_grid, c_.viscosity));
        }
    }

    const PseudospectralSystem& spectralSystem() const { return system_; }

    ContinuationObservation observe(const OptimizationState& state, double time) const {
        ContinuationObservation d;
        d.time = time;
        d.energy = system_.energy(state);
        d.enstrophy = system_.enstrophy(state);
        d.palinstrophy = system_.palinstrophy(state);
        d.h_half = system_.criticalHOneHalf(state);
        d.k_rms = std::sqrt(d.enstrophy / d.energy);
        d.cutoff_fraction = system_.cutoffShellEnergyFraction(state);
        d.divergence_defect = system_.divergenceDefect(state);
        d.reality_defect = system_.realityDefect(state);
        {
            const auto nonlinear = c_.fftw ? fftw_.rightHandSide(state, false)
                                           : system_.rightHandSide(state, false);
            d.stretching = system_.enstrophyDerivative(state, nonlinear);
            d.nonlinear_energy_residual = std::abs(system_.energyDerivative(state, nonlinear));
        }
        d.viscous_destruction = 2.0 * c_.viscosity * d.palinstrophy;
        d.net_enstrophy_rate = d.stretching - d.viscous_destruction;
        d.production_to_dissipation = d.stretching / d.viscous_destruction;
        const auto samples = sampler_
            ? sampler_->diagnostics(liftOptimizationState(system_, state, *sampler_grid_))
            : fftw_.diagnostics(state);
        d.l3_sample = samples.critical_l3_sample;
        d.vorticity_sample = samples.sampled_vorticity_max;
        d.l3_dense = d.l3_sample;
        d.vorticity_dense = d.vorticity_sample;
        if (dense_) {
            const auto dense_samples = dense_->diagnostics(liftOptimizationState(system_, state, *dense_grid_));
            d.l3_dense = dense_samples.critical_l3_sample;
            d.vorticity_dense = dense_samples.sampled_vorticity_max;
        }
        const auto& modes = system_.gridModes();
        for (std::size_t i = 0; i < state.size(); ++i) {
            d.vorticity_fourier_bound += norm(cross(modes[i], state[i]));
        }
        for (double value : continuation_detail::observationValues(d)) {
            if (!std::isfinite(value)) throw std::runtime_error("Nonfinite continuation diagnostic");
        }
        return d;
    }

    ContinuationCheckpoint initialize(const OptimizationState& state) const {
        requireOptimizationState(system_, state);
        ContinuationCheckpoint p;
        p.configuration = c_;
        p.state = state;
        p.initial = p.latest = observe(state, 0.0);
        p.peak_cutoff_fraction = p.latest.cutoff_fraction;
        p.cutoff_stopped = p.peak_cutoff_fraction > c_.cutoff_limit;
        continuation_detail::validateCheckpoint(p);
        return p;
    }

    // Endpoints must be observations on the global clock, so splitting a run
    // never adds a hidden small step or changes when a diagnostic is measured.
    void advance(ContinuationCheckpoint& p, double final_time,
                 const std::function<void(const ContinuationCheckpoint&)>& observer) const {
        continuation_detail::validateCheckpoint(p);
        if (!continuation_detail::sameConfiguration(c_, p.configuration)) {
            throw std::invalid_argument("Continuation cannot change checkpoint configuration");
        }
        const std::uint64_t final_index = continuation_detail::observationIndex(
            final_time, c_.observation_interval);
        if (p.cutoff_stopped || final_index <= p.completed_observations) {
            throw std::invalid_argument("Continuation needs a later endpoint and a passing checkpoint");
        }
        double time = p.latest.time;
        while (p.completed_observations < final_index) {
            const double target = (p.completed_observations + 1) * c_.observation_interval;
            double previous_energy = system_.energy(p.state);
            while (time < target) {
                const double cap = std::min(c_.maximum_time_step, target - time);
                double dt, cfl, viscous;
                if (c_.fftw) {
                    const auto step = fftw_.chooseAdaptiveTimeStep(p.state, cap, c_.target_cfl, c_.diffusion_safety);
                    dt = step.time_step; cfl = step.advective_cfl_upper_bound; viscous = step.viscous_stability_number;
                    fftw_.stepRungeKutta4Compact(p.state, dt);
                } else {
                    const auto step = system_.chooseAdaptiveTimeStep(p.state, cap, c_.target_cfl, c_.diffusion_safety);
                    dt = step.time_step; cfl = step.advective_cfl_upper_bound; viscous = step.viscous_stability_number;
                    system_.stepRungeKutta4(p.state, dt);
                }
                if (time + dt == time) throw std::runtime_error("Continuation time step underflow");
                time = dt == target - time ? target : time + dt;
                ++p.step;
                p.minimum_time_step = p.minimum_time_step == 0.0 ? dt : std::min(p.minimum_time_step, dt);
                p.maximum_time_step = std::max(p.maximum_time_step, dt);
                p.maximum_cfl = std::max(p.maximum_cfl, cfl);
                p.maximum_viscous_number = std::max(p.maximum_viscous_number, viscous);
                const double energy = system_.energy(p.state);
                if (!std::isfinite(energy) || energy <= 0.0 || energy > previous_energy * (1.0 + 1e-10)) {
                    throw std::runtime_error("Continuation violated unforced energy decay; last good checkpoint retained");
                }
                previous_energy = energy;
                const double cutoff = system_.cutoffShellEnergyFraction(p.state);
                if (!std::isfinite(cutoff)) throw std::runtime_error("Nonfinite cutoff monitor");
                p.peak_cutoff_fraction = std::max(p.peak_cutoff_fraction, cutoff);
                p.cutoff_stopped = p.peak_cutoff_fraction > c_.cutoff_limit;
                if (p.cutoff_stopped) break;
            }
            const auto observation = observe(p.state, time);
            if (observation.divergence_defect > 1e-9 || observation.reality_defect > 1e-9) {
                throw std::runtime_error("Continuation lost Fourier reality or incompressibility");
            }
            p.bkm_sampled_integral += 0.5 * (time - p.latest.time) *
                (p.latest.vorticity_sample + observation.vorticity_sample);
            p.latest = observation;
            if (time == target) ++p.completed_observations;
            observer(p);
            if (p.cutoff_stopped) return;
        }
    }

private:
    ContinuationConfiguration c_;
    PseudospectralSystem system_;
    FftwReferenceSystem fftw_;
    std::unique_ptr<PseudospectralSystem> sampler_grid_;
    std::unique_ptr<FftwReferenceSystem> sampler_;
    std::unique_ptr<PseudospectralSystem> dense_grid_;
    std::unique_ptr<FftwReferenceSystem> dense_;
};

// Same endian-stable doubles, checksum, exact-length checks and atomic rename
// discipline as the older checkpoint implementation; distinct version/magic.
inline std::vector<unsigned char> encodeContinuationCheckpoint(const ContinuationCheckpoint& p) {
    continuation_detail::validateCheckpoint(p);
    using namespace checkpoint_detail;
    std::vector<unsigned char> bytes;
    bytes.reserve(512 + p.state.size() * 48);
    const auto& c = p.configuration;
    for (int value : {c.grid, c.cutoff, c.sampling_grid, c.dense_sampling_grid}) appendUnsigned64(bytes, value);
    appendUnsigned64(bytes, c.fftw ? 1 : 0);
    for (double value : {c.viscosity, c.maximum_time_step, c.target_cfl, c.diffusion_safety,
                         c.observation_interval, c.cutoff_limit}) appendDouble(bytes, value);
    appendUnsigned64(bytes, p.step);
    appendUnsigned64(bytes, p.completed_observations);
    appendUnsigned64(bytes, p.cutoff_stopped ? 1 : 0);
    for (double value : {p.peak_cutoff_fraction, p.minimum_time_step, p.maximum_time_step,
                         p.maximum_cfl, p.maximum_viscous_number, p.bkm_sampled_integral}) appendDouble(bytes, value);
    continuation_detail::appendObservation(bytes, p.initial);
    continuation_detail::appendObservation(bytes, p.latest);
    appendUnsigned64(bytes, p.state.size());
    for (const auto& value : p.state) {
        for (double part : {value.x.real(), value.x.imag(), value.y.real(), value.y.imag(), value.z.real(), value.z.imag()}) {
            appendDouble(bytes, part);
        }
    }
    return bytes;
}

inline void saveContinuationCheckpoint(const std::string& path, const ContinuationCheckpoint& p) {
    if (path.empty()) throw std::invalid_argument("Checkpoint path is empty");
    const auto bytes = encodeContinuationCheckpoint(p);
    std::vector<unsigned char> header;
    checkpoint_detail::appendUnsigned64(header, bytes.size());
    checkpoint_detail::appendUnsigned64(header, checkpoint_detail::checksum(bytes));
    const std::string temporary = path + ".tmp";
    std::ofstream output(temporary.c_str(), std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot open checkpoint output: " + temporary);
    output.write("NSCONT1\n", 8);
    output.write(reinterpret_cast<const char*>(header.data()), header.size());
    output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    output.close();
    if (!output || std::rename(temporary.c_str(), path.c_str()) != 0) {
        std::remove(temporary.c_str());
        throw std::runtime_error("Cannot atomically save continuation checkpoint: " + path);
    }
}

inline ContinuationCheckpoint loadContinuationCheckpoint(const std::string& path) {
    using namespace checkpoint_detail;
    std::ifstream input(path.c_str(), std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("Cannot open continuation checkpoint: " + path);
    const std::streamoff file_size = input.tellg();
    if (file_size < 24 || file_size > INT64_C(16) * 1024 * 1024 * 1024) {
        throw std::runtime_error("Invalid continuation checkpoint length");
    }
    input.seekg(0);
    char magic[8];
    input.read(magic, 8);
    if (!input || std::memcmp(magic, "NSCONT1\n", 8) != 0) {
        throw std::runtime_error("Invalid continuation checkpoint magic/version");
    }
    std::vector<unsigned char> header(16);
    input.read(reinterpret_cast<char*>(header.data()), 16);
    std::size_t offset = 0;
    const auto length = readUnsigned64(header, offset);
    const auto expected_checksum = readUnsigned64(header, offset);
    if (length != static_cast<std::uint64_t>(file_size - 24)) {
        throw std::runtime_error("Truncated or trailing continuation checkpoint bytes");
    }
    std::vector<unsigned char> bytes(static_cast<std::size_t>(length));
    input.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
    if (!input || checksum(bytes) != expected_checksum) throw std::runtime_error("Continuation checksum mismatch");
    ContinuationCheckpoint p;
    auto& c = p.configuration;
    offset = 0;
    int* integers[] = {&c.grid, &c.cutoff, &c.sampling_grid, &c.dense_sampling_grid};
    for (int* value : integers) {
        const auto raw = readUnsigned64(bytes, offset);
        if (raw > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) throw std::runtime_error("Oversized grid metadata");
        *value = static_cast<int>(raw);
    }
    const auto backend = readUnsigned64(bytes, offset);
    if (backend > 1) throw std::runtime_error("Invalid continuation backend");
    c.fftw = backend == 1;
    double* config_values[] = {&c.viscosity, &c.maximum_time_step, &c.target_cfl,
        &c.diffusion_safety, &c.observation_interval, &c.cutoff_limit};
    for (double* value : config_values) *value = readDouble(bytes, offset);
    continuation_detail::validateConfiguration(c);
    p.step = readUnsigned64(bytes, offset);
    p.completed_observations = readUnsigned64(bytes, offset);
    const auto stopped = readUnsigned64(bytes, offset);
    if (stopped > 1) throw std::runtime_error("Invalid continuation stop flag");
    p.cutoff_stopped = stopped == 1;
    double* progress_values[] = {&p.peak_cutoff_fraction, &p.minimum_time_step,
        &p.maximum_time_step, &p.maximum_cfl, &p.maximum_viscous_number, &p.bkm_sampled_integral};
    for (double* value : progress_values) *value = readDouble(bytes, offset);
    p.initial = continuation_detail::readObservation(bytes, offset);
    p.latest = continuation_detail::readObservation(bytes, offset);
    const auto state_size = readUnsigned64(bytes, offset);
    if (state_size != expectedStateSize(c.grid) || state_size != (bytes.size() - offset) / 48 ||
        (bytes.size() - offset) % 48 != 0) throw std::runtime_error("Invalid continuation state length");
    p.state.resize(static_cast<std::size_t>(state_size));
    for (auto& value : p.state) {
        double parts[6];
        for (double& part : parts) part = readDouble(bytes, offset);
        value = ComplexVector(Complex(parts[0], parts[1]), Complex(parts[2], parts[3]), Complex(parts[4], parts[5]));
    }
    continuation_detail::validateCheckpoint(p);
    return p;
}

}  // namespace ns_cascade
#endif
