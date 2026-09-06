#ifndef NS_CASCADE_CHECKPOINT_HPP
#define NS_CASCADE_CHECKPOINT_HPP

#include "ns_cascade/spectral_profile.hpp"

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ns_cascade {

struct CheckpointConfiguration {
    int grid_size = 0;
    int cutoff = 0;
    double viscosity = 0.0;
    InitialCondition initial_condition = InitialCondition::TaylorGreen;
    VortexTubeParameters tube_parameters;
    double initial_energy = 0.0;
    bool adaptive = true;
    double maximum_time_step = 0.0;
    double target_cfl = 0.0;
    double diffusion_safety = 0.0;
    int profile_bin_count = 0;
    double profile_maximum_coordinate = 0.0;
};

struct CheckpointProgress {
    std::uint64_t step = 0;
    double time = 0.0;
    double bkm_sampled_integral = 0.0;
    double previous_vorticity_max = 0.0;
    double previous_diagnostic_time = 0.0;
    double initial_critical_l3 = 0.0;
    double initial_critical_h_half = 0.0;
    double initial_sampled_vorticity = 0.0;
    double peak_high_shell_fraction = 0.0;
    double peak_critical_l3 = 0.0;
    double peak_critical_h_half = 0.0;
    double peak_sampled_vorticity = 0.0;
    double maximum_production_to_dissipation = 0.0;
    double peak_forward_flux = 0.0;
    double minimum_time_step = 0.0;
    double maximum_time_step = 0.0;
    double maximum_cfl_bound = 0.0;
    double maximum_viscous_number = 0.0;
};

struct SimulationCheckpoint {
    CheckpointConfiguration configuration;
    CheckpointProgress progress;
    SpectrumProfile initial_profile;
    SpectrumProfile previous_profile;
    PseudospectralSystem::State state;
};

namespace checkpoint_detail {

inline void appendUnsigned64(std::vector<unsigned char>& bytes,
                             std::uint64_t value) {
    for (int shift = 0; shift < 64; shift += 8) {
        bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xffU));
    }
}

inline void appendDouble(std::vector<unsigned char>& bytes, double value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t),
                  "Checkpoint format requires 64-bit IEEE-style doubles");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(value));
    appendUnsigned64(bytes, bits);
}

inline std::uint64_t readUnsigned64(const std::vector<unsigned char>& bytes,
                                    std::size_t& offset) {
    if (offset > bytes.size() || bytes.size() - offset < 8) {
        throw std::runtime_error("Checkpoint payload ended unexpectedly");
    }
    std::uint64_t value = 0;
    for (int shift = 0; shift < 64; shift += 8) {
        value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
    }
    return value;
}

inline double readDouble(const std::vector<unsigned char>& bytes,
                         std::size_t& offset) {
    const std::uint64_t bits = readUnsigned64(bytes, offset);
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

inline std::uint64_t checksum(const std::vector<unsigned char>& bytes) {
    std::uint64_t hash = UINT64_C(1469598103934665603);
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        hash ^= static_cast<std::uint64_t>(bytes[i]);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

inline void appendProfile(std::vector<unsigned char>& bytes,
                          const SpectrumProfile& profile) {
    appendDouble(bytes, profile.characteristic_wavenumber);
    appendDouble(bytes, profile.rescaled_mean);
    appendDouble(bytes, profile.rescaled_standard_deviation);
    appendDouble(bytes, profile.maximum_rescaled_coordinate);
    appendDouble(bytes, profile.analyticity_radius_estimate);
    appendDouble(bytes, profile.analyticity_fit_r_squared);
    appendUnsigned64(bytes,
                     static_cast<std::uint64_t>(profile.analyticity_fit_shells));
    appendUnsigned64(bytes, profile.analyticity_fit_valid ? 1U : 0U);
    appendUnsigned64(bytes,
                     static_cast<std::uint64_t>(profile.energy_fractions.size()));
    for (std::size_t i = 0; i < profile.energy_fractions.size(); ++i) {
        appendDouble(bytes, profile.energy_fractions[i]);
    }
}

inline SpectrumProfile readProfile(const std::vector<unsigned char>& bytes,
                                   std::size_t& offset) {
    SpectrumProfile profile;
    profile.characteristic_wavenumber = readDouble(bytes, offset);
    profile.rescaled_mean = readDouble(bytes, offset);
    profile.rescaled_standard_deviation = readDouble(bytes, offset);
    profile.maximum_rescaled_coordinate = readDouble(bytes, offset);
    profile.analyticity_radius_estimate = readDouble(bytes, offset);
    profile.analyticity_fit_r_squared = readDouble(bytes, offset);
    const std::uint64_t fit_shells = readUnsigned64(bytes, offset);
    const std::uint64_t fit_valid = readUnsigned64(bytes, offset);
    const std::uint64_t bin_count = readUnsigned64(bytes, offset);
    if (fit_shells > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        fit_valid > 1U || bin_count < 5U || bin_count > 4097U) {
        throw std::runtime_error("Checkpoint contains invalid profile metadata");
    }
    profile.analyticity_fit_shells = static_cast<int>(fit_shells);
    profile.analyticity_fit_valid = fit_valid == 1U;
    profile.energy_fractions.resize(static_cast<std::size_t>(bin_count));
    for (std::size_t i = 0; i < profile.energy_fractions.size(); ++i) {
        profile.energy_fractions[i] = readDouble(bytes, offset);
    }
    return profile;
}

inline void requireFinite(double value, const char* label) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string("Checkpoint has non-finite ") + label);
    }
}

inline void validateProfile(const SpectrumProfile& profile,
                            int configured_finite_bins,
                            double configured_maximum_coordinate) {
    if (profile.energy_fractions.size() !=
        static_cast<std::size_t>(configured_finite_bins + 1) ||
        profile.maximum_rescaled_coordinate != configured_maximum_coordinate) {
        throw std::runtime_error("Checkpoint spectrum profile shape is inconsistent");
    }
    requireFinite(profile.characteristic_wavenumber,
                  "profile characteristic wavenumber");
    requireFinite(profile.rescaled_mean, "profile mean");
    requireFinite(profile.rescaled_standard_deviation,
                  "profile standard deviation");
    requireFinite(profile.analyticity_radius_estimate,
                  "analyticity-radius estimate");
    requireFinite(profile.analyticity_fit_r_squared, "analyticity-fit R squared");
    double sum = 0.0;
    for (std::size_t i = 0; i < profile.energy_fractions.size(); ++i) {
        const double value = profile.energy_fractions[i];
        requireFinite(value, "profile energy fraction");
        if (value < 0.0) {
            throw std::runtime_error(
                "Checkpoint profile has a negative energy fraction");
        }
        sum += value;
    }
    if (std::abs(sum - 1.0) > 1e-10) {
        throw std::runtime_error(
            "Checkpoint profile energy fractions do not sum to one");
    }
}

inline std::size_t expectedStateSize(int grid_size) {
    if (grid_size < 8 || (grid_size & (grid_size - 1)) != 0) {
        throw std::runtime_error("Checkpoint FFT grid size is invalid");
    }
    const std::size_t grid = static_cast<std::size_t>(grid_size);
    if (grid > std::numeric_limits<std::size_t>::max() / grid ||
        grid * grid > std::numeric_limits<std::size_t>::max() / grid) {
        throw std::runtime_error("Checkpoint FFT grid size overflows memory size");
    }
    return grid * grid * grid;
}

inline void validateCheckpoint(const SimulationCheckpoint& checkpoint) {
    const CheckpointConfiguration& configuration = checkpoint.configuration;
    const CheckpointProgress& progress = checkpoint.progress;
    const double pi = 3.1415926535897932384626433832795;
    const double two_pi = 2.0 * pi;
    const int maximum_cutoff = (configuration.grid_size - 1) / 3;
    if (configuration.cutoff < 1 || configuration.cutoff > maximum_cutoff) {
        throw std::runtime_error("Checkpoint Fourier cutoff is invalid");
    }
    requireFinite(configuration.viscosity, "viscosity");
    requireFinite(configuration.initial_energy, "initial energy");
    requireFinite(configuration.maximum_time_step, "maximum time step");
    requireFinite(configuration.target_cfl, "CFL target");
    requireFinite(configuration.diffusion_safety, "diffusion safety");
    requireFinite(configuration.profile_maximum_coordinate,
                  "profile maximum coordinate");
    requireFinite(configuration.tube_parameters.core_radius,
                  "vortex-tube core radius");
    requireFinite(configuration.tube_parameters.separation,
                  "vortex-tube separation");
    requireFinite(configuration.tube_parameters.bend_amplitude,
                  "vortex-tube bend amplitude");
    if (configuration.viscosity < 0.0 || configuration.initial_energy <= 0.0 ||
        configuration.maximum_time_step <= 0.0 || configuration.target_cfl <= 0.0 ||
        configuration.diffusion_safety <= 0.0 ||
        configuration.profile_bin_count < 4 ||
        configuration.profile_bin_count > 4096 ||
        configuration.profile_maximum_coordinate <= 0.0) {
        throw std::runtime_error("Checkpoint simulation configuration is invalid");
    }
    if (configuration.tube_parameters.core_radius <= 0.0 ||
        configuration.tube_parameters.core_radius > pi ||
        configuration.tube_parameters.separation <= 0.0 ||
        configuration.tube_parameters.separation >= two_pi ||
        configuration.tube_parameters.bend_amplitude < 0.0 ||
        configuration.tube_parameters.bend_amplitude >= pi ||
        configuration.tube_parameters.axial_wavenumber < 1 ||
        configuration.tube_parameters.axial_wavenumber >
            configuration.cutoff) {
        throw std::runtime_error(
            "Checkpoint vortex-tube parameters are invalid");
    }
    const int condition = static_cast<int>(configuration.initial_condition);
    if (condition < static_cast<int>(InitialCondition::Deterministic) ||
        condition > static_cast<int>(InitialCondition::VortexTubes)) {
        throw std::runtime_error("Checkpoint initial-condition identifier is invalid");
    }
    if (checkpoint.state.size() != expectedStateSize(configuration.grid_size)) {
        throw std::runtime_error("Checkpoint state length does not match its grid");
    }
    const double* progress_values[] = {
        &progress.time,
        &progress.bkm_sampled_integral,
        &progress.previous_vorticity_max,
        &progress.previous_diagnostic_time,
        &progress.initial_critical_l3,
        &progress.initial_critical_h_half,
        &progress.initial_sampled_vorticity,
        &progress.peak_high_shell_fraction,
        &progress.peak_critical_l3,
        &progress.peak_critical_h_half,
        &progress.peak_sampled_vorticity,
        &progress.maximum_production_to_dissipation,
        &progress.peak_forward_flux,
        &progress.minimum_time_step,
        &progress.maximum_time_step,
        &progress.maximum_cfl_bound,
        &progress.maximum_viscous_number};
    for (std::size_t i = 0;
         i < sizeof(progress_values) / sizeof(progress_values[0]);
         ++i) {
        requireFinite(*progress_values[i], "progress value");
    }
    if (progress.time < 0.0 || progress.previous_diagnostic_time < 0.0 ||
        progress.previous_diagnostic_time > progress.time ||
        progress.bkm_sampled_integral < 0.0) {
        throw std::runtime_error("Checkpoint time progress is inconsistent");
    }
    validateProfile(checkpoint.initial_profile,
                    configuration.profile_bin_count,
                    configuration.profile_maximum_coordinate);
    validateProfile(checkpoint.previous_profile,
                    configuration.profile_bin_count,
                    configuration.profile_maximum_coordinate);
    for (std::size_t i = 0; i < checkpoint.state.size(); ++i) {
        const ComplexVector& value = checkpoint.state[i];
        const double parts[] = {value.x.real(), value.x.imag(),
                                value.y.real(), value.y.imag(),
                                value.z.real(), value.z.imag()};
        for (std::size_t part = 0; part < 6; ++part) {
            requireFinite(parts[part], "Fourier coefficient");
        }
    }
}

inline std::vector<unsigned char> encodePayload(
    const SimulationCheckpoint& checkpoint) {
    validateCheckpoint(checkpoint);
    std::vector<unsigned char> bytes;
    bytes.reserve(512 + checkpoint.state.size() * 6 * sizeof(double));
    const CheckpointConfiguration& c = checkpoint.configuration;
    const CheckpointProgress& p = checkpoint.progress;
    appendUnsigned64(bytes, static_cast<std::uint64_t>(c.grid_size));
    appendUnsigned64(bytes, static_cast<std::uint64_t>(c.cutoff));
    appendDouble(bytes, c.viscosity);
    appendUnsigned64(bytes, static_cast<std::uint64_t>(c.initial_condition));
    appendDouble(bytes, c.tube_parameters.core_radius);
    appendDouble(bytes, c.tube_parameters.separation);
    appendDouble(bytes, c.tube_parameters.bend_amplitude);
    appendUnsigned64(bytes,
                     static_cast<std::uint64_t>(c.tube_parameters.axial_wavenumber));
    appendDouble(bytes, c.initial_energy);
    appendUnsigned64(bytes, c.adaptive ? 1U : 0U);
    appendDouble(bytes, c.maximum_time_step);
    appendDouble(bytes, c.target_cfl);
    appendDouble(bytes, c.diffusion_safety);
    appendUnsigned64(bytes, static_cast<std::uint64_t>(c.profile_bin_count));
    appendDouble(bytes, c.profile_maximum_coordinate);

    appendUnsigned64(bytes, p.step);
    appendDouble(bytes, p.time);
    appendDouble(bytes, p.bkm_sampled_integral);
    appendDouble(bytes, p.previous_vorticity_max);
    appendDouble(bytes, p.previous_diagnostic_time);
    appendDouble(bytes, p.initial_critical_l3);
    appendDouble(bytes, p.initial_critical_h_half);
    appendDouble(bytes, p.initial_sampled_vorticity);
    appendDouble(bytes, p.peak_high_shell_fraction);
    appendDouble(bytes, p.peak_critical_l3);
    appendDouble(bytes, p.peak_critical_h_half);
    appendDouble(bytes, p.peak_sampled_vorticity);
    appendDouble(bytes, p.maximum_production_to_dissipation);
    appendDouble(bytes, p.peak_forward_flux);
    appendDouble(bytes, p.minimum_time_step);
    appendDouble(bytes, p.maximum_time_step);
    appendDouble(bytes, p.maximum_cfl_bound);
    appendDouble(bytes, p.maximum_viscous_number);
    appendProfile(bytes, checkpoint.initial_profile);
    appendProfile(bytes, checkpoint.previous_profile);

    appendUnsigned64(bytes,
                     static_cast<std::uint64_t>(checkpoint.state.size()));
    for (std::size_t i = 0; i < checkpoint.state.size(); ++i) {
        const ComplexVector& value = checkpoint.state[i];
        appendDouble(bytes, value.x.real());
        appendDouble(bytes, value.x.imag());
        appendDouble(bytes, value.y.real());
        appendDouble(bytes, value.y.imag());
        appendDouble(bytes, value.z.real());
        appendDouble(bytes, value.z.imag());
    }
    return bytes;
}

inline SimulationCheckpoint decodePayload(
    const std::vector<unsigned char>& bytes) {
    SimulationCheckpoint checkpoint;
    std::size_t offset = 0;
    CheckpointConfiguration& c = checkpoint.configuration;
    CheckpointProgress& p = checkpoint.progress;
    const std::uint64_t grid_size = readUnsigned64(bytes, offset);
    const std::uint64_t cutoff = readUnsigned64(bytes, offset);
    if (grid_size > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        cutoff > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Checkpoint grid metadata is too large");
    }
    c.grid_size = static_cast<int>(grid_size);
    c.cutoff = static_cast<int>(cutoff);
    c.viscosity = readDouble(bytes, offset);
    const std::uint64_t condition = readUnsigned64(bytes, offset);
    if (condition > static_cast<std::uint64_t>(InitialCondition::VortexTubes)) {
        throw std::runtime_error("Checkpoint initial condition is invalid");
    }
    c.initial_condition = static_cast<InitialCondition>(condition);
    c.tube_parameters.core_radius = readDouble(bytes, offset);
    c.tube_parameters.separation = readDouble(bytes, offset);
    c.tube_parameters.bend_amplitude = readDouble(bytes, offset);
    const std::uint64_t axial_wavenumber = readUnsigned64(bytes, offset);
    if (axial_wavenumber >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Checkpoint axial wavenumber is too large");
    }
    c.tube_parameters.axial_wavenumber = static_cast<int>(axial_wavenumber);
    c.initial_energy = readDouble(bytes, offset);
    const std::uint64_t adaptive = readUnsigned64(bytes, offset);
    if (adaptive > 1U) throw std::runtime_error("Checkpoint adaptive flag is invalid");
    c.adaptive = adaptive == 1U;
    c.maximum_time_step = readDouble(bytes, offset);
    c.target_cfl = readDouble(bytes, offset);
    c.diffusion_safety = readDouble(bytes, offset);
    const std::uint64_t profile_bin_count = readUnsigned64(bytes, offset);
    if (profile_bin_count >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Checkpoint profile bin count is too large");
    }
    c.profile_bin_count = static_cast<int>(profile_bin_count);
    c.profile_maximum_coordinate = readDouble(bytes, offset);

    p.step = readUnsigned64(bytes, offset);
    p.time = readDouble(bytes, offset);
    p.bkm_sampled_integral = readDouble(bytes, offset);
    p.previous_vorticity_max = readDouble(bytes, offset);
    p.previous_diagnostic_time = readDouble(bytes, offset);
    p.initial_critical_l3 = readDouble(bytes, offset);
    p.initial_critical_h_half = readDouble(bytes, offset);
    p.initial_sampled_vorticity = readDouble(bytes, offset);
    p.peak_high_shell_fraction = readDouble(bytes, offset);
    p.peak_critical_l3 = readDouble(bytes, offset);
    p.peak_critical_h_half = readDouble(bytes, offset);
    p.peak_sampled_vorticity = readDouble(bytes, offset);
    p.maximum_production_to_dissipation = readDouble(bytes, offset);
    p.peak_forward_flux = readDouble(bytes, offset);
    p.minimum_time_step = readDouble(bytes, offset);
    p.maximum_time_step = readDouble(bytes, offset);
    p.maximum_cfl_bound = readDouble(bytes, offset);
    p.maximum_viscous_number = readDouble(bytes, offset);
    checkpoint.initial_profile = readProfile(bytes, offset);
    checkpoint.previous_profile = readProfile(bytes, offset);

    const std::uint64_t state_size = readUnsigned64(bytes, offset);
    const std::size_t expected_size = expectedStateSize(c.grid_size);
    if (state_size != static_cast<std::uint64_t>(expected_size)) {
        throw std::runtime_error("Checkpoint state size is inconsistent");
    }
    checkpoint.state.resize(expected_size);
    for (std::size_t i = 0; i < checkpoint.state.size(); ++i) {
        const double x_real = readDouble(bytes, offset);
        const double x_imaginary = readDouble(bytes, offset);
        const double y_real = readDouble(bytes, offset);
        const double y_imaginary = readDouble(bytes, offset);
        const double z_real = readDouble(bytes, offset);
        const double z_imaginary = readDouble(bytes, offset);
        const Complex x(x_real, x_imaginary);
        const Complex y(y_real, y_imaginary);
        const Complex z(z_real, z_imaginary);
        checkpoint.state[i] = ComplexVector(x, y, z);
    }
    if (offset != bytes.size()) {
        throw std::runtime_error("Checkpoint contains unexpected trailing payload");
    }
    validateCheckpoint(checkpoint);
    return checkpoint;
}

}  // namespace checkpoint_detail

inline void saveSimulationCheckpoint(const std::string& path,
                                     const SimulationCheckpoint& checkpoint) {
    if (path.empty()) throw std::invalid_argument("Checkpoint path cannot be empty");
    const std::vector<unsigned char> payload =
        checkpoint_detail::encodePayload(checkpoint);
    const std::string temporary_path = path + ".tmp";
    std::ofstream output(temporary_path.c_str(),
                         std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Cannot open temporary checkpoint: " +
                                 temporary_path);
    }
    const char magic[8] = {'N', 'S', 'C', 'C', 'H', 'K', '2', '\n'};
    output.write(magic, sizeof(magic));
    std::vector<unsigned char> header;
    checkpoint_detail::appendUnsigned64(
        header, static_cast<std::uint64_t>(payload.size()));
    checkpoint_detail::appendUnsigned64(
        header, checkpoint_detail::checksum(payload));
    output.write(reinterpret_cast<const char*>(&header[0]),
                 static_cast<std::streamsize>(header.size()));
    if (!payload.empty()) {
        output.write(reinterpret_cast<const char*>(&payload[0]),
                     static_cast<std::streamsize>(payload.size()));
    }
    output.close();
    if (!output) {
        std::remove(temporary_path.c_str());
        throw std::runtime_error("Failed while writing checkpoint: " + path);
    }
    if (std::rename(temporary_path.c_str(), path.c_str()) != 0) {
        const int saved_errno = errno;
        std::remove(temporary_path.c_str());
        throw std::runtime_error(
            "Cannot atomically publish checkpoint " + path + ": errno " +
            std::to_string(saved_errno));
    }
}

inline SimulationCheckpoint loadSimulationCheckpoint(const std::string& path) {
    if (path.empty()) throw std::invalid_argument("Checkpoint path cannot be empty");
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open checkpoint: " + path);
    char magic[8];
    input.read(magic, sizeof(magic));
    const char expected_magic[8] = {'N', 'S', 'C', 'C', 'H', 'K', '2', '\n'};
    if (!input || std::memcmp(magic, expected_magic, sizeof(magic)) != 0) {
        throw std::runtime_error("Checkpoint magic/version is invalid");
    }
    unsigned char raw_header[16];
    input.read(reinterpret_cast<char*>(raw_header), sizeof(raw_header));
    if (!input) throw std::runtime_error("Checkpoint header is truncated");
    const std::vector<unsigned char> header(raw_header, raw_header + 16);
    std::size_t header_offset = 0;
    const std::uint64_t payload_size =
        checkpoint_detail::readUnsigned64(header, header_offset);
    const std::uint64_t expected_checksum =
        checkpoint_detail::readUnsigned64(header, header_offset);
    const std::uint64_t maximum_payload_size = UINT64_C(16) * 1024 * 1024 * 1024;
    if (payload_size > maximum_payload_size ||
        payload_size > static_cast<std::uint64_t>(
                           std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("Checkpoint payload size is unreasonable");
    }
    std::vector<unsigned char> payload(static_cast<std::size_t>(payload_size));
    if (!payload.empty()) {
        input.read(reinterpret_cast<char*>(&payload[0]),
                   static_cast<std::streamsize>(payload.size()));
    }
    if (!input) throw std::runtime_error("Checkpoint payload is truncated");
    char trailing = 0;
    if (input.read(&trailing, 1)) {
        throw std::runtime_error("Checkpoint file has trailing bytes");
    }
    if (checkpoint_detail::checksum(payload) != expected_checksum) {
        throw std::runtime_error("Checkpoint checksum mismatch");
    }
    return checkpoint_detail::decodePayload(payload);
}

}  // namespace ns_cascade

#endif
