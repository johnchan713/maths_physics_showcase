#ifndef NS_CASCADE_SPECTRAL_PROFILE_HPP
#define NS_CASCADE_SPECTRAL_PROFILE_HPP

#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ns_cascade {

struct SpectrumProfile {
    double characteristic_wavenumber = 0.0;
    double rescaled_mean = 0.0;
    double rescaled_standard_deviation = 0.0;
    double maximum_rescaled_coordinate = 0.0;
    double analyticity_radius_estimate = 0.0;
    double analyticity_fit_r_squared = 0.0;
    int analyticity_fit_shells = 0;
    bool analyticity_fit_valid = false;
    std::vector<double> energy_fractions;
};

struct SpectrumProfileChange {
    double l1_distance = 0.0;
    double overlap = 1.0;
    double absolute_log_scale_change = 0.0;
    double l1_per_log_scale_change = 0.0;
    bool scale_normalized_drift_valid = false;
};

inline SpectrumProfile rescaledSpectrumProfile(
    const PseudospectralSystem& system,
    const PseudospectralSystem::State& state,
    int finite_bin_count = 32,
    double maximum_rescaled_coordinate = 4.0) {
    if (state.size() != system.gridPointCount()) {
        throw std::invalid_argument("Spectrum profile state has the wrong grid size");
    }
    if (finite_bin_count < 4 || finite_bin_count > 4096) {
        throw std::invalid_argument(
            "Spectrum profile requires between 4 and 4096 finite bins");
    }
    if (!std::isfinite(maximum_rescaled_coordinate) ||
        maximum_rescaled_coordinate <= 0.0) {
        throw std::invalid_argument(
            "Maximum rescaled spectral coordinate must be finite and positive");
    }

    SpectrumProfile profile;
    profile.maximum_rescaled_coordinate = maximum_rescaled_coordinate;
    profile.energy_fractions.assign(
        static_cast<std::size_t>(finite_bin_count + 1), 0.0);

    const double total_energy = system.energy(state);
    if (!std::isfinite(total_energy) || total_energy <= 0.0) {
        throw std::invalid_argument(
            "Spectrum profile requires a finite state with positive energy");
    }
    profile.characteristic_wavenumber =
        std::sqrt(system.enstrophy(state) / total_energy);
    if (!std::isfinite(profile.characteristic_wavenumber) ||
        profile.characteristic_wavenumber <= 0.0) {
        throw std::runtime_error("Characteristic wavenumber is not positive");
    }

    const int radial_shell_count = static_cast<int>(std::ceil(
        std::sqrt(3.0) * static_cast<double>(system.cutoff()) - 1e-12));
    std::vector<double> radial_shell_energy(
        static_cast<std::size_t>(radial_shell_count + 1), 0.0);
    double rescaled_first_moment = 0.0;
    double rescaled_second_moment = 0.0;
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t index = 0; index < state.size(); ++index) {
        const int wave_squared = modes[index].normSquared();
        if (wave_squared == 0) continue;
        const double mode_energy = 0.5 * normSquared(state[index]);
        if (mode_energy == 0.0) continue;
        const double wave_magnitude =
            std::sqrt(static_cast<double>(wave_squared));
        const double rescaled_wave =
            wave_magnitude / profile.characteristic_wavenumber;
        const double fraction = mode_energy / total_energy;
        rescaled_first_moment += rescaled_wave * fraction;
        rescaled_second_moment += rescaled_wave * rescaled_wave * fraction;

        // Deposit onto neighbouring bin centres rather than using a hard
        // histogram edge.  This cloud-in-cell projection makes the profile
        // distance continuous when k_rms moves by a tiny amount.
        const double bin_width =
            maximum_rescaled_coordinate / finite_bin_count;
        const double bin_position = rescaled_wave / bin_width - 0.5;
        if (bin_position <= 0.0) {
            profile.energy_fractions[0] += fraction;
        } else if (bin_position >= finite_bin_count) {
            profile.energy_fractions[
                static_cast<std::size_t>(finite_bin_count)] += fraction;
        } else {
            const int lower_bin = static_cast<int>(std::floor(bin_position));
            const double upper_weight = bin_position - lower_bin;
            profile.energy_fractions[static_cast<std::size_t>(lower_bin)] +=
                fraction * (1.0 - upper_weight);
            profile.energy_fractions[
                static_cast<std::size_t>(lower_bin + 1)] +=
                fraction * upper_weight;
        }

        const int shell = std::max(
            1, static_cast<int>(std::ceil(wave_magnitude - 1e-12)));
        radial_shell_energy[static_cast<std::size_t>(shell)] += mode_energy;
    }
    profile.rescaled_mean = rescaled_first_moment;
    profile.rescaled_standard_deviation = std::sqrt(std::max(
        0.0,
        rescaled_second_moment -
            rescaled_first_moment * rescaled_first_moment));

    // For an analytic velocity field, a radial energy tail often has the form
    // E(k) approximately k^a exp(-2 delta k).  A straight-line fit of log E(k)
    // over the upper half of the isotropic (non-corner) range gives the rough
    // diagnostic delta=-slope/2.  The fit is deliberately reported with its
    // point count and R^2 because it is not a rigorous analyticity bound.
    const int first_fit_shell = std::max(2, (system.cutoff() + 1) / 2);
    const int last_fit_shell = system.cutoff();
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xx = 0.0;
    double sum_xy = 0.0;
    std::vector<double> fit_x;
    std::vector<double> fit_y;
    for (int shell = first_fit_shell; shell <= last_fit_shell; ++shell) {
        const double shell_energy =
            radial_shell_energy[static_cast<std::size_t>(shell)];
        if (shell_energy <= total_energy * 1e-30) continue;
        const double x = static_cast<double>(shell) - 0.5;
        const double y = std::log(shell_energy / total_energy);
        fit_x.push_back(x);
        fit_y.push_back(y);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }
    profile.analyticity_fit_shells = static_cast<int>(fit_x.size());
    if (fit_x.size() >= 4) {
        const double count = static_cast<double>(fit_x.size());
        const double denominator = count * sum_xx - sum_x * sum_x;
        if (denominator > 0.0) {
            const double slope =
                (count * sum_xy - sum_x * sum_y) / denominator;
            const double intercept = (sum_y - slope * sum_x) / count;
            double residual_sum = 0.0;
            double total_sum = 0.0;
            const double mean_y = sum_y / count;
            for (std::size_t i = 0; i < fit_x.size(); ++i) {
                const double residual =
                    fit_y[i] - (intercept + slope * fit_x[i]);
                residual_sum += residual * residual;
                const double centered = fit_y[i] - mean_y;
                total_sum += centered * centered;
            }
            profile.analyticity_fit_r_squared =
                total_sum == 0.0 ? 1.0 : 1.0 - residual_sum / total_sum;
            if (slope < 0.0 &&
                std::isfinite(profile.analyticity_fit_r_squared)) {
                profile.analyticity_radius_estimate = -0.5 * slope;
                profile.analyticity_fit_valid = true;
            }
        }
    }
    return profile;
}

inline double spectrumProfileL1Distance(const SpectrumProfile& left,
                                        const SpectrumProfile& right) {
    if (left.energy_fractions.size() != right.energy_fractions.size() ||
        left.energy_fractions.empty() ||
        left.maximum_rescaled_coordinate !=
            right.maximum_rescaled_coordinate) {
        throw std::invalid_argument(
            "Spectrum profiles use incompatible rescaled bins");
    }
    double distance = 0.0;
    for (std::size_t i = 0; i < left.energy_fractions.size(); ++i) {
        distance += std::abs(left.energy_fractions[i] -
                             right.energy_fractions[i]);
    }
    return distance;
}

inline double spectrumProfileOverlap(const SpectrumProfile& left,
                                     const SpectrumProfile& right) {
    return std::max(0.0, 1.0 - 0.5 * spectrumProfileL1Distance(left, right));
}

inline SpectrumProfileChange spectrumProfileChange(
    const SpectrumProfile& current,
    const SpectrumProfile& reference) {
    if (!std::isfinite(current.characteristic_wavenumber) ||
        !std::isfinite(reference.characteristic_wavenumber) ||
        current.characteristic_wavenumber <= 0.0 ||
        reference.characteristic_wavenumber <= 0.0) {
        throw std::invalid_argument(
            "Spectrum-profile scale comparison requires positive finite "
            "characteristic wavenumbers");
    }
    SpectrumProfileChange change;
    change.l1_distance = spectrumProfileL1Distance(current, reference);
    change.overlap = std::max(0.0, 1.0 - 0.5 * change.l1_distance);
    change.absolute_log_scale_change = std::abs(std::log(
        current.characteristic_wavenumber /
        reference.characteristic_wavenumber));
    const double roundoff_floor =
        64.0 * std::numeric_limits<double>::epsilon();
    if (change.absolute_log_scale_change > roundoff_floor) {
        change.l1_per_log_scale_change =
            change.l1_distance / change.absolute_log_scale_change;
        change.scale_normalized_drift_valid =
            std::isfinite(change.l1_per_log_scale_change);
    }
    return change;
}

}  // namespace ns_cascade

#endif
