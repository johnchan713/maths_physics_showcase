#ifndef NS_CASCADE_CANDIDATE_SCORE_HPP
#define NS_CASCADE_CANDIDATE_SCORE_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ns_cascade {

// All inputs are dimensionless.  Keeping the score calculation separate from
// the solver makes the search policy deterministic, auditable, and directly
// testable without running an expensive three-dimensional trajectory.
struct SearchScoreEvidence {
    double peak_l3_ratio = 1.0;
    double peak_h_half_ratio = 1.0;
    double final_l3_ratio = 1.0;
    double peak_vorticity_ratio = 1.0;
    double peak_cutoff_fraction = 0.0;
    double cutoff_fraction_threshold = 0.01;
    int profile_scale_windows = 0;
    double first_profile_drift = 0.0;
    double latest_profile_drift = 0.0;
    double minimum_profile_drift = 0.0;
    double profile_drift_threshold = 1.0;
};

struct SearchScore {
    double critical_growth_reward = 0.0;
    double supporting_critical_reward = 0.0;
    double vorticity_tiebreak = 0.0;
    double final_l3_tiebreak = 0.0;
    double profile_drift_cost = 0.0;
    double profile_trend_reward = 0.0;
    double profile_rebound_cost = 0.0;
    double cutoff_cost = 0.0;
    double total = 0.0;
    bool profile_stationarity_improving = false;
};

inline bool finitePositiveScoreValue(double value) {
    return std::isfinite(value) && value > 0.0;
}

inline double clampScoreValue(double value, double lower, double upper) {
    return std::max(lower, std::min(value, upper));
}

inline double safeNonnegativeScoreRatio(double numerator,
                                        double positive_denominator) {
    if (numerator == 0.0) return 0.0;
    const double quotient = numerator / positive_denominator;
    return std::isfinite(quotient)
               ? quotient
               : std::numeric_limits<double>::max();
}

inline double safeRootMeanSquare(double left, double right) {
    const double scale = std::max(left, right);
    if (scale == 0.0) return 0.0;
    const double normalized_left = left / scale;
    const double normalized_right = right / scale;
    return scale * std::sqrt(
                       0.5 * (normalized_left * normalized_left +
                              normalized_right * normalized_right));
}

inline void validateSearchScoreEvidence(const SearchScoreEvidence& evidence) {
    if (!finitePositiveScoreValue(evidence.peak_l3_ratio) ||
        !finitePositiveScoreValue(evidence.peak_h_half_ratio) ||
        !finitePositiveScoreValue(evidence.final_l3_ratio) ||
        !finitePositiveScoreValue(evidence.peak_vorticity_ratio)) {
        throw std::invalid_argument(
            "Candidate score requires positive finite diagnostic ratios");
    }
    if (!std::isfinite(evidence.peak_cutoff_fraction) ||
        evidence.peak_cutoff_fraction < 0.0 ||
        !finitePositiveScoreValue(evidence.cutoff_fraction_threshold) ||
        !finitePositiveScoreValue(evidence.profile_drift_threshold) ||
        evidence.profile_scale_windows < 0) {
        throw std::invalid_argument("Candidate score has invalid gate metadata");
    }
    if (evidence.profile_scale_windows > 0 &&
        (!std::isfinite(evidence.latest_profile_drift) ||
         evidence.latest_profile_drift < 0.0 ||
         !std::isfinite(evidence.minimum_profile_drift) ||
         evidence.minimum_profile_drift < 0.0 ||
         evidence.minimum_profile_drift > evidence.latest_profile_drift)) {
        throw std::invalid_argument("Profile drift range is invalid");
    }
    if (evidence.profile_scale_windows > 1 &&
        (!std::isfinite(evidence.first_profile_drift) ||
         evidence.first_profile_drift < 0.0)) {
        throw std::invalid_argument("First profile drift is invalid");
    }
}

inline SearchScore scoreSingleResolution(
    const SearchScoreEvidence& evidence) {
    validateSearchScoreEvidence(evidence);
    SearchScore score;

    const double l3_log = std::log(evidence.peak_l3_ratio);
    const double h_half_log = std::log(evidence.peak_h_half_ratio);
    score.critical_growth_reward = std::max(l3_log, h_half_log);
    score.supporting_critical_reward = 0.25 * std::min(l3_log, h_half_log);
    score.vorticity_tiebreak =
        0.02 * std::log(evidence.peak_vorticity_ratio);
    score.final_l3_tiebreak = 0.01 * std::log(evidence.final_l3_ratio);

    // A candidate with no completed forward scale window has not supplied
    // stationarity evidence.  Give it a finite missing-evidence cost rather
    // than manufacturing a zero drift or rejecting short smoke runs.
    const double normalized_drift = evidence.profile_scale_windows > 0
                                        ? safeNonnegativeScoreRatio(
                                              evidence.latest_profile_drift,
                                              evidence.profile_drift_threshold)
                                        : 4.0;
    score.profile_drift_cost = 0.10 * std::log1p(normalized_drift);

    if (evidence.profile_scale_windows > 1) {
        const double drift_floor = 1e-15;
        const double trend = std::log(
            std::max(evidence.first_profile_drift, drift_floor) /
            std::max(evidence.latest_profile_drift, drift_floor));
        score.profile_trend_reward =
            0.05 * clampScoreValue(trend, -2.0, 2.0);
        const double rebound = std::log(
            std::max(evidence.latest_profile_drift, drift_floor) /
            std::max(evidence.minimum_profile_drift, drift_floor));
        score.profile_rebound_cost =
            0.05 * clampScoreValue(rebound, 0.0, 2.0);
        score.profile_stationarity_improving =
            evidence.latest_profile_drift <= evidence.first_profile_drift &&
            evidence.latest_profile_drift <=
                evidence.profile_drift_threshold &&
            evidence.latest_profile_drift <=
                1.10 * evidence.minimum_profile_drift;
    }

    const double normalized_cutoff =
        safeNonnegativeScoreRatio(evidence.peak_cutoff_fraction,
                                  evidence.cutoff_fraction_threshold);
    score.cutoff_cost = 0.03 * std::log1p(normalized_cutoff);
    score.total = score.critical_growth_reward +
                  score.supporting_critical_reward +
                  score.vorticity_tiebreak + score.final_l3_tiebreak +
                  score.profile_trend_reward - score.profile_drift_cost -
                  score.profile_rebound_cost - score.cutoff_cost;
    return score;
}

struct PairedScoreEvidence {
    SearchScoreEvidence coarse;
    SearchScoreEvidence fine;
    bool preliminary_cross_resolution_ok = false;
    double l3_ratio_relative_difference = 0.0;
    double h_half_ratio_relative_difference = 0.0;
    double characteristic_scale_relative_difference = 0.0;
    double profile_drift_relative_difference = 1.0;
    bool profile_drift_comparison_valid = false;
    bool coarse_profile_window_recent = false;
    bool fine_profile_window_recent = false;
    double coarse_characteristic_growth = 1.0;
    double fine_characteristic_growth = 1.0;
    double minimum_characteristic_growth = 1.10;
    double minimum_critical_growth = 1.005;
};

struct PairedSearchScore {
    double multi_resolution_cutoff_cost = 0.0;
    double critical_convergence_cost = 0.0;
    double scale_convergence_cost = 0.0;
    double profile_convergence_cost = 0.0;
    double total = 0.0;
    bool refinement_eligible = false;
};

inline PairedSearchScore scoreResolutionPair(
    const PairedScoreEvidence& evidence) {
    const SearchScore coarse_score = scoreSingleResolution(evidence.coarse);
    const SearchScore fine_score = scoreSingleResolution(evidence.fine);
    if (!std::isfinite(evidence.l3_ratio_relative_difference) ||
        evidence.l3_ratio_relative_difference < 0.0 ||
        !std::isfinite(evidence.h_half_ratio_relative_difference) ||
        evidence.h_half_ratio_relative_difference < 0.0 ||
        !std::isfinite(evidence.characteristic_scale_relative_difference) ||
        evidence.characteristic_scale_relative_difference < 0.0 ||
        !finitePositiveScoreValue(evidence.coarse_characteristic_growth) ||
        !finitePositiveScoreValue(evidence.fine_characteristic_growth) ||
        !finitePositiveScoreValue(evidence.minimum_characteristic_growth) ||
        !finitePositiveScoreValue(evidence.minimum_critical_growth)) {
        throw std::invalid_argument("Resolution-pair score evidence is invalid");
    }
    if (evidence.profile_drift_comparison_valid &&
        (!std::isfinite(evidence.profile_drift_relative_difference) ||
         evidence.profile_drift_relative_difference < 0.0)) {
        throw std::invalid_argument("Profile-drift disagreement is invalid");
    }

    PairedSearchScore score;
    const double threshold = evidence.coarse.cutoff_fraction_threshold;
    if (threshold != evidence.fine.cutoff_fraction_threshold) {
        throw std::invalid_argument(
            "Coarse and fine cutoff thresholds must be identical");
    }
    const double coarse_cutoff = safeNonnegativeScoreRatio(
        evidence.coarse.peak_cutoff_fraction, threshold);
    const double fine_cutoff = safeNonnegativeScoreRatio(
        evidence.fine.peak_cutoff_fraction, threshold);
    const double rms_cutoff = safeRootMeanSquare(coarse_cutoff, fine_cutoff);
    score.multi_resolution_cutoff_cost = 0.04 * std::log1p(rms_cutoff);
    score.critical_convergence_cost =
        0.25 * (evidence.l3_ratio_relative_difference +
                evidence.h_half_ratio_relative_difference);
    score.scale_convergence_cost =
        0.10 * evidence.characteristic_scale_relative_difference;
    score.profile_convergence_cost =
        0.05 * (evidence.profile_drift_comparison_valid
                    ? evidence.profile_drift_relative_difference
                    : 1.0);

    // The weaker resolution controls the joint score.  One impressive grid
    // cannot conceal a poor or under-resolved companion run.
    score.total = std::min(coarse_score.total, fine_score.total) -
                  score.multi_resolution_cutoff_cost -
                  score.critical_convergence_cost -
                  score.scale_convergence_cost -
                  score.profile_convergence_cost;

    const bool l3_growth_on_both =
        evidence.coarse.peak_l3_ratio >= evidence.minimum_critical_growth &&
        evidence.fine.peak_l3_ratio >= evidence.minimum_critical_growth;
    const bool h_half_growth_on_both =
        evidence.coarse.peak_h_half_ratio >= evidence.minimum_critical_growth &&
        evidence.fine.peak_h_half_ratio >= evidence.minimum_critical_growth;
    score.refinement_eligible =
        evidence.preliminary_cross_resolution_ok &&
        (l3_growth_on_both || h_half_growth_on_both) &&
        coarse_score.profile_stationarity_improving &&
        fine_score.profile_stationarity_improving &&
        evidence.profile_drift_comparison_valid &&
        evidence.coarse_profile_window_recent &&
        evidence.fine_profile_window_recent &&
        evidence.profile_drift_relative_difference <= 0.50 &&
        evidence.coarse_characteristic_growth >=
            evidence.minimum_characteristic_growth &&
        evidence.fine_characteristic_growth >=
            evidence.minimum_characteristic_growth &&
        evidence.characteristic_scale_relative_difference <= 0.10;
    return score;
}

}  // namespace ns_cascade

#endif
