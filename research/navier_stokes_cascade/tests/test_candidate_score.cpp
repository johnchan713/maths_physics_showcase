#include "ns_cascade/candidate_score.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void expectNear(double actual,
                double expected,
                double tolerance,
                const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

ns_cascade::SearchScoreEvidence usefulEvidence() {
    ns_cascade::SearchScoreEvidence evidence;
    evidence.peak_l3_ratio = 1.01;
    evidence.peak_h_half_ratio = 1.08;
    evidence.final_l3_ratio = 0.96;
    evidence.peak_vorticity_ratio = 1.8;
    evidence.peak_cutoff_fraction = 0.001;
    evidence.cutoff_fraction_threshold = 0.01;
    evidence.profile_scale_windows = 3;
    evidence.first_profile_drift = 0.8;
    evidence.latest_profile_drift = 0.4;
    evidence.minimum_profile_drift = 0.4;
    evidence.profile_drift_threshold = 1.0;
    return evidence;
}

void testSingleResolutionRewardsStationarity() {
    const ns_cascade::SearchScoreEvidence improving = usefulEvidence();
    const ns_cascade::SearchScore improving_score =
        ns_cascade::scoreSingleResolution(improving);
    expect(improving_score.profile_stationarity_improving,
           "Improving low drift was not recognized");
    expect(improving_score.profile_trend_reward > 0.0,
           "Decreasing drift did not earn a positive trend term");

    ns_cascade::SearchScoreEvidence worsening = improving;
    worsening.latest_profile_drift = 1.2;
    const ns_cascade::SearchScore worsening_score =
        ns_cascade::scoreSingleResolution(worsening);
    expect(!worsening_score.profile_stationarity_improving,
           "Large increasing drift incorrectly passed stationarity");
    expect(improving_score.total > worsening_score.total,
           "The score did not prefer a more stationary profile");

    ns_cascade::SearchScoreEvidence missing = improving;
    missing.profile_scale_windows = 0;
    missing.first_profile_drift = 0.0;
    missing.latest_profile_drift = 0.0;
    missing.minimum_profile_drift = 0.0;
    const ns_cascade::SearchScore missing_score =
        ns_cascade::scoreSingleResolution(missing);
    expect(!missing_score.profile_stationarity_improving,
           "Missing profile evidence incorrectly passed stationarity");
    expect(missing_score.profile_drift_cost > improving_score.profile_drift_cost,
           "Missing profile evidence was treated as zero drift");

    ns_cascade::SearchScoreEvidence rebounding = improving;
    rebounding.minimum_profile_drift = 0.2;
    const ns_cascade::SearchScore rebounding_score =
        ns_cascade::scoreSingleResolution(rebounding);
    expect(!rebounding_score.profile_stationarity_improving,
           "A materially rebounding final drift passed stationarity");
    expect(rebounding_score.profile_rebound_cost > 0.0 &&
               improving_score.total > rebounding_score.total,
           "A late profile-drift rebound was not penalized");
}

void testSingleResolutionPenalizesCutoffLoading() {
    const ns_cascade::SearchScoreEvidence clean = usefulEvidence();
    ns_cascade::SearchScoreEvidence loaded = clean;
    loaded.peak_cutoff_fraction = 0.009;
    const ns_cascade::SearchScore clean_score =
        ns_cascade::scoreSingleResolution(clean);
    const ns_cascade::SearchScore loaded_score =
        ns_cascade::scoreSingleResolution(loaded);
    expect(clean_score.cutoff_cost < loaded_score.cutoff_cost,
           "Cutoff cost did not increase with cutoff loading");
    expect(clean_score.total > loaded_score.total,
           "The score did not prefer the cleaner spectrum");
}

ns_cascade::PairedScoreEvidence usefulPair() {
    ns_cascade::PairedScoreEvidence pair;
    pair.coarse = usefulEvidence();
    pair.fine = usefulEvidence();
    pair.fine.peak_cutoff_fraction = 0.0001;
    pair.preliminary_cross_resolution_ok = true;
    pair.l3_ratio_relative_difference = 0.002;
    pair.h_half_ratio_relative_difference = 0.003;
    pair.characteristic_scale_relative_difference = 0.02;
    pair.profile_drift_relative_difference = 0.10;
    pair.profile_drift_comparison_valid = true;
    pair.coarse_profile_window_recent = true;
    pair.fine_profile_window_recent = true;
    pair.coarse_characteristic_growth = 1.20;
    pair.fine_characteristic_growth = 1.22;
    pair.minimum_characteristic_growth = 1.10;
    pair.minimum_critical_growth = 1.005;
    return pair;
}

void testPairGateAndJointCutoffCost() {
    const ns_cascade::PairedScoreEvidence eligible = usefulPair();
    const ns_cascade::PairedSearchScore eligible_score =
        ns_cascade::scoreResolutionPair(eligible);
    expect(eligible_score.refinement_eligible,
           "A converged growing stationary pair failed the refinement gate");

    ns_cascade::PairedScoreEvidence coarse_loaded = eligible;
    coarse_loaded.coarse.peak_cutoff_fraction = 0.009;
    const ns_cascade::PairedSearchScore coarse_loaded_score =
        ns_cascade::scoreResolutionPair(coarse_loaded);
    ns_cascade::PairedScoreEvidence fine_loaded = eligible;
    fine_loaded.fine.peak_cutoff_fraction = 0.009;
    const ns_cascade::PairedSearchScore fine_loaded_score =
        ns_cascade::scoreResolutionPair(fine_loaded);
    expect(eligible_score.multi_resolution_cutoff_cost <
               coarse_loaded_score.multi_resolution_cutoff_cost &&
               eligible_score.multi_resolution_cutoff_cost <
                   fine_loaded_score.multi_resolution_cutoff_cost,
           "Joint cutoff cost ignored a coarse or fine resolution");
    expect(eligible_score.total > coarse_loaded_score.total &&
               eligible_score.total > fine_loaded_score.total,
           "Pair score did not prefer lower cutoff loading on both grids");

    ns_cascade::PairedScoreEvidence disagreeing = eligible;
    disagreeing.profile_drift_relative_difference = 0.51;
    const ns_cascade::PairedSearchScore disagreeing_score =
        ns_cascade::scoreResolutionPair(disagreeing);
    expect(!disagreeing_score.refinement_eligible,
           "Disagreeing profile drifts incorrectly passed refinement");

    ns_cascade::PairedScoreEvidence weak_scale = eligible;
    weak_scale.coarse_characteristic_growth = 1.09;
    const ns_cascade::PairedSearchScore weak_scale_score =
        ns_cascade::scoreResolutionPair(weak_scale);
    expect(!weak_scale_score.refinement_eligible,
           "Insufficient spectral-scale movement incorrectly passed refinement");

    ns_cascade::PairedScoreEvidence stale = eligible;
    stale.fine_profile_window_recent = false;
    const ns_cascade::PairedSearchScore stale_score =
        ns_cascade::scoreResolutionPair(stale);
    expect(!stale_score.refinement_eligible,
           "A stale final scale window incorrectly passed refinement");
}

void testInvalidEvidenceIsRejected() {
    ns_cascade::SearchScoreEvidence invalid = usefulEvidence();
    invalid.peak_h_half_ratio = 0.0;
    bool rejected = false;
    try {
        ns_cascade::scoreSingleResolution(invalid);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    expect(rejected, "A zero diagnostic ratio was accepted by the score");

    const ns_cascade::SearchScoreEvidence evidence = usefulEvidence();
    const ns_cascade::SearchScore score =
        ns_cascade::scoreSingleResolution(evidence);
    expectNear(score.total,
               score.critical_growth_reward +
                   score.supporting_critical_reward +
                   score.vorticity_tiebreak + score.final_l3_tiebreak +
                   score.profile_trend_reward - score.profile_drift_cost -
                   score.profile_rebound_cost - score.cutoff_cost,
               1e-15,
               "Reported score components do not sum to the total");

    ns_cascade::SearchScoreEvidence extreme = evidence;
    extreme.cutoff_fraction_threshold =
        std::numeric_limits<double>::denorm_min();
    const ns_cascade::SearchScore extreme_score =
        ns_cascade::scoreSingleResolution(extreme);
    expect(std::isfinite(extreme_score.total),
           "An extreme positive threshold overflowed the score");
}

}  // namespace

int main() {
    try {
        testSingleResolutionRewardsStationarity();
        testSingleResolutionPenalizesCutoffLoading();
        testPairGateAndJointCutoffCost();
        testInvalidEvidenceIsRejected();
        std::cout << "All candidate-score tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
