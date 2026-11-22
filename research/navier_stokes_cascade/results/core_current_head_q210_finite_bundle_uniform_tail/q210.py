"""Q210: unify the finite primitive bundle and isolate all-band uniformity."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q210-FINITE-BUNDLE-UNIFORM-TAIL"
PROVIDERS = (
    "core_current_head_q209_cutoff_role_separation",
    "core_current_head_q207_residual_constant_frontier",
    "core_current_head_q192_frame_graph_jet_transfer",
    "weighted_stress_foundation",
    "leading_frame_foundation",
    "core_same_band_gap_reconstruction",
    "core_current_head_q189_same_profile_envelope_gap_audit",
    "core_direct_compact_slow_provider_reconstruction",
    "core_current_head_q165_finite_map_frontier_audit",
    "core_current_head_q184_generated_source_log_loss_correction",
    "core_current_head_q196_initial_generated_source",
)
EXPECTED = {
    PROVIDERS[0]: "Q209-background-and-stage-cutoff-roles-separated-Q208-identity-blocker-corrected",
    PROVIDERS[1]: "Q207-fixed-stage-complete-residual-constant-reconstructed-conditionally-primitive-seminorm-frontier-isolated",
    PROVIDERS[2]: "Q192-auxiliary-frame-jets-closed-and-physical-graph-loss-corrected",
    PROVIDERS[3]: "actual-weighted-first-stress-and-background-foundation-joined-relative-to-exact-core-and-finite-lemmas",
    PROVIDERS[4]: "current-leading-frame-floors-from-reviewed-core-and-same-root-regional-profile",
    PROVIDERS[5]: "same-band-continuum-gap-contract-reconstructed-on-finalized-leading-profile",
    PROVIDERS[6]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
    PROVIDERS[7]: "direct-compact-slow-operators-proved-with-actual-source-supremum-frontier",
    PROVIDERS[8]: "Q165-current-head-finite-map-and-weighted-source-chain-audited",
    PROVIDERS[9]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
    PROVIDERS[10]: "Q196-initial-generated-source-stage-and-log-quotient-corrected",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q210 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def finite_bundle(order=36479350301, stage=602, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q210 identity required")
    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("nonnegative integer order required")
    if isinstance(stage, bool) or not isinstance(stage, int) or stage < 0:
        raise ValueError("nonnegative integer stage required")
    return {
        "order": order,
        "stage": stage,
        "jet_index_set_finite": True,
        "stage_prefix_finite": True,
        "per_jet_contract": (
            "M_I=sup_(ell>=ell*) epsilon^-B S^-d_I edge_I |D^I(source/P)| < infinity"
        ),
        "unified_majorant": "A_R=max_(|I|<=R) M_I < infinity",
        "gap_bundle": "g_R=min_(g in finite stage-prefix gaps) g > 0",
        "operator_bundle": "K_R=max_(K in finite requested compact/slow archive) K < infinity",
        "finite_maximum_needs_numerical_materialization": False,
        "conditional_on_named_continuum_contracts": True,
        "contracts_independently_reproved_here": False,
    }


def all_band_obstruction(degree=100, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q210 identity required")
    if isinstance(degree, bool) or not isinstance(degree, int) or degree < 0:
        raise ValueError("nonnegative integer degree required")
    n = 8 * (degree + 1)
    ell = n ** 4
    lower = n * (n - 4 * degree)
    return {
        "degree": degree,
        "subsequence": "ell=n^4",
        "smooth_source": "F_ell=epsilon_ell^B exp(sqrt(ell)) phi, phi in C_c^infinity(Omega)",
        "endpoint_jets_vanish": True,
        "each_band_seminorm_finite": True,
        "absolute_source_tends_to_zero_for_positive_B": True,
        "normalized_log": "sqrt(ell)-2d log(ell)=n^2-8d log(n)",
        "elementary_bound": "log(n)<=n/2 gives normalized_log>=n(n-4d)",
        "witness_n": n,
        "witness_ell": ell,
        "witness_lower_bound": lower,
        "lower_bound_positive": lower > 0,
        "normalized_polynomial_supremum_finite": False,
        "counterexample_claimed_for_actual_selected_source": False,
    }


def corrected_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q210 identity required")
    return {
        "Q209_numeric_A_R_is_first_obligation": False,
        "conditional_finite_bundle_follows_from_existing_contracts": True,
        "first_independent_inequality": (
            "sup_(ell>=ell*) epsilon_ell^-B_j S^-d_(j,I) "
            "sup_Omega zeta^-1/2 y_a^r y_b^s |D^I(f_(j,m)/P)| < infinity"
        ),
        "must_hold_for": [
            "same selected profile", "same roots and frame", "complete pressure",
            "radial and axial velocity", "all six covariances", "all supported rows"
        ],
        "Q192_role": "propagates a supplied primitive polynomial-in-S majorant through a finite graph",
        "compact_slow_role": "maps the weighted source space boundedly but cannot prove its own input",
        "per_band_smoothness_suffices": False,
        "tests_or_hashes_prove_uniform_tail": False,
        "unconditional_residual_or_blowup": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q210 identity required")
    status = statuses()
    bundle = finite_bundle()
    obstruction = all_band_obstruction()
    frontier = corrected_frontier()
    checks = {
        "provider_statuses_match": all(status.values()),
        "finite_provider_rows_unify_without_evaluation": (
            bundle["jet_index_set_finite"] and bundle["stage_prefix_finite"]
            and not bundle["finite_maximum_needs_numerical_materialization"]
            and bundle["conditional_on_named_continuum_contracts"]
            and not bundle["contracts_independently_reproved_here"]
        ),
        "smooth_bandwise_countermodel_blocks_invalid_inference": (
            obstruction["endpoint_jets_vanish"]
            and obstruction["each_band_seminorm_finite"]
            and obstruction["absolute_source_tends_to_zero_for_positive_B"]
            and obstruction["lower_bound_positive"]
            and not obstruction["normalized_polynomial_supremum_finite"]
            and not obstruction["counterexample_claimed_for_actual_selected_source"]
        ),
        "frontier_reduced_to_uniform_tail_input": (
            not frontier["Q209_numeric_A_R_is_first_obligation"]
            and frontier["conditional_finite_bundle_follows_from_existing_contracts"]
            and not frontier["per_band_smoothness_suffices"]
            and not frontier["tests_or_hashes_prove_uniform_tail"]
        ),
    }
    return {
        "status": "Q210-finite-primitive-bundle-unified-conditionally-uniform-all-band-source-frontier-isolated",
        "checks": checks,
        "scope": {
            "stage_602_bundle_finite_relative_to_named_contracts": True,
            "numerical_stage_602_bundle_evaluated": False,
            "uniform_actual_source_tail_proved_independently": False,
            "counterexample_to_actual_selected_profile_claimed": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
        "next": (
            "derive the polynomial-in-S all-band primitive source bound from the explicit selected "
            "profile and first generated source, then propagate it through one sourced/generic cycle"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q210 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "bundle": finite_bundle(),
        "countermodels": {
            "d0": all_band_obstruction(0),
            "d100": all_band_obstruction(100),
        },
        "frontier": corrected_frontier(),
        "review": review(),
    }
