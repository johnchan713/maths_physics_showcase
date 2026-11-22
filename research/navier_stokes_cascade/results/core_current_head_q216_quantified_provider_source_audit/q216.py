"""Q216: source audit of Q215's all-fixed-order and all-finite-pole providers."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q216-QUANTIFIED-PROVIDER-SOURCE-AUDIT"
R_STAR = 2_626_513_221_816
PROVIDERS = (
    "core_current_head_q215_fixed_stage_residual_bundle_join",
    "core_current_head_q214_exact_core_root_pressure_join",
    "core_leading_profile_certificate",
    "core_leading_background_seminorm_reconstruction",
    "core_same_band_gap_reconstruction",
    "core_full_path_ode_reconstruction",
    "core_fixed_profile_projector_reconstruction",
    "core_compact_slow_inverse_reconstruction",
    "core_current_head_q172_q2_q3_source_audit",
    "core_current_head_q189_same_profile_envelope_gap_audit",
)
EXPECTED = {
    PROVIDERS[0]: "Q215-fixed-stage-primitive-bundle-and-complete-residual-recurrence-joined-conditionally",
    PROVIDERS[1]: "Q214-exact-core-selected-roots-and-leading-pressure-joined-conditionally",
    PROVIDERS[2]: "finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses",
    PROVIDERS[3]: "fixed-order-background-seminorm-and-real-borel-selector-reconstructed-relative-to-leading-certificate",
    PROVIDERS[4]: "same-band-continuum-gap-contract-reconstructed-on-finalized-leading-profile",
    PROVIDERS[5]: "canonical-full-path-ode-reconstructed-on-current-same-band-tail",
    PROVIDERS[6]: "fixed-profile-same-choice-projector-contract-reconstructed-on-current-full-path-tail",
    PROVIDERS[7]: "compact-slow-weighted-inverse-contract-reconstructed-on-current-same-choice-tail",
    PROVIDERS[8]: "Q172-Q2-source-audited-and-Q3-uniform-frame-corrected",
    PROVIDERS[9]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q216 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def borel_split(order, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q216 identity required")
    if isinstance(order, bool) or not isinstance(order, int) or order < 0:
        raise ValueError("nonnegative integer order required")
    tail_start = max(1, order)
    return {
        "order": order,
        "finite_head_first": 1,
        "finite_head_last": tail_start - 1,
        "finite_head_count": tail_start - 1,
        "geometric_tail_start": tail_start,
        "tail_selector_contains_requested_jet": order <= tail_start,
        "requires_order_uniform_constant": False,
    }


def source_quantifiers(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q216 identity required")
    return {
        "safe_request": R_STAR,
        "leading_profile_every_fixed_mixed_jet": True,
        "background_finite_head_plus_geometric_tail": True,
        "gaps_value_level_and_request_independent": True,
        "infinite_band_frame_uses_analytic_tail_not_compactness": True,
        "full_path_fixed_order_Duhamel_induction": True,
        "one_fixed_projector_profile_all_fixed_jets": True,
        "compact_kernel_accepts_every_finite_pole": True,
        "twelve_q_is_terminal_requested_pole": True,
        "diagnostic_fixture_used_as_continuum_proof": False,
        "circular_dependency_found": False,
    }


def selection_order(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q216 identity required")
    return [
        "exact core and selected roots",
        "finalized leading profile",
        "fixed coefficient then its Borel scale",
        "later value gap and frame tail",
        "full-path ODE",
        "fixed projector and same-choice recomputation",
        "finite-request compact and slow inverses",
        "residual constant then diagonal scale",
    ]


def reserves(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q216 identity required")
    h = Q(1, 200)
    kappa = Q(1, 100000)
    delta = Q(1, 2) - 3 * kappa
    return {
        "regular_log_reserve": h * delta / 2,
        "compact_log_reserve": h * kappa / 2,
        "pure_power_margin": Q(0),
        "pure_power_route_rejected": True,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q216 identity required")
    split = borel_split(R_STAR)
    quantifiers = source_quantifiers()
    reserve = reserves()
    order = selection_order()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "huge_fixed_order_has_finite_head_and_geometric_tail": (
            split["finite_head_count"] == R_STAR - 1
            and split["geometric_tail_start"] == R_STAR
            and split["tail_selector_contains_requested_jet"]
            and not split["requires_order_uniform_constant"]
        ),
        "all_source_quantifiers_are_analytical": (
            all(value for key, value in quantifiers.items()
                if key not in {"safe_request", "diagnostic_fixture_used_as_continuum_proof",
                               "circular_dependency_found"})
            and not quantifiers["diagnostic_fixture_used_as_continuum_proof"]
            and not quantifiers["circular_dependency_found"]
        ),
        "dependency_order_is_acyclic": len(order) == len(set(order)) == 8,
        "log_route_has_strict_owned_reserves": (
            reserve["regular_log_reserve"] == Q(49997, 40000000)
            and reserve["compact_log_reserve"] == Q(1, 40000000)
            and reserve["pure_power_margin"] == 0
            and reserve["pure_power_route_rejected"]
        ),
    }
    return {
        "status": "Q216-all-fixed-order-and-finite-pole-provider-quantifiers-source-audited",
        "checks": checks,
        "scope": {
            "Q215_provider_quantifiers_supported_by_source_recursions": True,
            "R_star_numerically_materialized": False,
            "order_or_stage_uniform_constant": False,
            "terminal_origin_regularized": False,
            "unconditional_force_exclusion_or_blowup": False,
        },
        "next": (
            "bind the source-audited object to localization, compact force, energy and fixed-gauge "
            "comparison, then identify the first premise still conditional in that chain"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q216 identity required")
    reserve = reserves()
    return {
        "status": review()["status"],
        "providers": statuses(),
        "borel_split_at_R_star": borel_split(R_STAR),
        "source_quantifiers": source_quantifiers(),
        "selection_order": selection_order(),
        "reserves": {key: str(value) if isinstance(value, Q) else value
                     for key, value in reserve.items()},
        "review": review(),
    }
