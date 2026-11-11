"""Reconcile Q164--Q170 with the retained continuum-provider frontier."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q171-current-head-provider-discharge-frontier-audit"
PROVIDERS = (
    "core_current_head_q170_protected_core_limit_frontier_audit",
    "core_current_head_q169_pressure_comparison_frontier_audit",
    "core_current_head_q168_localization_force_frontier_audit",
    "core_current_head_q167_corner_endpoint_frontier_audit",
    "core_current_head_q166_complete_residual_frontier_audit",
    "core_current_head_q165_finite_map_frontier_audit",
    "core_current_head_q164_core_fixed_point_pressure_audit",
    "comparison_foundation",
    "core_current_head_primitive_frontier_audit",
    "core_current_head_q2_regional_gap_join",
    "core_current_head_q3_full_frame_join",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q171 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q170-current-head-protected-core-and-limit-order-audited-conditionally",
        PROVIDERS[1]: "Q169-current-head-fixed-gauge-pressure-comparison-audited-conditionally",
        PROVIDERS[2]: "Q168-current-head-localization-force-energy-audited-conditionally",
        PROVIDERS[3]: "Q167-current-head-cartesian-corner-and-punctured-endpoint-audited-conditionally",
        PROVIDERS[4]: "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
        PROVIDERS[5]: "Q165-current-head-finite-map-and-weighted-source-chain-audited",
        PROVIDERS[6]: "Q164-all-index-core-contraction-and-leading-pressure-audited",
        PROVIDERS[7]: "current-source-compact-pressure-comparison-relative-to-retained-core-and-finite-analytic-lemmas",
        PROVIDERS[8]: "current-head-five-contract-frontier-corrected-to-lower-core-root-regional-frame-antecedents",
        PROVIDERS[9]: "current-head-Q2-regional-continuum-factored-edges-and-same-band-gap-bound",
        PROVIDERS[10]: "current-head-Q3-full-frame-whole-path-and-phase-domain-bound",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q171 provider status mismatch: " + name)

    paths = {
        number: f"results/core_current_head_q{number}_{suffix}/evidence.json"
        for number, suffix in (
            (169, "pressure_comparison_frontier_audit"),
            (168, "localization_force_frontier_audit"),
            (167, "corner_endpoint_frontier_audit"),
            (166, "complete_residual_frontier_audit"),
            (165, "finite_map_frontier_audit"),
            (164, "core_fixed_point_pressure_audit"),
        )
    }
    chain = {}
    for parent, child in zip(range(170, 164, -1), range(169, 163, -1)):
        parent_record = data[PROVIDERS[170 - parent]]
        chain[f"Q{parent}_to_Q{child}"] = (
            parent_record["dependency_evidence_sha256"][paths[child]]
            == sha(PROJECT / paths[child])
        )

    comparison = data[PROVIDERS[7]]["flags"]
    q164 = data[PROVIDERS[6]]["symbolic"]
    q165 = data[PROVIDERS[5]]["symbolic"]
    q2_scope = data[PROVIDERS[9]]["review"]["scope"]
    q3_scope = data[PROVIDERS[10]]["review"]["scope"]
    return {
        "identity_chain": chain,
        "all_identity_links_hold": all(chain.values()),
        "historical_comparison_flag_all_dependencies_discharged": comparison[
            "all_retained_analytic_dependencies_discharged"
        ],
        "historical_comparison_flag_core_and_finite_maps_still_dependencies": comparison[
            "retained_exact_core_and_finite_existence_inverse_gap_lemmas_still_dependencies"
        ],
        "Q164_later_closes_current_head_core_pressure_frontier": q164[
            "U1_core_contraction_and_leading_pressure_closed"
        ],
        "Q165_later_closes_current_head_finite_map_frontier": q165[
            "U2_finite_map_and_weighted_source_chain_closed"
        ],
        "Q2_bound_only_relative_to_named_providers": (
            q2_scope["Q2_regional_continuum_and_factored_edge_chain_bound_on_that_object"]
            and not q2_scope["every_provider_theorem_independently_reproved_here"]
        ),
        "Q3_bound_only_relative_to_named_providers": (
            q3_scope["Q3_full_frame_primitive_path_and_phase_domain_bound_on_same_later_tail"]
            and not q3_scope["every_provider_theorem_independently_reproved_here"]
        ),
    }


def successor_discharge(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q171 identity required")
    return {
        "comparison_foundation_historical_scope": (
            "accurate at its own head; historical evidence is not rewritten"
        ),
        "exact_core_and_leading_pressure": "later-current-head-audit-Q164",
        "finite_regional_frame_projector_inverse_chain": "later-current-head-audit-Q165",
        "complete_local_residual": "later-conditional-audit-Q166",
        "corner_and_punctured_endpoint": "later-conditional-audit-Q167",
        "localization_compact_force_and_energy": "later-conditional-audit-Q168",
        "fixed_gauge_pressure_and_comparison": "later-conditional-audit-Q169",
        "protected_core_and_limit_order": "later-conditional-audit-Q170",
        "old_comparison_core_finite_map_flag_is_current_frontier": False,
        "later_audits_make_every_named_provider_unconditional": False,
    }


def selector_contract(margins=(Q(1, 2), Q(1, 4)),
                      constants=(Q(4), Q(2)), h=Q(1, 200),
                      identity=IDENTITY):
    """Exact finite-selector implication for h=1/200.

    Raising 2^(2h(ell-2)) >= C/r to the 100th power avoids a
    floating logarithm.  It proves finiteness only after every actual
    continuum margin and constant has been supplied.
    """
    if identity != IDENTITY:
        raise ValueError("Q171 identity required")
    h = Q(h)
    if h != Q(1, 200):
        raise ValueError("the retained h=1/200 is required")
    try:
        margins = tuple(Q(value) for value in margins)
        constants = tuple(Q(value) for value in constants)
    except (OverflowError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("all selector inputs must be finite real numbers") from exc
    if not margins or len(margins) != len(constants):
        raise ValueError("paired nonempty margins and constants required")
    if any(value <= 0 for value in margins):
        raise ValueError("all actual continuum margins must be positive")
    if any(value <= 0 for value in constants):
        raise ValueError("all actual continuum constants must be positive and finite")
    ratios = tuple(constant / margin for constant, margin in zip(constants, margins))
    targets = tuple(ratio ** 100 for ratio in ratios)
    ell = 2
    target = max(targets)
    while Q(2) ** (ell - 2) < target:
        ell += 1
    return {
        "h": h,
        "m_base": min(margins),
        "ratios": ratios,
        "raised_power": 100,
        "ell": ell,
        "exact_selector_holds": all(
            Q(2) ** (ell - 2) >= item for item in targets
        ),
        "previous_integer_fails_when_ell_is_nontrivial": (
            ell == 2 or Q(2) ** (ell - 3) < target
        ),
    }


def first_open_premise(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q171 identity required")
    primitive = evidence(PROVIDERS[8])["review"]
    q2_scope = evidence(PROVIDERS[9])["review"]["scope"]
    q3_scope = evidence(PROVIDERS[10])["review"]["scope"]
    return {
        "first_not_independently_reproved_current_head_premise": (
            "on the selected Q164 core profile, independently establish the actual "
            "Q2 regional continuum margins and Q3 finite whole-path/frame constants"
        ),
        "required_margin": (
            "m_base=min(f,t,d,omega,m_phi,c_a,b_0,lambda_min,b_min^2)>0"
        ),
        "required_constants": "0<C_i<infinity for every value, relative-stress, chart, length and frame row",
        "blocking_selector": "2h(ell-2)>=max_i log2(C_i/r_i)",
        "logical_obstruction": primitive["blocking_inequality"]["logical_status"],
        "Q2_provider_chain_bound_but_not_fully_reproved_in_join": (
            q2_scope["Q2_regional_continuum_and_factored_edge_chain_bound_on_that_object"]
            and not q2_scope["every_provider_theorem_independently_reproved_here"]
        ),
        "Q3_provider_chain_bound_but_not_fully_reproved_in_join": (
            q3_scope["Q3_full_frame_primitive_path_and_phase_domain_bound_on_same_later_tail"]
            and not q3_scope["every_provider_theorem_independently_reproved_here"]
        ),
        "smaller_q_alone_can_prove_missing_positivity_or_finiteness": False,
        "positive_finite_inputs_imply_a_finite_selector": selector_contract()[
            "exact_selector_holds"
        ],
        "provider_claims_declared_false": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q171 identity required")
    return (
        {"id": "U1-U7", "status": "successor-audits-reconciled-relative-to-named-providers"},
        {"id": "U8", "status": "closed-by-Q171-provider-frontier-reconciliation"},
        {"id": "U9", "status": "next-independent-Q2-Q3-continuum-provider-audit"},
    )


def record():
    return {
        "status": "Q171-current-head-provider-frontier-reconciled-with-selector-obstruction",
        "correction": (
            "comparison_foundation's old core/finite-map dependency flag is historical, "
            "not the Q170 frontier: Q164 and Q165 later audit those interfaces on the "
            "current-head object, while Q166-Q170 close the downstream implications only "
            "relative to their named continuum providers"
        ),
        "established": (
            "the Q164-Q170 successor chain is byte-identical and acyclic; every current-head "
            "interface through comparison and the protected core is joined, and positive "
            "finite Q2/Q3 inputs would give a finite selector by an exact exponent-100 argument"
        ),
        "obstruction": (
            "the current-head joins do not independently reprove every Q2 regional continuum "
            "infimum or Q3 whole-path/frame constant on the selected profile; until "
            "m_base>0 and all C_i are finite are independently validated, "
            "2h(ell-2)>=max_i log2(C_i/r_i) does not define a finite ell"
        ),
        "scope": (
            "provider-frontier and exact selector audit; it neither falsifies the named "
            "provider theorems nor upgrades tests and frozen hashes into continuum proofs"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
        "next": (
            "independently audit the actual Q2 ten-region continuum infima and factored "
            "flat-edge coefficients, then the Q3 whole-path/frame constants, on the exact "
            "Q164 selected profile before rerunning the common selector and downstream chain"
        ),
    }
