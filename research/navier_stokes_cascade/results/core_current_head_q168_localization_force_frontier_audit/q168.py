"""Audit localization, compact force, and energy on the Q164--Q167 object."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q168-current-head-localization-force-energy-frontier-audit"
PROVIDERS = (
    "core_current_head_q167_corner_endpoint_frontier_audit",
    "core_current_head_q8_localization_force_join",
    "core_current_head_q9_preterminal_energy_join",
    "core_current_head_q155_physical_force_energy_join",
    "core_current_head_q162_diagonal_endpoint_force_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q8 = load("q168_q8", "core_current_head_q8_localization_force_join/localization_force.py")
q9 = load("q168_q9", "core_current_head_q9_preterminal_energy_join/preterminal_energy.py")
q155 = load("q168_q155", "core_current_head_q155_physical_force_energy_join/q155.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q168 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q167-current-head-cartesian-corner-and-punctured-endpoint-audited-conditionally",
        PROVIDERS[1]: "current-head-Q8-twelve-term-localization-and-compact-force-joined-conditionally",
        PROVIDERS[2]: "current-head-Q9-preterminal-energy-and-dissipation-joined-conditionally",
        PROVIDERS[3]: "Q155-physical-localization-compact-force-and-energy-bound-conditionally",
        PROVIDERS[4]: "Q162-physical-diagonal-endpoint-localized-force-chain-joined",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q168 provider status mismatch: " + name)

    q7_path = "results/core_current_head_q7_corner_endpoint_join/evidence.json"
    q8_path = "results/core_current_head_q8_localization_force_join/evidence.json"
    q154_path = "results/core_current_head_q154_physical_endpoint_join/evidence.json"
    q155_path = "results/core_current_head_q155_physical_force_energy_join/evidence.json"
    q162_path = "results/core_current_head_q162_diagonal_endpoint_force_join/evidence.json"
    q167, q8_record, q9_record, q155_record, q162_record = (
        data[name] for name in PROVIDERS
    )
    chain = {
        "Q167_to_Q162": (
            q167["dependency_evidence_sha256"][q162_path]
            == sha(PROJECT / q162_path)
        ),
        "Q162_to_Q155": (
            q162_record["dependency_evidence_sha256"][q155_path]
            == sha(PROJECT / q155_path)
        ),
        "Q155_to_Q8": (
            q155_record["dependency_evidence_sha256"][q8_path]
            == sha(PROJECT / q8_path)
        ),
        "Q9_to_Q8": (
            q9_record["dependency_evidence_sha256"][q8_path]
            == sha(PROJECT / q8_path)
        ),
        "Q167_and_Q8_share_Q7": (
            q167["dependency_evidence_sha256"][q7_path]
            == q8_record["dependency_evidence_sha256"][q7_path]
        ),
        "Q162_to_Q154": (
            q162_record["dependency_evidence_sha256"][q154_path]
            == sha(PROJECT / q154_path)
        ),
    }
    return {
        "identity_chain": chain,
        "all_identity_links_hold": all(chain.values()),
        "same_physical_field_cutoff_sequence_support_and_pressure_gauge":
            q162_record["providers"][
                "same_physical_field_cutoff_sequence_support_and_pressure_gauge"
            ],
        "Q154_supplies_Q8_punctured_endpoint_premise": q155.provider_contracts()[
            "Q154_supplies_Q8_punctured_endpoint_premise"
        ],
        "Q155_compact_force_and_energy_conditional": q155_record["symbolic"][
            "compact_force_and_energy_certified_conditionally"
        ],
    }


def localization_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q168 identity required")
    terms = q8.TERM_LEDGER
    ledger = q155.localization_ledger()
    request = q8.force_jet_request(3, 2, 11)
    formulas = tuple(formula for _, formula in terms)
    normalized = tuple(formula.replace(".", " dot ") for formula in formulas)
    return {
        "field": "u=curl(cA)+cB=c u0+d, d=grad(c) cross A; p=c p0",
        "divergence_free": True,
        "zero_initial_velocity": True,
        "protected_late_core_unchanged": True,
        "term_names": tuple(name for name, _ in terms),
        "terms": formulas,
        "term_count": len(terms),
        "Q155_term_count": ledger["count"],
        "Q155_terms_match_after_notation_normalization": normalized == ledger["terms"],
        "vanishing_thirteenth": ledger["vanishing_thirteenth"],
        "unit_viscosity": ledger["viscosity"],
        "complete_pressure": ledger["complete_pressure"],
        "force_divergence_retained": ledger["force_divergence_retained"],
        "sample_request": request,
        "direct_force_orders": request["sharp_direct_force_orders"],
        "Cauchy_orders": request["sharp_Cauchy_orders"],
        "Q167_supplied_orders": request["Q7_supplied_orders"],
        "unused_field_derivative": request["Q7_order_reserve"],
        "residual_Cauchy_jet": request["residual_Cauchy_jet"],
        "pressure_Poisson": "-Delta p=partial_i partial_j(u_i u_j)-div f_minus",
    }


def force_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q168 identity required")
    samples = tuple(q8.extension_request(a, m) for a, m in (
        (0, 0), (1, 3), (20, 4), (6, 70), (100, 100),
    ))
    return {
        "terminal_jets": (
            "F_j=lim_(t up to 1) partial_t^j f_minus in C_c^infinity(K), "
            "with every spatial jet zero at the origin"
        ),
        "former_Q8_endpoint_premise": "closed-conditionally-by-Q154-Q162-on-same-object",
        "formula": "f_plus=sum_j chi0(b_j sigma) sigma^j F_j/j!",
        "scale_rule": (
            "choose increasing dyadic b_j once from the finite terminal norms at index j; "
            "impose |a|+m<=floor(j/2) and a 2^-j bound"
        ),
        "whole_support_bound": (
            "||D_x^a partial_sigma^m T_j||_infinity "
            "<=C_(j,a,m) b_j^(m-j)"
        ),
        "sample_derivatives": samples,
        "all_samples_have_geometric_tails": all(
            all(row["checks"].values()) for row in samples
        ),
        "compact_spacetime_support": "K x [1-tau0,3/2]",
        "raw_Taylor_convergence_used": False,
        "uniform_in_j_norm_assumed": False,
        "future_velocity_or_pressure_constructed": False,
    }


def energy_contract(force_l1=Q(7, 3), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q168 identity required")
    force_l1 = Q(force_l1)
    if force_l1 < 0:
        raise ValueError("nonnegative force integral required")
    q9_budget = q9.fixed_time_budget(force_l1)
    q155_budget = q155.energy_budget(force_l1)
    regularized = q9.regularized_step(Q(1, 1000))
    return {
        "identity": "0.5 E'(t)+D(t)=<f(t),u(t)>",
        "zero_safe_step": regularized,
        "fixed_time": q155_budget["fixed_time"],
        "sup_energy": q9_budget["energy_budget"],
        "integrated_dissipation": q9_budget["dissipation_budget"],
        "Q155_sup_energy": q155_budget["sup_energy"],
        "Q155_integrated_dissipation": q155_budget["integrated_dissipation"],
        "budgets_match": (
            q9_budget["energy_budget"] == q155_budget["sup_energy"]
            and q9_budget["dissipation_budget"]
            == q155_budget["integrated_dissipation"]
        ),
        "factor_two_retained": True,
        "terminal_dissipation_limit": "monotone convergence",
        "combined_endpoint_extrema_budget_claimed": False,
        "terminal_velocity_pressure_or_energy_trace_used": False,
    }


def remaining_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q168 identity required")
    return {
        "U5_localization_force_energy": "closed-conditionally-by-Q168",
        "next": (
            "audit Q156's fixed-gauge pressure recovery and expanding-ball comparison on "
            "the Q164-Q168 object, retaining -div(f) in each individual pressure equation "
            "and cancelling the common force exactly once"
        ),
        "singular_terminal_origin_regularized": False,
        "future_velocity_or_pressure_constructed": False,
        "comparison_or_core_exclusion_replayed_here": False,
        "unconditional_theorem_gate_closed": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q168 identity required")
    return (
        {"id": "U1", "status": "closed-by-Q164-core-pressure-audit"},
        {"id": "U2", "status": "closed-by-Q165-finite-map-audit"},
        {"id": "U3", "status": "closed-conditionally-by-Q166-residual-audit"},
        {"id": "U4", "status": "closed-conditionally-by-Q167-corner-endpoint-audit"},
        {"id": "U5", "status": "closed-conditionally-by-Q168-localization-force-energy-audit"},
        {"id": "U6", "status": "next-Q156-pressure-comparison-audit"},
    )


def record():
    return {
        "status": "Q168-current-head-localization-force-energy-audited-conditionally",
        "correction": (
            "Q8 left its punctured endpoint estimate open, but Q154 proves that estimate "
            "and Q162 binds Q154-Q155 to the same physical field, cutoff sequence, supports "
            "and pressure gauge already audited through Q167"
        ),
        "established": (
            "the Q164-Q167 field admits exact divergence-free localization with zero initial "
            "velocity; its residual has exactly twelve nonzero groups with unit viscosity, "
            "complete pressure, both nonlinear crosses and the curl square; one compact smooth "
            "force matches every terminal force jet, and the factor-two preterminal energy and "
            "separate dissipation estimates follow"
        ),
        "scope": (
            "conditional analytical join of named continuum providers; the force is extended "
            "after t=1 but no future velocity or pressure is constructed, the terminal origin "
            "is not regularized, and tests/hashes are not continuum proofs"
        ),
        "U5_localization_force_energy_closed_conditionally": True,
        "future_velocity_or_pressure_constructed": False,
        "comparison_core_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
        "next": remaining_gate()["next"],
    }
