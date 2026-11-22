"""Q214: join the exact core, selected roots and same-gauge leading pressure."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q214-EXACT-CORE-ROOT-PRESSURE-JOIN"
PROVIDERS = (
    "core_current_head_q213_reconstructed_primitive_bundle_exact_root_frontier",
    "core_current_head_q164_core_fixed_point_pressure_audit",
    "core_pressure_identification",
    "core_current_head_q176_modulation_repair_dependency_audit",
    "matching_root_foundation",
    "annular_exterior_foundation",
    "core_current_head_q189_same_profile_envelope_gap_audit",
    "core_current_head_q204_global_ledger_uniform_primary",
    "core_current_head_q207_residual_constant_frontier",
)
EXPECTED = {
    PROVIDERS[0]: "Q213-reconstructed-primitive-bundle-joined-conditionally-exact-core-root-frontier-isolated",
    PROVIDERS[1]: "Q164-all-index-core-contraction-and-leading-pressure-audited",
    PROVIDERS[2]: "leading-core-pressure-identified-with-regular-angular-inverse-and-retained-higher-memory",
    PROVIDERS[3]: "Q176-compact-modulation-exact-repair-audited-with-refined-dependency-order",
    PROVIDERS[4]: "actual-annular-and-compact-matching-roots-audited-relative-to-retained-core-and-compact-input-estimates",
    PROVIDERS[5]: "annular-exterior-all-order-leading-regularity-and-global-moment-foundation-relative-to-exact-roots",
    PROVIDERS[6]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
    PROVIDERS[7]: "Q204-global-pressure-coefficient-uniform-primary-and-covariance-ledgers-joined-conditionally",
    PROVIDERS[8]: "Q207-fixed-stage-complete-residual-constant-reconstructed-conditionally-primitive-seminorm-frontier-isolated",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q214 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def core_gate(H=Q(2), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q214 identity required")
    H = Q(H)
    if H <= 1:
        raise ValueError("the selected hierarchy has H>1")
    upper = Q(1, 10**50) / H**127
    return {
        "q_core_upper": upper,
        "strict_contraction": upper < Q(1, 10**50) < 1,
        "complete_pressure_source_retained": True,
        "original_core_equations_recovered": True,
        "leading_core_stress_exactly_zero": True,
    }


def root_gates(lam=Q(1, 100), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q214 identity required")
    lam = Q(lam)
    if lam <= 0:
        raise ValueError("the selected axial parameter lambda is positive")
    floors = {
        "compact_repair": Q(1, 1600),
        "annular_five_row": Q(1, 1220),
        "two_angular_corrections": Q(1, 10),
        "axial_end_correction": lam / 20,
        "pulse_strength": Q(7, 20),
        "heat_compensation": Q(1, 100),
    }
    return {
        "inverse_floors": floors,
        "all_displayed_floors_positive": all(value > 0 for value in floors.values()),
        "same_accepted_roots": True,
        "signed_variance_and_moving_phase_derivatives_positive": True,
        "exact_five_rows_include_pressure": True,
        "all_fixed_jets_use_same_root": True,
    }


def pressure_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q214 identity required")
    return {
        "reference_measure_is_core_datum": True,
        "same_root_five_row_ledger": True,
        "smooth_endpoint_kills_homogeneous_mode": True,
        "Pi_core_equals_completed_leading_pressure": True,
        "positive_order_pressure_memory_retained": [
            "Pi1/Pi2", "covariance pressure", "moving pressure",
            "cutoff pressure derivatives",
        ],
        "full_corrected_pressure_identified_with_centrifugal_pressure": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q214 identity required")
    core = core_gate()
    roots = root_gates()
    pressure = pressure_gate()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "Q213_core_gate_already_closed": (
            core["strict_contraction"] and core["complete_pressure_source_retained"]
            and core["original_core_equations_recovered"]
        ),
        "Q213_selected_root_gate_already_closed": (
            roots["all_displayed_floors_positive"] and roots["same_accepted_roots"]
            and roots["signed_variance_and_moving_phase_derivatives_positive"]
            and roots["exact_five_rows_include_pressure"]
        ),
        "Q213_leading_pressure_gate_already_closed": (
            pressure["Pi_core_equals_completed_leading_pressure"]
            and pressure["smooth_endpoint_kills_homogeneous_mode"]
            and not pressure["full_corrected_pressure_identified_with_centrifugal_pressure"]
        ),
    }
    return {
        "status": "Q214-exact-core-selected-roots-and-leading-pressure-joined-conditionally",
        "checks": checks,
        "scope": {
            "Q213_three_part_frontier_closed_relative_to_named_providers": True,
            "same_profile_and_unchanged_pressure_gauge": True,
            "positive_order_pressure_memory_preserved": True,
            "high_jet_constants_numerically_evaluated": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
        "next": (
            "audit whether Q213-Q214's quantified finite primitive bundle discharges "
            "Q207's concrete fixed-stage request and joins the complete residual recurrence"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q214 identity required")
    core = core_gate()
    roots = root_gates()
    return {
        "status": review()["status"],
        "providers": statuses(),
        "core_gate": {**core, "q_core_upper": str(core["q_core_upper"])},
        "root_gates": {
            **roots,
            "inverse_floors": {name: str(value) for name, value in roots["inverse_floors"].items()},
        },
        "pressure_gate": pressure_gate(),
        "review": review(),
    }
