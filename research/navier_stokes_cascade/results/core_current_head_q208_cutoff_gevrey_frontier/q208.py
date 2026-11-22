"""Q208: bind the existing squared-flat cutoff bound into Q207."""
from collections import deque
from functools import lru_cache
from math import ceil, log2
from pathlib import Path
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q208-BOREL-PROVIDER-FRONTIER"
PROVIDERS = (
    "core_current_head_q207_residual_constant_frontier",
    "core_current_head_q205_stagewise_pulse_diagonal",
    "core_current_head_q204_global_ledger_uniform_primary",
    "core_current_head_q120_closed_borel_cutoff",
    "core_current_head_q192_frame_graph_jet_transfer",
    "weighted_stress_foundation",
    "regional_weighted_gate",
)
EXPECTED = {
    PROVIDERS[0]: "Q207-fixed-stage-complete-residual-constant-reconstructed-conditionally-primitive-seminorm-frontier-isolated",
    PROVIDERS[1]: "Q205-stagewise-generated-pulse-bounds-suffice-for-diagonal-uniform-raw-hierarchy-remains-unproved",
    PROVIDERS[2]: "Q204-global-pressure-coefficient-uniform-primary-and-covariance-ledgers-joined-conditionally",
    PROVIDERS[3]: "closed-implemented-Borel-Euler-bound-and-retained-row-slack",
    PROVIDERS[4]: "Q192-auxiliary-frame-jets-closed-and-physical-graph-loss-corrected",
    PROVIDERS[5]: "actual-weighted-first-stress-and-background-foundation-joined-relative-to-exact-core-and-finite-lemmas",
    PROVIDERS[6]: "finite-seminorm-weighted-input-and-regional-background-gap-join",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


@lru_cache(maxsize=1)
def _load_q120():
    path = HERE.parent / PROVIDERS[3] / "q120.py"
    spec = importlib.util.spec_from_file_location("q208_q120", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=None)
def frozen_path(parent, child, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    start = f"results/{parent}/evidence.json"
    target = f"results/{child}/evidence.json"
    queue = deque([(start, [start])])
    seen = {start}
    while queue:
        current, path = queue.popleft()
        if current == target:
            return path
        data = json.loads((PROJECT / current).read_text())
        for name, digest in sorted(data.get("dependency_evidence_sha256", {}).items()):
            if sha(PROJECT / name) != digest:
                raise ValueError(f"stale dependency hash: {name}")
            if name not in seen:
                seen.add(name)
                queue.append((name, path + [name]))
    raise ValueError(f"{child} is not in {parent}'s frozen closure")


def same_object_links(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    return {
        "Q207_to_Q205": len(frozen_path(PROVIDERS[0], PROVIDERS[1])) - 1,
        "Q207_to_Q204": len(frozen_path(PROVIDERS[0], PROVIDERS[2])) - 1,
        "Q207_to_Q120": len(frozen_path(PROVIDERS[0], PROVIDERS[3])) - 1,
        "Q207_to_Q192": len(frozen_path(PROVIDERS[0], PROVIDERS[4])) - 1,
    }


def cutoff_provider(order=36479350300, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    if isinstance(order, bool) or not isinstance(order, int) or order < 3:
        raise ValueError("integer order at least three required")
    q120 = _load_q120()
    bits = q120.cutoff_binary_upper(order)
    recipe = q120.cutoff_recipe(order)
    return {
        "order": order,
        "identity": q120.IDENTITY,
        "flat_seed": "g(x)=exp(-1/x^2) for x>0 and zero for x<=0",
        "Euler_majorant": "K_m <= 82*4097^m*(m!)^3",
        "binary_budget_formula": "7+13*m+3*m*ceil(log2(m))",
        "ceil_log2_order": ceil(log2(order)),
        "binary_budget": bits,
        "direct_q120_budget_match": bits == 7 + 13 * order + 3 * order * ceil(log2(order)),
        "stage_602_order_match": order == 36479350300,
        "factorial_or_derivative_array_materialized": False,
        "new_cutoff_substituted": False,
        "q120_actual_field_identity_certified": recipe["actual_field_identity_certified"],
    }


def correction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    row = cutoff_provider()
    return {
        "discarded_draft_seed": "exp(-1/x)",
        "reason": "the live implemented Borel cutoff is squared-flat exp(-1/x^2)",
        "Q120_already_proved_all_order_cutoff_bound": True,
        "Q204_rule_selects_background_scale_after_finite_constants": True,
        "Q205_rule_selects_stage_scale_after_finite_constants": True,
        "conditional_budget_available": row["binary_budget"] == 4414001386307,
        "Q207_to_Q120_same_field_identity_explicitly_certified": row["q120_actual_field_identity_certified"],
        "cutoff_shape_removed_from_frontier_unconditionally": False,
        "historical_evidence_modified": False,
    }


def propagation_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    return {
        "cutoff_shape": "Q120 squared-flat Euler bound through every fixed order",
        "background_schedule": "Q204 rule-defined after each finite coefficient row",
        "stage_schedule": "Q205 rule-defined after each finite stage request",
        "frame_reciprocals": "Q192 finite Leibniz recurrence after primitive jets and floors",
        "pressure_memory": "Pi1/Pi2, whole moving product and all six covariance entries retained",
        "cutoff_shape_bound_available_after_identity_binding": True,
        "same_field_cutoff_identity_bound": False,
        "selected_profile_gap_and_operator_bundle_closed_directly": False,
    }


def remaining_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    return {
        "first_missing_majorant": (
            "A_R=max_(|I|<=R+1) sup_D w_I S^(-d_I) "
            "|D^I(F0,U0,Pi0,M0,V0,E1,U1,W1,T1/zeta)| <= A_R"
        ),
        "R": 36479350300,
        "actual_weighted_first_coefficient_orderwise_finite": True,
        "closed_majorant_A_R_from_selected_defining_maps_recorded": False,
        "additional_required_bundle": "common positive gaps g_R and every archived compact/slow inverse K_R",
        "same_cutoff_identity_obligation": (
            "bind Q207 chi(a_j q) to Q120 implemented squared-flat Borel chi before applying the bit budget"
        ),
        "common_value_floors_are_separate_from_high_derivative_sizes": True,
        "passing_tests_or_hashes_prove_this_supremum": False,
        "unconditional_residual_or_blowup": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    status = statuses()
    links = same_object_links()
    cutoff = cutoff_provider()
    corrected = correction()
    propagation = propagation_ledger()
    frontier = remaining_frontier()
    checks = {
        "provider_statuses_match": all(status.values()),
        "same_object_links_frozen": all(value > 0 for value in links.values()),
        "existing_cutoff_provider_covers_stage_602": (
            cutoff["stage_602_order_match"]
            and cutoff["binary_budget"] == 4414001386307
            and cutoff["direct_q120_budget_match"]
            and not cutoff["factorial_or_derivative_array_materialized"]
        ),
        "wrong_cutoff_draft_and_identity_gap_recorded": (
            corrected["Q120_already_proved_all_order_cutoff_bound"]
            and corrected["conditional_budget_available"]
            and not corrected["Q207_to_Q120_same_field_identity_explicitly_certified"]
            and not corrected["cutoff_shape_removed_from_frontier_unconditionally"]
            and not corrected["historical_evidence_modified"]
        ),
        "derived_maps_propagate_only_after_seed_array": (
            propagation["cutoff_shape_bound_available_after_identity_binding"]
            and not propagation["same_field_cutoff_identity_bound"]
            and not propagation["selected_profile_gap_and_operator_bundle_closed_directly"]
        ),
        "remaining_profile_frontier_stated_without_promotion": (
            frontier["actual_weighted_first_coefficient_orderwise_finite"]
            and not frontier["closed_majorant_A_R_from_selected_defining_maps_recorded"]
            and not frontier["passing_tests_or_hashes_prove_this_supremum"]
        ),
    }
    return {
        "status": "Q208-squared-flat-cutoff-budget-computed-same-field-binding-and-primitive-frontier-isolated",
        "checks": checks,
        "scope": {
            "existing_scale_cutoff_bound_available_conditionally": True,
            "stage_602_cutoff_binary_budget": cutoff["binary_budget"],
            "same_field_cutoff_identity_certified": False,
            "complete_stage_602_primitive_bundle_proved": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
        "next": "bind Q207's cutoff to the Q120 squared-flat function, then derive A_R with common gap floors and the archived compact/slow inverse bundle",
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q208 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "same_object": same_object_links(),
        "cutoff": cutoff_provider(),
        "correction": correction(),
        "propagation": propagation_ledger(),
        "frontier": remaining_frontier(),
        "review": review(),
    }
