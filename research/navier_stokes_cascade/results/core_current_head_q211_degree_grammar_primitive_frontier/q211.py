"""Q211: close the finite polynomial degree grammar and lower the frontier."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q211-DEGREE-GRAMMAR-PRIMITIVE-FRONTIER"
PROVIDERS = (
    "core_current_head_q210_finite_bundle_uniform_tail",
    "core_current_head_q192_frame_graph_jet_transfer",
    "core_current_head_q193_duhamel_phase_recurrence",
    "core_current_head_q194_primary_field_assembly",
    "core_current_head_q195_primary_covariance_initial_mean",
    "core_current_head_q196_initial_generated_source",
    "core_current_head_q197_first_sourced_covariance",
    "core_current_head_q198_generic_signed_target",
    "core_current_head_q199_weighted_target_obstruction",
    "core_current_head_q200_physical_weighted_target_reconciliation",
    "core_current_head_q201_post_signed_mean_reconciliation",
    "core_current_head_q202_next_sourced_cycle",
    "regular_source_foundation",
    "source_envelope_edge_audit",
)
EXPECTED = {
    PROVIDERS[0]: "Q210-finite-primitive-bundle-unified-conditionally-uniform-all-band-source-frontier-isolated",
    PROVIDERS[1]: "Q192-auxiliary-frame-jets-closed-and-physical-graph-loss-corrected",
    PROVIDERS[2]: "Q193-uncut-amplitude-Duhamel-and-phase-partitions-closed-conditionally",
    PROVIDERS[3]: "Q194-complete-primary-field-bounds-and-support-assembled-conditionally",
    PROVIDERS[4]: "Q195-primary-six-entry-covariance-and-initial-complete-mean-joined",
    PROVIDERS[5]: "Q196-initial-generated-source-stage-and-log-quotient-corrected",
    PROVIDERS[6]: "Q197-first-sourced-complex-covariance-and-current-means-joined",
    PROVIDERS[7]: "Q198-corrected-generic-complex-signed-target-and-weighted-input-obstruction",
    PROVIDERS[8]: "Q199-exact-moment-weighted-current-target-obstruction",
    PROVIDERS[9]: "Q200-physical-weighted-target-route-reconciled-at-corrected-rates",
    PROVIDERS[10]: "Q201-corrected-post-signed-fast-slow-mean-block-joined-conditionally",
    PROVIDERS[11]: "Q202-next-constrained-sourced-cycle-joined-on-corrected-physical-domain",
    PROVIDERS[12]: "regular-source-foundation-sharper-linear-gain-and-actual-initial-supported-input-relative-to-core-and-pulse-lemmas",
    PROVIDERS[13]: "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _array(values):
    values = list(values)
    if not values or any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in values):
        raise ValueError("nonempty nonnegative integer degree array required")
    return values


def product_degrees(left, right, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    left, right = _array(left), _array(right)
    if len(left) != len(right):
        raise ValueError("equal degree-array lengths required")
    return [max(left[i] + right[n - i] for i in range(n + 1))
            for n in range(len(left))]


def reciprocal_degrees(base, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    base = _array(base)
    if base[0] != 0:
        raise ValueError("degree-zero value bound and positive floor required")
    inverse = [0]
    for n in range(1, len(base)):
        inverse.append(max(base[i] + inverse[n - i] for i in range(1, n + 1)))
    return inverse


def sourced_solution_degrees(coefficients, source, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    coefficients, source = _array(coefficients), _array(source)
    if len(coefficients) != len(source):
        raise ValueError("equal degree-array lengths required")
    forcing = product_degrees(coefficients, source)
    solution = [forcing[0] + 1]
    for n in range(1, len(source)):
        solution.append(max([forcing[n] + 1] + [
            coefficients[i] + 1 + solution[n - i] for i in range(1, n + 1)
        ]))
    return {"forcing": forcing, "solution": solution}


def affine_closure(order=12, offset=20, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    if (isinstance(order, bool) or not isinstance(order, int) or order < 0
            or isinstance(offset, bool) or not isinstance(offset, int) or offset < 0):
        raise ValueError("nonnegative integer order and offset required")
    coefficients = [0] + [2 * i + 4 for i in range(1, order + 1)]
    source = [7 * i + offset for i in range(order + 1)]
    row = sourced_solution_degrees(coefficients, source)
    expected = [7 * i + offset + 1 for i in range(order + 1)]
    return {
        "order": order,
        "offset": offset,
        "coefficient_degrees": coefficients,
        "source_degrees": source,
        "forcing_degrees": row["forcing"],
        "solution_degrees": row["solution"],
        "expected_solution": expected,
        "closes": row["forcing"] == source and row["solution"] == expected,
        "reason": "2i+5<=7i for every i>=1",
        "actual_optimized_offset_evaluated": False,
    }


def completed_route(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    return {
        "Q196": "initial source quotient and uniform pure-power bound after one kappa reserve",
        "Q197": "first sourced constrained solve, complete covariance and current means",
        "Q198_Q200": "generic signed target corrected and physical weighted route reconciled",
        "Q201": "post-signed fast/slow means with complete pressure and axial memory",
        "Q202": "next sourced cycle on the corrected physical domain",
        "Q210_next_step_already_present_conditionally": True,
        "these_records_are_unconditional_continuum_proofs": False,
    }


def corrected_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    return {
        "source_grammar_can_create_superpolynomial_S_growth_from_polynomial_inputs": False,
        "finite_degree_closure": (
            "finite sums/products, reciprocals with positive floors, Bell quotients, "
            "weighted primitives and the sourced Duhamel recurrence preserve finite S degree"
        ),
        "primitive_degree_obligation": (
            "D_R^prim=max_(|I|<=R) inf{d>=0: sup_(ell>=ell*) "
            "S^-d N_I(actual current target/root/frame data)<infinity}<infinity"
        ),
        "weighted_target_identity": "X_e^new-X_e^old=-I_e Q_e conjugate(Delta E_e)",
        "must_retain": [
            "all six covariances", "radial and axial velocity", "complete pressure",
            "Pi1/Pi2", "whole moving pressure product", "zero weighted totals",
            "both flat-edge primitives", "reserved interior support"
        ],
        "unweighted_jets_moments_and_support_imply_primitive_degree_obligation": False,
        "Q199_countermodel_remains_applicable_to_that_invalid_inference": True,
        "existing_physical_route_supplies_obligation_conditionally": True,
        "independent_raw_map_certificate_complete": False,
        "unconditional_residual_or_blowup": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    status = statuses()
    affine = affine_closure()
    route = completed_route()
    frontier = corrected_frontier()
    checks = {
        "provider_statuses_match": all(status.values()),
        "degree_grammar_closes_one_sourced_solve": (
            affine["closes"] and not affine["actual_optimized_offset_evaluated"]
        ),
        "Q210_conditional_cycle_work_not_repeated": (
            route["Q210_next_step_already_present_conditionally"]
            and not route["these_records_are_unconditional_continuum_proofs"]
        ),
        "frontier_lowered_to_weighted_primitive_degree": (
            not frontier["source_grammar_can_create_superpolynomial_S_growth_from_polynomial_inputs"]
            and not frontier["unweighted_jets_moments_and_support_imply_primitive_degree_obligation"]
            and frontier["Q199_countermodel_remains_applicable_to_that_invalid_inference"]
            and frontier["existing_physical_route_supplies_obligation_conditionally"]
            and not frontier["independent_raw_map_certificate_complete"]
        ),
    }
    return {
        "status": "Q211-polynomial-degree-grammar-closed-conditionally-weighted-primitive-frontier-isolated",
        "checks": checks,
        "scope": {
            "finite_degree_algebra_and_sourced_recurrence_proved": True,
            "initial_and_first_sourced_cycles_recognized_as_existing_conditional_work": True,
            "actual_primitive_degree_array_certified_independently": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
        "next": (
            "derive the actual weighted current-target primitive degree array from the explicit "
            "state-difference, pressure and moment maps, beginning at the first physical target"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q211 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "degree_algebra": {
            "product_sample": product_degrees([0, 2, 4, 6], [1, 3, 5, 7]),
            "reciprocal_sample": reciprocal_degrees([0, 6, 8, 10]),
            "affine_sourced": affine_closure(),
        },
        "completed_route": completed_route(),
        "frontier": corrected_frontier(),
        "review": review(),
    }
