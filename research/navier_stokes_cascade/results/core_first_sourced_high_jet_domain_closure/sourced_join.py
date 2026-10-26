"""Compose the first-sourced compact graph with the certified first-real provider."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sourced = load(
    "conditional_first_sourced_tail",
    HERE.parent / "core_first_sourced_physical_tail_closure/closure.py",
)
current = load(
    "certified_first_real_provider",
    HERE.parent / "core_first_real_high_jet_domain_closure/high_jet_join.py",
)


PRIMITIVE_BINDINGS = {
    "analytic_background_and_Borel_coefficients": "core_leading_background_borel_join",
    "same_band_phase_normal_root_and_value_cone_gaps": "core_same_band_gap_join",
    "canonical_forced_ODE_and_full_path_cutoffs": "first_sourced_foundation",
    "fixed_rho_projectors_and_complete_pressure_maps": "core_fixed_profile_projector_join",
    "fixed_compact_temporal_and_reserved_slow_inverse_maps": "core_compact_slow_inverse_domain_join",
    "whole_support_endpoint_and_seam_extensions": "core_finite_support_closure",
}


def primitive_join():
    if set(sourced.PRIMITIVES) != set(PRIMITIVE_BINDINGS):
        raise ValueError("First-sourced primitive names changed")
    return {
        name: {
            "provider": PRIMITIVE_BINDINGS[name],
            "joined_on_current_certified_tail": True,
        }
        for name in sourced.PRIMITIVES
    }


def _retained_kappa(kappa):
    k = Q(kappa)
    if k != Q(1, 100000):
        raise ValueError("This join freezes the retained kappa=1/100000 graph")
    return k


def normalized_join(n=0, kappa=Q(1, 100000)):
    k = _retained_kappa(kappa)
    request = sourced.normalized_request(n, k)
    interface = request["largest_first_real_ceiling"]
    provider = current.normalized_join(interface, k)
    B, D, C = Q(1) - k, Q(3, 2) - 2 * k, Q(3, 2) - k
    achieved_A0 = B + request["integrations"]["A0"] * k
    achieved_A1 = D + request["integrations"]["A1"] * k
    return {
        "output_order": n,
        "kappa": k,
        "A0_source_order": request["raw_sources"]["A0"],
        "A1_source_order": request["raw_sources"]["A1"],
        "first_real_interface_ceiling": interface,
        "fully_expanded_primitive_ceiling": provider["combined_graph_ceiling"],
        "A0_source_rate": B,
        "A1_source_rate": D,
        "target_rate": C,
        "A0_tail_rate": achieved_A0,
        "A1_tail_rate": achieved_A1,
        "A0_margin": achieved_A0 - C,
        "A1_margin": achieved_A1 - C,
        "radial_inequality": (
            f"N_{request['raw_sources']['A0']}(Qcirc gr_S) "
            f"<= epsilon^(1-kappa) P_{request['raw_sources']['A0']}(S+K_all)"
        ),
        "axial_inequality": (
            f"N_{request['raw_sources']['A1']}(Qcirc Ez_postwave) "
            f"<= epsilon^(3/2-2*kappa) Q_{request['raw_sources']['A1']}(S+K_all)"
        ),
        "finite_and_acyclic": (
            sourced.fixed_order_induction(n=n)["requested_ceiling_is_finite"]
            and provider["all_requests_finite"]
            and not sourced.fixed_order_induction(n=n)["desired_A0_output_used_earlier"]
            and not provider["desired_compact_outputs_used_as_inputs"]
        ),
        "numerical_materialization_required": False,
    }


def physical_join(m=0, N=0, h=Q(1, 200), kappa=Q(1, 100000)):
    k = _retained_kappa(kappa)
    request = sourced.physical_request(m, N, h, k)
    interface = request["largest_first_real_ceiling"]
    provider = current.normalized_join(interface, k)
    return {
        "physical_jet": request["output_order"],
        "requested_power": request["requested_q_power"],
        "A0_source_order": request["raw_sources"]["A0"],
        "A1_source_order": request["raw_sources"]["A1"],
        "first_real_interface_ceiling": interface,
        "fully_expanded_primitive_ceiling": provider["combined_graph_ceiling"],
        "A0_margin": request["A0_physical_exponent"] - request["requested_q_power"],
        "A1_margin": request["A1_physical_exponent"] - request["requested_q_power"],
        "statement": f"|Dphysical^{m} Tcompact_S| <= C_{{{m},{N}}} q^{N}",
        "finite_graph": isinstance(provider["combined_graph_ceiling"], int),
        "constant_may_depend_on_fixed_request_and_profiles": True,
        "uniform_in_all_orders_stages_or_scales": False,
    }


def residual_obstruction(m=0, N=0, h=Q(1, 200), kappa=Q(1, 100000), eta=1):
    k = _retained_kappa(kappa)
    row = sourced.residual_interface(m, N, h, k, eta)
    return {
        "regular_order": row["regular_order"],
        "regular_physical_exponent": row["regular_physical_exponent"],
        "current_regular_decay_certified": row["current_regular_decay_certified"],
        "sufficient_later_order_threshold": row["sufficient_later_order_threshold"],
        "required_inequality": "h*a > N + Lambda + 2*m + eta",
        "complete_residual_decay_proved": row["complete_residual_decay_proved"],
    }


def review():
    primitives = primitive_join()
    row = normalized_join()
    physical = physical_join()
    obstruction = residual_obstruction()
    k = row["kappa"]
    checks = {
        "six_first_sourced_primitive_contracts_are_bound": len(primitives) == 6,
        "historical_interface_counts_are_retained": (
            row["first_real_interface_ceiling"] == 300021
            and physical["first_real_interface_ceiling"] == 60400023
        ),
        "interface_is_composed_with_current_provider": (
            row["fully_expanded_primitive_ceiling"] == 700030
            and physical["fully_expanded_primitive_ceiling"] == 60800032
        ),
        "normalized_raw_sources_are_exact": (
            row["A0_source_order"] == 100005 and row["A1_source_order"] == 7
        ),
        "both_normalized_compact_tails_have_strict_kappa_margin": (
            row["A0_tail_rate"] == Q(3, 2)
            and row["A1_tail_rate"] == Q(3, 2)
            and row["A0_margin"] == k
            and row["A1_margin"] == k
        ),
        "composed_graph_is_finite_and_acyclic": row["finite_and_acyclic"],
        "physical_raw_sources_are_exact": (
            physical["A0_source_order"] == 60200007
            and physical["A1_source_order"] == 60100009
        ),
        "physical_fixed_request_margins_are_strict": (
            physical["finite_graph"]
            and physical["A0_margin"] == Q(1, 20000000)
            and physical["A1_margin"] == Q(1, 20000000)
        ),
        "physical_quantifier_is_not_uniform": not physical["uniform_in_all_orders_stages_or_scales"],
        "regular_physical_exponent_remains_negative": (
            obstruction["regular_physical_exponent"] == Q(-15025001, 10000000)
            and not obstruction["current_regular_decay_certified"]
            and not obstruction["complete_residual_decay_proved"]
        ),
    }
    return {
        "checks": checks,
        "primitive_join": primitives,
        "normalized_join": {key: str(value) for key, value in row.items()},
        "physical_example": {key: str(value) for key, value in physical.items()},
        "residual_obstruction": {key: str(value) for key, value in obstruction.items()},
        "scope": (
            "analytical fixed-request composition on the certified later tail; "
            "no numerical giant-jet evaluation or complete residual decay"
        ),
    }
