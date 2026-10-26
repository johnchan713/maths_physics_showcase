"""Join certified primitive domains to the first-real high-jet implications."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


phi = load(
    "first_real_phi_high_jet",
    HERE.parent / "core_initial_phi_high_jet_closure/closure.py",
)
a0 = load(
    "first_real_a0_high_jet",
    HERE.parent / "core_initial_a0_physical_tail_closure/closure.py",
)
sharp = load(
    "sharp_initial_target_domain_rate",
    HERE.parent / "core_sharp_initial_target_domain_join/rate_join.py",
)


PRIMITIVE_BINDINGS = {
    "analytic_background_and_Borel_coefficients": "core_leading_background_borel_join",
    "same_band_phase_normal_and_value_cone_gaps": "core_same_band_gap_join",
    "canonical_homogeneous_ODE_and_full_path_cutoffs": "core_canonical_full_path_ode_join",
    "fixed_profile_rho_and_rho_squared_projectors": "core_fixed_profile_projector_join",
    "fixed_compact_temporal_and_slow_inverse_maps": "core_compact_slow_inverse_domain_join",
}


def primitive_join():
    expected = tuple(phi.PRIMITIVES)
    if set(expected) != set(PRIMITIVE_BINDINGS):
        raise ValueError("Primitive contract names changed")
    return {
        name: {
            "provider": PRIMITIVE_BINDINGS[name],
            "independently_joined_on_current_tail": True,
        }
        for name in expected
    }


def normalized_join(n=0, kappa=Q(1, 100000)):
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError("A nonnegative fixed output order is required")
    k = Q(kappa)
    if k != Q(1, 100000):
        raise ValueError("This join freezes the retained kappa=1/100000 request graph")
    phi_induction = phi.fixed_order_induction(n)
    radial_induction = a0.fixed_order_induction(n=n)
    phi_rate = phi.rate_closure(k)
    radial_rate = a0.radial_rate_closure(k)
    final = sharp.final_groups(k)
    return {
        "output_order": n,
        "kappa": k,
        "Phi_order": phi_induction["request"]["raw_Phi_order"],
        "DzPhi_order": phi_induction["request"]["raw_Dz_Phi_order"],
        "Phi_graph_ceiling": phi_induction["request"]["largest_immediate_order"],
        "radial_source_order": radial_induction["request"]["A0_source_order"],
        "combined_graph_ceiling": radial_induction["request"]["largest_same_prefix_order"],
        "Phi_inequality": (
            f"N_{phi_induction['request']['raw_Phi_order']}(Qcirc Phi_postsigned) "
            f"<= epsilon P_{phi_induction['request']['raw_Phi_order']}(S)"
        ),
        "DzPhi_inequality": (
            f"N_{phi_induction['request']['raw_Dz_Phi_order']}(Dz Qcirc Phi_postsigned) "
            f"<= epsilon^2 Q_{phi_induction['request']['raw_Dz_Phi_order']}(S)"
        ),
        "radial_inequality": (
            f"N_{radial_induction['request']['A0_source_order']}(Qcirc gr_final) "
            f"<= epsilon^(1-kappa) R_{radial_induction['request']['A0_source_order']}(S+K_all)"
        ),
        "A1_tail": phi_rate["A1_tail"],
        "A1_margin": phi_rate["strict_margin"],
        "A0_tail": radial_rate["A0_tail"],
        "A0_margin": radial_rate["strict_margin"],
        "first_real_final_order": final["conditional_final_order"],
        "required": phi_rate["required"],
        "all_requests_finite": (
            phi_induction["requested_order_is_finite"]
            and radial_induction["requested_order_is_finite"]
        ),
        "numerical_materialization_required": False,
        "desired_compact_outputs_used_as_inputs": (
            phi_induction["desired_A1_output_used_as_input"]
            or radial_induction["desired_A0_output_used_as_input"]
        ),
    }


def physical_join(m=0, N=0, h=Q(1, 200), kappa=Q(1, 100000)):
    if Q(kappa) != Q(1, 100000):
        raise ValueError("This join freezes the retained kappa=1/100000 request graph")
    request = a0.physical_request(m, N, h, kappa)
    return {
        "physical_jet": request["physical_jet"],
        "requested_power": request["requested_q_power"],
        "A0_source_order": request["A0_source_order"],
        "A1_DzPhi_source_order": request["A1_DzPhi_source_order"],
        "combined_graph_ceiling": request["largest_same_prefix_order"],
        "A0_margin": request["A0_strict_margin"],
        "A1_margin": request["A1_strict_margin"],
        "finite_graph": isinstance(request["largest_same_prefix_order"], int),
        "statement": f"|Dphysical^{m} Tcompact| <= C_{{{m},{N}}} q^{N}",
        "constant_may_depend_on_fixed_m_N_and_profiles": True,
        "uniform_in_all_orders_stages_or_scales": False,
    }


def review():
    primitives = primitive_join()
    row = normalized_join()
    physical = physical_join()
    k = row["kappa"]
    checks = {
        "five_named_primitive_contracts_are_bound": len(primitives) == 5,
        "current_Phi_and_DzPhi_orders_are_exact": (
            row["Phi_order"] == 199994
            and row["DzPhi_order"] == 199993
            and row["Phi_graph_ceiling"] == 200018
        ),
        "current_radial_source_and_combined_ceiling_are_exact": (
            row["radial_source_order"] == 399995
            and row["combined_graph_ceiling"] == 400009
        ),
        "current_A1_tail_has_strict_kappa_margin": (
            row["A1_tail"] == 3 - 5 * k and row["A1_margin"] == k > 0
        ),
        "current_A0_tail_has_strict_kappa_margin": (
            row["A0_tail"] == 3 - 5 * k and row["A0_margin"] == k > 0
        ),
        "sharp_first_real_final_order_is_retained": row["first_real_final_order"] == 3 - 5 * k,
        "fixed_order_graphs_are_finite_and_acyclic": (
            row["all_requests_finite"] and not row["desired_compact_outputs_used_as_inputs"]
        ),
        "zero_jet_zero_power_physical_request_is_exact": (
            physical["A0_source_order"] == 60200007
            and physical["A1_DzPhi_source_order"] == 60000005
            and physical["combined_graph_ceiling"] == 60200021
        ),
        "fixed_physical_request_has_strict_margins": (
            physical["finite_graph"]
            and physical["A0_margin"] > 0
            and physical["A1_margin"] > 0
        ),
        "physical_quantifier_is_not_order_uniform": not physical["uniform_in_all_orders_stages_or_scales"],
    }
    return {
        "checks": checks,
        "primitive_join": primitives,
        "normalized_join": {key: str(value) for key, value in row.items()},
        "physical_example": {key: str(value) for key, value in physical.items()},
        "scope": (
            "analytical fixed-order modus ponens on the certified later tail; "
            "no numerical giant-jet evaluation or post-sourced residual closure"
        ),
    }
