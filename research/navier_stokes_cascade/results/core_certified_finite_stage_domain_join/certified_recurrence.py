"""Join the finite-stage compact recurrence to certified primitive domains."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


legacy = load(
    "conditional_finite_stage_recurrence",
    HERE.parent / "core_finite_stage_compact_recurrence/recurrence.py",
)
first_sourced = load(
    "certified_first_sourced_provider",
    HERE.parent / "core_first_sourced_high_jet_domain_closure/sourced_join.py",
)
first_real = load(
    "certified_first_real_provider_for_stages",
    HERE.parent / "core_first_real_high_jet_domain_closure/high_jet_join.py",
)


def _retained_kappa(kappa):
    k = Q(kappa)
    if k != Q(1, 100000):
        raise ValueError("This join freezes the retained kappa=1/100000 recurrence")
    return k


def primitive_contracts():
    contracts = first_sourced.primitive_join()
    if len(contracts) != 6:
        raise ValueError("First-sourced primitive contract set changed")
    return contracts


def expand_first_real_interface(interface_order, kappa=Q(1, 100000)):
    k = _retained_kappa(kappa)
    provider = first_real.normalized_join(interface_order, k)
    return {
        "first_real_interface_ceiling": interface_order,
        "fully_expanded_primitive_ceiling": provider["combined_graph_ceiling"],
        "provider_offset": provider["combined_graph_ceiling"] - interface_order,
        "provider_requests_are_finite": provider["all_requests_finite"],
        "provider_is_acyclic": not provider["desired_compact_outputs_used_as_inputs"],
    }


def state_join(j, n=0, kappa=Q(1, 100000)):
    k = _retained_kappa(kappa)
    conditional = legacy.state_ceiling(j, n, k)
    expanded = expand_first_real_interface(conditional["first_real_ceiling"], k)
    return {
        "stage": conditional["stage"],
        "output_order": conditional["output_order"],
        "first_real_interface_ceiling": expanded["first_real_interface_ceiling"],
        "fully_expanded_primitive_ceiling": expanded["fully_expanded_primitive_ceiling"],
        "provider_offset": expanded["provider_offset"],
        "layers": conditional["layers"],
        "finite_and_acyclic": (
            conditional["stage_decreases_on_every_recursive_edge"]
            and expanded["provider_requests_are_finite"]
            and expanded["provider_is_acyclic"]
        ),
    }


def physical_join(j, m=0, N=0, h=Q(1, 200), kappa=Q(1, 100000)):
    k = _retained_kappa(kappa)
    conditional = legacy.physical_request(j, m, N, h, k)
    expanded = expand_first_real_interface(conditional["first_real_ceiling"], k)
    return {
        "stage": conditional["stage"],
        "physical_jet": conditional["physical_jet"],
        "requested_q_power": conditional["requested_q_power"],
        "raw_sources": conditional["raw_sources"],
        "strict_margins": conditional["strict_margins"],
        "first_real_interface_ceiling": expanded["first_real_interface_ceiling"],
        "fully_expanded_primitive_ceiling": expanded["fully_expanded_primitive_ceiling"],
        "provider_offset": expanded["provider_offset"],
        "finite_and_acyclic": expanded["provider_requests_are_finite"] and expanded["provider_is_acyclic"],
        "statement": f"|Dphysical^{m} Tcompact_[{j}]| <= C_{{{j},{m},{N}}} q^{N}",
        "constant_may_depend_on_stage_jet_power_and_profiles": True,
        "uniform_in_stage_order_or_scale": False,
    }


def pole_support_join(j, m=0, N=0, h=Q(1, 200),
                      kappa=Q(1, 100000), original_pole=1,
                      retained_wave_pole=30):
    """Archive operators at each distinct pole after primitive expansion."""
    k = _retained_kappa(kappa)
    p = Q(original_pole)
    if p < 0:
        raise ValueError("nonnegative original pole required")
    physical = physical_join(j, m, N, h, k)
    requests = [
        ("physical", physical["fully_expanded_primitive_ceiling"]),
    ] + [
        (f"state_{i}", state_join(i, 0, k)["fully_expanded_primitive_ceiling"])
        for i in range(j + 1)
    ]
    retained = Q(retained_wave_pole)
    by_request = []
    distinct = []
    for label, L in requests:
        q = max(retained, 6 * p, 6 * L + 9, 6 * L + 12 * p + 6)
        by_request.append({"request": label, "primitive_order": L, "common_pole": q})
        if q not in distinct:
            distinct.append(q)
        retained = q
    inherited = legacy.pole_support_ledger(j, m, N, h, k, p, retained_wave_pole)
    return {
        "stage": j,
        "requests": requests,
        "common_pole_by_request": by_request,
        "archived_distinct_poles": distinct,
        "operator_archive": [
            {"KT1": 6*q, "KA1": 6*q, "Kp": (6*q, 12*q),
             "Kmom": (6*q, 12*q), "oriented": 6*q}
            for q in distinct
        ],
        "common_argument": (
            "S plus every archived K_6q_12q and K_Ie_6q; "
            "no distinct-pole monotonicity"
        ),
        "harmonics_at_most": inherited["harmonics_at_most"],
        "same_closed_supports_preserved": inherited["same_closed_supports_preserved_if_inherited_support_contracts_hold"],
        "finite_for_fixed_stage": all(isinstance(L, int) for _, L in requests),
        "actual_constants_or_high_jets_materialized": False,
    }


def diagonal_obstruction(N=0, m=0, h=Q(1, 200),
                         kappa=Q(1, 100000), reserve=1):
    k = _retained_kappa(kappa)
    threshold = legacy.first_stage_for_power(N, m, h, k, reserve)
    margin = threshold["achieved_power"] - N
    return {
        "first_fixed_stage": threshold["stage"],
        "preceding_stage": threshold["preceding_stage"],
        "fixed_stage_power_margin": margin,
        "preceding_power_margin": threshold["preceding_power"] - N,
        "fixed_stage_big_O_available": margin >= 0,
        "q_dependent_diagonal_required_inequality": (
            "log C_{j(q),m} <= (h*B_{j(q)}-loss_m-N)*log(1/q)"
        ),
        "stage_constant_growth_bound_available": False,
        "one_field_all_targets_diagonal_proved": False,
    }


def review():
    contracts = primitive_contracts()
    s0, s1 = state_join(0), state_join(1)
    p0, p1 = physical_join(0), physical_join(1)
    pole = pole_support_join(2)
    obstruction = diagonal_obstruction()
    checks = {
        "six_certified_primitive_contracts_are_retained": (
            len(contracts) == 6
            and all(row["joined_on_current_certified_tail"] for row in contracts.values())
        ),
        "state_interface_is_not_mislabeled_primitive": (
            s0["first_real_interface_ceiling"] == 500038
            and s0["fully_expanded_primitive_ceiling"] == 900047
            and s1["first_real_interface_ceiling"] == 1000067
            and s1["fully_expanded_primitive_ceiling"] == 1400076
        ),
        "provider_offset_is_applied_at_every_stage": (
            s0["provider_offset"] == 400009
            and s1["provider_offset"] == 400009
        ),
        "physical_stage_requests_are_fully_expanded": (
            p0["first_real_interface_ceiling"] == 60500042
            and p0["fully_expanded_primitive_ceiling"] == 60900051
            and p1["first_real_interface_ceiling"] == 60900077
            and p1["fully_expanded_primitive_ceiling"] == 61300086
        ),
        "physical_compact_margins_remain_strict": (
            min(p0["strict_margins"].values()) == Q(1, 20000000)
            and min(p1["strict_margins"].values()) == Q(1, 20000000)
        ),
        "fixed_stage_graphs_are_finite_and_acyclic": (
            s0["finite_and_acyclic"] and s1["finite_and_acyclic"]
            and p0["finite_and_acyclic"] and p1["finite_and_acyclic"]
        ),
        "common_pole_uses_expanded_primitive_requests": (
            pole["requests"][0] == ("physical", 61900109)
            and pole["archived_distinct_poles"] == [Q(371400672)]
        ),
        "support_and_harmonic_ledger_remains_finite": (
            pole["finite_for_fixed_stage"]
            and pole["same_closed_supports_preserved"]
            and pole["harmonics_at_most"] == 16
            and not pole["actual_constants_or_high_jets_materialized"]
        ),
        "fixed_stage_residual_threshold_is_retained": (
            obstruction["first_fixed_stage"] == 1002
            and obstruction["fixed_stage_power_margin"] >= 0
            and obstruction["preceding_power_margin"] < 0
        ),
        "diagonal_constant_growth_obstruction_is_explicit": (
            not obstruction["stage_constant_growth_bound_available"]
            and not obstruction["one_field_all_targets_diagonal_proved"]
        ),
    }
    return {
        "checks": checks,
        "primitive_contracts": contracts,
        "state_zero": {key: str(value) for key, value in s0.items()},
        "state_one": {key: str(value) for key, value in s1.items()},
        "physical_stage_zero": {key: str(value) for key, value in p0.items()},
        "physical_stage_one": {key: str(value) for key, value in p1.items()},
        "stage_two_pole_support": {key: str(value) for key, value in pole.items()},
        "diagonal_obstruction": {key: str(value) for key, value in obstruction.items()},
        "scope": (
            "analytical fixed-stage composition on certified primitive domains; "
            "no numerical giant-jet evaluation, stage-uniform diagonal or blow-up"
        ),
    }
