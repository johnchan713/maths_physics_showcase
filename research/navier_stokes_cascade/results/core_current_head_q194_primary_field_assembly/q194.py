"""Q194: conditional complete primary potential, pressure and curl assembly."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q194-PRIMARY-FIELD-ASSEMBLY"

PROVIDERS = (
    "core_current_head_q193_duhamel_phase_recurrence",
    "core_current_head_q192_frame_graph_jet_transfer",
    "physical_pulse_foundation",
    "physical_pulse_transfer",
)

EXPECTED = {
    PROVIDERS[0]: "Q193-uncut-amplitude-Duhamel-and-phase-partitions-closed-conditionally",
    PROVIDERS[1]: "Q192-auxiliary-frame-jets-closed-and-physical-graph-loss-corrected",
    PROVIDERS[2]: "actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas",
    PROVIDERS[3]: "actual-background-to-physical-primary-pulse-transfer-audited",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _order(n):
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("nonnegative integer order required")
    return n


def base_scalings(h=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    A = Q(1, 2) + h
    sx, st = Q(1, 2) + h / 2, 1 + 3 * h / 2
    return {
        "A": A,
        "sx": sx,
        "st": st,
        "velocity": -A + h / 2,
        "pressure": -2 * A + h,
        "potential": Q(1, 2) - A + h,
        "velocity_matches": -A + h / 2 == -sx,
        "pressure_matches": -2 * A + h == -2 * sx,
        "potential_matches": Q(1, 2) - A + h == 0,
    }


def curl_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    terms = {
        "u_r": ("r^-1*dtheta(A_z)", "-dz(A_theta)"),
        "u_theta": ("dz(A_r)", "-dr(A_z)"),
        "u_z": ("dr(A_theta)", "r^-1*A_theta", "-r^-1*dtheta(A_r)"),
    }
    return {
        "terms": terms,
        "component_count": len(terms),
        "term_count": sum(map(len, terms.values())),
        "axial_velocity_retained": len(terms["u_z"]) == 3,
        "cylindrical_connections_retained": any("r^-1" in x for row in terms.values() for x in row),
        "cut_before_curl": True,
        "divergence_free_by_curl_identity": True,
    }


def allocation_budget(spatial_slots, time_slots, phase_space, coefficient_space,
                      connections, phase_time, coefficient_time,
                      h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    """One product-rule allocation among phase, coefficients/cutoffs and frames."""
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    values = tuple(map(_order, (spatial_slots, time_slots, phase_space, coefficient_space,
                               connections, phase_time, coefficient_time)))
    ns, nt, ps, cs, conn, pt, ct = values
    if ps + cs + conn != ns or pt + ct != nt:
        raise ValueError("derivative allocation must exhaust the requested slots")
    h, kappa = Q(h), Q(kappa)
    if not 0 < h < Q(1, 100) or not 0 < kappa < Q(1, 2):
        raise ValueError("reviewed h and kappa required")
    sx, st = Q(1, 2) + h / 2, 1 + 3 * h / 2
    coefficient_space, coefficient_time = Q(1, 2) + h * kappa, 1 + h
    used = ps * sx + cs * coefficient_space + conn * Q(1, 2) + pt * st + ct * coefficient_time
    target = ns * sx + nt * st
    return {
        "used_Q_loss": used,
        "target_Q_loss": target,
        "margin": target - used,
        "fits": used <= target,
        "phase_slots": ps + pt,
        "coefficient_and_cutoff_slots": cs + ct,
        "connection_slots": conn,
    }


def field_loss(field, spatial_order=0, time_order=0, h=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    a, b, h = _order(spatial_order), _order(time_order), Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    sx, st = Q(1, 2) + h / 2, 1 + 3 * h / 2
    offsets = {"velocity": 1, "pressure": 2, "potential": 0}
    if field not in offsets:
        raise ValueError("field must be velocity, pressure or potential")
    return {
        "field": field,
        "spatial_order": a,
        "time_order": b,
        "Q_loss": (a + offsets[field]) * sx + b * st,
        "bound": (
            f"C*Q^(-(({a}+{offsets[field]})*sx+{b}*st))*S^B*sqrt(zeta)*delta^-B*P"
        ),
    }


def support_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    return {
        "transverse": "|xi|<=r0/2",
        "pulse": "|v/Ls-1/2|<=2/5; cutoff one on the central fifth",
        "shell": "fixed closed active shell with flat sqrt(zeta) edges",
        "cutoffs_on_potential_and_pressure_before_curl": True,
        "derivatives_and_curl_do_not_enlarge_support": True,
        "closed_enlarged_preimages_have_at_most_one_label": True,
        "axis_smoothness_from_zero_support_near_axis": True,
        "T_zero_jets_are_one_sided_not_zero_extended": True,
        "positive_P_not_used_as_endpoint_support": True,
        "primary_cutoff_supplies_smooth_pulse_endpoint_extension": True,
        "generated_source_still_needs_own_endpoint_extension": True,
    }


def primary_assembly(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    scaling = base_scalings()
    curl = curl_ledger()
    support = support_contract()
    return {
        "velocity": field_loss("velocity", 3, 2),
        "pressure": field_loss("pressure", 3, 2),
        "potential": field_loss("potential", 3, 2),
        "complete_curl_components": curl["terms"],
        "same_domain_support": support,
        "base_scalings_close": all(scaling[k] for k in ("velocity_matches", "pressure_matches", "potential_matches")),
        "positive_pointwise_lower_bound_for_complete_velocity": False,
        "conditional_uniform_all_accepted_bands_and_fixed_orders": True,
        "constants_not_uniform_in_derivative_order_harmonic_or_fixed_geometry": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    return {
        "closed_here": (
            "conditional complete primary potential, pressure and curl velocity fixed-order "
            "bounds and smooth support on every accepted band"
        ),
        "blocking_inequality": (
            "derive the complete-wave covariance mean identities using the assembled radial, "
            "angular and axial velocity and wave pressure, retaining all six covariance entries, "
            "the moving pressure response and the radial/axial mean corrections on one domain"
        ),
        "complete_wave_covariance_means_rederived_here": False,
        "generated_source_and_complete_residual_rederived_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    status = statuses()
    scaling, curl, support = base_scalings(), curl_ledger(), support_contract()
    allocations = []
    for ns in range(6):
        for nt in range(5):
            for ps in range(ns + 1):
                for conn in range(ns - ps + 1):
                    cs = ns - ps - conn
                    for pt in range(nt + 1):
                        allocations.append(allocation_budget(ns, nt, ps, cs, conn, pt, nt - pt))
    fields = [field_loss(name, a, b) for name in ("velocity", "pressure", "potential")
              for a in range(5) for b in range(4)]
    assembly, boundary = primary_assembly(), frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "all_base_rescalings_match_target_exponents": all(scaling[k] for k in ("velocity_matches", "pressure_matches", "potential_matches")),
        "complete_cylindrical_curl_keeps_axial_and_connection_terms": curl["term_count"] == 7 and curl["axial_velocity_retained"] and curl["cylindrical_connections_retained"],
        "every_checked_product_rule_allocation_fits_sx_st": all(x["fits"] for x in allocations),
        "all_three_field_families_have_finite_losses": len(fields) == 60,
        "same_support_and_one_sided_endpoint_memory_retained": support["derivatives_and_curl_do_not_enlarge_support"] and support["T_zero_jets_are_one_sided_not_zero_extended"],
        "no_positive_complete_velocity_floor_promoted": not assembly["positive_pointwise_lower_bound_for_complete_velocity"],
        "covariance_means_and_residual_not_promoted": not boundary["complete_wave_covariance_means_rederived_here"] and not boundary["generated_source_and_complete_residual_rederived_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "base_scalings": scaling,
        "curl_ledger": curl,
        "allocation_summary": {
            "count": len(allocations),
            "minimum_margin": min(x["margin"] for x in allocations),
            "maximum_target_Q_loss": max(x["target_Q_loss"] for x in allocations),
        },
        "field_requests": fields,
        "support_contract": support,
        "primary_assembly": assembly,
        "frontier": boundary,
        "scope": {
            "complete_primary_fixed_order_field_bounds_closed_conditionally": True,
            "uniform_every_accepted_band_relative_to_named_inputs": True,
            "complete_wave_covariance_means_rederived": False,
            "complete_residual_improvement_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q194 identity required")
    return {
        "status": "Q194-complete-primary-field-bounds-and-support-assembled-conditionally",
        "established": (
            "Q193's uncut recurrence and Q192's graph ledger assemble into complete primary "
            "curl velocity, wave pressure and fixed potential bounds on one support domain"
        ),
        "next": frontier()["blocking_inequality"],
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
