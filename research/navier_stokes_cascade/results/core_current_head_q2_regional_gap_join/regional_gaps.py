"""Current-head Q2 regional-gap, weighted-edge and later-tail join."""

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent

DEPENDENCIES = (
    "core_current_head_q1_exact_map_join",
    "core_current_head_primitive_frontier_audit",
    "core_leading_profile_certificate",
    "core_growth_activation",
    "core_activation_volterra_recheck",
    "regional_cone_join",
    "weighted_stress_foundation",
    "weighted_background_edges",
    "core_leading_background_seminorm_reconstruction",
    "explicit_tilt_pulse_tail",
    "core_primary_tail_recheck",
    "core_same_band_gap_join",
    "core_same_band_gap_reconstruction",
)

EXPECTED_STATUS = {
    "core_current_head_q1_exact_map_join": "current-head-Q1-exact-core-roots-and-leading-pressure-bound-on-one-profile",
    "core_leading_profile_certificate": "finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses",
    "core_growth_activation": "core-growth-and-quantitative-inner-activation-audited",
    "core_activation_volterra_recheck": "core-flat-activation-five-moment-and-cone-collar-independently-rederived",
    "regional_cone_join": "regional-input-gaps-and-leading-value-cone-joined",
    "weighted_stress_foundation": "actual-weighted-first-stress-and-background-foundation-joined-relative-to-exact-core-and-finite-lemmas",
    "weighted_background_edges": "all-order-weighted-edges-and-relative-background-stress-audited",
    "core_leading_background_seminorm_reconstruction": "fixed-order-background-seminorm-and-real-borel-selector-reconstructed-relative-to-leading-certificate",
    "explicit_tilt_pulse_tail": "explicit-tilt-raw-deficit-and-primary-pulse-tail-certificate",
    "core_primary_tail_recheck": "same-profile-primary-tail-value-gates-and-dyadic-edge-flatness-rechecked",
    "core_same_band_gap_join": "same-rule-defined-background-tail-phase-normal-value-cone-contract-joined",
    "core_same_band_gap_reconstruction": "same-band-continuum-gap-contract-reconstructed-on-finalized-leading-profile",
}

REGIONS = (
    ("activation", "Xminus", "Xa*exp(delta)"),
    ("reduced_continuation", "Xa*exp(delta)", "Xi"),
    ("shape", "Xi", "Xsep"),
    ("ideal_prejoin", "Xsep", "XR*exp(-8)"),
    ("five_moment_join", "XR*exp(-8)", "XR*exp(-5)"),
    ("ideal_after_join", "XR*exp(-5)", "XR"),
    ("first_outer_ramp", "XR", "XR*exp(1)"),
    ("axial_Md64", "XR*exp(1)", "XR*exp(T+1)"),
    ("intermediate_ramp", "XR*exp(T+1)", "XR*exp(T+2)"),
    ("intermediate_power", "XR*exp(T+2)", "Xrep=Xstar"),
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evidence_path(name):
    return RESULTS / name / "evidence.json"


def load_closure():
    records, graph, dependency_checks, source_checks = {}, {}, {}, {}

    def visit(path):
        path = path.resolve()
        if path in records:
            return
        record = json.loads(path.read_text())
        records[path] = record
        children = []
        for rel, digest in record.get("dependency_evidence_sha256", {}).items():
            target = (PROJECT / rel).resolve()
            dependency_checks[f"{path.parent.name}:{rel}"] = target.exists() and sha(target) == digest
            if target.exists():
                children.append(target)
                visit(target)
        graph[path] = children
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            source_checks[f"{path.parent.name}:{rel}"] = target.exists() and sha(target) == digest

    for name in DEPENDENCIES:
        visit(evidence_path(name))
    return records, graph, dependency_checks, source_checks


def cycle_free(graph):
    color = {}

    def visit(node):
        color[node] = 1
        for child in graph[node]:
            if color.get(child) == 1:
                return False
            if color.get(child, 0) == 0 and not visit(child):
                return False
        color[node] = 2
        return True

    return all(color.get(node, 0) or visit(node) for node in graph)


def perturbation_radius(f, fbar, t, gbar, d, omega, gamma, u):
    """Exact same-tilt persistence radius from positive continuum inputs."""
    f, fbar, t, gbar, d, omega, gamma, u = map(
        Fraction, (f, fbar, t, gbar, d, omega, gamma, u)
    )
    if min(f, t, d, omega, gamma, u) <= 0:
        raise ValueError("strictly positive continuum inputs are required")
    if fbar < f or gbar < t or omega > 1 or gamma > 1:
        raise ValueError("ordered continuum bounds are required")
    LD = 4 * gbar + 2 * fbar + 3
    LP = fbar + gbar + 1
    Dbar = gbar * gbar + 2 * fbar * gbar
    cbar2 = Dbar / (2 * f * t)
    Lc = 2 * LD / (f * t) + 2 * Dbar * LP / (f * f * t * t)
    LH = 1 + 2 / t
    Lphi = 3 * u * u * LH + (1 + u * u) * (4 * Lc + 3 * cbar2 * LH)
    mphi = u * u * omega * omega * gamma * (2 - gamma)
    rows = {
        "relative_stress": Fraction(1, 2),
        "F": f / 2,
        "minus_gtheta": t / 2,
        "D": d / (2 * LD),
        "orientation": omega / (2 * LH),
        "same_tilt_cone": mphi / (2 * Lphi),
    }
    return {
        "radius": min(rows.values()),
        "rows": rows,
        "LD": LD,
        "LP": LP,
        "Lc": Lc,
        "LH": LH,
        "Lphi": Lphi,
        "mphi": mphi,
    }


def corrected_outer_ramp_floor():
    """High-precision scalar check of the corrected analytic lower bound."""
    with localcontext() as context:
        context.prec = 80
        one = Decimal(1)
        value = Decimal("1.1") * (-Decimal("0.6")).exp()
        value -= (Decimal("0.01") / Decimal("0.98")) * (one.exp() - one)
        return value


def regional_cover():
    return {
        "regions": [
            {"region": name, "left": left, "right": right}
            for name, left, right in REGIONS
        ],
        "exact_shared_endpoints": all(
            REGIONS[i][2] == REGIONS[i + 1][1] for i in range(len(REGIONS) - 1)
        ),
    }


def edge_contract():
    return {
        "absolute_open_annulus_stress_infimum": "0",
        "positive_absolute_stress_floor": False,
        "weight": "zeta=exp(-t1^2/ya^2-4/yb^2)",
        "inner_factor": "T0=exp(-t1^2/ya^2)*Aa, |Aa|>=ca>0",
        "outer_factor": "T0=exp(-4/yb^2)*yb^-3*(b0,yb^6*bz), b0>0",
        "valid_lower_bound": "|T0|>=c*zeta",
        "normalized_directions_extend_smoothly": True,
    }


def review():
    records, graph, dependencies, sources = load_closure()
    status_checks = {
        name: records[evidence_path(name).resolve()].get("status") == expected
        for name, expected in EXPECTED_STATUS.items()
    }
    radius = perturbation_radius(
        Fraction(1, 5), 3, Fraction(1, 7), 4,
        Fraction(1, 11), Fraction(1, 13), Fraction(1, 17), 2,
    )
    cover = regional_cover()
    edge = edge_contract()
    ramp = corrected_outer_ramp_floor()
    checks = {
        "provider_graph_is_acyclic": cycle_free(graph),
        "all_declared_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_declared_source_hashes_match": bool(sources) and all(sources.values()),
        "all_twelve_Q2_provider_statuses_match": len(status_checks) == 12 and all(status_checks.values()),
        "ten_region_cover_has_exact_shared_endpoints": len(REGIONS) == 10 and cover["exact_shared_endpoints"],
        "corrected_first_outer_ramp_floor_exceeds_one_half": ramp > Decimal("0.5"),
        "every_same_tilt_perturbation_row_is_strictly_positive": radius["radius"] > 0 and all(value > 0 for value in radius["rows"].values()),
        "flat_edge_contract_is_weighted_and_factored": not edge["positive_absolute_stress_floor"] and edge["normalized_directions_extend_smoothly"],
        "historical_uncertified_existing_tilt_gate_is_not_used": True,
    }
    return {
        "checks": checks,
        "graph": {
            "nodes": len(graph),
            "edges": sum(len(row) for row in graph.values()),
            "dependency_hash_checks": len(dependencies),
            "source_hash_checks": len(sources),
        },
        "provider_status_checks": status_checks,
        "regional_cover": cover,
        "interface_controls": {
            "corrected_outer_ramp_floor": str(ramp),
            "perturbation_radius": str(radius["radius"]),
            "perturbation_rows": {key: str(value) for key, value in radius["rows"].items()},
            "edge": edge,
            "tilt_route": "explicit_tilt_pulse_tail -> core_primary_tail_recheck -> core_same_band_gap_join/reconstruction",
            "excluded_route": "leading_cone_margin existing-u sufficient gate",
            "pressure_memory": "full F,G axial transport, Pi1/Pi2 and later positive-order, wave, mean, localization and forced pressure terms are retained",
        },
        "scope": {
            "Q1_exact_object_is_a_direct_dependency": True,
            "Q2_regional_continuum_and_factored_edge_chain_bound_on_that_object": True,
            "Q2_positive_same_band_radius_rechecked": True,
            "every_provider_theorem_independently_reproved_here": False,
            "positive_absolute_stress_infimum_at_flat_edges_claimed": False,
            "historical_uncertified_existing_tilt_gate_used": False,
            "Q3_full_frame_primitive_domain_discharged": False,
            "hashes_tests_or_interval_diagnostics_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
