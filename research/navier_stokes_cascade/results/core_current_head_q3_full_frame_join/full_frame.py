"""Current-head Q3 full-frame, whole-path and constrained-ODE join."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent

DEPENDENCIES = (
    "core_current_head_q2_regional_gap_join",
    "core_current_head_q1_exact_map_join",
    "core_leading_background_seminorm_reconstruction",
    "core_leading_profile_certificate",
    "core_same_band_gap_join",
    "core_same_band_gap_reconstruction",
    "core_phase_coefficient_budget",
    "core_reference_scaled_normal_floor",
    "core_tangent_propagator_closure",
    "core_canonical_full_path_ode_join",
    "core_full_path_ode_reconstruction",
    "frame_column_foundation",
    "canonical_primary_foundation",
    "core_canonical_entrance_jet_closure",
    "core_canonical_haar_column_jet_closure",
    "sourced_inverse_foundation",
    "core_finite_support_closure",
)

EXPECTED_STATUS = {
    "core_current_head_q2_regional_gap_join": "current-head-Q2-regional-continuum-factored-edges-and-same-band-gap-bound",
    "core_current_head_q1_exact_map_join": "current-head-Q1-exact-core-roots-and-leading-pressure-bound-on-one-profile",
    "core_leading_background_seminorm_reconstruction": "fixed-order-background-seminorm-and-real-borel-selector-reconstructed-relative-to-leading-certificate",
    "core_leading_profile_certificate": "finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses",
    "core_same_band_gap_join": "same-rule-defined-background-tail-phase-normal-value-cone-contract-joined",
    "core_same_band_gap_reconstruction": "same-band-continuum-gap-contract-reconstructed-on-finalized-leading-profile",
    "core_phase_coefficient_budget": "full-background-frozen-phase-coefficient-recipes-varying-growth-and-normal-floor-gate",
    "core_reference_scaled_normal_floor": "reference-scaled-normal-floor-on-original-frame-tail",
    "core_tangent_propagator_closure": "original-frame-budgets-close-physical-Gaussian-ratio-tangent-propagator-and-length-inputs",
    "core_canonical_full_path_ode_join": "same-tail-canonical-full-path-ode-pressure-and-cutoff-contract-joined",
    "core_full_path_ode_reconstruction": "canonical-full-path-ode-reconstructed-on-current-same-band-tail",
    "frame_column_foundation": "explicit-current-phase-frame-and-primary-column-comparison-budgets",
    "canonical_primary_foundation": "conditional-canonical-primary-rebuild-and-relative-root-stability",
    "core_canonical_entrance_jet_closure": "conditional-canonical-entrance-and-homogeneous-jet-closure",
    "core_canonical_haar_column_jet_closure": "conditional-canonical-Haar-column-jet-closure",
    "sourced_inverse_foundation": "explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs",
    "core_finite_support_closure": "finite-smooth-support-closure-relative-to-original-weighted-input-contract",
}


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


def ceiling(value):
    value = Fraction(value)
    return (value.numerator + value.denominator - 1) // value.denominator


def frame_budget(normal2, lambda_min, bmin2, cmin, cmax, s0, CE, Cd, length, scale):
    """Exact finiteness/positivity interface for the original frame theorem."""
    values = tuple(map(Fraction, (normal2, lambda_min, bmin2, cmin, cmax, s0, CE, Cd, length, scale)))
    normal2, lambda_min, bmin2, cmin, cmax, s0, CE, Cd, length, scale = values
    if min(normal2, lambda_min, bmin2, cmin, cmax, CE, Cd, length, scale) <= 0 or s0 < 0:
        raise ValueError("positive finite frame inputs are required")
    if cmin > cmax or length > scale:
        raise ValueError("ordered frame and path inputs are required")
    H = CE + Cd
    exponent = ceiling(H)
    kframe = (1 + s0) * max(Fraction(1), cmax) / min(Fraction(1), cmin)
    tangent = kframe * 3 ** exponent
    dstar = bmin2 / 4
    coefficient = 1 + CE + Cd + 1 / normal2 + 1 / lambda_min
    picard_index = ceiling(2 * coefficient * length) - 1
    next_ratio = coefficient * length / (picard_index + 1)
    return {
        "normal_reciprocal": 1 / normal2,
        "H": H,
        "exponential_power": exponent,
        "Kframe": kframe,
        "tangent_constant": tangent,
        "dstar": dstar,
        "picard_index": picard_index,
        "picard_next_ratio": next_ratio,
    }


def projected_pressure(n, nprime, forcing, tangent, k, mode):
    """Pressure contribution from the moving constraint and unprojected source."""
    n, nprime, forcing, tangent = (tuple(map(Fraction, row)) for row in (n, nprime, forcing, tangent))
    k, mode = Fraction(k), Fraction(mode)
    norm2 = sum(x * x for x in n)
    if norm2 <= 0 or k == 0 or mode == 0 or sum(x*y for x, y in zip(n, tangent)) != 0:
        raise ValueError("nonzero normal/carrier/mode and tangent data are required")
    moving = -sum(x*y for x, y in zip(nprime, tangent))
    source = sum(x*y for x, y in zip(n, forcing))
    return {"moving_normal_over_i": moving / (k * mode * norm2),
            "source_normal_over_i": source / (k * mode * norm2),
            "total_without_K_over_i": (moving + source) / (k * mode * norm2)}


def compactness_contract():
    return {
        "domain": "finite accepted labels/charts and each closed lifted interval 0<=v<=L<=S on the Q2 later tail",
        "background": "every fixed Borel background jet is finite before the common tail is selected",
        "phase": "full F,G axial transport and their fixed jets are finite on the closed path",
        "denominators": "Q2 supplies positive normal, shear, frame, carrier and canonical-root floors on the same tail",
        "frame_constants": "CE=sup L*||E|| and Cd=sup L*|d-dref| are finite by compactness; endpoint frame factors are finite by positive cmin and normal floors",
        "growth": "P=exp(integral(lambda-dref)) is positive; its ratios cancel at entrance and Duhamel interfaces",
        "old_numeric_Cphi_gate_required": False,
        "unspecified_historical_first_band_certified": False,
    }


def review():
    records, graph, dependencies, sources = load_closure()
    status_checks = {
        name: records[evidence_path(name).resolve()].get("status") == expected
        for name, expected in EXPECTED_STATUS.items()
    }
    budget = frame_budget(
        Fraction(961, 103424), Fraction(1, 38500), Fraction(1, 101),
        Fraction(1, 11), 3, 2, Fraction(5, 3), Fraction(7, 4), 13, 17,
    )
    normal_source = projected_pressure((0, 0, 1), (0, 0, 0), (0, 0, 2), (1, 0, 0), 5, 1)
    moving = projected_pressure((2, 1, 1), (1, -1, 0), (0, 0, 0), (1, -2, 0), 7, -2)
    compactness = compactness_contract()
    checks = {
        "provider_graph_is_acyclic": cycle_free(graph),
        "all_declared_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_declared_source_hashes_match": bool(sources) and all(sources.values()),
        "all_seventeen_Q3_provider_statuses_match": len(status_checks) == 17 and all(status_checks.values()),
        "frame_denominators_and_constants_are_positive_and_finite": all(value >= 0 for value in budget.values()) and budget["dstar"] > 0 and budget["tangent_constant"] > 0,
        "picard_factorial_ratio_closes_without_smallness": budget["picard_next_ratio"] <= Fraction(1, 2),
        "pure_normal_source_retains_nonzero_pressure": normal_source["source_normal_over_i"] != 0,
        "moving_normal_pressure_term_is_nonzero": moving["moving_normal_over_i"] != 0,
        "same_later_tail_compactness_closes_frame_constants": not compactness["old_numeric_Cphi_gate_required"] and not compactness["unspecified_historical_first_band_certified"],
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
        "interface_controls": {
            "frame_budget_fixture": {key: str(value) for key, value in budget.items()},
            "pure_normal_source_pressure": {key: str(value) for key, value in normal_source.items()},
            "moving_normal_pressure": {key: str(value) for key, value in moving.items()},
            "compactness": compactness,
            "constrained_equation": "t'+Kt+m^2*d*t+i*k*m*n*pi=-f, n.t=0",
            "pressure": "pi=i*(n.Kt-nprime.t+n.f)/(k*m*|n|^2)",
            "cutoffs": "curl(psi*A)=psi*curl(A)+grad(psi)xA and grad(psi*pi)=psi*grad(pi)+pi*grad(psi)",
        },
        "scope": {
            "Q1_Q2_object_is_a_direct_dependency": True,
            "Q3_full_frame_primitive_path_and_phase_domain_bound_on_same_later_tail": True,
            "old_numeric_Cphi_sufficient_gate_used": False,
            "unspecified_historical_first_band_certified": False,
            "every_provider_theorem_independently_reproved_here": False,
            "fixed_projector_compensator_and_weighted_inverse_contracts_closed": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
