"""Audit fixed-gauge pressure recovery and comparison on the Q168 object."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import sys


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q169-current-head-pressure-comparison-frontier-audit"
PROVIDERS = (
    "core_current_head_q168_localization_force_frontier_audit",
    "core_current_head_q156_physical_pressure_comparison",
    "core_current_head_q163_comparison_core_exclusion_join",
    "comparison_foundation",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q156 = load(
    "q169_q156",
    "core_current_head_q156_physical_pressure_comparison/q156.py",
)
q163 = load(
    "q169_q163",
    "core_current_head_q163_comparison_core_exclusion_join/q163.py",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q169 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q168-current-head-localization-force-energy-audited-conditionally",
        PROVIDERS[1]: "Q156-physical-fixed-gauge-pressure-and-expanding-ball-comparison-conditionally",
        PROVIDERS[2]: "Q163-same-object-comparison-and-core-exclusion-joined-conditionally",
        PROVIDERS[3]: "current-source-compact-pressure-comparison-relative-to-retained-core-and-finite-analytic-lemmas",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q169 provider status mismatch: " + name)

    paths = {
        "q167": "results/core_current_head_q167_corner_endpoint_frontier_audit/evidence.json",
        "q166": "results/core_current_head_q166_complete_residual_frontier_audit/evidence.json",
        "q165": "results/core_current_head_q165_finite_map_frontier_audit/evidence.json",
        "q164": "results/core_current_head_q164_core_fixed_point_pressure_audit/evidence.json",
        "q163": "results/core_current_head_q163_comparison_core_exclusion_join/evidence.json",
        "q162": "results/core_current_head_q162_diagonal_endpoint_force_join/evidence.json",
        "q156": "results/core_current_head_q156_physical_pressure_comparison/evidence.json",
        "q155": "results/core_current_head_q155_physical_force_energy_join/evidence.json",
        "comparison": "results/comparison_foundation/evidence.json",
    }
    q168, q156_record, q163_record, comparison = (data[name] for name in PROVIDERS)
    q167 = evidence("core_current_head_q167_corner_endpoint_frontier_audit")
    q166 = evidence("core_current_head_q166_complete_residual_frontier_audit")
    q165 = evidence("core_current_head_q165_finite_map_frontier_audit")
    q164 = evidence("core_current_head_q164_core_fixed_point_pressure_audit")
    chain = {
        "Q168_to_Q167": q168["dependency_evidence_sha256"][paths["q167"]] == sha(PROJECT / paths["q167"]),
        "Q167_to_Q166": q167["dependency_evidence_sha256"][paths["q166"]] == sha(PROJECT / paths["q166"]),
        "Q166_to_Q165": q166["dependency_evidence_sha256"][paths["q165"]] == sha(PROJECT / paths["q165"]),
        "Q165_to_Q164": q165["dependency_evidence_sha256"][paths["q164"]] == sha(PROJECT / paths["q164"]),
        "Q164_to_Q163": q164["dependency_evidence_sha256"][paths["q163"]] == sha(PROJECT / paths["q163"]),
        "Q163_to_Q156": q163_record["dependency_evidence_sha256"][paths["q156"]] == sha(PROJECT / paths["q156"]),
        "Q163_and_Q168_share_Q162": (
            q163_record["dependency_evidence_sha256"][paths["q162"]]
            == q168["dependency_evidence_sha256"][paths["q162"]]
        ),
        "Q156_and_Q168_share_Q155": (
            q156_record["dependency_evidence_sha256"][paths["q155"]]
            == q168["dependency_evidence_sha256"][paths["q155"]]
        ),
        "Q156_uses_frozen_comparison_foundation": (
            q156_record["dependency_evidence_sha256"][paths["comparison"]]
            == sha(PROJECT / paths["comparison"])
        ),
    }
    flags = comparison["flags"]
    return {
        "identity_chain": chain,
        "all_identity_links_hold": all(chain.values()),
        "same_Q162_field_force_pressure_gauge_and_localization": q163_record[
            "providers"
        ]["same_Q162_field_force_pressure_gauge_and_localization"],
        "Q168_localization_force_energy_closed_conditionally": q168["symbolic"][
            "U5_localization_force_energy_closed_conditionally"
        ],
        "comparison_pressure_gradient_recovery_provider": flags[
            "pressure_gradient_recovered_without_assumed_competitor_pressure_decay"
        ],
        "comparison_R_minus_two_provider": flags[
            "full_expanding_ball_energy_remainder_improved_from_four_thirds_to_two"
        ],
    }


def pressure_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q169 identity required")
    row = q156.pressure_ledger()
    return {
        "individual_pressures": row["individual"],
        "difference_equation": row["difference"],
        "stress": row["stress"],
        "force_cancellations": row["force_cancellations"],
        "constructed_pressure_gauge": row["constructed_gauge"],
        "riesz_pressure_space": row["riesz_pressure_space"],
        "averaged_gradient_space": row["averaged_gradient_space"],
        "recovery": row["recovery"],
        "recovery_logic": (
            "for every compact time test, momentum puts grad(pi-pi_star) in H^-3; "
            "it is harmonic, so its weighted-L2 Fourier transform is supported at {0} "
            "and therefore vanishes"
        ),
        "spatial_constant_gauge_flux": "zero by div(w)=0",
        "finite_energy_essential": row["finite_energy_essential"],
        "finite_energy_counterexample": (
            "without L2 control, v=t e1, P=-x1 and u=p=0 have the same zero datum "
            "and force but a nonzero harmonic pressure gradient"
        ),
        "competitor_pressure_decay_assumed": row["competitor_pressure_decay_assumed"],
        "Poisson_equation_alone_removes_harmonic_gradient": False,
    }


def commutator_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q169 identity required")
    row = q156.commutator_ledger(Q(4, 3), 8)
    weak = q156.commutator_ledger(Q(4, 3), 7)
    near = Q(1, 3 - 2 * Q(4, 3))
    far = Q(1, 3 * Q(4, 3) - 3)
    return {
        "kernel": row["kernel"],
        "exponent": row["p"],
        "near_radial_power": row["near_radial_power"],
        "far_radial_power": row["far_radial_power"],
        "near_coefficient": near,
        "far_coefficient": far,
        "radial_coefficient": row["radial_coefficient"],
        "norm_decay": row["norm_decay"],
        "cutoff_power": row["cutoff_power"],
        "cutoff_sufficient": row["cutoff_sufficient"],
        "seventh_power_sufficient": weak["cutoff_sufficient"],
        "weight_condition": "k/2-1>=3k/8",
    }


def radius_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q169 identity required")
    row = q156.radius_ledger()
    comparison = q156.comparison_ledger()
    mixed_decay = Q(3) + Q(1) - Q(3, 2)
    return {
        "cutoff": row["cutoff"],
        "support_split": row["support_split"],
        "mixed_pressure_point_decay": Q(3),
        "gradient_cutoff_decay": Q(1),
        "annular_L1_volume_cost": Q(3, 2),
        "compact_mixed_flux_decay": mixed_decay,
        "quadratic_flux": row["quadratic_flux"],
        "terms": row["terms"],
        "cutoff_laplacian_decay": row["other_decays"]["cutoff_laplacian"],
        "worst_decay": row["worst_decay"],
        "obsolete_unsplit_decay": row["obsolete_unsplit_decay"],
        "localized_inequality": comparison["localized"],
        "gronwall": comparison["gronwall"],
        "limit_order": comparison["limit_order"],
        "conclusion": comparison["conclusion"],
        "competitor_derivative_bound_assumed": comparison[
            "competitor_derivative_bound_assumed"
        ],
        "competitor_dissipation_bound_assumed": comparison[
            "competitor_dissipation_bound_assumed"
        ],
        "competitor_boundary_flux_bound_assumed": comparison[
            "competitor_boundary_flux_bound_assumed"
        ],
        "terminal_uniform_constants_assumed": comparison[
            "terminal_uniform_constants_assumed"
        ],
    }


def remaining_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q169 identity required")
    return {
        "U6_fixed_gauge_pressure_and_comparison": "closed-conditionally-by-Q169",
        "next": (
            "audit Q157/Q163's protected physical core path, localization plateau and legal "
            "fixed-T then radius then terminal limit order on the Q164-Q169 object"
        ),
        "comparison_only_requires_fixed_T_constants": True,
        "protected_core_or_terminal_exclusion_replayed_here": False,
        "future_velocity_or_pressure_constructed": False,
        "unconditional_theorem_gate_closed": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q169 identity required")
    return (
        {"id": "U1", "status": "closed-by-Q164-core-pressure-audit"},
        {"id": "U2", "status": "closed-by-Q165-finite-map-audit"},
        {"id": "U3", "status": "closed-conditionally-by-Q166-residual-audit"},
        {"id": "U4", "status": "closed-conditionally-by-Q167-corner-endpoint-audit"},
        {"id": "U5", "status": "closed-conditionally-by-Q168-localization-force-energy-audit"},
        {"id": "U6", "status": "closed-conditionally-by-Q169-pressure-comparison-audit"},
        {"id": "U7", "status": "next-Q157-Q163-protected-core-limit-audit"},
    )


def record():
    return {
        "status": "Q169-current-head-fixed-gauge-pressure-comparison-audited-conditionally",
        "correction": (
            "Q156 already lies in the Q164 ancestry, but Q168's forward frontier did not bind "
            "its comparison explicitly; Q169 verifies the same Q155/Q162 field and rederives "
            "the pressure and radius estimates without adding competitor hypotheses"
        ),
        "established": (
            "for every fixed T<1, any locally smooth same-force zero-datum competitor with a "
            "uniform L2 bound on [0,T] has the same pressure gradient modulo a spatial constant "
            "and equals the Q168 velocity after the expanding-ball radius tends to infinity"
        ),
        "scope": (
            "conditional analytical join of named continuum providers; no competitor pressure "
            "decay, derivative, dissipation or boundary-flux bound is assumed, constants need "
            "not be uniform as T approaches one, and no terminal exclusion is replayed here"
        ),
        "U6_pressure_comparison_closed_conditionally": True,
        "protected_core_exclusion_or_blowup_certified": False,
        "future_velocity_or_pressure_constructed": False,
        "historical_frozen_records_modified": False,
        "next": remaining_gate()["next"],
    }
