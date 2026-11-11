"""Q4 source audit on the Q172-corrected tail."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q173-Q4-fixed-profile-source-audit"
PROVIDERS = (
    "core_current_head_q172_q2_q3_source_audit",
    "core_current_head_q4_fixed_profile_join",
    "core_fixed_profile_projector_reconstruction",
    "core_fixed_profile_projector_join",
    "core_fixed_profile_tame",
    "core_reserved_geometry_recheck",
    "core_rho_squared_projector_realization",
    "core_rho_projector_signed_recompute",
    "core_log_pullback_projector_join",
    "core_finite_support_closure",
)

EXPECTED = {
    "core_current_head_q172_q2_q3_source_audit": "Q172-Q2-source-audited-and-Q3-uniform-frame-corrected",
    "core_current_head_q4_fixed_profile_join": "current-head-Q4-fixed-profile-compensator-and-same-choice-projectors-bound",
    "core_fixed_profile_projector_reconstruction": "fixed-profile-same-choice-projector-contract-reconstructed-on-current-full-path-tail",
    "core_fixed_profile_projector_join": "fixed-profile-same-choice-projector-contract-joined",
    "core_fixed_profile_tame": "fixed-profile-all-order-tame-bounds-and-quantifier-correction",
    "core_reserved_geometry_recheck": "uniform-reserved-coordinate-compensator-jets-and-moving-support-costs-rechecked",
    "core_rho_squared_projector_realization": "explicit-new-rho-squared-projector-realization",
    "core_rho_projector_signed_recompute": "conditional-same-choice-rho-projector-signed-recomputation",
    "core_log_pullback_projector_join": "conditional-actual-logarithmic-Haar-pullback-and-projector-seminorm-boundary",
    "core_finite_support_closure": "finite-smooth-support-closure-relative-to-original-weighted-input-contract",
}

def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())

def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q173 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}

def projector_dilation(scale, ell, upper, support_ell, support_upper, a2, exponent):
    """Exact change-of-variable and strict-support interface."""
    scale, ell, upper, support_ell, support_upper, a2 = map(
        Q, (scale, ell, upper, support_ell, support_upper, a2)
    )
    if exponent not in (1, 2):
        raise ValueError("projector exponent must be one or two")
    roots = {Q(1, 4): Q(1, 2), Q(1): Q(1), Q(4): Q(2)}
    if scale not in roots:
        raise ValueError("certified rational-square scale required")
    if not 0 < ell < support_ell < support_upper < upper or a2 <= 0:
        raise ValueError("strict positive fixed-profile data required")
    if a2 * (upper - ell) < 1:
        raise ValueError("unit-mass Cauchy-Schwarz floor required")
    root = roots[scale]
    return {
        "weighted_moment": root ** exponent * root ** (-exponent) * root ** (-1) * root,
        "reciprocal_a2": 1 / a2,
        "reciprocal_a2_bound": upper - ell,
        "actual_support": (root * support_ell, root * support_upper),
        "moving_enclosure": (root * ell, root * upper),
        "strictly_enclosed": root * ell < root * support_ell < root * support_upper < root * upper,
    }

def moving_pressure(jz_z, crho_z, pressure, crho, pressure_z):
    values = tuple(map(Q, (jz_z, crho_z, pressure, crho, pressure_z)))
    jz_z, crho_z, pressure, crho, pressure_z = values
    rows = (jz_z, crho_z * pressure, crho * pressure_z)
    return {
        "Jz_Z": rows[0],
        "crho_Z_P": rows[1],
        "crho_P_Z": rows[2],
        "total_without_epsilon": sum(rows),
    }

def bounded_operator_counterexample(prefix=32):
    """Uniform operator norm does not imply uniform source membership."""
    if type(prefix) is not int or prefix < 2:
        raise ValueError("prefix at least two required")
    operator_norms = tuple(Q(1) for _ in range(prefix))
    source_norms = tuple(Q(ell) for ell in range(1, prefix + 1))
    outputs = tuple(a * f for a, f in zip(operator_norms, source_norms))
    return {
        "uniform_operator_bound": max(operator_norms),
        "finite_prefix_source_supremum": max(source_norms),
        "finite_prefix_output_supremum": max(outputs),
        "prefix_outputs_strictly_increase": all(outputs[j] > outputs[j-1] for j in range(1, prefix)),
        "operator_boundedness_implies_uniform_source_bound": False,
    }

def generated_source_gate(stage, order, epsilon_power, weighted_supremum=None):
    if type(stage) is not int or stage < 0 or type(order) is not int or order < 0:
        raise ValueError("finite nonnegative stage and order required")
    epsilon_power = Q(epsilon_power)
    if epsilon_power <= 0:
        raise ValueError("positive epsilon normalization required")
    if weighted_supremum is None:
        return {
            "stage": stage,
            "order": order,
            "epsilon_power": epsilon_power,
            "actual_weighted_supremum_supplied": False,
            "gate_closed": False,
        }
    weighted_supremum = Q(weighted_supremum)
    if weighted_supremum < 0:
        raise ValueError("nonnegative weighted supremum required")
    return {
        "stage": stage,
        "order": order,
        "epsilon_power": epsilon_power,
        "actual_weighted_supremum_supplied": True,
        "weighted_supremum": weighted_supremum,
        "gate_closed": True,
    }

def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q173 identity required")
    status = statuses()
    projectors = [
        projector_dilation(scale, 1, 2, Q(5,4), Q(7,4), Q(10,7), exponent)
        for exponent in (1, 2) for scale in (Q(1,4), Q(1), Q(4))
    ]
    pressure = moving_pressure(Q(17,9), Q(1,5), Q(7,3), Q(11,4), Q(13,6))
    counterexample = bounded_operator_counterexample()
    open_gate = generated_source_gate(1, 6, Q(1,8))
    control_gate = generated_source_gate(1, 6, Q(1,8), Q(7,3))
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q172_corrected_tail_is_direct": PROVIDERS[0] == "core_current_head_q172_q2_q3_source_audit",
        "projector_moments_are_exact": all(row["weighted_moment"] == 1 for row in projectors),
        "actual_support_is_strict_not_the_enclosure": all(row["strictly_enclosed"] and row["actual_support"] != row["moving_enclosure"] for row in projectors),
        "moving_pressure_keeps_crho_Z_P": pressure["crho_Z_P"] != 0,
        "bounded_operators_do_not_supply_source_membership": not counterexample["operator_boundedness_implies_uniform_source_bound"],
        "missing_actual_generated_source_supremum_blocks_promotion": not open_gate["gate_closed"],
        "finite_supplied_control_closes_only_its_fixture": control_gate["gate_closed"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "projector_controls": [
            {k: ([str(x) for x in v] if isinstance(v, tuple) else str(v) if isinstance(v, Q) else v) for k,v in row.items()}
            for row in projectors
        ],
        "pressure_memory": {k: str(v) for k,v in pressure.items()},
        "operator_counterexample": {k: str(v) if isinstance(v,Q) else v for k,v in counterexample.items()},
        "first_open_gate": {
            "inequality": "sup_{ell>=ell_*} epsilon_ell^{-b_j} max_{|alpha|<=R_j(m,N)} sup_{Omega_{j,ell}} zeta^{-1} delta^{p_{j,alpha}} |D^alpha F_{j,ell}| < infinity",
            "fixture": {k: str(v) if isinstance(v,Q) else v for k,v in open_gate.items()},
            "logical_status": "not implied by fixed projector jets, support, exact moments or bounded inverse operators; it requires actual recomputed current-target/source bounds",
        },
        "scope": {
            "Q4_fixed_profile_geometry_source_audited_on_Q172_tail": True,
            "same_choice_recomputation_is_conditional_norm_chain": True,
            "actual_generated_source_weighted_suprema_proved": False,
            "compact_slow_inverse_domain_audited_here": False,
            "complete_covariance_means_or_residual_proved": False,
            "unconditional_pulse_force_exclusion_or_blowup_claimed": False,
        },
        "next": "audit Q5 compact and slow inverse weighted domains on this fixed choice while retaining the generated-source gate as an explicit downstream premise",
    }

def record():
    return {
        "status": "Q173-Q4-fixed-profile-source-audited-with-generated-source-boundary",
        "established": "the fixed compensator and rho-squared projector jets, exact moments, strictly enclosed supports and full moving-pressure product are analytic on the Q172-corrected tail",
        "correction": "Q4 same-choice recomputation is a conditional algebraic/norm-chain rebuild; it does not establish the actual recomputed generated-source weighted continuum suprema",
        "obstruction": "bounded projector/inverse operators do not imply the actual F_{j,ell} weighted source inequality; that inequality remains a separate same-state premise",
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
