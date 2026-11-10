"""Q156 fixed-gauge pressure recovery and expanding-ball comparison ledger."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q156-physical-pressure-expanding-ball-comparison"
PROVIDERS = (
    "core_current_head_q155_physical_force_energy_join",
    "comparison_foundation",
    "core_current_head_q10_pressure_comparison_join",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def pressure_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q156 identity required")
    return {
        "constructed_gauge": "p=c p0 with Q155 complete fixed pressure memory",
        "individual": (
            "-Delta p=partial_i partial_j(u_i u_j)-div f; "
            "-Delta P=partial_i partial_j(v_i v_j)-div f"
        ),
        "difference": (
            "w_t+div g=Delta w-grad pi; -Delta pi=partial_i partial_j g_ij"
        ),
        "stress": "g=w tensor w+w tensor u+u tensor w",
        "recovery": (
            "grad pi=grad sum_ij R_i R_j g_ij in space-time distributions"
        ),
        "competitor_constant_gauge": (
            "a time-dependent spatial constant remains and has zero cutoff flux"
        ),
        "force_cancellations": 1,
        "stress_space": "L_infinity([0,T];L1)",
        "riesz_pressure_space": "L_infinity([0,T];H^-2)",
        "averaged_gradient_space": "H^-3",
        "finite_energy_essential": True,
        "competitor_pressure_decay_assumed": False,
    }


def commutator_ledger(p=Q(4, 3), cutoff_power=8, identity=IDENTITY):
    p = Q(p)
    if identity != IDENTITY or not 1 < p < Q(3, 2):
        raise ValueError("Q156 identity and 1<p<3/2 required")
    near_power = 2 - 2 * p
    far_power = 2 - 3 * p
    coefficient = Q(1, 3 - 2 * p) + Q(1, 3 * p - 3)
    norm_decay = 3 - 3 / p
    cutoff_sufficient = Q(cutoff_power, 2) - 1 >= Q(3 * cutoff_power, 8)
    return {
        "p": p,
        "kernel": "|z|^-3 min(|z|/R,1)",
        "near_radial_power": near_power,
        "far_radial_power": far_power,
        "radial_coefficient": coefficient,
        "norm_decay": norm_decay,
        "cutoff_power": cutoff_power,
        "cutoff_sufficient": cutoff_sufficient,
    }


def young_decay(dissipation_power, radius_power):
    a, b = Q(dissipation_power), Q(radius_power)
    if not 0 < a < 2 or b <= 0:
        raise ValueError("Require 0<a<2 and b>0")
    return 2 * b / (2 - a)


def radius_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q156 identity required")
    rows = (
        ("quadratic pressure and transport", Q(3, 2), Q(1)),
        ("quadratic pressure commutator", Q(3, 4), Q(7, 4)),
    )
    terms = [
        {
            "name": name,
            "dissipation_power": a,
            "radius_power": b,
            "young_decay": young_decay(a, b),
            "sobolev_offset_decay": a + b,
        }
        for name, a, b in rows
    ]
    other = {
        "compact_mixed_pressure": Q(5, 2),
        "cutoff_laplacian": Q(2),
    }
    decays = [item["young_decay"] for item in terms]
    decays += [item["sobolev_offset_decay"] for item in terms]
    decays += list(other.values())
    return {
        "cutoff": "chi_R=phi_R^8",
        "support_split": "g=w tensor w+h; h=w tensor u+u tensor w; supp h subset supp u",
        "compact_mixed_flux": "O(R^-5/2)",
        "quadratic_flux": "C_T(R^-1 B_R^(3/2)+R^-7/4 B_R^(3/4))",
        "terms": terms,
        "other_decays": other,
        "worst_decay": min(decays),
        "obsolete_unsplit_decay": young_decay(Q(1, 2), Q(1)),
    }


def comparison_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q156 identity required")
    return {
        "competitor": (
            "locally smooth v,P with the same f and zero datum and "
            "sup_[0,T] ||v||_2 finite for every T<1"
        ),
        "localized": (
            "0.5 E_R'+0.5 D_R^2<=||grad u||_infinity E_R+C_T R^-2"
        ),
        "gronwall": "E_R(t)<=2 C_T T exp(2 T L_T) R^-2",
        "limit_order": "fix T<1; recover pressure; send R to infinity",
        "conclusion": "v=u on R^3 x [0,T] for every T<1",
        "competitor_derivative_bound_assumed": False,
        "competitor_boundary_flux_bound_assumed": False,
        "competitor_dissipation_bound_assumed": False,
        "terminal_uniform_constants_assumed": False,
    }


def core_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q156 identity required")
    return {
        "required": (
            "there exist e0,A,tau0>0 such that x_tau lies in the Q155 cutoff "
            "plateau and |u_theta(x_tau,1-tau)|>=(e0/2) tau^-A for 0<tau<tau0"
        ),
        "Q155_frozen_interface_binds_quantitative_path": False,
        "global_smooth_competitor_exclusion_claimed": False,
        "blowup_claimed": False,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q156 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    if data[PROVIDERS[0]]["status"] != \
            "Q155-physical-localization-compact-force-and-energy-bound-conditionally":
        raise ValueError("Q155 physical force and energy parent required")
    if data[PROVIDERS[1]]["status"] != \
            "current-source-compact-pressure-comparison-relative-to-retained-core-and-finite-analytic-lemmas":
        raise ValueError("Comparison foundation required")
    if data[PROVIDERS[2]]["status"] != \
            "current-head-Q10-pressure-and-expanding-ball-comparison-joined-conditionally":
        raise ValueError("Frozen Q10 comparison replay required")
    flags = data[PROVIDERS[1]]["flags"]
    needed = (
        "common_force_divergence_cancels_exactly_once_in_difference",
        "compact_mixed_stress_pressure_estimated_outside_its_support",
        "full_expanding_ball_energy_remainder_improved_from_four_thirds_to_two",
        "general_comparison_lemma_has_explicit_ordinary_hypotheses",
        "pressure_gradient_recovered_without_assumed_competitor_pressure_decay",
    )
    if not all(flags.get(key) for key in needed):
        raise ValueError("Comparison foundation contracts missing")
    return {
        "Q155_same_field_force_energy": True,
        "refined_comparison_lemma": True,
        "historical_Q10_frozen": True,
        "quantitative_same_physical_core_path": False,
    }


def record():
    return {
        "status": "Q156-physical-fixed-gauge-pressure-and-expanding-ball-comparison-conditionally",
        "established": (
            "for the exact Q155 field, complete fixed-gauge pressure and force, every locally smooth "
            "same-force zero-datum competitor with fixed-interval finite L2 norm equals u for all t<1"
        ),
        "correction": (
            "the compact mixed-pressure split removes the older R^-4/3 remainder and leaves the "
            "coarse cutoff-Laplacian R^-2 remainder; pressure is recovered only up to a harmless "
            "time-dependent spatial constant"
        ),
        "obstruction": core_gate()["required"],
        "scope": (
            "conditional analytical pressure/comparison join; no competitor pressure, derivative, "
            "dissipation or boundary-flux bound is assumed, and no exclusion or blow-up is claimed"
        ),
        "next": (
            "bind a quantitative protected-core lower bound, with constants and path, to the exact "
            "Q153-Q155 physical field and localization plateau before attempting terminal exclusion"
        ),
        "preterminal_comparison_certified_conditionally": True,
        "global_competitor_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
