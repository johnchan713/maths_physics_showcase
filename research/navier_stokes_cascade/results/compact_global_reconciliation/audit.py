#!/usr/bin/env python3
"""Reconcile the explicit cutoff with compact repair and global pressure memory."""
import argparse
import hashlib
import json
from copy import deepcopy
from pathlib import Path

from join import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "b295e7c7179e638f4f2122b5477a8edbab7d7d03"
STATUS = "explicit-cutoff-to-same-root-compact-cone-and-global-pressure-reconciliation"
DEPENDENCIES = (
    "explicit_variance_cutoff",
    "mixed_shape_cone_bridge",
    "joint_fourth_loop_gamma",
    "fifth_loop_input_gate",
    "final_fifth_moment_interface",
    "regional_cone_join",
    "global_pressure_moment_foundation",
    "coefficient_support_foundation",
)

FLAGS = {
    "explicit_cutoff_discharges_mixed_third_and_joint_fourth_loop_gates": True,
    "fifth_gate_discharges_F5_and_complete_compact_fourth_map_budget": True,
    "same_N_same_root_repaired_state_and_inner_outer_value_cone_joined": True,
    "all_five_global_moment_rows_and_complete_pressure_memory_retained": True,
    "same_gauge_exterior_leading_stresses_joined": True,
    "finite_coefficient_support_and_reserved_mean_velocity_consequence_reached": True,
    "independent_core_native_regional_or_root_foundation_review": False,
    "global_FINAL_all_order_leading_profile_or_all_physical_scale_join": False,
    "uniform_all_family_pulses_complete_covariance_axial_pressure_means_residual_or_blowup": False,
}

REMAINING = [
    "discharge retained exact-core/native/regional/root premises rather than merely replaying their conditional estimates",
    "finish actual higher-coefficient and all-scale pulse-family defining/existence/inverse/gap/value/derivative/support estimates",
    "join complete-wave covariance means with radial and axial velocity and exceptional pressure tails",
    "prove the complete residual including viscosity, every cutoff and every interaction, then the distinct realization/comparison obligations",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protocol():
    return json.loads((HERE / "protocol.json").read_text())


def dependencies():
    records = {}
    statuses = {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        data = json.loads(path.read_text())
        records[str(path.relative_to(PROJECT))] = sha(path)
        statuses[name] = data["status"]
    return {"frozen_records": records, "statuses": statuses}


def sources():
    files = [HERE / n for n in ("README.md", "join.py", "audit.py", "test_audit.py", "protocol.json", "requirements.txt", ".gitignore")]
    files += [HERE.parent / n / "README.md" for n in DEPENDENCIES]
    files += [HERE.parent / n / "evidence.json" for n in DEPENDENCIES]
    return {str(p.relative_to(PROJECT)): sha(p) for p in files}


def dependency_claims():
    def load(name):
        return json.loads((HERE.parent / name / "evidence.json").read_text())
    e = load("explicit_variance_cutoff")["flags"]
    m = load("mixed_shape_cone_bridge")["flags"]
    j = load("joint_fourth_loop_gamma")["flags"]
    f = load("final_fifth_moment_interface")["flags"]
    r = load("regional_cone_join")["flags"]
    g = load("global_pressure_moment_foundation")["flags"]
    c = load("coefficient_support_foundation")["flags"]
    return {
        "selected_Z4_Z5_gates_closed": e["fourth_target_J_and_reduced_fifth_cutoff_gate_discharged_under_retained_foundations"],
        "same_N_exact_repair_bridge_available": m["same_N_exact_five_row_root_and_conditional_compact_value_cone_join_derived"],
        "joint_fourth_compact_field_interface_available": j["same_frequency_exact_repair_compact_FIELD_J4_less_than_K_squared_CONDITIONAL"],
        "complete_compact_pressure_flux_interface_available": f["complete_compact_joint_fourth_pressure_flux_and_ordinary_conversion_interface_derived"],
        "regional_same_frequency_value_join_available": r["same_finite_frequency_leading_value_cone_joined"],
        "global_five_moments_and_pressure_memory_available": g["all_five_global_leading_moment_identities_joined"] and g["full_pressure_memory_and_eta_derivative_retained"],
        "finite_coefficient_support_consequence_available": c["every_finite_coefficient_support_and_reserved_mean_velocity_joined"],
    }


def generate():
    rev = deepcopy(review())
    claims = dependency_claims()
    gates = dict(rev["checks"])
    gates.update(claims)
    return {
        "status": STATUS,
        "parent_commit": PARENT,
        "protocol": protocol(),
        "flags": dict(FLAGS),
        "reviews": rev,
        "gates": gates,
        "inherited_provenance": dependencies(),
        "source_hashes": sources(),
        "remaining_obligations": list(REMAINING),
    }


def validate(record):
    expected = generate()
    if json.dumps(record, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise ValueError("Compact/global reconciliation record changed")
    if not all(expected["gates"].values()):
        raise ValueError("A reconciliation gate failed")
    if any(type(v) is not bool for v in expected["flags"].values()):
        raise ValueError("Scope flags must be booleans")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use fresh output; frozen evidence is never overwritten")
    result = generate()
    validate(result)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"{STATUS}: {len(result['gates'])} checks; compact cone and global pressure stream reconciled relative to retained premises.")


if __name__ == "__main__":
    main()
