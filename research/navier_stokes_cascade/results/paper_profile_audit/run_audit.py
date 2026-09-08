#!/usr/bin/env python3
"""Fail-closed, source-bound audit of explicit building blocks, not a proof."""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy

from profiles import (ManufacturedProfile, coordinates, core_coefficients,
                      core_comparison, cylindrical_and_stress, field,
                      finite_difference_residual, full_residual, heat_derivatives,
                      heat_exterior, heat_reference, normalized_error, physical_point)

HERE = Path(__file__).resolve().parent
SOURCES = ("protocol.json", "profiles.py", "run_audit.py", "test_profiles.py")
STATUS = "building-block-audit-passed"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def specifications(protocol):
    """The expected check set and limits are independent of reported values."""
    algebra = protocol["maximum_algebraic_identity_error"]
    spec = {name: ("at-most", algebra) for name in (
        "coordinate-map", "coordinate-derivatives", "divergence", "axis-divergence",
        "centrifugal-balance", "cartesian-cylindrical", "integrated-stress")}
    spec["axis-finite"] = ("exact", 1.)
    for i, _ in enumerate(protocol["finite_difference"]["similarity_points"]):
        spec[f"fd-{i}-finest-error"] = ("at-most", protocol["finite_difference"]["maximum_normalized_error"])
        spec[f"fd-{i}-reduction"] = ("at-least", protocol["finite_difference"]["minimum_coarse_to_fine_error_reduction"])
        spec[f"fd-{i}-decreasing"] = ("exact", 1.)
    spec["heat-reference"] = ("at-most", protocol["heat"]["maximum_reference_error"])
    for name in ("heat-ode", "heat-full-residual", "heat-curvature-identity"):
        spec[name] = ("at-most", protocol["heat"]["maximum_identity_error"])
    spec["heat-positive-decreasing"] = ("exact", 1.)
    spec["core-reference"] = ("at-most", protocol["core_comparison"]["maximum_reference_error"])
    for name in ("core-exact-recurrence", "core-exact-positive-lower-bound"):
        spec[name] = ("exact", 1.)
    for name in protocol["required_negative_controls"]:
        spec["detect-" + name] = ("at-least", protocol["negative_control_minimum_normalized_gap"])
    return spec


def validate(report, protocol):
    """Reject missing/extra checks, altered limits, nonfinite values and promotion."""
    def finite_tree(value):
        if isinstance(value, dict):
            for child in value.values():
                finite_tree(child)
        elif isinstance(value, list):
            for child in value:
                finite_tree(child)
        elif isinstance(value, float) and not math.isfinite(value):
            raise ValueError("Nonfinite raw evidence")
    finite_tree(report)
    fixture_axes = [protocol["fixtures"][key] for key in ("q", "eta", "X", "theta")]
    if [r["similarity_point"] for r in report["manufactured_samples"]] != [list(p) for p in itertools.product(*fixture_axes)]:
        raise ValueError("Missing or changed manufactured samples")
    if report["axis_samples"] != len(protocol["fixtures"]["q"])*len(protocol["fixtures"]["eta"]):
        raise ValueError("Missing axis samples")
    if [(r["h"], r["Z"]) for r in report["heat_profile"]] != list(itertools.product(protocol["heat"]["h"], protocol["heat"]["Z"])):
        raise ValueError("Missing or changed heat profile samples")
    if [(r["h"], r["point"]) for r in report["heat_exterior"]] != list(itertools.product(protocol["heat"]["h"], protocol["heat"]["physical_points"])):
        raise ValueError("Missing or changed heat exterior samples")
    if [r["z"] for r in report["core_comparison"]] != protocol["core_comparison"]["z"]:
        raise ValueError("Missing or changed scalar comparison samples")
    if [r["similarity_point"] for r in report["finite_difference"]] != protocol["finite_difference"]["similarity_points"]:
        raise ValueError("Missing finite-difference samples")
    if any(r["relative_steps"] != protocol["finite_difference"]["relative_steps"] for r in report["finite_difference"]):
        raise ValueError("Changed finite-difference steps")
    checks, spec = report["checks"], specifications(protocol)
    if set(checks) != set(spec):
        raise ValueError("Missing or unexpected audit checks")
    failures = []
    for name, (kind, limit) in spec.items():
        check = checks[name]
        value = check["value"]
        if check["kind"] != kind or check["limit"] != limit:
            raise ValueError("Changed gate: " + name)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("Nonfinite or nonnumeric gate: " + name)
        passed = ((kind == "at-most" and 0 <= value <= limit) or
                  (kind == "at-least" and value >= limit) or
                  (kind == "exact" and value == limit))
        if not passed:
            failures.append(name)
    expected_status = "building-block-audit-failed" if failures else STATUS
    if (report["status"] != expected_status or report["full_proof_verified"] is not False or
            report["new_candidate"] is not False or
            report["prior_unforced_status"] != "unchanged-resolution-failure"):
        raise ValueError("Invalid scientific status or unsupported promotion")
    if report["unverified_obligations"] != protocol["unverified_obligations"]:
        raise ValueError("Unverified obligations changed or omitted")
    if failures:
        raise ValueError("Audit gates failed: " + ", ".join(failures))


def audit(protocol):
    if protocol["fixtures"]["viscosity"] != 1.0:
        raise ValueError("This audit implements viscosity one only")
    p = ManufacturedProfile(protocol["fixtures"])
    h = protocol["fixtures"]["h"]
    A, D = .5+h, .5-h
    spec = specifications(protocol)
    checks = {name: {"kind": kind, "limit": limit, "value": None}
              for name, (kind, limit) in spec.items()}
    maxima = {name: 0. for name in ("coordinate-map", "coordinate-derivatives", "divergence",
                                   "centrifugal-balance", "cartesian-cylindrical", "integrated-stress")}
    fixture_rows, axial_gaps, sign_gaps = [], [], []
    axes = [protocol["fixtures"][key] for key in ("q", "eta", "X", "theta")]
    for q0, eta0, X0, theta in itertools.product(*axes):
        point = physical_point(q0, eta0, X0, theta, h)
        u, pressure, (q, X, eta, _) = field(point, h, p)
        residual, scale, _ = full_residual(u, pressure)
        cylinder = cylindrical_and_stress(point, h, p)
        qmap = normalized_error(np.array([q.v, X.v, eta.v]), np.array([q0, X0, eta0]))
        d, L = 1-eta.v**2, 1-2*h*eta.v**2
        derivative_gap = 0.
        # q^b F is tested at positive and negative exponents; jets use implicit q.
        for b in (-2*A, -A, 1.):
            scalar = q**b*p.F(X, eta)
            f = p.F(X.v, eta.v)
            fx, fe = [p.F.derivative(i)(X.v, eta.v) for i in (0, 1)]
            time_formula = q.v**(b-1)/L*(-b*f + D*eta.v*fe + X.v*fx)
            z_formula = q.v**(b-D)/L*(2*b*eta.v*f + d*fe - 2*eta.v*X.v*fx)
            derivative_gap = max(derivative_gap, normalized_error(
                np.array([scalar.g[3], scalar.g[2]]), np.array([time_formula, z_formula])))
        divergence = abs(sum(u[i].g[i] for i in range(3)))/max(1., *(abs(u[i].g[i]) for i in range(3)))
        radius = np.linalg.norm(point[:2])
        er = point[:2]/radius
        ut = -er[1]*u[0].v+er[0]*u[1].v
        centrifugal = normalized_error(float(er @ pressure.g[:2]), float(ut*ut/radius))
        row = {"similarity_point": [q0, eta0, X0, theta], "physical_point": point.tolist(),
               "full_cartesian_residual": residual.tolist(), "term_scale": scale,
               "leading_tangential_residual": cylinder["leading_tangential"].tolist(),
               "axial_second_derivatives": cylinder["axial_diffusion"].tolist(),
               "full_residual_over_term_scale": float(np.max(np.abs(residual)))/scale}
        values = {"coordinate-map": qmap, "coordinate-derivatives": derivative_gap,
                  "divergence": divergence, "centrifugal-balance": centrifugal,
                  "cartesian-cylindrical": normalized_error(residual, cylinder["cartesian"]),
                  "integrated-stress": normalized_error(cylinder["leading_tangential"], cylinder["stress_rhs"])}
        for name, value in values.items():
            maxima[name] = max(maxima[name], value)
        row["identity_errors"] = values
        fixture_rows.append(row)
        axial_gaps.append(normalized_error(cylinder["cylindrical"][1:], cylinder["stress_rhs"]))
        sign_gaps.append(normalized_error(cylinder["leading_tangential"], -cylinder["stress_rhs"]))
    for name, value in maxima.items():
        checks[name]["value"] = value
    axis_gaps, axis_finite = [], True
    for q0, eta0 in itertools.product(protocol["fixtures"]["q"], protocol["fixtures"]["eta"]):
        u, pressure, _ = field(physical_point(q0, eta0, 0., 0., h), h, p)
        for value in [*u, pressure]:
            axis_finite &= bool(np.isfinite(value.v) and np.isfinite(value.g).all() and np.isfinite(value.H).all())
        axis_gaps.append(abs(sum(u[i].g[i] for i in range(3)))/max(1., *(abs(u[i].g[i]) for i in range(3))))
    checks["axis-divergence"]["value"] = max(axis_gaps)
    checks["axis-finite"]["value"] = float(axis_finite)
    fd_rows = []
    for i, similarity in enumerate(protocol["finite_difference"]["similarity_points"]):
        point = physical_point(*similarity, h)
        u, pressure, _ = field(point, h, p)
        residual, scale, _ = full_residual(u, pressure)
        errors = []
        for step in protocol["finite_difference"]["relative_steps"]:
            independent = finite_difference_residual(point, h, p, step)
            errors.append(float(np.max(np.abs(independent-residual)))/scale)
        reduction = errors[0]/max(errors[-1], np.finfo(float).tiny)
        checks[f"fd-{i}-finest-error"]["value"] = errors[-1]
        checks[f"fd-{i}-reduction"]["value"] = reduction
        checks[f"fd-{i}-decreasing"]["value"] = float(all(a > b for a, b in zip(errors, errors[1:])))
        fd_rows.append({"similarity_point": similarity, "relative_steps": protocol["finite_difference"]["relative_steps"],
                        "normalized_errors": errors, "coarse_to_fine_reduction": reduction})
    heat_rows, exterior_rows = [], []
    reference_gaps, ode_gaps, full_gaps, curvature_gaps, omitted_curvature = [], [], [], [], []
    positive_decreasing = True
    for hh, Z in itertools.product(protocol["heat"]["h"], protocol["heat"]["Z"]):
        derivatives, estimates = heat_derivatives(hh, Z, max(protocol["heat"]["derivative_orders"]))
        references = [heat_reference(hh, Z, j) for j in protocol["heat"]["derivative_orders"]]
        reference_gaps.append(normalized_error(np.array(derivatives), np.array(references)))
        ode_terms = [Z*Z*derivatives[2], (1+2*(1+hh)*Z)*derivatives[1], hh*(1+hh)*derivatives[0]]
        ode_gaps.append(abs(sum(ode_terms))/max(1., *(abs(v) for v in ode_terms)))
        positive_decreasing &= derivatives[0] > 0 and derivatives[1] < 0
        heat_rows.append({"h": hh, "Z": Z, "derivatives": derivatives, "special_function_reference": references,
                          "quad_error_estimates_not_certified_bounds": estimates, "ode_normalized_error": ode_gaps[-1]})
    for hh, point in itertools.product(protocol["heat"]["h"], protocol["heat"]["physical_points"]):
        result = heat_exterior(point, hh)
        full_gaps.append(float(np.max(np.abs(result["residual"]))/result["scale"]))
        curvature_gaps.append(normalized_error(result["without_swirl_curvature"], result["expected_curvature_error"]))
        omitted_curvature.append(abs(result["without_swirl_curvature"])/result["scale"])
        positive_decreasing &= result["K"] > 0 and result["Kr"] < 0
        exterior_rows.append({"h": hh, "point": point, "full_residual": result["residual"].tolist(),
                              "term_scale": result["scale"], "K": result["K"], "Kr": result["Kr"],
                              "normalized_error": full_gaps[-1],
                              "omitted_curvature_normalized_gap": omitted_curvature[-1]})
    for name, value in {"heat-reference": max(reference_gaps), "heat-ode": max(ode_gaps),
                        "heat-full-residual": max(full_gaps), "heat-curvature-identity": max(curvature_gaps),
                        "heat-positive-decreasing": float(positive_decreasing)}.items():
        checks[name]["value"] = value
    coefficients = core_coefficients(protocol["core_comparison"]["terms"])
    recurrence = all(2*(n+1)*(n+2)*coefficients[n+1]+coefficients[n] == 0
                     for n in range(len(coefficients)-1))
    tmax = Fraction(41, 20)
    lower = 1-tmax/2+tmax*tmax/12-tmax**3/144
    # P'(t)=-((t-4)^2+8)/48, and alternating terms decrease from n=1.
    positive_lower = lower == Fraction(305719, 1152000) and lower > Fraction(265, 1000) and tmax/6 < 1
    core_rows, core_gaps = [], []
    for z in protocol["core_comparison"]["z"]:
        value, reference, next_term = core_comparison(z, len(coefficients))
        core_gaps.append(normalized_error(value, reference))
        core_rows.append({"z": z, "series": value, "bessel_reference": reference,
                          "alternating_truncation_bound_excludes_roundoff": next_term})
    checks["core-reference"]["value"] = max(core_gaps)
    checks["core-exact-recurrence"]["value"] = float(recurrence)
    checks["core-exact-positive-lower-bound"]["value"] = float(positive_lower)
    for name, value in {"omitted-axial-viscosity": max(axial_gaps), "reversed-stress-sign": max(sign_gaps),
                        "omitted-swirl-curvature": max(omitted_curvature)}.items():
        checks["detect-"+name]["value"] = value
    failed = any(c["value"] is None or not math.isfinite(c["value"]) or
                 (c["kind"] == "at-most" and not 0 <= c["value"] <= c["limit"]) or
                 (c["kind"] == "at-least" and c["value"] < c["limit"]) or
                 (c["kind"] == "exact" and c["value"] != c["limit"]) for c in checks.values())
    return {"status": "building-block-audit-failed" if failed else STATUS,
            "full_proof_verified": False, "new_candidate": False,
            "prior_unforced_status": "unchanged-resolution-failure",
            "unverified_obligations": protocol["unverified_obligations"],
            "checks": checks, "manufactured_samples": fixture_rows,
            "axis_samples": len(axis_gaps), "finite_difference": fd_rows,
            "heat_profile": heat_rows, "heat_exterior": exterior_rows,
            "core_comparison": core_rows, "core_lower_bound": str(lower),
            "illustrative_scaling_not_measured_error": {
                "h": h, "tau_power_2h": {str(tau): tau**(2*h) for tau in (.05, 1e-6)},
                "log10_tau_for_factor_point_one": -1/(2*h)}}


def verify_record(path, protocol, source_hashes):
    record = json.loads(Path(path).read_text())
    if record["source_sha256"] != source_hashes:
        raise ValueError("Recorded evidence is not bound to the current source files")
    if record["paper"] != protocol["paper"]:
        raise ValueError("Paper provenance mismatch")
    validate(record, protocol)
    return record


def compare_reproduction(record, fresh, protocol):
    """Compare physical values, not machine-dependent cancellation-error ratios."""
    fields = {"manufactured_samples": ("full_cartesian_residual", "term_scale", "leading_tangential_residual", "axial_second_derivatives"),
              "heat_profile": ("derivatives", "special_function_reference"),
              "heat_exterior": ("full_residual", "term_scale", "K", "Kr"),
              "core_comparison": ("series", "bessel_reference")}
    tolerance = protocol["heat"]["maximum_reference_error"]
    for group, keys in fields.items():
        if len(record[group]) != len(fresh[group]):
            raise ValueError("Reproduction sample count mismatch")
        for old, new in zip(record[group], fresh[group]):
            for key in keys:
                if normalized_error(old[key], new[key]) > tolerance:
                    raise ValueError("Recorded physical value was not reproduced: " + group + "/" + key)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New JSON path; existing files are refused")
    parser.add_argument("--paper", type=Path, help="Optional local PDF whose bytes must match the frozen SHA256")
    parser.add_argument("--verify-record", type=Path, help="Check existing evidence's source binding and gates before independently recomputing")
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("Output already exists; use a fresh path")
    protocol = json.loads((HERE/"protocol.json").read_text())
    hashes = {name: sha256(HERE/name) for name in SOURCES}
    if args.paper and sha256(args.paper) != protocol["paper"]["sha256"]:
        raise ValueError("PDF bytes do not match the frozen paper")
    record = verify_record(args.verify_record, protocol, hashes) if args.verify_record else None
    report = audit(protocol)
    report.update({"source_sha256": hashes, "paper": protocol["paper"],
                   "paper_bytes_verified_in_this_run": bool(args.paper),
                   "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}})
    # Preserve failed finite results too. allow_nan=False prevents invalid JSON evidence.
    with args.output.open("x") as output:
        json.dump(report, output, indent=2, sort_keys=True, allow_nan=False)
        output.write("\n")
    validate(report, protocol)
    if record is not None:
        compare_reproduction(record, report, protocol)
    print(json.dumps({"status": report["status"], "checks": len(report["checks"]),
                      "manufactured_samples": len(report["manufactured_samples"]),
                      "output": str(args.output), "full_proof_verified": False}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
