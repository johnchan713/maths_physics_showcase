#!/usr/bin/env python3
"""Separate shared-mode and extra-mode vorticity differences without evolving."""
import importlib.util
import json
from pathlib import Path
import tempfile

import numpy as np

from probe import ROOT, RESEARCH, save, sha


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


reader = module("peak_budget_reader", RESEARCH / "results/frozen_continuations/analyze_critical_budget.py")
restorer = module("peak_restorer", RESEARCH / "results/frozen_continuations_resolution/restore_checkpoint.py")


def frequencies(n):
    return np.rint(np.fft.fftfreq(n) * n).astype(int)


def require_field(field):
    n = field.shape[0]
    if field.shape != (n, n, n, 3) or n < 8 or n & (n - 1) or not np.isfinite(field).all():
        raise ValueError("Need a finite cubic power-of-two Fourier vector field")
    return n


def low_pass(field, cutoff):
    n = require_field(field)
    if cutoff < 0 or 3 * cutoff >= n:
        raise ValueError("Invalid dealiased filter cutoff")
    f = np.abs(frequencies(n))
    keep = (f[:, None, None] <= cutoff) & (f[None, :, None] <= cutoff) & (f[None, None, :] <= cutoff)
    result = field.copy()
    result[~keep] = 0
    return result


def embed(field, n):
    old = require_field(field)
    if n < old or n % old or n & (n - 1):
        raise ValueError("Embedding must preserve all modes on a nested grid")
    result = np.zeros((n, n, n, 3), dtype=complex)
    indices = frequencies(old) % n
    result[np.ix_(indices, indices, indices)] = field
    return result


def spectral_norms(field):
    n = require_field(field)
    f = frequencies(n)
    k2 = f[:, None, None]**2 + f[None, :, None]**2 + f[None, None, :]**2
    power = np.sum(np.abs(field)**2, axis=-1)
    return {"velocity_l2_squared": float(power.sum()),
            "h_half_squared": float(np.sum(np.sqrt(k2) * power)),
            "vorticity_l2_squared": float(np.sum(k2 * power))}


def sample_curl(field, sampling_grid):
    n = require_field(field)
    m = sampling_grid
    if m < n or m % n or m & (m - 1):
        raise ValueError("Sampling must use a nested power-of-two grid")
    f = frequencies(n)
    k = (f[:, None, None], f[None, :, None], np.arange(n // 2 + 1)[None, None, :])
    positive = field[:, :, :n // 2 + 1, :]
    locations = np.ix_(f % m, f % m, np.arange(n // 2 + 1))
    result = np.empty((m, m, m, 3), dtype=float)
    for c in range(3):
        a, b = (c + 1) % 3, (c + 2) % 3
        spectrum = np.zeros((m, m, m // 2 + 1), dtype=complex)
        spectrum[locations] = 1j * (k[a] * positive[..., b] - k[b] * positive[..., a])
        result[..., c] = np.fft.irfftn(spectrum, s=(m, m, m), axes=(0, 1, 2)) * m**3
    return result


def peak_sum(*terms):
    """Maximum of a sum of vectors, never the sum of their separate maxima."""
    squared = np.zeros(terms[0].shape[:-1])
    for c in range(3):
        values = sum(term[..., c] for term in terms)
        squared += values * values
    index = tuple(int(i) for i in np.unravel_index(np.argmax(squared), squared.shape))
    return {"maximum": float(np.sqrt(squared[index])), "index": list(index)}


def direct_curl(field, points, sampling_grid):
    """Independent Fourier sums at selected points; no inverse FFT is used."""
    n = require_field(field)
    f = frequencies(n)
    k = (f[:, None, None], f[None, :, None], f[None, None, :])
    phases = [[np.exp(2j * np.pi * f * index / sampling_grid) for index in p] for p in points]
    result = np.empty((len(points), 3), dtype=complex)
    for c in range(3):
        a, b = (c + 1) % 3, (c + 2) % 3
        omega = 1j * (k[a] * field[..., b] - k[b] * field[..., a])
        for j, phase in enumerate(phases):
            result[j, c] = np.einsum("i,j,k,ijk->", *phase, omega, optimize=True)
    return result


def projection_fractions(shared, high):
    shared, high = np.asarray(shared), np.asarray(high)
    total = shared + high
    squared = float(np.dot(total, total))
    if squared == 0:
        return None
    return {"shared": float(np.dot(shared, total) / squared),
            "extra": float(np.dot(high, total) / squared)}


def relative(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def compare_fields(coarse, fine, cutoff, sampling_grid, sweep):
    n = require_field(fine)
    lifted = embed(coarse, n)
    common = low_pass(fine, cutoff)
    extra = fine - common
    shared_error = common - lifted
    total_error = fine - lifted
    fourier_identity = float(np.max(np.abs(total_error - shared_error - extra)))
    norms = {key: spectral_norms(field) for key, field in
             (("coarse", lifted), ("fine", fine), ("extra", extra),
              ("shared_error", shared_error), ("total_error", total_error))}
    splits = {key: relative(norms["total_error"][key], norms["extra"][key] + norms["shared_error"][key])
              for key in norms["fine"]}
    del lifted, common, total_error
    fields = {"coarse": coarse, "fine": fine, "shared_error": shared_error, "extra": extra}
    physical = {key: sample_curl(field, sampling_grid) for key, field in fields.items()}
    peaks = {key: peak_sum(value) for key, value in physical.items()}
    identity = max(float(np.max(np.abs(physical["fine"][..., c] - physical["coarse"][..., c] -
                    physical["shared_error"][..., c] - physical["extra"][..., c]))) for c in range(3))
    identity /= max(peaks["fine"]["maximum"], peaks["coarse"]["maximum"], 1e-300)
    peaks["total_error"] = peak_sum(physical["shared_error"], physical["extra"])
    peaks["filtered_fine"] = peak_sum(physical["coarse"], physical["shared_error"])
    points = [peaks[key]["index"] for key in ("coarse", "fine", "total_error")]
    direct_errors = {}
    for key, field in fields.items():
        direct = direct_curl(field, points, sampling_grid)
        samples = np.array([physical[key][tuple(index)] for index in points])
        direct_errors[key] = float(np.max(np.abs(direct - samples))) / max(peaks["fine"]["maximum"], 1e-300)
    pointwise = []
    for label, point in zip(("coarse_peak", "fine_peak", "maximum_vector_difference"), points):
        vectors = {key: value[tuple(point)].tolist() for key, value in physical.items()}
        pointwise.append({"location": label, "index": point,
                          "coordinates": [float(2 * np.pi * p / sampling_grid) for p in point],
                          "vectors": vectors,
                          "projected_error_fractions": projection_fractions(vectors["shared_error"], vectors["extra"])})
    del physical
    filters = []
    for k in sweep:
        restricted = low_pass(fine, k)
        filtered_norms = spectral_norms(restricted)
        if k == cutoff:
            peak = peaks["filtered_fine"]
        elif k == (n - 1) // 3:
            peak = peaks["fine"]
        else:
            samples = sample_curl(restricted, sampling_grid)
            peak = peak_sum(samples)
            del samples
        filters.append({"cutoff": k, **peak,
                        "relative_maximum_gap_to_full": relative(peak["maximum"], peaks["fine"]["maximum"]),
                        "removed_energy_fraction": max(0., 1 - filtered_norms["velocity_l2_squared"] / norms["fine"]["velocity_l2_squared"])})
    return {"norms": norms, "peaks": peaks, "pointwise": pointwise, "filtered_fine": filters,
            "high_mode_fractions_of_fine_squared_norms": {key: norms["extra"][key] / norms["fine"][key] for key in norms["fine"]},
            "shared_mode_fractions_of_squared_error": {key: norms["shared_error"][key] / max(norms["total_error"][key], 1e-300) for key in norms["fine"]},
            "fourier_identity_maximum_coefficient_error": fourier_identity,
            "spectral_split_relative_errors": splits,
            "sampled_curl_identity_relative_error": identity,
            "direct_fourier_relative_errors": direct_errors,
            "sampled_vorticity_relative_gap": relative(peaks["coarse"]["maximum"], peaks["fine"]["maximum"])}


def main():
    protocol = json.loads((ROOT / "protocol.json").read_text())
    output = ROOT / "diagnosis.json"
    if output.exists():
        raise ValueError("Refusing to replace a completed or partial diagnosis")
    for name, digest in protocol["source_sha256"].items():
        if sha(RESEARCH / name) != digest:
            raise ValueError("Frozen dependency changed: " + name)
    initial_path = RESEARCH / protocol["initial_state"]
    if sha(initial_path) != protocol["initial_state_sha256"]:
        raise ValueError("Frozen initial field changed")
    settings = protocol["snapshot_audit"]
    initial_samples = sample_curl(reader.load_initial(initial_path, 64), settings["sampling_grid"])
    baseline = peak_sum(initial_samples)
    del initial_samples
    report = {"status": "analysis-in-progress", "protocol_sha256": sha(ROOT / "protocol.json"),
              "analysis_sha256": sha(Path(__file__)), "initial_state_sha256": sha(initial_path),
              "numpy_version": np.__version__, "initial_vorticity_maximum": baseline,
              "cases": [], "original_vorticity_limit": .10,
              "scope": "Snapshot diagnosis only. Filtering does not replay a trajectory or certify missing PDE modes.",
              "projection_note": "Signed contributions along the vector error at one fixed point; not fractions of the difference between separate global maxima. They may be negative or greater than one."}
    save(output, report)
    with tempfile.TemporaryDirectory(prefix="peak-states-", dir=ROOT) as directory:
        for time in settings["times"]:
            fields, provenance, observations = {}, {}, {}
            for item in [c for c in protocol["checkpoints"] if c["time"] == time]:
                index = RESEARCH / item["index"]
                if sha(index) != item["index_sha256"]:
                    raise ValueError("Input archive index changed")
                path = Path(directory) / ("n%d-t%03d.chk" % (item["grid"], round(time * 1000)))
                restored = restorer.restore(index, path)
                if restored["uncompressed_sha256"] != item["state_sha256"]:
                    raise ValueError("Input evolved state changed")
                field, viscosity, observed, layout = reader.load_checkpoint(path, item["state_sha256"])
                if observed[0] != time or layout["grid"] != item["grid"] or viscosity != .02:
                    raise ValueError("Unexpected saved state")
                fields[item["grid"]] = field
                observations[item["grid"]] = observed
                provenance[str(item["grid"])] = item
                path.unlink()
            result = compare_fields(fields[64], fields[128], settings["shared_cutoff"],
                                    settings["sampling_grid"], settings["cutoff_sweep"])
            errors = {}
            for n, name in ((64, "coarse"), (128, "fine")):
                norms, observed = result["norms"][name], observations[n]
                errors[name] = {"energy": relative(.5 * norms["velocity_l2_squared"], observed[1]),
                                "enstrophy": relative(.5 * norms["vorticity_l2_squared"], observed[2]),
                                "h_half": relative(norms["h_half_squared"]**.5, observed[4])}
            residuals = [result["fourier_identity_maximum_coefficient_error"],
                         result["sampled_curl_identity_relative_error"],
                         *result["spectral_split_relative_errors"].values(),
                         *result["direct_fourier_relative_errors"].values(),
                         *(value for group in errors.values() for value in group.values())]
            if max(residuals) > settings["identity_tolerance"] or not all(np.isfinite(residuals)):
                raise ValueError("Snapshot identity or independent evaluation failed")
            case = {"time": time, "sampling_grid": settings["sampling_grid"], "sources": provenance,
                    "snapshot_relative_errors": errors, **result,
                    "vorticity_amplification_ratios": {key: result["peaks"][key]["maximum"] / baseline["maximum"] for key in ("coarse", "fine", "filtered_fine")},
                    "original_endpoint_vorticity_gate_failed": result["sampled_vorticity_relative_gap"] > .10}
            report["cases"].append(case)
            save(output, report)
            print("Analyzed", time, "gap", case["sampled_vorticity_relative_gap"],
                  "extra-mode enstrophy fraction", case["high_mode_fractions_of_fine_squared_norms"]["vorticity_l2_squared"], flush=True)
            del fields, field
    if not report["cases"][-1]["original_endpoint_vorticity_gate_failed"]:
        raise ValueError("Diagnosis unexpectedly changed the established .14 failure")
    report["status"] = "diagnosis-complete-original-resolution-failure-retained"
    save(output, report)
    print(report["status"], flush=True)


if __name__ == "__main__":
    main()
