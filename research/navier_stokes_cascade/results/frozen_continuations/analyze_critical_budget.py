#!/usr/bin/env python3
"""Independently evaluate the critical H1/2 budget using NumPy FFTs.

This is snapshot analysis, not another trajectory solver or a rigorous error
bound. Run after the campaign. NumPy is the only extra analysis dependency.
"""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import struct

import numpy as np


def budget(coefficients, viscosity):
    n = coefficients.shape[0]
    frequencies = np.rint(np.fft.fftfreq(n) * n)
    k = (frequencies[:, None, None], frequencies[None, :, None], frequencies[None, None, :])
    squared = k[0] ** 2 + k[1] ** 2 + k[2] ** 2
    radius = np.sqrt(squared)
    velocity = np.empty(coefficients.shape, dtype=float)
    vorticity = np.empty(coefficients.shape, dtype=float)
    for component in range(3):
        velocity[..., component] = (np.fft.ifftn(coefficients[..., component]) * n**3).real
        a, b = (component + 1) % 3, (component + 2) % 3
        omega_hat = 1j * (k[a] * coefficients[..., b] - k[b] * coefficients[..., a])
        vorticity[..., component] = (np.fft.ifftn(omega_hat) * n**3).real
    sampled_maximum = float(np.sqrt(np.sum(vorticity**2, axis=-1)).max())
    nonlinear = np.cross(velocity, vorticity)
    del velocity, vorticity
    power = np.zeros(squared.shape)
    transfer = np.zeros(squared.shape)
    for component in range(3):
        transformed = np.fft.fftn(nonlinear[..., component]) / n**3
        u = coefficients[..., component]
        power += np.abs(u) ** 2
        # u is solenoidal: its inner product with the removed pressure term
        # is zero, so no shared Leray-projector implementation is needed here.
        transfer += 2.0 * (u.conj() * transformed).real
    q = float(np.sum(radius * power))
    production = float(np.sum(radius * transfer))
    destruction = float(2 * viscosity * np.sum(radius**3 * power))
    return {
        "energy": float(0.5 * power.sum()),
        "enstrophy": float(0.5 * np.sum(squared * power)),
        "h_half": q**0.5,
        "h_half_squared_nonlinear_production": production,
        "h_half_squared_viscous_destruction": destruction,
        "h_half_squared_net_rate": production - destruction,
        "h_half_logarithmic_rate": (production - destruction) / (2 * q),
        "critical_production_to_dissipation": production / destruction,
        "enstrophy_stretching": float(0.5 * np.sum(squared * transfer)),
        "enstrophy_viscous_destruction": float(viscosity * np.sum(squared**2 * power)),
        "sampled_vorticity_native": sampled_maximum,
        "nonlinear_energy_residual": float(abs(0.5 * transfer.sum())),
    }


def load_checkpoint(path, expected_sha256):
    raw = gzip.decompress(path.read_bytes())
    # The campaign manifest records the exact bytes emitted by the validated
    # C++ checkpoint writer, including its internal checksum. SHA-256 also
    # detects any changed header, metadata, coefficient or trailing byte here.
    if hashlib.sha256(raw).hexdigest() != expected_sha256 or raw[:8] != b"NSCONT1\n":
        raise ValueError("Checkpoint identity or format changed")
    n, cutoff, sampling, dense, backend = struct.unpack_from("<5Q", raw, 24)
    length = struct.unpack_from("<Q", raw, 8)[0]
    state_size = struct.unpack_from("<Q", raw, 488)[0]
    if length != len(raw) - 24 or state_size != n**3 or len(raw) != 496 + 48 * state_size:
        raise ValueError("Checkpoint length is inconsistent")
    if 3 * cutoff >= n or sampling != n or backend != 1:
        raise ValueError("Budget comparison needs the dealiased native FFTW sampling grid")
    viscosity = struct.unpack_from("<d", raw, 64)[0]
    observed = struct.unpack_from("<19d", raw, 336)
    field = np.frombuffer(raw, dtype="<c16", offset=496).reshape(n, n, n, 3)
    if not np.isfinite(field).all():
        raise ValueError("Nonfinite Fourier coefficients")
    return field, viscosity, observed, {"grid": n, "cutoff": cutoff, "dense_sampling_grid": dense}


def load_initial(path, n):
    field = np.zeros((n, n, n, 3), dtype=complex)
    with path.open() as source:
        for row in csv.DictReader(source):
            index = tuple(int(row[key]) % n for key in ("kx", "ky", "kz"))
            field[index] = [complex(float(row[c + "_real"]), float(row[c + "_imag"]))
                            for c in ("ux", "uy", "uz")]
    return field


def test_exact_shear():
    field = np.zeros((8, 8, 8, 3), dtype=complex)
    field[0, 1, 0, 0] = field[0, -1, 0, 0] = 0.5
    result = budget(field, 0.02)
    assert abs(result["energy"] - 0.25) < 1e-14
    assert abs(result["h_half_logarithmic_rate"] + 0.02) < 1e-14
    assert abs(result["h_half_squared_nonlinear_production"]) < 1e-14
    assert abs(result["sampled_vorticity_native"] - 1) < 1e-14


def main():
    test_exact_shear()
    output = Path(__file__).resolve().parent
    research = output.parents[1]
    manifest = json.loads((output / "manifest.json").read_text())
    report = {
        "equation": "Q=||u||_Hhalf^2; Q'=2 Re sum |k| conj(u_k).N_k - 2 nu sum |k|^3 |u_k|^2",
        "measure": "Normalized periodic volume (2*pi)^-3 dx",
        "numpy_version": np.__version__,
        "analysis_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "analytic_shear_test": "passed",
        "interpretation": "Instantaneous finite-mode critical-norm budget; sampled spatial maxima and floating-point values are not continuum error bounds.",
        "candidates": [],
    }
    for candidate in manifest["candidates"]:
        checkpoint = next(item for item in candidate["checkpoints"] if item["grid"] == 64)
        field, nu, observed, metadata = load_checkpoint(output / checkpoint["path"], checkpoint["uncompressed_sha256"])
        result = budget(field, nu)
        comparisons = {"energy": 1, "enstrophy": 2, "h_half": 4,
                       "sampled_vorticity_native": 6, "enstrophy_stretching": 11,
                       "enstrophy_viscous_destruction": 12}
        errors = {name: abs(result[name] - observed[index]) / max(abs(observed[index]), 1e-300)
                  for name, index in comparisons.items()}
        if max(errors.values()) > 1e-10:
            raise ValueError("Independent snapshot budget failed agreement with FFTW: " + str(errors))
        initial = budget(load_initial(research / candidate["initial_state"], 64), nu)
        report["candidates"].append({"name": candidate["name"], "time": observed[0], **metadata,
                                      "checkpoint_sha256": checkpoint["uncompressed_sha256"],
                                      "initial": initial, "final": result,
                                      "numpy_vs_fftw_relative_errors": errors,
                                      "pair_status": candidate["stages"][-1]["status"]})
    (output / "critical_budget.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for item in report["candidates"]:
        print(item["name"], "initial/final critical production-to-dissipation",
              item["initial"]["critical_production_to_dissipation"], item["final"]["critical_production_to_dissipation"],
              "final d(log H)/dt", item["final"]["h_half_logarithmic_rate"], flush=True)


if __name__ == "__main__":
    main()
