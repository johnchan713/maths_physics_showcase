# Candidate-search evidence

These records support numerical finite-amplification experiments, not a PDE
singularity claim. The gates and equations are documented in the
[research README](../README.md#robust-amplification-search).

## Verified seven-seed pilot

[`robust_amplification_verified/manifest.json`](robust_amplification_verified/manifest.json)
is the authoritative record. It retains every command, implementation and
executable hash, selected coefficient hash, rejected gate and outcome. Each
run directory contains `trace.csv`, `evidence.csv`, `state.csv` and `run.log`.
The manifest is updated after each run, so partial output alone must never be
interpreted as a completed validation.

Controls were `E=10`, `nu=0.02`, initial bandwidth `K=3` (684 real variables),
evolution grids `N=16/32`, one adjoint iteration per seed, discovery `T=0.04`
and held-out `T=0.06`. Two deterministic starts were tested for each of the
packet, tube and orthogonal-bundle seeds, plus the previous smooth-maximum
packet checkpoint. The seven discoveries and nine held-out replays all
completed. Three seeds reached the non-dominated frontier and all three
passed the predeclared **finite-amplification** gates, including separate
independent FFTW trajectories. All three still failed the profile route.

The ranges below cover all six held-out trajectories per selected state:
both grids, the smaller timestep rerun, and the denser-sampling rerun.

| Frozen optimized seed | Final H1/2 ratio range | Final sampled L3 ratio range | Minimum scale ratio | Minimum late production/destruction | Worst per-step cutoff fraction |
|---|---:|---:|---:|---:|---:|
| Previous packet checkpoint | 1.056956669–1.057465591 | 0.994691171–0.994904529 | 1.158313104 | 12.347309 | 0.005913063 |
| Fresh packet, start 0 | 1.051314839–1.051559324 | 0.993269956–0.993583449 | 1.156928619 | 14.733717 | 0.005197318 |
| Tube, start 0 | 1.015182390–1.015588460 | 0.986194232–0.986450344 | 1.058132856 | 5.305355 | 0.003571197 |

All fields remain smooth truncated solutions and all final L3 ratios are
below one. Sampling an increasing critical norm at these times is not an
unbounded-growth result. The cutoff fractions are empirical triage evidence,
not Fourier-tail bounds.

The fine-grid timestep check increased accepted steps from `201 -> 408`
(checkpoint), `169 -> 344` (fresh packet), and `149 -> 304` (tube). Every
coarse reference used 120 steps and its timestep check used 240. The denser
spatial grid was `64^3` instead of `32^3`; fixed-time snapshots increased from
8 to 16. Extra shortened steps at snapshot boundaries explain that rerun's
different accepted-step counts. Endpoint comparisons, rather than sampled
peak times, govern the common-grid agreement tests.

Independent `N=32, T=0.06` FFTW/RK4 comparisons had peak relative whole-state
differences `6.5340e-16`, `7.0420e-16` and `5.2246e-16` for the checkpoint,
fresh-packet and tube states respectively. The largest adjoint/difference
slope discrepancy in discovery was `1.4423e-5`, below the `0.01` gate.

The leading state is also copied to
[`wave_k3_amplification_t004_screening.csv`](../candidates/wave_k3_amplification_t004_screening.csv).
Its SHA-256 is
`5f3dc4449abd58f5ae5e1aca3443b2a0fa7ad06029f9dadae939ff2b9a3b51d8`.
It is byte-identical to the verified campaign's selected checkpoint state.
Optimization family names describe the seed, not a geometric constraint on
the final Fourier coefficients.

## Longer and larger-resolution extension

The same frozen leading state is evaluated separately at `T=0.08` on
`N=32/64`, with maximum timestep `0.000125` and common spatial samples on
`N=64`. Commands, raw data and the independent `N=64` FFTW comparison live in
[`robust_amplification_extended`](robust_amplification_extended).
`run_extension.py` reproduces those commands with supplied executable paths;
use a new output directory. This additional spatial and horizon check does
not replace the pilot's separate timestep and sampling perturbations, which
were performed at `N=16/32, T=0.06`.

| Extension measurement at `T=0.08` | `N=32` | `N=64` |
|---|---:|---:|
| Final H1/2 / initial | 1.079044840 | 1.079088316 |
| Final sampled L3 / initial | 0.986383199 | 0.986594453 |
| Final enstrophy / initial | 1.493309729 | 1.495045195 |
| Final sampled vorticity / initial, common `64^3` samples | 1.506446678 | 1.604798100 |
| Characteristic wavenumber / initial | 1.233429178 | 1.234147383 |
| Minimum sampled late production/destruction | 9.648878953 | 9.500616478 |
| Peak per-step cutoff fraction | 0.000276205785 | 0.000000671045 |
| Latest rescaled-spectrum drift | 7.691962193 | 7.684392492 |
| Accepted steps | 640 | 667 |

The H1/2 ratio differs by `0.00403%` across these grids, but sampled vorticity
still differs by `6.13%`. The latter is an important remaining convergence
limitation even though it passes this workflow's 10% preliminary gate.
Tiny cutoff energy does not guarantee a comparably tiny pointwise vorticity
error. Both profile drifts remain well above one, and L3 still decreases.
This is a reproducible finite-cascade candidate for further investigation,
not evidence that a singularity will form.

The independent `N=64, T=0.08` FFTW/RK4 comparison passed with 666 shared
safe steps. Peak relative whole-state disagreement was `1.4092e-15`; peak
scaled diagnostic disagreement was `1.9137e-15`. The paired run takes one
extra step because it also lands on fixed evidence-snapshot times. Both
reach the same final physical time. See `fftw.csv`, `fftw.log` and the zero
return codes in `commands.json` for the recorded check.

The next useful experiment is to continue the frozen finalists while
monitoring late critical growth and reducing the remaining pointwise
vorticity discrepancy. Longer continuations should use the existing
checksummed restart path. A proof would additionally require rigorous
control of the unresolved Fourier tail and of the limiting PDE mechanism;
none of these finite runs supplies that control.
