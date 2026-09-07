# 128-grid resolution follow-up for the late-growth candidate

This separate experiment evolves the frozen late-growth initial field on
`N=128` through `T=.10`, comparing it with the previously verified `N=64`
trajectory. Its purpose is to resolve the near-threshold pointwise and
stretching disagreements of the earlier `32/64` comparison. No coefficients
are optimized or fitted on these held-out times.

**Status: passed the preliminary finite-amplification resolution and endpoint
sampling gates at both `.08` and `.10`.** The unchanged initial field completes
1,570 accepted steps on `N=128`, using about 65 minutes of solver time.
This is finite-growth evidence, not a proof of singularity or regularity.

## Completed evolution

All ratios below use the same initial field. L3 and vorticity use the common
physical sampling grid `128`; cutoff entries are energy fractions, not percentages.

| Quantity | T=.08, N=64 | T=.08, N=128 | T=.10, N=64 | T=.10, N=128 |
|---|---:|---:|---:|---:|
| H1/2 / initial | 1.119653950 | 1.119654159 | 1.145519113 | 1.145521605 |
| L3 / initial | .989370332 | .989372552 | .976621199 | .976652999 |
| Sampled vorticity / initial | 1.620565790 | 1.619422088 | 1.894252862 | 1.977433827 |
| Enstrophy / initial | 1.750114425 | 1.750136251 | 2.024450897 | 2.024996510 |
| Characteristic wavenumber / initial | 1.337956661 | 1.337965018 | 1.444901276 | 1.445096554 |
| Critical logarithmic growth rate | 1.255145462 | 1.255178823 | 1.026694166 | 1.026857962 |
| Enstrophy stretching / viscous destruction | 8.974793935 | 8.967917446 | 5.782149132 | 5.714409983 |
| Peak per-step cutoff fraction | 1.459355e-6 | 1.174260e-10 | 3.214034e-5 | 1.775897e-8 |
| Accepted steps | 662 | 1,130 | 880 | 1,570 |

The fine H1/2 norm grows `14.55216%` by `.10` and another `2.31031%` over
`.08--.10`, with no sampled late reversal. L3 falls `2.33470%` below its
initial value. These are separate critical norms; growth in one is not growth
in both.

## How far the numerical checks are from their limits

The gaps below are the largest relative discrepancies at common trajectory
observations through each horizon. Stretching uses the final `.02` window.
All percentages use `abs(a-b)/max(abs(a),abs(b))`.

| Comparison | Through .08 | Through .10 | Existing limit |
|---|---:|---:|---:|
| H1/2 ratio gap | .00001871% | .00021756% | 2% |
| L3 ratio gap | .00022436% | .00325603% | 2% |
| Enstrophy ratio gap | .00124708% | .02694389% | 5% |
| Characteristic wavenumber ratio gap | .00062461% | .01351319% | 5% |
| Sampled vorticity ratio gap | .17840260% | 4.20651068% | 10% |
| Late stretching gap | .04373024% | .72570452% | 10% |
| Horizon times maximum critical-rate gap | 2.668927e-6 | 1.637967e-5 | .002 |

The earlier `32/64` vorticity and stretching gaps at `.10` were `9.08869%`
and `9.42021%`. The finer comparison reduces both. However, its vorticity
gap rises sharply between `.08` and `.10`, so the near-perfect critical-norm
agreement must not be generalized to pointwise accuracy.

[`assessment.json`](assessment.json) preserves every declared gate and its
outcome. These limits control empirical screening; they are not proof targets
or a rigorous error bound.

## Denser physical sampling

Independent NumPy FFTs evaluate each saved field on physical grids `128`
and `256`, with the initial maximum recomputed on each sampling grid.
The initial maxima are `82.55311350` and `82.60944619` respectively.

| Evolved field | Vorticity ratio, 128 samples | Vorticity ratio, 256 samples | Raw maximum shift | Ratio shift |
|---|---:|---:|---:|---:|
| N=64, T=.08 | 1.620565790 | 1.622427070 | .182835% | .114722% |
| N=128, T=.08 | 1.619422088 | 1.621106344 | .172016% | .103895% |
| N=64, T=.10 | 1.894252862 | 1.893361747 | .021158% | .047043% |
| N=128, T=.10 | 1.977433827 | 1.981861400 | .291444% | .223405% |

All raw and normalized shifts are below the declared `2%` limits. The
cross-evolution vorticity-ratio gap on `256` samples is `.081404%` at `.08`
and `4.465481%` at `.10`, below the `10%` gate. Denser sampling therefore
does not explain away the remaining evolution-grid discrepancy. The raw
maxima increase on the nested grids; a normalized amplification ratio can
still decrease because its initial maximum also increases.

## Critical production and diffusion

| Time | N=128 H1/2 / initial | Relative growth rate | Critical production / viscous destruction |
|---|---:|---:|---:|
| .04 | 1.056829130 | 1.569394287 | 9.535736802 |
| .05 | 1.073393583 | 1.533054127 | 8.542737153 |
| .06 | 1.089596690 | 1.458609733 | 7.447596979 |
| .07 | 1.105088693 | 1.362427083 | 6.362403760 |
| .08 | 1.119654159 | 1.255178823 | 5.351202040 |
| .09 | 1.133161398 | 1.142499598 | 4.448974776 |
| .10 | 1.145521605 | 1.026857962 | 3.674157991 |

Independent NumPy snapshot diagnostics agree with the evolving solver to
at most `8.64e-13` relative error, below the `1e-10` check. Nonlinear critical
production remains larger than viscous destruction, but their ratio and the
relative growth rate decrease throughout the recorded budget window.

## Frozen protocol and provenance

[`protocol.json`](protocol.json) was recorded before the run. Initial energy
is `10`, Fourier bandwidth is `3`, and viscosity is `.02`. The coefficient
SHA-256 is `893ef61450d670936dc851d28ab3cf7c1aa315c6ff23934646e74137c40ac6ef`.
The new trajectory uses cutoff `42`, maximum timestep `.000125`, CFL target
`.4`, and the same `.01` physical observation clock as the earlier run.
The comparison uses the common physical sampling grid `128`.

The current executable, all recorded headers, the unchanged continuation
driver and the initial file were checked against the previous validation's
hashes before launch. `reference128/manifest.json` records exact commands,
source and executable hashes, accepted output identities, restart boundaries,
and independent NumPy critical-budget comparisons at `.04,.05,...,.10`.

The fine trajectory begins at time zero from the original coefficients.
Zero-padding an evolved coarse field would omit the finer modes' preceding
nonlinear interactions and is not the experiment performed here.

## What is checked

The historical `assess.py:paired` function compares every common observation
through `.08` and `.10`. It retains the original cutoff, critical-norm,
enstrophy, scale, vorticity, late-stretching and late-critical-rate thresholds.
All limits are empirical screening criteria, not bounds on the PDE error.

[`analyze.py`](analyze.py) also evaluates both saved endpoint fields on
physical grids `128` and `256`. It uses independent NumPy real FFTs and the
same initial-field normalization for both evolution resolutions. A denser
sampling grid evaluates the retained Fourier polynomial more densely; it
does not add missing evolution modes or certify the continuous maximum.
These extra sampling tests are at `.08` and `.10`, not at every time.

The previous smaller-timestep and full radix2/FFTW trajectory comparisons
are at `N=64`. They are reused as prior validation, not counted as new
experiments. No new full `N=128` timestep or independent trajectory study
is claimed. The independent NumPy budget checks here concern saved
snapshots, not a second complete fine-grid evolution.

## Interpretation and the next decision

The critical budget evaluates the instantaneous relative growth rate
`gamma_H = (P_H - D_H) / (2 H^2)`, where `P_H` is nonlinear production
of the squared H1/2 norm and `D_H` is its viscous destruction. Positive
`gamma_H` records growth at the sampled instant. A decreasing positive
rate records slowing relative growth; it neither proves a future turnover
nor excludes a later acceleration.

Passing a spatial comparison supports the finite numerical measurements.
It does not provide a bound on the omitted PDE modes, a continuous
vorticity maximum, or a time-independent critical-norm bound. There is no
finite numerical target value that completes a blow-up proof.

With the complete `.10` checks passed, the next bounded experiment is an
unchanged-field continuation to `.12` on both grids, tracking late H1/2 growth,
L3, the production--diffusion balance and the same resolution gates. The rising
pointwise discrepancy deserves particular attention. No longer trajectory is
part of this archive.

## Reproduction

From the repository root, use a new output directory:

```sh
python3 research/navier_stokes_cascade/results/late_growth_continuation/run.py \
  --solver ./build/research/navier_stokes_cascade/navier_stokes_continue \
  --grid 128 --output-dir new-resolution-run/reference128
```

`--solver` selects the previously built FFTW continuation executable.
`--grid 128` requests the finer evolution; its default retained cutoff is
`42`. `--output-dir` protects existing evidence by requiring a new directory.
The driver freezes the initial file, viscosity, timestep controls and
observation clock, and saves the actual evolved states at `.08` and `.10`.

The archive's `analyze.py` compares its `reference128` output with the
preceding campaign's `reference64` data. `archive.py` splits the compressed
states into bounded Git blobs and verifies the combined compressed stream.
It may pack the completed `.08` target while the `.10` trajectory continues;
that does not mark the full experiment complete. The original gzip files
are ignored local copies; the published indices and parts restore them.
`audit.py` checks source and run identities, recomputes the empirical verdict
from recorded evidence, and verifies restoration of both complete states.
Both evolved states restore exactly: each contains `100,663,792` bytes and
is preserved in ten ordered compressed parts. An archive-integrity pass is
reported separately from the scientific outcome in [`audit.json`](audit.json).

## Evolved restart states

Each final archive index describes an ordered gzip stream, individual part
hashes, the combined compressed hash and the restored-state hash. Restore
a chosen state into a fresh local file with the existing generic reader:

```sh
python3 research/navier_stokes_cascade/results/frozen_continuations_resolution/restore_checkpoint.py \
  --index research/navier_stokes_cascade/results/late_growth_resolution128/reference128/t100_checkpoint.json \
  --output fine128-restored.chk
```

`--index` selects the `.10` state; use `t080_checkpoint.json` for `.08`.
`--output` must not already exist. The restored file contains the full
evolved Fourier field, original baselines, clock and timestep settings.
Continue with `--restart` and a later absolute final time on that clock.
The archived checkpoints support reproduction; they are not proof certificates.
