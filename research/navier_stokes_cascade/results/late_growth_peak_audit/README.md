# Why the .14 peak comparison fails

This is a diagnosis and bounded cost probe, **not a new candidate or a new
long-time evolution**. The initial coefficients, evolution code, and original
10% scalar-peak agreement gate remain unchanged. The preceding
[64/128 comparison at .14](../late_growth_t014/README.md) still fails.

The main finding is that modes absent from the coarse evolution carry very
little kinetic energy but materially affect vorticity. At `.14`, removing
those modes from the saved fine field lowers its sampled maximum from
`229.882` to `194.051`, a `15.59%` reduction. They contain only `0.02687%` of
the fine field's energy. The values of modes retained by both evolutions
also differ: simply padding an evolved coarse state would not reproduce a
genuine fine-grid trajectory.

## Frozen inputs and interpretation

[protocol.json](protocol.json) was frozen before the snapshot analysis and
128/256 cost probe. It identifies the six archived states at `.12`, `.13`,
and `.14`, their hashes, the source hashes, and the resource/time limits.
The initial file is
[`wave_k3_late_rate_t004_screening.csv`](../../candidates/wave_k3_late_rate_t004_screening.csv),
SHA-256 `893ef61450d670936dc851d28ab3cf7c1aa315c6ff23934646e74137c40ac6ef`.
Its initial energy is 10, Fourier bandwidth 3, and viscosity `.02`.

For the saved 64- and 128-grid velocity fields, we compute

$$
\omega_{128}-\omega_{64}
=\nabla\times(P_{21}u_{128}-u_{64})
+\nabla\times((I-P_{21})u_{128}).
$$

Here `P21` is the **cube** cutoff `max(|kx|, |ky|, |kz|) <= 21`, not a
radial cutoff. The first term is the shared-mode error; the second contains
modes retained only by the fine evolution, whose cube cutoff is 42.
Their integrated squared errors add by disjoint Fourier support. Their
pointwise maxima do not add: vectors can cancel, and maxima can move.
All physical-space maxima below use the same `256^3` sampling lattice.
They are not certified continuous maxima or bounds on unresolved PDE modes.

## A scalar peak is not the full vorticity field

| Time | Coarse maximum | Fine maximum | Relative scalar-peak gap | Maximum vector difference / fine peak |
|---|---:|---:|---:|---:|
| .12 | 185.1382 | 189.1651 | 2.129% | 19.637% |
| .13 | 198.7221 | 203.8719 | 2.526% | 22.374% |
| .14 | 200.7923 | 229.8819 | 12.654% | 21.339% |

The scalar gap compares the magnitudes of two possibly different peaks,
using the larger magnitude as denominator. The last column compares vectors
at matching spatial points, maximizes their difference, and divides by the
fine-grid peak. Similar peak magnitudes did not establish uniformly similar
sampled vector fields, even at the earlier times.

The vector comparison is a **post-hoc diagnostic**, not a retroactive change
to the earlier pass/fail rules. Future search protocols should include it
prospectively, with an explicit tolerance and fixed observation times.

## Small energy does not imply a small derivative contribution

| Time | Extra-mode energy / fine energy | Extra-mode enstrophy / fine enstrophy | Shared-mode portion of squared vorticity error |
|---|---:|---:|---:|
| .12 | 0.010989% | 0.865400% | 33.2941% |
| .13 | 0.016724% | 1.277439% | 34.4217% |
| .14 | 0.026872% | 1.876547% | 34.8008% |

Curl weights Fourier coefficients by wavenumber, making vorticity more
sensitive to high modes than energy. The last column also shows why missing
modes are not the whole discrepancy: the shared modes have evolved
differently. These integrated error portions are not portions of the
difference between the two scalar maxima.

At the fine peak at `.14`, the extra-mode vector contributes `117.42%`
along the direction of the vector error, while the shared-mode error
contributes `-17.42%`. This is signed projection and partial cancellation,
not an attribution of the scalar peak gap. Full vectors and locations are
recorded in [diagnosis.json](diagnosis.json).

### Filtering the same fine field

| Time | Cutoff 21 | Cutoff 28 | Cutoff 35 | Cutoff 42 (full) |
|---|---:|---:|---:|---:|
| .12 | 184.7623 | 182.3791 | 186.7649 | 189.1651 |
| .13 | 195.2411 | 207.0053 | 204.9953 | 203.8719 |
| .14 | 194.0505 | 222.1375 | 230.3497 | 229.8819 |

These are post-hoc filtered sampled maxima, not four independently evolved
resolutions. Nonmonotonicity is possible because modes reinforce or cancel
and peak locations shift. In particular, the close cutoff-35/42 values at
`.14` do **not** bound the influence of modes beyond 42 in a finer evolution.

## Measured cost of the next supported grid

The unchanged FFTW RK4 path took two matching safe steps from the original
initial coefficients at each grid, with `dt=1e-5`, ending at `.00002`.
This is early-time capacity/agreement evidence only; no 256-grid trajectory
through `.14` was run. Padding this bandwidth-3 **initial** field is exact;
padding a previously evolved coarse field would be a different experiment.

| Measurement | 128-grid | 256-grid |
|---|---:|---:|
| Cube cutoff | 42 | 85 |
| First RK4 step | 2.4371 s | 98.3407 s |
| Second RK4 step | 2.4671 s | 92.9607 s |

The mean step-time ratio is `39.01`, the whole-state relative L2 difference
is `1.13e-18`, and the maximum coefficient difference is `3.47e-18`.
The combined probe completed in `199.11 s` with peak child RSS `5.052 GiB`,
inside its `300 s` wall timeout and `8 GiB` address-space cap. Fourier
reality, divergence, energy, and cutoff checks passed. The memory report
does not establish what a full production run with all diagnostics needs.

A rough cost model scales the previous 128-grid step count of 2749 by
`85/42`, assuming comparable velocity bounds. Combining the resulting
5563 modeled steps with the two measured fine-step times gives
**143.7–152.0 hours of stepping**, roughly six days. Even a hypothetical
1120-step run at the maximum allowed `dt=.000125` would take 28.9–30.6 hours
at those measured step costs. These are scenarios, not completion promises
or confidence intervals. They exclude initialization, observations,
checkpoint I/O, load variation, and differences in the actual fine trajectory.

Raw results: [probe.json](probe.json), [probe.stdout](probe.stdout),
[probe.log](probe.log), and [build provenance](build.json).
No threading or backend optimization was introduced in this experiment.

## Checks and reproducibility

- Thirteen analytic/adversarial tests cover Fourier curl normalization,
  direct sums, embedding, filters, Parseval splitting, vector cancellation,
  signed projections, and invalid/nonfinite inputs.
- The short `32/64` control passes (`5.23e-14` relative state gap). An
  initially attempted `16/32` smoke failed because generated modes are lost
  at the coarser cutoff; it is retained as a discovered negative control
  (`6.03e-6` gap). An invalid `16/64` pair is rejected before evolution.
- Strict-warning compilation and the three controls under standalone UBSan
  pass. The combined ASan/UBSan run is **incomplete**: LeakSanitizer cannot
  inspect this runner's process information, and attempts to disable its
  leak phase did not resolve that environment failure. No ASan-clean or
  leak-clean claim is made. See [verification.log](verification.log) and
  [UBSan control evidence](ubsan_controls.json).
- [audit.py](audit.py) verifies source/state hashes, recomputes Fourier
  norms (including independent `|k cross u|^2` sums), verifies recorded point
  vectors by direct Fourier sums without inverse FFTs, and compares the
  recorded scalar maxima with the prior assessments. It does **not** repeat
  every dense maximum search or certify continuous extrema.

From the repository root, with NumPy and FFTW development files installed:

```sh
python3 research/navier_stokes_cascade/results/late_growth_peak_audit/test_analysis.py -v
c++ -std=c++11 -O3 -Wall -Wextra -Wpedantic -Werror \
  -I research/navier_stokes_cascade/include \
  research/navier_stokes_cascade/results/late_growth_peak_audit/benchmark.cpp \
  -lfftw3 -o navier_stokes_peak_probe
python3 research/navier_stokes_cascade/results/late_growth_peak_audit/test_probe.py \
  --binary ./navier_stokes_peak_probe --output navier_stokes_peak_probe_controls.json
python3 research/navier_stokes_cascade/results/late_growth_peak_audit/audit.py \
  --output navier_stokes_peak_audit_ci.json
```

Use fresh output paths. CI runs these checks, not the costly 128/256 probe.
This local build linked the installed versioned FFTW library because the
unversioned `-lfftw3` development linker name was unavailable.
`probe.py` and `analyze.py` refuse to replace the completed experiment records;
any full rerun should use a separate working copy preserving these originals.

## Next research step

1. Profile and improve the higher-resolution evolution backend; compare its
   states, RK4 steps, invariants, and checkpoint behavior with the unchanged
   implementation before using any performance gain scientifically.
2. Freeze a prospective vector-vorticity agreement diagnostic alongside the
   existing critical-growth, cutoff, timestep, and scalar-peak checks. Do not
   select its threshold merely to pass this candidate.
3. Once a monitored full replay is feasible, evolve the **same original
   coefficients** at 256 through `.14` and compare with 128 at fixed times.
   Do not continue the failed 64/128 pair to a longer horizon first.

The prior finite H1/2 gain of `18.24%`, decreasing L3, and slowing relative
critical growth are unchanged. This diagnosis improves numerical screening;
it supplies neither an unresolved Fourier-tail bound nor a proof of a PDE
singularity or global regularity. There is no meaningful percentage-to-proof.
