# Candidate-search evidence

These records support numerical finite-amplification experiments, not a PDE
singularity claim. The gates and equations are documented in the
[research README](../README.md#robust-amplification-search).

## Separate paper reproduction track

The [paper building-block audit](paper_profile_audit/README.md) checks
explicit similarity identities, a scalar comparison function and the
viscous exterior from the user-supplied OpenAI manuscript. The source-bound
record passes 29 tests and 25 numerical/algebraic gates, including deliberate
missing-term controls. It does not reproduce the nonlinear matched profile,
verify the complete proof, or promote a search candidate. The unforced
research results below are unchanged.

## Late-growth objective comparison

The [peak-error decomposition and bounded cost probe](late_growth_peak_audit/README.md)
explain why low energy contamination does not rescue the failed `.14` peak
comparison. Modes absent from the coarse grid carry only `0.02687%` of fine
energy, but removing them lowers the fine sampled vorticity peak by `15.59%`.
The maximum sampled vector difference is `21.34%` of the fine peak, exceeding
the `12.65%` scalar-peak gap; shared Fourier modes differ as well. Two tiny
128/256 steps pass, but measured 256-grid step costs are `93–98 s`, making
a full same-field replay a multi-day task under a rough cost model. No full
256-grid `.14` evolution or new candidate is claimed. The next prerequisite
is a faster, independently validated high-resolution backend, with a
prospective vector-field agreement check in future search protocols.

The [frozen .12-to-.14 continuation](late_growth_t014/README.md) passes all
declared checks at `.13` but fails the local-vorticity agreement gate at `.14`.
The `64/128` discrepancy is `12.39%` on common physical samples and `12.65%`
on doubled samples, above the unchanged `10%` limit. Fine H1/2 grows `18.24%`,
but L3 falls `5.23%` and relative critical growth slows to `.557`. All four
new states are archived. Further time promotion of this grid pair pauses;
the next requirement is a finer evolution comparison at this same horizon.
The failed numerical gate neither proves nor rules out a PDE singularity.

The [frozen .10-to-.12 continuation](late_growth_t012/README.md) passes the
unchanged `64/128` and endpoint sampling gates. Fine-grid H1/2 grows `16.65%`,
while L3 falls `3.77%` and relative critical growth slows to `.792`. The worst
common-sampling vorticity gap is `7.19%` at `.11`; the denser endpoint gap at
`.12` is `2.13%`. All four new evolved states are archived with verified
restoration and independent snapshot budgets. These remain preliminary
finite-amplification results, not a proof or certified PDE error bound.

The [frozen late-growth continuation](late_growth_continuation/README.md)
extends the same field through held-out `.08` and `.10`. Fine-grid H1/2
growth reaches 14.55%, but its relative growth rate keeps falling and L3
decreases by 2.34%. The `32/64` spatial gates pass narrowly at `.10`, with
vorticity and late-stretching discrepancies of 9.09% and 9.42%. The report
and archive audit record the independent-solver and smaller-step checks.
The [completed 64/128 follow-up](late_growth_resolution128/README.md) passes
all retained spatial and additional endpoint sampling gates. Through `.10`,
the H1/2 ratio gap is `0.00022%`, the common-sampling vorticity gap is `4.21%`
and the late-stretching gap is `0.73%`. Fine-grid H1/2 still grows `14.55%`,
while L3 falls `2.33%` and the relative H1/2 growth rate declines to `1.027`.
Actual fine-grid states at `.08` and `.10` are archived with verified restoration.

The [matched late-growth pilot](late_growth_pilot/README.md) integrates the
critical-rate source into the full discrete adjoint. Four matched starts
with two iterations each improve their weakest sampled late growth by
1.7--2.8% against endpoint-objective controls. These are small finite
improvements. Neither low-grid arm passes all held-out checks; the raw
failures and the separate frozen-leader resolution follow-up are retained.
That follow-up passes all declared `32/64, T=.06` checks, with 8.96%
H1/2 growth and independent full-state disagreement `1.41e-15`. The
logarithmic growth rate already slows, and L3 remains slightly below its
initial value. A finite-amplification label is not a proof claim.

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

## Checkpointed frozen-field continuations

The bounded `32/64` continuation campaign is complete; see
[`frozen_continuations/README.md`](frozen_continuations/README.md).
The fresh packet and tube pass the continuation-stage gates through `T=.12`
with about 11.04% and 5.76% fine-grid H1/2 growth. The prior packet exceeds
the vorticity convergence tolerance at `.10`; the other two do so by `.16`.
All three still grow H1/2 at those endpoints, so a resolution limit stops
inference before any demonstrated critical-norm turnover. Their actual
evolved checkpoints are saved, along with an independent leading-field
`N=64, T=.10` replay and a separately computed critical production–diffusion
budget for designing the next search objective.

The subsequent [64/128 check through T=.10](frozen_continuations_resolution/README.md)
passes the preliminary continuation gates for the leading frozen field.
H1/2 grows 9.91%, sampled vorticity reaches about 2.02 times its initial
value, and the vorticity discrepancy decreases from 16.23% to 4.63% after
raising the evolution resolution. A new `N=64` timestep refinement and
independent trajectory pass. The final 128-grid state is archived and can
be restored without rerunning the beginning. L3 still falls, and the
remaining pointwise uncertainty is not a rigorous continuum error bound.

A proof would additionally require rigorous control of the unresolved
Fourier tail and of the limiting PDE mechanism; none of these finite runs
supplies that control.
