# Frozen late-growth continuation from .10 to .12

This bounded experiment continues the unchanged bandwidth-three initial field
from its independently checked, published `.10` evolved states on grids `64`
and `128`. The coefficients are not reoptimized on these held-out times.

**Status: both legs and all declared continuation and endpoint sampling gates
pass through `.12`.** The fine extension takes 534 new accepted steps, about
22 minutes of solver time; the coarse extension takes 261 new steps. All four
new evolved states are archived and restore to their original hashes.

The outcome remains finite amplification: H1/2 grows `16.65492%` by `.12`,
but L3 falls `3.77084%`, and relative H1/2 growth keeps slowing. This is not
a singularity or regularity proof.

## Completed measurements

All ratios use the original initial field. L3 and vorticity in this table
use the common physical sampling grid `128`. Cutoff entries are energy
fractions, not percentages; step counts include the preceding trajectory.

| Quantity | T=.11, N=64 | T=.11, N=128 | T=.12, N=64 | T=.12, N=128 |
|---|---:|---:|---:|---:|
| H1/2 / initial | 1.156664301 | 1.156667545 | 1.166558012 | 1.166549211 |
| L3 / initial | .969492876 | .969569973 | .962174724 | .962291572 |
| Sampled vorticity / initial | 1.986598862 | 2.140456632 | 2.239796270 | 2.283198617 |
| Enstrophy / initial | 2.172313677 | 2.173971230 | 2.324911447 | 2.328371906 |
| Characteristic wavenumber / initial | 1.500170121 | 1.500744653 | 1.555799118 | 1.556963262 |
| Critical logarithmic growth rate | .909765402 | .909564528 | .794072353 | .791847344 |
| Critical production / viscous destruction | 3.050139386 | 3.032888773 | 2.547326514 | 2.518860450 |
| Enstrophy stretching / viscous destruction | 4.657278429 | 4.532300968 | 3.828670482 | 3.648347290 |
| Peak per-step cutoff energy fraction | 8.226740e-5 | 1.039116e-7 | 1.350872e-4 | 4.453897e-7 |
| Accepted steps | 1,005 | 1,824 | 1,141 | 2,104 |

The fine norm gains another `1.83564%` over `.10--.12`, above the retained
`0.1%` late-window growth gate, with no sampled late reversal. Its relative
growth rate falls from `1.026857962` at `.10` to `.791847344` at `.12`.
The critical production-to-diffusion ratio falls from `3.674157991` to
`2.518860450`: net critical production is still positive at these samples,
but its margin over viscous loss is shrinking. This budget concerns H1/2
squared and is distinct from the enstrophy stretching budget.

## Resolution gaps versus the unchanged screening limits

The following are maximum relative gaps at common observations through the
specified horizon. Stretching uses the last `.02` window. Percentage gaps
are `100 * abs(a-b)/max(abs(a),abs(b))`.

| Check | Through .11 | Through .12 | Limit |
|---|---:|---:|---:|
| H1/2 ratio gap | .00028049% | .00075441% | 2% |
| L3 ratio gap | .00795167% | .01214272% | 2% |
| Enstrophy ratio gap | .07624540% | .14862144% | 5% |
| Characteristic wavenumber ratio gap | .03828311% | .07477016% | 5% |
| Sampled vorticity ratio gap | 7.18808160% | 7.18808160% | 10% |
| Late stretching gap | 1.76110337% | 2.93035970% | 10% |
| Horizon times maximum critical-rate gap | 2.209610e-5 | 2.670010e-4 | .002 |

The worst common-sampling vorticity gap occurs at `.11`; it is not the
endpoint gap at `.12`. An endpoint that agrees more closely cannot erase an
earlier discrepancy. Conversely, nearly identical critical norms do not
guarantee equally accurate local vorticity peaks. The largest native-to-common
sampling shift over the joined trajectories remains `1.345230%`, below `2%`.

These are numerical screening limits, not target values for a proof.
[`assessment.json`](assessment.json) records all gates and their outcomes.

## Endpoint physical sampling

The initial sampled maxima are `82.55311350` on `128` samples and
`82.60944619` on `256` samples, shared by both evolution grids.

| Evolved field | Vorticity ratio, 128 samples | Vorticity ratio, 256 samples | Raw maximum shift | Ratio shift |
|---|---:|---:|---:|---:|
| N=64, T=.11 | 1.986598862 | 1.993557108 | .416990% | .349037% |
| N=128, T=.11 | 2.140456632 | 2.141273167 | .106299% | .038133% |
| N=64, T=.12 | 2.239796270 | 2.241126209 | .127494% | .059342% |
| N=128, T=.12 | 2.283198617 | 2.289872340 | .359438% | .291445% |

Every raw and normalized shift is below `2%`. On `256` physical samples,
the cross-evolution vorticity-ratio gap is `6.898515%` at `.11` and `2.128771%`
at `.12`, both below `10%`. The final fine-grid sampled maximum is
`189.165085824`, or `2.289872340` times its initial sampled value. Denser
sampling evaluates the same retained Fourier field; it is not a `256`-grid
evolution or a certificate for its continuous maximum.

## Frozen protocol and provenance

[`protocol.json`](protocol.json) was fixed before the runs. It keeps initial
energy `10`, viscosity `.02`, maximum timestep `.000125`, CFL target `.4`,
observation interval `.01` and per-step cutoff-energy limit `.008`. Retained
Fourier cutoffs are `21` and `42`. Both grids use `128` physical samples for
their common-sampling diagnostics; native sampling remains `64` or `128`.

The initial coefficient SHA-256 remains
`893ef61450d670936dc851d28ab3cf7c1aa315c6ff23934646e74137c40ac6ef`.
The parent checkpoint is commit `96829fb2cd74d3c39699a6866fb601ed81596843`.
Each grid restores its own complete evolved `.10` field, with its original
normalization, configuration, accepted-step count, observation clock and
accumulated diagnostics. No evolved coarse field is padded into a finer grid.

The executable and all previously recorded solver sources were verified
unchanged. The `64` source is a published whole-gzip checkpoint; the `128`
source is a published multipart gzip archive. Both restore to the original
raw-state hashes and pass independent NumPy critical-budget checks before
evolution starts. Restart scientific settings cannot be overridden by the CLI.

## Measurements and limits

The earlier `assess.py:paired` implementation is reused without changing its
thresholds. It evaluates the complete joined histories through `.11` and
`.12`, retaining the original late `.02` window and critical-budget samples
over the last half of each horizon. The existing `2%` native-to-common
vorticity sampling gate is retained as well as the cross-evolution gates.

Independent NumPy FFTs check the new critical budgets and sample both saved
endpoint fields on physical grids `128` and `256`. Raw vorticity maxima and
their initial-normalized ratios must each shift by at most `2%`; the
cross-evolution vorticity-ratio gap on `256` samples must be at most `10%`.
These endpoint checks do not certify the continuous spatial maximum or
values between observation times.

The earlier full independent-solver and smaller-timestep tests reach `.10`
on `N=64`. This extension does not claim those full-trajectory tests beyond
`.10`, nor a new independent `N=128` evolution. Its additional independent
checks concern restored and newly saved snapshots and their physical samples.

All thresholds are finite-amplification screening criteria. A failed
resolution gate pauses promotion of this grid pair; it does not rule out the
initial field as a possible mechanism. A passing gate does not establish an
infinite-dimensional PDE error bound or a blow-up proof.

The next bounded time target is `.14`, retaining the same initial field and
the existing gates. Its purpose would be to determine whether late critical
growth persists or turns over. The present evidence does not justify
extrapolating the measured vorticity or growth rates toward infinity. Additional
independent-solver and timestep validation is needed for stronger claims about
the new interval than these preliminary screening results.

## Reproduction and archive verification

From the repository root, use fresh output directories:

```sh
python3 research/navier_stokes_cascade/results/late_growth_t012/extend.py \
  --solver ./build/research/navier_stokes_cascade/navier_stokes_continue \
  --grid 64 --output-dir new-t012-run/reference64

python3 research/navier_stokes_cascade/results/late_growth_t012/extend.py \
  --solver ./build/research/navier_stokes_cascade/navier_stokes_continue \
  --grid 128 --output-dir new-t012-run/reference128
```

The recorded driver intentionally requires the previously verified executable
hash. An independently rebuilt binary must be separately validated and its
new provenance recorded; it cannot silently replace the archived executable.
In a fresh environment, establish a separately verified preceding campaign
with that build before declaring a corresponding new extension protocol.
The fresh directories protect existing evidence. Output manifests record
exact commands, source and executable hashes, restart links, run identities,
independent budgets and saved-state identities.

`analyze.py` operates on this archive's `reference64` and `reference128`
directories, joining them to the explicitly identified earlier histories.
`archive.py` can pack already completed targets while another leg runs; that
does not mark the overall experiment finished. Original gzip outputs are
ignored local copies; published parts and indices are sufficient to restore
the full states. `audit.py` reports archive integrity separately from the
scientific pass or failure, and `test_tools.py` checks rejection of malformed
clocks, failed scientific gates and corrupted archive parts.

All 15 helper regression checks pass; see [`verification.log`](verification.log).
The restored `.10` fields and all four new snapshots pass the independent
`1e-10` relative diagnostic check. The final [`audit.json`](audit.json) verifies
source and restart provenance, assessment arithmetic, and restoration of all
four states. The two coarse states each restore `12,583,408` bytes; the two
fine states each restore `100,663,792` bytes. Archive integrity and the
scientific outcome are reported separately.

Restore a saved fine-grid state to a fresh file using the existing reader:

```sh
python3 research/navier_stokes_cascade/results/frozen_continuations_resolution/restore_checkpoint.py \
  --index research/navier_stokes_cascade/results/late_growth_t012/reference128/t120_checkpoint.json \
  --output fine128-t120-restored.chk
```

Use `reference64` for the coarse field or `t110_checkpoint.json` for `.11`.
The output file must not already exist. Checkpoints are reproduction inputs,
not proof certificates. No trajectory beyond `.12` is part of this archive.
