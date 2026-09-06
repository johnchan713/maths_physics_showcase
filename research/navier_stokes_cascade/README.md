# Navier-Stokes cascade research scaffold

This directory is an **experimental research scaffold, not a claimed proof or
disproof** of the Navier-Stokes Millennium Prize problem. It gives the
repository a reproducible place to test one concrete direction: whether
finite-stage transfer of energy toward high frequencies can be organized into
an accumulating cascade, or whether viscosity and geometric depletion prevent
that transfer.

## Exact scope

The Clay problem concerns smooth, divergence-free initial data for the
three-dimensional incompressible equations

```text
partial_t u + (u . grad)u = -grad p + nu Delta u,
div u = 0,
u(0) = u_0,
```

on R^3 or the periodic three-torus, with positive viscosity. A resolution must
either prove global smoothness for every allowed datum or rigorously construct
one allowed smooth datum whose solution breaks down in finite time. Weak
non-uniqueness from singular data, blow-up for a modified/averaged equation,
or a large finite numerical value does not settle that statement.

## What is implemented

Two independent numerical paths now evolve a mean-zero, real, divergence-free
velocity field on the 2-pi periodic torus using the Fourier cube
`|k_x|, |k_y|, |k_z| <= K`:

- `GalerkinSystem` evaluates every Fourier triad directly. Its quadratic cost
  makes it slow, but its simple formula is the low-cutoff correctness oracle.
- `PseudospectralSystem` uses an in-repository radix-2 three-dimensional FFT,
  the rotational nonlinearity `P[u x curl(u)]`, and strict 2/3 de-aliasing. Its
  cost scales approximately as `N^3 log(N)` and permits materially larger runs.

For each retained non-zero mode, both backends compute

```text
d u_k / dt = -nu |k|^2 u_k
             - i P_k sum_(p+q=k) (u_p . q) u_q,
```

where `P_k` is the Leray projector. Time stepping is classical RK4. The FFT
path can choose each step from the conservative Fourier bound

```text
U_F = sum_k |u_k| >= ||u||_infinity,
dt  = min(dt_max, CFL/(sqrt(3) K U_F), S_nu/(3 nu K^2)).
```

The defaults are `CFL=0.4` for a single FFT run (`0.35` in the search) and
`S_nu=2`, inside the negative-real-axis stability interval of RK4. The final
step is shortened to end at the same physical time at every resolution. The
CSV records every sampled accepted step and its conservative advective and
viscous stability numbers. The direct convolution is intentionally small and
auditable; it is not intended as a production turbulence solver.

Long FFT runs can atomically replace a checksummed, versioned binary
checkpoint. It contains the full complex Fourier state, scientific
configuration, physical time, accepted-step count, accumulated diagnostics,
and spectrum-profile references. Restart restores those values rather than
reconstructing the initial condition. A split fixed-step run and an
uninterrupted run are required by the tests and CI to produce bit-for-bit
identical final checkpoints. Each scheduled checkpoint is also a diagnostic
step, so its accumulated telemetry and profile reference correspond to the
stored state time. The checksum detects accidental corruption; it is not a
proof certificate or a guarantee against storage failure.

An optional test target uses FFTW3 as an implementation-independent transform
oracle. It compares forward and inverse three-dimensional transforms and
reconstructs the complete dealiased Navier-Stokes right-hand side through an
FFTW path before comparing every retained coefficient. This is a useful check
against a shared FFT bug, although it is not yet a second full time-evolution
codebase.

For an `N^3` FFT grid, the code enforces

```text
K <= floor((N - 1) / 3),
```

so `3K < N`. Therefore a wrapped quadratic product cannot alias back into the
retained cube. Modes outside the retained cube and the mean mode are never
evolved. The Fourier normalization is `u(x) = sum_k u_k exp(i k.x)`.

Four smooth, exactly divergence-free starting fields are available:

- `taylor-green` (default), a standard vortex-interaction benchmark;
- `abc`, a Beltrami flow whose projected nonlinear term vanishes, used here as
  a negative control: a cascade detector should report no nonlinear flux;
- `deterministic`, the original reproducible low-mode mixture;
- `vortex-tubes` (FFT only), a parameterized counter-rotating pair with core
  radius, separation, helical bend amplitude, and axial wavenumber.

The Taylor-Green field is

```text
u = (sin(x) cos(y) cos(z), -cos(x) sin(y) cos(z), 0),
```

and the equal-parameter ABC field is

```text
u = (sin(z) + cos(y), sin(x) + cos(z), sin(y) + cos(x)).
```

Both are normalized to the requested kinetic energy without changing their
shape.

The vortex pair starts from the smooth periodic vector potential

```text
A = (0, 0, G_1 - G_2),
G_j = exp(-(2(1-cos(x-c_jx)) + 2(1-cos(y-c_jy)))/(2 a^2)),
u = curl(A).
```

The two centres follow oppositely displaced helices in `z`. Setting the bend
to zero gives a z-invariant two-dimensional control; a non-zero bend populates
three-dimensional Fourier modes. The sampled field is truncated at the safe
cutoff, Leray projected, and energy normalized. This is a reproducible smooth
test family, not a claim that it resembles a singular profile.

The CSV diagnostics include:

- normalized kinetic energy, enstrophy, and palinstrophy;
- the scale-critical homogeneous H1/2 Fourier norm and a sampled
  scale-critical L3 velocity norm;
- sampled maximum vorticity and its time integral (a BKM-inspired diagnostic);
- a rigorous-for-the-truncated-polynomial Fourier upper bound on maximum
  vorticity;
- spectral centroid and energy fraction touching the cutoff shell;
- energy, nonlinear transfer, viscous loss, and forward flux for every radial
  Fourier shell;
- nonlinear enstrophy production, viscous enstrophy destruction, their net
  rate, and their ratio;
- a normalized rescaled energy spectrum, its characteristic wavenumber,
  profile distances, and a heuristic exponential-tail fit;
- divergence, Fourier-reality, and semi-discrete energy-balance defects.

With the repository's Fourier normalization, the critical Sobolev diagnostic
and enstrophy budget are

```text
||u||_(Hdot 1/2) = (sum_k |k| |u_k|^2)^(1/2),
Omega              = (1/2) sum_k |k|^2 |u_k|^2,
d Omega / dt       = S - 2 nu P,
S                  = sum_k |k|^2 Re(conj(u_k) . N_k),
P                  = (1/2) sum_k |k|^4 |u_k|^2.
```

Thus `S/(2 nu P) > 1` means nonlinear vortex stretching is instantaneously
creating enstrophy faster than viscosity destroys it. For `nu=0`, the ratio is
reported as zero because its denominator vanishes. Neither a ratio above one
nor finite growth of a critical norm implies blow-up.

The rescaled-spectrum diagnostic uses

```text
k_rms = sqrt(enstrophy / energy),
xi    = |k| / k_rms.
```

Each mode's fraction of total energy is deposited continuously onto
neighbouring `xi`-bin centres, with a separate overflow node. For two samples
the code reports their L1 distance, overlap `1-L1/2`, absolute change in
`log(k_rms)`, and

```text
shape drift = L1 distance / |change in log(k_rms)|.
```

This normalization prevents an artificially short final sampling interval
from looking like profile convergence. An exactly preserved rescaled shape
has zero drift; a credible self-similar candidate should make the drift
decrease toward zero as resolution grows. The command-line threshold
`shape drift <= 1` is only a deliberately conservative triage heuristic, not
a theorem.

The tail diagnostic fits `log E(k)` linearly over the upper half of the
isotropic retained shells and reports `delta=-slope/2`, the fit point count,
R-squared, and the rescaled product `delta*k_rms`. An exactly self-similar
exponential tail would keep that product constant. The fit is motivated by
`E(k) ~ k^a exp(-2 delta k)` but ignores the algebraic prefactor and provides
no rigorous lower bound on analyticity radius.

For the unit-width shell `S_j = {k : j-1 < |k| <= j}`, the code records

```text
E_j  = (1/2) sum_(k in S_j) |u_k|^2,
T_j  = sum_(k in S_j) Re(conj(u_k) . N_k),
D_j  = nu sum_(k in S_j) |k|^2 |u_k|^2,
Pi_j = -sum_(m <= j) T_m.
```

`Pi_j > 0` means the modes inside radius `j` are losing energy through the
nonlinear term to finer modes. The identities `sum_j T_j = 0`,
`sum_j D_j = 2 nu * enstrophy`, and `Pi_last = 0` are tested. Positive flux is
necessary for a forward cascade; it is nowhere near sufficient for blow-up.

The tests check the structural identities before any experiment is interpreted:
Leray projection, real/divergence-free invariants, nonlinear energy
cancellation, shell accounting, viscous energy dissipation, the exact
semi-discrete enstrophy budget, an analytic Taylor-Green H1/2 value, the ABC
negative control, and short-run numerical stability.
The FFT-specific suite also checks a complex 3D transform round trip, rejects
unsafe cutoffs, compares every nonlinear Fourier coefficient against direct
convolution, repeats that comparison with all modes populated at the maximum
safe cutoff, and compares complete short trajectories. It now also verifies
that straight tubes have no non-zero axial modes, bent tubes do, both remain
real and divergence-free, invalid geometry is rejected, and adaptive steps
obey both requested stability bounds. The profile tests check normalization,
amplitude invariance, an exact discrete scale shift with zero shape drift, and
metric ranges. Checkpoint
tests cover exact round trips, split/uninterrupted trajectory identity, and
checksum-corruption rejection. When FFTW3 is present, a separate target checks
the internal transforms and the full nonlinear right-hand side against FFTW.
The candidate-score suite separately checks that lower profile drift and lower
cutoff loading improve a score, that both grid levels contribute to the paired
cutoff cost, and that missing, worsening, or cross-resolution-inconsistent
profile evidence cannot pass the refinement gate.
They also encode the elementary guardrail that zero gradient does not imply
zero field value, a mistake found in some purported proofs.

## Build and run

From the repository root:

```bash
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/research/navier_stokes_cascade/navier_stokes_cascade
```

To make the independent transform oracle mandatory, install FFTW3 and add
`-DNS_CASCADE_REQUIRE_FFTW_REFERENCE=ON` to the configure command. The branch
workflow does this on Ubuntu so a missing reference library fails CI rather
than silently skipping the check.

Run the faster FFT backend on a `32^3` grid, retaining the strict safe cutoff
`K=10`:

```bash
./build/research/navier_stokes_cascade/navier_stokes_fft \
  --grid 32 --cutoff 0 --initial-condition taylor-green \
  --dt 0.005 --final-time 0.1 --cfl 0.4 --diagnostic-every 20 \
  --output fft-n32.csv --shell-output fft-n32-shells.csv \
  --profile-output fft-n32-profile.csv
```

Here `--dt` is the maximum allowed adaptive step. Add `--fixed-dt` to recover
fixed-step evolution; if `--final-time` is omitted, the end time is
`--steps * --dt`. The executable writes a time-series CSV and a separate
long-form shell CSV, and a long-form rescaled-spectrum CSV.

Run one member of the smooth vortex-tube family with:

```bash
./build/research/navier_stokes_cascade/navier_stokes_fft \
  --grid 32 --initial-condition vortex-tubes \
  --tube-core 0.7 --tube-separation 1.8 \
  --tube-bend 0.3 --tube-axial-mode 2 \
  --viscosity 0.02 --dt 0.005 --final-time 0.1 \
  --profile-bins 64 --profile-max-xi 4 \
  --output tube.csv --shell-output tube-shells.csv \
  --profile-output tube-profile.csv \
  --checkpoint-output tube.chk --checkpoint-every 100
```

Resume the exact stored state to a new absolute final time with:

```bash
./build/research/navier_stokes_cascade/navier_stokes_fft \
  --restart tube.chk --final-time 0.2 --diagnostic-every 20 \
  --output tube-part-2.csv --shell-output tube-part-2-shells.csv \
  --profile-output tube-part-2-profile.csv \
  --checkpoint-output tube.chk --checkpoint-every 100
```

Scientific options cannot be mixed with `--restart`; grid, cutoff, viscosity,
initial data, timestep policy, and profile binning all come from the
checkpoint. Output paths and the new absolute final time remain selectable.

For a slightly larger exploratory run:

```bash
./build/research/navier_stokes_cascade/navier_stokes_cascade \
  --initial-condition taylor-green --cutoff 3 --viscosity 0.02 \
  --dt 0.0002 --steps 1000 --diagnostic-every 25 \
  --output cascade-k3.csv --shell-output cascade-k3-shells.csv
```

Run the same physical experiment at several cutoffs and at both `dt` and
`dt/2` with one command:

```bash
./build/research/navier_stokes_cascade/navier_stokes_convergence \
  --initial-condition taylor-green --cutoffs 2,3 \
  --dt 0.0005 --steps 60 --diagnostic-every 20 \
  --output navier_stokes_convergence.csv
```

Use the same comparison program with the FFT backend:

```bash
./build/research/navier_stokes_cascade/navier_stokes_convergence \
  --backend fft --grids 16,32 --initial-condition taylor-green \
  --dt 0.0005 --steps 60 --diagnostic-every 20 \
  --output navier_stokes_fft_convergence.csv
```

Search a small parameter grid at `16^3`, rank only after applying cutoff and
constraint gates, and rerun the top three candidates at `32^3`:

```bash
./build/research/navier_stokes_cascade/navier_stokes_search \
  --coarse-grid 16 --fine-grid 32 \
  --cores 0.55,0.70 --separations 1.2,1.8 \
  --bends 0,0.30 --axial-modes 1,2 \
  --energies 1,4,10 --viscosity 0.02 \
  --dt 0.005 --final-time 0.1 --top 3 \
  --profile-bins 64 --profile-scale-window 0.025 \
  --profile-drift-threshold 1 \
  --output navier_stokes_candidate_search.csv
```

The straight cases are deduplicated because their axial wavenumber has no
effect. `--energy E` remains the single-energy shorthand; `--energies` includes
energy in the candidate grid. The CSV records initial/final/peak critical L3
and H1/2, sampled vorticity, enstrophy, palinstrophy, enstrophy
production/destruction, maximum positive shell flux, cutoff contamination,
constraint defects, accepted-step statistics, fixed-scale profile drift,
score components, and coarse/fine differences. Profile drift is sampled
whenever `log(k_rms)` advances by `--profile-scale-window`, independently of
`--diagnostic-every`; changing CSV verbosity therefore cannot change the scale
windows. Only forward movement to finer scales completes a window. The cutoff
fraction used by the resolution gate and score is also checked after every
accepted step, so a short contamination spike cannot hide between CSV rows.

For one resolution, let `G3`, `Gh`, and `Gw` be the peak L3, H1/2, and sampled
vorticity ratios; `F3` the final L3 ratio; `C/C*` the peak cutoff fraction
normalized by its threshold; and `D/D*` the latest profile drift normalized by
its threshold. The auditable single-grid score is

```text
S_N = max(log G3, log Gh) + 0.25 min(log G3, log Gh)
      + 0.02 log Gw + 0.01 log F3
      + 0.05 clamp(log(D_first/D_latest), -2, 2)
      - 0.05 clamp(log(D_latest/D_min), 0, 2)
      - 0.10 log(1 + D/D*) - 0.03 log(1 + C/C*).
```

The trend and rebound terms are zero until two scale windows exist. The
rebound term prevents an early low drift followed by a deteriorating final
window from looking stationary. No completed window is assigned the finite
missing-evidence cost `D/D*=4`, rather than being mistaken for zero drift.
Coarse finalists are rerun on the fine grid and receive

```text
C_pair = sqrt(((C_coarse/C*)^2 + (C_fine/C*)^2) / 2),
S_pair = min(S_coarse, S_fine) - 0.04 log(1 + C_pair)
         - 0.25 (relative L3 difference + relative H1/2 difference)
         - 0.10 relative k_rms-growth difference
         - 0.05 relative profile-drift difference.
```

Using the weaker single-grid score prevents one impressive resolution from
hiding its companion. A row becomes `refinement_eligible` only when it passes
the existing convergence gate, the same critical norm grows materially on
both grids, `k_rms` grows by at least 10% on both, each grid has at least two
profile windows whose latest drift is no larger than its first or than one,
no more than 10% above the best previous drift, each latest window occurs in
the final quarter of the run, the two drifts agree within 50%, and scale
growth agrees within 10%. These
thresholds and weights are deterministic triage choices, not probabilities or
mathematical implications.

The direct comparison uses one common spatial sampling grid for every cutoff.
The FFT comparison samples on each native grid; therefore its L3-grid
difference includes quadrature error as well as solution error. The H1/2
diagnostic is computed directly from Fourier coefficients. The CSV contains
total and shell energies, both critical-norm ratios, the enstrophy budget,
nonlinear transfer, viscous dissipation, forward flux, constraint defects, and
cutoff-shell checks.
The terminal summary reports both timestep differences and adjacent-resolution
differences, including common-shell flux disagreement.

### Verified reference run

The default Taylor-Green comparison (`nu=0.05`, `dt=0.0005`, final time
`0.03`) produced this small-cutoff baseline during implementation:

| Check | Result |
|---|---:|
| K=2 peak cutoff-shell energy fraction | 4.44501e-4 |
| K=3 peak cutoff-shell energy fraction | 3.72479e-7 |
| K=2 versus K=3 final-energy difference | 1.11183e-9 |
| K=2 versus K=3 final-L3 difference | 3.83754e-5 |
| K=2 versus K=3 common-shell flux difference | 1.67587e-3 |
| Maximum CSV shell-energy sum residual | 5.551e-16 |
| Maximum total nonlinear-transfer residual | 3.614e-18 |

The `dt` versus `dt/2` differences were at floating-point round-off for this
short run. These figures verify the implementation baseline only. They do not
show growing critical norms, an infinite cascade, or singular behaviour; so
far, the monster has merely submitted tidy paperwork.

The corresponding FFT comparison at `N=16, K=5` and `N=32, K=10` gave:

| Check | Result |
|---|---:|
| N=16 peak cutoff-shell energy fraction | 5.46049e-13 |
| N=32 peak cutoff-shell energy fraction | 1.17288e-26 |
| Final-energy resolution difference | 3.36079e-16 |
| Final-L3 resolution difference | 1.40368e-4 |
| Common-shell flux difference | 7.35405e-12 |
| N=32 final critical-L3 ratio | 0.995387 |

The run has a clean, converged finite forward flux, but its critical L3 norm
decreases by about 0.46%. That is evidence the measurement machinery works;
it is not a blow-up candidate.

### Vortex-tube search result

An extended 12-configuration sweep at `nu=0.02`, initial energy one, and
`t=0.1` found no critical-L3 growth. The only finalist whose growth/no-growth
signals and flux passed the preliminary resolution gate had
`core=0.7`, `separation=1.8`, `bend=0.3`, and `axial mode=2`. A dedicated
`32^3 -> 64^3` rerun produced:

| Check at t=0.1 | N=32, K=10 | N=64, K=21 |
|---|---:|---:|
| Peak critical-L3 / initial | 1.000000 | 1.000000 |
| Final critical-L3 / initial | 0.988878748 | 0.988879563 |
| Peak critical-H1/2 / initial | 1.000000 | 1.000000 |
| Final critical-H1/2 / initial | 0.989693218 | 0.989693227 |
| Peak sampled vorticity / initial | 1.03280068 | 1.04312237 |
| Final enstrophy / initial | 0.981387172 | 0.981387447 |
| Final nonlinear enstrophy production | 1.22370216 | 1.22377821 |
| Final production / viscous destruction | 0.809589276 | 0.809629929 |
| Maximum positive forward flux | 0.050832943 | 0.050832879 |
| Peak cutoff-shell energy fraction | 4.06726e-7 | 4.56699e-14 |
| Accepted adaptive steps | 59 | 124 |
| Maximum conservative CFL bound | 0.40 | 0.40 |

The final-L3 ratios agree to about `8.2e-7` relative and the peak forward flux
to about `1.3e-6` relative. The final H1/2 ratios differ by less than `1e-8`
relative, and the production/destruction ratios by about `5e-5` relative. The
sampled-vorticity growth differs by about 0.99% relative, partly because the
maximum is sampled on each native grid.

Extending the same candidate to `t=0.5` at `N=32` gave peak sampled-vorticity
growth `1.12599`, enstrophy growth `1.09518`, final critical-L3 ratio
`0.932052`, final critical-H1/2 ratio `0.972161`, and cutoff fraction
`4.96984e-4`. Net enstrophy production first became positive at the sampled
time `t=0.1264`, but neither critical norm reversed its decline. The
corresponding `N=16` run was under-resolved (cutoff fraction `0.0200`), so the
long-time enstrophy growth still needs an `N=64` check.

The defensible interpretation is narrow: this family exhibits resolved local
vortex amplification and forward transfer over the short interval, while its
scale-critical L3 norm decreases. That is a useful configuration for studying
vortex stretching or geometric depletion, but it is presently evidence
against this particular run being a blow-up candidate—not a proof of global
regularity and not a disproof of Navier-Stokes.

### Critical-H1/2 stress test

Increasing the normalized initial energy changes the competition between the
nonlinear and viscous terms without changing the broad-tube geometry. A
`16^3 -> 32^3` sweep for `core=0.7`, `separation=1.8`, `bend=0.3`, `axial
mode=2`, `nu=0.02`, and `t=0.1` produced the following fine-grid values; all
six coarse/fine pairs passed the preliminary cross-resolution gate:

| Initial energy | Final H1/2 / initial | Final L3 / initial | Peak enstrophy / initial | Max production / destruction | Peak cutoff fraction |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.989693 | 0.988879 | 1.000000 | 0.809589 | 4.07e-7 |
| 2 | 0.991650 | 0.988170 | 1.000000 | 1.618339 | 1.83e-6 |
| 4 | 0.995715 | 0.986683 | 1.022627 | 3.068124 | 1.09e-5 |
| 6 | 0.999876 | 0.985121 | 1.052683 | 4.214082 | 3.41e-5 |
| 8 | 1.004060 | 0.983527 | 1.084048 | 5.068076 | 7.95e-5 |
| 10 | 1.008226 | 0.981948 | 1.116363 | 5.694614 | 1.54e-4 |

The transition from finite-time H1/2 decay to growth lies between energies six
and eight for this geometry and observation time. Notice that the enstrophy
budget becomes production-dominated before H1/2 grows; `production /
destruction > 1` is not itself a critical-norm criterion. At energy ten, an
`N=64, K=21` rerun gave H1/2 growth `1.008248`, final L3 ratio `0.981963`, peak
enstrophy ratio `1.117031`, and cutoff fraction `1.53e-7`. Its H1/2 growth ratio
differs from `N=32` by about `2.2e-5` relative.

This is a reproducible finite critical-norm growth signal, not a singularity
signal. It is less than one percent, ends at a fixed finite time, coexists with
decreasing L3, and comes from smooth finite-dimensional systems at every
resolution.

A sharper profile, first flagged but not interpreted on the under-resolved
coarse grid, strengthened the signal. For `core=0.55`, `separation=1.2`,
`bend=0.3`, `axial mode=2`, energy ten, `nu=0.02`, and `t=0.1`:

| Check | N=32, CFL 0.4 | N=32, CFL 0.2 | N=64, CFL 0.4 |
|---|---:|---:|---:|
| Peak/final H1/2 / initial | 1.042902961 | 1.042902961 | 1.044373144 |
| Peak L3 / initial | 1.000000000 | 1.000000000 | 1.000000000 |
| Final L3 / initial | 0.962013680 | 0.962013680 | 0.961997639 |
| Peak/final enstrophy / initial | 1.409766731 | 1.409766731 | 1.450005700 |
| Peak sampled vorticity / initial | 1.396989578 | 1.396997100 | 1.674481303 |
| Maximum production / destruction | 6.755828833 | 6.755906152 | 6.742804414 |
| Maximum positive forward flux | 12.698342 | 12.698350 | 12.582115 |
| Peak cutoff-shell energy fraction | 5.02010e-3 | 5.02010e-3 | 7.68510e-5 |
| Accepted adaptive steps | 357 | 714 | 817 |

Halving the CFL target changes the final H1/2 value by only about `3e-13`
relative. From `N=32` to `N=64`, the H1/2 growth ratios differ by about
`1.4e-3` relative, the final L3 ratios by about `1.7e-5`, and the maximum
forward flux by about `9.2e-3`. The enstrophy ratios still differ by roughly
2.8%, while the native-grid sampled-vorticity ratios differ much more. The
critical norm and flux are substantially better converged than the pointwise
or more heavily derivative-weighted observables.

At both resolutions H1/2 first dips by about 0.14%, bottoms near `t=0.014`,
crosses its initial value near `t=0.0275`, and reaches its sampled maximum at
the final time. This is the strongest clue found by the scaffold so far: a
converged net finite-time increase of one scale-critical norm. It is still far
from a blow-up construction. The other monitored critical norm decreases,
H1/2 has not shown divergence or a stable rescaled profile, the observation
stops while H1/2 is still rising, and no unresolved Fourier-tail bound exists.

Checkpoint/restart then extended the `N=64` candidate from `t=0.1` to
`t=0.2`. The complete run used 2,457 accepted adaptive steps. A companion
`N=32` run crossed the one-percent cutoff gate near `t=0.1288` and ended with
2.86% cutoff loading, so its later values must not be interpreted. At the same
failure time the `N=64` cutoff fraction was only about 0.026%, and it remained
below the gate through `t=0.2`:

| Check through t=0.2 | N=32, K=10 | N=64, K=21 |
|---|---:|---:|
| Peak/final H1/2 / initial | 1.094806 | 1.102222 |
| Final L3 / initial | 0.917463 | 0.915652 |
| Final enstrophy / initial | 2.128676 | 2.525304 |
| Peak sampled vorticity / initial | 1.797454 | 4.009082 |
| Final characteristic wavenumber / initial | 1.533645 | 1.677045 |
| Peak cutoff-shell energy fraction | 2.86343e-2 | 1.85488e-3 |
| Final tail-fit delta | 0.137844 | 0.100243 |
| Resolution verdict | failed after t about 0.1288 | passed cutoff gate |

Only the `N=64` column passes the necessary cutoff criterion at the final time,
so the difference between the two columns is not a convergence error estimate
and no post-`t=0.1288` cross-resolution claim is made. On `N=64`, H1/2 was
still increasing at the final sample, but its incremental growth was slowing.
Sampled vorticity peaked at about four times its initial value near `t=0.1871`
and then fell to about 3.24 times initial by `t=0.2`. L3 decreased by about
8.43%.

Most importantly, the final rescaled-profile drift was about `10.5` per unit
change in `log(k_rms)`, well above the heuristic threshold one. The spectrum
therefore moved to finer scales without settling to a stationary rescaled
shape. The originally tempting raw distance to the immediately preceding row
was rejected because that row was separated by only seven final steps. The
normalized diagnostic closes that sampling-interval loophole. This is a
spectrally clean finite cascade episode on `N=64`, not a self-similar blow-up
profile converged across resolutions; the simple tail-fit `delta*k_rms`
product also fell from about 3.11 initially to 0.481 rather than remaining
constant. That regression is heuristic and does not alter the verdict by
itself.

### Profile-aware refinement result

The first profile-aware neighbourhood sweep varied `core=0.50,0.55,0.60`,
`separation=1.1,1.3`, and `bend=0.25,0.35` at axial mode two, energy ten,
`nu=0.02`, and `t=0.08`. All 12 candidates passed the `N=32` cutoff and
constraint gates. The paired score selected `core=0.50`, `separation=1.30`,
and `bend=0.35` for an `N=64` rerun:

| Check through t=0.08 | N=32, K=10 | N=64, K=21 |
|---|---:|---:|
| Peak H1/2 / initial | 1.044539613 | 1.046122325 |
| Final L3 / initial | 0.963589780 | 0.963279662 |
| Peak enstrophy / initial | 1.412584781 | 1.453402602 |
| Peak sampled vorticity / initial | 1.377488837 | 1.507441512 |
| Final characteristic wavenumber / initial | 1.210592209 | 1.228108183 |
| Completed 2.5% log-scale windows | 7 | 8 |
| First profile drift | 16.641357 | 16.635310 |
| Minimum profile drift | 10.980775 | 11.040019 |
| Latest profile drift | 10.980775 | 13.403208 |
| Peak cutoff-shell energy fraction | 8.00816e-3 | 7.91022e-5 |

The H1/2 growth ratios differ by about 0.15%, characteristic-scale growth by
about 1.43%, and latest drifts by about 18.1%; the pair passed the preliminary
cross-resolution gate. The fine run therefore confirms finite critical-norm
growth and spectral motion with very low cutoff loading. It does not confirm
profile stationarity: both latest drifts exceed the threshold one by more than
an order of magnitude, and the fine-grid drift rebounded after reaching its
minimum. Consequently zero candidates passed the strict refinement gate and
no longer run was promoted. This rejects this small neighbourhood at this
time horizon under the stated heuristic; it does not exclude other initial
data, parameters, or later behaviour.

Use `--help` for all parameters. The direct backend still grows quadratically
in the retained mode count; use it to audit small cases and the FFT backend to
explore larger ones.

The branch-scoped GitHub Actions workflow builds these CMake targets, runs the
direct, FFT, and FFTW-oracle tests, performs short direct and FFT convergence
comparisons plus a `16^3 -> 32^3` candidate-search smoke run, and verifies
bit-for-bit checkpoint/restart identity. It uploads the CSV products as
workflow artifacts.

## Interpretation guardrails

Every fixed Galerkin cutoff is a smooth finite-dimensional ODE, so it cannot by
itself demonstrate PDE singularity. In particular:

- growth that changes when `K`, the grid, or `dt` changes is a resolution
  artifact until proved otherwise;
- cutoff-shell energy above roughly one percent is reported as an
  under-resolution warning, not evidence of blow-up;
- `cutoff_shell_ok=true` is only one necessary resolution check, not a general
  certificate that the run is converged;
- the sampled vorticity maximum is a lower estimate of the truncated field's
  true maximum, while the Fourier sum is an often-loose upper bound;
- the sampled BKM integral is finite-run telemetry, not the hypothesis or
  conclusion of a theorem;
- finite growth of H1/2 or any other critical norm is a triage signal; only an
  appropriate unbounded or non-integrable limiting behaviour could support a
  blow-up argument;
- a small raw change between adjacent rescaled profiles is meaningless unless
  the accompanying time or spectral-scale change is controlled; the
  scale-normalized drift is still only a heuristic;
- the fitted exponential-tail slope is a floating-point regression over a
  short retained range, not a certified analyticity radius or Fourier-tail
  bound;
- ordinary floating point cannot certify inequalities needed by a proof.

## Research gates

The exact low-cutoff oracle, named benchmark fields, shell accounting, and the
cutoff/timestep comparison harness are implemented. The strictly dealiased FFT
backend, conservative adaptive timestep control, smooth parameterized
vortex-tube family, both critical-norm diagnostics, enstrophy budget,
energy-parameter sweep, saved single-run time series, and resolution-gated
search are implemented. The sharp critical-growth candidate has a short-time
`32^3 -> 64^3` check and a cutoff-clean `N=64` continuation to `t=0.2`.
Checkpoint/restart, the independent FFTW oracle, rescaled-spectrum output,
scale-normalized profile drift, and heuristic tail fitting are implemented.
The present sharp candidate fails the stationary-profile gate, so a `128^3`
run of exactly the same geometry is deprioritized. The candidate search now
uses fixed-forward-scale profile windows, explicit profile-drift and cutoff
costs, a conservative paired coarse/fine score, and a strict refinement gate.
Only candidates whose critical-norm growth, scale motion, profile stationarity,
and cross-resolution agreement all pass that gate may seed a narrower sweep.
The first 12-case `32^3 -> 64^3` neighbourhood search produced no survivor.
The next engineering milestone is therefore a second full FFTW-based evolution
path and trajectory-level comparison. After that independent check passes, the
search can broaden the initial-data family instead of spending larger grids on
a tube geometry that fails its profile gate.

A credible path from this scaffold to a theorem has several hard gates:

1. **Numerical credibility:** extend the established FFTW coefficient oracle
   into a second full trajectory implementation; reproduce standard
   benchmarks; run convergence studies across cutoff, time step, box, and
   precision; and search for stable rescaled profiles rather than isolated
   spikes.
2. **Analytic mechanism:** state a scale-by-scale transfer lemma that controls
   viscosity, nonlocal frequency interactions, pressure/Leray projection, and
   the time accumulated over infinitely many stages. Track a critical norm
   such as L3, not only supercritical quantities.
3. **Computer-assisted proof:** reformulate any candidate profile or cascade
   as a fixed-point/stability problem; use interval arithmetic; certify
   eigenpairs and nonlinear residuals; bound unresolved Fourier tails; and
   connect the construction back to admissible smooth initial data.
4. **Independent audit:** publish theorem statements, machine-readable
   certificates, resource requirements, and minimal verification code, then
   obtain expert review of every functional-analytic implication.

The opposite direction is equally useful: computations may suggest a
quantitative geometric-depletion or critical-norm estimate strong enough to
rule out an infinite cascade. Either outcome must ultimately become a rigorous
estimate, not a plot.

## Primary references

- C. Fefferman, *Existence and Smoothness of the Navier-Stokes Equation*, Clay
  Mathematics Institute problem statement:
  <https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf>
- T. Tao, *Finite time blowup for an averaged three-dimensional Navier-Stokes
  equation*: <https://arxiv.org/abs/1402.0290>
- S. Palasek, *Arbitrary norm growth in the 3D Navier-Stokes equations*:
  <https://arxiv.org/abs/2509.18595>
- T. Hou, Q. Wang, and D. Yang, computer-assisted weak non-uniqueness from
  singular initial data: <https://arxiv.org/abs/2509.25116>
- Z. Grujic, work on geometric depletion of vortex stretching:
  <https://arxiv.org/abs/2607.08866>
- L. Escauriaza, G. Seregin, and V. Sverak, the endpoint L3 regularity
  criterion:
  <https://www.pdmi.ras.ru/~seregin/Recent%20Publications/engESS3.pdf>
