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
- a sampled scale-critical L3 velocity norm;
- sampled maximum vorticity and its time integral (a BKM-inspired diagnostic);
- a rigorous-for-the-truncated-polynomial Fourier upper bound on maximum
  vorticity;
- spectral centroid and energy fraction touching the cutoff shell;
- energy, nonlinear transfer, viscous loss, and forward flux for every radial
  Fourier shell;
- divergence, Fourier-reality, and semi-discrete energy-balance defects.

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
cancellation, shell accounting, viscous energy dissipation, the ABC negative
control, and short-run numerical stability.
The FFT-specific suite also checks a complex 3D transform round trip, rejects
unsafe cutoffs, compares every nonlinear Fourier coefficient against direct
convolution, repeats that comparison with all modes populated at the maximum
safe cutoff, and compares complete short trajectories. It now also verifies
that straight tubes have no non-zero axial modes, bent tubes do, both remain
real and divergence-free, invalid geometry is rejected, and adaptive steps
obey both requested stability bounds.
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

Run the faster FFT backend on a `32^3` grid, retaining the strict safe cutoff
`K=10`:

```bash
./build/research/navier_stokes_cascade/navier_stokes_fft \
  --grid 32 --cutoff 0 --initial-condition taylor-green \
  --dt 0.005 --final-time 0.1 --cfl 0.4 --diagnostic-every 20 \
  --output fft-n32.csv --shell-output fft-n32-shells.csv
```

Here `--dt` is the maximum allowed adaptive step. Add `--fixed-dt` to recover
fixed-step evolution; if `--final-time` is omitted, the end time is
`--steps * --dt`. The executable writes a time-series CSV and a separate
long-form shell CSV.

Run one member of the smooth vortex-tube family with:

```bash
./build/research/navier_stokes_cascade/navier_stokes_fft \
  --grid 32 --initial-condition vortex-tubes \
  --tube-core 0.7 --tube-separation 1.8 \
  --tube-bend 0.3 --tube-axial-mode 2 \
  --viscosity 0.02 --dt 0.005 --final-time 0.1 \
  --output tube.csv --shell-output tube-shells.csv
```

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
  --viscosity 0.02 --dt 0.005 --final-time 0.1 --top 3 \
  --output navier_stokes_candidate_search.csv
```

The straight cases are deduplicated because their axial wavenumber has no
effect. The CSV records initial/final/peak critical L3, sampled vorticity,
enstrophy, palinstrophy, maximum positive shell flux, cutoff contamination,
constraint defects, accepted-step statistics, and coarse/fine differences.
An apparent growth signal must occur on both grids to pass the preliminary
cross-resolution gate. The ranking score uses critical-L3 growth first with a
small sampled-vorticity and final-L3 tie-break; it is merely a deterministic
triage rule.

The direct comparison uses one common spatial sampling grid for every cutoff.
The FFT comparison samples on each native grid; therefore its L3-grid
difference includes quadrature error as well as solution error. The CSV
contains total and shell energies, critical-L3 ratio, nonlinear transfer,
viscous dissipation, forward flux, constraint defects, and cutoff-shell checks.
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
| Peak sampled vorticity / initial | 1.03280068 | 1.04312237 |
| Final enstrophy / initial | 0.981387172 | 0.981387447 |
| Maximum positive forward flux | 0.050832943 | 0.050832879 |
| Peak cutoff-shell energy fraction | 4.06726e-7 | 4.56699e-14 |
| Accepted adaptive steps | 68 | 141 |
| Maximum conservative CFL bound | 0.35 | 0.35 |

The final-L3 ratios agree to about `8.2e-7` relative and the peak forward flux
to about `1.3e-6` relative. The sampled-vorticity growth differs by about 0.99%
relative, partly because the maximum is sampled on each native grid.

Extending the same candidate to `t=0.5` at `N=32` gave peak sampled-vorticity
growth `1.12599`, enstrophy growth `1.09518`, final critical-L3 ratio
`0.932052`, and cutoff fraction `4.96984e-4`. The corresponding `N=16` run was
under-resolved (cutoff fraction `0.0200`), so the long-time growth still needs
an `N=64` check.

The defensible interpretation is narrow: this family exhibits resolved local
vortex amplification and forward transfer over the short interval, while its
scale-critical L3 norm decreases. That is a useful configuration for studying
vortex stretching or geometric depletion, but it is presently evidence
against this particular run being a blow-up candidate—not a proof of global
regularity and not a disproof of Navier-Stokes.

Use `--help` for all parameters. The direct backend still grows quadratically
in the retained mode count; use it to audit small cases and the FFT backend to
explore larger ones.

The branch-scoped GitHub Actions workflow builds these CMake targets, runs the
direct and FFT invariant tests, performs short direct and FFT convergence
comparisons plus a `16^3 -> 32^3` candidate-search smoke run, and uploads all
three CSV products as workflow artifacts.

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
- ordinary floating point cannot certify inequalities needed by a proof.

## Research gates

The exact low-cutoff oracle, named benchmark fields, shell accounting, and the
cutoff/timestep comparison harness are implemented. The strictly dealiased FFT
backend, conservative adaptive timestep control, smooth parameterized
vortex-tube family, and resolution-gated search are implemented. The next
engineering milestone is an independent FFT-library oracle, saved time-series
search traces, and `64^3` long-time validation of the stretching candidate.

A credible path from this scaffold to a theorem has several hard gates:

1. **Numerical credibility:** cross-check this small radix-2 implementation
   against an established FFT library; reproduce standard benchmarks; run
   convergence studies across cutoff, time step, box, and precision; and
   search for stable rescaled profiles rather than isolated spikes.
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
