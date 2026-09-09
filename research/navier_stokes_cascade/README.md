# Navier-Stokes cascade research scaffold

This directory is an **experimental research scaffold, not a claimed proof or
disproof** of the Navier-Stokes Millennium Prize problem. It gives the
repository a reproducible place to test one concrete direction: whether
finite-stage transfer of energy toward high frequencies can be organized into
an accumulating cascade, or whether viscosity and geometric depletion prevent
that transfer.

## Exact scope

The existing numerical search studies smooth, divergence-free initial data
for the unforced three-dimensional incompressible equations

```text
partial_t u + (u . grad)u = -grad p + nu Delta u,
div u = 0,
u(0) = u_0,
```

on the periodic three-torus, with positive viscosity. The full
[Clay formulation](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf)
also covers R^3 and allows suitably smooth forcing in its breakdown
alternatives C and D. A resolution must prove one of its precise alternatives;
the existing unforced numerical search is only one research route. Weak
non-uniqueness from singular data, blow-up for a modified/averaged equation,
or a large finite numerical value does not settle that statement.

The search now has two distinct routes. The historical `profile` route looks
for a relatively stationary **rescaled shell-energy spectrum**. The newer
`amplification` route looks for sustained finite critical-norm growth and
vortex stretching that survive held-out numerical checks. Profile rejection
does not rule out other singularity mechanisms. Conversely, a stationary
shell spectrum is weaker than a self-similar velocity field and proves
neither singularity nor regularity. Historical statements below about a
candidate being rejected or not worth refining refer to that earlier route
and its chosen compute budget.

A separate [paper building-block audit](results/paper_profile_audit/README.md)
checks explicit similarity identities, a scalar comparison function and a
viscous exterior from the user-supplied OpenAI manuscript. It retains the full
momentum residual and deliberately tests omitted terms. Its 29 regression
tests and 25 audit gates pass, but the nonlinear matched profile, correction
construction and complete proof remain unverified. No new numerical
candidate is promoted, and the earlier `.14` resolution failure is unchanged.

The subsequent [nonlinear axis pilot](results/nonlinear_axis_pilot/README.md)
constructs finite Taylor approximations to the manuscript's nonlinear inner
equations using an explicitly labelled analytic pressure seed. Thirteen
numerical gates pass, including degree/precision refinements and a full
momentum-residual crosscheck. The seed is not the completed outer schedule's
pressure; the annular connection and five-moment matching remain unbuilt.
Direct attachment to the heat exterior fails, and the full base-flow residual
is not small. This is reusable local construction machinery, not a matched
paper profile, a proof, or a new candidate.

The [outer pressure pilot](results/outer_pressure_pilot/README.md) now builds
the complete unedited swirl schedule from Appendix A.2, evaluates its axis
pressure integral, solves the two pressure-preserving angular bump equations,
and couples that datum to the nonlinear inner equations. Tiny tail terms and
an insufficient-precision failure are retained explicitly. The axial pulse's
M/J/S closure, global cone thresholds, heat compensation and five-moment axis
annulus remain open; this checkpoint does not promote a global solution.

The [axial stress audit](results/axial_stress_audit/README.md) reproduces a
stress-cone failure for the pilot's `Md=4` choice and obtains positive bounds
for both axial stress ratios at `Md=64`, over the complete axial interval
and all angles. An analytic reduction and outward interval arithmetic retain
the finite-parameter remainders and the small-angle memory term. The old
`lambda=1e-5` also fails a later necessary cone condition by more than `1e27`;
retuning that parameter is the next prerequisite. No complete outer cone,
matched profile, or blowup result is claimed.

## What is implemented

Three numerical paths now evolve a mean-zero, real, divergence-free velocity
field on the 2-pi periodic torus using the Fourier cube
`|k_x|, |k_y|, |k_z| <= K`:

- `GalerkinSystem` evaluates every Fourier triad directly. Its quadratic cost
  makes it slow, but its simple formula is the low-cutoff correctness oracle.
- `PseudospectralSystem` uses an in-repository radix-2 three-dimensional FFT,
  the rotational nonlinearity `P[u x curl(u)]`, and strict 2/3 de-aliasing. Its
  cost scales approximately as `N^3 log(N)` and permits materially larger runs.
- `FftwReferenceSystem` is an optional independent evolution oracle. It rebuilds
  the Fourier grid, uses FFTW3 transforms, writes out its own Leray projection,
  computes its own adaptive bound and diagnostics, and assembles all four RK4
  stages independently.

For each retained non-zero mode, all three paths compute

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

The optional FFTW3 target now compares forward and inverse transforms, the
complete dealiased Navier-Stokes right-hand side, independently computed
physical and spectral diagnostics, adaptive timestep bounds, and complete RK4
trajectories. The comparison executable gives both evolvers the exact same
initial Fourier coefficients and advances them on a common timestep equal to
the safer of their two independently proposed steps. It checks the full state
after every accepted step and returns failure when state, diagnostic, reality,
or divergence tolerances are exceeded. This isolates evolution-code errors,
but it is not fully independent mathematics: both paths solve the same
finite-dimensional Fourier model in double precision and intentionally share
the initial coefficients and elementary complex-vector types.

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
cutoff, Leray projected, and energy normalized.

The search and independent-comparison tools also expose an
`orthogonal-bundle` family. It superposes three copies of the same smooth
counter-rotating pair, with vector potentials aligned with the `x`, `y`, and
`z` axes:

```text
A = (w F_x, w F_y, F_z),
u = curl(A).
```

Here each `F_axis` is a periodic Gaussian difference whose two centres bend
along that axis; the `x` and `y` helices use opposite phase offsets. Thus the
field contains six mutually oriented tubes, remains divergence-free by
construction, and reduces exactly to the original `z` pair at `w=0`. The
orthogonal weight and phase are genuine geometry parameters, not amplitude
rescalings; total energy is normalized only after the three pairs are
combined. Both families are reproducible smooth probes, not claims that they
resemble a singular profile.

The third family is a localized interacting wave-packet triad. Three real
periodic Gaussian-envelope vector potentials carry the integer wavevectors

```text
k1 = (m,m,0),  k2 = (-m,0,m),  k3 = (0,-m,-m),  k1+k2+k3 = 0,
A = sum_j w_j G(x) cos(k_j . x + phi_j) a_j,
u = curl(A).
```

The exact triad relation permits immediate quadratic interaction, while the
envelope localizes the packets in physical space and broadens each carrier in
frequency. The curl, real sampling, safe truncation, Leray projection, and
energy normalization preserve periodicity, Fourier reality, and numerical
incompressibility. Width, carrier, secondary-packet weight, and phase are
independent search coordinates. This is motivated by the central role of
Fourier triads in Navier-Stokes energy transfer; it does not assume that one
triad can sustain an infinite cascade.

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
that straight tubes have no non-zero axial modes, bent tubes do, and the
six-tube bundle is energy-normalized, real, divergence-free, cutoff-clean,
distinct from the pair, and populated in all three velocity components. A
zero orthogonal weight must recover the original pair, invalid geometry is
rejected, and adaptive steps obey both requested stability bounds. The profile
tests check normalization,
amplitude invariance, an exact discrete scale shift with zero shape drift, and
metric ranges. Checkpoint
tests cover exact round trips, split/uninterrupted trajectory identity, and
checksum-corruption rejection. When FFTW3 is present, a separate target checks
the internal transforms, grid ordering, full nonlinear right-hand side,
including a state that populates every retained maximum-cutoff mode, adaptive
bound, independently sampled diagnostics, and a nonlinear 40-step
orthogonal-bundle trajectory against the FFTW evolution path.
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

To make the independent evolution oracle mandatory, install FFTW3 and add
`-DNS_CASCADE_REQUIRE_FFTW_REFERENCE=ON` to the configure command. The branch
workflow does this on Ubuntu so a missing reference library fails CI rather
than silently skipping the check.

Compare the main solver with the independent FFTW trajectory on the strongest
profile-rejected candidate using:

```bash
./build/research/navier_stokes_cascade/navier_stokes_fftw_compare \
  --grid 32 --initial-condition vortex-tubes \
  --tube-core 0.50 --tube-separation 1.30 \
  --tube-bend 0.35 --tube-axial-mode 2 \
  --energy 10 --viscosity 0.02 --dt 0.005 --final-time 0.08 \
  --cfl 0.35 --diagnostic-every 50 \
  --output navier_stokes_fftw_comparison.csv
```

Select the six-tube geometry with `--vortex-family orthogonal-bundle`; use
`--orthogonal-weight` and `--phase-offset` to set its relative pair strength
and helical phase. The ordinary pair remains the default.

The CSV contains both copies of every key diagnostic, the maximum coefficient
difference, relative full-state error, and constraint defects. The state error
is checked after every step even when a row is not written. By default the
command fails if the state error or scaled diagnostic difference exceeds
`1e-9`; both gates are configurable. `--fixed-dt` is available for controlled
timestep studies.

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
  --families pair,orthogonal-bundle,wave-packets \
  --cores 0.55,0.70 --separations 1.2,1.8 \
  --bends 0,0.30 --axial-modes 1,2 \
  --orthogonal-weights 0.5,1.0 --phase-offsets 0,1.0471975512 \
  --packet-widths 0.8,1.1 --carrier-modes 1,2 \
  --packet-weights 0.75,1.0 --packet-phases 0,1.0471975512 \
  --energies 1,4,10 --viscosity 0.02 \
  --dt 0.005 --final-time 0.1 --top 3 \
  --profile-bins 64 --profile-scale-window 0.025 \
  --profile-drift-threshold 1 \
  --output navier_stokes_candidate_search.csv
```

The straight cases are deduplicated because their axial wavenumber and phase
have no effect. Pair candidates do not multiply over bundle-only parameters.
Wave-packet candidates likewise do not multiply over tube-only parameters.
`--energy E` remains the single-energy shorthand; `--energies` includes energy
in the candidate grid. The CSV records the family and bundle geometry together
with initial/final/peak critical L3
and H1/2, sampled vorticity, enstrophy, palinstrophy, enstrophy
production/destruction, maximum positive shell flux, cutoff contamination,
constraint defects, accepted-step statistics, fixed-scale profile drift,
score components, and coarse/fine differences. Profile drift is sampled
whenever `log(k_rms)` advances by `--profile-scale-window`, independently of
`--diagnostic-every`; changing CSV verbosity therefore cannot change the scale
windows. Only forward movement to finer scales completes a window. The cutoff
fraction used by the resolution gate and score is also checked after every
accepted step, so a short contamination spike cannot hide between CSV rows.

Optimize the three continuous wave-packet coordinates at a fixed carrier with
checked central differences and a cutoff-gated backtracking line search:

```bash
./build/research/navier_stokes_cascade/navier_stokes_optimize \
  --grid 16 --carrier-mode 1 \
  --initial-width 1.1 --initial-weight 0.75 \
  --initial-phase 1.0471975511965976 \
  --energy 10 --viscosity 0.02 --final-time 0.08 \
  --iterations 2 --gradient-step 0.015 \
  --cutoff-threshold 0.01 \
  --output navier_stokes_packet_optimization.csv
```

The optimizer uses normalized bounded coordinates. Each derivative is computed
at step `h` and `h/2`, then Richardson extrapolated; excessive disagreement
rejects the entire gradient. Its smooth finite-time objective rewards terminal
H1/2 and L3, forward characteristic-scale motion, and penalizes smooth
scale-normalized spectral-profile change plus terminal cutoff loading. Peak
cutoff loading remains a hard acceptance gate. This is a parameter-space
finite-difference optimizer, not yet a discrete adjoint or a proof tool.

Check the tangent-linear equation and every differentiated RK4 stage against a
sequence of centered directional differences:

```bash
./build/research/navier_stokes_cascade/navier_stokes_tangent_check \
  --grid 16 --dt 0.0005 --final-time 0.01 \
  --epsilons 0.001,0.0005,0.00025 \
  --output navier_stokes_tangent_check.csv
```

This check uses fixed timesteps, so it differentiates the RK4 flow map rather
than the nonsmooth minimum inside adaptive timestep selection. The perturbation
is divergence-free, real, and projected tangent to the fixed-energy sphere.
For `F(u)=P[u x curl(u)]-nu A u`, the implementation evolves
`DF(u)v=P[v x curl(u)+u x curl(v)]-nu A v` directly through the FFT backend.

Check the analytical adjoint of that tangent equation and the exact reverse
pass through all four fixed RK4 stages against both the forward tangent and
centered finite differences of a terminal critical-norm objective:

```bash
./build/research/navier_stokes_cascade/navier_stokes_adjoint_check \
  --grid 16 --dt 0.0005 --final-time 0.01 \
  --epsilons 0.004,0.002,0.001 \
  --output navier_stokes_adjoint_check.csv
```

On the real, solenoidal retained state space, the implemented RHS adjoint is

```text
DF(u)^* lambda = P[curl(u) x lambda + curl(lambda x u)]
                 - nu A lambda.
```

The reverse RK4 routine differentiates the discrete fixed-step map itself; it
does not differentiate adaptive step selection. The terminal objective in the
standalone check is `J(u)=0.5 ||u||_(H1/2)^2`, whose Fourier gradient is
`lambda_k=|k|u_k`.

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

### Independent trajectory result

The selected profile-rejected case was then replayed to `t=0.08` on an
`N=32, K=10` grid through the internal FFT and independent FFTW evolution
paths. Each of the 381 accepted steps used the smaller of the two separately
computed CFL bounds:

| Check | Internal FFT | Independent FFTW |
|---|---:|---:|
| Final H1/2 / initial | 1.044539612676717 | 1.044539612676717 |
| Final L3 / initial | 0.9635897800235352 | 0.9635897800235353 |
| Final enstrophy / initial | 1.412584781446576 | 1.412584781446577 |
| Final characteristic wavenumber / initial | 1.210592208763490 | 1.210592208763490 |

Across all accepted steps, the peak relative Fourier-state difference was
`1.10305e-15`. Across sampled diagnostics, the peak scaled difference was
`6.64882e-16`; the largest divergence and reality defects were `1.24127e-15`
and `3.95053e-16`. Thus the finite 4.45% H1/2 growth is not explained by a bug
specific to the in-repository FFT or its RK4 assembly. This does not rescue the
candidate: its large and rebounding rescaled-profile drift is a property of the
agreed trajectory. Agreement between two floating-point solvers also cannot
rule out common equation, truncation, modelling, or finite-resolution errors.
A separate 15-step `N=64, K=21` smoke evolution to `t=0.002` also passed, with
peak relative state difference `5.81702e-17`.

### Orthogonal-bundle pilot

The first broader-family pilot used the six-tube construction motivated by
reconnection-rich extreme-flow computations. It varied
`core=0.55,0.70`, `separation=1.2,1.6`, `bend=0.20,0.35`, and orthogonal
weight `0.5,1.0`, with axial mode two, phase offset `pi/3`, energy ten,
`nu=0.02`, and `t=0.08`. This was a bounded `16^3 -> 32^3` triage sweep, not
an adjoint optimization or an exhaustive search.

Only one of the 16 candidates passed the one-percent cutoff gate already at
`N=16`: `core=0.70`, `separation=1.60`, `bend=0.20`, and weight `0.5`.
Its `N=32` rerun passed the preliminary cross-resolution gate but not the
profile-refinement gate:

| Check through t=0.08 | Cross-resolved candidate | Most active N=32 finalist |
|---|---:|---:|
| Core / separation / bend / weight | 0.70 / 1.60 / 0.20 / 0.5 | 0.55 / 1.20 / 0.35 / 1.0 |
| Peak H1/2 / initial | 1.021059685 | 1.045863570 |
| Peak sampled vorticity / initial | 1.256117215 | 1.497489057 |
| Peak enstrophy / initial | 1.161895002 | 1.361950240 |
| Final characteristic wavenumber / initial | 1.086385801 | 1.185241615 |
| Latest rescaled-profile drift | 8.581635012 | 7.366613616 |
| Peak cutoff-shell energy fraction at N=32 | 8.90323e-5 | 2.89824e-3 |

The active narrow-core case was under-resolved at `N=16` (cutoff fraction
about `5.25e-2`), so it did not pass the cross-resolution gate even though its
`N=32` trajectory was cutoff-clean. Replaying that full `N=32, K=10`
trajectory through `t=0.08` with the internal FFT and independent FFTW
evolvers took 482 shared adaptive steps. The peak relative state difference
was `8.21632e-16`, the peak scaled diagnostic difference was `1.38476e-15`,
and both paths produced the ratios in the right-hand column.

The broader geometry therefore generated stronger finite vorticity growth
than the cross-resolved case, but neither candidate approached a stationary
rescaled spectrum: both drifts remained more than seven times the threshold,
and the active case's drift rebounded from a minimum near `3.07`. Zero of 16
candidates passed the strict refinement gate, so no `N=64` promotion was
justified. In short: more tubes made more drama, but not the missing
self-similar mechanism.

### Interacting wave-packet pilot

A separate 16-case pilot varied envelope width `0.8,1.1`, carrier mode `1,2`,
secondary-packet weight `0.75,1.0`, and phase `0,pi/3`, at energy ten,
`nu=0.02`, and `t=0.08`. It used the same `16^3 -> 32^3` cutoff,
cross-resolution, and profile gates as the vortex searches.

Three coarse candidates passed the one-percent cutoff gate. The best resolved
`N=32` case used width `1.1`, carrier one, weight `0.75`, and phase `pi/3`:

| Check through t=0.08 | N=16 | N=32 |
|---|---:|---:|
| Peak H1/2 / initial | 1.039069558 | 1.039894160 |
| Peak sampled vorticity / initial | 1.143555544 | 1.360050385 |
| Peak enstrophy / initial | 1.282393587 | 1.296462862 |
| Final characteristic wavenumber / initial | 1.140217559 | 1.146475609 |
| Latest rescaled-profile drift | 7.537323850 | 7.658901048 |
| Peak cutoff-shell energy fraction | 7.64993e-3 | 5.18432e-5 |

The scale and H1/2 changes agree fairly closely, but sampled-vorticity growth
and cutoff history do not satisfy the conservative cross-resolution gate.
More importantly, profile drift remains over seven times the threshold. An
independent FFTW replay of the full `N=32, K=10` trajectory used 219 shared
adaptive steps and matched the internal evolution to a peak relative state
difference of `7.60244e-16`; peak scaled diagnostic disagreement was
`6.76855e-16`. Thus the finite amplification is reproduced, but zero of 16
candidates qualifies for refinement and no `N=64` promotion is justified.

### Checked-gradient packet optimization

Starting from the resolved width `1.1`, weight `0.75`, phase `pi/3`, carrier-one
packet, a two-step `N=16, t=0.08` run used normalized perturbation `h=0.015`.
The first gradient's coarse/refined disagreement was `0.0891`; the second's was
`0.00584`. Both backtracking steps improved the evaluated objective without
crossing the one-percent peak-cutoff gate:

| Measurement | Initial | Optimized N=16 | Optimized N=32 |
|---|---:|---:|---:|
| Width / weight / phase | 1.100 / 0.750 / 1.047 | 1.088 / 0.758 / 0.299 | same |
| Smooth objective | -0.02027 | -0.01335 | 0.01475 |
| Peak H1/2 / initial | 1.03907 | 1.04744 | 1.04820 |
| Characteristic wavenumber / initial | 1.14022 | 1.16200 | 1.16849 |
| Smooth endpoint profile drift | 1.93271 | 1.91310 | 1.90936 |
| Strict windowed profile drift | 7.53732 | 7.46301 | 7.62088 |
| Peak cutoff-shell energy fraction | 7.64993e-3 | 9.59527e-3 | 9.26224e-5 |

The objective improvement transfers to the fine grid, but the strict profile
drift remains far above one and sampled-vorticity growth still disagrees across
resolutions. The optimized `N=32` trajectory therefore fails the existing
promotion gates. Its 221-step independent FFTW replay passed with peak relative
state disagreement `8.40746e-16`. This is evidence that optimization found a
slightly stronger finite cascade, not evidence of singular behaviour.

### Tangent-linear trajectory verification

The first fixed-step tangent check used the optimized packet on `N=16, K=5`
through `t=0.01`. A separate smooth solenoidal packet was orthogonally
projected onto the fixed-energy tangent space. Centered finite differences
converged to the directly evolved tangent as follows:

| Epsilon | Relative tangent error | Observed order |
|---:|---:|---:|
| 1.0e-3 | 2.99731e-10 | - |
| 5.0e-4 | 7.49884e-11 | 1.99893 |
| 2.5e-4 | 1.88739e-11 | 1.99027 |

The primal state produced inside differentiated RK4 was bit-identical to the
ordinary RK4 state. Final tangent divergence and Fourier-reality defects were
`4.44306e-17` and `1.24321e-16`. The nearly quadratic error decrease is the
expected centered-difference verification of the tangent implementation.

### Reverse discrete-adjoint verification

The corresponding fixed-step reverse check used the same optimized `N=16,
K=5` base trajectory through `t=0.01`. The RHS and complete-trajectory duality
errors, normalized by the associated vector-norm products, were
`1.90094e-17` and `3.74769e-17`. Forward tangent and reverse adjoint evaluation
of the terminal directional derivative both gave `-0.665160097044`.

Centered differences of `0.5 ||u(t)||_(H1/2)^2` converged to that derivative:

| Epsilon | Relative gradient error | Observed order |
|---:|---:|---:|
| 4.0e-3 | 2.88144e-8 | - |
| 2.0e-3 | 7.20551e-9 | 1.99962 |
| 1.0e-3 | 1.78960e-9 | 2.00947 |

The differentiated and ordinary primal states were bit-identical. The initial
adjoint's divergence and Fourier-reality defects were `8.89911e-16` and
`8.01300e-16`. These identities verify the reverse derivative of this
finite-dimensional fixed-step solver; they do not validate a PDE singularity
claim or remove the need for resolution and independent-evolution gates.

### Constrained many-mode Fourier optimization

The verified reverse RK4 pass now drives an initial-state optimizer rather
than only checking one derivative. Its variables are every non-zero Fourier
coefficient in the cube `|k_i| <= K_seed`, restricted exactly to real,
divergence-free fields and to a fixed kinetic-energy sphere. The resulting
real dimension is

```text
2 ((2 K_seed + 1)^3 - 1),
```

so `K_seed=2` gives 248 variables and `K_seed=3` gives 684. This follows the
large-scale adjoint-optimization strategy used by Kang, Yun, and Protas, while
retaining this project's critical-norm and resolution gates. The smooth
differentiated objective is

```text
J = log(H1/2(T) / H1/2(0))
    + 0.15 log(k_rms(T) / k_rms(0))
    - 0.04 log(1 + (E_cutoff(T) / E(T)) / 0.01)
    - 0.05 P_endpoint
    - 0.05 P_path.
```

`P_endpoint` compares nine smooth Gaussian measurements of the initial and
final spectra in the coordinate `log(|k|/k_rms)`. Each measurement is divided
by total energy, and the penalty is the mean squared change of their
regularized logarithms. It is therefore insensitive to amplitude and to a
pure shift of spectral scale. For snapshot penalties
`p_j = P_shape(u(0), u(j T / S))`, the default path term is now

```text
P_path = tau log((1 / S) sum_j exp(p_j / tau)),
S = 8, tau = 0.01.
```

This normalized log-sum-exp lies between the arithmetic mean and the hard
maximum. Its derivative gives the worst-changing snapshots the largest
softmax weights while remaining smooth enough for the discrete adjoint. The
integrator lands exactly on the fixed physical times, independently of
diagnostic cadence. The analytic gradient includes the state dependence of
`k_rms` and injects every weighted snapshot derivative at the corresponding
point in the reverse RK4 sweep. Centered differences test the shape derivative
both directly and through the full path. `--profile-shape-weight`,
`--profile-path-weight`, `--profile-path-samples`, and
`--profile-path-temperature` expose the weights and schedule. Historical
arithmetic-mean runs remain reproducible with
`--profile-path-aggregation mean`. These differentiable quantities only guide
optimization. The independent cloud-in-cell L1 drift remains the promotion
gate.

The terminal gradient and intermediate snapshot gradients are reversed
through every stored fixed-step RK4 stage, then projected onto the low-band,
real-solenoidal tangent space. Trial states follow energy-sphere geodesics.
Every iteration checks the adjoint slope by a centered geodesic difference,
requires Armijo improvement in `J`, requires an improvement in the existing
profile-aware search score, and preserves the hard cutoff and invariant
gates. Each accepted coarse state is also run on the fine grid. The saved
output is the best paired state encountered, not merely the last coarse
iterate; this prevents coarse-grid overfitting from replacing a better
candidate. Full-precision coefficient CSVs use one shared strict reader and
can be reloaded exactly by both the optimizer and independent FFTW trajectory
tool.

`--starts N` adds reproducible, hash-generated tangent directions around the
same base state and places them at `--start-angle` on the fixed-energy sphere.
Start zero is always the unperturbed control. The trace records `start_index`,
and `--start-offset I` can replay one numbered basin without rerunning earlier
starts. A loaded `--state-input` can now serve as the center of the same
deterministic local-start construction; exact replay uses one start at offset
zero. Every rejected line-search state is recorded as `line-trial`, making a
surrogate-versus-hard-score conflict auditable.

For example:

```bash
./build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --grid 16 --fine-grid 32 --seed-bandwidth 3 \
  --initial-family wave-packets --energy 10 --viscosity 0.02 \
  --dt 0.0005 --final-time 0.08 --iterations 3 \
  --starts 4 --start-angle 0.35 \
  --profile-shape-weight 0.25 --profile-path-weight 0.25 \
  --profile-path-samples 8 --profile-path-temperature 0.01 \
  --output navier_stokes_state_optimization.csv \
  --state-output navier_stokes_optimized_state.csv
```

The bounded three-family pilot selected the widened wave-packet branch:

Its exact full-precision coefficients are preserved in
`candidates/wave_k3_t008_screening.csv`; the candidate-directory warning is
part of that artifact's interpretation.

| Measurement | Coarse `N=16, K=5` | Fine `N=32, K=10` |
|---|---:|---:|
| Peak H1/2 / initial | 1.090335 | 1.091131 |
| Peak L3 / initial | 1.001521 | 1.001642 |
| Final L3 / initial | 0.989143 | 0.990212 |
| Peak sampled vorticity / initial | 1.104702 | 1.117222 |
| Final `k_rms` / initial | 1.252101 | 1.257258 |
| Latest profile drift (fixed-threshold audit) | 6.301960 | 6.284336 |
| Peak cutoff-shell fraction | 0.0074097 | 0.00004135 |

All three accepted slopes had relative adjoint/difference errors below
`1.5e-5`. Halving the coarse timestep to `0.00025` changed the reported H1/2
ratio only in the eleventh decimal place and `k_rms` ratio in the tenth. A
shape-aware step reduced the smooth endpoint penalty, but the fine trajectory
reached one additional forward-scale window and its latest L1 drift rebounded
to `9.18`; paired selection correctly retained the earlier state. The
`K_seed=3` orthogonal-bundle seed exceeded the one-percent coarse cutoff gate
and was rejected.

This is a stronger, reproducible finite-cascade candidate than the earlier
hand-parameterized packet, but it still fails decisively: the profile drift is
over six times the stationarity threshold, L3 falls by the final time, and the
coarse cutoff margin is modest. It is not eligible for a larger-grid or FFTW
promotion. It therefore motivated multi-start optimization of the same
constrained space with an explicitly differentiable profile-shape cost, not a
longer run of this rejected state.

### Smooth-profile deterministic multi-start result

That next search is now implemented. A bounded four-start `K_seed=2` pilot at
`T=0.04`, using shape weight `0.25`, accepted two checked steps from every
start. The largest adjoint/difference slope error was `1.42e-5`. Fixed coarse
and adaptive fine trajectories selected deterministic start 3. A continuation
with shape weight `0.25`, followed by a profile-dominant weight `1.0`, each
accepted one further step and then failed to find a jointly improving line
step. The best profile-oriented state is preserved exactly in
`candidates/wave_k2_profile_t004_screening.csv`.

| Measurement | Coarse `N=16, K=5` | Fine `N=32, K=10` |
|---|---:|---:|
| Peak/final H1/2 / initial | 1.013871 | 1.013876 |
| Peak L3 / initial | 1.000000 | 1.000000 |
| Final L3 / initial | 0.997357 | 0.997320 |
| Peak sampled vorticity / initial | 1.140758 | 1.166421 |
| Final `k_rms` / initial | 1.040772 | 1.040805 |
| Smooth endpoint shape penalty | 0.012229 | 0.012370 |
| Latest windowed L1 drift | 5.795320 | 5.838034 |
| Peak cutoff-shell fraction | 0.00014424 | 0.0000000202 |

The smooth penalty fell by about 71% from the unperturbed fine-grid baseline,
and the strict drift fell from about `11.12` to `5.84`. This is a real
improvement and agrees across resolutions, but it remains almost six times the
stationarity threshold. L3 also decreases. The state is therefore a screening
checkpoint, not a promoted blow-up candidate.

This pilot also exposed cumulative phase error in the profile-window
scheduler: each new window had previously been measured from the slightly
overshot prior crossing. Coarse and fine trajectories could consequently
report different window counts near an endpoint. Windows are now indexed by
integer multiples of `log(k_rms/k_rms(0))`; both the parameter search and the
state optimizer use the same tested rule. This changes no threshold and makes
late profile rebounds harder, not easier, to evade.

### Path-dependent `K_seed=3` screening result

The path-dependent extension and its staged search are complete. Four
deterministic 684-variable starts at `T=0.04`, with endpoint and path weights
both `0.25`, each accepted one checked ascent step. Start 2 won the paired
score, with fine-grid H1/2 growth `1.022052`, scale growth `1.064072`, cutoff
fraction `1.18e-6`, and strict drift `7.201165`. The largest finite-difference
adjoint error in this stage was `1.79e-4`, well below the `0.01` rejection
threshold.

The selected basin was continued at `T=0.06`. Balanced weights found no joint
improvement; increasing only the path weight to `1.0` accepted one step with
gradient error `1.32e-5`. That step moved the `T=0.08` coarse cutoff fraction
from `0.0102413` to `0.00997863` and reduced fine-grid drift from `4.746706`
to `4.605072`. A final full-horizon gradient check, using a smaller
energy-sphere angle so both probes remained inside the cutoff gate, had
relative error `2.58e-7`; no admissible line step improved both objectives.

The retained state is
`candidates/wave_k3_path_t008_screening.csv` (SHA-256
`cd73fe95e4a8c69f193dff3ebe5abd07a8716bf3739ae412a30b95b008b59491`).

| Measurement | Coarse `N=16, K=5` | Fine `N=32, K=10` |
|---|---:|---:|
| Peak/final H1/2 / initial | 1.054411 | 1.055435 |
| Peak L3 / initial | 1.000581 | 1.000560 |
| Final L3 / initial | 0.985774 | 0.985337 |
| Peak sampled vorticity / initial | 1.279198 | 1.488865 |
| Final `k_rms` / initial | 1.161960 | 1.169504 |
| Smooth endpoint shape penalty | 0.085033 | 0.100035 |
| Historical four-snapshot mean path penalty | 0.035026 | 0.040059 |
| Latest windowed L1 drift | 4.629944 | 4.605072 |
| Completed profile windows | 15 | 15 |
| Peak cutoff-shell fraction | 0.00997863 | 0.000122778 |

An exact `N=32, T=0.08` replay through the independent FFTW/RK4 evolver used
276 common adaptive steps. Its peak whole-state and diagnostic disagreements
from the internal evolver were `8.27e-16` and `8.54e-16`. The trajectory is
therefore reproducible inside the truncated numerical model.

This is the lowest strict drift found so far: about 21% below the earlier
`K_seed=2` profile candidate. It is still not promotable. Drift is 4.6 times
the threshold, L3 falls, the coarse cutoff margin is only `2.14e-5`, and the
coarse/fine vorticity-amplification difference is about 14%, above its 10%
agreement gate. Larger grids remain deliberately gated.

### Smooth-maximum local screening result

The eight-snapshot smooth maximum was tested at `T=0.06` around the retained
path checkpoint. Four adjoint checks spanning path weights `1`, `2`, and `4`
and temperatures `0.01` and `0.005` had relative errors between `8.81e-7` and
`5.07e-6`. No smooth-maximum ascent direction produced a line step that
improved both the differentiated objective and the strict profile-aware
score. The trace explains why: at angle `0.00375`, the coarse smooth path
penalty fell from `0.024036` to `0.022590`, but strict drift rose from
`6.166903` to `6.345465`. The derivative is correct; the smooth proxy is
locally misaligned with the hard cloud-in-cell drift.

Six deterministic perturbations of angle `0.00375` were therefore screened
around the exact checkpoint. Local start 3 slightly improved the paired score;
a half-angle six-start refinement and a new adjoint step found no further
improvement. The retained control is
`candidates/wave_k3_smoothmax_t008_screening.csv` (SHA-256
`75818222d22aaaeb0b92f1a417a0309400e810335712b1df9aef729884fcc7b1`).

With diagnostics evaluated on every accepted step, its `T=0.08` measurements
are:

| Measurement | Coarse `N=16, K=5` | Fine `N=32, K=10` |
|---|---:|---:|
| Peak/final H1/2 / initial | 1.054363 | 1.055388 |
| Peak L3 / initial | 1.000572 | 1.000550 |
| Final L3 / initial | 0.985719 | 0.985288 |
| Peak sampled vorticity / initial | 1.280092 | 1.489476 |
| Final `k_rms` / initial | 1.161872 | 1.169418 |
| Smooth endpoint shape penalty | 0.085031 | 0.100033 |
| Eight-snapshot smooth-maximum path penalty | 0.065240 | 0.079745 |
| Latest windowed L1 drift | 4.606862 | 4.574379 |
| Peak cutoff-shell fraction | 0.00997425 | 0.000122795 |

Against the prior checkpoint under the same dense schedule, the worse of the
two drift values fell by `0.50%` (`4.629944` to `4.606862`), vorticity-ratio
disagreement fell from `14.0808%` to `14.0576%`, and the coarse cutoff margin
grew by about `20.5%`. Halving the coarse timestep preserved the direction:
the old/new worst drifts were `4.610665` and `4.588704`. However, the fine-grid
drift alone worsened from `4.562908` to `4.574379`, critical growth decreased
slightly, and vorticity disagreement remains above the `10%` gate. This is a
secondary screening checkpoint, not a promoted candidate.

An independent `N=32, T=0.08` FFTW/RK4 replay used 276 common adaptive steps.
Peak whole-state and diagnostic disagreements were `7.85e-16` and `8.94e-16`.
The small numerical change is reproducible, but it remains nowhere near proof
evidence.

One possible further objective for the profile route is to compare
adjacent rescaled spectra and divide their smooth shape change by the
corresponding `log(k_rms)` advance, then take a smooth maximum over those
local rates. That mirrors the strict windowed drift much more closely than
comparing every snapshot with the initial spectrum. Its gradient must include
both endpoints of every adjacent pair and the scale-advance denominator.

## Robust amplification search

The repository already contains the relevant PDE structure. Adding a generic
heat or Poisson solver does not discover singular initial data; those methods
are useful for verification and for deriving constraints on a search.

| Existing component | Role in candidate discovery |
|---|---|
| `include/maths/pde_variational_methods.hpp` | Galerkin/Ritz and energy-method examples; its 1D boundary-value discretizations are not a replacement for the 3D evolution |
| `include/maths/pde_numerical_methods.hpp` and Fourier methods | Linear diffusion, stability and spectral checks; the exact decaying-shear test exercises that PDE limit |
| `include/ns_cascade/pseudospectral.hpp` | The actual dealiased 3D equation, spectral norms, vortex-stretching budget, and viscosity |
| `include/ns_cascade/state_optimizer.hpp` | Checked discrete adjoint and constrained optimization of Fourier coefficients at fixed energy and bandwidth |
| `scripts/robust_search.py` | Diverse starts, conservative multi-objective shortlist, and held-out trajectory validation |

The first two paths in the table are relative to the repository root; the
remaining paths are relative to this research directory. Adjoint searches for
large finite-time enstrophy growth have precedent in
[Kang, Yun and Protas](https://arxiv.org/abs/1909.00041). Their numerical search
fixes initial enstrophy and varies the horizon; this pilot instead fixes
energy, initial bandwidth, viscosity and horizon within each comparison.
It is not a reproduction of their optimization or a rigorous growth bound.

With volume-normalized integrals, define

```text
E = (1/2) integral |u|^2,    Z = (1/2) integral |curl u|^2,
P = (1/2) integral |grad curl u|^2,    S = (grad u + grad u^T)/2.
dE/dt = -2 nu Z,
dZ/dt = integral omega . S omega - 2 nu P.
```

The solver already measures the two terms in the second balance. Production
divided by destruction exceeding one means instantaneous positive net
enstrophy production. The new search checks it at every fixed snapshot in
the second half of the trajectory, along with continued H1/2 growth. Neither
condition is necessary or sufficient for blow-up; together they help separate
sustained finite amplification from a peak that is already decaying.

`navier_stokes_state_optimize --search-track amplification` disables the two
profile penalties and uses the existing differentiated objective

```text
J = log(H1/2(T)/H1/2(0)) + 0.15 log(k_rms(T)/k_rms(0))
    - 0.04 log(1 + cutoff_fraction(T)/0.01).
```

Line steps must improve this objective by the Armijo condition, pass the
adjoint/difference check and retain all existing coarse numerical gates.
The retained state maximizes the smaller coarse/fine objective, prioritizing
pairs that are valid at both resolutions. The profile score and
`refinement_eligible` column remain visible as **profile-route** diagnostics;
they cannot veto amplification discovery. Default `--search-track profile`
preserves the historical behavior. Two trace columns are appended:
`search_track` and `selection_score`.

`--evidence-output` exports budgets and physical samples at `j*T/S`, including
the initial state and the exact midpoint. Both evolutions' physical L3
quadrature and sampled vorticity maximum use one `--evidence-grid` (default:
the fine grid). Zero padding evaluates the same Fourier polynomial more
densely; it does not recover unresolved dynamics or certify a continuum
maximum. Spectral budgets and cutoff fractions still belong to each actual
evolving truncation. These rows do not depend on `--diagnostic-every`.

Run the default seven-seed pilot from the repository root:

```sh
python3 research/navier_stokes_cascade/scripts/robust_search.py \
  --optimizer build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --oracle build/research/navier_stokes_cascade/navier_stokes_fftw_compare \
  --checkpoint research/navier_stokes_cascade/candidates/wave_k3_smoothmax_t008_screening.csv \
  --output-dir robust-amplification-pilot
```

The output directory must be new, so stale files cannot pass as a new result.
The manifest records every command, exit status, rejection, executable and
source checksum, trace, fixed-time budget and selected coefficient file.
Defaults search the packet, tube and orthogonal-bundle families with two
deterministic sphere starts each, plus the optional saved packet. Each seed
receives one adjoint step at `T=0.04`, `E=10`, `nu=0.02`, initial `K=3` and
evolution grids `16/32`. This is a bounded search of 684 real variables per
seed, not coverage of all initial data.
Family names identify the starting seeds; optimization varies all allowed
low Fourier modes and does not preserve an exact tube or packet ansatz.

The discovery shortlist retains non-dominated trade-offs in conservative
H1/2 growth, late H1/2 growth, scale advance, late stretching, cutoff loading
and cross-resolution disagreement. Round-robin objectives and family
diversity choose up to three held-out probes; membership alone is not a pass.
The exported thresholds are predeclared experimental choices:

- per-step cutoff fraction at most `0.008`, leaving margin below the solver's
  `0.01` hard cutoff gate;
- coarse/fine differences at most 2% in final H1/2 and common-grid L3 ratios,
  10% in final sampled-vorticity ratio and late minimum production/destruction,
  and 5% in enstrophy and characteristic-scale ratios;
- both resolutions have final H1/2 ratio at least `1.005`, second-half ratio
  at least `1.001`, no sampled late reversal, scale ratio at least `1.05`, and
  production/destruction above one throughout the sampled late window;
- each selected **unchanged initial field** runs at a longer `T=0.06`, again
  with step caps no larger than half their original caps or half their
  observed mean steps (both trajectories must use at least 1.9 times as many
  steps), and again with twice the spatial
  sampling grid and twice as many fixed snapshots;
- perturbations change endpoint ratios by at most `0.1%` (2% for vorticity
  under denser spatial sampling, and for the late stretching minimum), and
  the minimum critical growth exceeds three times the observed H1/2-ratio
  spread across all six held-out trajectories;
- an independent fine-grid FFTW trajectory passes before any state is marked
  `validated-finite-amplification-shortlist`. Omitting the oracle leaves it
  `screening-only`.

The measured spread is an empirical error indicator, not a certified bound.
Passing these gates earns a larger-resolution investigation, not a claim
about the infinite-dimensional PDE. L3 decay remains visible even if H1/2
passes. The analytic shifted-shear regression checks the physical budget,
known diffusion, common-grid samples and rejection of dropped modes. Python
tests exercise incomplete/corrupt evidence, smooth decay, cutoff margin,
route separation, multi-objective trade-offs and growth smaller than error.

The first completed seven-seed pilot retained three finite-amplification
finalists. Their worst held-out H1/2 gains were 5.70%, 5.13% and 1.52% at
`T=0.06`; all three passed the denser samples, actual timestep refinement and
independent FFTW checks. None passed the profile route and all final L3
ratios were below one. The leading frozen packet subsequently reached
7.90–7.91% H1/2 growth at `T=0.08` on `N=32/64`. Its H1/2 ratio agreed to
0.00403%, while the sampled-vorticity ratio still differed by 6.13%.
See the [full experimental record](results/README.md) for the exact data,
checksums, limitations and reproduction commands.

The subsequent frozen-field continuation reaches `T=.10` on `N=64/128`
with 9.91% H1/2 growth and about 2.02x sampled vorticity. The pointwise
vorticity discrepancy decreases from 16.23% on `32/64` to 4.63% on `64/128`,
passing the preliminary gate while leaving a material uncertainty.
The `N=64` independent trajectory and a new timestep-refinement check pass;
L3 still decreases. The evolved 128-grid checkpoint is archived in
[the resolution record](results/frozen_continuations_resolution/README.md).
The same record links a critical production–diffusion budget and checked
local adjoint source for a proposed late-growth search objective.

### Late-window critical-growth optimization

The local critical-rate gradient is now integrated into the **full reverse
discrete RK4 trajectory**. `--growth-objective late-rate` replaces the endpoint
H1/2 reward with `T * Phi`, where `gamma=d log(||u||_(H1/2))/dt` includes
viscosity and `Phi=-tau log(mean(exp(-gamma/tau)))` over fixed late times.
The weakest sampled rates get the largest adjoint weights. The endpoint
scale reward and cutoff penalty remain; the new mode requires the
`amplification` track. The default `endpoint` mode remains reproducible.

```sh
python3 research/navier_stokes_cascade/scripts/robust_search.py \
  --optimizer ./build/research/navier_stokes_cascade/navier_stokes_state_optimize \
  --oracle ./build/research/navier_stokes_cascade/navier_stokes_fftw_compare \
  --growth-objective late-rate --late-window-start .5 \
  --late-rate-samples 5 --late-rate-temperature .1 \
  --starts 1 --iterations 2 --max-finalists 2 \
  --discovery-time .04 --holdout-time .06 --output-dir new-late-search
```

The window starts at `.5T` and includes both ends; temperature has units of
inverse time. These choices are fixed before discovery. The optimizer lands
on each observation time even when it splits a timestep. The clock is
independent of diagnostic printing and distinct from profile snapshots.
`--late-rate-output` exposes every selected rate and its adjoint weight;
the search driver retains it automatically. The doubled-sampling holdout
also doubles the number of late-window intervals.

Importantly, `Phi >= min(gamma)`: a positive smooth objective can hide a
negative sampled rate. The robust driver separately requires the **actual
sampled minimum** positive and compares all common rates between evolution
grids. The normalized rate discrepancy `T*max|gamma_coarse-gamma_fine|`
must be at most `.002`; the minimum's change under timestep/sampling
perturbation must be at most `.001/T`. These are empirical screening gates,
not continuous-time bounds or singularity criteria. Early decay is allowed,
so delayed-growth fields are not discarded solely on their initial rate.

The [matched four-start pilot](results/late_growth_pilot/README.md) finds
small improvements in weakest sampled late growth, but both low-resolution
search arms fail their complete held-out validation. Full-trajectory
finite differences test aligned and split observation clocks, the terminal
source, fixed-energy projection, and diagnostic-cadence independence. CLI
checks also preserve historical endpoint evidence and reject an over-tight
gradient check without accepting an update.

The leading new field passes a separate frozen `32/64, T=.06` follow-up:
about 8.96% H1/2 growth, positive sampled late rates, actual timestep
refinement, doubled spatial/time samples, and a 491-step independent FFTW
trajectory agreeing to `1.41e-15` in the whole state. L3 still decreases
slightly and the relative critical-growth rate slows near the endpoint.
See the pilot record for failures, source provenance and the exact new
initial coefficient checkpoint; this is not a singularity claim.

The [frozen-field continuation](results/late_growth_continuation/README.md)
now reaches held-out `.08` and `.10`, with fine-grid H1/2 gains of 11.97%
and 14.55%. A smaller-step replay and the 945-step independent trajectory
pass, but L3 decreases and relative H1/2 growth keeps slowing. The `32/64`
vorticity and late-stretching gaps approach the 10% screening limit at `.10`.
The subsequent [64/128 comparison](results/late_growth_resolution128/README.md)
passes the retained resolution and endpoint sampling checks. At `.10`, the
H1/2 ratio gap is `0.00022%`, while sampled vorticity and late stretching
differ by `4.21%` and `0.73%`. Denser physical sampling gives fine-grid
vorticity `1.982` times its initial value. The remaining pointwise discrepancy,
decreasing L3 and slowing relative critical growth remain material limits;
this result supports another bounded continuation, not a singularity claim.

The [next frozen continuation](results/late_growth_t012/README.md) now reaches
`.12` on both grids and passes the same preliminary gates. Fine-grid H1/2 grows
`16.65%`, L3 falls `3.77%`, and densely sampled vorticity reaches `2.290` times
its initial value. The worst common-sampling vorticity gap is `7.19%` at `.11`;
the denser-sampling endpoint gap at `.12` is `2.13%`. Relative critical growth
slows further to `.792`. The new evolved states, budget checks and complete
resolution/sampling assessment are preserved; finite growth is not a blow-up proof.

The [subsequent .14 continuation](results/late_growth_t014/README.md) reaches
both declared targets. All checks pass at `.13`, but the `64/128` pair fails
local-vorticity agreement at `.14`: `12.39%` on common physical samples and
`12.65%` on doubled samples, exceeding the fixed `10%` limit. The fine-grid
H1/2 gain is `18.24%`, L3 falls `5.23%`, and relative critical growth slows to
`.557`. The global norm gap is only `0.0118%`, and cutoff loading remains low;
neither rescues the local-peak disagreement. The archived states and failed
assessment pause further time promotion of this grid pair. A finer evolution
comparison at `.14` is needed before a longer continuation can be promoted.

The [snapshot decomposition and capacity probe](results/late_growth_peak_audit/README.md)
now identify a small-energy, large-vorticity contribution from modes absent
on the coarse grid. At `.14`, those modes carry only `0.02687%` of fine energy,
yet filtering them out lowers the sampled fine peak by `15.59%`. Shared-mode
error remains substantial. The maximum sampled vorticity-vector difference
is `21.34%` of the fine peak, compared with the `12.65%` scalar-peak gap.
This post-hoc diagnostic does not relabel the prior gates or promote the pair.
Two short 128/256 RK4 steps agree to `1.13e-18` in relative whole-state L2,
but the measured 256-grid cost is `93–98 s` per step. A rough same-horizon
replay model is about six days of stepping; no such full replay was launched.
The next prerequisite is a faster validated high-resolution backend, followed
by a same-initial-field `.14` comparison and prospective vector-field checks.

### Checkpointed continuation of a frozen Fourier field

The optional FFTW target `navier_stokes_continue` loads the optimized CSV
coefficients unchanged and saves the **evolved state**, rather than another
copy of the initial field. It observes on a fixed physical clock and atomically
replaces its checksummed checkpoint after every observation. For example:

```sh
./build/research/navier_stokes_cascade/navier_stokes_continue \
  --state-input research/navier_stokes_cascade/candidates/wave_k3_amplification_t004_screening.csv \
  --grid 64 --sampling-grid 64 --dense-sampling-grid 128 \
  --viscosity 0.02 --dt 0.000125 --observation-interval 0.01 \
  --final-time 0.08 --output continuation-08.csv --checkpoint-output current.chk
./build/research/navier_stokes_cascade/navier_stokes_continue \
  --restart current.chk --final-time 0.10 \
  --output continuation-10.csv --checkpoint-output current.chk
```

The restart restores the backend, viscosity, cutoff, adaptive timestep controls,
sampling grids, observation clock and initial normalization measurements.
Scientific overrides are rejected. Final times are absolute multiples of the
saved observation interval. With the same executable, FFTW build and
floating-point platform, splitting at those times produces the same
checkpoint bytes as an uninterrupted run. Evidence files contain the shared
boundary row, so discard that duplicate when joining segments. Existing
evidence and input files are protected against output aliases and overwrites.
The `NSCONT1` format is distinct from the benchmark solver's `NSCCHK2` format;
neither executable silently interprets the other's restart.

At every accepted step the continuation checks finite, positive, nonincreasing
unforced energy and cutoff-shell energy. Exceeding the default cutoff fraction
`0.008` saves the stopped state and returns exit code 2, including between
scheduled observations. Such a checkpoint cannot resume. Other numerical or
I/O failures return 1 and leave the last successfully published checkpoint.
Fourier constraints and the full enstrophy budget are checked at observations.
The CFL rule controls stability heuristically; it is not a local error estimator.

Primary and denser physical samples are measured on the **same** evolving
Fourier field. Padding changes neither its spectral budget nor its dynamics.
This separates quadrature/maximum sampling changes from evolution-cutoff
changes. The CSV reports both sampling grids, absolute values and ratios to
their own initial measurements. BKM integration uses the primary samples on
the frozen observation clock and remains sampled finite-run telemetry.

The FFTW continuation uses compact RK4 storage while preserving the original
stage arithmetic. The four-derivative reference path remains available to the
independent trajectory checker. Tests cover exact shear diffusion, dense
maximum sampling, interacting radix-2/FFTW evolution, coefficient-by-coefficient
compact/full RK4 identity, checkpoint corruption, output protection, unchanged
restart controls and an inter-observation cutoff stop. Clean Ubuntu CI also
runs ASan and UBSan over the continuation and restart tests.

`scripts/continue_candidates.py --solver PATH --output-dir NEW_DIRECTORY`
replays the three frozen finalists on `32/64`, checks stages at
`T=0.08,0.10,0.12,0.16`, and pauses each pair when its empirical growth or
convergence gates fail. It retains commands, hashes, all stage evidence and
compressed final checkpoints. A failed resolution comparison is a reason to
raise the cutoff or pause inference, not a rejection of all singularity
mechanisms for that initial field.

Use `--help` for all parameters. The direct backend still grows quadratically
in the retained mode count; use it to audit small cases and the FFT backend to
explore larger ones.

The branch-scoped GitHub Actions workflow builds these CMake targets, runs the
direct, FFT, full FFTW-trajectory, checked-gradient, tangent-linear, and
reverse-adjoint and Fourier-state optimizer tests; performs short convergence
and candidate-search smokes; verifies exact coefficient save/reload and
bit-for-bit checkpoint/restart identity; and uploads the CSV products as
workflow artifacts.

## Interpretation guardrails

At a fixed cutoff the Galerkin energy bound controls every coordinate of its
finite-dimensional ODE and prevents finite-time divergence. It cannot by
itself demonstrate PDE singularity. In particular:

- growth that changes when `K`, the grid, or `dt` changes is a resolution
  artifact until proved otherwise;
- cutoff-shell energy above roughly one percent is reported as an
  under-resolution warning, not evidence of blow-up;
- `cutoff_shell_ok=true` passes one empirical resolution gate, not a general
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
- close agreement between the internal FFT and FFTW evolvers rules out many
  implementation-specific mistakes, but not errors shared by the common
  truncated equation, floating-point model, or starting coefficients;
- ordinary floating point cannot certify inequalities needed by a proof.

## Research gates

The exact low-cutoff oracle, named benchmark fields, shell accounting, and the
cutoff/timestep comparison harness are implemented. The strictly dealiased FFT
backend, conservative adaptive timestep control, smooth parameterized
vortex-tube family, both critical-norm diagnostics, enstrophy budget,
energy-parameter sweep, saved single-run time series, and resolution-gated
search are implemented. The sharp critical-growth candidate has a short-time
`32^3 -> 64^3` check and a cutoff-clean `N=64` continuation to `t=0.2`.
Checkpoint/restart, the full independent FFTW evolution oracle,
rescaled-spectrum output, scale-normalized profile drift, and heuristic tail
fitting are implemented. The FFTW path reproduces the selected pair,
orthogonal-bundle, and interacting-wave-packet `N=32` candidate trajectories
to round-off through `t=0.08`.
The present sharp candidate fails the stationary-profile gate, so a `128^3`
run of exactly the same geometry is deprioritized. The candidate search now
uses fixed-forward-scale profile windows, explicit profile-drift and cutoff
costs, a conservative paired coarse/fine score, and a strict refinement gate.
For the profile route, candidates whose critical-norm growth, scale motion,
profile stationarity and cross-resolution agreement pass that gate may seed a
narrower sweep. The separate amplification route above does not impose
stationarity of the rescaled spectrum.
The first 12-case `32^3 -> 64^3` pair neighbourhood search produced no
survivor. The initial-data family has now been broadened to an exactly
divergence-free six-tube orthogonal bundle, with structural tests, search
parameters, and independent-trajectory replay. Its first 16-case
`16^3 -> 32^3` pilot also produced no survivor: the stronger finite growth
still came with large, rebounding profile drift.

The localized wave-packet family is now implemented with structural tests,
separate search coordinates, workflow coverage, and independent-trajectory
replay. Its first 16-case pilot also produced no survivor. The most active
resolved packet amplified sampled vorticity by 36% at `N=32`, but its profile
drift was `7.66` and the conservative cross-resolution gate rejected it.

The checked finite-difference parameter optimizer, full tangent-linear RK4
evolution, and reverse-mode discrete adjoint are now implemented. The adjoint
satisfies `<DF(u)v,lambda>=<v,DF(u)^*lambda>` at both RHS and full-trajectory
levels and agrees with the established finite differences. The packet
optimizer's first two-step run improved both coarse and fine objectives, but
the result failed the unchanged promotion gates. Larger grids for the rejected
geometry remain deprioritized.

The constrained many-mode optimizer over real, solenoidal Fourier
coefficients is now implemented with exact fixed-energy geodesic steps,
adjoint/difference slope checks, profile-aware line acceptance, fine-grid
selection, and exact coefficient save/reload. Its 684-variable pilot improved
the critical norm and forward scale motion with close timestep and
cross-resolution agreement, but profile drift stalled near `6.3`, so it also
failed the unchanged promotion gate. Deterministic multi-start search and a
differentiable, amplitude- and scale-normalized spectrum-shape term are now
implemented and checked through the full adjoint. The first four-start pilot
cut strict drift from about `11.12` to `5.84`, but the improvement plateaued
far above one and L3 still fell. Fixed-threshold profile-window scheduling now
prevents accumulated crossing overshoot from changing coarse/fine window
counts. The next useful adjoint extension is a path-dependent smooth shape
cost, followed by staged `K_seed=3` multi-start screening. Only a state whose
strict drift approaches one on both grids should reach the independent FFTW
and higher-grid gates.

A credible path from this scaffold to a theorem has several hard gates:

1. **Numerical credibility:** reproduce standard benchmarks with both complete
   trajectory implementations; run convergence studies across cutoff, time
   step, box, and precision; and search for stable rescaled profiles rather
   than isolated spikes.
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
- D. Ayala and B. Protas, *Extreme vortex states and the growth of enstrophy
  in three-dimensional incompressible flows*:
  <https://arxiv.org/abs/1605.05742>
- D. Kang, D. Yun, and B. Protas, *Maximum amplification of enstrophy in
  three-dimensional Navier-Stokes flows*:
  <https://arxiv.org/abs/1909.00041>
- R. Suaza Jaque and O. Velasco Fuentes, *Reconnection of orthogonal
  cylindrical vortices*:
  <https://doi.org/10.1016/j.euromechflu.2016.11.001>
- F. Waleffe, *The nature of triad interactions in homogeneous turbulence*:
  <https://pubs.aip.org/aip/pof/article/4/2/350/402478/The-nature-of-triad-interactions-in-homogeneous>
- H. K. Moffatt and Y. Kimura, *Towards a finite-time singularity of the
  Navier-Stokes equations. Part 3. Maximal vorticity amplification*:
  <https://doi.org/10.1017/jfm.2023.472>
- S. Palasek, *Arbitrary norm growth in the 3D Navier-Stokes equations*:
  <https://arxiv.org/abs/2509.18595>
- T. Hou, Q. Wang, and D. Yang, computer-assisted weak non-uniqueness from
  singular initial data: <https://arxiv.org/abs/2509.25116>
- Z. Grujic, work on geometric depletion of vortex stretching:
  <https://arxiv.org/abs/2607.08866>
- L. Escauriaza, G. Seregin, and V. Sverak, the endpoint L3 regularity
  criterion:
  <https://www.pdmi.ras.ru/~seregin/Recent%20Publications/engESS3.pdf>
