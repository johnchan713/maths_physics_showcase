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

`GalerkinSystem` evolves a mean-zero, real, divergence-free velocity field on
the 2-pi periodic torus using the Fourier cube
`|k_x|, |k_y|, |k_z| <= K`. For each retained non-zero mode,

```text
d u_k / dt = -nu |k|^2 u_k
             - i P_k sum_(p+q=k) (u_p . q) u_q,
```

where `P_k` is the Leray projector. Time stepping is classical RK4. The direct
convolution is intentionally small and auditable; it is not intended as a
production turbulence solver.

The CSV diagnostics include:

- normalized kinetic energy, enstrophy, and palinstrophy;
- a sampled scale-critical L3 velocity norm;
- sampled maximum vorticity and its time integral (a BKM-inspired diagnostic);
- a rigorous-for-the-truncated-polynomial Fourier upper bound on maximum
  vorticity;
- spectral centroid and energy fraction touching the cutoff shell;
- divergence, Fourier-reality, and semi-discrete energy-balance defects.

The tests check the structural identities before any experiment is interpreted:
Leray projection, real/divergence-free invariants, nonlinear energy
cancellation, viscous energy dissipation, and short-run numerical stability.
They also encode the elementary guardrail that zero gradient does not imply
zero field value, a mistake found in some purported proofs.

## Build and run

From the repository root:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
./build/research/navier_stokes_cascade/navier_stokes_cascade
```

For a slightly larger exploratory run:

```bash
./build/research/navier_stokes_cascade/navier_stokes_cascade \
  --cutoff 3 --viscosity 0.02 --dt 0.0002 --steps 1000 \
  --diagnostic-every 25 --output cascade-k3.csv
```

Use `--help` for all parameters. Runtime grows quickly because the current
quadratic convolution is direct rather than FFT-based.

## Interpretation guardrails

Every fixed Galerkin cutoff is a smooth finite-dimensional ODE, so it cannot by
itself demonstrate PDE singularity. In particular:

- growth that changes when `K`, the grid, or `dt` changes is a resolution
  artifact until proved otherwise;
- cutoff-shell energy above roughly one percent is reported as an
  under-resolution warning, not evidence of blow-up;
- the sampled vorticity maximum is a lower estimate of the truncated field's
  true maximum, while the Fourier sum is an often-loose upper bound;
- the sampled BKM integral is finite-run telemetry, not the hypothesis or
  conclusion of a theorem;
- ordinary floating point cannot certify inequalities needed by a proof.

## Research gates

A credible path from this scaffold to a theorem has several hard gates:

1. **Numerical credibility:** replace direct convolution with a dealiased 3D
   FFT solver; reproduce standard benchmarks; run convergence studies across
   cutoff, time step, box, and precision; and search for stable rescaled
   profiles rather than isolated spikes.
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
