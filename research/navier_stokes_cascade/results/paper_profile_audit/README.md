# Separate paper building-block audit

Status: **building-block-audit-passed**. This is not a reproduction or
verification of the complete blowup construction, and supplies no new
candidate to the existing unforced search. Its `.14` resolution failure
remains unchanged.

## Frozen input and scope

The input is OpenAI's [Finite Time Blowup for Navier-Stokes](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
166 pages, SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The same hash is in [protocol.json](protocol.json); local PDF bytes were
verified when generating [evidence.json](evidence.json).
Equation pages 25, 26, 27, 138 and 147 were visually inspected. The PDF is
referenced, not redistributed. Changed input bytes require a new audit.

Smooth forcing is allowed in Clay breakdown alternatives C and D, subject to
their hypotheses; its presence alone does not disqualify a construction.
See the [official statement, page 2](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf#page=2).
We do not claim to have checked all those hypotheses for this manuscript.

The protocol fixed fixtures, sample locations, main tolerances and negative
controls before the reported audit. The parent research checkpoint is
`3386bd8223d56b2c7dac1fee0f8a998ddd93f78d`.

| Implemented object | Source locator | What was checked |
|---|---|---|
| Similarity map and derivatives | (4.2), (4.5) | Implicit-root derivatives against profile-coordinate formulas |
| Divergence-free ansatz and integrated stress | (4.7), (4.11), (4.12), (4.16) | Manufactured polynomials, not the asserted nonlinear profile |
| Scalar inner comparison function | (B.11)-(B.13) | Rational recurrence, positivity lower bound, Bessel comparison |
| Viscous swirl exterior | (A.32)-(A.37) | Integral and special-function evaluations; full momentum residual away from the axis |

## Results

The evidence records Python, NumPy 2.3.5 and SciPy 1.17.0. There are 29
regression tests and 25 audit gates: 54 manufactured samples, 9 axis checks,
two three-step finite-difference refinements, 10 heat profile samples with
four derivatives each, and 6 exterior momentum checks.

| Measurement | Worst recorded value, rounded | Frozen gate |
|---|---:|---:|
| Cartesian versus cylindrical residual discrepancy | `9.32e-16` | `<=1e-10` |
| Integrated-stress identity discrepancy | `4.79e-15` | `<=1e-10` |
| Finest physical finite-difference discrepancy | `2.72e-9` | `<=1e-6` |
| Coarse-to-fine error reduction | `251.27x` minimum | `>=4x` |
| Heat integral versus special-function discrepancy | `9.04e-13` | `<=1e-9` |
| Exterior momentum residual / term scale | `8.52e-16` | `<=1e-9` |
| Scalar series versus Bessel reference | `1.11e-16` | `<=1e-12` |

Identity discrepancies divide the maximum absolute difference by the larger
of one and the compared magnitudes. Full-residual and finite-difference
errors use the largest individual time, advection, viscosity or pressure
term, floored at one. These are numerical conventions, not certified bounds.

The manufactured fields have **nonzero** normalized full residuals from
`0.852` to `2.753`. A passing identity does not make them zero-force
Navier-Stokes solutions. The exterior is singular at the axis and independent
of the axial coordinate; alone it is not a regular, finite-energy
whole-space candidate.

## Checker design and limitations

`profiles.py` computes Cartesian gradients and Hessians with derivative
arithmetic. A separate cylindrical assembly retains curvature and axial
viscosity. Fourth-order physical-coordinate stencils use only scalar field
values. The analytic assemblies share the implicit coordinate map and basic
arithmetic; this is not a fully independent implementation of the manuscript.

At viscosity one our convention is

\[
R=\partial_tu+(u\cdot\nabla)u+\nabla p-\Delta u.
\]

This is the force required by a prescribed field. Computing it does not show
it extends smoothly through a proposed singular time. Leading tangential
stress identities omit axial diffusion; the full checker does not. Deliberate
omission of axial viscosity, omission of swirl curvature, and reversal of
the stress sign produce maximum normalized gaps `1.754`, `0.781` and `2.000`.
Each exceeds the detection threshold `1e-4`.

The tests reject missing checks, changed tolerances, nonfinite evidence,
altered physical samples and unsupported proof/candidate labels. The scalar
comparison's rational lower bound is `305719/1152000 > .265` on the declared
interval. Its alternating-series truncation bound excludes floating-point
roundoff; a finite polynomial truncation is not the exact ODE solution.

The exterior integral is checked against SciPy's confluent hypergeometric
function on the frozen grid. A separate radial pressure integral and finite
difference check the pressure gradient. That test shares the integral
evaluator for the heat profile: nesting the special-function reference inside
adaptive pressure quadrature was too noisy for this derivative tolerance.
Quadrature warnings are errors, and error estimates are not rigorous bounds.

A practical scaling inference: the nominal axial/radial diffusion factor is
`tau^(2h)`. At the illustrative fixture `h=.005`, it is about `.9705` at
`tau=.05` and `.8710` at `tau=1e-6`; reaching `.1` requires `tau=1e-100`.
This is a scale factor, not a measured residual ratio or a parameter choice
certified for the complete construction. It motivates retaining the full
equation at computationally accessible scales.

## Reproduce without overwriting evidence

From the repository root, with NumPy and SciPy installed:

```sh
python3 research/navier_stokes_cascade/results/paper_profile_audit/test_profiles.py -v
python3 research/navier_stokes_cascade/results/paper_profile_audit/run_audit.py \
  --verify-record research/navier_stokes_cascade/results/paper_profile_audit/evidence.json \
  --output /tmp/paper-profile-new-run.json
```

Choose an unused output path. Add `--paper /path/to/navier-stokes.pdf` to
verify the PDF bytes. Otherwise the output explicitly records that the PDF
was not re-read. CI does not fetch a mutable external PDF: it checks source
hashes, recomputes every gate and compares physical values with the stored
record. Cancellation-error ratios are independently gated, not required to
be bit-for-bit equal. Dependency versions are recorded per run.

The earlier damaged local checkpoint chunk is neither an input nor a changed
publication artifact. Historical solvers and archived remote states are
untouched by this track.

## Next milestone and open obligations

Next: construct an explicit numerical nonlinear inner/annular profile under
documented parameter choices, test matching conditions, and feed it into
this full-residual checker. Small leading stress alone must not promote a
search candidate.

The protocol's obligation list remains open: valid ordered parameter bounds;
nonlinear profile and five-moment matching; admissible stress and oscillatory
corrections; stage-uniform residual improvement; infinite correction sums
and all required derivatives; localization, smooth force extension and
uniqueness comparison. This milestone verifies none of that list in full.
It is a tested foundation, not a percentage-of-proof estimate.
