# Nonlinear axis pilot: local checks pass, matching remains open

This experiment advances from manufactured test fields to **finite Taylor
approximations of the nonlinear inner equations (B.15)**. It does not construct
the complete inner/annular/outer profile of Theorem 4.6. The distinction is
essential: a local leading-order solution is not the full Navier-Stokes flow.

Status: `nonlinear-axis-pilot-passed-direct-join-failed`. No new candidate or
proof is promoted. The historical unforced `.14` resolution failure is unchanged.

## Input and deliberate scope limit

Source: [OpenAI, Finite Time Blowup for Navier-Stokes](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The equation pages, particularly 26-28 and 144-148, were read and the printed
nonlinear formulas visually checked. Parent commit: `e48a121`.

The pressure datum in the actual manuscript is determined by the completed
outer schedule through (A.21), **before** solving the axis problem. That
schedule and its quantitative parameter thresholds have not been constructed
here. Instead, this pilot explicitly fixes

\[
h=0.005,\quad j_0=0.05,\quad \sigma_*=0.002,\quad
\Lambda=65536,\qquad \Pi_0(\eta)=-\frac{13}{(1+\eta^2)^2}.
\]

This pressure is analytic, negative, even, and increasing for positive eta.
Those properties do not identify it with the paper's outer pressure. In
particular, this experiment does not verify the ordered global parameter
choices in (A.6)/(B.40) or the contraction thresholds in Proposition B.2.

The same finite parameters are used in every refinement. The protocol fixed
degrees 12, 18, 24; 60/90 decimal-digit arithmetic; eight eta centers; and
five radial samples from Y=0 through 4.1. Additional physical-check settings
were fixed before their evaluation. The first mixed-Cartesian stencil failed
through absorption of tiny terms; its replacement and retained negative
control are documented below. No scientific acceptance threshold was relaxed.

## What the construction actually solves

Following (B.1)-(B.3), set U*=4 eta+j0 and construct H*, W*, Z*, zeta*, chi.
Write Y=Lambda X and use the normalized unknowns in (B.12):

\[
\phi=\phi_*\Phi,\qquad U=U_*+u/\Lambda,\qquad
\Pi=\Pi_0+p/\Lambda,\qquad p_Y=g^2\Phi^2,\quad g=\phi_*/C.
\]

`inner.py` expands Phi, u and p as polynomials in Y with Taylor coefficients
in eta around each center. Their initial coefficients are Phi(0)=1, u(0)=0
and p(0)=0. At radial order n, the equations on page 147 give

\[
\Phi_{n+1}=
\frac{[-\chi\Phi+R_1/\Lambda]_n}{2(n+1)(n+2)},\qquad
u_{n+1}=\frac{[-Z_*/L+R_2/\Lambda]_n}{2(n+1)^2},\qquad
p_{n+1}=\frac{[g^2\Phi^2]_n}{n+1}.
\]

Here [.]_n means the coefficient of Y^n, not a sampled residual. The nonlinear
products and eta derivatives in both remainders are retained. One radial
step consumes one available eta derivative; the triangular Taylor budget
leaves three eta derivatives even at the final radial degree. No radial ODE
is solved independently at each eta while silently dropping eta transport.

The check is assembled separately from the **original** source equations
(4.9) and (B.15), not from the recurrence's R1/R2 expressions. Independently
re-solving at nearby eta centers also tests the propagated eta derivatives.
These checks do not bound the gaps between the eight centers or the infinite
Taylor tail. Stored local coefficient jets are not one certified global
profile over [-1,1].

## Amplitude and arithmetic

The azimuthal axis datum can be exponentially large, so C cannot be selected
by looking only at its real-axis maximum. This pilot uses a conservative
complex-neighborhood bound, implementing one requirement of (B.16).

For rho=sigma/[4(D+4+48+4j0)], the rectangle |Re eta|<=1+rho,
|Im eta|<=rho lies in |eta|<2. The change in H* along a vertical segment
is at most sigma/4, so both |H* +/- i sigma| are at least 3 sigma/4.
Consequently |zeta*| is bounded by

\[
M=4(1+8h)\{2D+5(8+j_0)\}/\sigma_*^2.
\]

A path from zero inside this rectangle has length below two. Choosing
log C=Lambda(2M+1) therefore gives |g|<=exp(-Lambda) there. This conservative
choice is not an optimized practical amplitude, nor a proof that every
threshold C0(Lambda) has been met.

The code retains g with mpmath's arbitrary-exponent numbers, and records log g;
g is **not set to zero**. Ordinary double conversion would return zero.
Pressure-increment and moment checks divide out the known tiny factors before
comparison or quadrature. A unit absolute tolerance is not allowed to make
the tiny pressure-balance equation pass automatically.

## Recorded numerical results

See [evidence.json](evidence.json) and [protocol.json](protocol.json). The
source-bound run passed all **13 numerical gates**.
All **36 regression tests** pass. The prior building-block audit was also
rerun unchanged: 29 tests and 25 audit gates still pass.

| Check | Worst recorded value, rounded | Limit |
|---|---:|---:|
| Original leading equations, including relative pressure balance | 2.08e-28 | <=1e-12 |
| Degree 18 versus 24, values and selected derivatives | 7.66e-27 | <=1e-12 |
| 60 versus 90 digits | 9.82e-49 | <=1e-30 |
| Minimum sampled Phi | 0.271048 | >=0.25 |
| Independently re-solved eta derivatives | 1.77e-18 | <=1e-10 |
| Full physical finite-difference discrepancy | 2.23e-12 | <=1e-6 |
| Coarse-to-fine physical stencil error reduction | 256.00x | >=4x |
| Five normalized moments versus independent quadrature | 3.75e-91 | <=1e-40 |

The errors are floating-point numerical comparisons, not interval bounds.
Very small local equation residuals are **not** a percentage of progress to
a global proof. The degree comparison includes values, first/second radial
derivatives, and the first three eta derivatives of Phi and u. The explicit
pressure balance is checked separately.

Omitting eta transport gives a normalized gap of 1.0; omitting the axial
pressure terms gives a gap of about 0.965. The nonlinear corrections are
nonzero: the result is not just the scalar Bessel comparison profile.

## Full equations and the retained stencil failure

`physical.py` maps the computed fields to Cartesian space and retains all
time, advection, pressure and viscosity terms, including axial diffusion.
At Y=2 and q=.2, the sampled component-normalized full residuals are:

| eta | Radial | Azimuthal | Axial |
|---|---:|---:|---:|
| Zero of H* (about -0.0111233) | 1.05372 | 1.00000 | 0.24529 |
| 0 | 1.05370 | 1.00000 | 0.24852 |
| .5 | 0.52549 | 1.00009 | 0.43021 |

Each component is divided by the largest magnitude of its individual
equation terms, with no unit floor. The full residual is substantial. Only
the leading tangential residual, which removes axial diffusion as specified
in (4.12), is small. The omitted higher-order construction is still needed
to establish a smoothly extending force.

The initially attempted Cartesian finite-difference check sampled mixed
radial/swirl components away from theta=0. The swirl is so tiny that adding
it to a finite radial component loses it even at 90-digit precision. The
resulting relative azimuthal error was enormous (above 10^(2.4e12)); this
is numerical absorption, **not physical blowup**. That implementation remains
as `cartesian_finite_difference` and a regression test requires its failure.

The accepted crosscheck instead uses scalar cylindrical components on the
theta=0 ray, with independent physical-coordinate stencils and the full
curvature terms. It agrees with Cartesian derivative arithmetic, and its
error decreases approximately fourth-order. Both implementations share the
finite input profile and implicit coordinate map; they are not independent
reproductions of the complete manuscript.

## Matching: calculated data, not a completed connection

All five inner moments (M,I,J,S,Cp) in (4.15) are computed from the nonlinear
polynomial. Their normalized values at eta=.5, Y=4 are independently checked
by quadrature. This verifies their evaluation; **it does not match them to
an exterior**.

All eight tested endpoints fail necessary tests for direct attachment to
the heat exterior. For example, at eta=.5:

- U is about 2.05075, whereas the heat exterior has U=0.
- The normalized M mismatch is about 1.00018; the target construction has
  zero total axial moment.
- The logarithmic swirl slope is about -1.19512, outside the heat exterior's
  necessary interval [-0.505,-0.5]. Fitting its amplitude cannot fix its slope.

This expected failure is **not a flaw in the paper**: it explicitly inserts
an annulus to change the flow and preserve the moments. This pilot has not
implemented the B.5-B.8 collar, the five bump corrections, or their stress
cone conditions. It also has not constructed the A.2 outer schedule that
must provide the correct axis pressure in the first place.

## Reproduce and review

Install `requirements.txt`, then from the repository root:

```sh
python3 research/navier_stokes_cascade/results/nonlinear_axis_pilot/test_pilot.py -v
python3 research/navier_stokes_cascade/results/nonlinear_axis_pilot/run.py \
  --verify-record research/navier_stokes_cascade/results/nonlinear_axis_pilot/evidence.json \
  --paper /path/to/navier-stokes.pdf \
  --output /path/to/unused/nonlinear-axis-recheck.json
```

Use an unused output path. Without `--paper`, the run explicitly records that
it did not re-read the PDF bytes. CI recomputes the pilot against the pinned
source hashes and checks the stored evidence; it does not download a mutable
external PDF. Numerical source inputs and tolerances are checked separately
from the unmatched scientific status. The original solvers and archives
are not modified.

Next: construct and test the actual outer schedule and its pressure integral,
then rerun this inner solver with that datum. Only then build the annular
connection with five-moment restoration and stress-cone margins. The present
module is reusable for that work but is not the completed matched profile.
