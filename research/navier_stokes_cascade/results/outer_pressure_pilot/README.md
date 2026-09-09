# Scheduled outer pressure coupled to the inner equations

This checkpoint replaces the earlier freely chosen pressure seed with a
numerical evaluation of the complete unedited swirl schedule in Appendix A.2.
It also solves the two actual angular correction equations (A.11), and uses
the resulting pressure datum in the nonlinear inner equations (B.15).

**The global matched profile remains incomplete.** In particular, computing
the swirl schedule does not solve the axial pulse's three moment conditions,
verify the outer stress cone, or attach the regular axis through an annulus.
The original unforced candidate's `.14` resolution failure is unchanged.

Source: [the supplied Navier-Stokes manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Printed pages 128-134 were read; the schedule, correction and pressure
formulas on pages 129, 130 and 133 were checked visually. Parent checkpoint:
`bbd3e397fbe0306d4e474d6f736548ce7a17bc06`.

## Why the pressure can be constructed first

Write `y=log(X/XR)` and let E be the swirl profile. Pressure normalized to
vanish at radial infinity gives the axis datum

\[
\Pi_0(\eta)=-\frac12\int_{-\infty}^{\infty}E(y,\eta)^2\,dy.
\]

The axial pulse changes U, not E. The two angular bumps in (A.11) are chosen
to change the angular moment while preserving this pressure integral.
Consequently Lemma A.5 explicitly defines the pressure with those bumps
omitted. This is the datum constructed here; no substitute pressure shape
is selected to make the inner equations easier.

The terminal waiting length is determined using the angular-moment target
in (A.11). This would be an assumption if that correction were not supplied.
Here its actual discrepancy is computed and both bump equations are solved
and independently integrated. Later heat and axis modifications still have
to preserve the same pressure datum.

## The finite schedule

`schedule.py` implements the specified flat smooth step and all thirteen
stages after the ideal reference interval: the first slope transition,
axial reduction, intermediate slope transition, reserved power interval,
pulse swirl, parameter interpolation, angular correction interval, two
exterior slope transitions, steep decay hold, terminal wait, terminal
cutoff, and infinite power tail.

The finite pilot fixes

\[
\begin{aligned}
M_d&=4,&T_d&=e^4+10, &\log P_*&=T_d+1,\\
\lambda&=10^{-5},&h&=e^{-2T_d},&T_f&=200,\\
c_o&=10^{-3},&j_0&=.05,&\sigma_*&=.002,\\
\Lambda&=65536P_*^2.
\end{aligned}
\]

Lowercase lambda controls the outer schedule. Uppercase Lambda controls the
inner radial scale. They are different parameters.

These choices satisfy the explicit inequalities `P*>exp(Td)` and
`0<h<min(lambda,exp(-Td))`, and leave room for the four reserved patches.
**They do not certify the unspecified sufficiently large/small thresholds.**
In particular, `Md=4` is a finite computational choice, not a demonstrated
small bound on the axial-transition derivative. Raising XR afterward cannot
repair a failed prerequisite at an earlier stage.

## How small contributions are treated

No physical radius of size `exp(y)` is needed. The code stores log amplitudes
and integrates each stage relative to its own starting amplitude. Constant
slopes, the ideal region and the infinite power tail have elementary
integrals. Smooth transitions use positive composite Gauss-Legendre rules
with orders 16, 24 and 32 on four panels.

The pressure takes the form

\[
\frac{\Pi_0}{P_*^2}
=-\frac12\sum_j w_j(1+\eta^2)^{-2\theta_j},
\qquad w_j>0,\quad 0\leq\theta_j\leq1.
\]

Each stage retains its separate exponential scale. The ideal region
contributes exactly `-(5/2) f^2`, with `f=(1+eta^2)^(-1)`. The remaining
positive masses make pressure strictly more negative. The representation
also explains evenness and the sign of its first derivative: differentiating
each summand gives a positive multiple of eta. Its analytic continuation
uses a branch of `log(1+eta^2)` away from the poles, as in Lemma A.5.

Late stages are much smaller than the precision of the combined pressure.
They are recorded separately, with a uniform bound on their relative
contribution over `[-1,1]`. A zero tail fails validation and reproduction.
This bookkeeping does not turn numerical quadrature errors into rigorous
interval bounds.

The pressure Taylor coefficients come from the differential identity

\[
(1+\eta^2)\frac{d}{d\eta}(1+\eta^2)^{-b}
=-2b\eta(1+\eta^2)^{-b}.
\]

Thus parameter derivatives are passed to the inner solver; they are not
discarded by solving unrelated radial equations at sampled eta values.

## The actual two-bump correction

Set `beta=1-lambda` and `rI=I/(XH)`. On an unchanged power interval,

\[
\delta=r_I-1/\beta,\qquad \delta'=-\beta\delta.
\]

The code propagates delta separately. Subtracting two rounded values near
`1/beta` after the long pulse would incorrectly erase its nonzero memory.
The interpolation contribution is integrated with a positive, exponentially
weighted kernel. This gives the actual angular discrepancy at the first
correction center.

The correction is `E -> E(1+c1 beta1+c2 beta2)`, with the two disjoint bumps
of width .3 and separation 2 specified on page 130. With radial weights
included, the angular equation is linear, while the pressure equation is
linear plus quadratic in the coefficients. Both are solved by the iteration
in Lemma A.2. The report includes the matrix, quadratic terms, coefficients,
residuals normalized by the size of the edit, and a conservative smallness
estimate. These estimates concern this correction only.

Independent tanh-sinh integration checks the edited moment functionals.
A direct integral of H checks the actual interpolation discrepancy at
eta=.5. Its negligible initial-memory term has an explicit bound. Ignoring
the pressure equation is retained as a failing control.

The first four-panel calculation of the concentrated interpolation kernel
failed this independent check: its relative discrepancy was about 8.96e-11,
above the fixed 1e-12 limit. The accepted calculation uses eight panels for
that kernel. The original failure is recomputed and retained alongside it;
the scientific limit is unchanged.

This is a **two-moment outer correction**, not the five-moment restoration
needed for the still-unbuilt axis annulus.

## Coupling and the rejected arithmetic setting

`coupling.py` expands the original nonlinear W, H, U and Pi products, instead
of copying the previous pilot's R1/R2 coefficient assembly. On the old seed,
the two formulations agree. The shared low-level Taylor operations and
coordinate map are explicitly reused; this is not an independent
implementation of the entire manuscript.

The first 60-digit probe resolved normalized fields but failed to resolve
the absolute value of `log(g)`, whose magnitude is about `5.14e69`. A small
relative error in that logarithm becomes an enormous multiplicative error
after exponentiation. The failed setting is retained in the evidence.

The accepted refinements use 120 and 160 digits and compare
`abs(expm1(log_g_low-log_g_high))`, alongside the field values and derivatives.
The physical finite-difference check uses 180 digits because taking two
differences across the very small physical scales requires extra precision.
No scientific residual tolerance was relaxed.

The full Cartesian derivative calculation includes time, advection,
pressure, radial viscosity and axial viscosity. The scheduled pressure
adds its analytic axial pressure gradient to the previous field evaluator.
Component-separated physical finite differences check that gradient and
the full residual independently of derivative arithmetic. They share the
finite input profile and coordinate map.

## Evidence and reproduction

All **27 regression tests** pass. The refined run passes **22 numerical gates**.
The measurements below are
floating-point comparisons on the declared samples, not certified bounds
on an infinite series or a continuum of eta values.

| Check | Worst recorded value, rounded | Limit |
|---|---:|---:|
| Pressure stages, quadrature refinement | 8.12e-17 | 1e-14 |
| Independent schedule integrals | 3.22e-16 | 1e-12 |
| Independently integrated edited moments | 1.66e-25 | 1e-12 |
| Original nonlinear leading equations | 2.07e-28 | 1e-10 |
| Radial degrees 18 versus 24 | 7.76e-28 | 1e-10 |
| 120 versus 160 digits, including amplitude | 1.45e-50 | 1e-30 |
| Minimum sampled normalized swirl Phi | 0.271104 | at least 0.25 |
| Full physical finite-difference discrepancy | 2.20e-12 | 1e-6 |

The computed pressure is dominated by
`-3.31462273001433255 P*^2/(1+eta^2)^2`. The separately recorded remaining
stages have a uniform relative mass bound below `9.86e-564924`, subject to
the numerical evaluation of their integrals. That tiny positive remainder
is not a reason to omit its construction or erase it from the record.

The largest relative angular edit is bounded by approximately `4.34e-175`.
The four-panel interpolation calculation had relative error `8.96e-11`;
eight panels reduce it to `8.89e-26`. Independent edited-moment integration
then has worst error `1.66e-25`. The two retained failures are detected, and
the accepted refinements keep the original scientific tolerances.

At eta=.5, Y=2 and q=.2, the full momentum residual divided componentwise
by the largest individual term is approximately `(1, 1, 1.81e-57)` in
radial, azimuthal and axial order. Thus the radial and azimuthal full
residuals remain substantial. This is not a new Navier-Stokes solution.
The physical stencil's coarse-to-fine error decreases about 256-fold.

`protocol.json` fixes the samples, numerical refinements, limits and scope.
`evidence.json` contains every stage's normalized pressure jets, the two-bump
solutions, inner radial refinements, full residuals and rejected arithmetic
control. `audit.py` binds the evidence to the source hashes, optionally checks
the supplied PDF bytes, recomputes the experiment, and rejects unsupported
scientific promotion. Prior audit modules and their evidence are unchanged.

Install `requirements.txt`, then run from the repository root:

```sh
python3 research/navier_stokes_cascade/results/outer_pressure_pilot/test_outer.py -v
python3 research/navier_stokes_cascade/results/outer_pressure_pilot/audit.py \
  --verify-record research/navier_stokes_cascade/results/outer_pressure_pilot/evidence.json \
  --paper /path/to/navier-stokes.pdf \
  --output /path/to/unused/outer-pressure-recheck.json
```

Existing output files are never overwritten. CI runs the same reproduction
without downloading a mutable external PDF; that distinction is recorded.

## Remaining work

The next outer construction must choose the axial pulse and its end bumps
to close M, J and S, then check stress-cone margins over the whole schedule.
Those checks may require new finite parameter choices and recomputation of
this pressure. The heat replacement needs its own three-moment compensation.
The axis annulus then needs five-moment restoration, followed by uniform
analytic estimates and the all-order momentum correction and forcing
construction. Small sampled leading residuals do not settle any of those
remaining requirements.
