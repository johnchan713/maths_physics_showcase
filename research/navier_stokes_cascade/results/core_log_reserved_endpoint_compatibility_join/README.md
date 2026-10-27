# Log-reserved endpoint compatibility join

Parent: `254ed16c56eb8b2e707e2bb689c48f331dca659e`.
Status: `log-reserved-corner-and-punctured-endpoint-inputs-joined-conditionally`.

The parent transfers every fixed jet of the complete six-group residual to
Cartesian terminal-corner flatness.  The inherited punctured endpoint
checkpoint separately supplies compatible one-sided jets of the complete
potential, direct swirl, velocity and fixed-gauge pressure away from the
singular origin.  This checkpoint joins those two inputs on the same field
and records the exact finite derivative requests needed by localization.

The result is conditional on the frozen all-order field and pressure
contracts.  Exact-arithmetic tests check the derivative ledger, stage
quantifiers and provenance; they do not numerically construct the comparison
state or prove the continuum estimates.

## 1. Keep the two endpoint roles distinct

For each punctured terminal compact `|x|>=d>0`, the inherited enlarged atlas
and axial/equatorial cover give, conditionally,

\[
 \sup_{E_+\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p)|<\infty                 \tag{1}
\]

at every fixed order.  Here `A` is the complete potential, `B` is the
separate axisymmetric direct swirl, and `p` retains its full fixed gauge,
both exceptional pressure tails, all changing covariance pressures and the
axial/moving-moment memory.  The heat exterior supplies its exact mixed
swirl and pressure jets.

At the excluded origin, the new parent supplies a different fact:

\[
 |D_x^\alpha\partial_t^j\mathcal R_1(u_*,p_*)|
 \le C_{\alpha,j,N}
 (\tau+|z|^{1/D})^{N}                                   \tag{2}
\]

for every fixed `alpha,j,N`.  Equation (2) regularizes residual and later
force traces.  It does not assert that velocity or pressure is regular at
the terminal origin.  Equation (1) cannot replace (2), and pathwise residual
flatness cannot replace its compact-`X` uniformity.

## 2. Exact derivative ledger

Let a requested force derivative contain `a` spatial and `j` time
derivatives and put `m=a+j`.  Differentiating the exact localization formula
at order `m` needs

\[
 A:m+3,\qquad B:m+2,\qquad p:m+1.                       \tag{3}
\]

Uniform one-sided Cauchy convergence reserves one further time derivative:

\[
 A:m+4,\qquad B:m+3,\qquad p:m+2.                       \tag{4}
\]

These are exactly the `potential_for_future_residual_Cauchy`,
`direct_for_future_residual_Cauchy`, and
`pressure_for_future_residual_Cauchy` orders in the punctured endpoint
ledger.  For the residual at the corner, the corresponding Cauchy request is
the fixed finite jet `m+1` in (2).

Thus the extra time derivative does not consume a stage-independent reserve.
For each fixed `(a,j,N)`, it merely raises the finite request order and may
raise the finite comparison stage.  The field, profile, gauge, labels,
cutoff sequence and localization cylinder stay unchanged.  A
request-dependent estimate is not a request-dependent field.

For example, at `a=j=N=0`, the punctured request is

\[
 A:4,\qquad B:3,\qquad p:2,                             \tag{5}
\]

and the corner request uses residual jet one.  Higher examples in the frozen
evidence verify the same offsets through total mixed order 200.

## 3. One compatible closed-time input package

The retained cover is now:

* punctured compacts: uniform one-sided `A,B,p` jets and their next time
  derivatives;
* the origin: all fixed jets of the complete residual are Cartesian flat;
* the heat exterior: exact mixed swirl and fixed-gauge pressure jets;
* overlaps: the inherited axial/equatorial cover inside the enlarged closed
  atlas.

All entries refer to the same locally finite field and selected cutoff
sequence.  The residual side retains time, unit viscosity, pressure, both
interaction orientations and the correction square.  The field side retains
complete pressure memory and the direct swirl separately from the curl.

This package is precisely the input needed to replay the twelve-term
axisymmetric localization identity.  That identity and the resulting force
extension are deliberately left to the next checkpoint.

## 4. Scope and next obligation

This checkpoint closes the derivative-order and quantifier join between the
punctured endpoint and the corrected corner residual, conditional on their
frozen continuum contracts.  It does not prove terminal regularity of
velocity or pressure at the origin, localize the field, construct an
admissible force, replay energy or comparison, or prove blow-up.

The next obligation is to replay

\[
 u=\operatorname{curl}(cA)+cB,
 \qquad p=cp_0,                                         \tag{6}
\]

with all twelve residual terms, unit viscosity, the complete pressure cutoff
term, every cutoff derivative, both nonlinear crosses and the curl square.
Only after smooth terminal force jets are proved should the future Borel
force extension, energy estimate and comparison be advanced.

Reproduce from the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_log_reserved_endpoint_compatibility_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_log_reserved_endpoint_compatibility_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_log_reserved_endpoint_compatibility_join/evidence.json \
  --output core_log_reserved_endpoint_compatibility_join_fresh.json
```
