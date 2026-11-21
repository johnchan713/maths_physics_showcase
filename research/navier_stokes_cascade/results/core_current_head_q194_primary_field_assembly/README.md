# Q194: complete primary field bounds and support assembly

Parent: Q193, `c58b5f2d5add7e2da6ed366f6f965e9f473ff836`.
The live branch, Q192--Q193, the physical-pulse transfer/foundation records,
their frozen evidence, and active CI were read first. Q189's selected profile,
full pressure memory, supports, frequency and unchanged `C^-20` activation
width remain fixed. Q189 CI has completed successfully; later runs were still
in progress when this checkpoint began.

Q193 closed the uncut homogeneous Duhamel and phase-partition recurrences.
This checkpoint places the actual cutoffs on potentials and pressures,
retains every cylindrical curl component and derives the complete primary
velocity, pressure and potential bounds on the same domain. The result remains
conditional on the named profile and frame inputs.

## 1. Exact zero-order rescalings

With

\[
 A={1\over2}+h,\qquad s_x={1\over2}+{h\over2},
 \qquad s_t=1+{3h\over2},
\]

the primary coefficient powers give

\[
 -A+{h\over2}=-s_x,qquad
 -2A+h=-2s_x,qquad
 {1\over2}-A+h=0.                                  \tag{1}
\]

These are respectively the complete velocity, pressure and fixed vector
potential zero-order powers before requested derivatives are added.

## 2. Cut before taking the complete curl

The transverse, pulse and shell cutoffs are applied to the vector potential
and wave pressure before differentiation. In cylindrical components,

\[
 u_r={1\over r}\partial_\theta A_z-\partial_zA_\theta,
\]
\[
 u_\theta=\partial_zA_r-\partial_rA_z,
\]
\[
 u_z=\partial_rA_\theta+{A_\theta\over r}
      -{1\over r}\partial_\theta A_r.              \tag{2}
\]

Thus (2) retains all three velocity components, seven component terms, the
axial velocity and both explicit `1/r` connections. Since the support stays
on a fixed annulus with `r=Q^(1/2)R` and `R>=R_*>0`, each connection costs at
most `Q^-1/2`, no worse than one `s_x` spatial slot. Repeated Cartesian/frame
conversion differentiates these factors but gives the same per-slot bound.

## 3. Every product-rule allocation fits the selected budget

A spatial derivative may hit the phase, a coefficient/cutoff, or a
cylindrical/frame connection. Their respective worst losses are

\[
 s_x,qquad {1\over2}+h\kappa,qquad {1\over2}.
\]

A time derivative may hit the phase or a coefficient/cutoff, costing

\[
 s_t,qquad 1+h.
\]

Because `0<kappa<1/2`,

\[
 {1\over2}+h\kappa<s_x,qquad {1\over2}<s_x,
 \qquad 1+h<s_t.                                   \tag{3}
\]

In any Leibniz term the derivative slots are merely distributed among these
possibilities. Summing (3) over the allocation proves that `a` spatial and
`b` time derivatives cost at most `a s_x+b s_t`. Q193 supplies the finite
polynomial `S` loss for derivatives of the uncut amplitude; cutoff and compact
frame derivatives only increase that finite power.

## 4. Conditional complete primary bounds

Combining (1)--(3), for every fixed `(a,b)` and every accepted band,

\[
 |\nabla^a\partial_t^b u_p|
 \le C_{a,b}Q^{-((a+1)s_x+bs_t)}S^{B_{a,b}}
       \sqrt\zeta\,\delta^{-B_{a,b}}P,              \tag{4}
\]

\[
 |\nabla^a\partial_t^b p_p|
 \le C_{a,b}Q^{-((a+2)s_x+bs_t)}S^{B_{a,b}}
       \sqrt\zeta\,\delta^{-B_{a,b}}P,              \tag{5}
\]

\[
 |\nabla^a\partial_t^b {\cal A}_p|
 \le C_{a,b}Q^{-(as_x+bs_t)}S^{B_{a,b}}
       \sqrt\zeta\,\delta^{-B_{a,b}}P.              \tag{6}
\]

Constants may depend on fixed derivative/harmonic order and the frozen
geometry. They are uniform in accepted band, label and copy. Equations
(4)--(6) are conditional analytical estimates, not absolute smallness as
`Q` tends to zero and not a numerical PDE verification.

The pressure denominator and moving-normal term remain those audited in
Q191--Q192. The source-normal term vanishes only for this homogeneous primary
pulse; it remains mandatory for sourced inverses.

## 5. One support domain and endpoint memory

Retain the actual supports

\[
 |\xi|\le r_0/2,qquad |v/L_s-1/2|\le2/5,
\]

with the pulse cutoff equal to one on the central fifth, and the fixed closed
active shell. Differentiation and curl do not enlarge these supports. Closed
enlarged preimages retain at most one active label; local finiteness holds for
`q>0`. Smoothness at the axis follows from zero support near the axis, not
from extending fractional radial powers through it.

The flat shell factor `sqrt(zeta)` and the explicit transverse/pulse cutoffs
give smooth zero extension at their edges. The positive Gaussian envelope
`P` is not used as endpoint support. At `T=0` the chart jets remain one-sided;
no temporal zero extension is asserted there. Generated sources still require
their own full-harmonic pulse-endpoint extension.

The lower bound for the uncut radial amplitude is not promoted to a positive
pointwise bound for the oscillating cutoff velocity.

## 6. Remaining boundary

The next obligation is the complete-wave covariance/mean calculation using
the assembled radial, angular and axial velocity and wave pressure. It must
retain all six covariance entries, the moving pressure response, and the
radial/axial mean corrections on one common domain. Generated sources and the
complete viscous/cutoff/interaction residual remain later obligations.

No unconditional pulse hierarchy, complete residual improvement, admissible
force, comparison result or Navier--Stokes blow-up is claimed. Tests verify
the exact component and exponent ledger; they are not a continuum proof.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q194_primary_field_assembly/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q194_primary_field_assembly/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q194_primary_field_assembly/evidence.json \
  --output q194-fresh.json
```
