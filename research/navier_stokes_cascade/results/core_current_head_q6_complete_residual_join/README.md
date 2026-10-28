# Current-head Q6 complete six-group residual join

Parent: `e8142c9e6a31242a860ea1403bd80bd8b194a0f7`.
Status: `current-head-Q6-complete-six-group-local-residual-bound-conditionally`.

Q1--Q5 bind one leading core/root object, its positive later-tail regional
gaps, complete frame/path domain, one fixed compensator/projector choice, and
the shifted compact and fixed slow inverse maps with their actual
log-renormalized sources.  This checkpoint binds that exact object directly
to the frozen log-reserved residual replay.  It preserves the physical
characteristic shift, complete covariance and pressure memory, actual axial
velocity, supports, gauges and locally finite cutoff sequence.

The result is analytical relative to the Q1--Q5 continuum contracts.  The
exact tests freeze the request, exponent and provenance interfaces; they do
not evaluate the enormous selected field or replace the continuum argument.

## 1. Same-object finite-state residual

At a fixed comparison stage retain the complete physical split

\[
 \mathcal R_{\rm phys}
 =\mathcal R_{\rm regular,current}
  +(-A_0Q^\circ g_{r,{\rm final}},0,+A_1Q^\circ E_{z,{\rm post}})
  +T_{\rm Gaussian,old}+T_{\psi,{\rm signed}}+T_{\rm base,ext}.       \tag{1}
\]

The regular term includes the current source, signed linear term, both
old/signed wave-cross orientations, signed square, old-mean interaction and
new fast and slow feedback.  Equation (1) keeps all six covariance entries,
`g_d-A1 g_d`, the complete `partial_Z(J_z+c_rho P)`, `c_(rho,Z)P`, the old
pressure and both `Pi1/Pi2` tails.  The potential-curl and full-pressure
cutoff derivatives remain in the complete representatives.

The Q5 and residual evidence records point to the identical frozen
`core_current_log_source_reconciliation` record.  Thus the source passed to
(1) is not a parallel reconstruction with merely matching formulas.

## 2. Four accounts in three disjoint logarithmic classes

Put

\[
 \delta={1\over2}-3\kappa>0,\qquad B_j=B_0+j\delta,
 \qquad \epsilon=q^h.                                  \tag{2}
\]

After the jet `m` and target `N` are fixed, the recurrence returns finite raw
requests `R_A0,R_A1`.  The one additional compact integration is invoked at

\[
             R_{A0}+2,\qquad R_{A1}+2.                 \tag{3}
\]

Each of `A0` and `A1` owns a separate `h kappa` gain, spends half on its own
finite source logarithm and retains half.  If `J0` is the old pure-exponent
threshold, take `J=J0+1`; half of the new `h delta` gain pays the regular
residual logarithm.  The first omitted cutoff representative separately
spends half of its gain on its coefficient, derivative and cutoff logarithms.
No account pays another class.  The elementary estimate

\[
 (1+|\log q|)^Pq^\sigma
 \le \exp\!\left({P^2\over4\sigma}\right),\qquad \sigma>0,             \tag{4}
\]

is applied only after the finite log degree and its owning positive reserve
have been fixed.

## 3. Exact six-group comparison

For the same locally finite field write
`(e,pi)=(u_*,p_*)-(U^[J],P^[J])`.  The exact difference is

\[
 \mathcal R(u_*,p_*)-\mathcal R(U^{[J]},P^{[J]})
 =e_t-\Delta e+\nabla\pi
 +(U^{[J]}\!\cdot\nabla)e+(e\!\cdot\nabla)U^{[J]}
 +(e\!\cdot\nabla)e .                                  \tag{5}
\]

Thus time, unit viscosity, pressure, old-new, new-old and the new square are
all present.  With

\[
 \tau={g_{J+1}\over2}-L_{m+2},\qquad H_{m+2}=L_{m+2}+1,                \tag{6}
\]

the cutoff-tail exponents are

\[
 \begin{array}{c|c}
 e_t,-\Delta e,\nabla\pi&\tau\\
 (U^{[J]}\cdot\nabla)e,(e\cdot\nabla)U^{[J]}&\tau-H_{m+2}\\
 (e\cdot\nabla)e&2\tau.
 \end{array}                                             \tag{7}
\]

The frozen threshold and the one-stage enlargement give, for every fixed
`m,N`, a finite `J` for which the log-corrected fixed-state exponent and all
six values in (7) are at least `N+1`.  Consequently

\[
             |\mathcal R(u_*,p_*)|_m=O(q^{N+1})          \tag{8}
\]

on the accepted later-tail chart, conditionally on the Q1--Q5 contracts.
Since `N` is arbitrary, the local residual is flat in every fixed jet.

## 4. Correction: disjoint margins cannot be added

The common exponent in (8) is the minimum of the fixed-state exponent and
the six cutoff-comparison exponents.  The regular `h delta/2`, the two compact
`h kappa/2` margins and the cutoff reserve have different owners and cannot
be summed into a stronger common margin.  For the default `m=N=0` ledger,
the comparison stage is `J=12817`; the limiting terms are the two nonlinear
cross orientations, with exponent

\[
 {20030773\over20000000}=1.00153865,                    \tag{9}
\]

while the log-corrected fixed-state exponent is much larger.  Equation (9)
is exact exponent arithmetic, not a numerical evaluation of 12,817 fields.

## 5. Scope and next obligation

This closes the complete six-group **local** residual replay on the exact
Q1--Q5 object, relative to the stated continuum contracts.  It does not prove
a uniform all-stage pulse bound, terminal velocity or pressure regularity,
an admissible compact force, global comparison, exclusion or blow-up.

The next justified obligation is to transfer this Q6 estimate to the closed
Cartesian terminal corner and join the punctured one-sided endpoint jets on
the same field.  That replay must retain the enlarged `A0/A1` requests and
must not infer terminal field regularity at the origin from residual
flatness.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q6_complete_residual_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q6_complete_residual_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q6_complete_residual_join/evidence.json \
  --output fresh.json
```
