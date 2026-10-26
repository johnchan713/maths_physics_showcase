# Certified-domain complete-potential diagonal

Parent: `9be61b33c2651c198f89604a76fd737e334d2c72`.
Status: `certified-domain-conditional-complete-potential-diagonal-join`.

The parent joins every fixed finite-stage compact graph to certified
primitive domains, but records the absence of a stage-growth estimate as an
obstruction to a `q`-dependent diagonal.  That statement is correct for the
naive substitution `j=j(q)` into an uncut fixed-stage estimate.  It is not an
obstruction to the already derived complete-potential cutoff diagonal.

This checkpoint corrects that scope and joins the certified fixed-stage
graphs to the locally finite construction.  The conclusion remains a
conditional analytical implication.  No infinite coefficient sequence,
large common-pole constant, or cutoff scale is numerically evaluated, and
the closed-time global forcing and comparison chain is not claimed here.

## 1. Certified finite states and complete representatives

Keep

\[
 \kappa=10^{-5},\qquad \delta=\tfrac12-3\kappa,
 \qquad B_j=\tfrac32-2\kappa+j\delta.                  \tag{1}
\]

The parent now supplies finite acyclic high-jet graphs on fully expanded
primitive domains for every fixed stage, jet, and requested power.  The
corrected reference values are

| Request | Primitive ceiling |
|---|---:|
| normalized stage 0 | 900,047 |
| normalized stage 1 | 1,400,076 |
| physical stage 0, `m=N=0`, `h=1/200` | 60,900,051 |

The illustrative stage-two physical pole archive uses `q=371400672`.
These are symbolic sufficient domains.  They are not allocated arrays or
sampled high derivatives.

For the actual cycle from stage `j-1` to `j`, retain the fixed operator
gauges and write

\[
 \Delta u_j=\operatorname{curl}A_j+b_j e_\theta,
 \qquad \partial_\theta b_j=0,
 \qquad \pi_j=p^{[j]}-p^{[j-1]}.                       \tag{2}
\]

`A_j` contains the complete sourced and signed wave potentials and the fast
and slow mean potentials.  Direct swirl remains outside the curl.  The
pressure difference telescopes every intermediate reconstruction, retaining
all six covariance entries, axial `A1` velocity, moving pressure moment, and
both same-base pressure tails.  Arbitrary stage-dependent potential or
pressure gauges are not allowed: cutting them would introduce
`grad chi x grad phi_j` or `c_j(t) grad chi`.

The inherited representative estimate is

\[
 |(A_j,b_j e_\theta,\pi_j)|_m
 \le C_{j,m}(1+|\log q|)^{P_{j,m}}
      q^{h\delta j-\ell_m},\qquad \ell_m=2A+2m,        \tag{3}
\]

with arbitrary finite `C_{j,m}` and `P_{j,m}` at each fixed stage.

## 2. Complete cut and scale selection

For one fixed smooth cutoff `chi`, define

\[
 V_{j,a}=\operatorname{curl}(\chi(aq)A_j)
          +\chi(aq)b_j e_\theta,
 \qquad P_{j,a}=\chi(aq)\pi_j.                          \tag{4}
\]

This cuts the complete potential, direct swirl, and pressure.  The velocity
is exactly divergence-free.  All cutoff-curl terms, including the radial
mean commutator, and all pressure cutoff gradients remain.

After reserving one curl, the exact implicit-coordinate derivative bounds
give

\[
 |(V_{j,a},P_{j,a})|_m
 \le \widetilde C_{j,m}(1+|\log q|)^{\widetilde P_{j,m}}
       q^{h\delta j-L_m},                               \tag{5}
\]

uniformly in `a` on the cutoff transition.

At stage `j`, collect the finitely many constants in (5) for jets through
`j`.  Put `gamma_j=h delta j/2`.  Once those constants are fixed, choose
`s_j=log a_j` strictly above

\[
 \max\left\{s_{j-1}+1,\log(q_*^{-1})+1,1,
 \frac{2(c+j)}{\gamma_j}+\frac{P^2}{\gamma_j^2}\right\} \tag{6}
\]

for every current pair `c>=log C`, `P>=0`.  The elementary inequality
`log(1+s)<=sqrt(s)` then yields, on the whole interval `0<q<=1/a_j`,

\[
 C(1+|\log q|)^Pq^{\gamma_j}\le2^{-j}.                 \tag{7}
\]

This choice exists for every finite list, however fast the constants grow.
It uses no stage-uniform bound and no convergence of the raw uncut series.

## 3. One locally finite field

Choose the sequence once and set

\[
 u^*=u^{[0]}+\sum_{j\ge1}V_{j,a_j},
 \qquad p^*=p^{[0]}+\sum_{j\ge1}P_{j,a_j}.              \tag{8}
\]

Since `a_j` grows strictly, only finitely many summands are active on every
compact subset of `q>0`.  Thus (8) is one smooth real divergence-free local
field.  The closed annular supports preserve the protected inner velocity
and heat exterior, conditional on the inherited support contracts.

For a fixed comparison stage `J` and `q<=1/(2a_J)`, the first `J` cutoffs
equal one and telescope in their fixed gauges.  The remaining locally finite
tail satisfies

\[
 |(u^*,p^*)-U^{[J]}|_m
 \le2^{-J}q^{h\delta(J+1)/2-L_m}.                       \tag{9}
\]

The comparison stage is chosen after fixing `m,N` and is held fixed as
`q` tends to zero.  There is no substitution of an uncontrolled `j(q)`
into the uncut estimate.

## 4. Correction of the previous obstruction

The parent recorded the sufficient inequality

\[
 \log C_{j(q),m}\le(a_{j(q)}-N)\log(1/q).              \tag{10}
\]

Equation (10) is indeed required if one tries to select an uncut stage
directly as a function of `q`.  Construction (8) does not do that.  It
absorbs each finite constant when choosing `a_j`, then compares to one fixed
stage.  Therefore the absence of a stage-growth estimate is not an
obstruction to the locally finite complete-potential diagonal.

For `(e,pi)=(u^*,p^*)-U^{[J]}`, the exact residual difference is

\[
 e_t-\Delta e+\nabla\pi+(U^{[J]}\!\cdot\nabla)e
 +(e\cdot\nabla)U^{[J]}+(e\cdot\nabla)e.               \tag{11}
\]

Time, viscosity, pressure, both nonlinear crosses, and the new square are
all present.  Choose one finite `J` satisfying the two positive-slope
exponent inequalities of the inherited diagonal theorem.  For the default
`h=1/200`, `m=N=0` budget, the conservative sufficient value is `J=12816`.
This is exact existence arithmetic, not a construction or numerical
evaluation of 12816 stages.  It yields conditionally

\[
 \forall m,N,\qquad |\mathcal R(u^*,p^*)|_m=O(q^N).    \tag{12}
\]

The next obligation is to join this same field, cutoff sequence, support,
and pressure gauge to the closed-time endpoint atlas and complete one-sided
jets.  Localization, zero initialization, compact smooth forcing,
whole-space energy, global comparison, and the core lower-bound argument
must then be checked on that same object.  No unconditional residual theorem,
admissible force, or blow-up claim is made.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_certified_complete_potential_diagonal_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_certified_complete_potential_diagonal_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_certified_complete_potential_diagonal_join/evidence.json \
  --output fresh.json
```
