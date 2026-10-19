# Sharp joint cutoff-shape feasibility and an eliminated scale-budget test

Status: `sharp-smooth-cutoff-capacity-and-eliminated-affine-budget-test`.
Parent: `355dce1e76bf0e6d97ec7da84aa06e0fddff5522`.

This checkpoint derives a new continuum inequality for the remaining
coefficient-cutoff input. Its two derivative upper bounds cannot be chosen
independently: a unit drop must fit inside the retained transition window.
We characterize exactly which pairs can bound SOME smooth shape, derive
the sharp unattained infimum of the viscous stress-row cost, and eliminate
the two shape variables from the existing four-row sufficient scale bound.

This is an **existence/obstruction test for that bound**, not identification
of the UNCHANGED cutoff or a proof of actual cone failure. No cutoff,
coefficient, leading field, root, frequency, tilt, scale, pressure gauge or
mean correction is replaced. Historical frozen evidence stays unchanged.
In particular the preceding rectangular corner (1,4) was explicitly only
necessary; it is not an admissible smooth shape. We sharpen that relaxation,
not silently reinterpret its old success flags as actual shape existence.

## 1. Exact joint feasibility of normalized derivative enclosures

Use the original keep/cut window in logarithmic X, with corrected

\[
 \Delta=t_c/16=\frac{C^{-20}}{160\sqrt{\log C}},\qquad
 K(s)=\kappa(X_{keep}e^{\Delta s}).
\]

The admissible class consists of smooth real K, equal to 1 for s<=0 and
0 for s>=1, with 0<=K<=1. Necessity below does **not** assume monotonicity.
Let M,N be positive finite UPPER ENCLOSURES of sup|K'|,sup|K''|. Set

\[
 W_{M,N}(s)=\min(M,Ns,N(1-s)),\quad 0\le s\le1,
\]
\[
 \mathcal A(M,N)=\int_0^1W_{M,N}(s)\,ds
 =\begin{cases}
 M-M^2/N,&2M\le N,\\
 N/4,&2M\ge N.
 \end{cases}                                                \tag{1}
\]

**Claim.** There exists an admissible smooth K with these two upper
enclosures if and only if

\[
 \boxed{\mathcal A(M,N)>1.}                                 \tag{2}
\]

Necessity: write v=-K'. The endpoint flatness gives v(0)=v(1)=0,
|v'|<=N, and hence v(s)<=W(s), while `int_0^1 v=1`. Smooth extension by
constants also gives v'(0)=0, so v(s)=o(s). For sufficiently small positive
s, W(s)=Ns and v(s)<Ns/2. The integrated comparison is therefore STRICT:
`1=int v<int W`. Negative v or oscillations cannot evade this upper bound.
The clipped-linear integration yields (1). This is an analytical argument
on the whole interval, not inference from finitely sampled derivatives.

Sufficiency, without using a nonsmooth tent as K: extend W by zero, put
`A=mathcal A>1`, `eps=(A-1)/(16A)<1/16`, and shrink its support:

\[
 V_\epsilon(s)=(1-2\epsilon)
 W\bigl((s-\epsilon)/(1-2\epsilon)\bigr).
\]

Its support is [eps,1-eps], sup<= (1-2eps)M, Lipschitz constant<=N,
and mass `m=(1-2eps)^2 A>1`. Convolve by ANY nonnegative unit-mass smooth
kernel supported in [-eps/2,eps/2]. Convolution preserves mass, nonnegativity
and those upper bounds. Divide by m to obtain a smooth v of mass one,
supported in [eps/2,1-eps/2], with

\[
 \|v\|_\infty\le(1-2\epsilon)M/m<M,\qquad
 \|v'\|_\infty\le N/m<N.                                  \tag{3}
\]

Then `K(s)=1-int_0^s v` is smooth, monotone and admissible, with the
required enclosures. This auxiliary existence proof does NOT select it as
the retained kappa. It is used to establish sharpness of the budget test.
All padding is INSIDE the same allowed window, not a change of that window.

Equivalently (2) requires M>1 and

\[
 N>\begin{cases}M^2/(M-1),&1<M\le2,\\4,&M\ge2.\end{cases}   \tag{4}
\]

Thus neither M=1 nor N=4 works at any finite companion bound. The old
corner has `mathcal A(1,4)=3/4`, not one. At M=2,N=4 the triangular area
is exactly one, but it still cannot be attained by a smooth flat-ended K.

Do not confuse upper enclosures with exact attained norms: a loose first
bound M=1000 and N=1216 is feasible by (2). Imposing N>=M^2 on arbitrary
upper bounds would wrongly reject it. The previously proved squared-flat
step already admits such loose enclosures; it is not newly assigned to the
retained field. Formula (2), not a relation between arbitrarily inflated
declared upper bounds, is the invariant statement.

## 2. Sharp, unattained lower cost in the complete stress row

The preceding complete bound has

\[
 R_w=A_0+A_1M+A_2N,\qquad A_1\ge0,\quad
 A_2=\frac{\nu c^3(c+2)J}{\tau\Delta^2 X_{keep}^2}>0.       \tag{5}
\]

All coefficients are those of the complete pressure/viscosity interface;
actual nu=1. For every feasible smooth-shape enclosure pair,

\[
 \boxed{R_w>A_0+A_1+2A_2+2\sqrt{A_2(A_1+A_2)}=:I_w.}       \tag{6}
\]

For 1<M<=2 write t=M-1 in (0,1]. By (4),

\[
 A_1M+A_2N>A_1+2A_2+(A_1+A_2)t+A_2/t
 \ge A_1+2A_2+2\sqrt{A_2(A_1+A_2)}.
\]

For M>=2, N>4 gives the result from
`2A1+4A2>=A1+2A2+2sqrt(A2(A1+A2))`, since the difference-of-squares
comparison reduces to A1^2>=0. The infimum in (6) is sharp: take
`t*=sqrt(A2/(A1+A2))`, `M*=1+t*`, and N approaching
`2+t*+1/t*` FROM ABOVE; (2)-(3) construct smooth shapes within these
budgets whose cost approaches the infimum. It is never attained.

For A1>0 this strictly improves the stress cost at the relaxed corner
(1,4). It is an infimum for Rw, NOT a claim that the full max of all four
rows has the same infimum or cannot attain a larger constant row.
For rational inputs, the necessary test B>Iw needs no floating square root:
with `beta=B-A0-A1-2A2`, require `beta>0` and
`beta^2>4A2(A1+A2)`. Passing it alone ignores the derivative and value rows.

## 3. Exact elimination of shape budgets in the four-row certificate

Retain the entire preceding max-affine expression

\[
 C_\star(M,N)=\max(1,R_v,D_0+D_1M,A_0+A_1M+A_2N),\quad D_1>0.
                                                               \tag{7}
\]

At an existing tilt/scale set `B=r_Kc*2^(2h(ell-2))` only AFTER its
same-tilt margin and other inputs are supplied. The following decision is
iff existence of SOME feasible smooth-shape enclosures satisfying (7)<=B,
with all displayed row coefficients held fixed. It is **not** an iff test
for the actual PDE cone, old cutoff or complete residual.

First require `1<=B`, `Rv<=B`. Monotonicity of the rows lets an existential
search restrict to 1<M<=2: if a feasible pair has M>2, then N>4 and the pair
(2,N) is also feasible and lowers the row bounds. This can produce a
DIFFERENT auxiliary shape; it is not a bound on the same old K.

Put

\[
 t_{cap}=\min(1,(B-D_0)/D_1-1),\quad
 P=A_1+A_2,\ Q=A_2,\ \beta=B-A_0-A_1-2A_2.                 \tag{8}
\]

If tcap<=0 or beta<=0, no feasible pair exists for this certificate.
Otherwise set

\[
 t_0=\min(t_{cap},\beta/(2P)),\qquad
 q_{min}=P t_0^2-\beta t_0+Q.                             \tag{9}
\]

The exact remaining test is

\[
 \boxed{q_{min}<0.}                                      \tag{10}
\]

Indeed t=M-1 must lie in (0,tcap]. Equations (4) and the stress row allow
some N iff `P t+A2/t<beta`, equivalently `P t^2-beta t+Q<0`.
The positive-vertex quadratic attains its minimum over this interval at
(9). Equality qmin=0 permits no smooth flat-ended shape; changing < to <=
would promote the nonsmooth limiting profile incorrectly.

When (10) holds there is a wholly rational enclosure witness:

\[
 M=1+t_0,\quad N_{min}=M^2/t_0,\quad
 N_{max}=(B-A_0-A_1M)/A_2,\quad N=(N_{min}+N_{max})/2.       \tag{11}
\]

It has `Nmin<N<Nmax`, capacity>1, and every original row<=B.
The construction (3) supplies an auxiliary smooth shape. No square root,
parameter scan or rounded optimization is needed in the decision/witness.

For example the exact auxiliary rows `Rd=M,Rw=3M+N`, constants one, pass
the old rectangular corner at B=8, yet (6) has infimum 9: NO smooth shape
can pass this bound at B=8 or 9. At B=10, (11) gives feasible budgets.
These are controls, not the actual enormous-profile constants.

## 4. Composition, historical scope and remaining obstruction

The implementation imports, without modifying, all four rows from the
preceding full corrected-width gate. A separate control checks its constant
stress row decreases if full w0/CF/CU inputs are erased; that erasure is a
negative control, not part of the result. Signed full-future pressure
moments, the separate pointwise source, five additive repairs, original
viscosity, tau only on [Xinner,Xc], both half-budgets and the +2 chart cost
remain unchanged. In particular the wrong tc=t1 identification is NOT
reintroduced. No reserved support or pressure memory is reset.

If (10) fails with valid full row bounds and an accepted r_Kc, it excludes
ALL smooth cutoffs from satisfying THIS sufficient upper-bound certificate
at that B. It does not prove actual stress/cone failure; upper bounds may
be too conservative. If it passes, only auxiliary shape-budget existence
has been established. The old kappa still needs identification and its own
upper enclosures. Applying a newly constructed shape would change the
coefficient and later realization and require a separately traced
construction, not merely relabelling this as acceptance of the old field.
The retained tilt margin and scale remain unverified too; implicit finite
constants and manufactured inputs are not their numerical certification.

This continuum cutoff theorem is independent of the Navier-Stokes regional
premises. Its application to the complete rows is conditional on those
premises. The source is the retained annular cutoff specification and
complete first-stress interface, with new elementary analysis (1)-(11);
no external claimed blow-up theorem is accepted as established.

Checks: exact integration on 25 rational capacity pairs; four strict
frontier/padding controls; 15 sharp weighted-cost approximation controls;
symbolic clipped-area, square-completion and capped-minimum identities;
boundary and full-interface compositions. Finite checks are algebra
diagnostics, not proofs of the continuum premises. Old source bindings and
evidence replays use original tolerances and remain frozen.

Next identify or bound the unchanged cutoff's actual derivative norms and
verify its same-tilt/retained-scale inputs; use (2) and (10) to reject
impossible proposed budgets before attempting actual acceptance. Uniform
physical pulse jets/supports/inverse gaps, complete wave-covariance means
including axial velocity/pressure, and the full viscous/cutoff/interaction
residual remain open. No unconditional pulse bound, residual improvement
or blow-up is proved.

```sh
python3 research/navier_stokes_cascade/results/cutoff_shape_feasibility/test_audit.py -v
python3 research/navier_stokes_cascade/results/cutoff_shape_feasibility/audit.py \
  --verify-record research/navier_stokes_cascade/results/cutoff_shape_feasibility/evidence.json \
  --output /tmp/fresh_cutoff_shape_feasibility.json
```
