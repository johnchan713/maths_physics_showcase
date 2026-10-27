# Closed-chart coefficient finiteness from the fixed primitive domains

Parent: `8b627a88a3d55259e9c35827192e58109b628968`.
Status: `closed-chart-coefficient-finiteness-joined-relative-to-certified-primitive-contracts`.

The parent corrects the cutoff-tail quantifier and leaves two coefficient
constants to be proved on the closed angular chart and compact `X` intervals.
This checkpoint traces those constants through the actual fixed-stage provider
graph.  For every fixed stage and derivative request, they follow by finite
composition from the already named primitive closed-domain contracts.  There
is no additional pole or compact-`X` concentration created by the finite
graph.

This is a conditional join.  It does not independently reprove the exact
leading profile, Borel background, global support geometry, or the incoming
estimates used by the root and compact maps.  Those remaining primitive
premises are stated precisely below.

## 1. Every primitive domain already contains the terminal poles

Use

\[
 T+\eta^2=1,\qquad |\eta|\le1,\qquad
 L=1-2h\eta^2\ge49/50.                                  \tag{1}
\]

The provider chain supplies the following domains:

| Primitive | Closed domain and nondegeneracy |
|---|---|
| implicit physical chart | `1/4<=q/Q<=4`, `T>=0`, `|eta|<=1`, `L>=49/50` |
| summed background | fixed compact radial enlargement with leading shear, growth and cone gaps |
| annular matching root | factorial `C^m([-1,1])` algebra; pointwise Jacobian inverse at most `1220` |
| compact repair root | the same closed angular interval; inverse at most `1600`, with `0<lambda<=1/100` |
| pulse ODE | closed pulse interval, positive spectral gap `mu`, positive radial amplitude and fixed slow boxes |
| compact and slow maps | fixed smooth projectors and closed reserved support hulls |
| support geometry | fixed compact `X` annuli and separated closed label supports |

All endpoint derivatives are one-sided.  None of these formulas divides by
`T` or `1-eta^2`.  Differentiating reciprocal factors uses positive floors;
for example, derivatives of `1/L` obey a finite recurrence at every fixed
order because `L>=49/50`.

## 2. Finite-composition lemma for one fixed stage

Fix a stage `J`, jet `m`, and requested power.  The certified recurrence
expands the state into a finite acyclic graph with a finite primitive ceiling,
finite harmonic set, finitely many archived poles, and the same closed
supports.  Its nodes use only finite sums and products, derivatives,
compositions, fixed-profile integrals, ODE solution operators, implicit roots,
and inverses on the positive floors in the table.

Induct on a topological ordering of that graph.  Sums and products preserve
finite `C^m` seminorms.  The chain and Leibniz rules use finitely many lower
seminorms.  Reciprocal differentiation follows from

\[
 D^n(f^{-1})=-f^{-1}\sum_{k=1}^n {n\choose k}D^kf\,
                                      D^{n-k}(f^{-1}),    \tag{2}
\]

with `|f|` bounded below.  The root identities use their bounded inverse
Jacobians and lower root jets.  Duhamel differentiation uses the finite closed
pulse interval and bounded propagator.  Hence every node has a finite
closed-domain seminorm for this fixed request.  Constants may be enormous and
may depend arbitrarily on `J,m`, the fixed profiles, and the requested power.

This is an analytical compactness argument.  The positive-polynomial controls
in `compactness.py` verify only the finite recurrence and request graph; their
numerical values are not actual field constants.

## 3. Why compact-`X` uniformity follows

Every normalized correction profile is supported in the same fixed compact
`X` annuli.  Physical scaling changes `q` and the band coordinates, while the
profile variable remains `X=r^2/(2q)` on those fixed supports.  The finite
provider graph neither introduces an `X/q` argument nor moves the normalized
support toward `X=0` as `q` tends to zero.  Products of distinct labels vanish
by the closed separation contract; same-label operations preserve the fixed
support.

Consequently the earlier model

\[
 (X/q)e^{-(X/q)^2}
\]

is a valid warning against pointwise inference but is excluded by the actual
fixed-profile support contract.  On any requested compact `K`, the finite
graph is uniformly bounded on `K` intersected with its fixed supports.  Outside
the correction annuli, the unchanged summed background and exact heat exterior
apply.

## 4. The two coefficient conclusions

The finite-stage complete residual is a polynomial in the state, its fixed
derivatives, viscosity one, and the complete fixed-gauge pressure.  The graph
therefore gives, relative to the primitive contracts,

\[
 C^{\rm stage}_{J,m,K}:=
 \sup_{\substack{0<q\le q_0,\ X\in K,\ |\eta|\le1\\T=1-\eta^2}}
 q^{-(hB_J-K_m)}|\mathcal R(U^{[J]})|_m<\infty.         \tag{3}
\]

The increment graph similarly gives a finite logarithmic power and

\[
 C^{\rm rep}_{j,m,K}:=
 \sup_{\substack{0<q\le q_0,\ X\in K,\ |\eta|\le1\\T=1-\eta^2}}
 {q^{-(g_j-L_m)}|Z_j|_m\over(1+|\log q|)^{P_{j,m}}}<\infty. \tag{4}
\]

Complete covariance, axial velocity, moving pressure moments, both exceptional
pressure memories, compact `A0/A1` pieces, viscosity, cutoff curls, nonlinear
crosses and squares remain in these graphs.  Equation (4), the parent's
half-gain absorption, and the previous Cartesian transfer close the corner
estimate relative to the same primitive contracts.

## 5. Exact remaining premise boundary

Graph finiteness cannot prove that its starting domain is nondegenerate.  The
independent foundation must establish on the unchanged common tail

\[
 \Delta_{\rm prim}:=
 \min\{\inf F,\inf|g|,\inf(-N_\theta),\inf\lambda^2,
          \inf(1-r_{\rm cone}),\Delta_{\rm support},
          \Delta_{\rm root\ input},\Delta_{\rm compact\ input}\}>0, \tag{5}
\]

together with finite all-fixed-order closed-chart seminorms of the exact
leading/background and input profiles.  The root inverse bounds then follow
from their proved contraction estimates, and every later constant follows by
the finite-composition induction above.

The cited primitive checkpoints assert (5) relative to their exact-core,
shape-entry, compact-input and support hypotheses.  This checkpoint does not
replace an independent proof of those inputs with graph arithmetic or finite
diagnostics.  Therefore compact forcing, the complete residual conclusion,
global comparison and blow-up remain conditional.

The next obligation is to audit (5) and the associated seminorms directly from
the defining leading profile, Borel background, support construction and root
inputs.  No numerical evaluation of giant scales is required for existence,
but every positive floor and support identification must be derived on the
actual unchanged field.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_closed_chart_coefficient_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_closed_chart_coefficient_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_closed_chart_coefficient_join/evidence.json \
  --output core_closed_chart_coefficient_join_replay.json
```
