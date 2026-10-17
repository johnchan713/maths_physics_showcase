# Actual annular and compact matching-root foundation

Parent: `655f9510308d39d001415d176bf259e085baf638`.
Status: `actual-annular-and-compact-matching-roots-audited-relative-to-retained-core-and-compact-input-estimates`.

This checkpoint supplies quantitative inverse, target and implicit-root estimates
for the **existing** five-row annular and compact-repair maps. It identifies the
same roots selected by the previous construction. It is relative to the retained
exact-core/shape-entry and actual compact-input estimates, whose own analytical
dependencies remain. It does not assert a complete unconditional construction,
residual improvement, Navier–Stokes solution or blow-up.

The reference axial and scalar-strength roots were handled in
[the parent checkpoint](../axial_root_foundation/README.md). Historical evidence
is unchanged. The new evidence binds 18 direct sources and replays 1,069
inherited source bindings. Rational inequalities below are continuum estimates;
finite quadrature and polynomial controls are diagnostics, not a continuum proof
or formal certificate. The supplied manuscript and its SHA-256 remain those
recorded by the parent. No parameter, frequency, support, pressure gauge, label,
cutoff sequence or original width `C^-20` is changed.

## 1. Norm and root lemma

For a scalar function of eta use

\[
 |f|_k=\sum_{j=0}^k\sup_{[-1,1]}|\partial_\eta^jf|/j!,\qquad
 |c|_k=\max_i|c_i|_k.
\]

These are product algebra norms. Endpoints use one-sided derivatives. For
`Bc+Q_eta(c,c)=z`, if `||B^-1||<=beta`, the bilinear map has algebra norm at
most q, and `|z|_k<=delta`, the ball `r=2 beta delta` is invariant and the map
is a contraction when `8 beta^2 q delta<1`. Its Lipschitz bound is
`4 beta^2 q delta<1/2`. The pointwise root Jacobian has inverse norm at most
`2 beta` by the Neumann series. This argument uses exact integrals throughout.

## 2. Actual annular matrix and eta-dependent quadratic term

Use the existing `MomentMap` in
[axis_matching_audit/moments.py](../axis_matching_audit/moments.py).
Set `x=X/X_R`, `y=log x`, `x_c=exp(-6)`, `t=y+6`, and `f=(1+eta^2)^-1`.
The ideal fields are `E=P f x^.1`, `U=4eta`, and `H=sqrt(2X) E`.
The five normalized rows, with hatted moments scaled by the corresponding
powers of `X_R`, are exactly

\[
 \frac{\Delta\widehat M}{x_c},\quad
 \frac{\Delta\widehat J-4\eta\Delta\widehat I}{\sqrt2 P f x_c^{1.6}},\quad
 \frac{\Delta\widehat I}{\sqrt2 P f x_c^{1.6}},\quad
 -\frac{\Delta\widehat S-8\eta\Delta\widehat M}{P^2f^2x_c^{1.2}},\quad
 \frac{\Delta C_p}{P^2f^2x_c^{.2}}.
\]

The two additive U bumps have centers `.3,.7`; the three relative E bumps
have centers `.15,.5,.85`, all width `.08`. A bump here is
`beta_j= sigma'((t-center)/.08+1/2)`, of exact mass `.08`, supremum at most
64, and `beta_j^2<=64 beta_j`. The five supports are pairwise disjoint.
The two constant linear blocks have rates `(1,1.6)` and `(1.6,1.2,.2)`.
For equally spaced centers their entries are `B_kj=m_k n_k^j`,
`n_k=exp(spacing*rate_k)`, and
`m_k>=.08 exp((first_center-.04)*rate_k)`.

The coefficients of the Lagrange polynomial
`product_(l!=k)(x-n_l)/(m_k product_(l!=k)(n_k-n_l))` give the inverse
columns. Taking absolute coefficients, enclosing exponentials with positive
Taylor series and a geometric tail, gives row upper bounds

\[
 (75.570696,44.195911),\qquad
 (417.298607,607.488714,214.267257).
\]

Thus **beta=610** is valid without quadrature. `matching_bounds.py` recomputes
outward rational bounds rather than importing these decimal conclusions.

Disjointness deletes the mixed U/E term in J; it does not delete either square.
The only nonlinear rows are

\[
 Q_3(c,d)=\tfrac12\sum_E\!\int e^{1.2t}\beta_j^2 c_jd_j
       -\zeta(\eta)\sum_U\!\int e^t\beta_j^2 c_jd_j,
 \quad Q_4(c,d)=\tfrac12\sum_E\!\int e^{.2t}\beta_j^2 c_jd_j,
\]
\[
 \zeta=e^{1.2}(1+\eta^2)^2/P^2,\qquad P\ge16.
\]

Let `E0=32*.24*exp(1.6)`, `U0=64*.16*exp(1)`,
`gamma=exp(1.2)/256`. The factorial normalized coefficient-jet bounds are
`q_0<=E0+4 U0 gamma` and `q_i<=U0 gamma*(8,8,4,1)_i` for i=1..4.
Their sum is below `47.065<50`. Hence **q=50** works in every finite factorial
eta algebra; this is not an assertion that high-order targets are small.

## 3. Actual incoming annular target and same root

The retained [shape-entry estimate](../shape_join_pressure/README.md) gives
factorial C2 normalized incoming rows below `epsilon/4` at the separation
point and `|g|_2<=g_b=3epsilon/100`, `epsilon=10^-16`.
The core logarithmic choice `log C=10000B+100 log(1/epsilon)`, `B>=2^16`,
makes its displayed ledger bound `2*10^6 exp(412B) C^-2<epsilon/4`:
indeed `8*10^6 epsilon^199 exp(-19588B)<1`.
These are the existing current-source input inequalities, not manufactured
incoming axis values from a numerical experiment.

Between separation and y=-8 the five additional bounds are

\[
 (e^{-2}g_b,\tfrac58e^{-16/5}g_b,0,
       \tfrac{20}{256}e^{-4/5}g_b^2,0).
\]

Restoring U from `4eta+g` to `4eta` on -8<y<-7 adds at most `g_b` to
each of the first two rows and `g_b^2` to the energy row. The exact identity
`U^2-(4eta)^2-8eta(U-4eta)=g^2` keeps the cancellation and all angular
derivatives. The normalized zeta C2 bound is below one. Normalization
derivatives are already included (`|f^-1|_2=5`, `|f^-2|_2=20`). Consequently

\[
 \delta<.286\epsilon,\quad
 r=2(610)(.286\epsilon)=348.92\epsilon<349\epsilon,
 \quad 8(610)^2(50)(.286\epsilon)=4.256824\,10^{-9}<1.
\]

The new C2 root is inside both the previous C2 ball `572epsilon` and the
previous C1 uniqueness ball `6000epsilon`; it therefore identifies the same
accepted root, rather than selecting a different branch. Its pointwise inverse
is at most 1220. The existing annular exterior and cone estimates, with their
larger accepted coefficient envelopes, continue to apply relative to their
retained inputs. This does not discharge an independent global pressure lemma.

## 4. Actual compact-repair matrix and target

Use [the existing Repair](../stress_realization_audit/repair.py) on the first
reserved patch, with `y=log(X/rho_p)` in (0,5), `U0=0` and
`E0=e(eta) exp((-1/2-lambda)y)`. Both corrections use the **same original** e:
`deltaU=e sum u_j beta_j`, `deltaE=e sum v_j beta_j`.
The U centers are 1,2; E centers .5,1.5,2.5; width .1. Here
`beta_j=sigma'((y-center)/.1+1/2)/.1` has unit mass, bound 640 and
`beta_j^2<=640 beta_j`. These are not the annular bump normalizations.

Normalize M,I,J,S,Cp respectively by
`rho_p e`, `sqrt(2)rho_p^1.5 e`, `sqrt(2)rho_p^1.5 e^2`,
`rho_p e^2`, `e^2`. In row order `(m,(m-j)/lambda,i,-s,c_p)`
the U weights are `exp(y)` and `exp(y)g_lambda(y)`, where

\[
 g_\lambda(y)=(1-e^{-\lambda y})/\lambda,\qquad
 g'_\lambda(y)=e^{-\lambda y}>0,\quad 0<\lambda\le.01.
\]

The two supports have separation at least .9. Their determinant is at least
`.9 exp(.95+1.95-.0205)`. Bounding the inverse numerator by
`3.05 exp(2.05)` gives U inverse below 2.
For the E block the rates are `(1.5,.5-lambda,-.5-lambda)`.
For each triple of integration points put `x_j=exp(y_j)` and factor
`x_j^(-.5-lambda)` from columns. The remaining rows are
`(x_j^(2+lambda),x_j,1)`. Their divided difference is at least one since
`(x^(2+lambda))''>=2` for x>=1. Multilinearity and positive unit masses give
the determinant lower bound

\[
 d_E=e^{-.51(4.65)}(e^{1.45}-e^{.55})
             (e^{2.45}-e^{.55})(e^{2.45}-e^{1.55}).
\]

Each cofactor is at most `2 exp(7.65)`, hence
`||B_E^-1||<=6 exp(7.65)/d_E<787.562<800`.
Both square contributions remain in energy; the E square remains in pressure.
Their bilinear algebra bound is
`640[2 exp(2.05)+1.5 exp(2.55)]<22237.734<25000`.
Thus **beta=800,q=25000** are valid. B and Q are independent of eta only
after the complete original-e normalization above.

Use the full **sum** C1 norm from section 1. The older repair used
`max(sup|f|,sup|f_eta|)`, whose product bound costs a factor two. The new sum
norm controls that older maximum norm; it must not be silently identified with
it. The actual loop discrepancy bound, with all normalization derivatives and
the divided-difference loss, is

\[
 \delta\le2A m/\lambda=64A^5H/(\lambda N)
 \le64A^6H/N\le D/N,\quad D=H^{16}.
\]

Here the retained compact envelope supplies normalizers and their needed
derivatives at most A, `lambda^-1<=A`, and `m=32A^4H/N` in the full sum norm.
The original frequency is `N=1+floor(H^32)>H^32`, with `A=C^4096>=2^16`,
`G=exp(A^32)`, `H=exp(G^4)`. In particular H>A and
`8*800^2*25000 delta<2^37 H^-16<2^-219<1`.
The same root has full sum C1 norm at most **1600D/N**, pointwise inverse
at most 1600, and lies in the older `2000D/N` uniqueness ball.
The older acceptance inequalities remain valid for this same pointwise root;
no new frequency choice or extrapolation from small diagnostic frequencies is
used. The original denominator floors, gaps and repair-to-state estimates are
retained dependencies. In particular `1600H^-16<A^-128` follows immediately
by comparing logarithms with `log H=exp(4A^32)`.

Stabilizing B does not remove the `1/lambda` input loss. If ordinary M and J
are already rounded equal, their lost difference cannot be recovered by matrix
preconditioning; the analytical or stable integral ledger must retain it.
The continuous extension `g_0=y` is useful for bounding weights but does not
make the original independent M and J equations invertible at lambda=0.

## 5. Fixed jets and complete pressure memory

Differentiate the exact same root identity, using normalized jets
`d_n=sup|partial_eta^n c|/n!`, `z_n=sup|partial_eta^n z|/n!`:

\[
 d_n\le2\beta\left[z_n+
 \sum_{i+j+l=n\atop(i,j,l)\ne(0,n,0),(0,0,n)}q_i d_jd_l\right].
\]

For the annular map i ranges from zero to four. Omitting the i>0 terms would
incorrectly treat zeta as constant; the exact polynomial control detects this
error. For the repair only q0 occurs. This gives every fixed finite one-sided
jet from the pointwise inverse and preceding jets, without a new high-order
small-target contraction. The base root identity differentiated inductively
gives the actual jets, not a separate root chosen for each order.

At either patch end, the five exact moment differences vanish as functions of
eta, and hence so do all fixed eta jets. With the same subsequent fields this
preserves M,J,I,S and **the complete** `Pi_axis+Cp` at every later radius.
No pressure memory is reset at an interface. A four-row repair would instead
leave an arbitrary c(eta) in pressure; its Ns defect is
`4 Aphys eta c-(1-eta^2)c_eta`, which need not vanish even where c itself is
zero. Thus local root existence alone is not the complete global pressure
identification theorem.

The next obligation is an independent same-gauge trace of the complete core
backward-pressure identity through these roots and every later attachment.
Then the remaining finite-background and actual-pulse existence, inverse,
gap, size, derivative and support lemmas, complete-wave covariance mean
corrections (including axial velocity and pressure), and the full viscosity,
cutoff and interaction residual remain to be discharged before unconditional
theorem assessment.

## 6. Reproduction

Install `requirements.txt`, then run `python -m unittest -v test_audit.py` and
`python audit.py --verify-record evidence.json --output fresh.json` from this
directory. Output must be fresh; the CLI refuses to overwrite frozen evidence.
The checks include independent exact Catalan/implicit series through order 8,
actual existing-map quadrature roots on explicitly manufactured inputs,
re-solved angular derivative controls, ordinary-integrand repair controls,
rounded-input and missing-pressure negative controls, and immutable provenance.
They do not numerically evaluate the selected enormous scales or turn the
retained analytical dependencies into an unconditional proof.
