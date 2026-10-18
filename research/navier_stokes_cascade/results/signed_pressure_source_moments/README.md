# Full-support signed pressure-source moments and an axis-viscosity term

Parent: `a70349e4b48b41c08033921824951ce286766195`.
The activated-collar operator input is already explicit. This checkpoint
advances the still-open GLOBAL pressure-source input to the first-coefficient
repair. It derives a conservative radial identity and exact formulas for
the two signed source moments actually needed by that repair. They eliminate
global second radial derivatives from THIS moment bound, retain the complete
future axial support and expose a nonzero axis-viscosity contribution.

This sharpens, rather than invalidates, the earlier sufficient absolute-source
bound in [first_seminorm_bound](../first_seminorm_bound/README.md).
Its stronger integrals `int |w0|_1` and `int X|w0|_1` are not necessary for
the two signed repair discrepancies. The new estimate is conditional on the
same finalized leading profile, its regular axis, the exact functional total
`M_infinity=0` and compact axial/radial support from the current same-root
pressure/moment join. No historical evidence, gauge, width, root, frequency
or support is changed. It does not quantify the global angular input from
every original construction stage, nor establish the complete first-stress
bound, physical pulses, full residual improvement or blow-up.

## 1. Exact conservative identity, including the axis

Write `U=U0`, `M=int_0^X U`, `V=V0`, `v=V/X`, and
`A=1/2+h`, `D=1/2-h`, `d=1-eta^2`, `L=1-2h eta^2`.
The regular average formula defines v at the axis. The retained divergence
identity is

\[
 V={2\eta XU-2D\eta M-dM_\eta\over L},\qquad
 LV_X=2A\eta U+2\eta XU_X-dU_\eta.                    \tag{1}
\]

For general nonnegative viscosity nu (the construction keeps nu=1), the
FULL radial source from the retained six-variable system is

\[
 w_0={\Omega_0\over X}
 ={v+Xv_X+D\eta v_\eta\over L}
 +v(v/2+Xv_X)
 +{U[-2\eta(v+Xv_X)+dv_\eta]\over L}
 -2\nu(2v_X+Xv_{XX}).                                 \tag{2}
\]

Using `V=Xv`, (1), and the product rule gives the POINTWISE identity

\[
 \boxed{w_0=\partial_X\mathcal J+
 {D\eta v_\eta+d\partial_\eta(Uv)-2A\eta Uv\over L}
 +\tfrac12v^2,\qquad
 \mathcal J={(1-2\eta U)V\over L}+Xv^2-2\nu V_X.}     \tag{3}
\]

For example the mixed term is rewritten by
`d U_eta v=2A eta Uv+2eta U_X V-L V_X v`.
Also `v V_X=v^2/2+partial_X(Xv^2/2)`.
Together these identities retain both nonlinear orientations; discarding
one produces a wrong conservative flux. There is no division by d, so
eta=+/-1 are included with the retained one-sided angular jets.

At the axis `V(0)=0`, `Xv^2->0`, but `V_X(0)=v(0)` is generally NONZERO.
Thus `J(0)=-2nu v(0)`. Axis regularity alone does not erase this flux.

## 2. Full-support moments and exact finite-endpoint corrections

The current exact functional moment gives `M=0` after the leading axial
support Xv; there U,V and all their fixed jets vanish. Use that SAME Xv,
not Xcut, Xvel, I3 or I4. Define full-support signed integrals

\[
 a_j=\int_0^{X_v}X^jv\,dX,\quad
 b_j=\int_0^{X_v}X^jUv\,dX,\quad
 c_j=\int_0^{X_v}X^jv^2\,dX\quad(j=0,1).
\]

Integrating (3), and then X times (3), proves

\[
 \boxed{\mathcal I_0:=\int_0^{X_v}w_0\,dX
 =2\nu v(0)+{D\eta(a_0)_\eta+d(b_0)_\eta-2A\eta b_0\over L}
 +\tfrac12c_0,}                                      \tag{4}
\]
\[
 \boxed{\mathcal I_1:=\int_0^{X_v}Xw_0\,dX
 ={ -a_1+D\eta(a_1)_\eta+d(b_1)_\eta+2D\eta b_1\over L}
 -\tfrac12c_1.}                                      \tag{5}
\]

The first viscosity total is `-2nu[V_X]_0^Xv=2nu v(0)`.
The second is `-2nu[X V_X-V]_0^Xv=0`.
This cancellation is for the WEIGHTED TOTAL, not for the pointwise source,
its absolute integral, local pressure or the complete residual.
Compact support and uniform fixed jets justify integration/differentiation
in eta. These formulas retain all later-source memory algebraically.

For an arbitrary positive finite endpoint R, use the same a,b,c integrals
only through R. The right sides of (4),(5) then require respectively

\[
 B_0(R)={(1-2\eta U_R)V_R\over L}+{V_R^2\over R}
                         -2\nu(V_X)_R,               \tag{6}
\]
\[
 B_1(R)={R(1-2\eta U_R)V_R\over L}+V_R^2
                    -2\nu\{R(V_X)_R-V_R\}.           \tag{7}
\]

The axis term in (4) is still required. Before the later axial pulse,
`U=0,V=-mu` does not make (6),(7) zero. Neither a local zero-velocity
coefficient nor vanishing M at one angle permits truncating the future
source or erasing the eta derivatives of its moments. A nonzero final
functional M defect would also invalidate the compact-V endpoint premise;
the same-root global moment identity remains indispensable.

## 3. Explicit angular-norm bounds without global radial derivatives

Use `|f|_k=sum_(r<=k) sup|partial_eta^r f|/r!`, a product algebra.
It obeys `|f_eta|_1<=2|f|_2`, `|f_eta|_2<=3|f|_3`.
For real eta in [-1,1], let

`l=1-2h`, `R_L=1/l+4h/l^2 >= |1/L|_1`.

Suppose the actual leading fields satisfy
`|U(X,.)|_2<=u2`, `|v(X,.)|_2<=v2` for EVERY X in [0,Xv].
Then `|eta|_1=2`, `|d|_1=3` and (4),(5) give

\[
 |\mathcal I_0|_1\le
 4DR_L X_vv_2+(6+4A)R_L X_vu_2v_2
                   +\tfrac12X_vv_2^2+2\nu v_2=:Q_0, \tag{8}
\]
\[
 |\mathcal I_1|_1\le
 {R_L X_v^2\over2}\{(1+4D)v_2+(6+4D)u_2v_2\}
                       +\tfrac14X_v^2v_2^2=:Q_1.     \tag{9}
\]

One can remove the separate v2 input. Define the single actual global
zero-radial/third-angular seminorm

\[
 K_U=\sup_{0\le X\le X_v}|U(X,\cdot)|_3.              \tag{10}
\]

The regular average `Ubar=int_0^1 U(sX,eta)ds` has the same angular
bound. The exact reciprocal second-jet bound is

`|1/L|_2 <= 1/l+6h/l^2+16h^2/l^3 < 2` for `0<h<=.01`.

All positive terms in this majorant increase with h, so the entire
interval is bounded by its endpoint value `127600/117649<2`, not just
by two diagnostic samples.

Using `|eta|_2=2`, `|d|_2=4`, the numerator of v costs at most
`(4+4D+12)K_U`. Thus

\[
 u_2\le K_U,\qquad v_2\le(16+4D)|1/L|_2K_U\le36K_U.   \tag{11}
\]

Substitution of `u2=K_U,v2=36K_U` in (8),(9) is an explicit
full-support bound for the two signed pressure-source inputs. It requires
NO global X derivative of U or v. This is not a bound for
`int |w0|_1`, a pointwise w0 estimate, or the shifted axial source in the
compact stress equation. Those distinct obligations are not discharged.

The previous actual-collar estimate supplies the local contribution to
(10): if K is its common outer-strip value bound and delta its retained
strip, angular Cauchy gives
`|U|_3 <= K(1+delta^-1+delta^-2+delta^-3)` there.
But the rest of the SAME finalized profile, including the later axial
pulse, still needs a quantitative bound from its original stage parameters.
Finiteness from the retained smooth compact profile is not a numerical
evaluation or permission to claim the complete first coefficient closed.

## 4. Feed the signed totals to the unchanged exact linear repair

From the exact first-pressure equation, with the same zero axis gauge,

`Pi1_X=2F0 F1-w0/2`.

For the tentatively cut first velocity coefficients the signed source
contributions to m3 and m5 are therefore `-I0/2` and `+I1/2`.
If J bounds their collar C_eta^1 values and f,u bound the leading values,
the ordinary moment bounds can be sharpened to

\[
 (\mu_1,\mu_2,\mu_3,\mu_4,\mu_5)\le
 (X_{cut}J,X_{cut}^2J,2fX_{cut}J+Q_0/2,
 X_{cut}^2(u+f)J,2uX_{cut}J+fX_{cut}^2J+Q_1/2).         \tag{12}
\]

These are bounds AFTER full signed integration. The algebraic triangle
step is valid because the repair acts on integrated ordinary discrepancies,
not on their source density. There is no assertion `Q_j>=int X^j|w0|_1`.
Insert (12) into the parent's SAME inverse-800 repair formula; all
`1/lambda_patch`, patch radii and reciprocal-amplitude angular costs stay.
The later sources can still affect an earlier repair coefficient through
(4),(5). No incoming pressure constant is reset.

Neither the global zero repaired m3 nor (5)'s viscosity cancellation sets
local I4 pressure to zero: `Pi1=K1(eta)-mu^2/(4X)` and `Pi2=K2(eta)`
remain exactly as in coefficient_support_foundation. Complete covariance
means must still retain them and the axial velocity terms.

## 5. Checks, status and next obligation

The independent Fraction-only radial polynomial/angular Taylor calculus
checks (2) against (3) and integrates it directly, without quadrature.
Eighty cases cover two nonseparable fixtures, both angular endpoints,
two h, two viscosity values, and full versus truncated supports.
Both moments and their first angular derivatives give 320 exact scalar
identities. Negative controls detect the missing axis-viscosity term and
missing future-endpoint corrections. The compact polynomial fixtures have
enough fixed regularity for these identities but are NOT the actual
infinite leading profile or a proof of its existence. Eighteen separate
floating comparisons against the historical regular-source implementation
guard sign and convention changes. Passing tests are not a continuum proof.

An additional exact manufactured control separates signed and absolute
integrals: a certified coefficient-triangle bound for K_U feeds (9), while
the absolute signed integral through R=1/3 already exceeds Q1 at auxiliary
viscosity `10^30`. Hence the full absolute weighted integral exceeds Q1,
even though the full SIGNED moment satisfies it. This is an algebraic
counterexample to misusing (9) as an absolute bound, not an evaluation of
the construction, whose viscosity remains one.

Next quantify K_U for all later defining stages, then the separate global
leading/shifted-axial/cutoff/repair-bump inputs to Cstar. Verify the SAME
unchanged tilt and all scale/freezing/inverse gates before any actual-pulse
or complete residual claim. A smaller physical scale cannot by itself
verify the tilt condition. No old first-band or blow-up certificate is made.

```sh
python3 research/navier_stokes_cascade/results/signed_pressure_source_moments/test_audit.py -v
python3 research/navier_stokes_cascade/results/signed_pressure_source_moments/audit.py \
  --verify-record research/navier_stokes_cascade/results/signed_pressure_source_moments/evidence.json \
  --output /tmp/signed_pressure_source_moments_fresh.json
```
