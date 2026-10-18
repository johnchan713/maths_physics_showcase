# Explicit operator and full-source envelope on the actual activation collar

Parent: `1e910ccc6aa4da81eefb1d05e76e010943d0e91b`.
This checkpoint removes the unspecified LOCAL collar operator constant
`C_col` from [first_seminorm_bound](../first_seminorm_bound/README.md).
The new bound is an explicit finite expression in the retained core norm,
analytic radius, radial scale and unchanged activation width. It applies
to the ACTUAL activated functions, not just their analytic reference or
the old protected-core `C^16` estimate. The enormous actual parameters are
not numerically materialized. This is a conditional analytical estimate
relative to the retained exact core and regional construction, not an
independent verification of every inherited premise or a continuum proof
from passing tests. Frozen historical evidence is unchanged.

## 1. Domain, cutoff and reference bounds

Use the same `S>=2`, `0<rho<=1`, `Lambda>0`, `Xa=4/Lambda`,
`t=t1=C^-20`, `ell=tc/4<t/20`, `t<=.01`, and `delta=rho/(800S)`.
The actual collar is `0<=X<=Xa exp(ell)`, including the axis.
Later compact modulation begins at `y=tc/2`, outside this collar.
Put

\[
 x={4.1\over\Lambda},\quad b=1+x,\quad
 B=2(256S+9+S/\Lambda)(1+\Lambda)^2.                 \tag{1}
\]

Since `ell<.01`, `4 exp(ell)<400/99<4.1`, so `X<=x` and
`sqrt(X)<=b`. On the complex eta strip of radius `4delta`, the core
joint-series margin is at least `5/8`. Its explicit all-order majorant
and norms `256S`, `9+S/Lambda` give, for `r=0,1,2`,

`|partial_X^r F_r|, |partial_X^r U_r| <= B`.

Indeed `r!20^-r(8/5)^(r+1)<=2` for these three orders, and the radial
rescaling costs `Lambda^r<=(1+Lambda)^2`. The same majorant gives
`|Phi_Y|,|Phi_YY|<=S`; the retained zero-free neighborhood gives
`|1/Phi|<5`. These are complex bounds, not real positivity alone.
Consequently, for `a_r=-2Y Phi_Y/Phi`,

\[
 |a_r|\le50S,\qquad |\partial_y a_r|\le2000S^2.     \tag{2}
\]

For example the second bound follows by differentiating the exact
quotient, using `Y<=4.1` and `S>=2`; no bare angular derivative is used.

Keep the actual SQUARED flat cutoff

`sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2))`.

Its derivative is `2 sigma(1-sigma)[s^-3+(1-s)^-3]`.
For `0<s<=1/2`, `sigma<=exp(4-1/s^2)`, and
`s^-3 exp(-1/s^2)<=8 exp(-4)`. The other denominator costs at most 8.
Thus `|sigma'|<=32<64`; reflection covers the other half, and flat
endpoint limits are zero. For `e_a=(1-t)sigma(y/t)`, this proves
`0<=e_a<=1` and `|partial_y e_a|<=64/t`. A cutoff with `exp(-1/s)`
would be a different prescription and is not substituted here.

## 2. Quantify the actual activation, including second radial derivatives

Retain the exact forward integrals from activated_collar_foundation:

\[
 J={1\over2}\int_0^y e_a a_r\,dv,\quad
 K=\int_0^y e_a\partial_vU_r\,dv,\quad
 F_a=F_r e^J,\quad U_a=U_r-K.                       \tag{3}
\]

On the outer `4delta` strip define

\[
 j_1=25S,\quad j_2={1600S\over t}+1000S^2,\quad E=e^{25\ell S}.
\]

Then `|J|<=25ell S`, `|J_y|<=j1`, `|J_yy|<=j2`.
The exact product rule retains `2(F_r)_y J_y`, `J_yy` AND `J_y^2`.
The axial activation obeys `U_a,y=(1-e_a) U_r,y`, so its second
derivative retains `-(e_a)_y U_r,y`. Using
`partial_X=X^-1 partial_y` and
`partial_X^2=X^-2(partial_y^2-partial_y)` gives the following bounds:

| Actual field | Value | First X derivative | Second X derivative |
|---|---|---|---|
| `F_a` | `EB` | `EB(x+j1)/Xa` | `EB[2x+x^2+2xj1+j2+j1^2+j1]/Xa^2` |
| `U_a` | `B(1+ell x)` | `Bx/Xa` | `B[2x+x^2+64x/t]/Xa^2` |

The subtraction of the lower logarithmic derivative is essential.
A convenient common bound for both fields and all three radial orders is

\[
 \mathcal K=EB(1+\ell x)\left[1+{x+j_1\over X_a}
 +{2x+x^2+2xj_1+j_2+j_1^2+j_1+64x/t\over X_a^2}\right].       \tag{4}
\]

On `X<=Xa` use the unchanged core bound `B<=K`; the division by `Xa`
in (4) is used only to bound the activation side, never division by zero.
Flatness at `Xa` joins these actual bounds smoothly. This proof needs no
radial analyticity at the flat join. In particular the displayed `1/t`
cost cannot be dropped or made uniform as `C` varies.

For every real radial point, Cauchy on an eta disk of radius `delta`
inside `4delta` bounds the mixed derivatives `r<=2, k<=2` on `2delta` by

\[
 |\partial_X^r\partial_\eta^k(F_a,U_a)|\le
 \mathcal H:=\mathcal K(1+\delta^{-1}+2\delta^{-2}).            \tag{5}
\]

Differentiating radial averages under the integral gives the SAME bound
for `Ubar=int_0^1 U_a(sX,eta)ds`; the factor `s^r<=1` is retained.
This includes the regular axis values and the mixed `XX,eta` jet needed
below. No angular operator is declared bounded on the original Banach norm.

## 3. Full radial pressure source and the six-variable system

Write `A=1/2+h`, `D=1/2-h`, `d=1-eta^2`, `L=1-2h eta^2`,
`0<h<=.01`, and retain the construction viscosity equal to ONE.
On `2delta`, `|eta|<1.01`, `|d|<2.021`, `|L|>.97`.
The actual regular radial velocity is

\[
 v={2\eta U_a-2D\eta\overline U-d\overline U_\eta\over L}.
\]

Direct differentiation gives
`|v|,|v_X|,|v_XX|,|v_eta|<=12H`. For the angular derivative, the
numerator costs less than `10.071H`, plus the derivative of `1/L`;
the latter is explicitly included. The full regular radial source is

\[
 w_0={v+Xv_X+D\eta v_\eta\over L}
 +v(v/2+Xv_X)
 +{U_a[-2\eta(v+Xv_X)+dv_\eta]\over L}
 -2(2v_X+Xv_{XX}).                                           \tag{6}
\]

Every term, including radial viscosity and both quadratic terms, is kept.
Since `H>=1`, the triangle inequality yields

\[
 |w_0|\le512(1+x)\mathcal H^2=:\mathcal W.                    \tag{7}
\]

For the retained operator convention
`Z_k=(-2eta X partial_X+d partial_eta+2k eta)/L`,
the full shifted second axial derivative is `Z_(k-D) Z_k`.
The two relevant exponents are `k=-1-h,-1/2-h`, so both it and its
shift have modulus below 2. Its coefficient bounds are
`|a|<=3X, |b|<=3, |c|<=6`,
`|a_X|<=3, |a_eta|<=3X, |b_eta|<=3, |c_eta|<=6`.
Applying the product rule rather than freezing coefficients proves

\[
 |Z_{k-D}Z_k(F_a,U_a)|\le128(1+x)^2\mathcal H=:\mathcal Z.   \tag{8}
\]

In the exact retained system
`W_xi+diag(0,0,2,0,3,1)W/xi=A0 W+A1 W_eta+f`,
`W=(F1,U1,K1,Pi1,P,Q)`, the row-sum and source bounds are

\[
 \|A_0\|_\infty\le128(1+x+b)\mathcal H,\quad
 \|A_1\|_\infty\le32(1+x)\mathcal H,\quad
 |f|_\infty\le(b+3x)\mathcal W+2\mathcal Z.                  \tag{9}
\]

For reproducibility, rows five and six of `A0` cost respectively at most
`(73+24x+16b)H` and `(40+33x+16b)H`; the other rows are smaller.
The two nonzero `A1` rows cost at most `(17+10x)H` and `(12+10x)H`.
These estimates bound the actual expressions in
`first_correction_construction/system.py`, not a different model system.
The full source is `(0,0,0,-xi*w0,-2*CF,2eta X*w0/L-2*CU)`.
In particular the pressure entry and both shifted axial terms are nonzero
in general and are not erased. Combining (7)--(9) proves the explicit input

\[
 \boxed{C_{col}=4096(1+x+b)^3\mathcal H^2}.                  \tag{10}
\]

This bounds both matrix row norms and the vector sup norm on the common
complex collar. It is not claimed to be sharp, `C^16`, or small.

Pressure VALUE memory also remains: let
`Pcore=2(10P^2+80*256^3*S^2/Lambda)` from the exact pressure-family norm.
The full local pressure is
`Pi_a=Pi_r(Xa)+int_0^y X F_a^2 dv`, not a zero-incoming primitive.
Thus `|Pi_a|<=Pcore+ell*x*K^2` on the outer strip, including the core;
`|Pi_a,X|<=K^2`, `|Pi_a,XX|<=2K^2`.
Angular Cauchy factors in (5) apply as well. The inherited `Pi_axis`
is still part of `Pcore`, and the global backward-pressure identification
is not replaced by a local estimate. Although (6) needs no additive
pressure value, that fact is not permission to reset the datum.

## 4. Join to the existing first-coefficient estimate; remaining gate

Use exactly the previously proved sparse Picard estimate with
`sigma=2delta`, `Delta=d=delta/4`, and the endpoint enclosure `a_col<=b`:

\[
 z=2C_{col}b,\quad \beta={e z^2\over\Delta},\quad
 M=(z+\beta)e^\beta.                                        \tag{11}
\]

The strip condition `Delta+d<sigma` is strict. Substituting these finite
explicit inputs in the parent's `collar_jets` yields its value, first and
second radial coefficient bounds. The singular diagonal and axis values
remain in the retained integral system; (11) does not assume a radial
Taylor expansion through activation. For exact diagnostic arithmetic,
`exp(25ell S)<=2^ceil(50ell S)` follows from `log 2>1/2`;
the implementation stores the exponent separately, never materializing
the actual huge exponential or confusing floating underflow with zero.

This discharges the previously free LOCAL `C_col` and strip inputs,
conditional on the retained exact core. It does NOT yet close the
complete first-coefficient `Cstar`: its `W0,W1` are full future integrals
through `Xv`, not the collar. Bounds on the global leading jets, full
shifted axial sources, actual radial cutoff and exact repair bumps still
need to be substituted. The SAME original tilt must separately pass
`u^2>(1-epsilon_star)/epsilon_star`, or a sharper same-u test; a smaller
physical scale does not fix that gate. Only after those obligations may
one test

`2h(ell_scale-2)>=log2(Cstar/r_Kc)`

together with freezing, geometry, inverse and all other acceptance gates.
No historical first-band, unconditional actual-pulse bound, complete
residual gain or blow-up conclusion is obtained here. Complete-wave
covariance and mean pressure/axial-velocity obligations remain unchanged.

## Checks

Exact rational controls check cutoff constants, common-strip geometry,
width losses and exponent-separated arithmetic. Symbolic product and
logarithmic chain rules retain the lower derivative and squared-J term.
Nonnegative polynomial coefficients verify domination in (10) for all
`x,b>=0,H>=1`. Eighteen independent retained-system evaluations include
the axis and both angular endpoints; they are manufactured diagnostics,
not samples of the actual infinite core or proofs of its quantifiers.
Relevant parent tests and frozen evidence are replayed under their
original tolerances; historical records are not rewritten.

```sh
python3 research/navier_stokes_cascade/results/collar_operator_envelope/test_audit.py -v
python3 research/navier_stokes_cascade/results/collar_operator_envelope/audit.py \
  --verify-record research/navier_stokes_cascade/results/collar_operator_envelope/evidence.json \
  --output /tmp/collar_operator_envelope_fresh.json
```
