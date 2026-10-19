# Global axial radial jets, pointwise pressure source and shifted axial viscosity

Status: `global-axial-pointwise-source-and-shifted-viscosity-envelope`.
Parent: `dfcf62462a6a10da6313be53638f3ab2b25d9583`.

This checkpoint supplies explicit global bounds for the leading pointwise
source `w0=Omega0/X`, the axial radial derivatives, V0 and the complete
shifted axial viscosity of U0. These are new inputs to the retained first
coefficient/stress estimate, conditional on the same exact core, loop,
matching-root and reference-stage lemmas. No gigantic parameter is sampled.
The unresolved coefficient-cutoff shape remains unresolved; no replacement
cutoff is introduced. All original supports, roots, frequency, width,
pressure gauges and both mean-pressure memories are unchanged.

Two distinctions matter. Signed source moment bounds were already available,
but could not bound the pointwise or absolute source. We now derive a
separate pointwise bound. Also, blindly differentiating the radial-average
formula twice and then in eta asks for a fourth mixed U derivative. The
exact divergence identity removes that unnecessary derivative demand.

## 1. Raise the same compact modulation to mixed third order

Retain `A=C^4096`, `Ahi=A^256`, `G=exp(A^32)`, `H=exp(G^4)`,
`N=1+floor(H^32)`, `C>1e9`. In particular `Ahi<H`, `H>2^1024` and
`H^32<N<2H^32`. For example `256 log A<=256A<=A^2<A^32<G^4` proves
the first comparison. All comparisons are symbolic, not evaluations of H.

The [higher-jet checkpoint](../higher_jet_budget/README.md) bounds the
ORIGINAL fields in mixed J4 by Ahi and the joint periodic primitives
script-A,script-B in J3(y,eta,phi) by H, on the actual compact annulus
`[Xa exp(tc/2),rho_p exp(5)]`. The same modulation is

`E_N=E exp(script-A(y,eta,Ny)/N)`, `U_N=U+script-B(y,eta,Ny)/N`.

For either primitive P define the weighted truncated mixed norm

\[
\|f\|_{3,N}=\sum_{r+k\le3}\frac{N^{-r}}{r!k!}
                   \sup|\partial_y^r\partial_\eta^k f|.
\]

It is a product algebra. For `beta=P(y,eta,Ny)/N`, the chain rule gives

\[
\frac{N^{-r}}{r!k!}\partial_y^r\partial_\eta^k\beta
=\sum_{a+c=r}\frac{N^{-a-1}}{a!c!k!}
     (\partial_y^a\partial_\phi^c\partial_\eta^k P)(y,\eta,Ny).
\]

Summing the twenty joint terms proves `|beta|_(3,N)<=H/N<1`.
The Banach exponential estimate yields `|exp(beta)-1|_(3,N)<=3H/N`;
conversion back to unweighted J3 costs at most N3. Thus both modulated
native fields E,U have J3 norm at most `4Ahi H N^2`. No fourth loop
derivative is needed, and the third-order radial N2 loss is NOT called small.

The SAME final leading repair has coefficient C_eta^3 norm
`<=2000H^16/N`. Its original factor is `e(eta)=e_star/(1+eta^2)`;
the reference value budget gives e_star<=A and `|e|3<16A<Ahi`.
The five width-.1 unit-mass bumps need sigma derivatives through order four,
already bounded by 1e100 in higher_jet_budget. Their radial J3 norm is less
than `4e104<C^12`. Even summing five supports gives

`|delta E,delta U|J3<=5*2000*Ahi*C^12*H^16/N < H^20/N <1`.

Uniqueness keeps the previously selected root. Consequently

\[
                 |\widetilde E,\widetilde U|_{J_3}
                      <H^4N^2<H^{70}.                 \tag{1}
\]

This is a new mixed-order consequence of the retained inputs, not a relabeling
of the earlier ordinary C2 or purely angular C3 estimate.

For the axial source a sharper anisotropic sub-bound is useful. The same
chain rule gives `|beta_y|2<=H` and `|beta_yy|1<=2NH`. Adding the original
U J3 bound and the repaired jets (at most 1 and 2 respectively) proves

\[
|U_y|_2\le A_{hi}+H+1<H^2,\qquad
|U_{yy}|_1\le2A_{hi}+2NH+2<H^2N,
\]
\[
                         |U_{yy}-U_y|_1<H^3N.          \tag{2}
\]

Only the coefficient repair, not this leading repair, has the still-unspecified
inner cutoff from the preceding checkpoint. Its missing shape is not silently
used as an input to (1) or (2).

## 2. Include every later axial reference stage

After the compact patch, the reference U is 4eta, then
`4eta[1-sigma(log(1+y)/Md)]`, zero, the prescribed axial pulse, and zero
after the pulse. The changes to angular energy, heat flow and its earlier
compensation do not change this U. This is the same regional cover used by
[global_axial_angular_bound](../global_axial_angular_bound/README.md).

For the taper, the first three derivatives of `log(1+y)/Md`, y>=0,
have absolute bounds 1,1,2. The analytical step bounds 32,1216,220000
from the preceding checkpoint give a taper J3 norm below 1e7.

On the pulse

\[
U=e_b f(\eta)e^{-(1/2+\lambda)y}
       [Amp(\eta)R_0(\lambda y)+\sum_{i=1}^2c_i(\eta)\beta_i(y)].
\]

Here `|f|3<16`, `|Amp|3<2`, `|c_i|3<=4lambda^20<1`, eb<1, at the
SAME roots already audited. The exact R0 is
`phi_b(xi)[1-sigma(xi-10)]`, `phi_b=int_0^xi sigma(v/.02)dv`.
On its three pieces, `phi_b=.02 int sigma` before .02,
`phi_b=xi-.01` afterwards, and R0=0 after 11. Global raw derivative
bounds through three are `(11,353,13440,3040000)`, so `|R0|J3<1e7`.
Composition with lambda<=.01 cannot enlarge these bounds.

The two original width-.3 pulse bumps are **unnormalized** translates of
sigma-prime (mass .3), distinct from the width-.1 unit-mass repair bumps.
The retained fourth-step-derivative bound 1e100 gives their radial J3 norm
below 1e103. The exponential factor has J3 norm below 3 for y>=0. Thus the
complete pulse J3 norm is bounded by

`3*16*(2e7+2e103)<1e110<H`.

This includes both affine end bumps; no infinite interval is multiplied by
a uniform energy estimate. Every nonzero piece has flat compatible joins.
The zero pieces and constant 4eta satisfy the same bound. Therefore (2),
with its larger budgets H2,H3N, holds on ALL radii after the axis-side region.
These estimates are analytical consequences of the reference formulas,
not sampling of the actual tiny lambda pulse.

Source: supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(A.7)-(A.10), SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its proposed blow-up theorem is not treated as an established result.

## 3. Global ordinary derivative inputs without dividing at the axis

Let Kaxis be the same common radial-order-0,1,2 / angular-C3 axis-side
envelope from global_axial_angular_bound (the activated-collar proof at
`ell_axis=tc/2`, followed by angular Cauchy). Retain

\[
K_U=\max(K_{axis},2A_{hi},352),\quad
R_1=\max(K_{axis},H^2/X_a),\quad
R_2=\max(K_{axis},H^3N/X_a^2).                         \tag{3}
\]

Then on the whole leading axial support `[0,Xv]`,
`|U|3<=KU`, `|U_X|2<=R1`, `|U_XX|1<=R2`.
Use the axis envelope before the modulation start; only afterwards use
`U_X=U_y/X`, `U_XX=(U_yy-U_y)/X^2`, with X>=Xa. There is no division
by zero and no assertion of radial analyticity across the flat activation.

The same functional terminal moment M_infinity=0 makes U,V and their fixed
jets zero after the original Xv. None of these supports is replaced by I3,
I4, Xcut or a later heat endpoint. These bounds extend by zero for the
axial fields beyond Xv, even though swirl and other pressure values persist.

## 4. Lower the pressure-source derivative demand exactly

Write `M=int_0^X U`, `v=V/X`, `L=1-2h eta^2`, `d=1-eta^2`,
`A=1/2+h`, `D=1/2-h`. The regular radial-average formula for v gives
`|v|2<=36KU`; it does not set v(0) to zero. Differentiating the exact
divergence identity, not estimating v_XX_eta from its raw average, gives

\[
LV_X=2A\eta U+2\eta XU_X-dU_\eta,
\]
\[
LV_{XX}=2(A+1)\eta U_X+2\eta XU_{XX}-dU_{X\eta}.       \tag{4}
\]

Taking one more angular derivative requires only U_XX_eta and U_X_eta_eta,
both TOTAL mixed order three, not U_XX_eta_eta. This is an exact identity;
it neither neglects a term nor adds an angular regularity premise.

Let `RL=(1-2h)^-1+4h(1-2h)^-2`, `v2=36KU`, `0<h<=.01`. Product
norms (`|eta|1=2`, `|d|1=3`) prove

\[
B_X=R_L[(4A+6)K_U+4X_vR_1],\quad
B_{XX}=R_L[(4A+10)R_1+4X_vR_2].                       \tag{5}
\]

These bound `|V_X|1`, `|V_XX|1`; also `|V|1<=Xv v2`.
The COMPLETE pointwise source is exactly

\[
w_0={V_X+D\eta v_\eta\over L}
 +v(V_X-v/2)+{U[-2\eta V_X+dv_\eta]\over L}-2\nu V_{XX}.
\tag{6}
\]

Retaining viscosity (nu=1 for the construction), both quadratic terms
and all angular derivatives gives the GLOBAL C_eta^1 envelope

\[
\boxed{W=R_L(B_X+4Dv_2)+v_2(B_X+v_2/2)
          +K_UR_L(4B_X+6v_2)+2\nu B_{XX}.}            \tag{7}
\]

Hence `sup_X |w0|1<=W`, `int_0^Xv |w0|1<=Xv W`, and
`int_0^Xv X|w0|1<=Xv^2 W/2`. These NEW absolute bounds may be much larger
than the earlier signed Q0,Q1, which remain available as signed repair inputs.
No signed cancellation has been promoted to a pointwise estimate.

Pressure memory is still nonzero where U=0 before the later pulse. On the
reserved plateau `M=mu eta`, `V=-mu`, equations (4),(6) give
`w0=-mu^2/(2X^2)`, not zero. Local zero axial velocity is not a source
support argument. The pressure value datum is unchanged, and the source
is integrated through the full future Xv whenever an incoming moment
requires it. No global viscosity cancellation is used in (7).

## 5. Complete shifted axial second derivative

For `Z_k=a partial_X+b partial_eta+c_k`, where
`a=-2eta X/L`, `b=d/L`, `c_k=2k eta/L`, expand

\[
\begin{aligned}
Z_{k-D}Z_kU={}&a^2U_{XX}+2abU_{X\eta}+b^2U_{\eta\eta}\\
 &+[aa_X+ba_\eta+a(c_k+c_{k-D})]U_X\\
 &+[bb_\eta+b(c_k+c_{k-D})]U_\eta
   +[b(c_k)_\eta+c_{k-D}c_k]U.
\end{aligned}
\]

For the actual axial exponent k=-1/2-h, both exponents have modulus below
2. The same estimates also apply to the angular exponent if its field
inputs are separately bounded. Since `|1/L|2<2`, coefficient C_eta^1
bounds are `|a|<=8X`, `|b|<=8`, `|c_k|,|c_(k-D)|<=16`,
`|a_X|<=8`, `|a_eta|<=16X`, `|b_eta|<=16`, `|(c_k)_eta|<=32`.
In particular derivatives of coordinate coefficients are not frozen.
Use `|U_X_eta|1<=2R1`, `|U_eta_eta|1<=6KU`, `|U_eta|1<=2KU` to get

\[
\boxed{C_U=64X_v^2R_2+704X_vR_1+1664K_U
         \ \ge\sup_X|Z_{k-D}Z_kU|_1.}               \tag{8}
\]

The first-stress source interface now has explicit axial inputs
`u=KU`, `ux=R1`, `v=Xv v2`, `w=W`, `CU` from (8).
They apply through Xc by zero extension of the axial source where necessary.
The swirl inputs f,fx,CF still need their full later-reference/heat-domain
estimates; their compact counterparts follow from (1), not automatically
from an estimate on a different radial domain. The coefficient-cutoff shape,
full Cstar, unchanged-tilt and combined scale gates remain separate.

## Checks, scope and next obligation

Exact symbolic identities verify (4), logarithmic conversion and the full
shifted operator. Ninety independent rational polynomial controls include
the axis, eta=+-1 and three viscosities; they compare (6) with the original
nonconservative source in value and its first angular derivative. The
moderate hierarchy and reference-shape arithmetic are auxiliary controls,
not the selected enormous profile. Historical replays use original comparators.

This supplies conditional analytical leading-source bounds, not an
unconditional physical pulse estimate or a complete wave residual gain.
The next obligation is the remaining global swirl/shifted-swirl inputs,
the actual cutoff-shape budget and verification of the full same-field
cone/scale inequalities, while retaining complete covariance and mean
pressure terms. Full residual closure and blow-up remain unproved.

```sh
python3 research/navier_stokes_cascade/results/global_axial_source_envelope/test_audit.py -v
python3 research/navier_stokes_cascade/results/global_axial_source_envelope/audit.py \
  --verify-record research/navier_stokes_cascade/results/global_axial_source_envelope/evidence.json \
  --output /tmp/fresh_global_axial_source_envelope.json
```
