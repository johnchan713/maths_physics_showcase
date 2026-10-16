# Core growth and the quantitative activation collar at the frozen width

Status: `core-growth-and-quantitative-inner-activation-audited`.

This batch supplies the positive shear/growth input left open by the core
contraction audit and bounds the activation constant left abstract by the
inner flat-quotient audit. Two different mechanisms give growth: angular
shear where chi>.99, and pressure-driven axial shear where chi<=.99.
No division by chi is used. The same fixed reference datum and the same
Lambda,C and t1=kappa0=C^-20 are retained.

The result is local to the leading core and its first activation collar.
The unedited reference defines the axis pressure exactly; identifying the
completed global backward pressure and all coefficient supports still
requires the later exact moment identities. Compact modulation/repair,
unconditional actual-pulse bounds, complete residual improvement and
blow-up remain unverified. No enormous profile is sampled or substituted
for the analytical core, and this is not a proof-assistant certificate.

Parent: `3108ba7d550f074afc1199188609c898195b525d`.
Use the [core contraction](../core_contraction_audit/README.md) and
[inner activation](../inner_collar_activation/README.md) at their frozen
bytes. Source: the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(A.21), (B.17)--(B.31), (4.16), and (7.1), SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its final theorem is not a premise. All inherited evidence and GOAL.md
remain unchanged. The new constants below are derived rather than fitted.

## 1. Identify the reference axis datum and its indispensable lower part

The axis datum is DEFINED before the core by the unedited reference
schedule: Pi0=-(1/2)*integral E_id,sched(y,eta)^2 dy. Every stage has
`E=c(y)*f(eta)^theta(y)`, f=(1+eta^2)^-1, 0<=theta<=1. Thus

\[
 -\Pi_0/P^2=\int_0^1 f^{2\theta}\,d\mu(\theta),\qquad
 \mu\ge(5/2)\delta_1.                                      \tag{1}
\]

The atom follows by integrating the exact inner power
`E=P*f*exp(y/10)` on y<=0. For the first unit transition, c/P<=exp(.1)
and c(1)/P=exp(-.2), since the flat step has integral 1/2. Subsequent
native logarithmic slopes are nonpositive. Parameter interpolation adds
a nonpositive derivative of the factor 2^(-(1-theta)). The terminal cutoff
ratio is below two. Therefore, for y>=1,
`c/P<=2*exp(-.2)*exp(-(y-1)/2)`. Integrating the three ranges gives

\[
 \mu([0,1])\le(5+e^{.2}+4e^{-.4})/2<5.                     \tag{2}
\]

This identifies the datum used by the inner construction without evaluating
the enormous stage lengths. It verifies the pressure-family input to the
previous local core theorem. It does not verify that all subsequent edits
preserve the global pressure normalization; that is a different identity.
The lower atom in (1) is additional information, not a consequence of the
upper mass bound alone.

Write A=.5+h, D=.5-h, d=1-eta^2, L=1-2h*eta^2, j=1e-18,
sigma=j/100, U*=4eta+j, H*=Deta+dU*. The pressure contribution to Z* is
`-eta*P^2*Kp(eta)`, where

\[
 K_p=\int[4d\theta/(1+\eta^2)+4A]f^{2\theta}\,d\mu\ge5/4.
\]

The remaining velocity polynomial is exactly

\[
 R_v=-20\eta+(16+32A)\eta^3-(A+4)j+(4+16A)j\eta^2+2Aj^2\eta.
\]

For |eta|<=1, 0<h<=.01 and j<=.05, its magnitude is bounded by
`40*|eta|+5*j`. Put delta*=j*P^2/100. If |Z*|<=delta*, then

\[
 |\eta|/j\le\frac{1/100+5/256}{5/4-40/256}<1/30.             \tag{3}
\]

On that smaller interval,
`H*/j>=1-4.5/30-4j^2/30^3-j^2/30^2>.8`. Consequently
`chi=H*^2/(H*^2+sigma^2)>=6400/6401>.99`. Taking the contrapositive,

\[
                   \chi\le.99\ \Longrightarrow\ |Z_*|>\delta_* . \tag{4}
\]

This is a whole-interval separation bound, including the tiny neighborhood
of the zero of H*. The diagnostics resolve eta/j explicitly; binary64
sampling of eta would miss this region.

## 2. Positive angular source with explicit error bounds

Retain the previous core constants E=Eb, S,Fb,K,Lambda,M and put
`Ubar=E*S`, `e=E*Fb/Lambda`. These bound the normalized u jets and
the error Phi-f0(Y*chi), respectively, through two Y and six eta
derivatives on 0<=Y<=4.1. They are bounds for the analytic limit, not
errors inferred from a finite Taylor truncation. In particular
`1/4<Phi<2`, `e<1e-50`, `Ubar/Lambda<1e-50`.

The alternating series gives `|f0'|<=1/4`, f0>=.265 and f0<=1 on
[0,4.1]. The identities
`|H* chi_eta|<=2*53*chi`, `|chi_eta|<=106/sigma` imply

\[
 |D_Y\log f_0(Y\chi)|\le4\chi,\qquad
 |H_*\partial_\eta\log f_0(Y\chi)|\le410\chi.
\]

Use the reciprocal-difference formula with Phi>=1/4 to get

\[
 |D_Y\log\Phi-D_Y\log f_0|\le32e,\quad
 |H_*\partial_\eta\log\Phi|\le500\chi+2000e,\quad
 |\Phi_\eta/\Phi|\le500/\sigma.                             \tag{5}
\]

For the weighted difference in (5), the errors are at most
`4*64*e+4*410*e=1896e`; the 2000 bound leaves slack. These calculations
retain the vanishing H* factor before estimating a parameter derivative.

The exact source is
`Sq=-W*(1+D_Y log Phi)-h*(1-2eta U)-Hc*(Lambda*zeta*+Phi_eta/Phi)`.
Here `|W-W*|<=2*Ubar/Lambda`, |W|<=5, -W*>2.9, |U|<4.1 and
`|zeta*|<=sqrt(chi)/sigma` on the real interval. Substitution into this
source, with no missing pressure contribution, yields

\[
 S_q\ge\Lambda L\chi+2.8-520\chi-(\bar U/\sigma)\sqrt\chi
      -2160e-2\bar U/\Lambda-500\bar U/(\Lambda\sigma).
\]

Absorb 520chi in .02*Lambda*L*chi. Young's inequality gives

\[
 (\bar U/\sigma)\sqrt\chi\le.03\Lambda L\chi
                            +9\bar U^2/(\Lambda\sigma^2).
\]

The remaining constant loss is

\[
 \mathcal E=2160e+2\bar U/\Lambda+500\bar U/(\Lambda\sigma)
                                +9\bar U^2/(\Lambda\sigma^2). \tag{6}
\]

There is an analytical bound for ALL selected pressure scales. Set
H=K^(1/16)>1. The parameter definition gives
`Ubar<=H/1e6`, `E*Fb<=H/1e6`, and `sigma^-1<=M<=H/1e6`.
Hence

\[
 \mathcal E\le[(2162/10^6)+(500/10^{12})+(9/10^{24})]H^4/\Lambda
             <3\,10^{-47}<.1.                              \tag{7}
\]

Also Lambda>1e44 and L>=.98 make the linear absorption valid. Therefore

\[
                        S_q\ge2.5+.95\Lambda L\chi.        \tag{8}
\]

This bound remains valid at chi=0; a constant error estimate divided by chi
would not establish it.

## 3. Positive shear and the two exit-growth mechanisms

The exact core equations and regular axis give
`p1=a=-2Y*Phi_Y/Phi`, `ns=-2u_Y`. The positive forward integral gives

\[
 p_1/X=Q_s/L=\frac1L\int_0^1 v\frac{\Phi(vY)}{\Phi(Y)}S_q(vY)\,dv
                    \ge5/32>1/10.                          \tag{9}
\]

The same scalar estimates imply p1<=10 throughout the core. At Y=4 and
chi>.99, put z=4chi, t=z/2 in [1.98,2]. The decreasing alternating
polynomial
`1-t+t^2/4-t^3/36+t^4/576` is at most
`-75535511/400000000<-.188`. It bounds `f0+z*f0'` from above, so
the comparison shear exceeds 2.376. Its actual error is <=64e<.001.
Thus p1>2.3 on this sector.

On chi<=.99, (4) and `|ns-Z*/L|<=2e<delta*/2` give
`|ns|>=delta*/2`. At X=4/Lambda,
`E=sqrt(2X)*phi* Phi/C<=2*sqrt(2X)*exp(2Lambda M)/C`. Consequently

\[
 \frac{p_2^2}{p_1}\ge\frac{C^2\delta_*^2}{80\Lambda e^{4\Lambda M}}>2.3. \tag{10}
\]

The last inequality holds with the frozen log C=10000B+100log(1/epsilon),
B=1e20(1+Lambda+K)^64. Indeed M<Lambda, -2log(delta*)<100 for P>=16,
and Lambda>1e44; the logarithm of the right side is greater than
`20000*Lambda^64-4*Lambda^2-log Lambda-105>log(2.3)`.
There is no numerical evaluation of C or of an unnormalized tiny velocity.
Combining the sectors gives, for the exact core at its exit,

\[
 a_r\ge2/(5\Lambda)>0,\qquad
 v_r=a_r+q_r^2/a_r>2.3,\qquad q_r=-2U_y/E.                 \tag{11}
\]

These are the previously missing actual local core gaps for the fixed
reference datum. They are not a claim about the compact modulated annulus.

## 4. Bound the activation constants before choosing a smaller collar

Gauge the swirl by `Ebar=C*Er`, `Hbar=C*Hr` and `qbar=qr/C`.
Then `Er*qr=Ebar*qbar=-2Ur_y` has no C factor. The core estimates are
uniform in subsequent C: normalized Phi,u vary with C but their bounds
do not. Put J=1+Lambda+K, B=1e20*J^64 and R=exp(B/100).

Here is a common envelope on 4<=Y<=4.1. Norms include up to two eta
derivatives and one y derivative, wherever the displayed formula requires
it, with one additional eta derivative of a primitive when a stress formula
requires it. The six available eta jets suffice. A reciprocal is only
used on this collar, not at X=0.

| Basic quantities | Sufficient envelope before R |
|---|---|
| Phi,u jets; Phi reciprocal jets | J^6 |
| phi* and its reciprocal, parameter jets | exp(3J^2) |
| X, X^-1, L, L^-1, sqrt(2X) and reciprocal | J^2 |
| Ur, Ebar, Hbar; only Ebar and Hbar reciprocals | exp(4J^2) times J^6 |
| ar and ar reciprocal jets | J^40 |
| qbar and its y/eta jets | exp(5J^2) times J^6 |
| M, C*I, C*Jmoment, S, C^2*Cp, Pi0, W, Ns, C*NQ | exp(12J^2) times J^40 |

For example ar^-1<=5Lambda/2 by (9); differentiating a reciprocal up to
total order three gives at most four inverse factors and three ar-jet
factors. The bounds ar-jets<=J^10 therefore give J^40. The qbar_y
formula only needs u_YY, not a third Y derivative. For the pressure
primitive use `C^2*Cp_X=(phi* Phi)^2`, which is regular at the axis;
no bound on Ebar^-1 is integrated through zero. Forward integrals have
length X<1, and their differentiated integrands give the last row.
These elementary products are all below exp(20J^8), which is below R.
This explains the envelope and does not assume it for the completed flow.
No reciprocal of Ur is taken: the axial velocity can vanish on the collar.

Use the previous exact flat Volterra operator on 0<=s=y/t1<=1/4.
Its norm on values and eta jets is <=3/2, since 2<=g_flat<=6 in its
positive integral formula. For `t1*R^3<=.1`, exprel and the exponential
of the logarithmic field change have C1 norm <=2. The C1 norm here is
sup|f|+sup|f_eta|; it is a product algebra and the cutoffs have no eta
dependence. Apply the COMPLETE moment formulas with the same Pi0.
The following conservative counts give the needed coefficients:

| Quantity after removing the specified factor | Upper bound |
|---|---|
| Delta log E / (t1*ea*s^3) | R |
| Delta U / (t1*ea*s^3) | R^2 |
| C*Delta E, C*Delta H / (t1*ea*s^3) | R^4 |
| Delta M, C*Delta I, C*Delta J, Delta S, C^2*Delta Cp / (t1^2*ea*s^6), C1 norm | R^8 |
| Delta W, C*Delta NQ / (t1^2*ea*s^6) | R^10 |
| Delta p1 / (t1*ea*s^3) | R^18 |
| Delta Ns / (t1*ea*s^3) | R^16 |
| Delta p2 / (C*t1*ea*s^3) | R^24 |
| (Er/E-1) / (t1*ea*s^3) | R^2 |

For Qs, use the exact difference
`-Delta W+Delta NQ/(XH)-NQ,r*Delta H/(XHr H)`; for Ns, keep
`4A*eta*Delta Cp-d*Delta Cp_eta` as well as every U,W and M,S term.
The p2 difference is
`(X/L)*(Delta Ns/E-Ns,r*Delta E/(Er*E))`.
Thus only p2 and ts carry one positive C factor, while pressure differences
carry C^-2. The table bounds the surviving C^-2 terms; it does not delete
them. For the moment coefficients, the J row is the largest: its three
terms are bounded by R^5, R^3 and t1*R^6<=.1*R^3. Multiplication by
X<1 and the Volterra norm leaves a bound below R^8, also in C1. The
W and NQ differences use one parameter derivative of those moments and
one X reciprocal, giving R^10. The reciprocal-difference term in Qs
costs at most R^8 before harmless sums and X/L factors; R^18 suffices.
The Ns terms use R^16, and multiplication by Ebar^-1 and X/L costs at
most two further basic factors; R^24 suffices for p2/C. These counts
use no reciprocal of the pressure or axial velocity and no inverse
flat factor beyond the exact Volterra identities.

The cancellation-safe cone identities of the preceding activation note
then give, with coefficients Dc,Jc*,Vc,

\[
 (P_c-v_s)/e_a=v_r+t_1s^3D_c,\quad J_c/e_a=t_1s^3J_c^*,\quad
 (v_s-kv_r)/e_a=t_1s^3V_c,
\]
\[
 |D_c|,|J_c^*|,|V_c|\le K_{act}:=C^2R^{80}<C^3.             \tag{12}
\]

For example `(ts/C)*(Delta p2/(C*t1*ea*s^3))` costs at most R^27;
the other terms cost less than R^40. R^80 leaves room for all sums.
This is an explicit bound at the frozen amplitude, not just existence of
an unspecified compact supremum. K_act is distinct from the core K.

## 5. The existing width works, with scale-independent leading gaps

Set t1=kappa0=C^-20 as already prescribed. On the native core collar,
`|partial_y vr|<=C^2*R^8`, so (11) gives vr>=9/4 for 0<=y<=t1.
Also vr<=C^3 there. The same amplitude makes `t1*R^3<.1` and
`K_act*t1<C^-17`. No parameter is retuned after seeing a diagnostic.

Use the fixed smaller collar

\[
 s_c=1/(10\sqrt{\log C}),\qquad t_c=t_1s_c.                 \tag{13}
\]

It lies inside the native core and the activation interval. Since
`sigma(s)<=exp(16/9-1/s^2)` for s<=1/4,
`ea*(C^3+1)<=2*exp(16/9)*C^-97<1/8` on this collar. Equations (12)
give the strict bounds

\[
 v_s\ge17/8>2,\qquad (P_c-v_s)/e_a\ge9/8,
\]
\[
 \frac{(v_s-2)J_c^2}{2(P_c-v_s)^2}
       \le(32/81)C^{-31}s^6<1/2.                            \tag{14}
\]

Thus the pointwise pulse cone quotient is <=1/sqrt(2). Both shear and
growth have explicit positive profile bounds. With
`Fmin=exp(-2Lambda M)/(4C)`, one has

\[
 F\ge F_{min},\quad |g_0|\ge F_{min}/(5\Lambda),\quad
 \lambda_0^2=2aF^2(1-2/v_s)
             \ge e^{-4\Lambda M}/(680\Lambda C^2)>0.        \tag{15}
\]

The reference bound and the table give `|Bstress/F|<=2C^2` and
`sqrt(1+ts^2)<=2C^2` in T0=ea*Bstress. Hence orientation is at least
`-n dot N>=9/(32C^4)>0`, and
`|T0|>=9*Fmin*exp(-1/s^2)/(16C^2)`. The flat-quotient proof supplies
the smooth limiting direction and every fixed weighted derivative;
physical radial derivatives keep their t1^-j costs. Constants depend on
the frozen profile, not physical q, dyadic band or correction stage.
They may be astronomically small. Positive does not mean practical to sample.

This closes the local leading inner-collar growth/cone estimate for the
defined reference datum. It does not certify the summed background stress
without the later support/moment identities, and does not certify a pulse
on every part of the still unaudited compact annulus.

## 6. Checks and next task

The evidence includes exact rational scalar/separation inequalities,
high-precision eta/j diagnostics for a labelled atomic pressure fixture,
positive scalar series checks, and outward logarithmic budgets at 100 and
160 decimal digits. Exponentiating the nested parameter logarithms loses
about 33 relative digits at Md64; the 80/120-digit refinement trial only
agreed to about 47 digits. The stronger pair agrees past 55 digits and
both outward enclosures satisfy the required strict inequalities.
Interval gates use the full internal precision; printed endpoints are
rounded numerical diagnostics.
A separate diagnostic replays the complete forward moment
formulas of the activation checkpoint, keeping pressure and checking the
C-scaling. Manufactured measures are not substituted for the defining
schedule. The continuum arguments are sections 1--5, not finite samples.

Next: independently audit compact modulation and the exact final repair,
including pressure memory, then join their cone bounds to these inner and
the already checked outer estimates. Global coefficient/mean supports and
the unconditional total residual gate remain open.

```sh
python3 research/navier_stokes_cascade/results/core_growth_activation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_growth_activation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_growth_activation/evidence.json \
  --output core_growth_activation_recheck.json
```
