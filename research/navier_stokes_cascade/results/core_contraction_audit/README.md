# Nonlinear core contraction and exact cancellation of the leading stress

Status: `core-contraction-and-exact-stress-audited-for-pressure-family`.

This batch checks the infinite coefficient-space argument behind the inner
core, rather than only its final parameter arithmetic. It derives the
all-index product and integrated-derivative bounds, bounds every nonlinear
source including pressure, and recovers the exact original equations from
the Banach fixed point. For the stated exact analytic pressure family, the
selected parameter formulas give a contraction factor below `1e-50` and
both leading stress components vanish identically on the core.

The exact exterior pressure-family representation is an input to this
local theorem. Its identification with the finalized global construction
is not independently reverified here. Positive shear, the endpoint growth
alternative, the quantitative activation constant K_act of the previous batch,
compact modulation/repair and global moments remain separate obligations.
No unconditional actual-pulse estimate, complete residual gain, or blow-up
result is asserted. This is an analytical audit with reproducibility checks,
not a formal proof-assistant certificate or an independent peer review.

Parent: `5ac315d13e54797b3a8c692a3d465e41d55642ea`.
The source under audit is the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.9)--(4.13), (B.1)--(B.16), with SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its conclusion is not a premise. Earlier checkpoint bytes and GOAL.md
remain frozen. There is no new choice of N,h or a sampled pressure fit.

## 1. Exact pressure data and parameter order

Fix P>=16, 0<h<=.01, j=1e-18, sigma=1e-20 and an exact positive measure
mu on [0,1] with mass less than five. Assume the supplied datum is

\[
 \Pi_0(\eta)=-P^2\int_0^1(1+\eta^2)^{-2\theta}\,d\mu(\theta). \tag{1}
\]

The contraction below only uses this upper mass bound and analyticity;
the lower pressure contribution used to prove endpoint growth is not
inferred from it. Set A=.5+h, D=.5-h, d=1-eta^2, L=1-2h*eta^2 and

\[
 U_*=4\eta+j,\quad H_*=D\eta+dU_*,\quad
 W_*=1-dU_{*\eta}-2D\eta U_*,
\]
\[
 Z_*=-A(1-2\eta U_*)U_*-H_*U_{*\eta}-d\Pi_{0\eta}+4A\eta\Pi_0.
\]

Use r0=sigma/256, rho=r0/8, M=200/sigma^2 and Zb=1024*P^2/r0.
On `|Re eta|<1+r0`, `|Im eta|<r0`, one has |eta|<2, |d|<=5,
|L|>.9, |U*|<9, |H*|<64, |W*|<64. Since |H*'|<53, its vertical
variation is less than sigma/4. H* is real on the real axis, so both
`|H* +/- i*sigma|>=3*sigma/4`. Therefore

\[
 \chi=H_*^2/(H_*^2+\sigma^2),\quad
 \zeta_*=-LH_* /(H_*^2+\sigma^2),\quad |\chi|<3,\ |\zeta_*|\le M.
\]

The real part of 1+eta^2 exceeds 1-r0^2. A single logarithm branch in
(1) gives |Pi0|<=10P^2. Cauchy's inequality on the half-size rectangle
gives |Pi0_eta|<=20P^2/r0 and |Z*/L|<=Zb. For example the latter is
bounded by `1178+200*P^2/r0+160*P^2`, strictly below Zb.
These bounds apply to the whole measure family, not just three sampled P.

First fix these data, then choose Lambda below, and finally C. Define
`phi*=exp(Lambda*integral_0^eta zeta*)` on the simply connected rectangle,
and g=phi*/C. The integration path has length below two, so
`log C>=Lambda*(2M+1)` makes |g|<=exp(-Lambda)<=1 on the complex
rectangle. A bound only on the real axis would be insufficient for the
coefficient norm. This order does not assume C before Lambda is chosen.

The frozen amplitude formula also meets this condition: its
`B=1e20*(1+Lambda+K)^64` and `log C=10000*B+100*log(1/epsilon)`
dominate Lambda*(2M+1), since the K below exceeds M+1 and Lambda>2M+1.
No amplitude or enormous array of coefficients is numerically materialized.

## 2. Completeness and all-index operator bounds

For F=sum F_n(eta)Y^n, define

\[
 a_{nk}=\frac{20^{-n}\rho^{-k}k!\binom{n+k}{k}}{(n+1)^2(k+1)^2},\quad
 \|F\|_\rho=\sup_{n,k,\eta\in[-1,1]}|\partial_\eta^kF_n|/a_{nk}. \tag{2}
\]

A norm-Cauchy sequence converges uniformly at each coefficient and each
derivative order. The fundamental theorem of calculus identifies the
successive limits as derivatives of one coefficient function. Taking the
supremum in (2) then proves convergence in norm. Thus this space is complete;
it is not a finite Taylor space or a norm on field values alone.

For a product, the derivative binomial cancels the k! factors. The remaining
ratio is at most one, by counting subsets in two disjoint blocks:

\[
 \binom{i+j}{j}\binom{n-i+k-j}{k-j}\le\binom{n+k}{k}.         \tag{3}
\]

Splitting each squared-denominator convolution at its midpoint gives
`sum_{i=0}^n (n+1)^2/((i+1)^2*(n-i+1)^2)<=8*sum_{i>=1} i^-2`.
Apply it separately in n,k. The product constant is below 256; write A0=256.
An analytic degree-zero coefficient bounded by b on the half rectangle
has norm at most b, because `(k+1)^2/4^k<=1` after Cauchy differentiation.
This is why the fixed coefficient suprema of section 1 control every jet.

Write `D_Y=Y*partial_Y`, average F=sum F_n*Y^n/(n+1), and let Jnu
solve `Y*G_YY+nu*G_Y=F`, G(0)=0, for nu=1,2. Then

\[
 (J_\nu F)_{n+1}=F_n/[(n+1)(n+\nu)],
\]
\[
 a_{nk}/a_{n+1,k}=\frac{20(n+2)^2}{(n+1)(n+k+1)}\le80,
\]
\[
 a_{i,k+1}/a_{i+1,k}
 =\frac{20(i+2)^2(k+1)^2}{\rho(i+1)(k+2)^2}
 \le(80/\rho)(i+1).                                        \tag{4}
\]

Averaging has norm <=1; multiplication by Y, ordinary radial integration,
and Jnu have norm <=80. Parameter differentiation followed by radial
integration has norm <=80/rho. Bare parameter differentiation is not
bounded on this space and is never treated as such.

For `Jnu[(partial_eta F)*(D_Y G)]`, let the input radial degrees be i,j
and assign the one extra output degree to the differentiated F factor.
Equation (4) leaves the fraction

\[
 \frac{(i+1)j}{(i+j+1)(i+j+\nu)}\le1.                       \tag{5}
\]

The product convolution now uses degrees i+1,j, so (3) and the two sums
still apply. If F is averaged, its divisor i+1 improves (5). Omitting
either derivative also improves the bound. Additional undifferentiated
factors cost their usual algebra constants. A fixed analytic coefficient
is included as one such factor, with all its derivatives retained.
Consequently a J-integrated monomial with m unknown factors and a fixed
coefficient of norm b has bound `D0*A0^m*b` times their norms, where
`D0=80/rho`. This proof includes i=0, j=0 and every parameter order;
finite index diagnostics are only regression checks of these formulas.

## 3. The order-one inverse and every nonlinear term

Let `T=J2*(multiplication by chi)/2`. On sequences starting at radial
degree b, J2 has norm <=80/((b+1)(b+2)) and raises that degree by one.
Repeated application therefore gives

\[
 \|T^k\|\le\frac{30720^k}{k!(k+1)!},\quad
 \sum_{k\ge0}\|T^k\|\le e^{30720}=:V.                       \tag{6}
\]

Multiplication by chi preserves the lower radial degree, including after
parameter differentiation. The absolutely convergent series sum(-T)^k
is a two-sided inverse of 1+T: multiply its partial sums and let the
factorial remainder go to zero. Smallness of chi is not assumed.

Use Y=Lambda*X and
`phi=phi* Phi`, `U=U*+Lambda^-1 u`, `Pi=Pi0+Lambda^-1 p`,
`p=I(g^2*Phi^2)`. Here I is integration in Y. Put
`B=-2D*eta*average(u)-d*partial_eta average(u)`,
`W=W*+Lambda^-1 B`, `Hc=H*+Lambda^-1*d*u`.
Direct substitution into the original equations yields

\[
 2(Y\Phi_{YY}+2\Phi_Y)=-\chi\Phi+\Lambda^{-1}R_1,\quad
 2(Yu_{YY}+u_Y)=-Z_*/L+\Lambda^{-1}R_2,                     \tag{7}
\]
\[
 LR_1=[W+h(1-2\eta U)+du\zeta_*]\Phi+W D_Y\Phi+H_c\Phi_\eta,
\]
\[
 LR_2=[A(1-4\eta U_*)+dU_{*\eta}]u-2A\eta\Lambda^{-1}u^2
 +W D_Yu+H_*u_\eta+d\Lambda^{-1}uu_\eta
 -4A\eta p+dp_\eta-2\eta D_Yp.                             \tag{8}
\]

The mixed derivative products in (8) and R1 have the form controlled by
(5). All three pressure terms stay in the map. In particular the estimate
is on `Jnu*(b*p_eta)`; it does not assign a bounded norm to bare p_eta.
The explicit term sums, including 1/L<=2, are

| Source group | Sum of coefficient bounds | Used upper bound |
|---|---:|---:|
| R1 linear | 384.74 | 386 |
| R1 quadratic excluding zeta | 38.08 | 40 |
| R1 zeta quadratic | 10M | 10M |
| R2 linear | 442 | 442 |
| R2 quadratic excluding pressure | 32 | 32 |
| R2 pressure p,p_eta,Dp | 16+10+8 | 34 |

For example the R1 linear sum is `2*(64+.01*37)+128+128`.
The pressure norm is `||p||<=80*A0^3*S^2` when ||g||<=1 and
||Phi||<=S. A difference of its two Phi factors costs
`160*A0^3*S*||Delta Phi||`; no pressure amplitude is set to zero.

Take the center `Phi0=(1+T)^-1 1=f0(Y*chi)`,
`u0=-Y*Z*/(2L)`, and S=2+V+40Zb. Its closed unit ball in the maximum
pair norm has both component norms below S. Equations (4)--(8) give

\[
 M_1=D_0[386A_0S+(40+10M)A_0^2S^2],
\]
\[
 M_2=D_0[442A_0S+(32A_0^2+2720A_0^4)S^2].                   \tag{9}
\]

The 2720 term is `34*80`, with the additional algebra constant for its
fixed coefficient. Replacing one unknown factor at a time gives the
Lipschitz bounds M1'(S),M2'(S) in the maximum pair norm. Let
`Fb=V*(M1+M2)`, `Lb=V*(M1'+M2')`. The actual map is

\[
 (\Phi,u)\longmapsto
 (\Phi_0+(1+T)^{-1}J_2R_1/(2\Lambda),\ u_0+J_1R_2/(2\Lambda)). \tag{10}
\]

Its displacement from the center is <=Fb/Lambda and its Lipschitz
constant is <=Lb/Lambda. Discarded factors of 1/2 increase these bounds.

## 4. Uniform contraction, common analyticity, and the exact limit

On |Y|<=5, summing (2) and using the binomial generating function gives

\[
 |\partial_\eta^kF|\le\|F\|_\rho
 \frac{\rho^{-k}k!}{(k+1)^2}(1-5/20)^{-k-1}.                \tag{11}
\]

Cauchy's radial estimate on the disk of radius .9 gives the additional
factor `r!/.9^r` on [0,4.1]. For r<=2,k<=6 these constants are below
`Eb=10000*rho^-6`. The coefficient bounds also give a common eta radius:
any radius below `rho*(1-5/20)` works; use `3rho/8` with a safety margin.
Fixed radial derivatives have this same parameter neighborhood. This
proves common analyticity, not merely separate finite jet bounds.

Use the same frozen parameter formulas, with epsilon=1e-16:

\[
 K=\{10^6E_b(1+S+F_b+L_b+M+Z_b+P^2)\}^{16},\quad
 \Lambda=10^{12}K^8\epsilon^{-2}.                            \tag{12}
\]

There is a direct bound for the whole admissible family; it does not
depend on a scan of P values. Since K>1 and
`Eb*max(S,Fb,Lb)<=K^(1/16)/1e6`,

\[
 \frac{E_b\max(S,F_b,L_b)}\Lambda
 \le\frac{\epsilon^2}{10^{18}K^{127/16}}<10^{-50}.            \tag{13}
\]

Thus (10) maps the unit ball into itself and contracts by less than
1e-50. Completeness gives a unique fixed point in that ball. Iteration
from the center has tail at step m bounded by
`(Fb/Lambda)*(Lb/Lambda)^m/(1-Lb/Lambda)`. The limit has norm distance
<=Fb/Lambda from the center; it is not defined by the last coefficient
of a finite jet. Applying (11) proves all fixed derivative conclusions.

For real 0<=Y<=4.1, 0<=chi<=1 and
`f0(z)=sum_{n>=0}(-z/2)^n/(n!*(n+1)!)` exceeds .265 by its decreasing
alternating cubic lower bound. Equation (13) gives Phi>.264. Also
`||U(4/Lambda)-U*||C1<=Eb*S/Lambda<1e-50<epsilon/100`.
The solution is real by conjugation and uniqueness. The positive real
phi*,C and Phi make F=g*Phi positive; its lower bound is allowed to
depend on the already frozen finite parameters. No uniform lower bound
as C tends to infinity is claimed.

If a reciprocal of Phi is needed, shrink the complex neighborhood around
the real rectangle: compact positivity excludes zeros there. All operations
in (10) agree with the analytic equations on smaller common disks by
absolute convergence. Applying the radial differential
operators to the fixed point therefore gives (7) exactly. Undo the
normalization in (8):

\[
 -2L(X\phi_{XX}+2\phi_X)/\phi=S_q,\quad
 -2L(XU_{XX}+U_X)=S_n,\quad \Pi_X=\phi^2/C^2.                \tag{14}
\]

The last identity follows from `Pi=Pi0+Lambda^-1*p` and
`p_Y=g^2*Phi^2`. The holomorphic amplitude bound is essential here;
pressure and its eta and radial derivatives have not been removed.

## 5. Exact forward stress cancellation and axis regularity

Define candidates `Q=-2L*phi_X/phi`, `N=-2L*U_X`. They are regular at
X=0. Direct differentiation of (14), with `l=1+X*phi_X/phi`, gives

\[
 D_X Q+(1+l)Q=S_q,\qquad D_X N+N=S_n.                       \tag{15}
\]

These are the defining radial source equations. Their homogeneous
solutions are a parameter constant divided by XH and by X. Regularity
at the axis forces both constants to zero, so these candidates equal
the actual forward primitives Qs,Ns, with the original pressure datum.
There is no integration constant to estimate or fit numerically. Hence

\[
 p_{s,1}=XQ_s/L=-2X\phi_X/\phi=a,\quad
 p_{s,2}=XN_s/(LE)=-2XU_X/E=-b_s,
\]
\[
                         T_0=F[p_s-(a,-b_s)]\equiv0.       \tag{16}
\]

This removes the baseline-defect term of the previous activation theorem
for this exact core, relative to its stated pressure input. The infinitely
small activation factor is never used to divide an approximate core error.

F=phi/C, U and Pi are analytic in X,eta through X=0. Radial averaging
preserves that property, so the divergence-compatible expression for V0/X
has its regular axis value. In Cartesian coordinates both radial and swirl
velocities are r times analytic functions of r^2; the axial velocity is
analytic in r^2. This supplies local axis regularity, not a global solution
with controlled forcing. Later activation leaves X<=Xa unchanged.

## 6. Reproducibility, failure controls, and what remains

`core.py` supplies exact rational weight ratios and a padded bivariate
Taylor algebra. Finite coefficient iteration is triangular: radial
inversion raises degree. Its original angular/axial residuals, pressure
identity and forward stress primitives are computed separately. They
vanish through the recorded interior jet orders when pressure is retained.
Deleting pressure, only p_eta, or only Dp produces a nonzero original
residual. These finite jets use toy finite parameters; they do not evaluate
the giant schedule or prove convergence.

The numerical reviews also compare independently computed outward core
budgets with the frozen parameter certificate, test both sides of the
angular inverse, and retain nonzero endpoint integration constants as
failure controls. Finite index samples accompany the proof in section 2;
they are not an exhaustive replacement for its quantifiers. Likewise
the three numerical P scales do not replace the universal bound (13).

The new conclusion supplies the exact-core existence/cancellation step
FOR THE STATED PRESSURE FAMILY. It does not supply the positive gap
`ar>0`, the complementary angular/axial endpoint growth bound, or the
actual normalized activation constant K_act from the previous batch. In
particular, one still needs strictly positive actual minima of ar and
`ar+qr^2/ar-2` at Y=4, with the collar derivative bounds defining K_act.
Checking those bounds at the same fixed parameters is next. The identification of
the exact globally assembled pressure, compact cone repairs, mean-bump
supports and total residual improvement still need independent audits.

```sh
python3 research/navier_stokes_cascade/results/core_contraction_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_contraction_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_contraction_audit/evidence.json \
  --output core_contraction_audit_recheck.json
```
