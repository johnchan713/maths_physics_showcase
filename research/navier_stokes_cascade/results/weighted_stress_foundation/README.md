# Actual weighted first stress and the background foundation join

Parent: `a779312a3036c03ca7abea540b18e19785c9bc03`.
This independently checks the flat-edge argument left open by the
[annular/exterior foundation](../annular_exterior_foundation/README.md),
then supplies its actual leading inputs to the existing formal coefficient
and weighted-background maps. The result is relative to the retained exact
pressure-family core and finite existence/inverse/gap lemmas. Those remain
explicit analytical dependencies; neither sampled approximate roots nor
passing tests replace them. No final manuscript theorem is assumed.

The new conclusion is a smooth all-order relative FIRST-stress quotient on
the CLOSED active annulus, followed by one admissible diagonal schedule
with normalized background error O(q^(2h)) in every fixed profile/Euler jet.
It retains the exact inner t1^2 and outer 4 exponential constants, full
pressure memory, the shifted order-one axial viscosity and the exceptional
first two pressure tails. Historical sources, evidence, supports, width,
h and finite N are unchanged. This is not a proof-assistant certificate,
an unconditional actual-pulse bound, complete residual closure or blow-up.

Source: (4.7)--(4.16), (4.26)--(4.27), (5.2)--(5.22), (5.42)--(5.46),
(6.23) and the explicit tail (A.47)--(A.56) of the retained
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The preceding [weighted estimate](../weighted_background_edges/README.md)
is an implication of then-unverified leading inputs, not itself evidence
that they hold. Sections 1--4 below connect the defining maps to those inputs.

## 1. Actual inner factor, not an absolute residual tolerance

Keep Xa=4/Lambda, t1=C^-20, y_a=log(X/Xa), s=y_a/t1, kappa0=t1 and
e_a=(1-kappa0)sigma(s), with the original flat step sigma. The exact core
reference remains valid through a fixed smaller first collar: its taper
starts at y_a=t1, and the coefficient solution region ends before tc/4<t1/20.
Later modulation starts at tc/2 and cannot alter this collar. The
[actual activation](../activated_collar_foundation/README.md) supplies all
fixed jets there on a common eta strip. The exact core equations give
T_r=0 as a FUNCTION, not merely at selected radii or angles.

For g(s)=[exp(-1/s^2)+exp(-1/(1-s)^2)]^-1, define the flat Volterra map

\[
 {\int_0^s\sigma(u)u^m b(u,eta)du\over\sigma(s)}
 =s^{m+3}{1\over2g(s)}\int_0^\infty e^{-v}
 (1+s^2v)^{-(m+3)/2}g(u)b(u,eta)dv,
 \quad u={s\over\sqrt{1+s^2v}}.                    \tag{1}
\]

Every fixed derivative has an exponential times polynomial majorant.
Thus (1), including its nonzero endpoint b(0,eta)/2, is smooth jointly.
The actual shear equations give Delta log E and Delta U divisible by
t1 e_a s^3; the five moment excesses and full Delta Pi are divisible by
t1^2 e_a s^6. Original incoming traces cancel in these DIFFERENCES only.
Inserting them into the original forward stress map retains
4A eta Delta Pi-d Delta Pi_eta and the local reciprocal-field variation.
The latter is why the stress-vector excess has only t1 e_a s^3, not s^6.
It gives exactly

\[
 T_0=e^{-1/s^2}A_a(s,eta),\quad
 A_a=(1-kappa0)g(s)F_a\{(a_r,q_r)+t1 s^3 B_a\},
 \quad A_a(0,eta)=(1-kappa0)e F_r(a_r,q_r).        \tag{2}
\]

All coefficients in (2) have every fixed joint jet. F_r and a_r have
strict positive lower bounds from the retained core-growth estimate.
Compactness in eta therefore gives ONE fixed smaller collar where
|A_a|>=c_a>0. Shrinking a proof collar does not change the activation width.
The regional value/cone estimates treat the remaining middle. Since
exp(-1/s^2)=exp(-t1^2/y_a^2), its constant is exactly t1^2.

If the seed has a defect, the identity is instead
T_0=(F_a/F_r)T_r+e_a B. A tolerance |T_r|<=epsilon cannot control
T_r/e_a near zero for ANY epsilon>0. Nor can replacing t1^2 by another
positive constant be absorbed in a polynomial loss. These are exact
obstructions; the numerical seed is not substituted for the exact core.

## 2. Supply the actual coefficient and backward-support inputs

The parent now supplies smooth F0,U0,Pi0,M0/X,V0/X on every needed closed
rectangle, a common analytic FIRST collar, the actual positive I3 power
patch, and U0=M0=V0=0 beyond the completed axial edge Xv. These are precisely
the inputs, previously separated as assumptions, of
[background-jet foundation](../background_jet_foundation/README.md).
On the common inner collar the sparse derivative blocks satisfy
A1(xi)D0 A1(s)=A1(xi)D0 partial_eta A1(s)=0 even at different radii.
The retained factorial Picard bound, for each finite coefficient order,
therefore supplies the same-radius smooth inner coefficient. A finite
intersection of parameter strips suffices at each order; global heat
analyticity, uniform coefficient norms and a common strip for ALL orders
are not required.

Cut only the velocity fields at the prescribed keep/cut radii, reconstruct
M,V and pressure by their forward equations, then solve the actual five
linear I3 moment equations. Their fixed ordered positive bumps and
positive e3 give the retained normalized inverse <1000; the physical
lambda^-1 and e3^-1 losses remain. Current-order squares have order 2n,
not n. All incoming lower-order, pressure and viscosity terms are kept.
This supplies the smooth coefficient induction, not a small-target argument.
The third zero total removes pressure constants, not merely its derivative.

For clarity, the conservative integrals giving the backward stresses are
part of the join. At every positive order n the integrated time and
nonlinear angular terms use respectively the second and fourth moments;
the axial flux uses the first and fifth, including pressure by integration
by parts. Radial divergence boundaries vanish by axis regularity and the
exterior moment equations. Axial viscosity uses order n-1, not n. At n=1
its angular integral is the convergent LEADING swirl minus the stationary
power tail; the parent's all-order Gamma majorants justify differentiating
it. Its total is zero by R=0. The leading axial-flux integral is zero by
M_infinity=0. At later orders the corresponding previous-order compact
moments vanish. Thus BOTH weighted residual totals vanish at every order.
The negative forward stress integrals consequently equal positive backward
integrals, with no unknown integration constants.

The inner coefficient equations imply T_n=0 through
Xinner=Xa exp(tc/16), for every n>=1. With Xvel=rho3 exp(3), all positive
velocity coefficients and M_n,V_n vanish after Xvel. The exact support
and shifted radial-pressure source then give:

| Order | Pressure right edge | Stress right edge |
|---|---|---|
| n=1 | Xv | Xb |
| n=2 | Xv | Xv |
| n>=3 | Xvel | Xvel |

In particular n>=2 stresses have a COMMON compact interior support
K=[Xinner,Xv] inside (Xa,Xb). T1 alone reaches the outer flat edge.
This uses the coefficient equations and exact total moments, not U0=0.
On I4 the leading M0=mu eta,V0=-mu remain, and Pi1=K1-mu^2/(4X),
Pi2=K2. Neither is erased; potential cutoffs add the extra M_n term to
radial velocity, which vanishes on I4 because M_n does, not because U_n does.

## 3. Independently remove the actual outer flat factor

Set delta=log(Xb/X), r_delta=sqrt(2Xb)exp(-delta/2),
K_delta=c_infinity Xb^(-A)exp(A delta)H(2d exp(delta)/Xb),
A=1/2+h, D=1/2-h, L=1-2h eta^2, d=1-eta^2, rho=h/1000.
The exact tail multiplier has f=1-rho sigma(delta/2), so

\[
 f_y(3-u)=rho e^{-4/u^2}u^{-3}B(u),\quad
 B(u)={e^{-1/(1-u/2)^2}\over
 [e^{-4/u^2}+e^{-1/(1-u/2)^2}]^2}
 \left(8+{u^3\over(1-u/2)^3}\right),\quad B(0)=8e.  \tag{3}
\]

Write tH=-ZH'/H, Hrad=r^3 K, Bheat=2+2h-2tH. Physical axial
differentiation at q=1 gives alpha=4eta^2/L^2 and
beta=-2(d-2Deta^2)/L^2-8h eta^2 d/L^3. They are fixed during radial
integration. Retain BOTH derivatives in -K(alpha f_yy+beta f_y).
Integration by parts, with f_y=0 at the edge, gives

\[
 T_{1theta}={alpha rK\over2}f_y+
 {1\over2r^2}\int_y^3 Hrad(v)C(v)f_y(v)dv,\quad T_{1z}=0,
\]
\[
 C=alpha(1-h+tH)-beta
 ={2\over L^2}+{8h eta^2d\over L^3}+alpha tH>0.       \tag{4}
\]

The leading backward formulas are
T0theta=2K f_y/r+r^-2 int rK Bheat f_y+
(2Lr^2)^-1 int Hrad f_y and
T0z=eta/(Lr) int (r_v^2-r^2)K_v^2 f f_y.
Section 2 and the parent's exact global ledger now identify them with the
ACTUAL forward stress functions, including parameter derivatives.

For any smooth joint coefficient J the exact transform is

\[
 \int_0^delta e^{-c/u^2}u^{-3}J(u,eta)du
 ={e^{-c/delta^2}\over2c}\mathcal L_c J(delta,eta),
\quad \mathcal L_cJ=\int_0^\infty e^{-w}
 J\left({delta\over\sqrt{1+delta^2w/c}},eta\right)dw. \tag{5}
\]

All fixed mixed derivatives of (5) have polynomial-exponential majorants.
The parent supplies every heat jet at eta=+-1 by one-sided limits; no
heat integral is evaluated at negative Z and no analytic heat strip is used.
For normalized radial Taylor jets J=sum j_k delta^k, the exact endpoint
coefficients v_k of L_c J satisfy TWO equivalent finite formulas:

\[
 v_k=\sum_{n+2a=k}j_n { {-n/2}\choose a}{a!\over c^a},
 \qquad v_k=j_k-{k-2\over2c}v_{k-2}.               \tag{6}
\]

The second follows from delta^3 V'+2cV=2cJ; missing its radial Jacobian
would give V=J, which is wrong already at degree three. These are finite
jets with Taylor remainder, NOT claims of convergent infinite series.

For the axial factor use, with u=delta/sqrt(1+delta^2w/4),

\[
 {delta-u\over delta^3}={w/4\over
 \sqrt{1+delta^2w/4}(1+\sqrt{1+delta^2w/4})}.
\]

Also (r_u^2-r_delta^2)/delta^3 equals r_delta^2 times this expression
times exprel(delta-u). This removes the apparent division by delta^3
and has all fixed integrable derivatives. Equations (3)--(6) then give

\[
 T_0=e^{-4/delta^2}delta^{-3}(b_0,delta^6 b_z),\quad
 T_1=e^{-4/delta^2}delta^{-3}(b_1,0),                \tag{7}
\]

with jointly smooth coefficients, b0 bounded strictly below on a smaller
collar, and b1>=0. At delta=0, putting rb=sqrt(2Xb), Kb=K(0,eta),

\[
 b_0=16rho e Kb/rb,\quad b_z=eta rho e rb Kb^2/(8L),
 \quad b_1=4alpha rho e rb Kb.                     \tag{8}
\]

There is no division by eta or b1. Hence

\[
 W_1={T_1\over|T_0|}=
 \left({b_1\over\sqrt{b_0^2+delta^{12}b_z^2}},0\right),
 \quad W_1(0,eta)=(2eta^2Xb/L^2,0).               \tag{9}
\]

This is smooth even at eta=0. A new sharp diagnostic consequence of the
actual formulas is

\[
 W_{1theta}(delta,0)={Xb\over8}delta^3+O(delta^4),
 \quad partial_delta^3 W_{1theta}(0,0)=3Xb/4,
 \quad partial_eta^2 W_{1theta}(0,0)=4Xb.           \tag{10}
\]

Indeed alpha=0,beta=-2,C=2 at eta=0; the integral in (4) has normalized
limit rho e rb Kb, while b0 has limit 16rho e Kb/rb.
The first two delta jets vanish, but the third does NOT. Keeping only
alpha f_yy gives zero and fails (10). The cubic coefficient is independent
of heat amplitude and h; Xb is fixed and may be huge, not a small constant.

## 4. Close the weighted and relative estimates

Use EXACTLY zeta=exp(-t1^2/y_a^2-4/delta^2) and
delta_edge=min(1,y_a,delta). For E_c(y)=exp(-c/y^2),
D_y^j[E_c y^p]=E_c y^p P_j(1/y), with
P_(j+1)=(2cz^3+pz)P_j-z^2 P'_j and degree <=3j.
The other exponential has a smooth positive reciprocal on each opposite
collar. Equations (2),(7), the exact X/log-X Stirling conversion and compact
middle positivity consequently prove, for total profile order |I|<=m,

\[
 |T_0|>=c zeta,\quad
 |partial^I T_0|+|partial^I T_1|
 <=C_I zeta delta_edge^{-3m-3}.                   \tag{11}
\]

W1 is zero on the fixed inner collar by the actual coefficient equations,
has the smooth expression (9) at the outer edge, and is smooth in the
middle where |T0|>0. For n>=2, both Wn=Tn/|T0| and Tn/zeta are smooth
with the SAME interior support K. Their every fixed seminorm is finite;
no uniform-in-n bound or uncut series convergence is asserted.

One admissible choice of the still-free Borel scales now includes these
finite seminorms through n alongside ALL earlier physical/potential and
residual requirements. This specifies a common choice, not an alteration
of any numerically frozen scales. If M_n,m bounds the two quotient jets
through m, and K_m bounds Euler cutoff jets, take
C_n,m=max(1,2K_m(1+2nh)^m M_n,m). Choose log c_n greater than
1+max(log(2c_(n-1)),log(2/q0),max_(m<=n)(log C_n,m+n log2)/(nh),
the earlier finite requirements). For sigma=c_n q, Euler derivatives of
chi depend on sigma, not c_n. On its WHOLE support q<=1/c_n,

\[
 |D^I[chi(c_nq)q^{2nh}Wn]|,
 |D^I[chi(c_nq)q^{2nh}Tn/zeta]|<=2^{-n}q^{nh}
 \quad(|I|<=n,n>=2).                              \tag{12}
\]

D=(q partial_q at fixed X,eta,partial_X,partial_eta). For any fixed m,
split at J=max(m,2). The finite block is O_m(q^(2h)); the remaining
geometric sum is <=2^-J q^((J+1)h)<=2^-J q^(2h) for q<=1.
This proves on ONE q-domain, for every fixed jet,

\[
 D^I[(That-T0)/|T0|]=O_I(q^{2h}),\quad
 |D^I(That-T0)|<=C_I q^{2h}zeta delta_edge^{-3m-3}. \tag{13}
\]

That is the direct cutoff STRESS series. Velocity uses potential cutoffs,
including their extra radial M term. Its augmented residual is a separate
quantity, not declared zero by summing stress. Euler/frozen derivatives
retain Bell/factorial conversion constants.

The smooth leading direction n0 extends by (2),(7). Choose qstar>0 once
using only the value and first tangential-velocity jets and (13): make the
combined perturbation C qstar^(2h) less than the minimum of half a fixed
denominator-neighborhood radius, half the leading growth/oriented-cone
gap divided by its local Lipschitz bound, and 1/2. All these constants are
positive and finite on the closed annulus. Then
|That| is between |T0|/2 and 3|T0|/2, its direction has all fixed errors
O_I(q^(2h)), and growth and the oriented cone remain strict. Higher jets
change constants, not this domain. The radial velocity in q^A normalization
is still O(q^h), NOT O(q^(2h)) or zero. Pressure on I4 retains both
exceptional orders; no pressure equality follows from velocity support.

A generic weighted perturbation cannot establish (13): adding
-q^(2h)e^(-4/delta^2)delta^-4 to a positive leading
e^(-4/delta^2)delta^-3 reverses its sign when delta<q^(2h).
Exact common supports and (9), not an absolute tolerance, remove this
obstruction for the specified stress coefficients.

## 5. Checks, limits and next obligation

Independent exact rational kernel/ODE jets test (6) through order 15,
including mixed angular coefficients; actual finite Gamma POLE jets test
the factored outer algebra without presuming a convergent Taylor series.
Manufactured power-tail jets check (10), and independent actual Gamma
tail quadrature checks its eta=0 cubic limit at decreasing delta. A dropped
beta term, dropped radial Jacobian and wrong flat constants are negative
controls. Exact all-smaller-q budget inequalities test the single schedule.
None numerically evaluates the actual gigantic h,N,Xb or Borel constants.

The next obligation is to apply these now joined background inputs to
uniform ACTUAL physical-pulse values, derivatives and supports at EVERY
scale, including carrier rounding and edge/axis charts. Then re-audit
complete-wave covariance means with axial velocity and Pi1/Pi2, and the
complete residual with viscosity, all cutoffs and interactions. Earlier
conditional cycle/realization/comparison results are not promoted to
unconditional results by this foundation checkpoint.

```sh
python3 research/navier_stokes_cascade/results/weighted_stress_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/weighted_stress_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/weighted_stress_foundation/evidence.json \
  --output /tmp/weighted_stress_foundation_fresh.json
```
