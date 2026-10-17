# The actual flat-activation collar: common parameter strip and every fixed jet

This checkpoint supplies the first-collar estimate (5) left open in
[core_all_order_foundation](../core_all_order_foundation/README.md), for the
ACTUALLY activated leading functions, relative to the same exact
pressure-family core. It propagates one common complex eta neighborhood,
retains all incoming moments and pressure, and records every fixed-width
radial derivative loss. The unmodified seed is not substituted for the
activated profile. No parameters, supports or width C^-20 are changed.

Parent: `b1c66663d020a3920aafc58cb46aae421aeda8c1`. Defining inputs are the
infinite core and explicit zero-free neighborhood in the parent, the actual
shear prescription in [axis_core_attachment](../axis_core_attachment/README.md)
section 4 and [inner_collar_activation](../inner_collar_activation/README.md)
sections 1--3, and the frozen collar length in
[core_growth_activation](../core_growth_activation/README.md). These are
defining maps and local estimates, not an assumed final manuscript theorem.
Historical source, evidence and tolerances remain unchanged.

The exact positive-measure pressure datum, core contraction hypotheses,
fixed finite N,h,Lambda,C and the prior analytical estimates remain explicit
inputs. The completed global backward-pressure identity is NOT independently
re-proved here. Neither this local certificate nor its tests prove the full
background, actual pulse bounds, complete residual improvement or blow-up.

## 1. Verify which defining fields occur on this collar

Write Xa=4/Lambda, y=log(X/Xa), t1=kappa0=C^-20,
tc=t1/(10 sqrt(log C)) and ell=tc/4. The required background collar is
0<=X<=b^2=Xa exp(ell). Since C>=2 and log 2>=1/2,

ell/t1<1/20, ell<1/100 and 4 exp(ell)<400/99<4.1.

Thus the native reference on 0<=y<=ell is still the exact analytic core:
its taper starts at y=t1. Later modulation starts at y=tc/2, strictly
after this collar. The exact cutoff ordering remains

tc/8 < 3tc/16 < tc/4 < tc/2

for keep, cut, collar end and modulation start. This verifies the domain
against the actual radii, not a rounded log Xa. The parent core alone
describes the finalized leading functions only on Y<=4; on the present
right-hand part Y=4 exp(y)>4 the following actual activation is essential.

Let F_r=g Phi and U_r be that core reference, E_r=sqrt(2X)F_r, and
a_r=1-2 partial_y log E_r=-2Y Phi_Y/Phi. Retain the actual flat cutoff

e_a(y)=(1-kappa0)sigma(y/t1), kappa=1-e_a,

where sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), extended by
zero to s<=0. It has NO eta dependence. Define the two forward integrals

\[
J(y,eta)=\tfrac12\int_0^y e_a(v)a_r(v,eta)\,dv,\qquad
K(y,eta)=\int_0^y e_a(v)\partial_vU_r(v,eta)\,dv.\tag{1}
\]

The ACTUAL fields are exactly

\[
F_a=F_r e^J,\qquad E_a=E_r e^J,\qquad U_a=U_r-K.\tag{2}
\]

Indeed partial_y log E_a=(1-kappa a_r)/2 and
partial_y U_a=kappa partial_y U_r, exactly the frozen shear prescription.
There is no division by kappa, kappa0 or the flat number e_a.

## 2. A common eta strip without a small-perturbation assumption

The parent supplies a zero-free seed neighborhood with eta radius
delta0=rho/(100S) for every real 0<=Y<=4.1. Choose the SAME
delta=delta0/8=rho/(800S) for this entire argument. Work successively on
the strips of radius 4delta, 2delta and delta around [-1,1]. All fit strictly
inside the parent's neighborhood. In particular a_r and every fixed radial
derivative are holomorphic in eta and bounded there. This uses the parent's
explicit reciprocal bounds, not just real positivity of Phi.

For fixed real y, differentiation and forward integration in (1) preserve
holomorphy. Uniform compact bounds permit differentiation under these
integrals, including at eta=+/-1. Exponentiation in (2) introduces no zeros:

|1/F_a| <= |1/F_r| exp(ell sup|a_r|/2).

This argument does not require the complex perturbation to be small, and
does not retune t1 after seeing a norm. The bound may be enormous but is
finite for the frozen seed. Radial differentiation changes constants, not
the eta strip, because the cutoff depends only on y. The same is true of
every fixed mixed derivative by the product rule and the exponential chain
rule. On a smaller strip Cauchy gives k! delta^-k times the corresponding
outer-strip bound for the kth parameter derivative.

Flatness of e_a at y=0 makes J,K and every fixed mixed derivative vanish
there. Equations (2) therefore join the unchanged core C-infinity in X,
with the SAME eta strip. This is not radial analyticity: at a real point
where a_r>0, J>0 for y>0, so F_a differs from F_r despite having all the
same radial jets at the join. A Taylor series there cannot replace the
actual flat activation.

## 3. Full incoming memory and the actual pressure increment

Carry the exact core traces of M,I,J_m,S,Cp at Xa. The symbol J_m here
denotes the moment, not the logarithmic field integral J in (1).
The five forward densities in y are

| Moment | Actual derivative |
|---|---|
| M | X U_a |
| I | 2 X^2 F_a |
| J_m | 2 X^2 U_a F_a |
| S | X(U_a^2-X F_a^2) |
| Cp | X F_a^2 |

For each row, integrate from its ORIGINAL core trace at y=0. Keep also
the fixed axis datum Pi_axis, so the full pressure is

\[
Pi_a=Pi_{axis}+Cp_a
     =Pi_r(Xa,eta)+\int_0^y X(v)F_a(v,eta)^2\,dv.\tag{3}
\]

Using F_r instead of F_a in (3), or setting the incoming trace to zero,
changes the pressure. A derivative equation alone cannot determine its
eta-dependent integration constant. Pressure value AND parameter
derivatives are retained. The diagnostics separately retain Cp and Pi;
Pi-Cp is the same nonzero radially constant axis datum.

Each density has every fixed radial jet uniformly holomorphic on the
outer strip, so the moment integrals have that property too. The core
traces are holomorphic by the parent's whole-core radial integrals.
On X>=Xa, division by X is harmless with constants depending on Xa.
On the axis, M/X keeps its actual average formula int_0^1 U(sX,eta)ds,
not an arbitrary zero assignment. Finally

V/X=[2eta U-2Deta(M/X)-(1-eta^2)partial_eta(M/X)]/L,

where D=1/2-h and L=1-2h eta^2. This needs one additional parameter
derivative. For the retained tiny strip |eta|<1.01 and h<=.01,
|L|>=1-.02(1.01)^2>.97. Nested-strip Cauchy bounds therefore control
every fixed jet of V/X. No new axis singularity or inverse-kappa factor
is introduced.

Equation (3) is the local forward pressure with the same prescribed datum.
Agreement with the globally completed backward pressure still requires the
global exact moment identities; locality alone does not prove it.

## 4. All orders with honest width costs

To expose the costs, first freeze the seed and allow an auxiliary
0<t<=.01, 0<=kappa0<=1/2 and 0<=y<=t/4. This auxiliary interval includes
the diagnostic s=.24 and is wider than the actual collar's y<t1/20;
the reference still stays inside Y<=4.1. Write C_rk for finite constants
depending on that seed, strip and requested orders, but not t or kappa0.
The flat step gives |partial_y^r e_a|<=c_r t^-r. Thus (1) gives

|partial_y^r partial_eta^k J|+|partial_y^r partial_eta^k K|
 <= C_rk t^(1-r).

For r=0 this is the length of the integral; for r>=1 differentiate its
integrand r-1 times. In the rth derivative of exp J, a Bell monomial
with n>=1 factors of derivatives of J has total width power t^(n-r),
which is no worse than t^(1-r). At r=0 use
|exp J-1|<=|J| exp|J|. Consequently

\[
|\partial_y^r\partial_\eta^k(F_a-F_r,U_a-U_r)|
 \le C_{rk}t^{1-r}.\tag{4}
\]

Polynomial density differences have the same gain. Their incoming traces
CANCEL only because the original traces are kept. One additional Volterra
integration then proves, for all five moments and the full pressure,

\[
|\partial_y^r\partial_\eta^k\Delta(M,I,J_m,S,Cp,Pi)|
 \le C_{rk}t^{2-r}.\tag{5}
\]

These estimates include r=0,1 and every higher fixed order; no unbounded
bare eta-derivative operator is assumed. Parameter derivatives can instead
be obtained by Cauchy on the unchanged strips after the radial estimates.
They are estimates for the differences, not assertions that incoming
moments or pressure themselves are small.

Convert to physical profile radius with the exact operator identity

\[
\partial_X^r=X^{-r}\prod_{j=0}^{r-1}(\partial_y-j).\tag{6}
\]

The signed first-kind Stirling coefficients, including all lower orders,
are retained. For example the fourth derivative contains
X^-4(D_y^4-6D_y^3+11D_y^2-6D_y), not just X^-4 D_y^4.
Since X>=Xa, (4)--(6), the regular quotients and the unchanged core give
on one common delta>0, for every r,k,

\[
\sup_{0\le X\le Xa e^{tc/4},\;dist(eta,[-1,1])<delta}
\sum_{H\in\{F_0,U_0,Pi_0,M_0/X,V_0/X\}}
|\partial_X^r\partial_\eta^k H|\le B_{rk}t_1^{-r}<\infty.\tag{7}
\]

This supplies the parent's common-collar obligation for the ACTUAL fields,
including both coefficient cutoff radii. For the actual frozen width,
t1^-r=C^(20r). Constants are not uniform as C varies, and may grow with
derivative order. Neither the finite-order K_act bound nor C2/C6 closeness
was substituted for (7).

The parent's exact Cartesian similarity operators now apply to the whole
actual first collar, including the axis, giving the same leading velocity
and pressure bounds uniform in positive physical q<=1, with weight
w=(alpha_x+alpha_y)/2+(1/2-h)alpha_z+j. Their constants include these fixed
width costs. This is not an unconditional estimate for the full corrected
field or the actual pulses.

## 5. Checks and the next foundation obligation

The new manufactured implementation assembles radial Taylor jets through
order four from the actual shear ODEs, independently integrates all moment
excesses on a partial interval, and compares values with the older
cancellation-safe five-moment implementation. Refinement and comparisons
are RELATIVE TO THE SMALL EXCESSES, not normalized by the order-one state.
Forty-two complex-parameter Cauchy checks cover eta orders zero through six.
Cp, Pi and their nonzero incoming difference are kept separately.

At widths .01,.001,.0001 the fourth field derivative grows by about a
million, whereas its predicted t^3-scaled size and the t^2-scaled fourth
pressure difference stay within a factor of two. These are finite
diagnostics, not proofs of uniformity. The continuum derivations are
(1)--(7). Their toy seed is NOT the exact huge core and has a nonzero
baseline stress; no cone or zero complete residual is inferred from it.

Next independently audit the later reference continuation/shape and
compact modulation/repair for the SAME actual leading profile, with full
global pressure/moment identification and weighted edges. Then join the
verified collar input to the background inverse and actual pulse bounds.
The global construction, complete residual gain and blow-up remain
conditional/unverified; this local collar result does not discharge them.

```sh
python3 research/navier_stokes_cascade/results/activated_collar_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/activated_collar_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/activated_collar_foundation/evidence.json \
  --output /tmp/activated_collar_foundation_fresh.json
```

The record binds 14 current source/evidence files and 771 inherited source
checks. Historical evidence, tolerances and source bytes are preserved.
