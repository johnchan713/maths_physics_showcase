# Q30: bind the cutoff and slow-response support to the fixed choices

Parent: `cee654124843aeb2ce26234df8dae44c04366f3d` (Q29).
This checkpoint discharges the SUPPORT quantifier left in Q29 for the
current construction as defined by its fixed profiles. It also makes the
finite matrix denominators of the slow response explicit. It does not
verify an independently chosen historical implementation or evaluate the
selected huge profile seminorms. No cutoff, bump, root, h or band is changed.

## 1. The physical cutoff supplies one gap, not one gap per band

`compact_mean_corrections` (4) defines a fixed cutoff, zero near the inner
shell edge and one near the outer edge. Its physical gluing statement
specifies chi_m(r^2/(2q)); `compact_mean_foundation` section 1 retains that
same X-dependent cutoff under all chart changes. Thus write the EXISTING
choice as chi_hat(X). The shell endpoints Xa,Xb are fixed positive numbers.
The plateau definition supplies fixed c0,c1 satisfying

    Xa<c0<c1<Xb,
    chi_hat=0 on (Xa,c0), chi_hat=1 on (c1,Xb).         (1)

If necessary move the witnesses c0,c1 toward their respective endpoints
so they are ordered; this changes a proof interval, not the cutoff.
For xlog=log(X/Xa), Lx=log(Xb/Xa), set

    eta=min(log(c0/Xa),log(Xb/c1))>0.                  (2)

Both summands are positive and their sum is less than Lx, so eta<Lx/2.
Every positive-order derivative of the actual moving chart cutoff
chi_hat(R^2/(2s)) is supported where X is outside the two plateaus.
Consequently Q29's uniform gap condition holds with (2) in every band,
for all s in [1/4,4]. Differentiating a moving plateau creates no nonzero
jet at a point in its open zero/one region. Boundary jets follow by
smoothness. Q29's proof collar is min(1,Lx/4,eta).

An explicit outward lower bound, avoiding numerical logarithms, is

    eta >= min(2(c0-Xa)/(c0+Xa),2(Xb-c1)/(Xb+c1)).      (3)

This uses log r>=2(r-1)/(r+1) for r>=1. Numerical values of the witness
intervals are not required for the existence conclusion, but any proposed
numerical operator constant must use valid witnesses for the SAME cutoff.
A new per-band cutoff family would not be covered. Its uniform gap would
remain a separate condition even if each member had some plateau.

The fixed shape has finite derivatives on [c0,c1]. On the active positive
shell, R stays in one fixed positive integration hull and s stays in
[1/4,4]. The retained scale-coordinate jet bounds hold for every
0<h<1/100, including the exact root h. Finite real chain rules therefore
supply all cutoff and multiplier seminorms used in Q29, uniformly in the
running band. No arbitrary-profile uniformity or derivative-order
uniformity is asserted. Q3's earlier compactness argument for an infinite
family of ODE propagators is not used here: this is one fixed shape composed
with already bounded coordinate maps.

## 2. One convex hull covers all fixed slow profiles and their potential

Use x=r/sqrt(q), distinct from xlog. The fixed nonnegative base bump eta0
has support in [l0,u0], 0<l0<u0, and positive mass. The actual geometric
bumps are eta_j(x)=g^-j eta0(x/g^j), g=exp(dpatch)>1.
The angular block uses j=0,1,2; the axial block uses j=0,1. Hence the
common closed potential hull is

    K=[l0,g^2 u0].                                    (4)

The original five-map chooses these fixed bumps with their WHOLE convex
hull strictly inside the x-image of I4. This is stated explicitly in
`compact_mean_foundation` sections 2--3 and `core_slow_inverse_budget` (1).
The current same-profile I4 has 6<log(X/rho3)<11, by the coefficient-support
recheck. Since X=x^2/2 is increasing, (4) lies strictly inside the same
positive shell. This uses the already fixed choice, not permission to move
bumps in an independently specified historical field.

For Fg, the zero first moment int x Fg=0 implies

    Psi(x)=x^-1 int_left^x y Fg(y)dy

vanishes outside the axial convex hull. Inside gaps it need not vanish.
All finite derivatives have the same closed support and smooth zero
extension. The angular profiles FP,Fz and every slow velocity/potential
therefore have support in K. In a Q chart that support is sqrt(s)K, and
R^2/(2s)=x^2/2. The reserved-geometry provider's weight floor applies on
that actual moving hull, with one fixed positive distance to each edge.
The rectangular integration enclosure [l0/2,2g^2u0] is NOT its weight-floor
set. No extra q loss is introduced by this support statement.

## 3. All fixed response coefficients are finite at the actual lambda

Put v=exp(-2dpatch lambda), 0<v<1, and mu_p=int x^p eta0(x)dx>0.
For the actual fixed lambda>0, all these moments are finite and positive:
mu_p is between the bump mass times the minimum and maximum of x^p on
[l0,u0]. The normalized angular Vandermonde nodes are

    z1=g^2, z2=v/g^2, z3=v,
    |det Vtheta|=(g^2-v/g^2)(g^2-v)(v-v/g^2)>0,
    |det Vz|=g(1-v)>0.                                (5)

Include the original row factors (mu2,2mu_(-2-2lambda),-mu_(-2lambda))
and (mu1,mu_(1-2lambda)) at a=1. Their nonzero products multiply these
determinants. For either complete fixed matrix B, the explicit bound

    ||B^-1||_infinity <= max_i sum_j |cofactor_(j,i)(B)|/|det B|  (6)

is finite (indeed equality for this norm). Its entries give the fixed
response coefficients bP,bz,bg in the slow-inverse provider. Multiplying
by the finite bump derivative seminorms gives every FP,Fz,Fg derivative
norm. Psi derivatives follow from its compact integral and x>=l0>0.
More explicitly, for any response F=sum b_j eta_j and every fixed m,

    ||F^(m)|| <= ||eta0^(m)|| sum_j |b_j| g^(-(m+1)j).

For a common positive hull [l,u], take G_m>=||Fg^(m)|| and define

    P_0=u(u-l)G_0/l,
    P_(m+1)=G_m+sum_(i=0)^m binom(m,i) i! l^(-i-1) P_(m-i).

The integral formula bounds ||Psi|| by P_0. Differentiating
Psi'=Fg-Psi/x and using |(1/x)^(i)|<=i!l^(-i-1) proves inductively
||Psi^(m)||<=P_m. This includes every derivative in the gaps, with no
assumption that Psi vanishes there. These finite recurrences supply the
response-profile seminorm slots of the slow inverse directly from the
fixed base bump and matrix coefficients. No target-dependent root or
smallness condition enters this linear map.

The denominator 1-v is retained: this does not prove a bound uniform as
lambda tends to zero. The running band varies with the SAME positive
lambda, so its large reciprocal is a fixed constant. The subsequent
factor a(eta)^-1 still requires the retained background floor a>=a0>0
and its smooth jets. The slow-inverse finite Taylor formulas, with their
three distinct scale powers, then give the response constants. Those
background hypotheses are not proved merely by matrix invertibility.

## 4. Consequence and remaining obligation

For this fixed-choice branch, (1)--(4) supply Q29's actual support/gap
hypotheses, and (5)--(6) supply finite slow-profile coefficients and hulls.
Together with the fixed-profile coordinate and background-floor premises,
these justify band-uniform compact and slow operator constants at every
finite request in Q28's two archives. Their numerical materialization is
not required for a fixed-order analytical estimate. Their value cannot be
inferred from tests on arbitrary rational matrices.

The coefficient-support recheck identifies the full VELOCITY with the
leading reserved velocity on I4, but does not identify pressure there:
Pi1/Pi2 memory persists. All fresh pressure, c_rho P differentiation,
radial base flow, axial A1 repair, covariance and cutoff terms remain in
the source computation. No pressure equality follows from (4).

Next collect the remaining same-profile background, canonical value/jet,
weighted-stress and sharp signed-target premises of the first-source
estimate and resolve their actual bounds, rather than leaving independent
operator-support assumptions that are already consequences of the fixed
choices. Historical localization, all-scale realization, complete residual
improvement and blow-up remain unproved.

Exact tests check support scaling, plateau gap arithmetic, complete
Vandermonde row factors, inverse norms, compact potential support and the
failure of uniform axial conditioning at lambda=0. Frozen Q29 evidence
is replayed unchanged. These fixtures are not actual giant-parameter
suprema. Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` here.
