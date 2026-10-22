# Euler certificate for the implemented Borel cutoff and the unchanged schedule

Parent: `dbc9aeb1696c6a93b5128a5a04f95f76c528db2e`.
The preceding full radial reconstruction supplies the low native tuple
constants needed by Cphi, but retains the Borel cutoff Euler norm and the
normalized diagonal schedule. This checkpoint proves an analytic Euler
certificate for the EXACT cutoff already implemented in the background
curl checks, gives explicit constants at every fixed order, and joins its
second-order bound to the same seven-term phase comparison. It also adds a
verifier for supplied intervals of the existing cutoff scales; it never
selects replacement scales.

The elementary cutoff theorem is independent of the Navier--Stokes profile
premises. Its application to the summed field still requires that field to
use this same implemented chi, its full coefficient bounds and its retained
schedule. The historical mathematical specification permits a qualitative
smooth chi; flatness alone cannot identify it with a particular function.
No cutoff, profile, potential gauge, pressure memory, width C^-20, radius,
carrier, length budget or accepted band is changed here. All historical
files and frozen evidence remain unchanged. No unconditional pulse, complete
residual improvement or blow-up claim follows.

## 1. Identify the implemented function before assigning a bound

`background_summation/checks.py:cutoff_derivatives` uses exactly

    chi(s)=1-sigma(2s-1), 1/2<s<1,
    chi=1 for s<=1/2, chi=0 for s>=1,
    sigma(x)=g(x)/[g(x)+g(1-x)], g(x)=exp(-1/x^2),          (1)

with g extended by zero for x<=0. The imported squared-flat sigma is the
same `outer_pressure_pilot/schedule.py:step`, whose interior formula is
1/[1+exp(x^-2-(1-x)^-2)]. Its reflection branch represents the same exact
real function. AST checks bind these source expressions. Their floating
outputs are NOT used to certify a supremum.

This is the BOREL cutoff chi(a_n q), not the radial coefficient cutoff
kappa(X) or the modulation variance cutoff. Their arguments and norms
are different even when they use related smooth-step constructions.
If another unspecified qualitative chi was used to define a historical
field, (1) must be identified with it before the numerical bounds below
can be applied to that field. This checkpoint does not select a new chi.

Indeed the qualitative plateau/value/flatness conditions cannot supply a
universal K2: for N>=1 use the auxiliary monotone cutoffs
1-sigma((s-a_N)/w_N), w_N=1/(4N), a_N=3/4-w_N/2. Each transitions inside
[1/2,1], but |E chi_N(3/4)|=6/w_N=24N. They are an obstruction to an
invalid norm inference, not replacement cutoffs for the actual background.

## 2. Explicit all-fixed-order real derivative bounds

For u=1/x>=1 and 0<x<=1, differentiation gives exactly

    g^(j)(x)=P_j(u) exp(-u^2),
    P_0=1, P_(j+1)=-u^2 P_j'(u)+2u^3 P_j(u).              (2)

For a monomial u^d put r=ceil(d/2). Since u>=1 and
v^r exp(-v)<=r! for v>=0,

    A_j=sum_d |[u^d]P_j| ceil(d/2)!
       >=sup_(0<=x<=1)|g^(j)(x)|.                         (3)

All derivatives extend as zero at x=0: the polynomial times exp(-u^2)
tends to zero, using the same exponential estimate with a strictly higher
power. Thus this argument also proves the flat endpoint extension at every
fixed order; it assumes no complex neighborhood at the flat endpoint.

At least one of x,1-x is >=1/2, so H(x)=g(x)+g(1-x)>=exp(-4)>1/81.
The last inequality follows from e<3. Also |H^(j)|<=2A_j.
The RAW reciprocal/product rules therefore give the explicit finite bounds

    B_0=81,
    B_j=81 sum_(k=1)^j binom(j,k) 2A_k B_(j-k),
    S_j=sum_(k=0)^j binom(j,k) A_k B_(j-k).                (4)

These bound H^-1 and sigma respectively. For j=0 use the sharper
0<=sigma<=1. Reflection supplies the other flat endpoint. Equations
(2)--(4) prove finiteness with integer coefficients at every fixed order,
not a constant uniform in derivative order.

For the two needed low derivatives the existing direct proof in
`cutoff_repair_derivative_budget`, section 3, is much sharper. On x<=1/2,
z=x^-2>=4, logistic differentiation gives

    |sigma'|<=exp(4-z)4z^(3/2)<=32,
    |sigma''|<=exp(4-z)(12z^2+16z^3)<=1216.                (5)

Each power times exp(-z) decreases for z>=4. Reflection and flatness cover
the whole interval. These are analytic bounds for the SAME sigma in (1),
not sampled peaks. The code retains the independent coarse route (4) and
uses the smaller proved low bounds (5).

Let E=s partial_s and S(m,j) be the Stirling numbers of the second kind.
The exact operator identity and affine chain rule imply

    E^m chi=sum_(j=1)^m S(m,j)s^j partial_s^j chi,
    sup|E^m chi|<=T_m=sum_(j=1)^m S(m,j)2^j S_j,
    T_0=1, K_m=max_(0<=j<=m)T_j.                          (6)

On the transition, s<=1; on either plateau all positive derivatives vanish.
Thus K_m applies on the ENTIRE support, independently of a_n. In particular

    K_0<=1, K_1<=64, K_2<=64+4*1216=4928.                 (7)

The exact midpoint s=3/4 has chi=1/2, chi_s=-16, chi_ss=0,
and E chi=E^2 chi=-12. Therefore the actual padded through-two norm is
at least 12. Omitting the s chi_s contribution to E^2 chi would falsely
give zero there. Ordinary and Euler second derivatives are not interchangeable.

The preceding radial checkpoint deliberately used K2=3 as a MANUFACTURED
rational input, with actual-cutoff certification left false. It is not an
upper bound for (1). The new same-function substitution uses 4928 instead;
no historical fixture is relabeled as actual evidence or overwritten.

## 3. Verify supplied scales without choosing new ones

The parent requires, for each coefficient n and every m<=n,

    C_(n,m) a_n^(-nh)<=2^-n,
    C_(n,m)=max(1,2K_(m+1)(1+2nh)^m M_(n,m)).              (8)

The extra cutoff derivative m+1 remains because the velocity is the curl of
a cutoff potential. Full finalized coefficient/primitive norms M_(n,m),
including complete pressure and moment repairs, enter before summation.
No velocity or pressure datum is reset at a Borel collar or reserved annulus.

Suppose supplied rational intervals bound the EXISTING b_n=log2(a_n):
l_n<=b_n<=u_n. Let gamma_n>0 be a certified lower bound for nh, and
let c_(n,m)=ceil(log2 C_(n,m)), computed exactly from a padded rational
upper enclosure for C_(n,m). Sufficient, explicitly checkable gates are

    l_n>=u_(n-1)+1,
    l_n>=Bdomain+1, Bdomain>=log2(q0^-1),
    gamma_n l_n>=n+c_(n,m), every m=0,...,n.               (9)

Then a_n>=2a_(n-1), the common domain gate holds, and, for every
0<q<=a_n^-1,

    C_(n,m)q^(nh)<=C_(n,m)2^(-gamma_n l_n)<=2^-n.          (10)

This is an all-smaller-q inequality, not endpoint sampling. Constants
outside the support vanish by (1). If h or a_n are not rational, their
outward gamma/log intervals must be justified analytically before the
verifier is used; a caller assertion is not interval certification.

`verify_schedule` requires all normalized rows through n and consecutive
indices, checks the supplied intervals, and reports failures. It creates
no scales and does not enlarge a_n. A passing FINITE supplied prefix is not
an infinite schedule certificate. All original physical, potential,
weighted-stress and augmented-residual schedule constraints also remain;
(9) does not replace them.

The existing natural-log prescription in `background_jet_foundation`,
equation (11), proves (8) analytically if that SAME prescription was applied
to bounds dominating the constants in (8), together with all its original
requirements. A later, more conservative derivative upper bound need NOT
be dominated by an already fixed earlier scale. It must pass a sufficient
check such as (9), or the actual sharper bound must be proved. This code
neither retroactively identifies symbolic/diagnostic scales with the actual
infinite sequence nor silently reruns the old selection.

## 4. Exact dependence of the complete Cphi recipe on K2

Fix the parent's full native leading through-two tuple bounds, finalized
first through-one tuple bound M11 and the entire-chart geometry. The radial
checkpoint retains the same first pressure/viscosity/moment inverse. Its
Cphi dependence on a padded cutoff norm K is exactly

    Cphi(K)=D+Ecost max(1,a K),
    a=2(1+2h)M11,
    Ecost=5 Lhigh(pplus+qz)(CF_1+CG_1), D>=0.              (11)

The factor FIVE is essential: the radial-background term contains FOUR
times the correction norm, while the complete axial-transport term contains
ONE. D retains tilt, radial/azimuthal box variation, both rounding costs,
leading axial transport and the parent's half-tail. The code checks (11)
against the full seven-group formula, including the saturated aK<=1 branch
and the M11=0 case. It does not discard pressure memory or primitive inputs
to make the cutoff cost small.

For (1), 12<=K2<=4928. Thus, with every other input held fixed, the
UNCHANGED sufficient normal-floor test obeys the two rational checks

    necessary for this recipe:
       beta Lreq^2>=144 r0^2 Cphi(12)^2,
    sufficient cutoff enclosure for this recipe:
       beta Lreq^2>=144 r0^2 Cphi(4928)^2,
    beta=lambda_min/[4(1+U^2)^2].                          (12)

Failure of the first means no valid cutoff upper enclosure can pass THIS
recipe with the given other bounds and original Lreq. Failure of the second
alone can reflect conservative upper bounds. Neither proves failure of the
actual normal floor, stress cone or PDE construction. Passing either does
not verify any other selector gate or change Lreq/band. In particular a
finite Cphi is not itself proof of the original first-band inequality.

## 5. Evidence, remaining premises and next obligation

The exact controls bind the implemented source expressions, test rational
Gaussian and Euler jets, preserve the missing first term in E^2, verify
supplied-scale failure cases, and reproduce every group of (11). They use
manufactured coefficient and log-scale inputs. No enormous selected field
or actual infinite cutoff sequence is evaluated. The elementary derivations
(2)--(12), not finite tests, justify their continuum implications.

The audit binds six frozen evidence records and all their declared source
hashes. Reproduce with `python test_audit.py -v` and `python audit.py
--verify-record evidence.json --output /path/to/fresh.json`. Old evidence
cannot be overwritten. Numerical historical curl controls remain numerical
diagnostics, not interval proofs of derivative suprema.

Next identify the actual summed field's chi and normalized scale prescription
with the above certificate, retain all physical/weighted schedule gates,
verify the entire native chart coverage and actual leading gaps, and test
(12) using the original selector. The tangent propagator, full supports and
fixed maps must then be discharged before the first sourced complete-wave
covariance and fresh full-pressure/axial means are promoted. Uniform actual
pulse bounds, the full viscous/cutoff/interaction residual and blow-up remain
unproved.
