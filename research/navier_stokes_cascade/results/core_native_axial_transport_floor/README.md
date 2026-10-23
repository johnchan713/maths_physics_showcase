# Native axial transport and the unchanged-band normal floor

Parent: `9151274fbad956a8e2cabdf99b1e06af919ceabf`.
Its finite-support-closure CI step passed when inspected on 23 October
2026; the complete workflow was still running without a failed step.

The parent leaves the actual weighted inputs and original-band gaps open.
Here we reduce one of those value gates, the FULL tangential normal floor,
to a directional first-jet estimate for the two native LEADING profiles
E0,U0 and the SAME retained complete low-order Euler diagonal. We derive
the exact axial derivative before bounding it, preserve different powers
of S in the error, and prove the original slow gate from its retained
scalar band hypothesis without evaluating the tiny physical scales.

These are analytical identities and conditional sufficient estimates.
The actual leading first-jet enclosures, diagonal identity and reference
gaps are not numerically or unconditionally certified. No new band,
length, representative, cutoff, carrier, root, width C_setup^-20, pressure
gauge or physical field is selected. Historical evidence stays frozen.
The complete residual and physical decay orders are unchanged.

## 1. Differentiate the actual background in the native variables

Retain 0<h<1/100, A=1/2+h, D=1/2-h and the original physical chart

    s=q/Q, s-Z^2 s^(2h)=T, X=R^2/(2s), eta=Z s^-D,
    1/4<=s<=4, |eta|<=1, L=1-2h eta^2>49/50.               (1)

The original enlarged norm hull has Rmin<=R<=Rmax and therefore
`Rmin^2/8<=X<=2 Rmax^2`. It need not have positive growth/cone floors
away from the ORIGINAL active representative; only upper profile norms
are used there. At T=0 all statements use one-sided physical derivatives.

Let E(q,X,eta)=q^A utheta_B and U(q,X,eta)=q^A uz_B be the phase's
COMPLETE native background components, including every Borel coefficient
and cutoff.
Their actual phase coefficients are

    F=s^-A E/R, G=s^-A U,
    H=p F+pz G.                                             (2)

All label parameters, Q, carriers p,pz and representative choices stay
frozen. Euler_q means q partial_q at FIXED native X,eta. Implicit
differentiation, not freezing q inside a cutoff, gives

    partial_Z log s=2eta s^-D/L,
    partial_Z X=-2eta X s^-D/L,
    partial_Z eta=(1-eta^2)s^-D/L.                          (3)

Since A+D=1, define

    D_h=2eta(Euler_q-X partial_X-A)+(1-eta^2)partial_eta.

The EXACT identities on the entire physical norm hull are

    F_Z=s^-1 D_h E/(R L), G_Z=s^-1 D_h U/L,
    H_Z=s^-1[(p/R)D_h E+pz D_h U]/L.                       (4)

In particular Euler_q differentiates BOTH q^(2nh) and chi(c_n q).
Omitting it would lose a term of the complete axial phase transport.
The radial streamfunction primitive M is absent from (4) because these
are the angular/axial velocity components, not because M or its pressure
memory has been set to zero in the field. Radial cutoff-curl contributions
remain in the retained complete diagonal and all other phase equations.

## 2. The leading directional norm, including an exact heat cancellation

The LEADING E0,U0 are native q-independent. For either B=E0 or U0,
use actual outward first-jet enclosures on the full hull:

    |B|<=v, |B_X|<=r, |B_eta|<=e, X<=Xhi=2Rmax^2.

Put a=A v+Xhi r and b=e. The coefficient of |eta| is bounded by
`2a t+b(1-t^2)`, 0<=t<=1, whose exact maximum is

    D(a,b)=2a                    if a>=b,
           b+a^2/b              if 0<=a<b.                  (5)

Thus |D_h E0|<=DE and |D_h U0|<=DU follow directly from primitive
profile bounds, not a desired fine H_Z callback. Formula (5) is sharp
for independent component enclosures; their actual profile extrema need
not occur together. It never exceeds the looser 2a+b norm sum. No leading
second jet, first Borel coefficient, cutoff K2 or independent radial-M
norm is needed for THIS floor route. The same infinite Euler diagonal
remains a premise; it has not been justified from these first-jet bounds.

There is extra EXACT structure in the retained heat exterior. The actual
physical swirl is C r^(-1-2h) Hheat(4tau/r^2), and uz=0. Therefore its
native profiles are

    E0=C(2X)^(-A) Hheat(2(1-eta^2)/X), U0=0.                (6)

Substituting (6) into D_h cancels identically for ANY smooth Hheat:

    D_h E0=D_h U0=0.                                       (7)

This is the expected axial independence of the actual heat field, even
though its individual native X/eta jets are nonzero. A direct enclosure
of the directional combination may therefore improve (5) on that known
region. It does not authorize treating the inner/modulated profile as
heat, or setting an exterior/axial pressure constant to zero. Formula
(5) remains an available whole-hull bound without this refinement.

## 3. Bound the complete correction from the original Euler diagonal

The existing same-background diagonal supplies every positive coefficient
through native Euler/X/eta order one by `2^-n q^(nh)`, including n=1.
Summing the geometric majorants gives, for each such jet,

    |delta E|_1, |delta U|_1 <= q^h/(2-q^h) <= q^h.          (8)

This is conditional on that VERY SAME complete schedule, not a new
prescription. It reads only its E/U components for (4), while the full
radial streamfunction cutoff term stays in the retained field. For
X<=Xhi, the sharp elementary envelope (5) now has
a=1+A+Xhi, b=1, so

    |D_h delta E|,|D_h delta U| <=2(1+A+Xhi) q^h
                                      <=4(1+A+Xhi)epsilon. (9)

Here epsilon=Q^h and q^h=s^h epsilon<2epsilon. No generic analytic
three-coordinate Taylor majorant is needed for this single derivative.
All original quantitative coefficient/diagonal assumptions remain.

Keep the parent's outward carrier bounds pplus>=|p|, qz>=|pz|,
and set K=4/(1-2h), using the SAME whole-hull Rmin. Then

    |H_Z| <= H0+epsilon H1,
    H0=K[(pplus/Rmin)DE+qz DU],
    H1=4K(pplus/Rmin+qz)(1+A+Xhi).                          (10)

These are primitive-to-field bounds for the actual complete H_Z. They
are not numerical enclosures of the selected enormous leading profiles.

## 4. Retain the different S powers in the normal error

Use the SAME reference carrier identity and box comparison as the earlier
tangential-floor checkpoint. Let S=ell^2, Llow S<=Ls<=Lhigh S, Ubound>=u,
lambda_max bound the positive growth, gmin>0 bound reference shear, and
the original box diameter be Dbox S^-3. With B=max(1,lambda_max),

    qz=B[1+Ubound/(Llow gmin)], p0=Rmax qz, pplus=p0+1,
    ct=2B Ubound/(Llow gmin), ca=1/Rmin,
    cb=p0 Dbox/Rmin^2,
    beta=lambda_min/[4(1+Ubound^2)^2].                     (11)

The original nearest-NONZERO rounding still costs 1/k. Both tangential
tilt errors, box variation and full axial transport remain:

    |n_tan-Bs Kunit|
       <=ct/S+ca/k+cb/S^3+Lhigh S epsilon(H0+epsilon H1).

Under the original slow gate S^2(epsilon+epsilon^2+k^-1)<=1 this yields

    error <= A0/S+B0/S^2+C0/S^3,
    A0=ct+Lhigh H0, B0=ca, C0=cb+Lhigh H1.                 (12)

C0 here is just a positive floor-envelope coefficient, not the research
mean-row order 2-2kappa. The sufficient SAME-S floor is

    4[A0/S+B0/S^2+C0/S^3]^2 <= beta.                       (13)

Then |n_tan|>=sqrt(beta)/2 and |n|^-2<=4/beta by the original reference
calibration and triangle inequality. No individual angular or axial
component is asserted to be bounded below.

Collapsing (12) to `(A0+B0+C0)/S` remains valid but loses two powers of S
on the correction term and one on rounding. The retained-power test can
pass when this collapsed SAME positive program fails. This is not a
universal comparison with the older generic-chart route: those programs
read different leading norm enclosures. Where both are justified, either
sufficient route remains available. No failure is a PDE error lower bound.

All coefficients in (12) are nonnegative and fixed under their uniform
input enclosures. Its error decreases with S. Consequently checking (13)
at the actual first band suffices for later bands, once the original
slow gate and uniform input contracts apply. The implementation accepts
only a supplied S and never selects another scale or larger Lreq.

## 5. The original scalar band hypothesis supplies the slow gate

Use the ORIGINAL k=ceil(epsilon^-1/2), so k^-1<=sqrt(epsilon),
Q=2^-ell and epsilon=2^(-h ell). For all ell>0,

    S^2(epsilon+epsilon^2+k^-1) <=3 ell^4 2^(-h ell/2).

Since log 2>1/2, the positive exponential series at order 16 gives

    2^(-h ell/2) <=16! [4/(h ell)]^16.

If the SUPPLIED unchanged band satisfies its retained scalar requirement
ell>=512/h^2, then

    S^2(epsilon+epsilon^2+k^-1)
       <=3*16!*4^16/(h^16 ell^12)
       <=3*16!*4^16*h^8/512^12
       <3*16!/(2^76*100^8) <10^-24<1.                     (14)

The log inequality follows from integrating 1/t>1/2 on (1,2). Thus (14)
is an all-h estimate on 0<h<1/100 and every such later band, not a sample
at a diagnostic h. It uses the rounded carrier via an upper bound; it
does not falsely assume the stepwise exact ceil carrier expression is
monotone. No tiny exponential, giant first integer or new band is
materialized. If the actual first-band prescription is not identified
with ell>=512/h^2, this scalar identification remains open: the code
reports it and does not repair it by increasing ell.

## 6. Exact remaining gate and scope

The normal-floor route is now (13) with explicit native first-jet norms
of E0,U0, the original complete low-order diagonal, actual reference gaps
and the supplied band. Equation (14) removes an independent slow-error
test IF the original scalar band hypothesis is established. Original
representative membership, clearance and enlarged physical coverage stay
as previously proved/required; this result supplies no new negative-time
extension, cone floor on the norm hull or replacement growth reference.

The next obligation is to enclose those actual leading primitive jets and
reference gaps and identify the SAME diagonal/first-band prescription,
then evaluate (13) there. The separate target-freezing/root/cone gates,
weighted current target/source inputs, complete pressure/axial means,
full viscosity/cutoff/interactions and infinite forcing/growth obligations
remain. The full regular order is still E=2-3kappa; no absolute residual
decay, unconditional all-scale pulse bound or blow-up is proved.

Exact symbolic controls differentiate the implicit chart with Euler
dependence, check the heat cancellation, sharp native envelope and
independent full transport arithmetic. Rational norm examples distinguish
retained from collapsed powers and exercise supplied-scale guards.
They are not selected actual profile or physical-band enclosures.
Historical frozen records are bound without modification. Reproduce with
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH_PATH.json`.
