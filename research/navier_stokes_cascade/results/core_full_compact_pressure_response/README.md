# Full compact pressure and the changed-wave fresh-fast response

Parent `c82464277857a178873fa93df08712ad01d727ce`, verified remotely on
24 October 2026. Its same-choice recomputation CI step passed; run
40105092005 remained in progress without reported failure. Historical
frozen evidence, the chosen rho-squared projectors and old fields are
preserved. This checkpoint derives the FULL auxiliary compact T0 bound
from original raw profile inputs and joins it to the recomputed changed
covariance response. It does not construct the whole new fast/slow block.

## 1. Two averages and two different pressure operators

The wave covariance is averaged in the PHYSICAL angular carrier. Its
six entries may still depend on the auxiliary torus y. Physical angular
cosine/sine cross terms vanish in the first real covariance identity;
that does NOT say the resulting coefficients are auxiliary independent.
For example (1+cos y)cos theta has angular covariance
(1+cos y)^2/2, whose auxiliary Haar mean is 3/4 and whose complement
is nonzero. No omitted physical quadrature cross is restored here.

The preceding same-choice pressure_target_budget controls the Haar
pressure and projected Haar targets. On that domain, exact radial mass
zero makes T0 the ordinary two-edge primitive I0 and A0 vanishes.
The CURRENT full auxiliary-dependent radial source instead uses

    h=gr-rho P, P=integral Haar(gr), pm=T0 h.

Only integral Haar(h)=0 is assured. An ordinary left primitive can have
a right tail for each nonzero auxiliary mode. It cannot replace T0.
The exact radial identity remains D0 pm=gr-rho P-A0 h. The compact
A0 defect generally survives and is not made zero by a Haar moment.

## 2. Retain the actual shifted full compact operator

In U=R^d, on the same positive radial hull and frozen local band, use
the ORIGINAL fixed profiles

    f0=R^(1-d)/d, g0=1, k0=d R^(d-1),
    I_minus=integral_(a<0) (f0 h)(U+a,y+Ma vr) da,
    I_plus =integral_(a>0) (f0 h)(U+a,y+Ma vr) da,
    T0 h=(1-chi) I_minus-chi I_plus,
    A0 h=k0 chi_U (I_minus+I_plus).

The source has the original smooth zero extension. Hold the offset a
fixed when taking coefficient derivatives: y+Ma vr then has no U,Z,T
dependence. Thus no M power is introduced in these coefficient jets.
Graph/physical derivative losses remain separately in the row equations;
they are not removed by this observation.

Let N_n denote the maximum RAW five-variable coefficient derivative
through n. Supply the actual original profile seminorms for 1-chi,
chi, f0 and k0 chi_U. The last includes cutoff derivatives through n+1.
Let Jin bound dU/dx on the full log shell, LU its U integration length,
and gap the ACTUAL log plateau gap for this compact cutoff. These are
lower fixed-profile/hull inequalities, not current fine-residual bounds.

Using the original annular rational bounds Llo<=L<=Lhi, put
r=min(1,Llo/4,gap), q=max(0,ceil((p-3)/2)) and

    D(c)=sum_(j=0)^q binom(q,j)j!(r^2/c)^j/(2c).

The same sufficient symbolic edge constant is

    C_p=D(w^2)exp(4/(Llo-r)^2)
        +D(4)exp(w^2/(Llo-r)^2)
        +Lhi E_p exp((w^2+4)/r^2).

On the left collar chi=0 and the forward integral has the left-edge
bound; on the right chi=1 and the backward integral has the right-edge
bound. Derivatives of the opposite cutoff vanish there. On the middle
region the retained floor gives the third branch. Absolute values of
the shifted y input obey the SAME source bound, uniformly over y.
Consequently the raw weighted T0 coefficient is bounded by

    KT_raw=Jin C_p 2^(3n) [N_n(1-chi)+N_n(chi)] N_n(f0). (1)

This is the existing general compact-map proof specialized to e=0,
not an extension of the Haar zero-moment shortcut. All multipliers and
cutoff derivatives are retained. No smallness condition on M is added.

On the actual A0 output support, let
Zgap=exp((w^2+4)/gap^2). Unweighting the input by the original E_p
and using its whole-line integration length gives

    KA_raw=Zgap E_p 2^(3n) LU N_n(k0 chi_U) N_n(f0).       (2)

The output is compact on the cutoff support. It is not asserted zero,
flat to a chosen epsilon order, or free of its moving cutoff derivatives.
The stronger raw integration-by-parts estimate for A0 still needs its
actual high derivatives and strict epsilon margin.

## 3. Convert native normalized source and output jets explicitly

Supply Bin,Bout bounding the ORIGINAL raw same-weight R-to-U and U-to-R
coefficient jet transfers at order n. They concern the same physical
point, so they change neither zeta nor p. A normalized R through norm
H bounds the raw R maximum by n! H; the raw output maximum bounds the
normalized five-variable through norm by binom(n+5,5) times that maximum.
Thus multiply (1) and (2) by

    Cconvert=binom(n+5,5) n! Bin Bout.                    (3)

This deliberately safe conversion includes radial, slow and auxiliary
derivatives; the three-variable Haar conversion is insufficient here.
Fixed chart-transfer bounds remain retained profile premises. No source
derivative is lost, no extra pole is introduced, and giant exponentials
stay symbolic. Constants are not asserted uniform as n or p tends to
infinity. The physical width, plateau and accepted band are unchanged.

## 4. Join the new COMPLETE changed covariance to full pressure

FullCompactPressure requires the identical RecomputedSigned state and
derives (1)--(3); it cannot be replaced by the Haar-only maps dictionary.
For row output n it requests the new complete angular covariance bound
D at N=n+3 and full T0 coverage through n+1, at a pole covering D.
All six entries and signed squares remain in D; its auxiliary dependence
is retained rather than Haar averaged away.

Using the original graph and reciprocal-R coefficient recipes, let
d0=n+4, g bound the graph coefficient through N, and r bound R^-1.
The historical complete row-response functional gives

    Gr=(d0 g+d0+2r)D, Pd=Kmom Gr,
    Pm=KT0(Gr+rho Pd),
    Et=(d0 g+d0+2r)D,
    Ez=(d0 g+d0+r)D+d0 Pm,
    ||A0 h||<=KA0(Gr+rho Pd),
    ||Er_change||<=rho Pd+KA0(Gr+rho Pd).                (4)

Kmom and the moving rho norm are DERIVED from the SAME pressure/profile
container, on the new pole. Formula (4) retains the full compensated
pressure, not its Haar component. The radial response uses the exact
compact identity, retaining BOTH the compensator and A0. No fine radial
gain is inferred from its unintegrated operator norm. The required source derivatives fit
the N=n+3 covariance ledger: pressure through n+1 uses the corresponding
radial derivative of covariance; fixed coefficient norms use that larger
N. This conditional changed-wave contribution has order alpha-kappa
under the retained coefficient-class hypotheses.

This is NOT the whole current Qcirc fast source. Its remaining input
is the full old completed initial row, with its actual A0/A1 tails and
same-pressure cancellation. The sufficient next inequality is

    ||Qcirc E_postSigned||_n
      <=2 [||E_initial_completed||_n + Et_or_Ez(n)],

with the first term bounded at the required original order on this full
auxiliary domain. A Haar-row bound cannot replace it. Then the actual
temporal/compact fast inverse must use the source through m+4/m+5,
retaining gd-A1 gd, and the slow response must use CURRENT scalar defects.
The existing high-jet/domain ledger and zero strict fast-difference
margin remain; this low-order reproduction does not certify them.

Next join that original full-row cancellation and temporal inverse to
(4), then reconstruct postfast pressure/scalars and the slow correction.
Actual fixed raw profile/chart enclosures, root/value/support premises,
uniform generated fields, all viscous/cutoff/interacting residuals and
blow-up remain unresolved. No unconditional pulse or full mean-block
closure is claimed.

## Evidence

Exact diagnostics distinguish physical angular covariance from auxiliary
Haar averaging, verify the shifted compact differential identity with
nonzero A0, and check fixed-offset coefficient derivatives without M
loss. Marked manufactured lower norms evaluate the complete new signed
covariance and (4). Domain controls reject short pressure orders/poles
and Haar-only substitutes. Historical frozen sources are replayed;
these finite algebraic checks are not continuum enclosures.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json

For the nonzero-mode diagnostic, cos(V)>=1-V^2/2 on [0,1] gives
a strictly positive rational lower bound for the real part of its
Fourier moment against V^4(1-V)^4. With nonzero chi′ at the midpoint,
this proves the A0 defect is nonzero without a floating-point zero test.
