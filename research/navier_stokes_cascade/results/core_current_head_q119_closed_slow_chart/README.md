# Q119: closed real slow-chart factors for the full background

Parent: `b8ec2e08e8076c0e1249bf13d519b119da490e9f` (Q118).
On 3 November 2026 the inspected remote was Q118; research CI 41594347916
and field CI 41594347819 were running. Q117 field CI 41587075430 had passed.
The latest notes, frozen evidence and core_growth_activation scope were
read. The leading-core/first-collar result and original width C_setup^-20
are unchanged. This checkpoint proves a closed bound for the GEOMETRY
factor in the full-background jet conversion and supplies it to Q118's
current phase/source norm interface. Actual profile/tuple/schedule/source
inequalities remain inputs. Historical frozen evidence is preserved.

## 1. Start from the same real coordinate jets

Use core_background_finite_head_budget (6)--(8), with the original
s in [1/4,4], positive hull a<=R<=b, |eta|<=1, and fixed 0<h<1/100.
The source variables are log q, X=R^2/(2s), eta=Z s^-D. Euler q derivatives
are exactly log-q derivatives; they are not frozen q^j partial_q^j.
core_fixed_profile_tame supplies the normalized log-s jet majorant

    Y(z,t)=(1/3)[(1-32z)^-1(1-64t)^-1-1].

Collapse all coordinate increments to one positive variable t. For the
actual zero-constant increments, sufficient coefficientwise majorants are

    Y(t)=(1/3)[(1-32t)^-1(1-64t)^-1-1],
    Xd(t)=2[(b^2+2bt+t^2)exp(Y(t))-b^2],
    Ed(t)=(1+4t)exp(Y(t)/2)-1.                         (1)

These are the SAME coordinate expansions used in the old finite-sum
proof, retaining the radial increment and both implicit slow coordinates.
The constants 2 and 4 retain 1/(2s0)<=2 and s0^-D<=4. They are not new
coordinate maps or assumptions about the coefficient profiles.

For any fixed order m, bound the profile's raw derivatives through m in
(log q,X,eta) by M_m. Its finite Taylor composition is bounded through m by

    M_m sum_(i+j+k<=m) Y(t)^i Xd(t)^j Ed(t)^k/(i!j!k!). (2)

No profile derivatives above m enter a coefficient of degree <=m. The
positive infinite series exp(Y+Xd+Ed) can therefore majorize (2) through m
without asserting convergence of the PROFILE's Taylor series. Only the
already bounded coordinate majorants are evaluated below.

## 2. A fixed evaluation radius and explicit raw factors

Set

    tau=min(1,a)/[1024(1+b)^2].                         (3)

This is solely a generating-function evaluation point. It changes no
physical scale, cutoff width, accepted band or q-domain. In particular
tau<=1/1024 and tau<=a/2. For 0<=t<=1/1024, rational denominator bounds give

    Y(t)<=36t.

Indeed (1-32t)(1-64t)>=(31/32)(15/16), and the numerator of the product
inverse minus one is at most 96t. Thus Y(tau)<1/16. From the nonnegative
exponential series, exp(Y)<=2 and exp(Y)-1<=2Y there. Consequently

    Xd(tau)<=144 b^2 tau+8b tau+4tau^2
                  <=156(1+b)^2 tau<1/4,
    Ed(tau)<=44tau<1/16,
    Y(tau)+Xd(tau)+Ed(tau)<3/8<1.                      (4)

The composition series at tau is less than e<3. The original field
prefactors s^-A/R, s^-A and s^-1/2 have base bounds 4/R,4,2 and log
exponents bounded in absolute value by 1,1,1/2. The reciprocal radial
series at tau is

    sum_(j>=0) a^(-j-1) tau^j<=2/a.

Multiplying the majorants gives respective values at tau bounded by
48/a, 24 and 12. Positivity bounds each normalized coefficient of degree
j by that value times tau^-j. A raw mixed derivative of total order j
costs at most j!. Since tau<=1, every j<=m is thus covered by

    CF_m=(48/a) m! tau^-m,
    CG_m=24 m! tau^-m,
    CB_m=12 m! tau^-m.                                (5)

This is an all-fixed-order estimate, uniform in the running band under
the original fixed-hull/chart hypotheses. The dimension count is already
in the collapsed normalized coefficient sums. The profile multiplier M_m
remains arbitrary: (5) proves no analytic or Gevrey bound for the actual
profiles. It bounds only the coordinate-composition geometry.

## 3. Retain the complete tuple, head, tail and schedule

For the SAME full normalized tuple and fixed Borel cutoff/schedule, let

    C_(n,m)=max(1,2K_(m+1)(1+2nh)^m M_(n,m)),
    N=max(2,m), Bcorr_m=sum_(1<=n<N) C_(n,m)+2^(1-N).

Under the retained normalized diagonal C_(n,m) a_n^(-nh)<=2^-n for m<=n,
the old proof bounds every Euler/X/eta derivative of the correction by
Bcorr_m q^(2h). The leading tuple has its original bound Lead_m. Therefore
(5) gives the sufficient full slow-jet estimates

    F_m<=CF_m(Lead_m+Bcorr_m),
    G_m<=CG_m(Lead_m+Bcorr_m),
    btilde_m<=CB_m(Lead_m+Bcorr_m).                      (6)

Both radial tuple entries W_n^r and J_n^r remain, including the radial
cutoff-curl term. The complete Pi1/Pi2 histories and moment repairs are
inside the original tuple norms; they are not replaced by leading values.
The extra Borel Euler-cutoff derivative m+1 and coefficient ONE in the
finite head at low order remain indispensable. No a_n or cutoff is chosen
or modified by this estimate.

`chart_recipe` stores the factors in (5) as finite instructions. Its default
materialization helper rejects large factorial orders. The norm provider
can evaluate a finite supplied order explicitly; this run evaluates only
moderate manufactured controls, never the actual giant source graph.
No uniform native norm is inferred from the mere finiteness of a recipe.

## 4. Preserve the accepted value proof while using the new high jets

`ClosedBackground` subclasses the EXACT BackgroundBudget type used by Q118.
It retains its callbacks, h, Euler convention, geometry, finite head and
geometric tail. At orders ZERO, ONE and TWO it delegates to the original
budget without alteration. Its inherited comparison recipe uses only
these low orders, so the original seven-term Cphi and low value constants
are identical. No enlarged replacement constant is inserted into a frozen
normal-floor or length selector. Other selector contracts retain their
own original recipes and proof obligations.

At m>=3 the provider uses (5)--(6). These are independently proved upper
bounds, not asserted to dominate the old finite-sum constants. The resulting
sufficient constants need not be monotone across the order-two/order-three
boundary. Actual through-m norms are still bounded by (6) at each order;
no monotonicity of independently supplied norm bounds is assumed.

The subclass passes directly into Q118.BackgroundPhase. Its same-hull,
fixed-ceiling R->U->R conversions, varying reference-rate and linear slots
then feed Q117's current covariant inverse. The full Gaussian restoration,
pressure, potential, complete curl and covariance interfaces retain their
existing contracts and pressure memory. All additional current root/target,
normal value, frame, path, support and regular/cut pulse hypotheses remain.

## 5. Verification and next obligation

Exact rational checks verify (3)--(4) on hulls both below and above R=1.
Independent normalized series algebra computes the full composition and
prefactor series in (1)--(2), including all products, through order nine;
the raw coefficients obey (5). Other controls compare every original
low-order result and Cphi exactly, retain head/tail/cutoff traces at order
four, and run the current phase/complete response interface with the closed
budget. They use manufactured tuple inputs and are not canonical-field
norm enclosures. The analytical proof is (1)--(6); the tests do not prove
the native inputs. Frozen replay preserves Q118 and its evidence graph.

Next establish the actual leading/finalized tuple seminorms, original Borel
cutoff identity and schedule gates, and remaining path/pulse/source bounds
on the same proof domain. This checkpoint removes growing finite sums from
the slow-chart GEOMETRY factor; it does not supply those physical data.
The source-order-100005/high-pole radial obligation, full cone join,
coefficient/reserved-mean supports, unconditional all-scale pulse bounds,
complete viscous/cutoff/interaction residual improvement, infinite
realization and blow-up remain unproved.
