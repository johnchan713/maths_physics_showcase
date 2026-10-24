# Native Haar pressure and weighted target primitives

Parent `67f59dfe8743de3c092abbddb574a931f95b8f8a`, checked against the
remote branch on 24 October 2026. Its new compensator CI step passed;
run 40087871548 remained in progress with no reported failure.
This checkpoint discharges the Kp and oriented target-primitive inputs
of that conditional signed-pressure join. It uses the existing two-edge
kernel theorem, not a new primitive lemma. The new estimate works
DIRECTLY in native R,Z,T jets and then substitutes the complete new pole.
No profile, width, support, accepted band or pressure gauge changes.

## Exact domain and zero moment

Let h be an auxiliary-independent smooth zero-extended source on the
original moving annulus, whose union fits the same positive norm hull
[a,b]. Assume integral h(R,Z,T)dR=0 at EVERY slow point. For every pure
slow derivative this integral also vanishes: use fixed limits containing
the support union, with smooth flat extensions. Thus

    I0 h(R)=integral_left^R h=-integral_R^right h.

For pressure, h=Haar(f)-rho integral Haar(f), with the SAME unit-mass
rho and full current f. This is not h=f-rho integral Haar(f) for an
arbitrary oscillatory f. Global or auxiliary-averaged cancellation alone
does NOT make the latter pointwise compact. The parent's displayed
primitive obligation must be read on its Haar domain; this checkpoint
makes that restriction explicit and tests the excluded oscillatory case.
No A0/A1 oscillatory remainder or flat pressure tail is removed.

For targets the original e=1,2 Haar source has exact weighted moment
integral R^e h=0. Its shifted compact map reduces exactly to

    T_e h=R^-e I0(R^e h).

This uses the original auxiliary-independent zero-moment identity,
so the cutoff dependence cancels and A_e vanishes in this domain.
It is not a formula for a general source before Haar averaging.

## Direct native derivative bound

Use the existing x=log(R^2/(2s Xa)), L=log(Xb/Xa), width w and
zeta=exp(-w^2/x^2-4/(L-x)^2), delta=min(1,x,L-x). For every fixed
slow point dR=(R/2)dx. Apply the existing scalar two-edge theorem to
(R/2) partial_Z^j partial_T^k h, WITHOUT differentiating this logarithmic
pullback in Z,T. If the native source through norm is at most
H zeta delta^-p, its primitive pure-slow jets are bounded by

    (b/2) C_p H zeta delta^-p.

The source bound here is for derivatives at FIXED R; no motion derivative
is omitted. Differentiation under the fixed-limit native integral gives
those exact derivatives before the scalar kernel is applied pointwise.
The backward formula is available for each such source by exact mass
cancellation. This avoids any need for a differentiated chart transfer.

For a normalized jet with i>=1 radial derivatives,

    partial_R^i partial_Z^j partial_T^k I0h/(i!j!k!)
      =partial_R^(i-1) partial_Z^j partial_T^k h/(i!j!k!).

The extra normalized factor is 1/i<=1. Summing the source subfamilies
therefore yields, for the factorial-normalized full through norm,

    ||I0 h||_(n,weighted,p) <= [1+(b/2)C_p] ||h||_(n,weighted,p). (1)

There is no source derivative loss or additional edge pole. The source
must already have all native derivatives through n at this pole.
This does not upgrade the previous input enclosures.

## Exact finite symbolic kernel, then both target channels

Recover the SAME annular ratio from the existing target geometry, and
use its rational bounds Llo<=L<=Lhi. The API checks the compensator's
annular ratio agrees. Let r=min(1,Llo/4), q=max(0,ceil((p-3)/2)),

    D(c)=sum_(j=0)^q binom(q,j)j!(r^2/c)^j/(2c).

A sufficient C_p is the SUM of these three positive branches:

    D(w^2) exp(4/(Llo-r)^2),
    D(4) exp(w^2/(Llo-r)^2),
    Lhi E_p exp((w^2+4)/r^2).

The original theorem uses a maximum; the sum dominates it. Replacing L
by the lower denominator bound and upper length bound preserves every
inequality. The proof collar is smaller, but the physical support is not.
Huge exponentials stay symbolic. E_p is the inherited finite edge
absorption factorial recipe; no uniformity as n,p grow is asserted.

Native radial multiplication has exact through-order upper bounds

    A_e(n)=sum_(j=0)^min(e,n) binom(e,j)b^(e-j),
    B_e(n)=sum_(j=0)^n binom(e+j-1,j)a^(-e-j).

The factorial-normalized product inequality and (1) give

    ||T_e h||_n <= A_e(n)B_e(n)[1+(b/2)C_p] ||h||_n.       (2)

The sum of the two channel constants is a common sufficient primitive
coefficient for the preceding full pressure_target_budget. Multiplying
by R^e changes no endpoint pole on this positive fixed hull.

## Current complete pressure join and remaining obligation

NativeHaarPrimitives wraps the SAME MovingCompensator and RealCovariance
state. At the actual requested new order/pole it supplies (1), (2),
the parent's moving rho and c_rho recipes, and its all-weight Kmom.
Only the ACTUAL weighted projector bump remains a map callback. No
current pressure, covariance entry, signed square, scalar or mean is
reset; the whole epsilon partial_Z(Jz+c_rho P) is retained.

The actual fixed rho profile, exact moment/support identifications and
its finite seminorms remain premises. Most critically the actual sigma_e
projector family still needs a uniform bound for its full native jets
and moving support. The normalization-only counterfamily from
core_log_pullback_projector_join still blocks inferring this from mass
and interior support alone. Deriving rho's jets does not identify sigma_e.
Next justify that family, then rebuild fast/slow means from CURRENT
post-signed rows and recomputed scalar defects. Complete actual pulse
enclosures, full viscous/cutoff/interacting residual and blow-up remain
unproved. This checkpoint proves conditional operator estimates, not a
continuum realization or a residual gain.

## Evidence

Exact controls check the positive kernel branches against the retained
independent recipe, the radial multiplier jets, and native primitive
and weighted-target differential identities on endpoint-flat polynomial
diagnostics. A nonzero angular mass example rejects the invalid
Haar-only compensation extension. Manufactured norms check the complete
pressure substitution and whole current scalar memory. Historical frozen
sources are bound and replayed unchanged. Tests are finite diagnostics.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
