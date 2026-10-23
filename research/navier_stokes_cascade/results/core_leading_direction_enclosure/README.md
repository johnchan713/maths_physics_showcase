# Leading native directional enclosure on the original whole physical hull

Parent: `102f6fed681ccffc474887df75c99f40d16ed259`.
Status: `regional-leading-first-jets-joined-with-same-selector-margin-open`.

The parent reduced the tangential normal floor to native leading E0/U0
first jets, the complete first-order Euler diagonal and one same-band
inequality. This checkpoint supplies explicit analytical leading-direction
enclosures from the EXISTING regional bounds. In particular the first
logarithmic-radial derivative of the compact modulation has NO N loss.
It also identifies the low-order diagonal and the slow gate of the already
defined primary tail. The remaining normal-margin inequality is (9).
Its actual parameter substitution is not certified here.

These are implications of the reviewed core, zero-free collar, same-root
loop/repair and exterior defining-map lemmas. They are not an independent
unconditional verification of those foundations. Rational controls and
symbolic differentiations below do not evaluate the giant selected profile
or prove a continuum quantifier by sampling. All historical records remain
unchanged. No cutoff, width C^-20, frequency, carrier, tilt, root, pressure
datum, physical hull, accepted band or selector constant is reselected.

## 1. Explicit axis-side inputs, including the regular swirl conversion

Distinguish the core coefficient norm S_core from the physical S=ell^2.
Let w=tc/2<t/20, t=C^-20, delta=rho/(800 S_core), Xa=4/Lambda,
x=4.1/Lambda, b=1+x. Retain the zero-free core and collar estimates of
[collar_operator_envelope](../collar_operator_envelope/README.md), (1)--(4).
Their extension to THIS w was justified in
[global_axial_angular_bound](../global_axial_angular_bound/README.md), §4;
4 exp(w)<4.1 still holds. Write

    B=2(256 S_core+9+S_core/Lambda)(1+Lambda)^2,
    j1=25 S_core, j2=1600 S_core/t+1000 S_core^2,
    K=exp(25 w S_core) B(1+w x)
        *[1+(x+j1)/Xa
          +(2x+x^2+2xj1+j2+j1^2+j1+64x/t)/Xa^2].       (1)

For the actual regular swirl Freg and axial U0 on 0<=X<=Xa exp(w),
the common complex collar bound and angular Cauchy give

    |Freg|, |Freg_X|, |U0|, |U0_X| <= K,
    |Freg_eta|, |U0_eta| <= K/delta.                     (2)

No division by X at the axis is used. The NATIVE swirl is E0=sqrt(2X)Freg,
not Freg. With y=log X, E0_y=sqrt(2X)(Freg/2+X Freg_X).
Since X<=x and sqrt(2X)<=1+x=b, (2) bounds the POINTWISE first-jet sums by

    JE_axis=b K(3/2+x+delta^-1),
    JU_axis=K(1+x+delta^-1).                             (3)

E0 need not be smooth as a function of X AT X=0; all physical norm hulls
have Xlo=Rmin^2/8>0. Its logarithmic derivative has the displayed regular
limit. Ordinary X derivatives on a specified positive hull also follow from
|E0_X|<=b K(1/(2Xlo)+1), while |U0_X|<=K. The directional route below
avoids this unnecessary inverse-Xlo loss.

## 2. Compact first jets do not cost the large frequency

Keep Aenv=C^4096, Ahi=Aenv^256, G=exp(Aenv^32), H=exp(G^4),
and the SAME N=1+floor(H^32). In particular Ahi<H, H>2^1024 and
H^32<N<2H^32. The original fields have mixed J4 norm <=Ahi and the
joint periodic primitives P have J3(y,eta,phi) norm <H by
[higher_jet_budget](../higher_jet_budget/README.md), §§3--5.
The actual modulations are E exp(beta_A), U+beta_B, where

    beta=P(y,eta,Ny)/N,
    beta_y=P_y/N+P_phi, beta_eta=P_eta/N.                (4)

The joint factorial norm bounds the SUM of the four relevant primitive
seminorms, so |beta|J1<=H and |beta|<=H/N<1. Hence
|exp(beta)|J1<=3(1+H) and the modulated E,U J1 bounds are
3 Ahi(1+H), Ahi+H. This is a first-order calculation; differentiating
again gives N P_phiphi and does NOT have the same N-independent bound.

The same leading repair has J3 norm <H^20/N<1, as established in
[global_axial_source_envelope](../global_axial_source_envelope/README.md),
§1. It retains the original inverse-width bump costs, small root and full
incoming pressure/moment memory. Consequently, on the entire compact
modulation and first repair patch,

    |E0|J1 < H^4, |U0|J1 < H^2.                         (5)

For H>=16, the elementary inequalities
3 H(1+H)+1<H^4 and 2H+1<H^2 justify the padding. The parent H^70 mixed
third-order estimate remains valid; (5) is a sharper consequence at the
order actually needed here, not a correction to that estimate.

After the compact patch, the complete reference swirl has J3<H^2 by
[global_swirl_envelope](../global_swirl_envelope/README.md), §3. This
includes angular interpolation/edits, heat compensation, waiting, the
original terminal taper and genuine heat continuation. The axial taper,
both affine pulse end bumps and zero exterior have J3<H by
global_axial_source_envelope §2. Flat compatible joins and the unchanged
support cover therefore complete (3),(5) on EVERY positive native X.
They cover the parent's whole norm hull [Rmin^2/8,2Rmax^2], including
portions outside the active annulus; no spectral gap there is assumed.

For completeness, with the ORIGINAL Xlo=Rmin^2/8, separate whole-hull
primitive (value, X-derivative, eta-derivative) bounds are

    E0: (max(bK,H^4), max(bK[1+1/(2Xlo)],H^4/Xlo),
         max(bK/delta,H^4)),
    U0: (max(K,H^2), max(K,H^2/Xlo), max(K/delta,H^2)).

These supply the parent's ordinary first-jet input format explicitly.
The logarithmic route (6) gives a smaller directional estimate without
amplifying by the ratio Xhi/Xlo.

Care is needed with regional norms: max of regional SUMS of suprema need
not bound the global SUM of suprema. We use regional POINTWISE directional
bounds and then take their maximum. Equivalently, a global separate-jet
enclosure would take the maximum for EACH jet before summing it.

## 3. Directional estimates and existing reference gaps

Let Aphys=1/2+h and d=1-eta^2. The parent's exact native operator is

    D_h=2eta(Euler_q-X partial_X-Aphys)+d partial_eta.

The leading profiles have Euler_q E0=Euler_q U0=0. For any leading B,
|D_h B|<=2 (|B|+|B_y|+|B_eta|), since the absolute coefficients are
2Aphys|eta|<2, 2|eta|<=2, and d<=1. This pointwise dual-norm bound is
sharp using B=log(X/X0) at X0, eta=1. In the axis region the swirl
conversion allows the better bound

    M(a,c)=max_[0<=z<=1] [2az+c(1-z^2)]
          =2a if a>=c, otherwise c+a^2/c,
    DE=max{b K M(1+h+x,delta^-1), 2H^4},
    DU=max{K M(Aphys+x,delta^-1), 2H^2}.                 (6)

Then |D_h E0|<=DE, |D_h U0|<=DU on the whole hull. This uses primitive
field/loop bounds and fixed axis constants, not a supplied desired H_Z
bound. The actual heat exterior's exact cancellation D_h E0=D_h U0=0
from the parent is retained; replacing it by the larger bound (6) is safe.
No interior profile is replaced by heat and no pressure is reset.

Reference gaps are needed only at the ORIGINAL representatives in the
active annulus, as in core_physical_representative_closure. The same
[leading_frame_foundation](../leading_frame_foundation/README.md), (2)--(7),
supplies f_*, f^a, alpha, abar, bbar, nu and thus

    V=abar+bbar^2/alpha,
    g_*^2=f_*^2 alpha(2+nu),
    lambda_*^2=2 alpha f_*^2 nu/V,
    (lambda^*)^2=2 abar (f^a)^2.                         (7)

The analytic anchor and integrated log-swirl formula in that checkpoint
specify f_*,f^a; no diagnostic rational fixture replaces them. On the
1/4<=s<=4 chart, 1/8<s^(-1-h)<8. Explicit outward bounds without square
roots are therefore gmin=min(1,g_*^2/64),
lambda_min=min(1,lambda_*^2/64), lambda_max=max(1,64(lambda^*)^2).
For the original fixed u^2=2^19 G^22 use Ubound=1+u^2>=u.
These coarse choices change bounds, not the actual tilt, gaps or field.
Their enormous selected values are kept as symbolic positive expressions.

## 4. Same complete diagonal, same selector, and the remaining inequality

The normalized Euler prescription already in
[background_jet_foundation](../background_jet_foundation/README.md),
(9)--(11), and core_background_finite_head_budget (1)--(3), bounds the
COMPLETE tuple (E_n,U_n,W_n^r,J_n^r) with K_(m+1) cutoff jets and
requires C_(n,m) a_n^(-nh)<=2^-n for m<=n. Thus m=1 is included for
EVERY n>=1, including coefficient one. The whole-support geometric bound
used by the parent follows from that SAME prescription:

    |delta E,delta U|_(Euler,X,eta through1)
       <= q^h/(2-q^h) <= q^h.                           (8)

This identifies the analytical schedule condition; actual numerical a_n
were never materialized. An external schedule specified only by
unconverted frozen jets must still satisfy the stated conversion. We do
not retrofit new a_n. The full radial term
chi W_n^r-(Euler chi)J_n^r remains, although it is not an independent
leading input in this particular tangential floor.

Insert DE,DU from (6) into the parent's H0 and retain its H1, carrier,
rounding and box constants. Explicitly, with pplus,qz,ct,ca,cb as in the
parent, Kchart=4/(1-2h), Xhi=2Rmax^2,

    H0=Kchart[(pplus/Rmin)DE+qz DU],
    H1=4Kchart(pplus/Rmin+qz)(1+Aphys+Xhi),
    Aerr=ct+Lhigh H0, Berr=ca, Cerr=cb+Lhigh H1.

The existing core_primary_tail_recheck selector (2)--(3) ALREADY imposes
Cslow S^2(epsilon+epsilon^2+k^-1)<=tau, Cslow>=1, tau<=1, and
S>=Lreq/(2r0). It supplies the parent's slow gate directly on that tail.
The separate scalar ell>=512/h^2 bound need not be identified with this
selector. This does not identify the certified tail with an unknown old
first band or the previously localized force.

Put Sfloor=max(1,Lreq/(2r0)), using the EXISTING Lreq without modification.
The precise remaining sufficient comparison is

    4[Aerr/Sfloor+Berr/Sfloor^2+Cerr/Sfloor^3]^2 <= beta,
    beta=lambda_min/[4(1+Ubound^2)^2].                  (9)

Every coefficient is now a displayed primitive regional/reference/geometry
expression. All powers decrease with S, so (9) proves the normal floor
on the whole existing accepted tail if it passes. A failure of this
conservative comparison is not nonexistence of the actual normal floor.
The code reports failure and never increases Lreq or selects a later band.
Neither the inherited records nor this checkpoint proves actual substitution
in (9), including the actual Llow,Lhigh,Dbox and unchanged Lreq. Existence
of some sufficiently late half-line does not discharge it.

Next: compare these same-profile expressions against the existing selector's
Lreq and its target-freezing/root/cone budgets. Complete path-source norms,
wave-covariance axial/pressure mean updates and the full viscous, cutoff
and interaction residual remain separate obligations. No unconditional
pulse bound, complete residual improvement, infinite field or blow-up is
claimed.

## 5. Reproduction

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
The output must be fresh. Symbolic controls differentiate the actual
pullback, regular swirl conversion and native direction; exact controls
test the regional norm distinction, hierarchy padding, reference squared
floors and unchanged-selector pass/failure. They are analytical/algebraic
checks, not numerical measurements of the selected field.
