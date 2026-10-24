# Explicit rho-squared projector realization and uniform finite jets

Parent `b1a956753cf3e96203d151537466b3b153362d42`, checked against the remote
branch on 24 October 2026. Its first-real refresh CI step passed; run
40095916664 remained in progress without reported failure. Historical
frozen evidence and the existing signed-target objects are unchanged.

The inherited abstract sigma_e definitions give normalized interior bumps,
but no identified scale law. Earlier work proved those properties alone
insufficient for uniform derivatives. This checkpoint supplies an explicit
NEW admissible choice, using ONLY the retained fixed pressure profile.
It does not claim the historical unspecified sigma_e already equals this
choice, nor that the old signed target has been recomputed with it.

## 1. Exact nonnegative projector without a new cutoff profile

The original compensator has rho(R,Z,T)=s^-1/2 rhohat(R/sqrt(s)),
integral rhohat=1, with one fixed smooth compact profile supported in
0<l<=x<=u. Its smoothness, exact mass and support identifications remain
the retained premises. rhohat may be signed. Define the FIXED number

    A2=integral rhohat(x)^2 dx,

and explicitly choose, for e=1,2,

    sigma_e(R,Z,T)=R^-e sqrt(s) rho(R,Z,T)^2/A2
                 =s^(-(e+1)/2) f_e(R/sqrt(s)),
    f_e(x)=x^-e rhohat(x)^2/A2.                         (1)

Cauchy-Schwarz gives 1=(integral rhohat)^2 <=(u-l) A2.
Thus A2>0 and 1/A2<=u-l even for signed rhohat. No division by the
possibly zero first/second SIGNED moment of rhohat is used. The exact
change of variable gives

    integral R^e sigma_e dR
       =sqrt(s)/A2 integral rho^2 dR=1.                (2)

Both projectors are nonnegative. Their supports are contained in the
SAME actual moving pressure-profile hull [sqrt(s)l,sqrt(s)u], strictly
inside the retained annulus. Squaring and multiplying by powers of a
positive R preserve smooth zero extension. No new physical cutoff,
profile edge, flat width C_setup^-20, graph tilt or accepted band is
introduced. This is nonetheless a NEW projector CHOICE: its consequent
signed target and root cannot be copied from a different sigma family.

The associated Q_e f=f-sigma_e integral R^e f has exact zero weighted
moment and satisfies Q_e^2=Q_e by (2). These are pointwise slow-variable
identities; the bumps remain moving functions and their derivatives
cannot be discarded from later commutators.

## 2. Uniform native finite-jet coefficient from fixed rho data

Let H_n be the parent's derived UNWEIGHTED full factorial-normalized
R,Z,T through norm bound for rho. It is composed from the finite FIXED
rhohat seminorms and the actual implicit-coordinate estimate. Its value
is independent of physical band index. No signed target/covariance norm
is used to derive H_n; the rho container is only a geometry/profile
interface. The construction is therefore not circular even when that
container is associated with a preexisting signed-state object.

The original coordinate gives s in [1/4,4] and the log-s jet majorants
L_j. Both sqrt(s) and s^-1/2 have normalized total-degree coefficients
bounded by U_0=2,

    U_n=(1/(2n)) sum_(j=1)^n j L_j U_(n-j).

Let S_n=sum_(j=0)^n U_j. On the same positive whole radial norm hull
[a,b], the exact reciprocal-power through bound is

    B_e(n)=sum_(j=0)^n binom(e+j-1,j) a^(-e-j).

The full normalized through norm is submultiplicative. With A2 fixed,
the ordinary Leibniz rule applied to (1) therefore proves

    ||sigma_e||_(n,unweighted)
       <= (u-l) B_e(n) S_n H_n^2.                     (3)

This includes ALL radial/slow mixed derivatives, including derivatives
of sqrt(s) and of the moving rho. There is no derivative loss, no
new edge pole and no physical-scale factor. The constant grows with n
and the retained fixed-profile seminorms; uniformity in n is NOT claimed.
Only finite smooth jets are used, not analyticity of a compact bump.

On the ACTUAL moving support the previous geometric estimate gives
zeta>=exp(-E). Multiply (3) by exp(E) to obtain a weighted sigma_e
through norm. This floor is not asserted on a rectangular product hull.
Multiplying once by the existing logarithmic Haar pullback coefficient
PB_e(n) gives a bound for (R^(e+1)/2)sigma_e in the logarithmic source
chart. The code exposes native, weighted and logarithmic_bump separately
to prevent inserting the floor twice. The enormous floor stays symbolic.

## 3. What this closes and what still needs recomputation

For the choice (1), the formerly independent uniform moving projector
seminorm inequality is now reduced to the already stated FIXED rhohat
seminorms and actual coordinate/support premises. No normalization-only
argument is used. The old high-frequency counterfamily still disproves
such a claim for an arbitrary inherited sigma family; it does not apply
to this fixed-profile realization with (3).

RhoSquaredProjectors is a bound/definition object for (1). It deliberately
does not mutate the current ProjectedSource, Xi_e, SignedTarget or frozen
evidence, and does not assert identity with their historical sigma_e.
To use (1) in the continuing construction, instantiate the projected
current sources with its native norm and explicit field definition,
then RECOMPUTE Xi_e=-T_e(E_e-sigma_e M_e), the signed root/curl, complete
covariance and pressure/scalar updates on that SAME choice. The previous
conditional algebra applies once its required inputs hold; its numeric
or symbolic coefficients cannot simply be declared unchanged.

The first-real refresh still requests new-layer target jets through
400009 at output zero, in addition to the nested old requests. Formula
(3) supplies a finite recipe at every fixed requested order when the
original finite rhohat profile bounds are available; this checkpoint
does not evaluate those enormous jets or certify the original fixed
profile enclosures. Whole pressure memory, all six covariance entries,
actual compact A0/A1 maps, the full reserved potential hull and all
interactions remain obligations of that recomputation.

Next connect this explicit projector realization to a freshly recomputed
first signed target, preserving its identity through the full derivative
and pressure domains. Actual mean-block closure, unconditional actual
pulses, complete residual improvement and blow-up remain unproved.

## Evidence

Exact controls integrate a signed endpoint-flat polynomial diagnostic,
its square and its normalized weighted projectors at three physical
dilations. They check nonnegative normalization, idempotence, mixed jets
of the actual implicit coordinate, reciprocal-power multipliers and
single floor/pullback routing. The polynomial has only the finite endpoint
regularity requested in those tests and is not the retained smooth
rhohat. The analytical argument above supplies the stated all-finite-order
conclusion under the original smooth-profile premises. Historical frozen
records are bound and replayed unchanged; tests are not continuum proofs.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
