# Actual logarithmic Haar pullback and projector seminorm boundary

Parent: `e73aac0d9a14ad33eb171db92fbef6c203c576ae`.
On 24 October 2026 its complete initial Haar CI step passed; the full
run remained in progress without failures. This checkpoint derives the
missing R-to-logarithmic coefficient pullback from the ORIGINAL scale
chart and routes the actual projector's fixed-family seminorms. It
does not select a replacement bump or claim its uniform jets follow
from weighted mass normalization. Historical evidence, physical bands,
widths, supports, gauges and roots remain unchanged.

## 1. Finite jets of the actual coordinate map

Use the retained branch s-Z^2 s^(2h)=T, 1/4<=s<=4, |eta|<=1.
It has NO R dependence. Thus in x=log(X/Xa), X=R^2/(2s),

    R(x,Z,T)=sqrt(2Xa) exp(x/2) sqrt(s(Z,T)).

At a real center let log(s/s0) have nonconstant normalized slow jet
sum bounded by L(t), with exact-degree coefficients

    L_j=(1/3)sum_(i=0)^j 32^i 64^(j-i), j>=1.

These are the existing scale-chart bounds, not new analyticity
assumptions on a physical profile. Let

    E_n(b)=[t^n] exp(b[t+L(t)]), E_0(b)=1,
    E_n(b)=(b/n)sum_(j=1)^n j(delta_j1+L_j)E_(n-j)(b),
    C_j=2 delta_j1+Rmax E_j(1/2), j>=1, C_0=0,
    T_n=sum_(k=0)^n [t^n] C(t)^k.                     (1)

The two identity coordinates are Z,T. Haar averaging has already
removed torus dependence, so this is a THREE-variable composition,
not a reused five-variable U/R bound. Rmax is the SAME whole radial
norm-hull upper bound. Higher coordinate derivatives do not vanish.

If the normalized through-order R,Z,T norm of a scalar F is <=K,
its outer exact-degree Taylor coefficient sum is <=K at each degree.
Finite multivariate Taylor composition therefore bounds the output
exact-degree sum by K T_n. Only derivatives through n are requested.
The coefficient multiplier R^(e+1)/2 has normalized jets bounded by
Rmax^(e+1) E_i((e+1)/2)/2. Multiplication and summing degrees through n
give the explicit required constant

    PB_e,n=Rmax^(e+1)/2
          sum_(i+j<=n) E_i((e+1)/2) T_j, e=1,2.      (2)

This proves the parent's fixed PB inequality for the actual nonlinear
map on its retained common hull. All constants are finite, positive
and band independent at fixed order. The physical X,eta and hence
zeta,delta at corresponding points are identical. Derivatives of F
already include derivatives of its flat factors: no quotient by zeta
is differentiated and no inverse edge exponential is introduced.
Thus the row's epsilon order and pole2p remain unchanged. Physical
endpoint derivatives use the existing one-sided/zero-extension
premises. Neither a new tail nor a new radial hull is selected.

## 2. Route the SAME moving projector family, without inventing its jets

Let Sigma_e,n(S) bound the actual fixed family's normalized through-n
R,Z,T derivatives, uniformly over the accepted bands. Its physical
support is on the SAME declared moving interior hull x_phys=r/sqrt(q)
in [l,u], strictly inside the active annulus. Reuse the existing
reserved-geometry proof for that entire support to get zeta>=exp(-E).
Every differentiated bump remains on that moving support. Then

    B_e,n=exp(E) PB_e,n Sigma_e,n                      (3)

bounds the differentiated pulled-back projector (R^(e+1)/2)sigma_e
with full zeta weight and pole0. Exp(E) stays symbolic; this is not a
rectangular-hull weight floor. The code obtains PB from (2), uses the
same primary width, and substitutes (3) into the parent's projected
source recipe. It keeps BOTH moment channels and the WHOLE moving
c_rho P contribution. There is no new independent current-row,
projected-source or moment callback.

Sigma_e,n is a fixed-profile/family seminorm, NOT an already bounded
current PDE field. A derivation for a specifically recorded moving
sigma_e must supply it. The repository's abstract normalization and
support statements alone are insufficient to infer a band-uniform
seminorm. Formula (3) is conditional on this named lower input, not
an assertion that it has been evaluated or that a convenient standard
scaling is the inherited actual choice.

## 3. Precise obstruction to a normalization-only inference

Take any nonzero nonnegative C-infinity Phi compactly supported in
(1,2). Put m_e=int R^(e+1)Phi/int R^e Phi in (1,2), A_e=int R^e Phi,
and consider the AUXILIARY family

    sigma_e,N=Phi(R)/A_e
                  [1+(R-m_e)sin(NZ)/2], N=1,2,... .   (4)

It has the SAME support, is nonnegative (the bracket is >=1/2), and
int R^e sigma_e,N=1 exactly for every Z,N. Yet at Z=0,

    partial_Z sigma_e,N=N Phi(R)(R-m_e)/(2A_e),

whose supremum is N times a fixed positive number. Hence neither
mass1, positivity, fixed interior support nor smoothness at each band
proves a UNIFORM moving derivative estimate. The retained logarithmic chart keeps Z unchanged; at Z=0 its
R_Z is zero because s_Z=0, so the same derivative growth persists. This is an obstruction to
an inference, not an actual canonical projector or a construction
counterexample; (4) is never substituted for the inherited sigma_e.

The remaining exact inequality is

    sup_accepted_bands |sigma_e|_(R,Z,T,through n)
        <= Sigma_e,n(S),

for the actual chosen fixed-family formula and support. A scale-covariant
formula with fixed smooth shape would permit the existing composition
recipes to supply such constants, but it must be identified as the
actual choice rather than silently imposed.

## 4. Validation and next obligation

Independent finite expansions check (1), mixed three-variable chain
controls check (2), and distinct angular/axial multipliers are retained.
An explicit normalized polynomial bump checks every identity in (4)
at finite regularity; the continuum obstruction above uses a genuine
C-infinity bump. Passing these manufactured controls is not a continuum
bound on the selected giant field. Historical frozen evidence is preserved.

Next identify/bound the actual projector family jets on its original
moving hull, then apply the existing two-edge primitive, reciprocal
weight and native composition to these source coefficients and join the
ORIGINAL signed root. The implementation continues to mark the signed
primitive/root uncertified. Full background/phase/value/support gates,
global pressure moments, later complex covariances, all-scale realization
and complete viscosity/cutoff/interaction residual remain original
obligations. No unconditional pulse, residual improvement, infinite
realization or blow-up result is claimed.
