# Q131: angular swirl modulation and all five full-memory discrepancies

Parent: `5af6852b597b632db59e6f8f57fa692131f1b36d` (Q130).
Fresh fetch on 3 November 2026 agreed with this head. Q130 research CI
41699897115 and field CI 41699897210 were in progress; Q129 field CI
41689746383 had passed. Current notes, Q130 frozen evidence, Q129, Q64,
the all-order modulation map and the original loop bounds were read first.
Historical evidence, selected roots, phase origin, N, physical scales,
supports, pressure datum and C_setup^-20 width remain unchanged. The setup
activation result still concerns only the leading core/first collar.

Q128--Q130 give conditional closed angular bounds for the two SAME
zero-mean loop primitives. Here we propagate them through the exponential
swirl update, its inverse, the axial update, and the five primitive
DIFFERENCES before the later exact repair. This completes an angular
modulation map estimate, not an actual-profile certificate or the missing
final compact/cone estimate. Positive floors, full PREMOD caps, real loop
caps and the original value/frequency contract remain premises.

## 1. Retain the real value gate, not an exponentiated high norm

Hold the slow radius y fixed; the fast evaluation zeta=Ny is independent
of eta because N is the ORIGINAL selected fixed integer. Write f=Acal,
g=Bcal after the original moving-phase pullback. Q129/Q130 supply

    |f^(m)|<=MA LA^m(m!)^3, |g^(m)|<=MB LB^m(m!)^3,
    LA,LB>=1, m>=0.                              (1)

Separately retain the original sharper REAL value bound

    |f|,|g|<=V0<=N.                              (2)

The original |Acal|,|Bcal|<=H and N=1+floor(H^32), H>=1, imply (2):
N>H^32>=H. This is a check on the retained frequency/value contract, not
a new selector. Its actual premises must still hold on the whole interval.
Since e<3, (2) gives exp(|f|/N)<3. We never replace this by exp(MA/N)
or exp of a high derivative norm. The latter could be enormous without
invalidating the SAME value gate. No huge exponential is evaluated.

At a real eta center let t=MA/N and Lexp=LA(1+t). Normalized positive
Taylor increments of f/N have coefficients at most t LA^k(k!)^2.
Finite exponential composition gives, for m>=1,

    |(exp(f/N))^(m)|
      <=3 LA^m(m!)^3 sum_(j=1)^m t^j binom(m-1,j-1)/j!
      <=3t Lexp^m(m!)^3.                         (3)

Each product of positive index factorials is at most m!, and 1/j!<=1;
the last sum is <=t(1+t)^(m-1). At m=0, use exp(|f|/N)<3.
The same estimates hold for exp(-f/N), without division by f.
For the DIFFERENCE, |exp(+-f/N)-1|<=3V0/N at m=0, so with

    D=3 max(V0,MA)/N,
    |partial_eta^m(exp(+-f/N)-1)|
                       <=D Lexp^m(m!)^3, m>=0.  (4)

This retains a 1/N prefactor, while Lexp itself depends on MA/N.
It is not a uniform estimate in arbitrary frequency-dependent inputs.
All steps use finite real jets; there is no analytic eta continuation
through the flat variance cutoff.

## 2. Original field and inverse update, at every angular order

Supply uniform complex PREMOD caps e,i,u for E,E^-1,U on the same eta
disks of radius 0<rho<=1. These are actual same-field caps, not a
consequence of knowing a,b,p alone. A smaller common radius from Q128 is
permitted. Set R=rho^-1 and

    B=max(R+Lexp,LB),
    DE=eD, DI=iD, DU=MB/N.                       (5)

For the unchanged update E_N=E exp(f/N), U_N=U+g/N,

    |partial_eta^m Delta E|<=DE B^m(m!)^3,
    |partial_eta^m Delta(E^-1)|<=DI B^m(m!)^3,
    |partial_eta^m Delta U|<=DU B^m(m!)^3.        (6)

Cauchy bounds the PREMOD normalized coefficients by eR^k, iR^k, uR^k.
For a product, distribute the normalized derivatives and use
sum_(k=0)^m x^k z^(m-k)<=(x+z)^m, with k! (m-k)!<=m!.
Applying this to (3),(4) yields both (6) and

    |partial_eta^m E_N|<=3e max(1,t) B^m(m!)^3,
    |partial_eta^m E_N^-1|<=3i max(1,t) B^m(m!)^3,
    |partial_eta^m U_N|<=(u+DU) B^m(m!)^3.       (7)

The real value bounds remain E_N<=3e, |E_N^-1|<=3i,
|U_N|<=u+V0/N. If E has the original positive real floor e_*,
then E_N>=e_*/3. Positivity does not require differentiating the
exponential or accepting a new root. These coarse bounds do not replace
the sharper inherited low-order constants in the repair acceptance gate.

## 3. Five moment differences with the genuine incoming memory

Fix a real interval [x0,x1], x0>0, with endpoints independent of eta.
The SAME five X densities from Q67 are

    U, sqrt(2X)E, sqrt(2X)EU, U^2-E^2/2, E^2/(2X).

With r=max(1,2x1)>=sqrt(2x1), their modulated-minus-PREMOD differences
have raw angular bounds D_row*(2B)^m(m!)^3, where

    D_M=DU,
    D_I=r DE,
    D_J=r(e DU+u DE+DE DU),
    D_S=2u DU+DU^2+e DE+DE^2/2,
    D_Pi=(e DE+DE^2/2)/x0.                       (8)

These follow by expanding the actual signed densities before taking
absolute values. In particular D_Pi bounds a pressure-integrand difference,
not a velocity plateau or a terminal pressure reconstruction. Two-factor
products cost base 2B, also for the quadratic terms. The 1/N factors in
DE,DU remain explicit, including their quadratic interactions.

For EACH row supply its actual incoming discrepancy at x0 with raw bound
T_row L_row^m(m!)^3, T_row>=0,L_row>=1. All five are mandatory inputs.
With ell=x1-x0 and Lmem=max(2B,all L_row), the fundamental theorem gives

    |partial_eta^m Delta P_row(X)|
           <=(T_row+ell D_row)Lmem^m(m!)^3.      (9)

Integration endpoints are fixed in eta. All intervening supports and the
incoming pressure datum are included. Zero trace discrepancies may be
entered only when equality on the whole preceding prefix and the SAME
pressure gauge has been proved. Vanishing velocity modulation at x0 or a
zero value at one eta does not establish that equality. In particular no
incoming Pi is silently reset.

The forward pressure contains BOTH 4Aphys eta Pi and -(1-eta^2)Pi_eta.
To retain the latter, a further angular derivative of (9) is bounded by

    Dmem Lmem^(m+1)((m+1)!)^3
              <=Dmem Lmem (8Lmem)^m(m!)^3,      (10)

since m+1<=2^m for integer m>=0. At a finite request m this needs the
primitive inputs through m+1; no derivative has disappeared. For an
incoming perturbation c(eta)=eta, the original forward pressure changes
by -X/E at eta=0 despite c(0)=0. Thus a value-only incoming discrepancy
is insufficient for (10). Equations (9),(10) preserve the full pressure
memory needed by the eventual five-row repair.

These are pure angular estimates at fixed y. Radial differentiation of
f(y,eta,Ny)/N retains the original N^(r-1) field cost from Q64; it is
NOT generally O(1/N). The manufactured radial countercontrol checks this
separation. No physical support or pulse width has changed.

## 4. Checks and next obligation

Independent exponential-series coefficients check (3),(4), including the
m=0 difference, signed inverse exponent and the retained value gate.
Explicit modulated E/U functions check raw angular derivatives, including
E_N^-1. Direct five-density expansions and polynomial/exponential memory
controls check (8),(9); nonzero incoming pressure and its derivative are
retained. Shifted factorial controls check (10), and a radial oscillatory
example rejects an O(1/N) radial reinterpretation. Moderate exact bounds
are enclosed by binary budgets at order 100005 without factorial expansion.
These are manufactured diagnostics, not a selected-parameter evaluation or
continuum proof by testing. The audit replays Q130 and checks the frozen
source/evidence hash closure; fresh output must match byte for byte.

Next bind the actual PREMOD complex caps, real value/frequency gate and
Q128--Q130 premises over every original active region, then propagate
these SAME full-pressure discrepancies through the LATER accepted exact
five-row repair, retaining its original smallness/positive-field gates.
The final compact norm and cone join, global coefficient/reserved-mean
supports, uniform physical pulse estimates, complete wave-covariance
axial/pressure mean corrections and full viscous/cutoff/all-interaction
residual remain open. No unconditional pulse bound, complete residual
improvement or blow-up follows. Historical frozen evidence is unchanged.
