# Q76: all-order heat compensation on the first exterior power segment

Q75 stops at Xstar=Xrep. The next leading field is NOT the five-row
compact repair repeated at a new radius. It is the three-row **relative**
swirl correction driven by the terminal heat edit. We enclose its
coefficients at each angular order with the original h^4/Xstar factor,
then bound mixed leading jets from Xstar to the axial pulse entrance
Xp=Xstar exp(20). These are fixed-profile estimates under the retained
exact heat root and reference schedule, not an enclosure of every
physical pulse or of the ambient propagator.

## Differentiate the tail before normalizing the root

Use the [heat edit](../heat_exterior_audit/README.md) and the real
[Gamma derivative formula](../annular_exterior_foundation/README.md).
Beyond XK the pre-heat swirl Eref is angularly uniform and
|Eref(XK x)|<=2 eK x^(-1/2-h), x>=1. The multiplier is
1+chi(H(Z)-1), Z=2(1-eta²)/X, 0<=chi<=1 independent of eta.
Write ||v||_k=sum_(j<=k) sup|partial_eta^j v|/j!, summing from zero.
With 0<h<=1 and X>=XK>=1 put

    K_k=4+sum_(i=1)^k 6^i (i+1)!,    L_k=2K_k+K_k².

Then ||H-1||_k<=h K_k/X. Indeed the derivative bound
|H^(i)|<=(h)_i(1+h)_i gives |H^(i)|/i!<=h(i+1)! for i>=1.
The value follows by the mean value theorem since Z<=2/X.
For positive angular degrees compose the finite jet with
Delta Z=(-4eta z-2z²)/X. The positive polynomial (4z+2z²)^i
has coefficient sum 6^i, and X^(-i)<=X^-1. This gives a sum of
coefficientwise suprema, valid one-sided even at the poles. In
particular it assumes neither a convergent angular Taylor series
nor a complex evaluation of H at negative Z.
By multiplication the squared multiplier minus one has norm <=h L_k/X.

In the three ordinary terminal increments (new minus reference) integrate
these bounds in X. With the ORIGINAL X densities E²/(2X), -E²/2 and
sqrt(2X)E respectively, this gives

    ||Delta Cp||_k <= h L_k eK²/XK,
    ||Delta S ||_k <= 2h L_k eK²,
    ||Delta I ||_k <= 2 sqrt(2) K_k eK sqrt(XK).             (1)

The scalar integrals are int_1^infty x^(-3-2h), int x^(-2-2h)
and int x^(-1-h), respectively. In the last row h int x^(-1-h)=1;
discarding h or setting it to zero before integration invalidates
the bound. The majorants also justify differentiating the convergent
renormalized-I difference. The fixed axis pressure is common to both
profiles; (1) never bounds an incoming full pressure by a local edit.

## The same accepted heat root at high angular order

Set g=h^4/Xstar, f=(1+eta²)^(-1), e=e_star. The target divides the
signed increments (-Delta Cp,Delta S,-Delta I) by
(e² f², Xstar e² f², sqrt(2) Xstar^(3/2) e f). The retained schedule
gives eK sqrt(XK)/(e sqrt(Xstar))<=h^4 and XK>=Xstar>=1.
The angular algebra norms of f^-1 and f^-2 are <=5 and <=25
at EVERY finite order (they are polynomials). Thus (1) bounds the
target in component maximum angular norm by g Z_k, where

    Z_k=max(25 L_k,50 L_k,10 K_k).                           (2)

For example the pressure scale before f^-2 is <=h^9 L_k Xstar/XK²
and the energy scale <=2h^9 L_k/XK. The I row retains only h^4.
This is an upper bound at each order, not a contraction condition
imposed at every order.

The same constant matrix and quadratic map have inverse norm <=50
and bilinear norm <=100. The accepted real root has |c|<=6400g;
the frozen smallness test 8·50²·100·64g<1 gives the pointwise
Jacobian inverse <=100 by a Neumann series. Put
D_0=6400 and for j>=1 use

    D_j=100 (Z_j+100 sum_(l=1)^(j-1) D_l D_(j-l)).

Then sup|c^(j)|/j!<=g D_j. Differentiate the exact quadratic map;
the two endpoint terms belong to its Jacobian. In the remaining terms
use g²<=g. Existence of the fixed jets follows locally by that inverse
and the branch already selected (or by differentiation from the interior
and the one-sided smooth limits). This does not reselect any root or
infer complex angular analyticity of the heat compensation.

## Physical inputs through the reserved supports

Write t=log(X/Xstar), alpha=-1/2-lambda. By the original schedule
Tw=240T, Xstar=XR exp(241T-18), and
[Xp](../pulse_stress_audit/README.md)=XR exp(T+2+Tw).
For 0<=t<=20, U=0 and E=e f exp(alpha t)(1+sum c_j beta_j(t)).
Only the three disjoint bumps in the second reserved patch occur in
this LEADING field; the later positive-order coefficient supports
are not new leading bumps. Here beta=sigma'((t-center)/.3+1/2),
not the compact repair's width-.1 unit-mass bump. Put b_0=64,
b_i=162·3260^(i+1)((i+1)!)³ (.3)^(-i) for i>=1, using Q63.
Let F_k=2 sum_(j=0)^k 2^j, bounding ||f||_k by the Q68 disks.
For each r the raw log-radius derivative in the angular algebra obeys

    ||partial_t^r E||_k <= e F_k [ |alpha|^r +
       3g (sum_(j=0)^k D_j) sum_(i=0)^r binom(r,i)|alpha|^(r-i)b_i ].

This follows by differentiating the RELATIVE bump product. The second
summand bounds the correction alone. The real floor uses the original
small value bound (and actual disjoint supports), not this high-order
upper bound. Convert to physical X derivatives by Q70 (3) with
X0=Xstar. Thus the high leading background jets requested by Q60 have
explicit finite constants on this segment; this does not supply the
two-time propagator estimate. For Q58, w3=0 here; insert the above
bound and its next radial order in Q75 (1),(2) for the w4 change
(linear in the swirl correction). Its fixed-order size is O_k,r(g)
with the original geometric factors and e. The diagonal weights are
unchanged and the pressure-cutoff moment Cchi still moves with s.
On the remainder of this interval the shear weight uses the unedited
power alone. No new power of the physical band epsilon is obtained.

For the state/target carry the genuine five matched traces at Xrep
from Q68--Q74 and integrate Q67's densities through this patch.
In particular Pi(X)=Pi(Xstar)+int_(Xstar)^X E²/(2x) dx, with the
edited E inside the integral; U=0 does not suppress that contribution
or set M,J to zero. The three heat equations restore the totals only
in combination with the FUTURE terminal edit. They do not assert
that the compensation-only increments vanish at t=5.

Next enclose the axial pulse and the subsequent exterior stages to Xb
on the same physical hull, then verify the inverse and propagator gates.
The physical axial/radial tail estimates, uniform pulse bounds, complete
residual improvement and blow-up remain unproved.

Run test_audit.py and audit.py --verify-record evidence.json --output a
fresh path. The algebra checks and replay do not prove continuum existence.
