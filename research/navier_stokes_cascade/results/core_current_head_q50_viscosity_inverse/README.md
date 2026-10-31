# Q50: a quantitative necessary target estimate for strict viscosity gain

At Q49 the original axial remainder already fits the first-real row budget.
The [earlier quotient calculation](../core_initial_quotient_fast_inverse/README.md)
also extracts local pressure -Wrr and gives order two under its hypotheses.
Neither original-tail result estimates the later postwave tail: its distinct
potential includes -2 G_R Psi and full projected pressure
([axial pressure audit, (6)](../core_axial_potential_pressure_gain/README.md)).
Here we address the equality in the first-real fast viscosity difference.

The [actual-profile gate](../core_initial_actual_profile_gate/README.md)
proved injectivity modulo constants for the principal cross profile, but
nonvanishing alone says nothing about the size of the viscous term along
the accepted bands. We derive a reverse estimate uniform in the band and
in the size of the transverse zero set. This tests a proposed strict gain
without assuming a numerical lower bound for the physical target.

## Recovering the target from the leading viscosity

Freeze a slow point. Let Lambda denote the COMPLETE selected principal cross
in one of the two target channels, Xi its normalized auxiliary Haar mean, and
G=L^3 Nt^-1 Qcirc Lambda on the unit two-torus, with the same frozen directions
L=(1,-b).grad, Nt=(b,1).grad, b=sqrt(2)-1. Suppose Lambda has a zero on this
common torus. This is the transverse support premise for the summed profile
in the actual-profile gate. Zeros of different summands at different points
would not suffice. No fixed common zero across slow points or bands is needed.

For k nonzero, write l=k1-b k2 and t=b k1+k2. The inverse symbol has
modulus |t|/((2 pi)^2 |l|^3). In the shell K<=|k|_1<2K, K dyadic,
the quadratic norm identity l t=b[(k1+k2)^2-2k2^2] gives |l_h|>1/(3|h|_1)
for every nonzero integer h. Distinct l values in this shell, with zero
adjoined, thus have separation >1/(12K). Ordering each sign by distance
from zero gives sum |l|^-6 <=4(12K)^6. As |t|<=2K, the square-sum norm
of the inverse symbols on that shell is <=4*12^3 K^4 (discarding favorable
2 pi factors). Split the shell by the larger coordinate |ki|>=K/2.
Parseval for partial_yi^5 G bounds the square-sum norm of G's shell
coefficients by sqrt(2)(2/K)^5 J5, where J5=max_i sup|partial_yi^5 G|.
Cauchy--Schwarz leaves K^-1. Summing the geometric series gives

    ||Lambda-Xi||_infty <= C J5,    C=2*(4*12^3)*2*2^5.
    |Xi| <= C J5,                   ||Lambda||_infty <= 2 C J5.

The second line evaluates the first at a zero; the first line holds even
without support information. Constants are for raw auxiliary derivatives;
a normalized through-five bound costs at most 5! to recover J5. A weight
depending only on slow variables factors out at the same point. No slow
derivative of a moving zero location has been taken.

In the gate's notation the magnitude of the leading viscous component is
that of V=epsilon c_i^-1 M^3 j^3 G, with j=d R^(d-1) positive on the
original hull (an overall sign does not matter). All prefactors are constant
in y. Write j>=j0>0, M=c_i^rho epsilon^-kappa, 0<rho<1, and
1<=1/c_i<=Tg S. Since c_i^(1-3rho)<=(Tg S)^2, the necessary bound is

    |Xi| <= C j0^-3 (Tg S)^2 epsilon^(-1+3kappa) J5(V).

Thus a weighted bound on V at nu+eta, nu=3-6kappa, requires a
same-weight VALUE bound on Xi at alpha+eta, alpha=2-3kappa, up to a
polynomial in S. This is a converse to the historical sufficient route
using improved target jets. A value bound alone does not supply those jets.
If an estimate concerns the complete fast difference, first subtract
the other terms using their separately proved five-jet bounds; if those
have margin sigma, the necessary gain is min(eta,sigma). If another term
of the same order is present, the inequality tests only the isolated
component. No claim is made about cancellation between different residual
blocks. Neither target channel has been proved to attain its upper envelope.

## Why the jet order matters

A pointwise estimate with no derivative cost does not follow uniformly
from compact transverse support and nonzero mean alone. Take smooth
one-periodic h supported strictly inside one period, of mean one, and
Lambda_k(y)=h(k.y). Each profile has an open zero set and mean one.
On this subspace G=(l^3/t) h''(k.y), exactly. For the Pell sequence
a+b' sqrt(2)=(1+sqrt(2))^j, take k=(a-b',b'). Then |lt|=sqrt(2)-1,
and |t| is comparable to |k|_1 as j grows. Hence J_s(G)=O(|k|_1^(s-4))
for each fixed s; even the through-three norm tends to zero. These are
comparison profiles, not claims about physical pulse frequencies. We
do not resolve the endpoint at four derivatives or assert optimal constants.

Next test the value of the actual Q36 signed target, equivalently its
projected source primitive, against this necessary rate on the common
support. A matching target-jet upper bound or a quantitative sequence
violating the value bound would decide different sides of this obstruction.
The latter must beat the allowed polynomial and use the actual profiles.
Physical field enclosures, remaining weighted postslow rows and the full
viscous, cutoff and interaction residual still need verification.
