# Q68: angular traces on the native outer prefix

Q67 took incoming moment derivatives as inputs on each fixed interval.
Here those inputs and the pure angular state derivatives are enclosed on
the **pre-modulation** prefix from X_R to Xrep, relative to the retained
exact earlier annular match and native schedule. This does not bound the
earlier joining root or any later modulated field by the same formulas.

## Start at the matched trace, not at zero

Write R=X_R, y=log(X/R), P=Pstar, f(z)=1/(1+z²). The right end of the
[annular join](../axis_matching_audit/README.md) matches all five ideal
moments at y=-5. The intervening ideal segment agrees in the densities,
so at y=0 the incoming traces are exactly

    M=4R eta,  I=(5 sqrt(2)/8) P R^(3/2) f,  J=4 eta I,
    S=R(16 eta²-(5/12)P² f²),  Pi=Pi_axis+(5/2)P² f².       (1)

In the [axis datum](../core_growth_activation/README.md), subtract the
known atom (5/2)delta_1 from the pressure measure mu. The remainder nu
is nonnegative with mass <5/2, and (1) gives exactly
Pi(R)=-P² integral f^(2 theta) dnu. This is cancellation of a known
primitive against a known component of the datum, not a pressure reset.
It uses no future exact repair or backward terminal identity.

For each real center in [-1,1] use the closed complex eta disk of radius
1/2. There |z|<=3/2, Re(1+z²)>=3/4, |1+z²|<=13/4. The logarithm
with real values on the real axis is single valued and has modulus <4
there (bound real and imaginary parts separately). For 0<=theta<=1,
|f^theta|<=2, |f^(2 theta)|<=4, |1/f^theta|<=4.
The measure integral is holomorphic there by local domination and its
modulus at R is <=10P². This remains valid at real eta=±1; it is an
extension of the explicit native formulas, not of the Gamma heat tail.

## Explicit majorants on every fixed prefix [R,B], R<=B<=Xrep

Use the [native outer schedule](../annular_exterior_foundation/README.md)
before the moving loop and compact repair: E=c(y) f^theta(y), U=k(y)eta.
The [pre-modulation shear audit](../premodulation_mixed_shear/README.md)
bounds |(log E)_y| by 1000 on this prefix, |theta_y|<=64/Tf<=1
(theta is constant off the interpolation), |k'|<=4 and |k|<=4.
At eta=0 this bounds |(log c)'|, with c(0)=P. Thus with Y=log(B/R),
set C=P exp(1000Y), F=exp(1000Y)/P. Then c<=C, 1/c<=F, uniformly
on the interval. These deliberately coarse positive bounds involve no
unevaluated high seminorm. On the complex disks |E|<=2C, |1/E|<=4F,
|U|<=6. They hold across the fixed radial seams because the native
fields there have the same smooth continuation; no complex y is used.

Define H=(5 sqrt(2)/4)P R^(3/2), q=B-R. Disk suprema of all five
primitives are bounded by the following numbers:

    Mbar=6R+6q,                       Ibar=H+2C sqrt(2B) q,
    Jbar=6H+12C sqrt(2B) q,
    Sbar=R(36+5P²/3)+(36+2C²)q,     Pibar=10P²+2C² Y.     (2)

Indeed (1) gives the first terms using |f|<=2. Integrate Q67's
original densities U, sqrt(2X)E, sqrt(2X)EU, U²-E²/2 and E²/(2X)
forward with these same inputs. The last integral costs log(B/R),
not B-R. The endpoints are independent of eta, and the uniform disk
bounds on finite intervals justify holomorphic parameter integration.
For every j, |partial_eta^j moment|/j! <= 2^j times its entry in (2).
In particular Q67's incoming trace inputs at ANY cut in this prefix
are now explicit to order l+1 without sampling or using final-state
moment identities. One may set B to that cut for a tighter bound.

## Pure angular state outputs and the remaining radial requirement

For the state work on disks of radius 1/4 about the same real points.
Each first angular derivative of a primitive is <=4 times its bound
in (2), by Cauchy's formula on the remaining radius. Write D=1/2-h,
A=1/2+h, d=1-z², L=1-2h z²; |d|<=13/4 and 1/|L|<=2 for
0<h<=.01. For the exact w=X-2D z M-d M_eta use
W=B+(3D+13)Mbar. Q67's original p2 numerator gives

    N=6W+7D Mbar+(6h+13)Sbar+B(6A+13)Pibar,
    |partial_eta^j p2|/j! <= 4^j 8F N.                     (3)

This uses both X(4A z Pi) and -Xd Pi_eta. The native shear identities
a=1-2[(log c)'+theta_y log f], b_s=2 k' z/E give, on the same disks,

    |partial_eta^j a|/j! <= 4^j (1+2(1000+4)),
    |partial_eta^j b_s|/j! <= 4^j 48F.                     (4)

Equations (2)--(4) bound actual pre-modulation inputs, conditional only
on the stated inherited join and schedule; they need not be small.
They close the pure angular part of Q67 on this prefix. They do not
bound the high mixed radial derivatives of the radial ramps: those
still require quantitative derivatives of c,theta,k and their physical
widths. Earlier inner/shape/annular inputs and the normalized repair
defect must also be enclosed. Uniform actual pulses, the tail rates,
complete residual improvement and blow-up remain unproved.

Run test_audit.py and audit.py --verify-record evidence.json --output a fresh
path. Exact tests check the trace, Cauchy margins and forward integrals.
They are not continuum verification of the underlying existence lemmas.
