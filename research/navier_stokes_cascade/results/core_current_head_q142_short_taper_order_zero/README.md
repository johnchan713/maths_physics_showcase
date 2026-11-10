# Q142: selected order-zero short-taper caps and the regional obstruction

Parent: Q141, `3801ea92a4d114236b93a6084576a578a00910e4`.
The live branch, Q141 evidence, Q72, the prescribed continuation, the
axis/core attachment and current CI were read first. Historical evidence,
parameters, supports, roots, pressure gauge, frequency and the width
`t=C^-20` are unchanged.

Q141 left its order-zero inputs open. They can be supplied analytically on
the same short taper. Put `Xa=4/Lambda`, `0<=y<=2t`, and retain all five
Q72 incoming traces. Since `C>=2`, `exp(2t)<2`, so

    X<=8/Lambda,  X^-1<=Lambda/4,
    sqrt(2X)<=4,  1/sqrt(2X)<=Lambda.                    (1)

With `F0=2S`, `U0=2(9+S/Lambda)`, the common trace cap is the maximum of

    Xa U0, Xa^2 F0, Xa^2 F0 U0,
    Xa U0^2+Xa^2 F0^2/2, 20P^2+Xa F0^2.                (2)

The last entry is the full pressure trace. The `I,J` rows enter `p1r`; the
`Pi,Pi_eta` rows enter `nsr` and hence the actual axial field and `p2`.

For the reference logarithmic increment, integration by parts against the
monotone taper gives `|Jr|<=3(log(10S)+2)`. Thus
`exp(+-Jr)<729000 S^3`. For the actual field,

    Ja=-1/2 integral kappa p1r,

and the unchanged axis/core attachment proves `0<p1r<B/100` on the whole
continuation. The selected hierarchy has `tB<1`; hence
`exp(Ja)<=1` and `exp(-Ja)<3`. These are value caps. No derivative seminorm
is exponentiated.

The same attachment also proves `p1r>=X/10`, while the selected step has
`kappa>=t`. Therefore the actual short-taper shear has the explicit floor

    a=kappa p1r >= 2t/(5 Lambda)>0.                       (3)

This corrects the apparent floor gap in Q141. Equations (1)--(3), the five
traces and the two exponential caps feed the unchanged Q141 recurrence.
Astronomical caps such as `C exp(2B)` are recorded by integer binary
exponents; materializing them is neither necessary nor a proof.

The short taper is no longer the first missing order-zero premise. The next
global step needs one identity-bound finite cover using a common convention.
For every Q69--Q71/Q73/Q127 region `R`, it must establish

    H_R,Ay_R,Be_R<infinity, a_R>=a_*R>0,
    G_R>=max(2,|mu|,|p2|,|mu p2|), V0_R<=N,              (4)

with the unchanged integer `N`. Only then may one take the maximum of the
upper bounds and minimum of the floors, apply Q140--Q138, and join the later
full-pressure repair to the regional cone estimates. Existing regional
notes use different inputs and do not yet supply (4) as one actual-profile
record. That is a precise unproved inequality package, not evidence that it
fails.

The tests check exact rational caps, guard failures, the unchanged Q141 map,
binary representation and finite-cover max/min logic. They are diagnostics,
not a continuum proof. Uniform physical pulses, complete wave-covariance
mean corrections, coefficient and reserved-mean supports, the full residual
and blow-up remain unproved.
