# Q143: value-based mixed recurrence after the reference taper

Parent: Q142, `71424aa12fa04d125f835ff9ac63ff5845de2810`.
The live branch, Q142 evidence and CI, Q69--Q73, Q127, the prescribed
continuation and axis/core attachment were read first. Historical evidence,
width, roots, pressure gauge, supports and frequency are unchanged.

The finite-cover audit found a correction. Q73 supplies each requested
finite angular order on the earlier continuation, but its field cap exponentiates
`T_m^2`, where `T_m` already contains `C^(20m)(m!)^2`. That order-dependent
cap cannot be compressed into Q140's fixed mixed Gevrey-3 constants. This
is the same kind of artificial loss Q141 removed on the short taper; it is
not evidence that the actual continuation is non-Gevrey.

On the post-taper interval `B<=X<=Xi=110`, Q71 proves that the complete
reference sources have exact radial spans

    p1r in span{X^-1,1,X},   X nsr in span{1,X,X^2}.       (1)

Start instead from Q141's fixed mixed reference endpoint. Make its
`F,1/F,U` traces constant in `y`, carry its five actual reference moment
rows at `B`, and reintegrate the unchanged five densities up to `Xi`.
The complete source formula then gives a uniform mixed cap, with `I,J`
in `p1r` and full `Pi,Pi_eta` in `X nsr`. Applying Q71's exact
three-point inversion to that cap bounds radial derivatives by fixed bases
one and two. Thus no all-order cap is imported from Q73.

Use the unchanged prescriptions

    a=t(1-theta)p1r+(4/5)theta,
    (log F)_y=-a/2,  U_y=-(t/2) beta X nsr,               (2)

where `t=C^-20`. The Q63 steps are Gevrey-3 with a fixed base proportional
to `t^-1`; products and primitives preserve fixed bases. For the exponential
in (2), use its actual value. Before `X=100` its exponent is at most
`t Benv^2/200`; on `[100,110]`, `a<=4/5` and `log(11/10)<1/10`.
The selected hierarchy gives `t Benv^2<1/10`, so both multiplier and
reciprocal are below `3`. No derivative norm is exponentiated.
The integration length is certified without floating point: since
`B>=4/Lambda`, the guard `Xi Lambda/4<=2^length` implies
`log(Xi/B)<length`; the selected record uses `length=Benv=10`.

At the left endpoint the incoming rows are the five ACTUAL Q141 moments.
Reintegrating their unchanged densities gives new actual `M,I,J,S,Pi`;
they are not copied from the reference or reset. The full state formula
then yields fixed `G(H,Ay,Be)` envelopes for `a,b,p2,E`, including

    X{4 Aphys eta Pi-(1-eta^2)Pi_eta}/(L E).              (3)

The same real source inequalities give

    a>=2t/(5 Lambda)>0                                   (4)

throughout: before the final interpolation `p1r>=X/10`, and the
interpolation joins positive `t p1r` to `4/5`.

The join also checks the two Q142 endpoint identities
`a_floor=2t/(5 Lambda)` and `a_cap=Benv/100`. Thus its carried moments
cannot be paired with a different `Lambda` or `Benv`; equality of the
remaining inherited profile data is an explicit identity premise.

This closes the post-taper continuation region in Q142's convention. The
shape and earlier annular restoration still need a common mixed four-field
recurrence and order-zero gate binding before taking the finite-cover
maximum and applying Q140--Q138. Tests exercise the recurrence, pressure
memory, floor, old-cap obstruction and high-order binary budget. They are
diagnostics, not a continuum proof. Uniform pulses, complete covariance
means, residual improvement and blow-up remain unproved.
