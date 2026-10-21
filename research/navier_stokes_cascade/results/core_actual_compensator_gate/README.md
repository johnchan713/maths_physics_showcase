# Actual compensator data gate and conditional shape bounds

Parent: `2b3e99d09b7973831b4fa6d804d0a03270884639`. That checkpoint proved
that the principal pointwise axial pressure coefficient depends on the full
primitive of the difference between the potential-derived reference density
and the physical pressure compensator. This checkpoint traces the actual
compensator data available in the construction and gives the precise
inequality that remains blocked.

The result is an analytical obstruction to the present data ledger, not an
obstruction to choosing a suitable compensator. No physical profile is
changed. In particular, the explicit bump used by the old numerical review
is a diagnostic fixture with its own affine scale; it is not silently
promoted to the canonical physical profile.

## 1. What the construction actually fixes

The analytical mean-repair construction says to choose a fixed smooth compact
profile `rhohat` in the reserved interior with

    integral rhohat = 1,
    rho(R,Z,T)=s^(-1/2) rhohat(R/sqrt(s)).                 (1)

The later uniform-geometry checkpoint explicitly makes no nonnegativity
assumption. It records a support hull and treats profile norms and the second
moment as finite fixed constants. No formula for `rhohat`, no values for its
derivative seminorms and no equality with the reference density are frozen in
the analytical data. Therefore the phrase "same fixed compensator" provides
consistency between steps, but not a numerical value for the new shape norm.

Let

    Dhat(x)=integral_left^x [rhohat_*(y)-rhohat(y)] dy,
    B_k=||E^k Dhat||_infinity,  E=x partial_x.             (2)

The preceding checkpoint gives `B_0 <= ||rhohat_*-rhohat||_L1/2`, but the
current data do not give the right-hand side: signed unit mass does not
control L1. Nor do mass, support and the scalar second moment determine any
`B_k`.

## 2. Exact non-identifiability at fixed mass, support and second moment

This is not merely a missing numerical evaluation. Take a normalized smooth
compact base profile `rhohat_0` and a nonzero smooth `w` compactly supported
strictly inside the same hull. For every real N set

    rhohat_N = rhohat_0 + N w'''.                          (3)

All profiles are smooth and have the same support hull. Compact integration
by parts gives

    integral w'''=0,   integral x^2 w'''=0.               (4)

Thus every `rhohat_N` has the same unit mass and second moment. Nevertheless,

    I0(rhohat_0-rhohat_N)=-N w'',                          (5)

whose supremum grows linearly with `|N|`. Hence there is no finite
profile-independent bound for `B_0` from the declared mass, support and scalar
moment data. This also proves that the scalar pressure cancellation derived
earlier cannot select the full compensator shape.

Adding nonnegativity would control `B_0` crudely because both profiles would
have L1 mass one, but it still would not supply the derivative bounds. If `b`
is a nonnegative unit-mass bump on `[0,1]`, then

    rhohat_n(x)=n b(n(x-x0))                               (6)

is nonnegative, smooth, unit mass and supported inside one fixed hull for all
large n, while its height grows like n. Since

    E Dhat=x(rhohat_*-rhohat),                             (7)

`B_1` is unbounded. Higher profile jets grow faster. Thus even positivity,
mass and fixed support would not close the slow derivative budget without
quantitative shape data.

## 3. Exact dilation formula and explicit missing constants

With `x=R/sqrt(s)`, equation (1) makes the mismatch primitive especially
simple:

    D_rho(R,Z,T)=Dhat(x).                                  (8)

Write `L=log s`. Direct differentiation gives, for either slow variable y
and any pair y,v,

    partial_y D_rho = -(1/2)L_y E Dhat,
    partial_yv D_rho = (1/4)L_y L_v E^2 Dhat
                       -(1/2)L_yv E Dhat.                 (9)

On the inherited chart `1/4<=s<=4`, the certified Cauchy estimates imply

    |L_Z|<=128, |L_T|<=256,
    |L_ZZ|<=24576, |L_ZT|<=40960, |L_TT|<=98304.          (10)

Consequently

    |D_Z| <= 64 B_1,             |D_T| <= 128 B_1,
    |D_ZZ| <= 4096 B_2+12288 B_1,
    |D_ZT| <= 8192 B_2+20480 B_1,
    |D_TT| <=16384 B_2+49152 B_1.                         (11)

These are uniform in the physical band. At every fixed higher order the same
chain rule gives a finite polynomial in the certified `log s` jets times
`B_1,...,B_m`. Smoothness says those constants exist for a fixed selected
profile; it supplies no numerical ledger value or scale comparison.

For example, the order-zero pointwise axial mismatch term from the parent
obeys the concrete sufficient estimate

    ||epsilon^5 partial_Z[D_rho partial_Z^3(KJ)]||_infinity
      <= epsilon^5 [64 B_1 ||partial_Z^3(KJ)||_infinity
                    +B_0 ||partial_Z^4(KJ)||_infinity].   (12)

This is the next blocking inequality. Neither `B_0,B_1` nor the actual target
derivatives on its right-hand side are currently bounded by the frozen
canonical data. The earlier scalar moment identity only sets one integral of
`Dhat`; it cannot replace either supremum.

## 4. Consequences and remaining scope

The radial principal coefficient from the parent remains
`-epsilon^4 K rho`. Unit mass proves its radial L1 coefficient is at least
`epsilon^4 K`, independent of the missing upper bounds. This is an operator
coefficient statement, not a lower bound for the complete residual: other
linked-cycle terms could cancel it and must be evaluated together.

The next justified action is therefore to freeze the intended physical
`rhohat` (or explicitly supply `B_k` through every required fixed order), then
insert those values and actual current-target derivatives into (12), its
higher-order analogues, the moving-density commutators and the radial term.
Only after that can the proof join generated-source supports, viscosity,
cutoffs and all interactions. No unconditional pulse bound, complete residual
improvement or blow-up theorem is established here.

## 5. Evidence

Exact compact polynomial families verify (3)--(7), including unchanged mass
and second moment and linear growth of the shape primitive. A separate
nonnegative concentrating family verifies that positivity does not bound the
first shape jet. Rational arithmetic independently evaluates every constant
in (10)--(11). The parent pointwise pressure and persistent radial coefficient
checks are rerun and historical evidence is hash-bound. These finite checks
support the algebra; the continuum conclusions are the integration-by-parts
and chain-rule arguments above.

```sh
python3 research/navier_stokes_cascade/results/core_actual_compensator_gate/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_actual_compensator_gate/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_actual_compensator_gate/evidence.json \
  --output /tmp/core_actual_compensator_gate_fresh.json
```
