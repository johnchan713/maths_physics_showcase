# Fresh postfast pressure and slow response with symbolic constants

Parent: `91b233e043f2bb8f6efd6c81854547e2dc0d5ca0`.
On 24 October 2026 its prefast CI step passed; the full run remained
in progress without failures. This checkpoint completes a CONDITIONAL
primary -> first fast -> fresh scalars -> first slow size-budget join.
It removes independent current-source, fast-state and scalar-norm slots
from this join. It does not independently prove the remaining complete
background, fixed profile/operator, chart or value/support premises.
Historical frozen evidence and original physical parameters are unchanged.

## 1. Keep giant fixed constants in an exact positive coefficient ring

The existing `ExpSum` represents finite sums c exp(e), c,e>=0 rational.
`EP` now represents positive scale polynomials sum S^k C_k with each C_k
an ExpSum. Addition, multiplication and exact historical rational
coefficient copies need no exponential evaluation or floating conversion.
Closure under these finite operations proves every expression below is
a finite positive polynomial in S with fixed real coefficients. There
is no claim of uniformity in order, tilt, width or support gap.

The complete historical `full_mean` positive formula is reproduced over
this coefficient ring. Exact regression against the original rational
formula checks every output. This extension does not change the equations,
drop pressure memory, or identify a numerical diagnostic with a field.

## 2. Build fast-state bounds from the derived primary, not row callbacks

Fix a maximum slow-output order m and ONE common covariance/source pole
p>=12(m+12)+18, padded for any additional background input pole. Retain
the original matched common band and positive radial hull. In the
parent's notation let Et_j,Ez_j be the complete PREFast rows derived
from the charted primary, full Sigma, rho, Kp(p) and Kmom(p).
The graph/radius coefficients are derived as in that checkpoint.
For output j define J_j=binom(j+5,5), Gamma_j=24*2^(j+4), and
I_j>=c_i^-1/S. The symbolic actual weighted KT/KA constants give

    T_j=2 J_j Gamma_j I_j S (j+4)!
                              [Et_(j+4)+(1+KA_j)Ez_(j+4)],
    X_j=2 J_j Gamma_(j+1) I_j S KT_(j+1)(j+5)! Ez_(j+5).

T is a common sufficient bound for BOTH tangential and axial velocities;
it includes the actual axial gd-A1 gd. The factor2 is the zero-Haar
projection bound. The actual field orders are epsilon^mu for tangential
and axial velocities and epsilon^(mu+1) for the radial velocity,
mu=1-kappa. All carry zeta delta^-p. Neither the source factorial nor
the five-variable output conversion is omitted. No first fast field
norm is supplied independently.

The fixed callbacks supply only original background jets and fixed
operator/profile constants. The prefast source is rebuilt directly in
the extended coefficient ring so its Kp/Kmom can ALSO be symbolic.
This closes the prior rational-only diagnostic interface limitation.
The code does not materialize giant weight floors.

## 3. Recompute WHOLE pressure and scalars from the updated state

For fresh output n take state inputs through N=n+3:
use x=X_N, u=z=T_N as sufficient coefficients, the unchanged COMPLETE
primary covariance w_N, and the same full base b,V,G,Sigma. In
`core_signed_prefix_source_budget.full_mean` the exact positive recipe is

    srr=2b x+x^2+w, srz=b z+G x+x z+w, stt=2V u+u^2+w,
    Gr=time*x+(dr+r)srr+d*srz+r*stt+vis*x,
    P=Kmom Gr, Pm=Kp(Gr+Rho P),
    Jtheta=Kmom(G u+V z+z u+w),
    Jz=Kmom(2G z+z^2+w+Gr).

The original d=n+4, graph/time and viscosity coefficients and both
complete tangential row formulas are retained in code. These quantities
are recomputed AFTER the fast update; the prefast scalar norms are not
reused. All mean squares, cross products and all six covariance entries
remain. The pressure compensation uses the WHOLE current P and moving
rho, not a leading or projected surrogate. Pi1/Pi2 and additive flat
residuals remain subject to the same full-base subtraction as before.

Products of two fast fields have weight zeta^2 delta^-2p, bounded by
zeta delta^-2p. Pad single weighted fields to pole2p and use Kp and
Kmom at 2p for this fresh recomputation. The full conservative size
proof gives order mu for these rows/pressure/scalars: radial covariance
transport has order1-kappa, and tangential graph viscosity has order
mu+1-2kappa>=mu. No sharp Haar cancellation is used for unaveraged
rows, no mean product is assigned pole p by assertion, and no residual
improvement is inferred from this size estimate.

## 4. Apply the SAME slow response to FRESH scalar inputs

Let H,Kp_s,Kz_s,Kg_s,Kpsi_s be the already derived fixed slow reciprocal
and response-profile norms from `core_slow_inverse_budget`, with its
distinct physical exponents and actual fixed 1/lambda conditioning.
For total field jets, include all mixed radial profile derivative
bounds from its composition formula, not just slow-coordinate jets.
Using the FRESH scalars above, the existing response gives

    du_n=H_n(Kp_s,n P_n+Kz_s,n Jz_n),
    dz_n=H_n Kg_s,n Jtheta_n,
    dx_n=(n+1)H_(n+1)Kpsi_s,n+1 Jtheta_(n+1).

The radial FIELD retains the extra epsilon. The code uses the NEXT
fresh scalar row for dx; it does not reuse Jtheta_n. Reuse the ALREADY
proved `core_reserved_geometry_recheck.geometry` on the ENTIRE moving
potential hull [l,u] inside I4, including gaps filled by its primitive.
It supplies E=w^2/da^2+4/db^2 and zeta>=exp(-E) on that actual support.
Multiply each slow coefficient by exp(E) SYMBOLICALLY. All field jets
remain on the same moving hull under smooth zero extension, so this
gives full zeta weight with pole0, for differentiated fields rather
than derivatives of field/zeta. No weight floor is asserted on the
Cartesian union integration hull. The constructor retains the SAME
flat width as the primary; actual reserved support inclusion remains
an original geometry premise, not something manufactured by the API.

Finally ADD the slow response to the fast state; do not reset the old
means. The final size is tangential/axial epsilon^mu zeta delta^-p and
radial epsilon^(mu+1) zeta delta^-p. This is a conditional first-prefix
size bound, not a sharper full row or complete residual estimate.

For slow output m the safe covariance order is m+12: fresh scalar
order m needs fast fields through m+3; their source route needs primary
covariance through m+11; the radial slow response adds one order.
The original-root/background jet costs INSIDE that covariance recipe
remain. This is a finite requested-order ledger, not a minimality claim.

## 5. Validation and remaining boundary

Exact coefficient-ring and historical-row comparisons, fresh/stale
scalar controls, whole-pressure and A1 sensitivity, moving-hull rejection
and source-order traces test the implementation. Full pipeline traces
use a marked manufactured covariance subclass to keep diagnostic sizes
moderate; the actual lower primary checkpoint is separately replayed.
No selected giant continuum field is numerically evaluated.

The coarse slow bound above DOES NOT discharge the sharper initial Haar
contract, which uses slow tangential/axial order1 and radial order2.
Here mu=1-kappa<1. A bound epsilon^mu C cannot imply epsilon Cprime
uniformly in the band: the ratio is epsilon^-kappa, unbounded as
epsilon tends to zero. Before using that sharper signed-target route,
derive the better fresh scalar cancellations (or retain their existing
stronger hypotheses) so that the actual slow response satisfies

    |du|+|dgamma| <= epsilon zeta delta^-p U(S),
    |dbeta| <= epsilon^2 zeta delta^-p X(S).

This is an explicit blocking inequality, not a failure of the physical
construction. The old sharper initial scalar identities are not contradicted
by a weaker sufficient size estimate. They must be joined with their full
moving pressure and actual updated covariance, rather than inferred from
this coarse join. The finite rate control detects that invalid promotion.

Remaining inputs are the actual full background/phase/value gates,
original stress and fixed map/profile norms on the same hull, actual
cutoff plateaus, full pressure/moment identities and support premises.
Next carry this accumulated mean state into the original signed-target
and first sourced-wave bounds, retaining the whole moving c_rho P and
fresh covariance updates. Do not silently extend this initial real-primary
join to a later complex wave or promote coarse size to improved residual.
All physical-scale realization, viscosity, cutoff tails, nonlinear
interactions and complete residual closure remain open. Unconditional
pulse bounds, infinite realization and blow-up remain unproved.
