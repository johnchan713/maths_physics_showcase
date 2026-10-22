# Generated sourced covariance and complete pressure-moment bounds

Parent: `04aa29f2481efce6e6e9ce4bb1f3108341a83f5c`.
The preceding normal/frame estimates now feed the actual first-source
budget, its inverse, the SAME cut potential and pressure, the complete
cutoff curl, and all six covariance entries. This removes an independent
norm input for the NEW sourced covariance. Explicit scalar jet bounds
then follow for all three changed defects and the whole moving axial
pressure moment. These are conditional estimates from the retained
original primitive, old-wave and source/operator inputs, not actual-fluid
enclosures or a complete residual improvement.

The exact covariance/mean identities already exist in the frozen research.
The new contribution is their quantitative composition with the generated
inverse and the following correction: a next mean update cannot simply
reuse operator bounds certified only at the previous edge poles. In the
unweighted-fixed-coefficient route below its required poles are 6p and
12p, whereas the previous pressure calls were padded only for 2p and 4p.
This is a proof-domain mismatch, not failure or nonexistence of the maps.
Historical evidence and all actual profiles/supports remain unchanged.

## 1. A one-way generated-wave budget, including the cutoff curl

Work on one frozen label/common lifted chart with five coefficient
variables and the same phase, envelope P, normal and background. The
BEFORE prefix is exactly `core_signed_prefix_source_budget`: its actual
fresh fast/slow source obeys the conditional estimate at B=mu=1-kappa.
Its source callback requests only earlier field/wave/operator bounds.
The inverse uses `core_tangent_propagator_closure`, hence the original
4/bmin^2 denominator, derived T and length, with full nprime and n.f in
pressure. The real first source has signed harmonics {-2,-1,1,2}.
The common harmonic-l1 source bound bounds each coefficient; equal upper
bounds for conjugate modes do not identify their complex coefficients.

For raw coefficient jets X_j define the through-order factorial norm

    J_n(X) <= sum_(j=0..n) 5^j X_j/j!.

Conversely a through-order cutoff bound gives raw |D^alpha psi| <=
j! J_j(psi), |alpha|=j. `GeneratedWave` uses precisely this conversion
from the prefix's original temporal cutoff. It applies that SAME cutoff
to potential AND pressure using every Leibniz term.

Let t_m be the tangent inverse velocity and
A_m=i n cross t_m/(km|n|^2), with its epsilon^(B+1/2) normalization.
The full cut velocity coefficient is

    a_m=psi*t_m + curl_*(psi*A_m).

The phase part of curl equals psi*t_m by exact n.t_m=0. In cylindrical
components the remainder is

    (-Dz(psi Atheta), Dz(psi Ar)-Dr(psi Az),
      Dr(psi Atheta)+(psi Atheta)/R).

The amplitude has no separate physical-angular dependence. All graph and
cylindrical connection terms remain. Let g_n bound the full graph
coefficient norm in Dr, so Dr costs epsilon^-kappa g_n times a common
coefficient derivative; Dz=epsilon partial_Z. With r_n bounding R^-1,
and C_n=J_n(psi A_m), a safe remainder polynomial is

    R_m,n = 2(n+1)(g_n+1) C_(n+1) + r_n C_n,
    U_m,n = J_n(psi t_m) + R_m,n.                         (1)

The factor n+1 is the normalized derivative cost. This uses one MORE
potential jet; it does not drop curl-of-cutoff velocities. The principal
order is B and the remainder order B+1/2-kappa>B, so U_m,n bounds the
complete coefficient at order B. The cut wave pressure has order B+1/2
and retains the direct unprojected source term.

Sum (1) over all four signed modes in harmonic l1. No inverse output
frequency or unproved real-primary quadrature cancellation is used.
Original closed label separation removes distinct-label derivative
products; common-coordinate/aggregate constants must already be in the
retained BEFORE inputs. No new label-count assumption is introduced.

## 2. All six generated covariance entries

Let Uold_n bound the complete OLD wave at order 1/2, and Uadd_n the sum
above at order B. For D=<wnew tensor wnew-wold tensor wold>_theta,
the normalized coefficient norm of EVERY entry is bounded by

    D_n(S) = 2 Uold_n Uadd_n + Uadd_n^2.                  (2)

Products use the submultiplicative through-order norm before angular or
auxiliary averaging. Equal angular frequencies are combined FIRST; their
phases cancel in the covariance before coefficient differentiation. In
real quadratures the exact expression is

    (1/2) sum_(m>0) [a_m c_m^T+c_m a_m^T+b_m d_m^T+d_m b_m^T
                          +c_m c_m^T+d_m d_m^T].

In particular the second-harmonic square survives despite the absence of
an old second harmonic. Both quadratures, axial entries, principal/curl
crosses and curl squares are retained. A covariance DIFFERENCE need not
be positive semidefinite. This is a bound, not a two-flux substitute for
the complete tensor.

Because B>=1/2, C=B+1/2=3/2-kappa is the covariance order: 2B>=C.
If the incoming common pole is p and fixed coefficient/cutoff jets are
unweighted, the prior source has pole 3p. The inverse and (1) retain 3p;
the cross in (2) has pole 4p, the self-square 6p. Thus the full bound is

    |D^alpha D_ab|/(alpha! normalization)
       <= epsilon^C zeta delta^-6p P^2 D_n(S).           (3)

One finite p is chosen to cover all input jets requested for the chosen
output order. This is not a universal pole for every derivative order.
If a fixed coefficient has an additional pole, add it at its actual
product before applying this ledger. It is not silently included in 3p.

## 3. Full-axis pressure and moment response

Keep means fixed during the wave update. After auxiliary Haar averaging,
write the same D_ab for the averaged change. The exact full responses are

    dg=-D1^0 Drr-epsilon partial_Z Drz+Dtt/R,
    dP=integral (Dtt-Drr)/R-epsilon partial_Z integral Drz,
    dJtheta=integral R^2 Dtheta,z,
    dJz=integral R[Dzz-(Drr+Dtt)/2]
                    +(epsilon/2) partial_Z integral R^2 Drz. (4)

All integrals are on the complete positive radial axis, using compact
extension. The changed pressure is the SAME T0(dg-rho dP). Its actual
A0 remainder stays in the pointwise equation; its Haar mean is zero.
These are INCREMENTS: add them to the existing current defects and
pressure. They are not a reset of old means or an assumed bound for the
new total defects.
For c_rho=(integral R^2 rho)/2, exactly

    integral R bar(dp) = -(integral R^2 dg)/2+c_rho dP,
    axial residual moment change
       =epsilon partial_Z(dJz+c_rho dP).                 (5)

The compensator may extend beyond an individual wave's support; (5)
uses its FULL second moment. Both background Pi1/Pi2 tails cancel only
on subtraction of the SAME full base, not by discarding axial gradients.
The nonzero wave pressure has zero angular mean, not zero pointwise value.

For unweighted scalar bounds use the original fixed-width edge constant
E_(6p) from `core_primary_aggregate_atlas` and set V_n=E_(6p) D_n.
It bounds zeta delta^-6p because zeta<=sqrt(zeta), with P<=1. This does
not divide a scalar by a radial weight or prove a new weighted mean norm.
Let [a,d] contain ALL changed-wave supports on the full common chart,
0<a<d, and M_j=(d-a)max(a^j,d^j), j=-1,0,1,2. Use the SAME original
uniform frame/field hull; the generated API takes it from Rmin,Rmax.
The reserved I4 mean hull cannot replace that full covariance hull.

At order n, with c_n the unweighted through-order norm of the whole
moving c_rho, sufficient polynomials after factoring epsilon^C are

    p_n = 2 M_-1 V_n + (n+1) M_0 V_(n+1),
    q_n = 2 M_1 V_n + (n+1) M_2 V_(n+1),
    jt_n = M_2 V_n,
    jz_n = 2 M_1 V_n + (n+1) M_2 V_(n+1)/2,
    pressure_moment_n = q_n/2 + c_n p_n.                 (6)

The input axial epsilon in (4) is conservatively bounded by one here.
The residual moments have the EXTRA epsilon, hence order C+1 and bounds
(n+1)jt_(n+1) and (n+1)[jz_(n+1)+c_(n+1)p_(n+1)]. The entire moving
product is differentiated, including c_rho,Z dP. These statements request
covariance jets through n+2, potential/source jets through n+3, and,
with the prior safe prefix recipe, initial mean jets through n+18. No
output is reused to estimate its own source. All fixed-order constants
may depend on the original width, profiles and pole; no first band changes.

## 4. Correction before the NEXT fresh mean maps

The full wave change is now bounded, but post-wave means are not silently
certified using operators for a different input space. The next common
input pole must cover old means at 2p and the new complete wave at 3p;
3p suffices. Substituting that into the EXISTING finite-prefix ledger gives

| Next operation | Required pole |
| --- | ---: |
| Full covariance and prefast rows | 6p |
| KT1, KA1 fast-map input | 6p |
| Postfast full means | 6p |
| Recomputed postfast rows and moments | 12p |
| Kp, Kmom calls to pad | 6p AND 12p |
| Final means / following regular source | 6p / 9p |

Earlier Kp/Kmom constants were padded for 2p and 4p, not these domains.
For p>0 no finite old 4p input norm controls a generic 6p input:
on an inner collar, f=zeta*y_a^-6p is smoothly flat at the edge and

    f/(zeta*y_a^-4p)=y_a^-2p -> infinity.

Subtracting its finite integral times the same interior unit bump imposes
zero total without changing this divergent inner-edge ratio. The pressure
compensation constraint therefore does not fix the input-norm mismatch.

Thus those old proof domains cannot certify the next calls. This does
NOT prove the actual generated covariance saturates the worst pole, or
that the maps fail. It proves why their bounds must be reevaluated at the
required poles, or a sharper specific estimate supplied. The retained
shifted-primitive recipe exists at EVERY finite pole with the SAME width,
gap and supports; its kernel degree changes with the pole. This is an
operator-constant recalculation, not a new selector or smallness condition.
`operator_coverage` rejects the stale domains and records the exact missing
ones. This checkpoint does not claim that all actual new operator constants
or the ensuing fresh fast/slow field bounds have been supplied.

## 5. Ledger, validation and next obligation

Removed independent norm input: the NEW sourced complete covariance and
its unweighted pressure/moment response, now generated by (1)--(6).
Added analytical assumptions: none. Retained: valid original primitive
and BEFORE wave/source/fixed-map bounds on one full domain, original
root identity where required, full-path supports, and the same fixed
cutoff/width. Exact support propagation still uses those retained lemmas;
pointwise or endpoint source vanishing is not a replacement.

The exact controls check all six entries with both harmonics/quadratures,
a separate full-box graph/cylindrical curl bound, complete compact pressure
moments and a moving compensator. Generated positive-polynomial tests
trace the n+18 dependency and preserve full pressure. Pole-domain controls
detect reuse of old constants and verify the unchanged primitive recipe
at the required finite poles. Manufactured finite fixtures are not actual
fluid enclosures; the analytical estimates are the derivations above.

Next: reevaluate the SAME weighted mean operators at 6p and 12p, then
apply the fresh fast and slow maps and recompute full pressure/scalars
after EACH update. Keep old means, A0/A1 velocities, all interactions and
the separate additive tail ledger. The full signed target/pointwise
residual is not bounded by (6) alone. No unconditional pulse bound, full
residual improvement or blow-up is claimed.

```sh
python3 research/navier_stokes_cascade/results/core_sourced_covariance_moment_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_sourced_covariance_moment_join/audit.py \
 --verify-record research/navier_stokes_cascade/results/core_sourced_covariance_moment_join/evidence.json \
 --output /tmp/core_sourced_covariance_moment_join_fresh.json
```
