# Uniform reserved-hull geometry and explicit pressure-compensator jets

Parent: `4c082f5a3d343643dced20f22497badc33d695d8`. The dynamic-mean norm
budget remains conditional on its actual field norms. This checkpoint
discharges its scale-coordinate, reserved radial geometry and compensator
costs, with explicit constants independent of the band index. It does not
supply all generated-field or sourced-inverse input norms. Historical frozen
evidence, the common primary tail, fixed tilt and width C_setup^-20 are unchanged.

The implicit chart and identity c_rho=s*c_hat already appear in
`physical_pulse_foundation` and `compact_mean_corrections`. The new analytical
result is an explicit all-orders jet majorant, a certified reserved-support
weight comparison and their substitution into the dynamic budget. We do not
claim those identities as newly discovered or infer continuum bounds from tests.

## 1. A uniform complex neighborhood for the coordinate alone

On the same accepted physical chart,

    s-Z^2 s^(2h)=T, eta=Z s^(-D), D=1/2-h,
    1/4<=s<=4, |eta|<=1, 0<h<1/100.

Fix a physical point (s0,Z0,T0). Put s=s0(1+t),
p=(Z-Z0)s0^(-D), q=(T-T0)/s0. Its EXACT normalized equation is

    t=q+(eta0+p)^2(1+t)^(2h)-eta0^2.                         (1)

Use the branch of the power analytic on |t|<1, with value 1 at zero.
For complex |p|,|q|<=1/16 and |t|<=1/4, the derivative of the right side
in t has magnitude at most

    (6/5)(1/50)(4/3)=4/125 < 1.

Here |eta0+p|^2<6/5 and |1+t|^(2h-1)<=4/3. At t=0 the right side is
bounded by 3/16+1/256=49/256. Therefore the map sends the closed t-disk
strictly into itself: 49/256+1/125=6381/32000<1/4. Uniform contraction
gives one holomorphic solution, agreeing at the center with the selected
real branch. Its |t| bound is 1/4. The constants work also at eta0=+-1,
T0=0. This complex extension concerns ONLY the algebraic scale coordinate;
it does not analytically extend a heat profile or evaluate it at negative
heat time. Physical endpoint derivatives remain their one-sided limits.

Cauchy's formula on this bidisk gives normalized Taylor coefficients
|[p^i q^j]t|<=(1/4)16^(i+j). Since s0^(-D)<=2 and s0^-1<=4,

    |partial_Z^i partial_T^j s| <= i! j! 32^i 64^j,
                          i+j>=1.                         (2)

For order zero use |s|<=4. No Q or band index occurs. For any fixed real
power a and integer P>=|a|, the same argument uses
|s^a|<=(16/3)^P on the complex disk and proves

    |partial_Z^i partial_T^j s^a|/(i!j!)
       <= (16/3)^P 32^i 64^j,  i,j>=0.                    (3)

This covers all three retained q-to-Q target factors s^(2A),
s^(2A-3/2), s^(2A-1), A=1/2+h, with P=2, as well as the fixed physical
normalization powers. It does not remove the independent fixed-lambda
inverse cost or derivatives of the actual a(eta) profile.

## 2. The entire moving pressure product has explicit uniform constants

For the SAME physical unit-integral compensator, with fixed smooth compact
rhohat (no nonnegativity assumption needed),

    rho(R,Z,T)=s^-1/2 rhohat(R/sqrt(s)),
    c_rho=(1/2)integral R^2 rho = s*c_hat,
    c_hat=(1/2)integral x^2 rhohat(x) dx.                    (4)

The identity follows by change of variable; it is not an approximation of
pressure by its leading term. Let c*=|c_hat|. If its support is [l,u], the
always-valid bound is c*<=u^2 ||rhohat||_L1/2. Unit signed mass alone does
not imply ||rhohat||_L1=1. For the normalized slow coefficient norm used in
`core_dynamic_mean_recheck`, (2) yields

    |c_rho|_m <= C_m=c* [4+sum_(1<=i+j<=m) 32^i 64^j].     (5)

In particular this supplies every c_(m+1) in that checkpoint's complete
axial residual-moment estimate. With the conditional P/Jz budgets there,

    |epsilon partial_Z(Jz+c_rho P)|_m
      <=epsilon(m+1)[jz_(m+1)+C_(m+1) p_(m+1)].             (6)

This retains c_rho,Z P. At a real point L=1-2h eta^2,
s_T=1/L and s_Z=2 Z s^(2h)/L, so c_rho,Z is generally nonzero.
The budget (6) uses actual field input norms through m+3 as before; the
compensator part of those constants is now explicit for every order. Full
Pi1/Pi2 background pressure and all cutoff derivatives remain in the
incoming mean/base equations; no equality with leading pressure is inferred.

## 3. Moving support is not a rectangular weight-floor assertion

Let the fixed q-normalized slow potential hull be x in [l,u], 0<l<u,
strictly inside the actual reserved I4. Here x=r/sqrt(q) and X=x^2/2.
Its Q-chart support at each slow point is

    sqrt(s) l <= R <= sqrt(s) u.

Since 1/4<=s<=4, one common integration hull is [l/2,2u]. Thus the explicit
radial moment constants M_j in the preceding dynamic budget can be chosen
using a=l/2, d=2u, independently of the band. Derivatives of smooth zero
extensions retain the actual moving support. Potential primitives can fill
the gap between bumps, which is why [l,u] must be the whole convex hull.

Suppose Xa<l^2/2<=u^2/2<Xb for the SAME active annulus. Set

    r_a=l^2/(2Xa)>1, r_b=2Xb/u^2>1,
    d_a=2(r_a-1)/(r_a+1), d_b=2(r_b-1)/(r_b+1),
    E=t1^2/d_a^2+4/d_b^2,  t1=C_setup^-20.

For r>=1, log r>=2(r-1)/(r+1): subtract the right side and differentiate
to obtain (r-1)^2/[r(r+1)^2]>=0. Hence on the ACTUAL moving hull,

    zeta=exp(-t1^2/log(X/Xa)^2-4/log(Xb/X)^2) >= exp(-E).

Consequently, for any smooth field supported in that moving hull,

    sum_(i+j<=m) sup |partial_Z^i partial_T^j f|/(i!j! zeta)
       <= exp(E) |f|_(m,1).                               (7)

The weighted supremum is taken where zeta>0, with zero extension outside;
all field jets vanish off the actual moving hull. The weight divides the differentiated field, as in the retained coefficient
classes; (7) is not a claim about derivatives of f/zeta. E is recorded as a
logarithmic constant, without materializing exp(E). This avoids destroying
the actual width constant, which may be enormous in reciprocal estimates.
These are explicit, uniform-in-band geometry and comparison costs. Their
dependence on the actual fixed hull, width and background is retained.

There need NOT be a positive zeta floor on the Cartesian product
[l/2,2u] x Omega: an R admissible for one s can lie outside the active
annulus for another. The fields vanish at those off-support points. Using
the union hull for integration and the actual moving hull for (7) is
essential. If a covariance change has wider radial support, its hull and
edge behavior must be handled separately; the I4 floor is unavailable for
wave factors approaching the active shell edges. The oriented edge bounds
and the full source grammar remain required there.

## 4. Join to the current dynamic budget, with the remaining inequality

For a fresh slow update all six D_ab are zero. Substitute a=l/2,d=2u and
the explicit (5) into the preceding budget, taking omega=1. The other input
norms must still be valid unweighted current-field bounds through m+2;
the joint residual-moment budget uses them through m+3. This now supplies
the shell and c_rho costs for EVERY accepted band on the existing common
tail, with no further band threshold and no change of profile or tilt.

The remaining quantitative inequality is that the resulting
max(p_m,jtheta_m,jz_m) and (6) fit the proposed epsilon-order tolerances
with constants independent of the band. This checkpoint does not supply
those actual generated-field constants. The inherited slow inverse gives
their conditional orders and fixed-lambda factors; it must be joined to
the actual target norms. Full changed-wave supports, fast correlations,
all sourced-inverse admissibility bounds and the complete pointwise
viscous/cutoff residual remain open. No unconditional pulse bound,
complete residual gain or blow-up is claimed.

## 5. Evidence

Exact bivariate rational Taylor jets solve (1) degree by degree using the
nonzero divisor 1-2h eta0^2. They check the defining equation independently,
the explicit coefficient bounds, both polar endpoints and a tiny fixed h
without constructing a huge physical scale. Dropping the divisor fails a
first-derivative control. Independent radial polynomial integrals check (4)
on moving intervals. Exact support and logarithmic gap certificates test
(7); the rectangular-hull negative control is outside the active annulus.
These finite diagnostics do not prove all orders: sections 1--4 do that
for the stated coordinate/geometry conclusions, not for the remaining fields.

Run (standard-library diagnostics):

```sh
python3 research/navier_stokes_cascade/results/core_reserved_geometry_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_reserved_geometry_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_reserved_geometry_recheck/evidence.json \
  --output /tmp/core_reserved_geometry_recheck_fresh.json
```
