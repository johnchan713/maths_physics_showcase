# Q41: same-state postslow scalars and whole moving pressure moment

Parent: `02f3a428a5536e0515b6722806f271dbad1f45c4` (Q40).
Parent focused CI 41146882107 passed; broad CI 41146881998 was still
running at inspection on 31 October 2026. This checkpoint composes the
retained Haar scalar reset with Q39/Q40, without an independent outgoing
scalar callback. It is a conditional finite norm recipe, not evaluation
of the actual high-order field. Historical frozen evidence is unchanged.

## 1. Full old means for fast changes, Haar old means for slow changes

Q39 correctly retains FULL old means in the new fast scalar change.
Zero Haar of the new fast field does not kill its products with the old
fast field. That conclusion remains unchanged.

For an auxiliary-independent slow increment S, however,

    Haar(O S)=Haar(O) S.                                    (1)

The initial fast and new complete fast fields have zero auxiliary Haar;
the real signed wave changes no mean velocity. Thus the old perturbation
Haar mean before the new slow solve is exactly the retained INITIAL slow
mean, with tangential order 1 and radial order 2. It is not zero, and it
is not the full old fast+slow envelope with the weaker order 1-kappa.
The code retrieves it from the SAME initial slow-field constructor.
The new fast still includes its actual axial gd-A1 gd term; (1) uses its
exact zero-Haar property, not omission of that field from other terms.

The new slow increment is Q39's response to CURRENT postfast scalars,
which already contain the complete signed-wave change once. No old scalar
or pressure is reset before this solve. Exact prescribed swirl moments,
G=0 on the entire reserved potential hull, the five-map response identities
and flat zero extensions are premises of the linear cancellation below.

## 2. Retain all nonlinear scalar remainders after the linear solve

At scalar output n use old Haar and new slow norms through n+2, both
unweighted on the whole reserved potential hull. Write old components
as (ox,ou,oz), increments as (x,u,z), and the original radial base bound
as b (with its retained extra epsilon). Put

    S=2(b+ox)x+x^2,
    H=(b+ox)z+oz*x+x*z,
    T=2ou*u+u^2,
    K=oz*u+ou*z+u*z,
    Z=2oz*z+z^2,
    tau=(n+1)+(n+1)(n+2).

For a common unweighted radial moment factor M, the inherited positive
remainder functional gives

    P_after <= M[tau*x+S+T+(n+1)H],
    Jtheta_after <= M K,
    Jz_after <= M Z+(M/2)[tau*x+S+T+(n+1)H].                 (2)

These are bounds on the remainder AFTER the exact linear response, not
a rule setting pressure or scalars to zero. Old/current scalar defects
enter the new slow amplitudes, and their nonlinear effects remain in (2).
The original complete pressure and its archived components are not
replaced by a leading pressure or by only its latest change.

For the whole reserved potential hull [a,bh], the implementation uses

    M=(bh-a)(a^-1+1+bh+bh^2).                              (3)

This dominates weights -1,0,1,2. The hull is taken from the initial slow
budget, separately from the moving compensator/projector hull. Since
both slow fields are unweighted inputs, no second edge unweighting or
small-bump support restriction belongs in (3).

With alpha=2-3kappa, new slow tangential/radial orders alpha,alpha+1,
and old Haar orders 1,2, the possible worst orders in (2) are
alpha+1, 2alpha and alpha+2. For kappa=1/100000 the minimum is alpha+1.
Using full old means here would only give the weaker coarse alpha+1-kappa;
using Haar old means in the FAST update would be invalid. The distinction
is a consequence of (1), not a blanket license to discard correlations.

If G is not zero on the WHOLE potential hull, the separate linear scalar
leaks involving G*u and 2G*z are not removed by (2). Neither the wrapper
nor a finite test verifies the continuum reserved-hull premise.

## 3. Whole axial pressure moment and finite composition ledger

The resulting axial moment remains

    epsilon partial_Z(Jz_after+c_rho P_after),               (4)

including c_rho,Z P_after. For normalized derivative order n its sufficient
coefficient is

    (n+1)[Jz_after,n+1+C_rho,n+1 P_after,n+1].                (5)

C_rho is Q36's SAME full moving-profile norm. Formula (5) has conditional
order alpha+2; it does not freeze the compensator or use only a pressure
change. The angular scalar in (2) stays a distinct channel.

For scalar and whole axial outputs through m, the factory reserves:

| Layer | Maximum required order |
|---|---:|
| Postslow scalar for the axial derivative | m+1 |
| New slow fields | m+3 |
| New fast fields | m+5 |
| Current post-signed rows | m+10 |
| Q37 pressure ledger | m+11 |
| Q36 signed target | m+14 |
| Q36 fixed profile / primary covariance | m+27 |

The last entries compose the existing finite index bounds. They are NOT
expanded canonical primitive counts or substitutes for the larger A0/A1
compact-tail graph. The factory binds Q40's row and fast objects to Q39's
current scalar state, checks object identity, and rejects an anonymous
fast callback. It refuses insufficient derivative reserves before running
the norm program. Tail, compact/time and background inputs remain the
same explicit lower contracts as Q40.

## 4. Verification and remaining scope

Independent rational controls verify every term in (2). Angular averaging
controls distinguish the legitimate slow use of Haar from the illegitimate
fast shortcut. Symbolic differentiation retains both moving-pressure
product terms. Factory tests check the full ledger, same-state identity,
whole potential hull and rejection of a free fast callback without
materializing its enormous high-order coefficient programs. The existing
providers' scalar, pressure and refresh checks are rerun. These tests are
not a continuum proof or numerical evaluation of the selected field.

Next bind the original initial-tail and compact/time-map contracts across
the full required graph, and assemble the complete final row differences
and unaveraged source with these fresh means. The zero strict margin in
the first-real fast difference is still present; the scalar gain here does
not prove a stronger complete residual rate. Actual weighted high jets,
support/value premises, uniform pulse bounds, full viscosity/cutoff
interactions, residual improvement and blow-up remain unresolved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.
