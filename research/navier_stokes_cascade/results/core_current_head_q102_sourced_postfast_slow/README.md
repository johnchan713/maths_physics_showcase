# Q102: complete sourced postfast scalars and fresh slow response

Parent: `5be85f46c17ba97514553f01f795b538eef763a5` (Q101).
Its CI runs 41450637373 and 41450637089 were still running when inspected
on 2 November 2026. Q100 field-comparison run 41441649281 had passed.
This checkpoint composes Q101's complete sourced fast norm map with the
frozen full scalar-correlation identities and original slow five-map.
It derives conditional postfast scalar/slow coefficients on the CURRENT
state. These identities were already established in core_actual_scalar_gain_join
and core_preserved_haar_scalar_join; their current Q99/Q101 binding and
derivative/pole requests are the new contribution.

## 1. Add the sourced wave change exactly once

Let O be the complete incoming first-real FINAL mean velocity, F the new
Q101 fast increment (including its axial gd-A1 gd term), and W the full
changed Q97 covariance. During this fast mean update W is FIXED. Thus

    Theta_before_fast=Theta_old_postslow+delta Theta_wave,
    Theta_postfast=Theta_before_fast+delta Theta_fast.      (1)

The first line is Q99's SAME ChangedScalars object. There is no extra wave
term in the second line and no reset of the old pressure or scalar state.
The implementation adds positive norm coefficients, not physical signed
scalar values. Old means are derived from the SAME Q86 FinalMeans.fields;
they include all previous fast and slow components, not only their Haar.
No desired postfast-scalar or slow-field callback is introduced.

All three fast components have exact zero auxiliary Haar under the retained
compact/time inverse identities. This cancels their linear interactions
with auxiliary-independent b,V,G in the scalar integrals, and the Haar of
their linear time and viscous terms. It does not remove products with old
nonzero-Haar fields. In particular the old axial velocity and new axial
A1 correction must remain in the correlation tensor. Time and viscosity
remain in the full pointwise momentum rows; only their SCALAR CHANGES
have these exact cancellations.

## 2. Complete scalar polynomial and pressure moment

For normalized through-(n+1) full mean bounds O=(ox,ou,oz), F=(x,u,z), put

    S=2 ox x+x^2, T=2 ou u+u^2,
    H=ox z+oz x+x z,
    K=oz u+ou z+u z, Zf=2 oz z+z^2,
    M=(b-a)(a^-1+1+b+b^2).                              (2)

The fixed positive interval [a,b] is the WHOLE current radial/potential
hull, not the reserved bumps alone. The inherited sufficient scalar norm
coefficients are

    delta P_fast <= M[S+T+(n+1)H],
    delta Jtheta_fast <= M K,
    delta Jz_fast <= M Zf+(M/2)[S+T+(n+1)H].             (3)

All products are formed BEFORE Haar. To audit the signs in the exact
integral identities, let S,H,T now denote the actual Haar products and
write Dz=epsilon partial_Z. Full radial-source subtraction gives

    Haar(delta gr)=-partial_R S-S/R-Dz H+T/R,
    delta P_fast=int (T-S)/R dR-Dz int H dR,
    delta Jz_fast=int R Zf dR-int R(S+T) dR/2
                                      +Dz int R^2 H dR/2. (4)

These require full compact endpoint/mixed-trace matching for integration
by parts. Unmatched seams and exterior traces need the corrections already
specified in Q98; numerical zero samples or separate patch supports do not
prove those contracts. Absolute values in (3) replace signs in (4), and
drop only beneficial epsilon factors. Both direct and pressure-induced
axial terms are retained.

The current moving pressure moment still enters as the WHOLE product

    epsilon partial_Z(Jz_postfast+c_rho P_postfast),
    c_rho=(int R^2 rho dR)/2.                             (5)

The code bounds (5) through n by (n+1) times the through-(n+1) norm of
the entire sum, using the original c_rho(n+1) recipe (with safe 243 norm
padding). It does not freeze c_rho or delete c_rho,Z P. Complete old/base
pressure including Pi1/Pi2, Q97 covariance memory and all earlier scalar
correlations remain additive. This is not a new pointwise postfast T0
pressure estimate; that will be needed in subsequent final-row budgets.

## 3. Norm domains and conditional grades

At scalar output n, unweight the Q101 fast norm at its OWN exact returned
pole p, using Eedge(p,width). Apply a safe factor 243 for raw-max versus
five-variable factorial-through conversion. Q86 already pads its complete
old unweighted means. Only then form (2)--(3). The scalar products are
unweighted integrals; no monotonicity of weighted operator constants or
silent same-pole squaring is assumed. Both fields must obey the input
jet bounds on the same complete hull and closed physical cover.

For B=1-kappa, D=3/2-2kappa and C=3/2-kappa, old mean tangential orders
are B and new fast orders D (radial orders B+1,D+1). Therefore

    delta Theta_fast has sufficient order min(B+D,2D)=B+D,
    B+D=5/2-3kappa > C,
    Theta_postfast has order C,
    (5) has order C+1=5/2-kappa.                          (6)

Q99's old-postslow plus wave scalar change is used ONCE at grade C; it
dominates the smaller fast correction in (6). These are normalized
conditional grades for 0<kappa<1/4 and 0<epsilon<=1, not numerical continuum
enclosures, physical absolute-decay results or a complete residual gain.

## 4. Same original fresh slow five-map

Let a bound the new postfast scalars through n, and an through n+1. Read
the original budget.slow(n) and slow(n+1) coefficient programs, with exactly
the keys H,Kp,Kz,Kg,Kpsi. The complete original response is

    u_s=H_n(Kp_n a_P+Kz_n a_Jz),
    z_s=H_n Kg_n a_Jtheta,
    x_s=(n+1)H_(n+1) Kpsi_(n+1) an_Jtheta.              (7)

This gives conditional tangential/axial order C and radial order C+1.
The optional weighted bound multiplies the original whole-potential-hull
floor ONCE. Kp here is a scalar slow-response coefficient, not an ordinary
Haar primitive substituted for the full-torus pressure operator.

The original nondegenerate five-map, fixed compensator and reserved support
conditions, including EXACT G=0 on the ENTIRE potential hull, are analytical
premises. The code does not certify them from dictionary keys, successful
norm arithmetic or local bump support. They are unchanged from Q39 and
the frozen slow construction. These are slow INCREMENT bounds; old means
are not replaced by (7), and the incoming state object is not mutated.

## 5. Requests and next obligation

Postfast output n requires fast fields through n+1, Q101 rows through n+6,
changed covariance through n+8, Q94 mean-factory reserve n+11, conservative
target/profile n+25,n+38 and new cutoff shape through n+8. Slow output n
or (5) through n requires postfast n+1, hence covariance n+9, mean-factory
n+12, target/profile n+26,n+39 and new cutoff shape n+9. These counts cover
the newly composed calls; the original compact tails and nested high-jet
graphs remain unexpanded and conditional. No giant actual jet was evaluated.

Next derive the postslow reset using the CURRENT incoming Haar means,
which now contain both the original slow field and the first-real refreshed
slow increment. Reusing the first-real Q41.old_haar alone would omit that
previous slow increment. The new complete fast field has zero Haar, but
its full correlations have already been used in (3). Then recompute final
radial pressure, coefficients and complete mean/residual differences on
the new state. No final scalar reset or complete residual upgrade is claimed
here. Native/trace/support hypotheses, uniform actual pulse bounds,
viscous/cutoff/interaction residual closure, infinite realization and
blow-up remain unproved. Original fields, root, width C_setup^-20,
accepted domain and historical frozen evidence are unchanged.

Exact Fourier/polynomial controls subtract the full radial source including
time and vector viscosity and recover (4), including the pressure part of
Jz. They retain old/new products and detect freezing the moving pressure
factor. Manufactured norm leaves verify the current state, actual new fast
interface, next scalar in the radial response, exact pole and output reserves;
they are not divergence-free canonical continuum data. Run
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH.json`.
