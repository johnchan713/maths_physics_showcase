# Axial potential, pressure complement and the actual postwave A1 gain

Parent: `1999c84d2476b7b332a08bc108d158078f3a7f6c`.
Its angular-normal-form CI step passed when inspected on 23 October 2026;
the complete workflow was running without a failed step.

The parent derives order-one FINAL mean pressure from the final angular
row. Here we extend its pressure COMPLEMENT estimate to the SAME POSTWAVE
state before the latest fast/slow means, and use the fine current axial
row to bound the nonzero-Haar compact potential. This gives

    Qcirc Psi, Qcirc pm in M_1,
    A1 Qcirc Ez_postwave in M_2,                       (1)

under the retained current coarse field, derived fine row and whole
operator/support contracts. The actual A1 error improves from 2-kappa
to 2. The route in (1) needs NO additional sharper Haar-mean callback at
this new layer. Haar assumptions used upstream in the fine-row DAG
remain retained; they have not been independently validated. Separately,
the preserved finer Haar gamma history yields full Psi in M_1 and full
beta=-Dz Psi in M_2, instead of the old coarse 2-kappa radial-mean rate.
That history remains an actual-input obligation; it is not inferred
from the coarse field norms or fine current rows.

These are conditional analytical estimates, not validation of actual
all-scale norm, profile or support inequalities. The full regular source
stays E=2-3kappa and its physical exponent stays negative. No total
pressure estimate, unconditional pulse, complete residual improvement,
absolute decay, infinite field or blow-up is proved. Historical frozen
evidence, fields, roots, C_setup^-20, h, carriers and band are unchanged.

## 1. Start with the actual compact potential and old coarse pressure

At either ONE specified state, keep the finite sum of the actual fast
and slow potentials, including all A1 velocity corrections:

    gamma=D1 Psi, beta=-Dz Psi, T1 gamma=Psi.

Compact support here means the WHOLE potential hull, including filled
gaps. Divergence and velocity moments alone do not assert this actual
potential identification. Use the fixed-label operators of the parent,
and let

    Hrz=b gamma+G beta+beta gamma+Wrz,
    X=2G gamma+gamma^2+Wzz+pm.

pm is the COMPLETE same-state mean pressure, not merely its nonlocal
primitive or its last increment. The complete axial mean equation gives

    Ez=D1 Lambda+Dz X,
    Lambda=tstar Psi+Hrz-epsilon D_r gamma
                              -epsilon Dz^2 Psi-Sigma_z.       (2)

Both scalar axial viscosities remain. Since Sigma_z is independent of
the auxiliary torus, compact inversion on the full Qcirc Lambda gives

    T1 Qcirc(Ez-Dz X)=Qcirc Lambda,
    Qcirc Psi=c^-1 Nt^-1 A,
    A=T1 Qcirc(Ez-Dz X)+epsilon partial_T Qcirc Psi
              -Qcirc Hrz+epsilon D_r Qcirc gamma
                                      +epsilon Dz^2 Qcirc Psi. (3)

T1 is applied to Dz X as a WHOLE input. It has not been commuted through
Dz; derivatives of its moving cutoff remain in the original operator
recipe. Nt commutes with T1 on the fixed chart. No raw prefast row sum,
label renormalization, or multiplier-through-inverse step is used.
The complete difference Ez-Dz X has the required characteristic moment
because it is D1 Lambda; its two summands need not separately have zero
characteristic moments.

The RHS of (3) uses ONLY the OLD coarse bounds Psi=T1 gamma in M_B
and pm in M_B, B=1-kappa. No improved Psi or pressure is used to prove
itself. The source Ez is the already derived full same-state row:
E=2-3kappa POSTWAVE, C0=2-2kappa FINAL. The postwave recipe is old
fine full rows plus the complete new covariance row difference; it
precedes the latest means. All earlier source/A0/A1 dependencies remain.
Additive Gaussian/base/exterior residual tails are excluded from this
regular mean-row inverse and retain their assigned ledger.

From beta in M_(B+1), gamma in M_B, b=epsilon b0, W in M_1, and
the old pm in M_B, one has X in M_B and Hrz in M_1. The orders in A are

| Term | Order |
|---|---:|
| T1 Qcirc Ez | E postwave; C0 final |
| -T1 Dz Qcirc X | B+1 |
| epsilon partial_T Qcirc Psi | B+1 |
| -Qcirc Hrz | 1 |
| epsilon D_r Qcirc gamma | B+1-kappa=C0 |
| epsilon Dz^2 Qcirc Psi | B+3 |

All are at least one for the original 0<kappa<1/4. Original weighted
T1 and four-jet temporal bounds therefore prove Qcirc Psi in M_1.
`potential_envelope` records a positive norm program: output n reads
current fields/KT1 through n+6, Ez through n+4 and the OLD coarse pm
recipe through n+5. Its original conservative pressure recipe can read
fields through n+8. No new desired potential or fine pressure callback
replaces those lower programs. KT1 means the finite maximum of the
actual original constants at all requested orders on ONE domain, not
an assumed monotonicity of inverse norms.

## 2. Preserve the second viscous contribution carried by pressure

The parent's local pressure contains +epsilon D1 beta. Thus X in (2)
contains -epsilon D1 Dz Psi. If it is removed by defining

    Xflat=X-epsilon D1 beta,

the EXACT alternative is

    Ez=D1[tstar Psi+Hrz-epsilon D_r gamma
                          -2epsilon Dz^2 Psi-Sigma_z]+Dz Xflat. (4)

A single explicit Dz^2 Psi belongs with FULL X; the factor TWO belongs
with Xflat. Mixing them loses an axial viscous force. We use full X in
(2)--(3), and full pm below. Its derivative, local K0/K1, complete
covariance, moving rho P, Pi1/Pi2 and old pressure memory remain. The
exact controls detect this factor error and a dropped slow-time term.

## 3. The same POSTWAVE pressure complement needs no finer Haar means

The parent's angular normal form used Etheta_final in M_C0. Its proof
requires only Etheta at least E>1: the worst angular viscous remainder
already has order E. Hence it applies to the derived POSTWAVE angular
row at E, on its OWN 6q proof domain, with the old postwave means and
the complete updated covariance. FINAL pressure remains on 12q.

Use the parent's a=2V/R, a1=6V/R^2-2V_R/R, angular quantities Htheta,
Ktheta and its exact SAME-state pressure

    K1=K0-a Htheta,
    r1=a Haar(v)+v^2/R+Rrest-a1 Htheta-a Ktheta,
    pm=K1+T0(r1-rho P).

Haar commutes with T0, and rho P is auxiliary-independent even when
moving in Z,T. Therefore

    Qcirc pm=Qcirc K1+T0 Qcirc r1,
    Qcirc r1=Qcirc(v^2/R+Rrest)-a1 Htheta-a Ktheta.           (5)

The a Haar(v) term is absent ONLY from this complement. Both sources
in (5) have order one from the coarse means and the angular row alone.
Thus Qcirc pm is M_1 without any finer Haar(v) bound, moment norm or rho
norm. Full pressure still has the complete moving compensation and
variance; its Haar part has NOT been declared zero or order one here.
In particular Dz(rho P) is generally nonzero in the full axial equation.

`complement_envelope` computes only this channel. It needs neither a
finer Haar norm nor an incoming fine pressure estimate. The positive
formula is 2*(local+KT0*source) from the parent's angular program with
the a Haar(v) term suppressed. This bounds (5); its zero placeholders
are NOT estimates for full pressure or its moments. All other full-field
products are projected AFTER multiplication, with their factor two.

## 4. Keep G_R Psi before projecting the actual axial fast source

G is auxiliary-independent but radially variable. Exactly,

    2G gamma=D1(2G Psi)-2G_R Psi.

The complete compact support of G Psi, the commutation of D1 and Dz,
and A1 D1=0 imply

    A1 Qcirc Ez=A1 Dz Qcirc Phi,
    Phi=-2G_R Psi+gamma^2+Wzz+pm,
    Qcirc Phi=-2G_R Qcirc Psi+Qcirc(gamma^2+Wzz)+Qcirc pm. (6)

Equations (3),(5) prove Qcirc Phi in M_1 without sharper Haar data.
The gamma square has order 2B>1. G_R is retained on the WHOLE potential
hull; G=0 only on the reserved slow hull does not remove this global
term. Its R/U conversion G_R=d R^(d-1)G_U and all fixed jets remain.
The product is differentiated as a whole. A1 stays OUTSIDE Dz: a moving
cutoff has a nonzero [partial_Z,A1] on a general Phi. No division by
the cutoff at its zeros is used.

Apply (6) at the POSTWAVE state before the latest fast map. Its actual
fast-time identity is c Nt delta_gamma=-fz+A1 fz, fz=Qcirc Ez_postwave.
The positive A1 fz error, now in M_2, retains its actual A1 velocity in
all subsequent products. FINAL pressure or potentials are not inserted
into this earlier source. A FINAL version of (6) is a different identity
and is not a bound on the earlier fast source.

For raw coefficients the sufficient direct bound is

    N_m(A1 Dz Qcirc Phi)<=2epsilon KA1_m N_(m+1)(Qcirc Phi).

The extra factor two is conservative when a zero-Haar source is already
supplied. The original raw KA1, moving cutoff/profile constants and
weighted conversions remain. A factorial-through bound has the deferred
factor 2*6^m*(m+1)! KA1_m, not an arbitrary inverse norm. Thus the actual
postwave A1 error has order TWO, margin 2kappa over C0. This improves
the previous 2-kappa error, not the pointwise full axial row or regular E.

## 5. Full potential and radial mean still require the Haar history

Under the EXISTING finer same-field Haar gamma bound from the finite
preserved slow sum, compact inversion gives

    Haar Psi=T1 Haar gamma in M_1,
    Psi in M_1, beta=-Dz Psi in M_2.                         (7)

This keeps the whole moving potential before its derivative. The full
gamma remains M_B because D1 is a graph derivative. `beta_envelope`
uses the fine potential one order higher; output n consequently reads
fields through n+7, Ez through n+5 and coarse pm through n+6. Only
(7), not the A1 complement route, requires that sharper Haar history.

This is an essential extra input. On a reserved-hull slow plateau take
an auxiliary-independent increment delta_Psi=epsilon^B f(R), compact
nonzero f, delta_beta=delta_v=0 and delta_gamma=epsilon^B D1^0 f.
The divergence and both velocity moments are preserved. With G=0 and
b=epsilon b0(R), its homogeneous mean-row response is

    delta_gr=0, delta_Etheta=0,
    delta_Ez=epsilon^(B+1) D1^0[b0 D1^0 f-partial_R D1^0 f]. (8)

This is finer than both E and C0. Against coarse old means the extra
crosses are at least 2B+1-kappa>C0; a fixed smooth slow potential
extension also preserves fine response orders and exact divergence.
Yet epsilon^-1 Haar(delta_Psi)=epsilon^-kappa f grows exponentially
in the band index and cannot have a fixed S-polynomial bound. All
fixed jets of this manufactured increment obey the coarse class.
Thus fine current row response and the coarse constraints do not supply
the missing Haar-history inequality. This is a class-level homogeneous
comparison, not the actual prescribed slow map, covariance history or
a canonical-source counterexample. Its purpose is to keep (7)'s precise
additional gate visible, rather than promote it from the tests.
The comparison does not impose the full current scalar-defect constraints;
it does not rule out proving the missing Haar bound using those further
identities or the prescribed finite history.

## 6. Physical requests, common domains and remaining obligation

For the unchanged physical target t=(N+3/2+2h+2m)/h, or normalized
t=C0, let p=max(0,floor((t-2)/kappa)+1). The new source is Qcirc Phi,
through raw order m+1 if p=0, otherwise m+2p+4. At kappa=1/100000:

| Output m=0 request | New Qcirc Phi raw jet | Previous H1 raw jet |
|---|---:|---:|
| Normalized C0, direct p=0 | 1 | 1 |
| Physical h=1/200,N=0 | 60000006 | 60000008 |

These are DIFFERENT functions. Qcirc Phi at J reads current Ez/Etheta
through J+4, old coarse pressure through J+5, fields through J+6, and
the old conservative pressure field recipe through J+8. The last
physical illustration is 60000014, below the retained final A0 local
field request 60200007. That does not reduce a total original-profile
jet count: every earlier source/root/normal/A0/A1 graph remains intact.

Choose one original finite pole covering ALL new and old nested jets.
At postwave, old means have 2q, covariance 6q, and both the new potential
and pressure COMPLEMENT estimates can be embedded at 6q; fixed G_R
introduces no new pole under the retained unweighted multiplier contract.
FINAL estimates use the separately retained 12q operators. No monotonicity
between different-pole operator constants is inferred. Add any actual
fixed-multiplier pole before products, and retain all original KT1,
temporal, Kp, KA1, oriented primitive, raw profile, width/gap and edge
conversion costs at their actual orders. No giant factorial, selected
field or new derivative-dependent first band is evaluated.

The actual inequalities still needed are the SAME-stage weighted/raw
Ez, Etheta, Hrz, X and Qcirc Phi norms at these finite requests, all
whole compact potential/flux traces and fixed profile/G_R jets, and
the older coarse pressure/input identification. The sharper preserved
Haar gamma inequality is additionally needed for (7). Passing finite
fixtures supplies none of those canonical norm enclosures.

Next validate these reduced complement and support inputs before the
next generic signed cycle, retaining complete pressure and mean history.
The full regular order E, negative physical exponent, infinite diagonal,
localization/forcing regularity and velocity growth remain open.

Exact checks reconstruct the full axial equation, both viscosity forms,
variable-G projection and moving-cutoff commutators, and compare independent
epsilon budgets. Polynomial/Laurent fixtures have finite endpoint regularity
and are not canonical C-infinity fields or continuum PDE certificates.
Historical frozen evidence is bound without modification. Reproduce with
`python test_audit.py -v` and `python audit.py --verify-record evidence.json
--output FRESH_PATH.json`.
