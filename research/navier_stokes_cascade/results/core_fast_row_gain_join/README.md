# Join actual fast-row cancellation to the fresh scalar gain

Parent: `7228cde1b0281e46054e381a2e708080b1b346b9`.
The parent derived the postfast scalar-C and fresh slow-C bounds from the
existing coarse mean recipes. Here the missing **post-wave zero-Haar
fast-row D bound** follows from exact cancellation in the preceding fast
map, its complete pressure/viscosity differences, the intervening slow
increment, and the generated complete covariance change. The same weighted
fast map then gives an actual fresh fast increment at D. With the parent's
slow-C estimate this proves its sufficient mean-feedback inequality (8)
and the conditional NEXT REGULAR nonzero source order D.

This is relative to the retained original actual coefficient, all-fixed-order
profile, support and exact-map contracts. It does not independently validate
them, prove all angular-zero equations small, bound the entire additive tail
ledger, or prove complete residual improvement, an infinite field or blow-up.
Historical frozen evidence, fields, width C^-20 and accepted band are unchanged.

The abstract postfast remainder orders already appear in
`current_target_cycle_bridge` section 4. New here is their explicit full
three-component pressure-difference polynomial applied at the current coarse
B, the composition with THIS generated prefix, and the very large finite A1
dependency needed to make that implication quantitative. No fine old-row
norm, new row-rate callback or new accepted-band gate is assumed.

## 1. Subtract the full equations before estimating

Set B=1-kappa, C=B+1/2, D=C-kappa and H=B+1-2kappa=2-3kappa,
0<kappa<1/4. Write old complete means as (x,u,z), the actual increment
as (X,U,Z), and the same full base as (b,V,G). All tangential means have
the retained coarse order B and radial means B+1; b has an extra epsilon.
The covariance is fixed during a mean update. Put

    dSrr=2bX+2xX+X^2,
    dSrz=bZ+GX+xZ+Xz+XZ,
    dStt=2VU+2uU+U^2,
    dgr=-t*X-D1 dSrr-Dz dSrz+dStt/R
                                      +epsilon(Delta0-R^-2)X,
    dP=int bar(dgr), dp=T0(dgr-rho dP),
    dRt=bU+VX+xU+Xu+XU,
    dAt=GU+VZ+zU+Zu+ZU,
    dAz=2GZ+2zZ+Z^2+dp.                                (1)

The actual compact pressure operator T0 is linear. Thus dp is exactly the
difference of the SAME pressure reconstructions. It retains dP and the
moving compensator rho. Neither a local pressure truncation nor resetting
old means to zero is used. The fixed Sigma terms and full Pi1/Pi2 cancel
only because both equations subtract the identical complete base.

With t*=-epsilon partial_T+c Nt, Dr=partial_R+graph derivative,
Dz=epsilon partial_Z and De=Dr+e/R, exact tangential subtraction is

    dEtheta=c Nt U+Rtheta,
    Rtheta=-epsilon partial_T U+D2 dRt+Dz dAt
                                      -epsilon(Delta0-R^-2)U,
    dEz=c Nt Z+Rz,
    Rz=-epsilon partial_T Z+D1 dSrz+Dz dAz-epsilon Delta0 Z. (2)

Every old/new mean cross, new square, cylindrical connection, radial/time
pressure source and viscosity term is present. The wave/base are held fixed
ONLY in this difference; the next wave change is included separately below.
Actual A1 velocities belong to (X,U,Z) and remain in every product.

## 2. An explicit remainder polynomial at H

Use through-n factorial coefficient norms in the same five variables; fields
and fixed data through n+3 suffice. Put d=n+4, dr=g d,
time=d(1+c), vis=d^2(g^2+1)+r dr+r^2, where r bounds R^-1.
Form the positive versions of the products in (1). Sufficient polynomials are

    Gr=time X+(dr+r)dSrr+d dSrz+r dStt+vis X,
    Pd=Kmom Gr, Pm=Kp(Gr+rho Pd),
    Rt=bU+VX+xU+Xu+XU,
    At=GU+VZ+zU+Zu+ZU,
    Az=2GZ+2zZ+Z^2+Pm,
    Ftheta=dU+(dr+2r)Rt+d At+vis U,
    Fz=dZ+(dr+r)dSrz+d Az+vis Z.                       (3)

Gr/Pm are bounded at the coarse order B, not H. Their axial pressure
derivative in (2) gains an epsilon. The epsilon powers accompanying (3) give

| Group | Sufficient order |
|---|---:|
| Slow time of the increment | B+1 |
| Linear radial transport | B+1-kappa |
| Tangential viscosity, including two graph derivatives | H |
| Mean products with a radial component | 2B+1-kappa |
| Full pressure gradient | B+1 |

Thus both remainders in (2) are at H. This is a subtraction estimate: applying
an absolute norm bound separately to Eold and Enew would lose the cancellation.
It uses only the already-derived coarse old field and actual increment bounds.

For the old FAST step let ftheta=Q1 Etheta,old and fz=Q1 Ez,old,
Q1=I-Haar. The retained exact operator identities are

    c Nt U=-ftheta, c Nt Z=-fz+A1 fz.

Hence Q1 Eafterfast=(Q1 Rtheta, Q1 Rz+A1 fz). A1 fz is generally nonzero.
For the following old SLOW step, c Nt(U,Z)=0 exactly. Its complete row change
is again bounded by (3) at H, even using only its coarse order B. This does
not require any finer old Haar-row theorem, scalar reset rate, or G=0 support
assumption. The latter remains necessary for the parent's slow scalar gain.

## 3. Keep the exact A1 term and its finite derivative cost

For any fixed output order n choose

    q=max(1,floor((H-B)/kappa)+1), J=n+2q+3.             (4)

Then B+q kappa>H. The ORIGINAL compact-kernel formula (8) of
`core_fast_inverse_budget`, with Bn=N_n(k_e chi_U), AJ=N_J(f_e), gives

    N_n(A1 fz) <= 2^(n+J)L Bn AJ [8*3^q*2^(n+q+3)]
                                      *M^-q N_J(fz).    (5)

There is no additional temporal inverse or c^-1 in this error: commute c Nt
with the actual A1 operator before estimating. The zero-Haar input excludes
the nonflat zero-mode radial-moment obstruction. Full smooth compact extension
and the original all-fixed-order coefficient/profile estimates are essential.

Use M^-q<c_i^(-rho_g q)epsilon^(q kappa). If I(S) bounds c_i^-1, the
positive polynomial (1+I)^ceil(rho_g q) bounds that fixed power on the SAME
band. The full old prefast source gives N_J(fz)<=2 J! Eold,z,J(S), after
unweighting its finite pole by Eedge(2p,width). A factor 6^n safely converts
raw output jets to the five-variable factorial through norm. Restore the
output edge weight with the original fixed chi_U support floor; a sufficient
inverse floor is 3^ceil((width^2+4)/gap^2). These are retained fixed geometry
constants, not a new shell or rectangular floor near a moving edge.

`flat_recipe` records precisely the resulting finite positive product:

    2*6^n*8 * 2^(2n+J+q+3) * 3^q * J! * L * Bn * AJ
      *Eold,z,J *Eedge(2p,width)*output_floor_inverse
      *(1+I)^ceil(rho_g q).                             (6)

Bn and AJ are the original RAW kernel/profile seminorm polynomials. An
arbitrary operator norm KA1 does NOT bound their product; no such inference
is made. Those profile recipes and the old source recipe are already retained
lower inputs. Formula (6) is not a new outgoing row estimate assumed as input.
For each fixed requested output, choose the common original finite pole p
to cover ALL required jets, including J and its old-source dependencies.
This does not assert that a pole chosen for low-order data still covers those
high jets. The corresponding fixed higher-pole constants must be used.

At kappa=10^-5, q=99999 and J=n+200001. In particular row n=0 needs the
old prefast source through 200001 (and the old initial-mean inputs in that
source recipe through 200004). A fresh fast field of order zero reads rows
through 5, hence old prefast source through 200006 and old initial means
through 200009. The order-zero next source, using fast fields through 3,
can read old prefast source through 200009 and old initial means through
200012. These are finite coefficient-jet requirements, not physical-scale
exponents or a claim that those actual high jets have been independently enclosed.

The implementation keeps (6) as an exact deferred polynomial recipe. It
does not allocate the giant factorial or invent its coefficients. Moderate
instances are checked against the independent raw constant in (5). A guarded
materializer and improved-inverse callback prevent silently using an
uninstantiated recipe as a numerical norm. Materialization is an execution
choice, not a new analytical condition or a need for an infinite numerical list.

## 4. Add the actual sourced covariance and join the mean feedback

The complete sourced covariance CHANGE has order C by the preceding generated
inverse estimate, with all old/new crosses, squares, quadratures and cutoff curls.
For a common entry majorant Dcov through n+3, sufficient row-change polynomials are

    Gw=(dr+d+2r)Dcov, Pw=Kmom Gw, pmw=Kp(Gw+rho Pw),
    Ew,theta=(dr+d+2r)Dcov,
    Ew,z=(dr+d+r)Dcov+d pmw.                            (7)

These retain Drtheta, Dthetaz, Drz and Dzz as well as the radial pressure
entries. The full pressure in pmw includes the global compensated moment.
The largest graph loss is one kappa, so (7) has order C-kappa=D. Since
H-D=1/2-kappa>0, (3), the intervening old slow change, (6) and (7) give
the desired actual post-wave Q1 row bound at D. Q1 has norm at most two;
`regular_rows` keeps the summands and that projection factor explicitly.

Old mean products have pole at most 4p and complete sourced covariance at
most 6p. Their common fast-row pole is 6p. The original 2p/4p and newly
covered 6p/12p pressure/moment recipes are used on their proper domains.
The code uses the explicit larger embedding max(1,L)^ceil(12p) when joining
lower poles, rather than assuming operator constants increase with pole.
The same weighted fast operator at input pole 6p now gives actual tangential
increment D and radial increment D+1. Its A1 VELOCITY is retained by the
same (1+KA1) factor. The parent's independently derived slow increment is C.

Both mean/current-wave transport orientations therefore have order
min(D,C)=D. This proves the parent's precise mean-feedback inequality (8)
relative to the retained premises. All other wave-stage groups were already
bounded at D in `core_multiharmonic_mean_residual_budget`; their finite
positive-polynomial sum gives the NEXT REGULAR source at D, with pole 9p.
No arbitrary additive tail is put into that source or any sourced inverse.
Actual A1 velocities remain in the complete field products, while (5) bounds
their exact zero-angular fast-time remainder. The old exact temporal cutoff,
A0/A1 and base residual ledgers remain separate and unchanged.

This is a finite-stage conditional regular gain, not complete residual
improvement. Full angular-zero Haar targets/rows, actual input/support and
same-primary/frame validity, all tail estimates, physical derivative/cutoff
joins, forcing and infinite-field/velocity-growth arguments still need their
complete proofs. The current primitive/pulse enclosures are not made true by
the very high finite derivative count or by passing tests.

## 5. Obligation ledger and reproduction

Removed relative to the existing all-fixed-order actual-map/profile contracts:
the extra post-wave fast-row-D hypothesis, its fresh fast-D consequence, and
the mean-feedback-D gap for this finite regular-source join. Added analytical
assumptions or accepted-band gates: none. Explicitly exposed: (6)'s high-jet
and original raw-kernel dependencies. No historical evidence or actual field
is modified and no full residual, unconditional pulse or blow-up is certified.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
Controls check exact full nonlinear row subtraction, divergence, pressure and
viscosity omissions, independent compact A1 signs/endpoint phases, strict
whole-range exponent margins and moderate raw-kernel constants. The polynomial
pressure surrogate in the subtraction control is explicitly not an evaluation
of the selected graph inverse. These are algebra/provenance diagnostics; the
continuum implication is the derivation above with its stated retained premises.
