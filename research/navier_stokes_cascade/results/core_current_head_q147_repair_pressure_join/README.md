# Q147: same-modulation repair, cone join, and global pressure identity

Parent: Q146, `37aa64ffed83f55541bcb5f04eab7c806d419db1`.
The live branch, Q146 evidence and current CI, Q74, Q131/Q132, Q137/Q138,
the accepted repair neighborhood, regional cone join and global pressure
reconciliation were read first. Historical evidence, roots, supports,
pressure datum, frequency, scales and the `C^-20` width remain unchanged.

Q146 corrected a stage-only join: it did not identify Q138's modulation
with Q132's target memory. This checkpoint makes that identification. It
also resolves Q132's baseline data from the earlier exact match rather than
setting them to zero because a support ended. The same accepted repair then
joins the protected value cones and discharges the finite endpoint premise
in the global five-moment/pressure identity.

## 1. The baseline defects are zero functions for an exact reason

At the accepted earlier annular endpoint all five unmodulated primitive
differences vanish as functions of `eta`. The target profile and PREMOD
profile then use the same fields and hence the same five densities up to the
compact modulation. They retain the same `Pi_axis+Cp` gauge. Therefore each
Q132 baseline function satisfies

    d_r^base(eta)=0,  r in {M,I,J,S,Pi},                 (1)

and every angular derivative is zero. Thus valid Q132 pairs are
`(T_r,L_r)=(0,1)`. This is propagation of an exact functional identity from
the five-row annular match. It is not inferred from vanishing local velocity,
the end of a bump, a sampled angle, or a pressure reset.

The fixed-origin loop primitives vanish on the preserved entrance collar.
Consequently the modulated-minus-PREMOD moment differences themselves start
at zero there, although the incoming moments need not be zero.

## 2. Q138 and Q74 give the same five Q131 memory rows

Q146 supplies the same-profile complex caps needed to run the already derived
Q140--Q139--Q138 chain. Write its actual modulation as

    E_N=E exp(Acal/N),  U_N=U+Bcal/N.                    (2)

This is exactly Q74's field pair: the same selected integer, moving phase,
fixed-origin density-weighted means and zero-mean primitives. With bounds
`e,u,DE,DU` for `E,U,E_N-E,U_N-U` on `[x0,x1]`, direct expansion of the five
signed densities gives

    D_M=DU,
    D_I=r DE,
    D_J=r(e DU+u DE+DE DU),
    D_S=2u DU+DU^2+e DE+DE^2/2,
    D_Pi=(e DE+DE^2/2)/x0,                              (3)

where `r=max(1,2x1)`. Integrating from the zero entrance trace yields the
Q131 prefactors `(x1-x0)D_r`. Products use twice Q138's angular base. The
shifted derivative bound retains `Pi_eta`; no pressure memory is discarded.

Q74 independently expands the exact terminal increments through quadratic
order and, more decisively, defines the exact target from the full increment
`Delta`. In ordinary normalized rows `(m,i,j,s,cp)`, it is

    z=-T_c Delta=(-m,(j-m)/lambda,-i,s,-cp).             (4)

This is the signed negative of Q132's
`(m,(m-j)/lambda,i,-s,cp)`. Hence Q132's target is now bound to the same
Q138 object, with (1) as its complete baseline and (3) as its full incoming
memory. The inverse-`lambda` loss, pressure square, energy squares and all
five rows remain.

## 3. The original repair and value-cone join now apply

The low-order accepted-root proof used this same exact target. Its selected
finite `N=1+floor(H^32)`, inverse bound, quadratic bound and target estimate
put the root in the old uniqueness ball and inside the `A^-128/2` coefficient
radius. The repair-state estimate then preserves the positive swirl
neighborhood. Since (4), rather than a stage-compatible surrogate, is the
target, Q137 is the physical derivative envelope of those same coefficients.

At the repair exit all five defects vanish as functions of `eta`. The compact
gap remains at least `G^-4/2`, so the existing protected inner collar and
annular/exterior value estimates join to it with unchanged endpoints. This
closes the leading compact modulation/exact-repair/value-cone obligation.
It does not provide uniform physical pulse derivatives.

## 4. Global five-moment and pressure identities

After the repair exit the compared profiles have identical fields, including
the already audited common heat tail. Finite endpoint equality therefore
propagates to the global renormalized identities

    M_infty=R_infty=J_infty=S_infty=0,
    Pi_axis+Cp_infty=0.                                  (5)

With the unchanged axis datum, the completed leading pressure is

    Pi(X)=-integral_X^infty E(x)^2/(2x) dx.               (6)

Substitution of (5) into the complete exterior difference formulas makes the
forward and backward stresses agree. The fifth identity is independent: a
remaining `p(eta)` would contribute
`4 Aphys eta p-(1-eta^2)p_eta` to the axial numerator.

Equations (1)--(6) are analytical consequences of the cited continuum
identities. Exact-arithmetic tests check the density expansions, row signs,
normalizations, inverse-`lambda` term, pressure row, profile-cap binding and
frozen provenance. They are diagnostics, not a continuum proof by testing.

The next obligation is the finite positive-order coefficient and reserved-
mean support ledger on this finalized profile. Uniform actual-pulse size,
derivative and support bounds at every physical scale, complete covariance
means including axial velocity and pressure, the full viscosity/cutoff/all-
interaction residual and blow-up remain unproved.
