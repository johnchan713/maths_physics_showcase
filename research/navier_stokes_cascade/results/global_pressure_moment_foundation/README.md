# Current same-root global pressure and moment foundation

Parent: `29b1f5f281650f5fae082d5d4a90a1f0d66e9ae3`.
Status: `current-same-root-global-pressure-moments-and-leading-exterior-stresses-joined-relative-to-retained-inputs`.

This checkpoint closes the next **current-source** obligation: it traces the
same accepted angular, heat, axial, annular and compact-repair roots through
the original pressure gauge, proves that the two transformed five-row systems
are equivalent to the original five moment equations, and joins the resulting
zero global defects to the signed exterior leading stresses and the already
audited inner/compact/outer leading-value cone.  It is relative to the retained
exact-core and compact-input estimates.  It is not an unconditional completed
profile, physical-pulse theorem, complete Navier--Stokes residual estimate or
blow-up proof.

No parameter, width, frequency, support, root, label, cutoff sequence or gauge
is changed.  In particular `lambda=exp(-4T)>0`, the repair frequency remains
`N=1+floor(H^32)`, and the activation width remains `C^-20`.  The new evidence
replays the latest matching-root checkpoint and the frozen global-moment and
regional-cone derivations instead of modifying their historical records.

## 1. The transformed roots really are the original five roots

Write the original functional defect vector in the order

\[
             d=(M,R,J,S,C_p),                              \tag{1}
\]

where `R` is the renormalized heat moment

\[
 R=\int_0^\infty(H-H_{pow})\,dX
  =\lim_{X\to\infty}\left(I(X)-\frac{XH_{pow}(X)}{1-h}\right). \tag{2}
\]

At a finite join the implemented annular row order, after dropping only
strictly positive normalizers, is

\[
 T_a d=(M,J-4\eta I,I,-(S-8\eta M),C_p).                  \tag{3}
\]

Its exact inverse is

\[
 (r_1,r_2,r_3,r_4,r_5)\mapsto
 (r_1,r_3,r_2+4\eta r_3,8\eta r_1-r_4,r_5).              \tag{4}
\]

The future heat tail and its amplitude are identical on both sides of this
compact join, so a zero endpoint `I` difference is exactly a zero global `R`
difference.  Equations (3)--(4) show that the transformed zero set is the
original five-row zero set.  This is a triangular identity, not a numerical
conditioning assertion.

On the compact repair the implemented row order is

\[
 T_c d=(M,(M-J)/\lambda,I,-S,C_p),                        \tag{5}
\]

and, because the selected `lambda` is strictly positive,

\[
 (r_1,r_2,r_3,r_4,r_5)\mapsto
 (r_1,r_3,r_1-\lambda r_2,-r_4,r_5).                     \tag{6}
\]

Thus the same compact root restores `M,R,J,S,Cp`, including their functional
eta dependence and all fixed one-sided eta jets.  The old warning at
`lambda=0` remains essential: (5) is then not a definition of two independent
`M,J` equations.  The latest checkpoint's `1/lambda` target loss is retained.

## 2. Complete same-gauge pressure trace

Let `E_id` be the unedited exterior schedule and

\[
 C_{p,id}(\eta)=\int_0^\infty\frac{E_{id}^2}{2X}\,dX,
 \qquad \Pi_{axis}=-C_{p,id}.                              \tag{7}
\]

This is the datum used by the exact core.  Starting there, the exact edit
ledger is:

| Exact edit | Total action |
|---|---|
| angular roots | preserve total `Cp`; restore the angular heat condition |
| axial roots | preserve `I,Cp`; set `M=J=0`; the strength root sets complete `S=0` |
| heat edit and compensation | opposite changes in `R,S,Cp` |
| axis/shape replacement and annular root | opposite changes in all five rows by (3)--(4) |
| compact modulation and repair | opposite changes in all five rows by (5)--(6) |

Every statement is an equality of functions of eta.  The heat edit and its
separate compensation can carry nonzero pressure memory between their
supports; that memory is never reset.  With the original gauge retained,

\[
 M_\infty=R=J_\infty=S_\infty=0,
 \qquad \Pi_{axis}+C_{p,\infty}=0.                         \tag{8}
\]

Consequently the completed leading swirl has the exact backward pressure

\[
 \Pi(X)=\Pi_{axis}+\int_0^X\frac{E^2}{2x}\,dx
       =-\int_X^\infty\frac{E^2}{2x}\,dx.                \tag{9}
\]

The fifth row is indispensable.  If it is omitted, the remaining functional
pressure defect `p(eta)` contributes

\[
 X\{4(1/2+h)\eta p-(1-\eta^2)p_\eta\}/(1-2h\eta^2)       \tag{10}
\]

to `r T_z`.  A tiny defect, or a defect whose value vanishes at one angle,
does not give (9).

## 3. Join to exterior stresses and the leading cone

For two profiles with identical exterior local fields and defect vector
`(m,R,j,s,p)`, put `D=1/2-h`, `A=1/2+h`, `d=1-eta^2`,
`L=1-2h eta^2`, and `B_m=2D eta m+d m_eta`.  The exact differences are

\[
 \Delta V_0=-B_m/L,                                      \tag{11}
\]
\[
 \Delta(r^2T_\theta)=
 \{HB_m+(1-h)R-D\eta R_\eta-dj_\eta+2(h-D)\eta j\}/L,   \tag{12}
\]
\[
 \Delta(rT_z)=\{D(m-\eta m_\eta)+4h\eta s-ds_\eta
       +X(4A\eta p-dp_\eta)\}/L.                        \tag{13}
\]

Substitution of (8) makes all three expressions zero identically.  Hence the
forward leading stresses equal the backward stresses on the exterior, and
the frozen outer-edge sign/orientation calculation applies to this same
leading profile.  The matching checkpoint identifies the annular and compact
roots with the old accepted roots, while the regional-cone checkpoint uses
the same `N`, parameters and endpoints.  Its leading-value bounds therefore
join the protected inner collar, compact modulation/repair interval and
heat-compensated outer collar without reselecting a root or frequency.

This completes the current compact-modulation/exact-repair pressure and
global-moment stream.  It does **not** prove that the full positive-order
coefficients share untouched reserved mean supports, or that every actual
physical pulse and derivative obeys the required uniform estimates.

## 4. Checks and remaining limits

`identities.py` independently checks (3)--(6) with exact rationals, tests the
signed pressure ledger and backward-tail equality, and copies (11)--(13)
independently before comparing its zero-defect consequence with the frozen
global-moment ledger.  Negative controls reject `lambda=0`, omit the pressure
repair, retain an eta pressure slope with zero value, and verify the resulting
linear-in-`X` axial stress defect.  These are algebraic checks; the analytical
existence input is the replayed same-root chain, not a numerical residual.

The next obligation is to discharge the remaining finite-background
coefficient and untouched reserved-mean support lemmas, then the uniform
actual-pulse existence/inverse/gap/value/derivative/support estimates at every
physical scale.  Complete wave-covariance mean corrections must still retain
axial velocity and pressure, and the final residual estimate must include
viscosity, cutoffs and every interaction.  No unconditional pulse bound,
complete residual improvement or blow-up is claimed here.

Reproduce with `python -m unittest -v test_audit.py` and
`python audit.py --verify-record evidence.json --output fresh.json` in this
directory.  The output path must be fresh; frozen evidence is never replaced.
