# The same original tilt: raw-deficit radius and the joint freezing/scale gate

Status: `same-tilt-raw-deficit-radius-and-joint-freezing-scale-gate`.
Parent: `abfaf6c133a1b24a88c3eb55f54ac0546178a770`.

This checkpoint derives an explicit value-domain certificate directly from
the ORIGINAL fixed tilt's raw-cone deficit. It bypasses the earlier stronger
coarse sufficient tilt test, not the need to quantify the deficit for a
retained scale. It also makes the independent target-freezing cost explicit
in that same deficit. No tilt, carrier, cutoff, coefficient, root, support,
frequency, gauge, pressure datum or mean correction is reselected.

## 1. What the original tilt choice already guarantees

Use the complete leading stress direction n0 on the CLOSED active annulus,
with its smooth limiting directions at both flat edges. For the pointwise
leading shear frame N0,K0 and mixing parameter c0, define

\[
 H_0=-N_0\cdot n_0\ge\omega>0,\qquad
 R_{act}=\sup\frac{c_0^2(K_0\cdot n_0)^2}{H_0^2},\qquad
 \kappa_u=\frac{u^2}{1+u^2}.
\]

The manuscript's printed page 74, immediately before (7.2), fixes ONE u>0
with `sqrt(kappa_u)>sqrt(Ract)`, including the edge interpretation. The
retained actual-band and leading-margin checkpoints reproduce that choice.
Thus, conditional on that original admissibility and the retained leading
lemmas,

\[
 d_u:=\kappa_u-R_{act}>0.                                 \tag{1}
\]

This is already a qualitative positive SAME-TILT gap. It would be wrong to
infer that it disappears when the later coarse sufficient test fails.
Conversely, original strictness gives no supplied quantitative lower value
for d_u. That missing value still matters at the retained first scale.

The earlier explicit bound was Ract<=1-epsilon*. Its test
`(1+u^-2)(1-epsilon*)<1` is equivalent to
`kappa_u-(1-epsilon*)>0`. It is sufficient, not necessary. If it passes,
that raw difference is a valid lower bound on d_u. If it fails, the actual
original choice can remain admissible, with a smaller actual Ract. We do
not claim that the retained u satisfies the coarse test. This is a
refinement of that checkpoint's stated limitation, not an alteration of
its frozen evidence.

Band shear/swirl factors are common positive factors, so Ract, kappa_u and
d_u are unchanged by the quarter-to-four chart conversion. The magnitude
of the leading stress is not assumed nonzero at its flat edges, and no
pressure term is dropped from its direction.

## 2. A normalized cone polynomial with an explicit deficit radius

Keep the previous chart input floors
`f<=F0<=Fbar`, `-g0theta>=t>0`, `|g0|<=Gbar`,
`D0=|g0|^2+2F0*g0theta>=d>0`, with all band factors included.
Let e bound all three perturbations

\[
 \max(|F-F_0|,|g-g_0|,|r|)\le e,\qquad
 r=(\widehat T-T_0)/|T_0|.
\]

The quotient is the inherited smooth factored extension, NOT division of
numerical zero stresses. Reuse the full background-gap constants

\[
\begin{split}
 L_D&=4\bar G+2\bar F+3,&L_P&=\bar F+\bar G+1,\\
 \bar D&=\bar G^2+2\bar F\bar G,&
 \bar c^2&=\bar D/(2ft),\\
 L_c&=2L_D/(ft)+2\bar D L_P/(f^2t^2),&L_H&=1+2/t,\\
 r_b&=\min(1/2,f/2,t/2,d/(2L_D),\omega/(2L_H)).
\end{split}
\]

For the perturbed shear frame N,K put `H=-N.(n0+r)`, `J=K.(n0+r)`.
Instead of the unnormalized cone functional use

\[
 \Psi=\kappa_u H^2-c^2J^2,
 \qquad L_\Psi=3\kappa_u L_H+4L_c+3\bar c^2L_H.             \tag{2}
\]

This is exactly the old Phi divided by 1+u^2, so
`Lpsi=Lphi/(1+u^2)`. No pulse/frame is renormalized by doing this algebra.
The old discriminant, quotient and unit-vector difference arguments give,
for e<=rb,

\[
 |\Psi-\Psi_0|\le L_\Psi e,\quad
 \Psi_0=H_0^2(\kappa_u-\rho^2)\ge\omega^2 d_0,
\]

whenever `0<d0<=d_u` is a VALID raw-deficit lower bound. Therefore

\[
 \boxed{r_d=\min(r_b,\omega^2d_0/(2L_\Psi)),\quad e\le r_d}
 \quad\Longrightarrow\quad \Psi\ge\omega^2d_0/2>0.          \tag{3}
\]

All spectral/swirl floors in the original gap certificate remain valid:
F>=f/2, -gtheta>=t/2, D>=d/2, and
`lambda^2>=f*t*d/[4(Gbar+1)^2]`,
`c^2>=d/[4(Fbar+1)(Gbar+1)]`. No square root is differentiated to prove
the cone estimate. The same-point frame here is not a replacement for the
subsequent frozen representative or moving wave frame.

Since |n0+r|<=3/2 and H>=omega/2, the normalized stress direction obeys
`-N.nB>=omega/3`, with magnitude >=|T0|/2 in the active interior. There
is no nonzero stress floor asserted at its zero edges. Also H^2<=9/4, so

\[
 \frac{c^2J^2}{\kappa_u H^2}
 \le1-\frac{2\omega^2d_0}{9\kappa_u},\qquad
 \boxed{\gamma_B=\frac{\omega^2d_0}{9\kappa_u}>0}.           \tag{4}
\]

Taking a square root uses sqrt(1-x)<=1-x/2; normalizing the target cancels
from the ratio. The formula is safe since 0<d0<=kappa_u and omega<=1.
The old route gamma=d0/(2kappa_u) would spend part of the gap first:
its normalized starting floor is
`omega^2[d0-d0^2/(4kappa_u)]`. Thus the direct cone-radius term in (3)
is strictly larger, and the full radius is never smaller, for identical
supplied d0 and geometry. This comparison does not certify the actual d0.

## 3. Keep the separate representative-freezing gate

The primary-covariance/source-envelope foundations require a valid
`tau_freeze<=Cfreeze*S^-3`, `S=ell^2`, for the SAME relative covariance
coordinates. This is not automatically supplied by a pointwise cone
estimate. Feeding (4) into their retained padded half-margin test gives

\[
 \ell^6\ge8C_{freeze}/\gamma_B
       =\frac{72\kappa_u C_{freeze}}{\omega^2d_0},\qquad
 \gamma_{frozen}=\gamma_B/2.                              \tag{5}
\]

The actual Cfreeze must still be bounded from the inherited smooth limiting
direction and closed-box derivatives. We do not identify it with an
unrelated constant or claim it was numerically evaluated. Every pulse
length/inverse budget that depends on a frozen margin must consume
gamma_frozen, not an unrelated larger margin. Geometry, carrier/normal,
length, source and nonlinear-cycle acceptance gates remain independent.

## 4. Precise joint deficit/shape/scale inequality

Keep the complete corrected-width first-stress constant Cstar(M,N), with
the actual cutoff's valid derivative enclosures, all pressure/source and
shifted viscosity terms, and the shared Borel half-budgets. On a whole
`q<=4Q`, `Q=2^-ell` chart, put

\[
 E=C_\star(M,N)2^{-2h(\ell-2)},\qquad
 d_{req}=\max\left(\frac{2L_\Psi E}{\omega^2},
                 \frac{72\kappa_u C_{freeze}}{\omega^2\ell^6}\right).
                                                               \tag{6}
\]

A sufficient joint background-cone and padded-freezing condition is

\[
 \boxed{ E\le r_b,\qquad R_{act}+d_{req}\le\kappa_u.}       \tag{7}
\]

This is the precise QUANTITATIVE inequality still blocking acceptance of
the retained scale. Bare original admissibility states only Ract<kappa_u;
it does not compare its distance from the boundary with dreq. If dreq is
larger than kappa_u, no nonnegative raw supremum can pass this bound.
Failure is failure of this sufficient estimate, not of the actual PDE cone.

At a supplied d0, the sharp cutoff-shape program receives
`B=rd*2^(2h(ell-2))`, while (5) is checked SEPARATELY. Passing its auxiliary
shape-existence program still does not identify the unchanged kappa. All
rows and their full-memory dependence remain exactly those of the preceding
two checkpoints; no repair or pressure datum is reset.

For rational diagnostic inputs the implementation uses
`m=floor(2h(ell-2))`, `Eupper=Cstar/2^m>=E`, avoiding fractional-power
roundoff. Replacing E by this upper bound in (6)-(7) is conservative, not
an exact equivalence at a nonintegral exponent. The +2 cannot be dropped.

A conditional later-band witness for valid supplied M,N,d0,Cfreeze is

\[
 \ell\ge\max\left(2+\left\lceil
 \frac{\max(0,\lceil\log_2(C_\star/r_d)\rceil)}{2h}
 \right\rceil,
 \left\lceil\left(\frac{72\kappa_u C_{freeze}}{\omega^2d_0}
       \right)^{1/6}\right\rceil,2\right).                 \tag{8}
\]

The exact sixth-root ceiling is computed by integer comparisons. Later
bands preserve these two gates by monotonicity. This is NOT a new choice
of the old first band or a complete pulse/residual acceptance threshold.
For the actual enormous hierarchy (1)-(8) remain symbolic; no actual
profile or constant is numerically materialized.

## 5. Why strict admissibility alone has no uniform scale budget

An exact algebraic countercontrol keeps F=1, g=(-82/25,0), c^2=16/25,
D=2624/625>4, u=3/4 and kappa_u=9/25 FIXED. With
`co=(1-z^2)/(1+z^2)`, `si=2z/(1+z^2)`, take the unit target

\[
 n(z)=((4co+3si)/5,(-3co+4si)/5).
\]

At z=0 it lies on the cone boundary. For 0<z<=1/32, n(z) is strictly
inside and n(-z) strictly outside. Factoring the normalized polynomial gives

\[
 \Psi(n(z))=\frac{2z[24(1-z^2)-14z]}{25(1+z^2)^2}>0.
\]

Here H(z)>=4/5 and `0<d_u(z)=Psi/H^2<=3z`. Consequently the strictly
admissible deficit and radius in (3) tend to zero while all those spectral
data and the tilt stay fixed. Yet `|n(z)-n(-z)|=4z/(1+z^2)` also tends to
zero and crosses the cone. This proves why no positive uniform radius
follows just from the strict original tilt relation. These are manufactured
constant algebra data, NOT counterexamples to the actual leading profile
or its pressure identities. No actual stress direction is replaced.

## Evidence and remaining obligation

Ninety-six exact rational perturbations check the normalized polynomial,
Lipschitz estimate, spectral floors, unit-direction floor and same-u margin.
Three comparisons check exact agreement of geometry constants with the old
certificate and the improved cone budget. Four thin-margin controls include
arbitrarily small exact rational gaps; symbolic identities and a complete
full-pressure/viscosity, sharp-shape and independent-freezing composition are
also checked. These are algebra diagnostics, not continuum proofs of the
retained Navier-Stokes premises. Historical source bytes and original
evidence comparators are preserved.

Next supply the unchanged cutoff's valid M,N, an actual quantitative d0
(or a sharper upper bound for Ract), and the relevant Cfreeze, then verify
(7) plus all other retained acceptance gates at the existing scale. The
original choice supplies qualitative positivity only. Uniform physical
pulse size/derivative/support/inverse-gap estimates, complete covariance
means including axial velocity and pressure, and the full viscous/cutoff/
interaction residual remain open. No unconditional pulse estimate, complete
residual improvement or blow-up follows.

Source: manuscript (7.1), printed page 74 before (7.2), and retained
background-gap/leading-margin, primary-covariance and source-envelope
derivations. [Manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its proposed final theorem is not accepted as established.

```sh
python3 research/navier_stokes_cascade/results/same_tilt_deficit/test_audit.py -v
python3 research/navier_stokes_cascade/results/same_tilt_deficit/audit.py \
  --verify-record research/navier_stokes_cascade/results/same_tilt_deficit/evidence.json \
  --output /tmp/fresh_same_tilt_deficit.json
```
