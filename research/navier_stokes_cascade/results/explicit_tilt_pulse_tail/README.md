# Explicit tilt witness, raw deficit, and a certified primary-pulse tail

Parent: `9ff05e276d247ec23c15f08f22906234d7037c5a`.
Status: `explicit-tilt-raw-deficit-and-primary-pulse-tail-certificate`.

This checkpoint corrects a quantifier inherited by several scale-gate notes.
The construction never froze a numerical tilt `u` or a numerical first band
`ell0`: it said to choose them after finitely many strict inequalities were
known.  Consequently there is no historical integer whose acceptance can be
recovered from the frozen records.  It would be invalid to mark that absent
integer as verified.

There is, however, a legitimate constructive resolution.  We now specialize
the previously existential tilt using the proved global leading-cone bound,
obtain an explicit positive raw deficit, and give an exact integer formula for
a **separate certified later tail**.  This does not identify the specialization
with an earlier physical field or localized force.  The same core, coefficient
cutoff, `C^-20` width, carriers, supports, pressure gauges and full pressure
memory are retained.  Historical evidence is unchanged.

The result is analytical relative to the retained core/cone/frame/pulse
foundations.  Exact-arithmetic tests reproduce the parameter algebra; they are
not continuum supremum computations.  Complete generated source classes,
complete-wave covariance means, the nonlinear residual and blow-up remain
open.

## 1. The corrected first-band statement

The [leading-cone checkpoint](../leading_cone_margin/README.md) proves on the
closed active annulus

\[
 H_0\ge\omega_*={1\over640G^{13}},\qquad
 R_{act}:=\sup c_0^2{(K\cdot n_0)^2\over H_0^2}
 \le1-\varepsilon_*,\qquad
 \varepsilon_*={1\over2^{18}G^{22}}.                 \tag{1}
\]

The source then asks for one `u>0` with
`kappa_u=u^2/(1+u^2)>Ract`; it does not assign a value.  Likewise the pulse
foundation asks for one sufficiently large `ell0`; it does not assign an
integer.  Later notes calling these the “unchanged historical tilt/first
band” correctly avoided altering them, but the phrase obscured that the
parameters were still existential.  Qualitative strictness alone cannot
certify an unknown prior `u`, as the thin-margin countercontrols already show.

We therefore make two statements separately:

1. an unspecified earlier `u,ell0` remains unverified;
2. the explicit specialization below is certified and may seed a new tail
   theorem once its localization and force are propagated.

No inference from eventual acceptance back to an unknown first integer is
made.

## 2. An exact same-cone tilt and raw deficit

For any `0<epsilon<=1` satisfying `Ract<=1-epsilon`, choose

\[
 u^2={2\over\varepsilon},\qquad
 \kappa_u={2\over2+\varepsilon},\qquad
 d_0={\varepsilon\over2}.                              \tag{2}
\]

Then

\[
 \kappa_u-(1-\varepsilon)
 ={\varepsilon(1+\varepsilon)\over2+\varepsilon}
 =d_0+{\varepsilon^2\over2(2+\varepsilon)}>d_0.         \tag{3}
\]

Thus (1) makes `d0` a valid quantitative lower bound for the SAME cone
functional.  Equivalently,

\[
 (1+u^{-2})(1-\varepsilon)<1.
\]

For the actual hierarchy this specializes to

\[
 u^2=2^{19}G^{22},\qquad d_0={1\over2^{19}G^{22}}.       \tag{4}
\]

This is not a numerical approximation to the enormous hierarchy.  It is an
exact symbolic selection from the already proved bound.  The mixing factor
is large but finite.  The leading-frame freezing cost remains finite because
its only apparent large-`u` ratio is `(1+u)/u<=2` here.  Full axial and
pressure components of the target remain in `Ract` and the frame coordinates.

Let

\[
 \gamma_B={\omega_*^2d_0\over9\kappa_u}>0.              \tag{5}
\]

No positive stress-amplitude floor is asserted at either flat edge; (1)--(5)
are direction and normalized-cone statements using the inherited smooth edge
limits.

## 3. Exact selected-cutoff tail threshold

Use the explicit coefficient cutoff already selected in the construction,
with `M=729`, `N=2916`, and write

\[
 C_{coeff}=C_\star(729,2916).
\]

This contains the finalized `K0,5`, `K1,3*`, `H2,3*`, complete pressure
source moments, both shifted axial-viscosity moments and all repair/cutoff
jets.  No optimal hypothetical cutoff replaces it.

For positive retained constants `rb,Lpsi,Cmesh,Cpert,Dbox,dchart`, define

\[
 B_0=\min\left\{{r_b\over8},
 {\omega_*^2d_0\over16L_\Psi},
 {\gamma_B\over128C_{pert}}\right\},\quad
 m_0=\left\lceil\log_2{C_{coeff}\over B_0}\right\rceil_+,
                                                               \tag{6}
\]

and

\[
 \ell_{exp}=2+\left\lceil{m_0\over2h}\right\rceil,
 \quad
 \ell_{mesh}=\left\lceil
 \max\left\{1,{16C_{mesh}\over\gamma_B},
 {D_{box}\over d_{chart}}\right\}^{1/6}\right\rceil,
 \quad \ell_{tail}=\max(2,\ell_{exp},\ell_{mesh}).       \tag{7}
\]

All ceilings in the implementation use integer/Fraction comparisons.  Put
`Sfloor=2^floor(2h(ell-2))`.  For every integer `ell>=ell_tail`,

\[
 E={C_{coeff}\over S_{floor}}\le B_0,
 \qquad C_{mesh}\ell^{-6}\le\gamma_B/16,
 \qquad H_\ell={\gamma_B\over8}-C_{mesh}\ell^{-6}
 \ge\gamma_B/16.                                        \tag{8}
\]

Consequently

\[
 C_{coeff}\le S_{floor}\min\left\{{r_b\over8},
 {\omega_*^2d_0\over16L_\Psi},{H_\ell\over8C_{pert}}\right\},
 \qquad \ell^6\ge D_{box}/d_{chart}.                    \tag{9}
\]

This is the full selected-cutoff mixed current-frame certificate, with the
value-error term and independent chart gate retained.  Using `Sfloor` is
conservative relative to the physical `2^(2h(ell-2))`; it is not a floating
evaluation.  The other finite pulse length, carrier, normal, inverse and
background-domain thresholds are joined by taking their maximum with (7),
as in the existing value-domain and pulse foundations.  That maximum is
finite and independent of physical derivative order.

Equations (6)--(9) prove existence of a certified tail.  They do not certify
an absent old `ell0`, and they do not by themselves relocalize the old force.

## 4. What transfers to physical pulses

The finalized global coefficient bounds now supply the finite current-frame
jet inputs that the [physical-pulse foundation](../physical_pulse_foundation/README.md)
previously left conditional, while (2)--(9) supply its quantitative cone and
tail selection.  On the resulting common tail, the already derived primary
physical estimates therefore keep

\[
 s_x={1\over2}+{h\over2},\qquad
 s_t=1+{3h\over2}.                                      \tag{10}
\]

For every fixed spatial order `a` and time order `b`, uniformly over every
accepted scale, label and copy,

| Primary quantity | Bound |
|---|---|
| complete curl velocity | `C_ab Q^{-[(a+1)sx+b st]} S^B sqrt(zeta) delta^-B P` |
| wave pressure | `C_ab Q^{-[(a+2)sx+b st]} S^B sqrt(zeta) delta^-B P` |
| fixed vector potential | `C_ab Q^{-[a sx+b st]} S^B sqrt(zeta) delta^-B P` |

`K0,5`, `K1,3*` and `H2,3*` enlarge only fixed prefactors
and polynomial `S` powers.  They do not change (10), the fixed transverse
support `|xi|<=r0/2`, the pulse-cutoff support
`|v/Ls-1/2|<=2/5`, the exact 1250-color disjointness or smooth zero extension
of the cut potentials/pressures.

This closes the **primary pulse transfer on the new tail**, relative to its
named analytical foundations.  It does not prove that every source generated
by the nonlinear correction hierarchy satisfies the sourced-inverse
hypotheses.  Hence “all complete waves” or the complete residual must not be
inferred from this checkpoint.

## 5. Reproducibility and remaining obligation

Exact tests cover the algebraic slack in (3), missed-gate controls for an
unspecified admissible tilt, binary ceilings, sixth-root ceilings, monotonic
later-band gates, and all physical exponents through total order six.  A
manufactured rational fixture exercises (6)--(9); it is not the enormous
profile or a PDE sample.  Frozen dependency hashes and current source bytes
are recorded.

Next propagate this separately certified tail through localization and the
force, or supply evidence for a genuinely pre-existing first integer.  Then
audit the complete-wave covariance means—including axial velocity, curl
remainders and both exceptional pressure tails—and verify that every
generated source lies in the sourced inverse class before assembling
viscosity, all cutoffs and every nonlinear interaction in the complete
residual.  No unconditional residual improvement or blow-up follows.

Reproduce with:

```sh
python3 test_audit.py -v
python3 audit.py --verify-record evidence.json --output fresh.json
```

The output must be fresh; historical evidence is never overwritten.
