# Q186: correct the exact-core/finite-map frontier to the actual Q2 infima

Parent: Q185, `0653ebda03b82fdfb46dfda749bbde53d4babdd0`.
The live branch, current CI, Q164/Q165, the root foundations, Q1/Q2 and the
frozen provider evidence were read first. Historical evidence and every
physical parameter, pressure gauge, support, selected root, frequency and
`C^-20` activation width are unchanged.

Q185 correctly replayed the diagonal and comparison at the corrected source
rates, but its next obligation was too broad. It asked for independent
source-level validation of the exact core and all finite existence, inverse
and gap lemmas. At current head the exact core/leading-pressure calculation
has already been independently rederived by Q164, and the actual angular,
heat, axial, annular and compact matching maps have explicit source audits.
Q165 binds the finite Q2--Q5 chain to that object.

The first surviving independent-validation boundary is narrower. The Q2
same-band theorem is symbolic in the actual selected-profile regional infima.
Its executable freezes a manufactured rational fixture. That fixture checks
the algebra of the radius formula; it is not the actual ten-region profile.

## 1. Exact core and finite roots are not the current missing objects

Q164 rederives the all-index product and integrated-derivative bounds, all 21
normalized source terms, the complete three-term leading pressure source,
the invariant ball and the strict contraction below `10^-50`. It recovers
the original equations from the fixed point, removes the regular forward
homogeneous modes by axis regularity, and identifies the same-gauge completed
leading backward pressure. Later positive-order, wave, mean, localization and
forced pressure memory remains present.

The root audits separately trace the actual maps:

1. the reference angular and heat moment maps and terminal waiting root;
2. the axial two-moment map and complete scalar strength equation; and
3. the five-row annular and compact-repair maps, with inverse bounds below
   `610` and `800` respectively.

Q165 then binds Q2--Q5: the weighted regional gaps, full moving frame and
whole path, fixed compensator/projector profiles, and logarithmically
renormalized inverse/source domain. These conclusions remain relative to
their named continuum providers. They do, however, make “the exact core and
some unnamed finite map” an obsolete description of the next frontier.

## 2. What Q2 actually needs

For the selected leading profile define

\[
 f=\inf F_0,\quad t=\inf(-g_{0\theta}),\quad
 d=\inf D_0,\quad D_0=|g_0|^2+2F_0g_{0\theta},
\]

together with finite `Fbar=sup F0`, `Gbar=sup |g0|`, the orientation margin
`omega`, the same-tilt deficit parameter `gamma`, and the fixed positive tilt
`u`. Set

\[
 L_D=4\bar G+2\bar F+3,\qquad L_P=\bar F+\bar G+1,
\]

\[
 \bar D=\bar G^2+2\bar F\bar G,
 \quad L_c={2L_D\over ft}+{2\bar D L_P\over f^2t^2},
 \quad L_H=1+{2\over t},
\]

\[
 L_\phi=3u^2L_H+(1+u^2)(4L_c+3\bar c^2L_H),
 \quad m_\phi=u^2\omega^2\gamma(2-\gamma).
\]

The exact perturbation calculation gives

\[
 \rho_{\rm gap}=\min\left\{{1\over2},{f\over2},{t\over2},
 {d\over2L_D},{\omega\over2L_H},{m_\phi\over2L_\phi}\right\}>0. \tag{1}
\]

If the weighted background estimate is

\[
 \max\{|F-F_0|,|g-g_0|,
 |(\widehat T-T_0)/|T_0||\}\le C_{\rm val}q^{2h}, \tag{2}
\]

then the actual promotion condition is

\[
 \boxed{C_{\rm val}q^{2h}\le\rho_{\rm gap}
 (f,\bar F,t,\bar G,d,\omega,\gamma,u).}              \tag{3}
\]

At the two flat stress edges, (2) uses the exact factored quotient. It does
not divide by a positive absolute stress floor, because that floor is zero.
Complete `Pi1/Pi2`, axial velocity and later pressure memory stay in the
weighted construction.

For a dyadic chart `q<=4*2^-ell`, (3) follows from

\[
 2^{2h(\ell-2)}\ge {C_{\rm val}\over\rho_{\rm gap}}. \tag{4}
\]

Thus, once the actual constants in (1)--(2) are finite and positive, the
Archimedean property supplies a finite later `ell`. No numerical resolution
of the huge selected profile is needed for this existential conclusion.

## 3. The frozen rational radius is a diagnostic fixture

The Q2 executable calls (1) with

\[
 (f,\bar F,t,\bar G,d,\omega,\gamma,u)
 =(1/5,3,1/7,4,1/11,1/13,1/17,2),
\]

and obtains

\[
 \rho_{\rm fixture}={33\over387618781940}.           \tag{5}
\]

Equation (5) is exact and all six fixture rows are positive. None of its
eight inputs is identified with the actual selected-profile infimum in the
frozen executable. Substituting (5) for the right side of (3) would therefore
be unjustified.

The analytical regional provider does derive positive value margins over its
ten pieces and the background certificate derives (1)--(4) from such inputs.
But those checkpoints explicitly rely on replayed core, envelope, loop and
inverse lemmas and do not independently reprove every primitive provider.
Their conditional route remains valid. The missing independent-validation
work is to rederive the actual ten-region `f,t,d,omega,gamma` and weighted
`Cval` inputs on one selected profile and bind them to (3), not to rerun (5).

## Result and next obligation

This checkpoint corrects Q185's frontier without weakening its conditional
diagonal, localization, comparison or protected-core replay. The exact core
and named finite root maps already have source-level audits. The first
surviving independent-validation edge is the actual Q2 input set and the
precise inequality (3).

Q186 does not independently rederive those ten-region infima, promote the
fixture to an actual-profile certificate, construct a future solution, or
prove an unconditional pulse bound, residual improvement, admissible force,
competitor exclusion or Navier--Stokes blow-up.

Next independently rederive the actual selected ten-region
`f,t,d,omega,gamma` and weighted `Cval` constants and verify (3) with one
later-tail selector. Then continue through any surviving provider premise in
dependency order, retaining the corrected Q185 rate route and complete
pressure memory.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q186_exact_finite_frontier_correction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q186_exact_finite_frontier_correction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q186_exact_finite_frontier_correction/evidence.json \
  --output q186_frontier.json
```
