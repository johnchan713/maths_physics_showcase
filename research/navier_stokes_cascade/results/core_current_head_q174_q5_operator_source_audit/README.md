# Q174: Q5 operator audit and typed same-family source boundary

Parent: `afa4cd8d2ab1e23f434dbd25d0ce0de55866d445` (Q173).
Status: `Q174-Q5-operator-source-audited-with-family-and-typed-domain-boundary`.

Q173 established the fixed compensator/projector geometry on the corrected
Q2/Q3 tail and left generated-source membership open. This checkpoint audits
the shifted compact maps, zero-Haar temporal inverse and fixed slow response
against the later radial-weight and physical-source corrections. No field,
profile, root, support, gauge, band or frozen historical record is changed.

## 1. The operator estimates survive on their declared inputs

In `U=R^d`, retain the physical shift and the two necessary orientations

\[
 I_-(U,y)=\int_{-\infty}^0F(U+a,y+Ma v_r)\,da,
 \quad I_+(U,y)=\int_0^\infty F(U+a,y+Ma v_r)\,da,
\]
\[
 T_eF=g_e\{(1-\chi)I_- -\chi I_+\},\qquad
 A_eF=k_e\chi_U\int_{\mathbb R}F(U+a,y+Ma v_r)\,da. \tag{1}
\]

The inner and outer Gaussian factors absorb every fixed pole; the cutoff
derivative lies on an interior weight floor. On a nonzero torus mode,
integration by parts gives the recorded `M^{-p}` gain with the explicit
finite source request. The zero-Haar temporal inverse costs four torus
derivatives, and the fixed five-profile slow response costs only its fixed
finite profile and `lambda^{-1}` constants. These are analytical one-way
mapping estimates. They require the matching source norm and exact weighted
moments; they cannot prove those inputs.

The actual axial row remains `g_d-A_1g_d`. The slow response retains

\[
 \epsilon\partial_Z(J_z+c_\rho P)
 =\epsilon(J_{z,Z}+c_{\rho,Z}P+c_\rho P_Z),             \tag{2}
\]

old pressure, both `Pi1/Pi2` tails and all six covariance entries. The
complete potential support may fill the reserved convex hull.

## 2. Correction: the source domain is typed

Q5 stated a full-`zeta` bound for every generated source. Q18 gives an exact
obstruction. With

\[
 \zeta(x)=\exp(-w^2/x^2-4/(L-x)^2),
\]

the admissible half-weight function `sqrt(zeta)` satisfies, for every finite
pole `p`,

\[
 \zeta^{-1}x^p\sqrt\zeta\longrightarrow\infty
 \quad(x\downarrow0).                                  \tag{3}
\]

No finite logarithmic factor, derivative request or epsilon reserve repairs
this radial-weight mismatch. Q19 therefore routes nonzero harmonics through
the half-weight wave space and means/signed targets through the full-weight
space. Q20 requires jets of the Gaussian quotient on the whole lifted path;
Q21 transfers a proved extended-domain estimate to that path. Q26 proves
fixed-`q` domain closure, not source magnitude.

The later canonical chain does not close this gap. Q94 derives a conditional
first-real source polynomial while recording the native provider inequalities
as uncertified. Q117 propagates a supplied current source through the
covariant response while retaining the actual native entrance and high-pole
contracts as open. Those same-object records reinforce the typed source
frontier rather than supplying (4)--(5).

Thus Q173's pure-power display is a sufficient but non-exact gate. For every
fixed stage, jet and target, the exact canonical inputs are finite log degrees
and the two typed bounds

\[
 \boxed{\sup_{\ell\ge\ell_*}\epsilon_\ell^{-B_j}
 (1+|\log\epsilon_\ell|)^{-K_{j,R}}
 \max_{|I|\le R}\sup_{\widetilde\Omega_{j,\ell}}
 \zeta^{-1/2}\delta^{p_{j,I}}
 |D^I(f_{j,\ell}/P_{j,\ell})|<\infty,}                \tag{4}
\]
\[
 \boxed{\sup_{\ell\ge\ell_*}\epsilon_\ell^{-C_j}
 (1+|\log\epsilon_\ell|)^{-L_{j,R}}
 \max_{|I|\le R}\sup_{\Omega_{j,\ell}}
 \zeta^{-1}\delta^{q_{j,I}}|D^I\Xi_{j,\ell}|<\infty.} \tag{5}
\]

Equation (5) keeps both exact tangential moment identities. Equations
(4)--(5) retain the full endpoint extension, all six covariances, axial
velocity, moving pressure product and pressure memory.

## 3. The same-family obstruction

Q160--Q161 derive analogues of (4)--(5) by simultaneous induction for the
**inherited physical primary family**. That is a meaningful continuum theorem
on its stated object. It cannot be substituted into the current Q3--Q5
canonical chain: the canonical full-path provider explicitly records
`canonical_and_inherited_primary_families_identified = false`, and Q149
preserves that the inherited actual-mixing entrance is not the canonical
reference-mixing entrance.

This is an object-identity obstruction, not a claim that either family or
source estimate fails. Hash overlap, common background data, common
projectors and identical rate arithmetic do not identify the two pulse
solutions or their generated sources. Consequently the Q5 operator contract
is now audited on Q173's coefficients, while actual canonical source
membership in (4)--(5) remains open.

## 4. Result and next obligation

The shifted compact, temporal and slow maps preserve their typed weighted
domains at every fixed order, conditional on matching same-family inputs.
The blanket full-edge source statement is corrected, and physical-family
source bounds are not promoted to canonical-family bounds. Tests and frozen
hashes verify formulas and provenance; they are not continuum proofs of
(4)--(5).

Next either prove (4)--(5) for the canonical reference-mixing family, or
rebuild the current Q3--Q5 spine explicitly on the inherited physical family
with its own frame, root and same-choice projector identifications. Only then
replay Q6's complete residual with viscosity, cutoffs, pressure and all
interactions. No unconditional pulse, force, exclusion or blow-up conclusion
follows.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q174_q5_operator_source_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q174_q5_operator_source_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q174_q5_operator_source_audit/evidence.json \
  --output fresh.json
```
