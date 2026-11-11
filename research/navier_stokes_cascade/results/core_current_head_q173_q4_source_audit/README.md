# Q173: Q4 fixed-profile source audit and generated-source boundary

Parent: `e8f403601c85d2942f4cd5bb2d9c127aa5a35f7e` (Q172).
Status: `Q173-Q4-fixed-profile-source-audited-with-generated-source-boundary`.

Q172 corrected the uniform-frame argument on the accepted infinite tail. This
checkpoint audits the next fixed compensator/projector contract against that
corrected tail and the underlying profile, geometry, recomputation and support
providers. No selected field, width, gauge, pressure tail, profile or frozen
record is changed.

## 1. The fixed-profile geometry survives the Q3 correction

The compensator jets do not rely on Q3's invalid compactness sentence. The
reserved-geometry provider proves a band-independent complex neighborhood for
the implicit scale coordinate `s`, and hence explicit fixed-order `Z,T` jets.
For one fixed `rhohat in C_c^infinity((ell,u))`, ordinary finite chain rules
then give finite profile-dependent constants at every fixed order. There is no
claim uniform over profile choices, derivative orders or infinitely many
construction cycles.

Set

\[
 A_2=\int\widehat\rho^2,
 \qquad \sigma_e=R^{-e}\sqrt{s}\,\rho^2/A_2,
 \qquad e=1,2.                                         \tag{1}
\]

Unit mass and Cauchy--Schwarz give `A2>0` and
`A2^-1 <= u-ell`. With `R=sqrt(s)x`, the exact change of variable yields

\[
                 \int R^e\sigma_e\,dR=1.               \tag{2}
\]

If `K=supp(rhohat) compactly contained in (ell,u)`, then the actual support is
`sqrt(s)K`, strictly inside the moving interval. The interval and the common
hull are enclosures, not the support itself. The logarithmic weight floor is
used only on this actual support.

The complete moving axial pressure moment remains

\[
 \epsilon\partial_Z(J_z+c_\rho P)
 =\epsilon(J_{z,Z}+c_{\rho,Z}P+c_\rho P_Z),             \tag{3}
\]

with current full pressure and both `Pi1/Pi2` tails. The middle term is
generally nonzero.

## 2. Correct scope of the same-choice recomputation

The explicit projectors in (1) are a new admissible choice. The source chain
correctly rebuilds `h_e`, both edge primitives, signed amplitude/root input,
curl, all six covariance rows and pressure norm recipes after making that
choice. It does not reuse coefficients from an unspecified projector.

The source record is explicitly a **conditional norm-chain recomputation**.
It retains the actual initial Haar, target, root/value/cone and support
enclosures as premises. Q4 therefore closes the fixed-profile geometry,
moment, support and identity contract; its wording must not be promoted to an
actual continuum enclosure for every recomputed target or generated source.

This distinction matters even when every inverse operator is uniformly
bounded. The identity operators on scalar sources have norm one, while the
source family `F_ell=ell` and its outputs are unbounded. Operator boundedness
does not manufacture a source-family estimate.

## 3. First blocking inequality

For every requested finite stage and jet order, the next required same-state
input remains

\[
 \boxed{
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 \max_{|\alpha|\le R_j(m,N)}\sup_{\Omega_{j,\ell}}
 \zeta^{-1}\delta^{p_{j,\alpha}}
 |D^\alpha F_{j,\ell}|<\infty .}                       \tag{4}
\]

Fixed projector jets, exact moments, support closure and bounded inverse maps
do not prove (4). It requires the actual recomputed current-target/source
bounds with complete pressure, axial velocity, covariance and stage history.
The executable finite fixture verifies only that a supplied finite weighted
source constant can pass through the algebra; it does not supply (4).

## 4. Result, limits and next obligation

The Q4 fixed compensator/projector contract is now source-audited on the
Q172-corrected tail, with the conditional recomputation boundary made exact.
Historical evidence remains frozen. Tests and hashes corroborate the recorded
algebra and provenance; they are not a continuum proof of (4).

Next audit Q5's shifted compact maps, zero-Haar temporal inverse and fixed
slow response on these same-choice coefficients and complete weighted
domains. Keep (4) as an explicit downstream premise. Complete covariance and
fresh axial/full-pressure means, generated interactions, viscosity, cutoffs,
the residual, force, exclusion and blow-up remain unproved.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q173_q4_source_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q173_q4_source_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q173_q4_source_audit/evidence.json \
  --output fresh.json
```
