# Q177: reconcile selected inputs, weighted stress, and reserved supports

Parent: `d125a9aeec0eec74cea0ebbcc7d63da36825be74` (Q176).
Status: `Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled`.

Q176 validated the same finite-frequency modulation and exact five-row repair,
and corrected the historical order of the low-order acceptance and later
matching-root reaudit. Its remaining frontier was stated too broadly if read
as saying that the repository lacks selected-profile continuation and regional
input estimates. The frozen analytical providers already derive those inputs
conditionally. This checkpoint records that correction, audits their join to
the global pressure/moment identities, and checks the next weighted-stress and
coefficient/reserved-support providers at source level.

Historical evidence, fields, roots, frequency, pressure gauge, supports,
cutoffs and the `C^-20` width are unchanged. The provider chain still depends
on the retained exact pressure-family core and finite existence, inverse and
gap lemmas. Replaying its records is not an external certificate or a new
proof of every continuum theorem.

## 1. The selected input and frequency chain is already analytical

The actual source order supplies the following implications.

1. `prescribed_continuation_foundation` constructs the pre-modulation fields,
   five moments and full pressure from the exact core, retaining all fixed
   cutoff costs and one common parameter strip.
2. `compact_jet_envelope` proves that the selected unmodulated compact inputs,
   denominators and required mixed jets fit the envelope `A=C^4096`.
3. `loop_modulation_bounds` derives the fixed-phase modulation constant
   `Cstate=D=H^16`, including the moving inverse phase and the
   inverse-`lambda` target loss.
4. `repair_state_bounds` proves the positive repair neighborhood,
   `Ccorr=A^128`, radius `A^-128`, and accepts
   `N=1+floor(H^32)` relative to those inputs.
5. `regional_cone_join` derives the ten regional value bounds and joins the
   same finite `N` to the protected inner and compensated outer collars.

For the state, contraction and coefficient-radius thresholds the repair
source bounds their logarithms by

\[
 16G^4+193G<32G^4,
\]

which reduces to `193G<16G^4` and already holds for `G>=3`. The cone loss is
absorbed because

\[
 512\,4^{10}<G^{30},
\]

which holds for `G>=4`. The selected `G=exp(A^32)` is far larger. These are
finite symbolic inequalities; neither `N` nor the giant profile is rounded or
numerically materialized.

Thus Q176's statement that the selected continuation and regional suprema
remain unrecomputed is correct only as an independent-validation warning.
Within the repository, the named providers already form a conditional
analytical route. Their exact-core and finite-map premises remain explicit,
and none of their manufactured fixtures is substituted for the continuum
fields.

## 2. Exact global pressure and five moments remain on the same roots

The annular/exterior provider composes the whole compact input with the same
accepted maps and preserves the nonzero reserved leading memory
`M0=mu eta`, `V0=-mu`. The global pressure provider then proves that both
transformed row systems have the original zero set `(M,R,J,S,Cp)`, retains
`lambda>0`, and gives

\[
 M_\infty=R=J_\infty=S_\infty=0,
 \qquad \Pi_{axis}+C_{p,\infty}=0,
\]

so the completed leading pressure is the same-gauge backward tail. The
functional `eta` dependence and one-sided jets are retained. A four-row
repair would leave the axial pressure-stress defect and is not used.

## 3. Weighted first stress is smooth at the outer flat edge

The weighted-stress source factors the actual leading and first stresses as

\[
 T_0=e^{-4/\delta^2}\delta^{-3}(b_0,\delta^6b_z),
 \qquad T_1=e^{-4/\delta^2}\delta^{-3}(b_1,0),
\]

with smooth coefficients and `b0>0`. Therefore

\[
 W_1={T_1\over|T_0|}
 =\left({b_1\over\sqrt{b_0^2+\delta^{12}b_z^2}},0\right)
\]

is smooth without dividing by `eta` or `b1`. At `eta=0` the exact expansion
is

\[
 W_{1\theta}(\delta,0)={X_b\over8}\delta^3+O(\delta^4),
 \quad \partial_\delta^3W_{1\theta}(0,0)={3X_b\over4},
 \quad \partial_\eta^2W_{1\theta}(0,0)=4X_b.
\]

The nonzero cubic term retains the radial Jacobian and the `beta f_y` term.
The same provider joins every fixed weighted leading/first-stress edge and
chooses one admissible exact-real diagonal after the finite seminorms are
known. It does not assume a convergent heat Taylor series at the poles or
constants uniform in derivative order, `h`, `N`, or coefficient index.

## 4. Positive-order supports do not reuse the reserved mean patch

The coefficient source solves the actual current-order linear system, not a
new quadratic leading repair. Its two blocks obey

\[
 \|B_U^{-1}\|<2,
 \qquad \|B_E^{-1}\|<787.562,
 \qquad \|B^{-1}\|<800.
\]

Every fixed angular and radial jet is finite with the retained amplitude,
`lambda^-1`, radius and width costs. The support geometry is exact:

\[
 2.55=51/20<3<6<11.
\]

All I3 repair bumps end before logarithmic radius `3`; positive-order
velocity and its moments vanish there; the reserved I4 patch begins at `6`.
Consequently the later coefficient construction does not silently reuse the
leading repair support. The pressure is reconstructed rather than cut off and
keeps

\[
 \Pi_1=K_1(\eta)-{\mu^2\over4X},
 \qquad \Pi_2=K_2(\eta),
 \qquad \Pi_n=0\quad(n\ge3).
\]

`K1` and `K2` contain later source memory and are not set to zero. The leading
radial velocity also remains nonzero on I4. These terms must be carried into
all complete covariance and mean corrections.

## 5. Scope and next obligation

This checkpoint establishes the in-repository conditional route through the
selected compact inputs, finite frequency, cone, exact pressure/moments,
weighted first stress, and finite coefficient/reserved-velocity supports. It
also identifies the narrower unresolved premise: the exact core and finite
existence/inverse/gap providers have not been externally certified or
reproved here from first principles.

The next source-level obligation is `physical_pulse_foundation` and the Q149
same-profile join. That audit must retain the distinction between the uniform
primary pulse theorem and a sourced inverse whose input class is separately
assumed. Afterward, complete covariance/mean corrections must retain radial
and axial velocity, `Pi1/Pi2`, moving pressure and every cutoff derivative.
No unconditional pulse hierarchy, complete residual improvement, force
admissibility, exclusion result or blow-up follows here.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q177_weighted_support_source_reconciliation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q177_weighted_support_source_reconciliation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q177_weighted_support_source_reconciliation/evidence.json \
  --output fresh.json
```
