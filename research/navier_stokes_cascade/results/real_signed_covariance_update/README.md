# First phase-locked REAL signed covariance update

Parent: `e06495cda01638fe17e20f21f5032c606c806099`, the actual-primary
initial mean checkpoint. On 2026-10-16 its [CI run](https://github.com/johnchan713/maths_physics_showcase/actions/runs/38906537757)
was still running, with the initial means, physical pulses, compact means
and conditional cycle steps already successful. This checkpoint preserves
every historical source and evidence record and uses the same fixed Md64,
h, N, positive lambda, width `C^-20`, support geometry and primary waves.

Result using the retained analytical lemmas: the first signed update about
the ACTUAL unchanged real primary wave is well-defined with all fixed
weighted coefficient/physical jets and closed supports. For its target
order `alpha=2-3kappa`, the complete pair-covariance error and the updated
auxiliary-averaged tangential residual have order `3-6kappa`.
The signed increment square and every curl term are retained. This sharper
bound uses phase-locked real homogeneous amplitudes; it is not a generic
later complex-amplitude estimate or a complete residual improvement.

## 1. Join actual inputs to the fixed signed inverse

Use the all-later-band domain, discrete choices, fixed `Q=2^-ell`,
`epsilon=Q^h`, `S=ell^2`, `kappa=1/100000` and classes from the
[physical transfer](../physical_pulse_transfer/README.md). The
[parent](../complete_covariance_initial_means/README.md) has constructed
the actual mean state without changing the wave: complete primary curl
`U in W_(1/2)`, tangential means `v,gamma in M_(1-kappa)`, radial mean
`beta in M_(2-kappa)`, current pressure reconstructed from the full radial
source, and exact zero velocity moments. Its three compatibility defects
are in `S_2`, while both full tangential means are `M_alpha` modulo the
retained flat summands, for `alpha=2-3kappa`.

Here a bar denotes AUXILIARY Haar average; angular average is separate.
Exclude the independently retained flat BASE residual from the mean
increment equations. All actual compact reconstruction errors remain in
those equations and velocities. Their interior support and arbitrary
fixed-order flatness allow inclusion in the stated M classes without
pretending they vanish. The base flat error need not have shell support
and is not assigned a shell edge weight.

Define the real, auxiliary-independent pair target using the CURRENT
means and full pressure:

\[
\Xi_\theta=-T_2[\bar E_\theta-\sigma_\theta\epsilon\partial_ZJ_\theta],
\qquad
\Xi_z=-T_1[\bar E_z-\sigma_z\epsilon\partial_Z(J_z+c_\rho P)].
\tag{1}
\]

The exact integrated identities give zero weighted integrals of the two
sources. Thus the compact inverses are exact here. The moving product
`c_rho P` is differentiated as a whole. Current moments are `S_2`, so
their bump sources have order3; consequently `Xi in M_alpha`, with
smooth shell/slow extensions. Haar averaging and radial inversion may
fill the auxiliary torus, but (1) itself is auxiliary-independent.

On each actual label use the FIXED positive primary roots
`a_sigma=sqrt(y_sigma)` and actual pulse-column matrix H. Their audited
estimates, retaining the fixed positive radius costs, are

\[
|D^I H^{-1}|\le C_I S^{B_I},\quad
y_\sigma\ge cS^{1/2}\zeta,\quad
|D^I y_\sigma|\le C_I S^{B_I}\zeta\delta^{-B_I}.
\]

The first bound is shorthand for the more precise radius-dependent
`r0^(-3/2)` costs at the SAME fixed admissible radius, not a uniform
claim as the radius tends to zero. Set

\[
d=H^{-1}(\Xi/\epsilon),\quad
\delta a_\sigma={d_\sigma\over2a_\sigma},\quad
L\Xi=\sqrt\epsilon\sum_{\sigma=\pm}\delta a_\sigma b_\sigma.
\tag{2}
\]

The target estimates supply `D^I d=O(epsilon^(alpha-1) S^B zeta delta^-B)`.
Every derivative of `y^(-1/2)` is a sum of
`y^(-n-1/2)` times n differentiated y factors. After multiplication by a
derivative of d, its edge exponent is exactly
`(-n-1/2)+n+1=1/2` for EVERY n. Hence (2) is in `W_(alpha-1/2)` with
all fixed jets and smooth zero extension, using the exact inherited
exponential weight, not a modified edge constant. The actual-background,
cone, pulse and weighted-edge hypotheses of this fixed inverse have now
been supplied by the earlier checkpoints. Its continuum inverse/pulse
lemmas remain explicit analytical dependencies, not test-derived facts.

There is no sign or small-target condition. Never take a square root of
`y+d`, divide by a slow partition function, or replace the denominator by
an updated amplitude. A generic weighted target can have an unbounded
relative amplitude near an edge even while its signed wave extends
smoothly by zero; positivity of the UPDATED amplitudes is unnecessary.

The primary real-cosine factor 1/2 is already in H. Thus
`B(W0,L Xi)=epsilon H(2a*delta a)=Xi`. Physically assemble the same slow
partition on BOTH factors. Each chart target is `Q^(2A) xi_phys`, so
`sum Q^(-2A) eta_beta^2 Q^(2A) xi_phys=xi_phys` by `sum eta_beta^2=1`.
Closed separated supports remove distinct-label products. Bounded chart
overlap and polynomial S counts preserve the estimates across ALL
accepted later scales. This is not a sampled-band or numerical ODE claim.

## 2. Take the actual curl: a sharper INITIAL covariance identity

Cut the potential and associated wave pressure BEFORE the curl. Since
the homogeneous pulse, target and all cutoffs are real and angularly
independent in cylindrical coefficients, each label has exactly

\[
U_\beta=A_\beta\cos\vartheta_\beta+R_\beta\sin\vartheta_\beta,
\quad V_\beta=L_\beta\cos\vartheta_\beta+S_\beta\sin\vartheta_\beta,
\tag{3}
\]

with the SAME phase `vartheta_beta=k_beta Phi_beta`. All slow and cutoff
curl commutators are in R,S, not discarded. The classes are

\[
A\in W_{1/2},\quad R\in W_{1-\kappa},\quad
L,V\in W_{\alpha-1/2},\quad S\in W_{\alpha-\kappa}.
\tag{4}
\]

At each auxiliary point, the mixed cosine/sine angular average vanishes.
For the full six-entry tensor, define
`Lambda=1/2 sum bar(A tensor L+L tensor A)` and let `Delta W` be the
change in the doubly averaged complete wave tensor. Exactly,

\[
\Delta\bar W=\Lambda+\mathcal Q,\qquad
\mathcal Q=\tfrac12\sum\overline{
R\otimes S+S\otimes R+L\otimes L+S\otimes S}.
\tag{5}
\]

The two cross entries of Lambda are Xi; its OTHER entries need not be
zero. Equation (5) bounds every updated component, not just radial flux.
The error orders are respectively

| Retained contribution | Coefficient order |
|---|---:|
| old sine curl times new sine curl | alpha+1-2kappa |
| signed transverse square | 2alpha-1 |
| new sine curl square | 2alpha-2kappa |

Thus `nu=min(alpha+1-2kappa,2alpha-1)` for `alpha>=1` and
`Delta bar W-Lambda in M_nu`. In particular for `alpha=2-3kappa`,
`nu=3-6kappa`: the signed square is the limiting term, NOT zero.
The covariance gain is positive precisely when
`min(1-2kappa,alpha-1)>0`; at alpha=1 the square prevents a gain.

No additional graph loss is charged AFTER Haar averaging:
`bar D_r f=partial_R bar f`. Coefficient derivatives of (5) preserve
its orders; the graph losses already incurred inside R,S are retained.
The full tensor CHANGE is only `M_alpha` for alpha>=1, because Lambda
has that order. This fact matters in the radial and axial pressure terms.

The generic Corollary7.8 estimate is still valid and unchanged. It allows
arbitrary old complex harmonic coefficients and is weaker. Once a sourced
wave has introduced an old cosine-quadrature remainder E, a term
`B(E,L)` survives; (5)'s sharper bound is not licensed. All physical
velocities are real: "real" here means the phase-locked coefficient
structure (3), not merely taking a real part of general complex amplitudes.
PDF pp.84 and87, (7.31)--(7.42), were rendered and visually checked for
the factor2, squared partition and full actual-velocity remainder.
Source SHA256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.

## 3. Recompute the FULL pressure and all three defects

Keep the actual mean VELOCITY fixed in this substep, replace the wave by
U+V and recompute its full covariance, radial source and compact pressure.
Both background pressure tails remain in the same `p_B`; they cancel
only when subtracting that SAME exact base balance. The new nonzero wave
pressure has zero angular mean, not zero pointwise value. Neither fact
licenses dropping the new mean pressure change.

The exact changes are

\[
\Delta\bar g_r=-(\partial_R+1/R)\Delta\bar W_{rr}
 -\epsilon\partial_Z\Delta\bar W_{rz}+\Delta\bar W_{\theta\theta}/R,
\quad \Delta P=\int\Delta\bar g_r,
\]
\[
\partial_R\Delta\bar p_m=\Delta\bar g_r-\rho\Delta P,\quad
\int R\Delta\bar p_m=-\tfrac12\int R^2\Delta\bar g_r+c_\rho\Delta P,
\]
\[
\Delta J_\theta=\int R^2\Delta\bar W_{\theta z},\quad
\Delta J_z=\int R\Delta\bar W_{zz}-\tfrac12\int R^2\Delta\bar g_r.
\tag{6}
\]

All changes in (6) have order alpha; pointwise before Haar averaging
`Delta g_r,Delta p_m` have the conservative order `alpha-kappa`.
The complete compact pressure identity also retains its actual flat A0
remainder, with zero Haar average. Thus `epsilon partial_Z Delta bar p_m`
has order alpha+1. Keep the WHOLE `partial_Z(c_rho Delta P)` in the
integrated axial residual. In particular `W_zz`, `W_theta z`, and the
pressure moment are not inferred to vanish from leading polarization.

Using `(partial_R+e/R)Xi` to cancel the two old mean sources in (1),
the updated auxiliary-averaged tangential rows consist of

| Contribution | Order |
|---|---:|
| radial divergence of the covariance error (5) | nu |
| original retained compatibility bumps | 3 |
| new axial covariance flux | alpha+1 |
| changed averaged axial pressure | alpha+1 |

Consequently BOTH averaged tangential rows are `M_(3-6kappa)`.
But the full, unaveraged rows are conservatively only `M_(alpha-kappa)`;
the recomputed defects are `S_alpha`, NOT asserted to stay in `S_2`.
The velocity constraints persist because this step changes only angularly
nonzero waves. Pressure and every defect must be recomputed from U+V;
reusing the old ones would invalidate the next mean step.

## 4. This substep's generated wave sources and physical jets

Write `b=alpha-1/2` for the signed wave order (not the radial background).
Its homogeneous principal equation has no new fast source because
delta a is auxiliary-independent. The inherited full linear remainder
is `W_(b+1/2-3kappa)=W_(alpha-3kappa)`, including viscosity, normal
pressure, slow derivatives, cylindrical terms and potential cutoffs.
The homogeneous cutoff tail is retained separately flat, not deleted.

For a COMPLETE harmonic curl amplitude a, exact divergence gives
`ikm n.a=-div_*a`. It removes the apparent half-frequency loss in a
same-label wave--wave product, giving order `b+b'-kappa`. Distinct
labels have zero derivative products by CLOSED support separation; this
cancellation cannot be used for overlapping distinct phases. Mean--wave
transport DOES retain its half-frequency loss. The parent's actual mean
order is `mu=1-kappa`, with an extra radial power, so this update's
generated nonzero rows have orders

| Generated term | Order |
|---|---:|
| full linear remainder | alpha-3kappa |
| cross with complete old primary | alpha-kappa |
| cross with actual mean, including flat velocities | alpha-kappa |
| signed complete wave self-interaction | 2alpha-1-kappa |

Each nonzero term retains a wave factor, its closed slow/transverse/time
supports, its `sqrt(zeta)` edge weight and its Gaussian P envelope:
`P^2<=P`, and the stronger product edge weight is admissible. The original
potentials and pressures vanish on endpoint collars, so their derivatives
and products extend smoothly by zero. Finite harmonics remain finite
(here at most 2 times a primary angular multiplier). This verifies the
generated-source inputs of THIS substep, not every later sourced solve.
It does not cancel the OLD nonzero-angular residual left by initialization.

For space/time jet orders a,b, use the same losses `sx=1/2+h/2`,
`st=1+3h/2`. The uniform signed physical powers are

\[
u:\ -A+h(\alpha-1/2)-as_x-bs_t,\quad
\pi:\ -2A+h\alpha-as_x-bs_t,\quad
\mathcal A:\ 1/2-A+h\alpha-as_x-bs_t.
\tag{7}
\]

All labels, covers, Q and epsilon are held fixed during differentiation.
These are estimates with the inherited weights/envelopes and polynomial
S factors, not bare power equalities. The phase and support geometry
are unchanged, so the parent's all-fixed-order physical transfer proof
applies on every accepted later scale, with no new amplitude sign test.
Even the physical covariance error exponent `-2A+h nu` can be NEGATIVE;
normalized order gain is not absolute physical residual decay.

## 5. Evidence and remaining obligations

The new exact rational diagnostics verify the signed inverse/factor2,
squared partition scaling, half-edge-weight bookkeeping, all six updated
tensor entries, full curl divergence, variable-amplitude curl commutators,
the nonzero signed square, pressure reconstruction, both integrated flux
identities and the moving pressure bump. Controls detect old cosine errors,
missing axial pressure/covariance and false updated-amplitude positivity.
These are manufactured finite-regularity polynomial fixtures, not an
evaluation of the actual giant profile or an interval/proof-assistant
certificate. The all-orders derivations above use explicitly retained
analytical lemmas. The parent's current evidence and 610 inherited source
digests are replayed, not reinterpreted as a continuum proof.

Next: apply the fast inverse to the CURRENT zero-Haar residual after this
wave change, then the slow five-map to CURRENT recomputed defects; keep
all old/new products and pressure changes. After that, verify the later
sourced-wave and signed-covariance stages WITHOUT the special initial
quadrature shortcut, and join every finite cycle with viscosity/cutoffs.
The full nonlinear hierarchy, complete residual improvement, infinite
realization, forcing, energy and blow-up remain unverified.

```sh
python3 research/navier_stokes_cascade/results/real_signed_covariance_update/test_audit.py -v
python3 research/navier_stokes_cascade/results/real_signed_covariance_update/audit.py \
  --verify-record research/navier_stokes_cascade/results/real_signed_covariance_update/evidence.json \
  --output /tmp/real_signed_covariance_update_replay.json
```
