# Plan toward a verified finite-time blow-up

Updated: 2026-10-15, following research commit
`d80f6ff1b45b6111a400f04db463a9e3afd84ccf` and the full-potential-realization checkpoint.

The objective is a mathematical construction for the original three-dimensional
incompressible Navier–Stokes equations at positive viscosity. Our immediate
route is the supplied manuscript's forced construction. This is a research
strategy, not a prediction that the construction will succeed.

## Target and present position

We seek smooth divergence-free initial data u0, a force f with the required
smoothness and decay, and a solution u,p smooth for 0<=t<T, such that

\[
\partial_tu+(u\cdot\nabla)u-\nu\Delta u+\nabla p=f,
\quad\nabla\cdot u=0,\quad\nu>0,
\]
\[
\sup_{t<T}\|u(t)\|_2^2<\infty,
\qquad\limsup_{t\uparrow T}\|u(t)\|_\infty=\infty.
\]

The force must extend smoothly through T. It cannot supply its own singularity.
For the whole-space Clay alternative C, its derivatives and those of u0 must
meet the specified decay conditions. Finally we must justify why this local
singular solution rules out a global smooth bounded-energy solution with the
same data, rather than simply assert that implication.

The [official problem statement](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf)
allows forcing in the breakdown alternatives. Success on this route would
not automatically prove blow-up with f=0. The
[supplied manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
is the proposed construction to investigate, not a substitute for verifying
its arguments.

The repository contains analytic arguments and reproducible checks for a
joined leading profile, its derivative envelope, and oscillation errors.
Those results depend on preceding estimates. Manufactured diagnostic inputs
are not the actual global field, and the work has not received independent
mathematical review. No complete stress construction, PDE solution, or blow-up
has been verified. There is no defensible percentage of the proof completed.

## Ordered milestones and acceptance conditions

| Milestone | Work | Required result before advancing |
|---|---|---|
| 0. Audit the foundation | Recheck the actual core, joins, inherited parameter inequalities and loop estimates directly from their defining equations. Trace every quantitative assumption to a bound. | A dependency ledger with hypotheses, domain, derivative order, constants, source and proof for each required bound. No sampled fixture or success flag supplies a continuum estimate. |
| 1. Complete the five-moment repair | Propagate the correcting bumps through swirl positivity, shear, moments and integrated pressure. | An actual correction-to-state bound Ccorr and coefficient radius r, plus exact restoration of all five moment functions with the original axis pressure datum. |
| 2. Accept a finite oscillation frequency | Combine verified input, modulation and repair constants in the existing frequency criterion. | One finite N satisfying every required inequality. The corrected leading profile has all four admissibility gaps positive, preserves the protected axis region and matches the exterior. |
| 3. Build the full residual calculation | Compute all terms of R=dt u+(u.grad)u-nu Delta u+grad p in physical coordinates for the constructed field. | A derivative and scale ledger retaining viscosity, pressure, cutoffs, mixed derivatives, oscillatory cross terms and localization errors. Every uncancelled term has a proved bound. |
| 4. Construct and sum PDE corrections | Verify the wave and mean corrections, solvability conditions, derivative losses, residual improvement and convergence. | Convergence to a divergence-free field, smooth on every compact set before T, with the required energy control and residual regularity. Finite truncations alone do not pass. |
| 5. Establish admissible forcing | Localize and initialize the field while controlling the resulting residual, including behavior through T. | f is defined globally, meets the official smoothness/decay requirements for every derivative order, and comes from the exact equation. Continuity or finitely many bounded derivatives is insufficient. |
| 6. Prove actual divergence | Bound the complete corrected field below along points and times approaching T. Check data, energy and the global-solution contradiction. | A lower bound tending to infinity that survives all corrections, bounded energy, and a justified comparison with any putative global smooth solution. |

Milestone 0 runs alongside milestones 1–2. A later milestone cannot treat an
unverified earlier hypothesis as established. A failed argument sends the
affected construction back for revision; it does not prove the manuscript
impossible unless the obstruction itself is proved.

## Immediate work package: physical waves and the retained background stress

The [physical residual calculation](results/physical_residual_budget/README.md)
gives exact formulas for the unlocalized base and the first positive order.
The subsequent [higher-jet budget](results/higher_jet_budget/README.md) now
instantiates its derivative constant on the selected compact annulus:
`K=H^64`, conditional on the inherited analytic leading profile. It derives
the original C4 transfer and third-order loop, inverse-phase and final repair
bounds at the same N. This closes the compact background derivative gap;
it does not pass the full physical residual milestone, because wave tensor,
global coefficient extension and localization terms are still missing.

The executed [first-correction checkpoint](results/first_correction_construction/README.md)
now gives a conditional actual inner analytic coefficient on `X<=3/Lambda`,
strictly inside the protected leading core. Its six-variable system has an
explicit actual-core envelope `C_op=C^16`; a sparse Picard argument proves
convergence without requiring a tiny product of radius and operator bound.
The original incoming moment transfer extends to factorial C3, and the first
positive-order five-moment map has a linear inverse. A protected-core lower
bound survives this one correction at sufficiently small q. These are
conditional continuum arguments and separate manufactured algebra checks;
they are not a completed global first correction or blow-up proof.

The finalized compact C3 estimate includes the original C4 input that the
earlier checkpoint lacked. It also retains the positive N loss in second
radial derivatives. Its physical background residual upper bounds are still
singular as q tends to zero; bounded profile constants do not prove residual
decay.

The [annular first-correction checkpoint](results/annular_first_correction/README.md)
now continues the same coefficient beyond Xa into the preserved analytic-eta
collar, specifies the five actual incoming discrepancies, solves them linearly,
and closes the signed first tangential stress tails, conditionally on the
inherited leading construction. An important correction to the earlier plan:
I3 has U0=0 but NOT the earlier plan's extra M0=Omega0=0 premise. Its reference radial memory
gives V0=-mu and Omega0=-mu^2/(2X), mu>0. Zero memory is unnecessary for
the inverse; that source must instead be retained in the pressure and fifth
moment through the whole leading axial support. Velocities vanish beyond
the positive-order patch, while pressure can extend through Xv. The exact
first augmented-background remainder gains q^(2h), but still has the singular
sufficient bound C_m q^(-3/2+2h-m). This is not the full wave residual.

The [background-summation checkpoint](results/background_summation/README.md)
now supplies the conditional general-order induction, common higher stress
support, and fixed-derivative finite residual losses. It retains the radial
axial-viscosity term entering pressure at order two. A logarithmic cutoff
choice and the exact streamfunction curl product rule yield a smooth
divergence-free background for q>0, with a flat AUGMENTED error and a core
swirl lower bound surviving all background orders. The argument requires
finite coefficient bounds at each order, not uniform bounds in order or
convergence of the raw series. Actual coefficient values and cutoff scales
have not been numerically materialized. The bare momentum residual still
contains the annular stress divergence, so this does not give smooth forcing.

The [physical-wave interface](results/physical_wave_interface/README.md)
now supplies exact real curl and full covariance formulas, the complete
physical linearized residual, the uncancelled mean system and a conditional
moving-normal pulse propagator argument. For real primary cosine carriers,
the mean covariance error from the curl is quadratic in inverse frequency.
That exact cancellation does not bound the actual amplitude derivatives.
Two-coordinate physical fixtures and a separate reference pulse check the
algebra. Neither fixture is the actual summed-background band field. The
mean pressure calculation also retains axial defects and possible inner
pressure constants. Kinematic covariance matching alone can leave a
frequency-squared viscous residual; no actual full residual gain is claimed.

The [actual-band pulse checkpoint](results/actual_band_pulse_bounds/README.md)
now derives the phase/frame comparisons from the inherited background
bounds, with whole-pulse rounding costs and one common small-q domain.
An invariant ratio interval gives a positive radial pulse lower bound.
Fixed-order derivative estimates, Gaussian covariance concentration and
a quantitative cone inverse give conditional primary and signed edge
bounds. The full linear wave operator gains a positive epsilon power,
with pressure and viscosity retained. Exact real curls also give a relative
primary covariance gain. These statements depend on the prior continuum
estimates and support geometry; finite fixtures do not evaluate the actual
global field or its constants. The physical mean solve and complete
nonlinear residual cycle remain open. Covariance and linear gains alone
do not give smooth forcing or a complete-field blow-up proof.

The [compact-mean checkpoint](results/compact_mean_corrections/README.md)
now supplies the conditional radial-path, pressure, temporal and signed
five-equation mean maps. It derives the complete mean balances from every
covariance entry and retains the moving axial pressure-bump moment. Exact
recomputation of the same velocity yields a conditional compatibility-defect
gain, with radial memory, viscosity and the new velocity square retained.
Actual primary mean initialization, all wave/mean interactions and the
complete residual cycle are still open. The preceding pulse checkpoint's
32 tests and frozen evidence reproduce; its conditional scope is unchanged.

The [finite-cycle checkpoint](results/residual_cycle_closure/README.md) now
constructs the conditional initial state and closes all three residual
order inequalities through the complete four-step cycle. The proof uses
the full-curl longitudinal identity, retains every signed covariance
remainder and mean interaction, and tracks sources after auxiliary-time
inversion fills the torus. All-stage affine bounds give a 1/10 normalized
order gain, with a common domain and fixed physical jet losses. The
leading/profile/support foundation is still conditional, low-stage
physical residual bounds can remain singular, and the infinite full-field
realization, forcing, energy and actual blow-up remain open.

The subsequent [full-potential realization](results/full_potential_realization/README.md)
supplies the conditional infinite LOCAL field: sum complete potentials,
retain all cutoff curls and compare with a single fixed finite stage.
The complete local residual is flat in every physical jet, with no raw-series
convergence assumption or uniform stage constants. Annular support preserves
the inherited core growth, and the radial/axial energy integral is locally
finite since D-2h>0. These are conditional continuum implications plus
manufactured identity checks. The closed-time auxiliary atlas, admissible
global forcing, whole-space energy and solution-exclusion argument remain
unverified, as does independent review of the inherited foundation.

1. Review the higher-jet derivation and inherited foundation, retaining its
   exact compact domain and fixed N. Extend estimates when needed by the
   complete physical residual; do not treat singular upper bounds as decay.
2. Review the all-order background induction and asymptotic summation against
   their dependencies. Keep the common collar, the lower-order pressure
   sources, arbitrary coefficient growth and the curl-cutoff commutator.
3. Review the new conditional actual-band phase/pulse/linear estimates and
   their exact dependencies, including frozen integer choices, the common
   support geometry and temporal zero extension required for sourced linear
   residuals. Do not replace continuum inputs by the huge-band fixtures.
4. Verify the signed higher-stress corrections and weighted edge estimates
   at the wave interface. Inverse edge-distance powers remain; small q alone
   does not prove uniform cone stability at the flat stress boundaries.
5. Review the new conditional initialization and finite-cycle proof against
   the inherited profile, pulse and support premises. Preserve the exact
   radial pressure/flux identities, flat-term bookkeeping and fixed-stage
   quantifiers; its normalized gain is not a fixed-scale contraction.
6. Review the full potential sum and fixed-stage comparison, then check the
   closed-time auxiliary atlas and its one-sided potential jets away from
   the singular point. Localize and initialize the SAME field with every
   curl, pressure, viscosity and nonlinear term retained. Establish smooth
   compact forcing through T and the whole-space energy/comparison argument.

These steps require continuum estimates. Manufactured crosschecks reject
omitted terms and numerical cancellation errors but cannot supply those
estimates. The exact C2 counterexample exposes a missing bound, not an
irreparable obstruction to the construction.

## A six-hour research block

These are focused-work budgets, not promises that a mathematical argument
will close on schedule. Begin with the finalized leading profile at
research commit `c55ca1e6a19a452bb406841e180d7dcb0f827656`, retaining the fixed
N and all previous conditional dependencies.

| Time budget | Work and reviewable output | How it advances the blow-up argument |
|---|---|---|
| Hours 0–2 | Differentiate the actual moment, phase-map and repair equations to the required angular order. Record the final radial derivative costs and either an explicit third-angular moment bound and annular K, or the precise input preventing that bound. Include the regular axis representation. | Controls the viscous term missing from the current C2 estimate. This is required before the base residual budget applies to the actual field. |
| Hours 2–4 | Use the derived coupled system to construct first inner coefficients from the finalized leading profile, specify compatible axis traces, and attempt a local existence and derivative bound for F1,U1,Pi1. Retain V1's D+2h term and the radial pressure correction. If the first block remains open, trace the missing inputs instead of treating them as bounded. | A successful correction removes the next singular residual coefficient while preserving incompressibility and radial momentum balance. Formal coefficients alone do not establish existence. |
| Hours 4–5 | Set up the five-moment annular extension. Recompute the corrected residual with its quadratic terms and axial viscosity. If existence and bounds are available, prove the expected extra q^(2h) gain in the augmented residual after separating retained annular stress. At a protected core point, check q^(2h) abs(E1)<=E0/2 from actual bounds. Otherwise retain the corresponding open inequality. | Tests whether the correction is compatible with the exterior and reduces physical error while preserving the leading velocity growth in the finite corrected field. This is not yet a lower bound for the infinite corrected sum. |
| Hours 5–6 | Review the derivation against separate physical-coordinate calculations and deliberate omitted-term controls. Check inherited hashes, reproduce the evidence, and push a scoped checkpoint with the next unclosed inequality. | Makes a successful estimate reproducible and a failed estimate actionable. Test counts do not measure distance to a singularity proof. |

The end-of-session report must distinguish an identity, a manufactured
numerical check, a formal coefficient, and a bound for the actual continuum
profile. If a time box expires with an open argument, publish that gap and
continue useful work that does not assume it is solved. No numerical
threshold is to be relaxed to manufacture progress.

The useful outcome in six hours is a new justified estimate or a sharper
obstruction, with the first correction attempted where its inputs permit.
There is no defensible estimate in hours or percent for the complete proof.
Even success in every block leaves all-order background and wave corrections,
divergence-preserving summation, forcing smoothness through T, and the energy
and blow-up estimates for the same complete field to establish.

### Executed work and the next block

The first-correction checkpoint executes the core part of the historical block. It gives
the C3 transfer before modulation, identifies the extra C4 loop input, derives
and constructs a conditional actual first inner background coefficient, checks
the linear patch moments and supplies a one-correction core lower bound.
The higher-jet checkpoint now adds the conditional final compact derivative
estimate, including C4 inputs and all third-order inverse-phase terms. The
annular checkpoint now adds the conditional global first coefficient, corrects
the I3 radial-memory premise, closes the first signed stress tails and derives
the relative augmented-background remainder gain. At that checkpoint the
next tasks were higher-order background induction and physical waves.
None of these checkpoints completes a decaying full residual or blow-up proof.

The background-summation checkpoint subsequently completes that background
induction and asymptotic realization conditionally. The remaining first
interface is now the retained signed stress and the actual physical waves;
flat augmented background error does not cancel the stress by itself.

| Next focused budget | Concrete output | Why it matters |
|---|---|---|
| 0–2 hours | Check the closed-time label choices, atlas transitions and one-sided jets of the realized potentials on q>=c>0 | Flatness at the singular corner does not establish endpoint regularity elsewhere |
| 2–4 hours | Derive the exact spatial and temporal localization residual for the complete potentials, direct swirl and pressure | Preserves divergence while retaining every new curl, pressure and viscous term |
| 4–6 hours | Construct the smooth compact force extension through T and establish the resulting whole-space energy budget, or record the exact obstruction, then push | Connects the conditional local construction to admissible data and forcing; the global-solution contradiction is a further separate gate |

These are time budgets for research outputs, not completion deadlines or a
measure of distance to a blow-up proof. If either analytic input remains open,
retain it explicitly while progressing on independent equations.

## Accepted leading-profile frequency, conditional on the foundation

The [repair-state estimate](results/repair_state_bounds/README.md) bounds the
first reserved repair patch and gives Ccorr=A^128 and r=A^-128. With the
inherited beta=1000, qstar=10^5 and Cstate=D=H^16, its sufficient criterion is

\[
N>\max\left\{1,
\frac{2(C_{state}+2C_{corr}\beta D)}{\epsilon},
16\beta^2q_*D,\frac{4\beta D}{r}\right\}.
\]

N=1+floor(H^32) meets this criterion analytically, conditional on the earlier
continuum estimates. The epsilon in this criterion is the cone-state
tolerance, distinct from the physical expansion parameter q^(2h).
Large radial frequencies create large higher radial derivatives. Later
physical-scale choices must absorb those constants without changing N or
relying on circular parameter choices.

This work package should produce a derivation, exact equation checks,
outward checks where numerical constants are used, negative controls for
omitted terms, and a clearly scoped result. It does not require a new large
fluid simulation. Reproducibility supports the argument; it cannot replace it.

## The difficult part after the leading stress

Forcing and convergence are major proof obligations, not final formatting.
We need control of all spatial and temporal derivatives of the residual near
T, including mixed derivatives. An all-order argument must respect the
quantifiers: for every fixed derivative order there is a valid bound and a
convergent construction. It need not make infinitely many derivative norms
small with one unsupported finite estimate.

A shrinking core with growing leading velocity is a useful mechanism only
if the exact corrected velocity still grows. Obtain a lower bound at a
specified core point or sequence, and prove that the corrections cannot
cancel it there. Energy and force estimates must concern the same complete
field. Checking those properties on different surrogates is insufficient.

## Resource allocation and decision rules

Prioritize the actual derivative bounds, dependency review and coupled
correction construction. Spend numerical effort on specific
algebraic, integration or stability questions that can reject an erroneous
argument. Do not optimize test counts, extreme parameter values or isolated
vorticity peaks as a measure of progress.

Keep the historical unforced numerical search as a secondary route. Its
published .14 resolution failure remains unresolved. Before resuming, verify
archive integrity and obtain a faster independently checked backend capable
of replaying the same initial field at higher resolutions. A numerical
candidate must survive spatial/time refinement, common-point vector-field
comparisons and an independent solver before it receives more expensive
continuation. Even then it is evidence to investigate, not a proof of infinity.

If the forced construction encounters an irreparable obstruction, document
the exact failed implication and decide whether to modify that mechanism or
return to the unforced search. Do not quietly relax viscosity, regularity,
forcing or the target equation to preserve a success claim.

Report progress by the earliest unresolved milestone and its exact missing
inequality. Avoid calendar promises for a proof. Once an end-to-end argument
exists, seek independent mathematical scrutiny of its most fragile estimates
before claiming a resolution.

## Practical next decision

The subsequent [repair-state estimate](results/repair_state_bounds/README.md)
records a correction-to-state bound and accepts a finite frequency conditional
on the preceding continuum estimates. This supplies the analytic repair and
leading-cone acceptance arguments for milestones 1–2. Milestone 0 remains
open to independent scrutiny; those conditional results are not an end-to-end
proof or a numerically resolved singular flow.

The [physical residual budget](results/physical_residual_budget/README.md)
identifies every unlocalized base term and the first correction system. The
[higher-jet estimate](results/higher_jet_budget/README.md) supplies its actual
compact-background derivative bound conditionally, and the
[annular extension](results/annular_first_correction/README.md) now supplies
the conditional global first coefficient, matching moments and signed stress
tails. The [background summation](results/background_summation/README.md)
now supplies the conditional all-order background and flat augmented error.
The earliest new construction tasks are the actual physical wave covariance,
signed stress corrections and mean interactions. Milestone 3 remains open
for the full field: retained stress is still a source in the bare momentum
equation. Smooth forcing and a lower bound for the complete velocity remain
separate gates. A failure must revise the affected construction.
