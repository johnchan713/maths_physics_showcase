# Plan toward a verified finite-time blow-up

Updated: 2026-10-05, following research commit
`5a94f20dfc1de702a507cc9764bfc4a68ff01ab8` and the first-inner-correction
checkpoint.

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

## Immediate work package: actual residual jets and the first correction

The [physical residual calculation](results/physical_residual_budget/README.md)
now gives exact formulas for the unlocalized base and the first positive
order. It does not pass milestone 3: its constants are conditional on an
actual final-profile derivative bound. The earlier C2 bound cannot supply
the third angular moment derivative in radial axial diffusion.

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

The next derivative work is sharper: the loop's third slow derivative needs
stress-coordinate C3, and those coordinates already differentiate the moments.
Thus original moment C4 input is needed. The new unmodulated C3 transfer does
not close the finalized post-modulation bound. In parallel, the constructed
coefficient can be extended into the annulus only after checking the actual
I3 support/memory hypotheses, higher-order stress and all five tails.

1. Bound the finalized modulated/repaired fields through the radial and angular
   orders entering the exact physical residual. Retain phase derivatives of
   N log X, logarithmic-to-ordinary radial conversion, moving phase maps,
   bump widths and inverse normalizations. N remains fixed.
2. Derive the actual third angular moment bound entering Z_-D Z_0 V. State
   each norm and domain; C2 control of U alone is insufficient. Evaluate a
   valid annular K for the new component bounds, and bound the axis through
   the regular Cartesian coefficients F,U,v0,Pi.
3. Review the constructed inner F1,U1,Pi1 coefficient, its core operator
   envelope and sparse angular-loss bound. Carry it beyond `X=3/Lambda`
   while retaining V1's D+lambda factor and Pi1_X=2F0F1-Omega0/(2X).
4. Extend this coefficient into the annulus, restore all five positive-order
   compatibility moments and retain the required higher-order stress. Prove
   support and exterior properties for the same corrected fields.
5. Recompute the full residual with this coefficient, including quadratic
   products, axial diffusion of the correction and the full wave covariance
   when waves are introduced. Prove actual improvement in absolute physical
   units and after fixed derivative orders.
6. Establish the all-order solvability, derivative constants and summation
   argument. Cutoffs must act on divergence-preserving potentials; include
   their product-rule terms. Finite truncations do not give a smooth force
   through the singular time.

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

The first-correction checkpoint executes the core part of this block. It gives
the C3 transfer before modulation, identifies the extra C4 loop input, derives
and constructs a conditional actual first inner background coefficient, checks
the linear patch moments and supplies a one-correction core lower bound.
It does not complete the final-profile residual estimate or global extension.

| Next focused budget | Concrete output | Why it matters |
|---|---|---|
| 0–2 hours | Original angular C4 moment transfer; then the loop's third slow derivative and moving-phase costs, or a precise failed inequality | Supplies the missing input for the final third-angular viscous bound |
| 2–4 hours | Actual I3 support/memory ledger and extension of the convergent core coefficient, with five linear moment corrections and higher stress tails | Turns the local background coefficient into a compatible global coefficient |
| 4–6 hours | Full residual of that extended coefficient, with quadratic terms, its axial diffusion and the required physical wave interactions; reproduce and push | Tests an actual residual gain for the same field and records the next all-order obligation |

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
now identifies every unlocalized base term and derives the first correction
system. The earliest outstanding construction inequality is an actual
post-modulation derivative bound, including the third angular moment
derivative required by radial axial diffusion. Milestone 3 remains open.
After bounding those inputs, solve and extend the first background correction,
then verify all-order background summation, physical waves, mean corrections,
residual improvement and convergence. Smooth forcing and a lower bound for
the complete velocity remain separate gates. A failure must revise the
affected construction.
