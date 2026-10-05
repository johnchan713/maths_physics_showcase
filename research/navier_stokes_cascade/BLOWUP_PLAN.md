# Plan toward a verified finite-time blow-up

Planning checkpoint: 2026-10-05, based on research commit
`b822f89792dd408a6e55f33f2de737f08de6dfb9`.

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

## Immediate work package: the correction patch

Use the first reserved repair patch already selected in the construction:
X=rhop exp(y), 0<y<5, with input U=0 and
E=e(eta) exp[(-1/2-lambda)y]. Keep lambda strictly positive. Use the existing
two U bumps and three E bumps, with their disjoint supports and exact integrals.

1. Write the exact perturbed fields in terms of the five coefficient functions.
   Derive a coefficient neighborhood where E remains positive uniformly.
2. Differentiate the actual bumps. Retain their inverse-width factors and
   angular derivatives of e and of the coefficients. Bound the resulting
   changes in a and b with the perturbed E in the denominator.
3. Integrate all five moment changes, including quadratic terms. Bound their
   first angular derivatives, then propagate them through the original pressure
   formulas. Incoming modulation errors must remain in the comparison domain.
4. Prove Ccorr and r in the norms used by the acceptance criterion. Test the
   proposed Ccorr=A^128 and r=A^-128; replace them with derived bounds if needed.
   Their size alone is not a proof.
5. Combine with the exact small-root condition. Preserve the inverse-lambda
   loss when normalizing and combining the M/J rows. Establish exact matching
   of all five moments beyond the patch, not just agreement on a finite grid.
6. Review the frequency criterion against these bounds and the original cone
   gaps. Accept N only after all inputs are justified. Otherwise retain the
   failed condition and revise the construction.

The existing constants for this criterion are beta=1000 and qstar=10^5;
they also need their stated hypotheses checked. With verified constants, the
sufficient condition from the prior stress note is

\[
N>\max\left\{1,
\frac{2(C_{state}+2C_{corr}\beta D)}{\epsilon},
16\beta^2q_*D,\frac{4\beta D}{r}\right\}.
\]

The latest work bounds Cstate and D by H^16 conditional on the preceding
profile estimates. The proposed N=1+floor(H^32) is still unaccepted.
Large radial frequencies create large higher radial derivatives. Once N is
fixed, later physical-scale choices must absorb those constants without
changing N or relying on circular parameter choices.

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

Prioritize the analytic repair and dependency review, followed by the exact
residual and correction construction. Spend numerical effort on specific
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

The next useful result is either a valid repair estimate enabling the leading
stress construction, or a concrete failure that changes that construction.
That is the nearest measurable step toward the target. Full PDE correction,
smooth forcing and a blow-up lower bound remain beyond it.
