# External criticism check, 15 October 2026

This note responds to the reported criticism before advancing the construction.
It does not certify the manuscript, the upstream Lean project or this repo's
conditional construction. No change to a frozen mathematical checkpoint is
justified by the two specific comparisons checked here.

## Sources and what was actually checked

The 6 October preprint by Bastounis, Circelli and Hansen,
[Navier--Stokes lost in translation](https://arxiv.org/html/2610.08144v1),
compares the written argument with its formalization. Its authors expressly
withhold a verdict on the written proof's correctness. This distinction is
essential: a mismatch is not by itself a counterexample to an estimate.

The upstream default-branch head checked through GitHub on 15 October was
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, dated 10 September. This is
the revision cited in the criticism. We inspected these declarations:

* [SmoothFamilyTorusInverse.lean, inverse_finiteJets](https://github.com/openai/NavierStokesAndEuler/blob/f9e8bc5b38b6e212696e8a30e3e91517af887bbd/NavierStokes/SmoothFamilyTorusInverse.lean):
  its input jet order is n+5.
* [R3/PressureFlux.lean, exists_uniform_actual_pressure_flux_bound](https://github.com/openai/NavierStokesAndEuler/blob/f9e8bc5b38b6e212696e8a30e3e91517af887bbd/NavierStokes/R3/PressureFlux.lean):
  its displayed flux bound contains dissipationRoot/R as well as cutoffL6.

These source-level observations confirm the two reported statement
differences. We did NOT compile Lean, trace every upstream theorem dependency,
or check the formal definitions of the final Millennium problem statement.
The existing repo never supplied such a formal certificate: its evidence
flags continue to mark that status false.

## 1. Four versus five derivatives: direct mathematical check

The repo's [compact mean inverse](results/compact_mean_corrections/README.md),
section 4, uses output order m with m+4 torus input derivatives. Here is
an explicit justification, independent of the upstream formal declaration.

Let beta=sqrt(2)-1, N=beta partial_y1+partial_y2 and F smooth and periodic
on the unit two-torus, with zero Haar mean. For k=(a,b) nonzero, set
n=|a|+|b|. The product

\[
(\beta a+b)((b-a)-\sqrt2a)=(b-a)^2-2a^2
\]

is a nonzero integer: equality to zero would force an irrational rational
ratio unless a=b=0. The conjugate factor is at most (1+sqrt(2))*n in
absolute value. Thus |beta a+b|^-1<=(1+sqrt(2))*n.

At least one coordinate has |k_j|>=n/2. Integrating the Fourier coefficient
by parts m+4 times in that coordinate gives

\[
|\widehat F(k)|\le (\pi n)^{-(m+4)}\|F\|_{C^{m+4}}.
\]

Define the inverse coefficient by Fhat(k)/(2*pi*i*(beta a+b)) and set
its zero mode to zero. An output mixed derivative of order at most m has
factor at most (2*pi*n)^m. Each summand is bounded by C_m*n^-3 times
the displayed input norm. There are exactly 4n lattice points in the
l1 shell of radius n. Therefore the absolute sum is bounded by

\[
4C_m\sum_{n\ge1}n^{-2}\|F\|_{C^{m+4}}<\infty.
\]

This proves the requested C^m estimate. Taking arbitrarily high fixed m
proves smoothness; zero mean fixes uniqueness. Parameter derivatives commute
with this fixed Fourier multiplier. Constants may depend on m.

Conclusion: m+5 is a weaker sufficient bound, not evidence that m+4 is
false. The repo's m+4 statement need not be changed. The upstream n+5
theorem is not a formal certificate for this n+4 argument. If a later
stage elects to use n+5, the additional fixed input jet must be recorded;
the current all-fixed-jet classes allow it conditionally without a new
epsilon-power loss. This does not supply missing actual background jets.

## 2. Pressure flux: compare the actual derivations

Use the notation of [our comparison argument](results/global_solution_comparison/README.md):
W=sup_t||w||2, B_R=||phi_R^4 w||6 and
D_R=(integral phi_R^8|grad w|^2)^(1/2), for a fixed pre-singular interval.
The repo's proof is a direct L^(3/2) Riesz-transform argument, not an
invocation of the different upstream Lean flux bound.

The weighted stress main term has norm at most C_T(B_R+1). Its commutator
has L^(4/3) norm at most C_T R^(-3/4): the kernel bound is
|z|^-3 min(|z|/R,1), whose L^(4/3) norm scales as R^(-3/4).
The remaining test factor is R^-1 phi_R^3 w. Interpolation bounds its
L3 and L4 norms by R^-1 B_R^(1/2) W^(1/2) and
R^-1 B_R^(3/4) W^(1/4), respectively. Hence

\[
\left|\int\pi w\cdot\nabla\chi_R\right|
\le C_T\{R^{-1}(B_R^{3/2}+B_R^{1/2})
                +R^{-7/4}B_R^{3/4}\}.
\]

With B_R<=C(D_R+R^-1 W), Young's inequality absorbs the D_R powers
3/2, 1/2 and 3/4 into the dissipation square. The respective radius
remainders are R^-4, R^(-4/3) and R^(-14/5), all decaying.
Our earlier proof separately recovers the pressure gradient from the
equation and finite energy; no pressure decay or global gradient bound
for the competitor is silently introduced here.

Conclusion: the inspected formal flux statement is different, but no missing
term was found in this repo's displayed flux/absorption derivation. This
limited recheck does not certify the full comparison argument or the
conditional field on which its blow-up application depends.

## 3. Additional leads and requirements before advancing

The authors also link an [AI-assisted list of twelve potential comparisons](https://www.damtp.cam.ac.uk/research/afha/anders/Navier-Stokes_Experiment.pdf).
It explicitly warns that it has not been fully manually checked. Treat it
as an audit queue. Its items are not twelve established mathematical errors.
The following local safeguards were checked in the existing repo:

| Topic | Local status and remaining obligation |
|---|---|
| Localization potential | Our formula fixes a particular A and B and keeps grad(c) cross A. Changing to A-J, B+curl J changes the localized velocity by grad(c) cross J; equality before cutting cannot justify equality afterward. Keep the chosen representatives in the audit. |
| Domain margin | The localization note chooses tau0=qstar/16 and z0=(qstar/16)^D, obtaining q<=3qstar/8 on the cutoff support. Verify the same domain when assembling the actual field; do not import unrelated fixed upstream cutoffs. |
| Scale derivatives | The next background audit must distinguish (q partial_q)^j from q^j partial_q^j. They are triangularly related, not equal. Transfer any diagonal summation constants explicitly. |
| Axis and reserved mean interval | Trace the actual axis extension and pressure-bump support directly from our defining profile and physical q-based bumps. No claim is made here to have completed that audit. |

Other listed differences in choices or intermediate bounds still require
case-by-case assessment if used. This check has not verified all twelve
upstream comparisons, nor proved that alternative constructions yield
identical fields. All actual profile, background and mean-inverse inputs
remain conditional.

## 4. Analytic forcing is a separate issue

Constantin, Ignatova and Vicol's [regularity paper](https://arxiv.org/html/2609.20803v1)
gives a restriction under its anisotropic bounds and exact-core-axisymmetry
hypotheses: a singular construction of that type cannot have forcing
identically zero on a neighborhood of the event, or the stipulated local
spatial analyticity of the force. This note has not independently verified
that paper's proof or every hypothesis for our eventual field.

Our route requires C-infinity forcing, not analytic or unforced dynamics.
Flatness of every force jet at one event does not mean the force is zero
on a neighborhood. The elementary smooth function exp(-1/x^2), extended
by zero at x=0, has all derivatives zero there and is positive at every
nonzero x. Similarly, our flat-residual claim must never be upgraded to
an exact unforced neighborhood without a new argument. The existing
localization note makes no such upgrade.

## Decision

Keep the mathematical formulas in the frozen checkpoints. Tighten the
evidence policy and the next audit: a compiling formalization is not a
substitute for statement correspondence, and successful fixtures are not
continuum proofs. Before another downstream advance, audit the actual
background-jet inequality (23) of the pulse replay, explicitly tracking
scale-derivative conventions, the common collar, fixed potentials, domain
margins, axis values and mean-bump support. A genuine failed estimate must
be corrected and its dependent claims reassessed; the current mismatch
evidence alone does not justify declaring the construction disproved.

Validation for this review: the unchanged compact-mean and global-comparison
suites pass 35 and 34 tests respectively. These checks reproduce their
existing evidence; they do not establish the upstream Lean compilation,
semantic correspondence or the complete construction. This commit changes
the review and research plan only.
