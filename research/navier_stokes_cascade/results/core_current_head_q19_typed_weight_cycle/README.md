# Q19: a typed radial-weight invariant for the finite cycle

Parent: `3ab3fa4cd0e5c85e8343ea7b10f7356c76c6f2d9` (Q18).
Q18 rejects the blanket full-edge source promotion. This checkpoint repairs
that **radial-weight induction**, conditionally on the retained mapping and
root hypotheses, using the exact angular-mean equations. It does not
reinstate Q5's blanket assertion or certify all actual coefficient bounds.
No field, root, pressure gauge, band, cutoff or width is changed.

## 1. Statement and hypotheses

Use E_theta from Q18: at each fixed requested jet, a bound by
C epsilon^a (1+|log epsilon|)^K zeta^theta d^-p, with finite constants,
log degree and poles. The epsilon order a is tracked separately, not
inferred from theta. The two types are:

* W: nonzero angular harmonics, radial weight sqrt(zeta), with the original
  Gaussian P and full lifted-path f/P and endpoint conditions;
* M: angular-independent mean increments, assigned compact stress and
  signed targets with radial weight zeta. Auxiliary torus dependence is
  permitted; angular mean and auxiliary Haar mean are different maps.

Assume the original primary waves are W, initial compact means and assigned
stress are M, and all fixed-order background multipliers have their retained
bounds (finite poles allowed). At every requested finite stage also retain:

1. The same original root a_c has a_c^2 >= c zeta/s0 and the corresponding
   all-jet reciprocal-root bounds. The signed map uses this root, not an
   updated amplitude. Matrix/frame denominators keep their stated gaps.
2. The harmonic inverse maps its original W source domain to W response,
   including complete potential/curl and moving/source-normal pressure.
   This is the `physical_pulse_transfer` section 3 contract, requiring the
   whole lifted-path normalized source jets, not a value norm alone.
3. The shifted T_e/A_e, zero-Haar temporal inverse, oriented target
   primitives, fixed projectors and reserved response map preserve M on
   the requested pole domains with finite log costs. Their exact moment
   hypotheses, smooth support extensions and common-chart identifications
   are retained. The fixed interior profiles restore zeta after a moment.

Then every finite sourced/signed cycle preserves W for nonzero harmonics
and M for mean increments and signed targets. No compact inverse on an
arbitrary half-weight input is needed. This is an implication from these
hypotheses, not their certification. In particular it does not independently
prove the epsilon orders in Q16 or the all-band constants at Q14's orders.

## 2. Why the mean rows have full weight

Use the exact cylindrical decomposition from `compact_mean_corrections`,
(2): u=B+m+w, with angular-independent B,m and <w>_theta=0.
All operators there have theta-independent coefficients and preserve
angular frequency. Therefore a linear term in w or its wave pressure has
zero angular average, including every derivative and cylindrical
connection. The complete curl still consists of nonzero angular harmonics.
This is not the incorrect assertion that auxiliary Haar removes every
harmonic; the angular projection is taken first.

In the mean equations the wave enters only through

    W_ab=<w_a w_b>_theta.

Each of its six entries has weight zeta P^2, hence M after P<=1. Linear
mean terms have weight zeta; quadratic mean terms have zeta^2; derivatives
keep these exponents with finite poles. Products with B do not remove the
weight. The separately assigned stress Sigma must itself satisfy M; no
unweighted base/exterior error is smuggled into this induction.

In particular the displayed radial source gr, Etheta and Ez in the exact
mean system are M once the compact mean pressure is M. There is no circular
pressure assumption: gr depends on means and covariance, not p_m. First
form gr, then Pscalar=int Haar gr, then

    p_m=T0(gr-rho Pscalar).

The source and compensator are M by hypothesis 3, hence so is p_m. Its
axial derivative then enters Ez. Keep the exact remainder
- rho Pscalar - A0(gr-rho Pscalar); zero Haar of A0 is not zero A0.
Similarly retain gamma=gd-A1 gd in every subsequent product.

The current mean rows and their weighted moments therefore feed the
fast/slow maps within M. The radial moment itself need not retain any
pointwise edge weight: its multiplication by a fixed interior profile is
what restores M. Slow potentials can fill the entire reserved hull; the
same oriented/compact bounds, rather than bump-by-bump support claims,
are required. Q17 supplies the polynomial projected-label count.

## 3. Signed target and covariance update

The actual current target is Xe=-Ie Qe Haar(Ee), with Qe the same moving
moment projector. Sections 1--2 put it in M, using the original oriented
primitive estimate. Derivatives of sigma_e, rho and c_rho remain in this
formula. The signed principal correction is H^-1 Xi/(2 epsilon a_c).
The inverse-root calculus loses half an edge weight, so M targets produce
half-weight signed amplitudes; multiplication by the original homogeneous
pulse supplies its retained Gaussian and endpoint factors. Those pulse
bounds are still required. The factor epsilon^-1 changes the epsilon ledger,
not the radial weight. All fixed reciprocal-root and H^-1 derivative
requests remain necessary.

For U the current complete wave and V the complete signed addition,

    Delta W=<U tensor V+V tensor U+V tensor V>_theta

has full weight because every term has two half-weight factors. This
holds for every covariance entry and both cross orientations. It also
holds in the refined decomposition Lambda+Ecov of
`core_generic_signed_covariance_join` (2); no bound for an error summand
is inferred by cancellation from a bound for U alone. Reconstruct all
pressure and scalar changes before the fresh mean maps. Those inputs
are again M by section 2. The selected principal cancellation is relevant
to epsilon gains, but is unnecessary for this radial-weight conclusion.

The entire moving axial moment epsilon partial_Z(Jz+c_rho Pscalar) is
retained. Total pressure includes old Pi1/Pi2 memory and the nonzero wave
pressure. Only the compact mean increment is M; the wave pressure remains
W. Even a purely normal source can have zero response velocity and nonzero
wave pressure. Thus the full fixed-gauge pressure difference is safely
assigned half weight, not full weight. The common base cancels only in
same-base differences, never by deleting its individual axial gradient.

## 4. Sourced step and induction

The exact regular nonzero source consists of the linear wave remainder,
complete wave-wave transport and both mean-wave orientations. Their radial
weights are respectively 1/2, 1 and at least 3/2. The sum is W, precisely
as in the regular-source foundation. A zero angular output of a quadratic
product goes to M; it is not sent through the nonzero harmonic inverse.
All actual cutoff-curl pieces stay in these products.

Under hypothesis 2 the sourced response, potential, complete curl and wave
pressure remain W. Its covariance changes are M. Section 2 then applies
to the fresh sourced means, and section 3 to the next signed step. Induction
over any finite number of these steps proves the radial invariant. Finite
sums/logarithmic degrees may grow with stage; no uniform-in-stage estimate
or convergence of an uncut infinite series is asserted.

Additive Gaussian endpoint and assigned base/exterior residual tails remain
outside the forward inverse, as required by the original grammar. The
argument assumes the retained supported-source/full-path contract; radial
weight counting cannot prove it. In particular it does not replace an f/P
bound with a bound for f or transport Gaussian values between different
points of a lifted path without the inverse estimate.

## 5. Consequences and remaining obligation

The radial part of Q18's obstruction is repaired by routing each source to
its correct space. Q5's full-weight formula remains invalid as a blanket
consequence; the corrected harmonic formula is Q18 (3), while the mean
and signed-target formulas retain full weight. Q15 treats both positive
weights and Q16/Q17 can use their unweighted outputs. Their conditional
physical exponents are not altered by this radial type correction.

Next verify the retained harmonic inverse's whole-path normalized-source
jets and its band-uniform constants for the actual recomputed W sources,
then connect the separate epsilon orders to every Q16 component at Q14's
fully expanded derivative/pole requests. Neither the executable type ledger
nor finiteness of that request graph proves these inequalities. Endpoint
closure, unconditional residual improvement, forcing and blow-up remain open.

Exact Laurent-polynomial tests distinguish angular projection from auxiliary
Haar, preserve opposite-frequency covariance, detect a missing cross term,
and exercise the weight invariant and its failed half-target premise.
Parent evidence and inherited source hashes are replayed. These diagnostics
are not continuum pulse enclosures.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q19_typed_weight_cycle/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q19_typed_weight_cycle/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q19_typed_weight_cycle/evidence.json --output FRESH.json
```
