# Certified fixed-stage high jets and the diagonal obstruction

Parent: `c4762504fef22a4683afbdc81f1961e0ce62fc26`.
Status: `certified-fixed-stage-high-jet-domain-join-with-diagonal-obstruction`.

The parent closes the complete first-sourced compact request graph on the
certified later tail.  The historical `core_finite_stage_compact_recurrence`
already supplies a finite backward recurrence for every later fixed stage,
but calls its terminal counts “first-real primitive orders.”  Those counts
actually stop at the **first-real interface**.  This checkpoint corrects the
label, composes every terminal request with the certified first-real provider,
and recomputes the common-pole archive from the fully expanded orders.

The result is an analytical fixed-stage composition.  No giant derivative,
factorial, operator constant, or thousand-stage field is evaluated.  Passing
tests checks exact arithmetic and frozen provenance; it is not a numerical or
independent continuum proof.

## 1. Retained finite-stage graph

Keep

\[
 \kappa=10^{-5},\qquad \mu=1-\kappa,
 \qquad B_j=\tfrac32-2\kappa+j\delta,
 \qquad \delta=\tfrac12-3\kappa>0.                 \tag{1}
\]

One later cycle consists of a sourced wave and fresh means followed by a
generic signed wave and fresh means.  The actual traced backward offsets are

\[
 d_s(j)=2p_{0s}(j)+18,\qquad
 d_g(j)=2p_{0g}(j)+17,                               \tag{2}
\]

with the strict integer counts frozen by the historical recurrence.  At
`j=0`,

\[
 d_s(0)=200018,\qquad d_g(0)=300011.                 \tag{3}
\]

All six covariance entries, complete pressure and moments, actual `A1` axial
velocity, fresh fast and slow means, unit-mass radial compensation, whole
supports, and growing finite harmonic set remain in those traces.  Every
recursive edge decreases the stage index, so the graph is finite and acyclic
for each fixed stage and output order.

## 2. Correct interface and primitive ceilings

Let `I_j(n)` denote the terminal first-real **interface** order in the
historical recurrence.  Its exact formula is

\[
 I_j(n)=n+500038+\sum_{i<j}\bigl(d_s(i)+d_g(i)\bigr). \tag{4}
\]

The certified first-real provider at output order `L` has conservative
primitive ceiling `L+400009`.  Hence the actual fully expanded order is

\[
 P_j(n)=I_j(n)+400009.                                \tag{5}
\]

For the first two normalized states,

| stage | first-real interface | primitive ceiling |
|---:|---:|---:|
| 0 | 500,038 | 900,047 |
| 1 | 1,000,067 | 1,400,076 |

Thus the inherited values 500038 and 1000067 were finite and useful, but
were not primitive ceilings.  Historical frozen evidence is preserved; the
correction appears only in this new checkpoint.

For the illustrative fixed physical request `m=N=0`, `h=1/200`, the same
composition gives

| stage | first-real interface | primitive ceiling |
|---:|---:|---:|
| 0 | 60,500,042 | 60,900,051 |
| 1 | 60,900,077 | 61,300,086 |

Both `A0` and `A1` branches retain strict physical exponent margin
`1/20000000`.  More generally, for every fixed `j,m,N`, the certified
all-fixed-order contracts and the finite recurrence imply

\[
 |D_{\rm physical}^{m}T_{{\rm compact},j}|
       \le C_{j,m,N}q^N.                               \tag{6}
\]

The constant may depend on the stage, jet, requested power, and fixed
profiles.  Equation (6) is not uniform in any of those indices.

## 3. Correct common-pole archive

The common-pole rule must be applied to the primitive order `P`, rather than
the smaller interface order `I`:

\[
 q_{\rm new}=\max(q_{\rm old},6p,6P+9,6P+12p+6).      \tag{7}
\]

At fixed stage, archive the `6q/12q` pressure and fast-map constants and the
separate oriented `6q` primitive constant for every distinct `q`.  The proof
argument contains their sum; it does not assume that a larger-pole operator
constant dominates an earlier one.

For the stage-two illustrative physical request, the fully expanded
primitive order is 61,900,109.  With original pole one, (7) gives

\[
 q=371400672.                                         \tag{8}
\]

The lower normalized-state requests introduce no additional distinct pole
in this example.  Their request-to-pole entries are still recorded, while
the operator archive stores the one distinct value once.  The harmonic
bound remains `M_2<=16`, and inherited closed-support extension remains
valid.  Neither the fields nor the constants at (8) are materialized.

## 4. Fixed-stage residual consequence and exact diagonal gap

The inherited finite-stage residual calculation retains viscosity, full
pressure gradients, cutoff commutators, old/new interactions, new squares,
and assigned Gaussian tails.  Conditional on those exact identities, a fixed
stage has

\[
 |\mathcal R(U^{[j]})|_m
 \le C_{j,m}q^{hB_j-(3/2+2h+2m+\eta)}.                \tag{9}
\]

For `h=1/200`, `m=N=0`, and reserve `eta=1`, the first stage whose exponent
is at least zero is `j=1002`; stage 1001 is still negative.  This is exact
rational threshold arithmetic.  It neither constructs nor numerically
evaluates 1002 fields.

The remaining obstruction is now quantitative.  If a proposed diagonal
uses `j=j(q)` and seeks a `q^N` bound, (9) requires

\[
 \boxed{\log C_{j(q),m}
 \le \bigl(hB_{j(q)}-(3/2+2h+2m+\eta)-N\bigr)
          \log(1/q).}                                  \tag{10}
\]

No available estimate bounds the stage growth of `C_{j,m}`, the archived
operator constants, or the corresponding complete-potential increments
strongly enough to prove (10).  The fixed-stage statement

\[
 \text{for every target there exists a finite conditional state}
\]

therefore cannot be exchanged for one `q`-dependent realized field.  This
is an insufficiency of the present upper bounds, not a lower-bound theorem
that every possible construction must fail.

The next obligation is to control the stage growth in (10) and converge one
complete-potential diagonal.  Only on that same field can the full viscosity,
cutoff, interaction, localization, force, energy, and core-growth estimates
be joined.  Uniform pulse bounds, a complete residual improvement, admissible
forcing, and blow-up remain unproved.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_certified_finite_stage_domain_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_certified_finite_stage_domain_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_certified_finite_stage_domain_join/evidence.json \
  --output fresh.json
```
