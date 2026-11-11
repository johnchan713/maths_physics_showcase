# Q171: provider-discharge frontier and selector obstruction

Parent: `bff4c7188acf7c2dd41046f64e0d5029c27bf3a8` (Q170).
Status: `Q171-current-head-provider-frontier-reconciled-with-selector-obstruction`.

This checkpoint reconciles the conditional Q164--Q170 current-head chain
with the older `comparison_foundation` scope flags and the lower provider
frontier. Historical evidence is preserved. The correction is about which
statement is current; it does not rewrite an old conditional record as an
unconditional theorem.

The conclusion is a provider-frontier audit and an exact implication for the
common selector. Tests verify provenance, logic and exact arithmetic. They do
not establish the continuum infima or whole-path bounds identified below.

## 1. What the later audits actually supersede

At its own head, `comparison_foundation` correctly recorded

- `all_retained_analytic_dependencies_discharged = false`, and
- `retained_exact_core_and_finite_existence_inverse_gap_lemmas_still_dependencies = true`.

Those bytes remain true descriptions of that historical checkpoint. They are
not, however, the current Q170 frontier. On the byte-identical successor
object:

1. Q164 independently rederives the all-index core operator estimates,
   contraction, original equations and leading same-gauge pressure identity.
2. Q165 binds the finite regional, frame, projector and log-reserved inverse
   chain to that core, relative to the named Q2--Q5 continuum providers.
3. Q166--Q170 successively replay the six-group residual, punctured endpoint,
   localization and force, fixed-gauge comparison, and protected core limit.

Thus the old local phrase “core and finite maps remain” is stale as an
interface frontier. The later audits do **not** make every theorem used by
the Q2/Q3 providers independently proved. Q165 itself says that it is a
same-object join relative to named continuum providers and that exact
fixtures and hashes are not proofs of actual continuum infima or selected
field values.

## 2. First premise not independently discharged by the joins

After Q164 fixes the exact core profile, the first independent promotion gate
is the actual Q2 regional and Q3 whole-path/frame continuum input on that same
profile. In the notation already frozen by the primitive-frontier audit, it
requires

\[
 m_{\rm base}=\min\{f,t,d,\omega,m_\phi,c_a,b_0,
                     \lambda_{\min},b_{\min}^2\}>0,                  \tag{1}
\]

together with

\[
 0<C_i<\infty                                                       \tag{2}
\]

for every value, relative-stress, chart, length and full-frame row. Here the
flat-edge requirement is the valid weighted/factored estimate; no positive
absolute stress infimum at a zero edge is inserted.

The Q2 and Q3 current-head records bind provider chains asserting (1)--(2),
but their joins explicitly do not independently reprove every theorem in
those chains. Q171 therefore neither declares those provider assertions
false nor promotes them to independently certified premises. It records
(1)--(2) as the first place where an unconditional promotion still needs a
continuum proof audit.

## 3. Exact selector implication and obstruction

Only after (1)--(2) may the shared later-tail selector choose an integer
`ell` satisfying

\[
 2h(\ell-2)\ge \max_i\log_2(C_i/r_i),
 \qquad 0<r_i\le m_{\rm base}.                                     \tag{3}
\]

For the unchanged `h=1/200`, (3) is equivalent to the exact rational
inequalities

\[
 2^{\ell-2}\ge (C_i/r_i)^{100}.                                   \tag{4}
\]

Hence positive rational margins and finite positive constants imply a finite
selector without a floating logarithm. The executable sample takes
`r=(1/2,1/4)` and `C=(4,2)`, so every ratio is `8` and the least selector is
exactly `ell=302`.

Conversely, if any actual `r_i` is zero or is not proved positive, or a
corresponding `C_i` is not proved finite, the right-hand side of (3) is not a
finite real number. Taking `q` smaller cannot create the missing positivity
or finiteness. This is the precise blocking inequality.

## 4. Scope and next obligation

Q171 establishes that the Q164--Q170 interface chain is coherent and that
the older comparison flag is not the current interface frontier. It does not
prove the actual Q2/Q3 continuum premises, a uniform actual-pulse estimate,
an unconditional complete residual theorem, force admissibility, competitor
exclusion or Navier--Stokes blow-up.

The next justified checkpoint is an independent source-level audit of the
actual Q2 ten-region continuum infima, orientation/cone margins and factored
flat-edge coefficients on the exact Q164 profile, followed by the Q3
whole-path/frame constants. It must either prove (1)--(2) from the displayed
continuum estimates or identify the first specific row that lacks such a
proof. Only then may (3) be rerun and the downstream chain promoted.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q171_provider_discharge_frontier_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q171_provider_discharge_frontier_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q171_provider_discharge_frontier_audit/evidence.json \
  --output q171-fresh.json
```
