# Complete dependency and quantifier audit

Parent: `527d758f8270763e4959f2c5e6b1a3e381414a43`.
Status: `complete-common-object-quantifier-audited-with-provider-frontier`.

The parent proves the protected-core lower bound for the complete-potential
diagonal relative to the frozen analytical providers.  This checkpoint audits
the resulting field/force/comparison theorem as one dependency graph.  It
checks object identity and the order of every limiting choice before deciding
whether the conditional conclusion may be promoted.

The audit finds no downstream field reset: the active construction retains
one completed leading profile, one rule-defined Borel background, one accepted
tail, fixed potential and pressure gauges, one complete-potential cutoff
sequence, one localization and one force.  The full transitive frozen graph
has 215 evidence nodes and 837 directed dependency edges and is acyclic.

This does **not** make the result unconditional.  The graph still imports its
all-order continuum estimates through provider contracts.  Hash agreement and
finite exact controls prove provenance and arithmetic, not those analytic
estimates.

## 1. The active object spine

The object-bearing chain is:

1. the completed leading profile and its rule-defined real Borel background;
2. the same-band phase, normal, value and cone tail;
3. the full-path homogeneous and sourced ODEs with moving-normal and
   unprojected pressure terms;
4. the fixed compensator/projector profiles and compact/slow weighted maps;
5. matching-state finite cycles with all six covariance entries, axial
   velocity, pressure moments and supports;
6. the locally finite complete-potential diagonal in fixed gauges;
7. the punctured endpoint atlas, one localization and one compact force;
8. the fixed-gauge expanding-ball comparison; and
9. the protected-core lower bound on that diagonal.

Each later item reaches the preceding item through the frozen dependency DAG.
Historical alternative realizations occur only as analytical providers; they
do not authorize replacing the object selected on this spine.

## 2. Quantifier correction

The fixed-stage records contain derivative ceilings, compact
integration-by-parts counts and common-pole operator-norm archives that depend
on a requested jet `m` and power `N`.  These are proof requests for the same
fixed smooth representatives.  They are not new profile, pulse or field
choices.  Consequently constants may depend on `j,m,N` without violating the
one-field requirement.

The valid order is

\[
 (u^*,p^*,\{a_j\})\ \hbox{fixed first};\qquad
 \forall m,N\ \exists J=J(m,N)<\infty\ \exists q_0>0\
 \forall 0<q<q_0.                                      \tag{1}
\]

The stage `J` is held fixed as `q` tends to zero.  For the conservative
default budget `h=1/200`, `m=N=0`, the existing inequalities give `J=12816`.
There is no uncontrolled substitution `J=J(q)`.

At the endpoint, every compact with `q>=c>0` sees finitely many stages before
the terminal limit.  The future force uses one Borel diagonal for all fixed
jets.  Finally, comparison has the noninterchangeable order

\[
 \boxed{T<1\ \hbox{fixed}\quad\longrightarrow\quad R\to\infty
        \quad\longrightarrow\quad\tau\downarrow0} .    \tag{2}
\]

Thus the request-dependent proof archives and the limit order create no new
same-field obstruction.

## 3. Exact remaining promotion condition

The remaining issue is upstream of (1), not in the comparison limit.  For the
same fixed representatives, the complete proof still needs, for every fixed
finite `j,m,N`,

\[
 \boxed{\mathcal C_{j,m,N}<\infty,\qquad
  |\mathcal R(U^{[j]})|_m
  \le \mathcal C_{j,m,N}
       q^{hB_j-(3/2+2h+2m+\eta)}} .                    \tag{3}
\]

Here `C_(j,m,N)` must be derived from the defining profile, roots and inverse
maps for the same `U^[j]`.  It must retain complete pressure, unit viscosity,
every cutoff, both old/new cross interactions and the new square.  No
stage-uniform bound is required.

The finite constants in (3) require the scale-compatible provider conditions

\[
 \Delta_{bg}>0,\qquad {3\over8012}>0,\qquad
 \vartheta_{ann}<1,\qquad\vartheta_{rep}<1,             \tag{4}
\]

together with finite fixed-order weighted seminorms of the profile, roots,
compact maps and slow maps on their unchanged normalized support hulls.  The
existing checkpoints assert (3)--(4) through frozen analytical contracts and
the later corner replay composes them correctly.  This checkpoint does not
reconstruct those continuum estimates from the raw defining maps, so it
cannot turn those assertions into an unconditional theorem.

Equation (3) is the precise frontier.  A failure would invalidate residual
flatness and smooth force gluing.  A proof for one request-dependent field
would also be insufficient; the constants must belong to the same fixed
representatives used in the diagonal.

## 4. Consequence and scope

Relative to (3)--(4), the already joined argument is coherent: it produces
one conditional smooth compact force with zero initial data and finite
preterminal energy, and the expanding-ball argument conditionally excludes a
global smooth locally finite-energy solution with that same force by using the
protected-core velocity lower bound.

No numerical diagnostic, test pass or hash match proves (3).  This checkpoint
therefore makes no unconditional pulse, complete-residual, admissible-force or
blow-up claim.  The next justified obligation is to reconstruct (3) directly
for the fixed all-order profile and matching finite-stage states, beginning
with the compact/slow inverse weighted estimates and propagating them through
the complete residual ledger.  Only after that reconstruction should the
force and comparison implication be replayed without a conditional provider
flag.

## 5. Reproduction

```sh
python3 research/navier_stokes_cascade/results/core_complete_dependency_quantifier_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_complete_dependency_quantifier_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_complete_dependency_quantifier_audit/evidence.json \
  --output fresh.json
```
