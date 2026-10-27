# Log-renormalized actual-source tail induction

Parent: `985189b8dae86e2df7a4d3bfe94a592d20bbaae7`.
Status: `log-renormalized-actual-source-tail-induction-proved-with-residual-replay-frontier`.

The parent reconstructs the shifted compact and fixed slow operators and
then states a pure-power normalized supremum for every actual generated
source as the next obligation.  That formulation conflicts with the earlier
logarithmic-slack correction on this branch: fixed-stage coefficients may
carry a finite power of `1+|log epsilon|`, so normalization by the exact
epsilon power alone need not be bounded.

This checkpoint corrects the source space, expands the actual initial source
and one sourced/generic cycle, and closes the corrected source induction at
every fixed finite stage.  Complete pressure memory, all six covariance
entries, the actual axial velocity and the unchanged supports and gauges are
retained.  The result is an analytical weighted-space calculation; the tests
freeze the formula and exponent ledgers but do not prove the continuum
estimates by computation.

## 1. Correction of the parent source frontier

The parent asks, schematically, for

\[
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 \max_{|\alpha|\le R}\sup_{\Omega_{j,\ell}}
 \zeta^{-1}\delta^{p_\alpha}|D^\alpha F_{j,\ell}|<\infty. \tag{1}
\]

Equation (1) is stronger than the established primitive domain.  If
`epsilon_l=exp(-ell)` and

\[
 F_\ell=\epsilon_\ell^{b}(1+|\log\epsilon_\ell|)\zeta,
\]

then the ratio in (1) is `1+ell` and diverges.  This example does not claim
that a canonical source equals this family.  It shows that finite
logarithmic degree does not imply (1), so (1) cannot be called the exact
remaining condition without separately proving cancellation of every log.

For each fixed stage and finite request the correct domain is

\[
 \boxed{
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 (1+|\log\epsilon_\ell|)^{-P_{j,R}}
 \max_{|\alpha|\le R}\sup_{\Omega_{j,\ell}}
 \zeta^{-1}\delta^{p_{j,\alpha}}
 |D^\alpha F_{j,\ell}|<\infty,}                         \tag{2}
\]

where `P_(j,R)` is finite and may depend on the fixed stage and request.
No uniformity in stage or derivative order is asserted.

The pure-power estimate is recovered whenever a positive exponent reserve
`sigma` is available.  With `s=|log epsilon|` and
`log(1+s)<=sqrt(s)`,

\[
 (1+s)^P e^{-\sigma s}
 \le \exp\!\left({P^2\over4\sigma}\right).             \tag{3}
\]

Thus an `epsilon^sigma` reserve absorbs every fixed logarithmic degree.  A
finite constant cannot absorb the log at `sigma=0`.

## 2. Weighted algebra and the actual initial source

The corrected primitive certificate supplies, at every fixed order, uniform
background gaps, normalized support separation, strict root contraction and
finite weighted profile/root/compact seminorms.  Its leaves obey (2) for
finite log degrees.  The following operations preserve that class:

* a finite sum takes the maximum log degree and pole;
* a product adds log degrees, epsilon exponents and poles; since
  `0<zeta<=1`, every product containing a weighted factor satisfies
  `zeta^r<=zeta` for `r>=1`;
* fixed coefficient differentiation is already included by raising `R` and
  the pole finitely;
* auxiliary Haar projection has norm at most two;
* radial moments are finite because every edge integral of
  `zeta delta^-p` is finite, while multiplying the resulting scalar by the
  fixed compensator/projector restores the pointwise weight;
* the oriented compact maps, zero-Haar inverse and five-profile slow response
  preserve the same weighted/logarithmic class by the parent's direct
  estimates.

The actual first prefix uses

\[
 F=F_{\rm linear}+F_{\rm wave-wave}+F_{\rm mean-wave}. \tag{4}
\]

The inherited exact rate ledger puts every group in (4) at least at
`mu=1-kappa`.  The pressure is recomputed after the actual fast map before
the current scalars enter the slow map.  Hence (4), all six covariance
entries, `g_d-A_1g_d`, and the whole
`partial_Z(J_z+c_rho P)` obey (2).  This uses the actual source formula and
does not insert a desired outgoing-source callback.

## 3. One sourced/generic cycle and finite-stage induction

For the next sourced block the retained exact rates are

\[
 B=1-\kappa,\quad C_0=2-2\kappa,\quad E=2-3\kappa.
\]

The actual postwave rows and fast fields have order `E`; recomputed
postfast scalars and slow fields have order `C0`; the complete regular
source remains `E`.  Its radial compact source has the coarse order `B`.
Every entry is a finite composition of the operations in Section 2, so it
has a finite log degree and satisfies (2).

The generic signed block uses

\[
 C={3\over2}-\kappa,\qquad D={3\over2}-2\kappa.
\]

Its complete covariance is order `C`, covariance error `2B`, pointwise
pressure change `D`, post-signed zero-Haar rows and new fast fields `D`,
fresh scalars and slow fields `C`, and complete regular source `D`.  The
covariance formula retains principal cross, old-error cross, complete-curl
cross and signed square.  The pressure/moment formulas retain all six
entries and the entire moving product.  Therefore this block also maps the
class (2) to itself.

The fixed-stage backward recurrence has finitely many nodes and decreases
the stage on every recursive edge.  Its derivative ceilings, pole archive
and harmonic set are finite for each fixed stage.  Induction over that graph
therefore proves (2) for every actual source requested at a fixed
`(j,m,N)`.  Constants and log degrees may grow arbitrarily with these fixed
indices; no infinite-stage uniformity is used.

## 4. Reserve for the compact tails

The existing compact integration count is strict.  To avoid spending an
unknown fraction of its margin, increase each requested count by one.  This
adds the explicit source gain `epsilon^kappa`, or `q^(h kappa)` after
`epsilon=q^h`.  Allocate half to (3) and retain half as strict margin:

\[
 q^{h\kappa}(1+|\log\epsilon|)^{P_{j,R}}
 \le C_{j,R}\,q^{h\kappa/2}.                            \tag{5}
\]

The extra integration requests two additional source derivatives and is
finite.  It changes no field, support, gauge, width or band.  The same rule
applies separately to `A0` and `A1`; the distinct operator constants remain
archived rather than assumed monotone.

## 5. Established scope and next obligation

Equations (2)--(5) correct and close the actual-source tail induction needed
by the direct compact/slow maps for every fixed stage.  They do not yet replay
the final residual estimate.  That replay must retain all six groups

\[
 \text{time},\quad\text{viscosity},\quad\text{pressure},\quad
 \text{old-new},\quad\text{new-old},\quad\text{new-square},              \tag{6}
\]

and assign positive reserve to every logarithmic coefficient, cutoff and
interaction term.  Only after that calculation should localization, compact
force, energy and fixed-gauge comparison be replayed.  No pure-power source
bound at zero reserve, unconditional complete residual improvement, compact
force theorem, comparison theorem or blow-up theorem is claimed here.

## 6. Reproduction

```sh
python3 research/navier_stokes_cascade/results/core_log_renormalized_actual_source_induction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_log_renormalized_actual_source_induction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_log_renormalized_actual_source_induction/evidence.json \
  --output fresh.json
```
