# Certified first-sourced high jets and compact tails

Parent: `592487a3f4acb055097a346636af3d7efc209533`.
Status: `certified-first-sourced-high-jet-and-compact-tail-domain-closure`.

The parent closes the first-real, pre-sourced `A0/A1` high-jet graph on the
certified later tail.  The next state adds the first sourced complex wave,
its complete covariance and pressure changes, and fresh fast and slow means.
The historical `core_first_sourced_physical_tail_closure` already proves the
corresponding finite acyclic implication, but leaves its primitive contracts
conditional.

This checkpoint joins that implication to the now-certified first-real
provider and the retained first-sourced construction.  It also corrects an
important bookkeeping ambiguity: the historical ceilings 300021 and
60400023 are orders requested at the **first-real interface**.  They are not
fully expanded primitive ceilings.  Composing them with the current
first-real fixed-order graph gives conservative primitive ceilings 700030
and 60800032, respectively.

The result is confined to the first sourced state.  It does not close the
later stage-indexed recurrence or the complete physical residual.

## 1. Keep one changed state and all pressure memory

Let `R` be the closed first-real state and let `v` be the first sourced cut
curl, with zero entrance data and the full forced pressure.  The changed
state `S` is constructed in the following acyclic order:

1. form the complete supported regular source of `R`;
2. solve the moving-normal constrained equation, retaining the unprojected
   normal source in pressure, and cut both potential and pressure;
3. form all six entries of

   \[
   \Delta W=\operatorname{Bil}(w_R,v)+\operatorname{Cov}(v); \tag{1}
   \]

4. reconstruct \(\Delta P,\Delta p_m,\Delta J_\theta,\Delta J_z\), including
   the whole moving product \(c_\rho\Delta P\) and old pressure memory;
5. form the postwave tangential rows and their `Qcirc` parts;
6. apply the fresh fast map, including the actual axial velocity
   `g_d-A1 g_d` computed from the pre-fast postwave row;
7. recompute pressure and all three scalar defects;
8. apply the same reserved five-map slow response; and
9. form the final radial source, unit-mass compensation, and only then apply
   `A0`.

The desired `A0` output is last and is never used to construct its own
source.  The postwave `A1` response precedes the new fast field but depends
only on the pre-fast row.  Thus the changed-state graph remains acyclic and
retains the exact compact residual

\[
 T_{\rm compact,S}
 =(-A_0Q^\circ g_{r,S},0,+A_1Q^\circ E_{z,\rm postwave}). \tag{2}
\]

Neither input in (2) may be replaced by its first-real predecessor.

## 2. Normalized high-jet closure

Write

\[
 B=1-\kappa,\qquad D=\tfrac32-2\kappa,
 \qquad C=\tfrac32-\kappa,qquad \kappa=10^{-5}.          \tag{3}
\]

At normalized output order zero, the retained strict compact counts give

| Branch | Integrations | Raw source order |
|---|---:|---:|
| radial `A0` | 50001 | 100005 |
| axial `A1` | 2 | 7 |

The historical changed-state graph requests the first-real provider through
order 300021.  The parent proves that the current first-real graph at output
order \(n\) has conservative combined ceiling \(n+400009\).  The fully
expanded primitive ceiling is therefore

\[
 300021+400009=700030.                                   \tag{4}
\]

No derivative, factorial or profile seminorm at (4) is numerically sampled.
The all-fixed-order continuum contracts and finite positive recurrences imply

\[
 \mathcal N_{100005}(Q^\circ g_{r,S})
 \le \epsilon^{1-\kappa}P_{100005}(S+K_{\rm all}),       \tag{5}
\]
\[
 \mathcal N_{7}(Q^\circ E_{z,\rm postwave})
 \le \epsilon^{3/2-2\kappa}Q_{7}(S+K_{\rm all}).         \tag{6}
\]

Both compact outputs reach order \(3/2\), hence exceed the required
\(C=3/2-\kappa\) by the strict margin \(\kappa\).  Equations (5)--(6) use the
changed covariance, pressure moments, axial velocity and refreshed means;
they do not mix states.

## 3. Fixed physical jets and powers

For every fixed physical jet \(m\), requested power \(q^N\), and retained
\(\epsilon=Q^h\), the same construction produces finite strict compact
counts.  For the illustrative fixed pair \(m=N=0\), \(h=1/200\),

| Branch | Raw source order |
|---|---:|
| radial `A0` | 60,200,007 |
| axial `A1` | 60,100,009 |

The first-real interface ceiling is 60,400,023.  Expanding it through the
current provider gives the conservative primitive ceiling

\[
 60400023+400009=60800032.                               \tag{7}
\]

Both physical exponent margins are exactly \(1/20000000\).  Therefore, for
each fixed pair,

\[
 |D_{\rm physical}^{m}T_{\rm compact,S}|
 \le C_{m,N}q^N.                                        \tag{8}
\]

The constant in (8) may depend on \(m,N\) and the fixed profiles.  This is
not one estimate uniform over all derivative orders, stages or scales.

## 4. The complete residual obstruction remains

The synchronized regular nonzero source at `S` has order
\(D=3/2-2\kappa\).  With \(\Lambda=3/2+2h\), its zero-jet physical exponent is

\[
 hD-\Lambda=-\frac{15025001}{10000000}<0.                \tag{9}
\]

Thus flatness of (2) does not prove decay of the complete residual.  After a
positive polynomial reserve \(\eta\), a later synchronized regular order
\(a\) can imply a requested \(q^N\) bound only when

\[
       ha>N+\Lambda+2m+\eta.                             \tag{10}
\]

Equation (9) is an insufficiency of the available upper bound, not a lower
bound for the residual.  The next obligation is the stage-indexed backward
high-jet and common-pole recurrence for every later sourced state, keeping
the growing harmonic set, complete pressure memory, all `A1` velocities and
whole supports.  Each compact graph must then be paired with its matching
regular state before a controlled physical diagonal is attempted.

Exact controls freeze the interface/primitive distinction, request indices,
rational margins and the negative regular exponent.  They test algebra and
provenance; the analytical conclusion follows from the cited finite-order
composition theorems.  Uniform pulse estimates, all later harmonics, the
complete viscosity/cutoff/interaction residual, forcing and blow-up remain
unproved.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_first_sourced_high_jet_domain_closure/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_first_sourced_high_jet_domain_closure/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_first_sourced_high_jet_domain_closure/evidence.json \
  --output fresh.json
```
