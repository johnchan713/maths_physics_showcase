# Goal: establish or reject a genuine Navier–Stokes blow-up construction

The target is a solution of the three-dimensional incompressible equations,
with positive viscosity and admissible smooth initial data, that becomes
unbounded at a finite time. A forced construction must also prove that its
force has the smoothness and decay required by the
[Clay statement](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf).
A singular force, an unresolved numerical spike, or divergence of an
approximate leading field does not meet that target.

Current priority: finish and critically verify the manuscript-based
construction. Each checkpoint must either close a named mathematical
obligation, construct a needed part of the field, or expose an obstruction.
Test counts and extremely small reference parameters are not measures of
distance to a proof. The historical unforced numerical search remains a
separate route; its `.14` resolution failure has not been repaired.

| Obligation | Current position |
|---|---|
| Reference outer schedule and pulse | Analytic moment and stress bounds established in the preceding audits. |
| Heat exterior, compensation and outer edge | Reference construction and bounds in `results/heat_exterior_audit`; this does not supply a regular axis. |
| Regular axis and annular attachment | Local Taylor machinery exists; a complete matched profile is unverified. |
| Full momentum residual and corrections | Unverified. A leading reference field is insufficient. |
| Smooth admissible force and actual singularity | Unverified. No blow-up result is claimed. |

The next decision is determined by the first unclosed obligation. A failed
condition must remain visible and must change the construction or its stated
scope before further claims are made.
