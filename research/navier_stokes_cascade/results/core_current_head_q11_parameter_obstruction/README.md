# Q11: parameter identity blocks the current-head promotion

Parent: `d1aea3715e41ae2db7440190d6a8f4a69f5ee733`.
Status: `same-object-parameter-interface-obstructed`.

This checkpoint corrects the claimed same-object Q6--Q10 join. It does not
change historical frozen records, rebuild the field, or disprove any
Navier--Stokes existence theorem. The obstruction is to the recorded proof
interface, not a counterexample to the intended continuum estimate.

## 1. Two different values of h

Q1's README and protocol fix `h=lambda^2` and retain the reference roots.
Its root_foundation provider fixes

\[
 T=e^{64}+10\ge64,\quad \lambda=e^{-4T},\quad
 h_{\rm root}=e^{-8T}\le e^{-512}<1/513<1/200.
\]

The strict bound uses only `exp(x)>1+x` for positive x. There is no floating
underflow, approximation of the giant parameter, or replacement of positive
h by zero.

Q6 `complete_residual.join` calls `current_residual.comparison` with five
positional arguments and no h argument. That provider defaults to `h=1/200`.
Q7 sets `H=Fraction(1,200)` and explicitly describes this as the selected Q6
object. Q8 repeats that parameter choice. In contrast, Q5's reserve example
at 1/200 can be read as a diagnostic of a parameterized estimate; the error
is promoting that example to the fixed Q1 object's parameter identity.

A shared source hash and positive sample exponents cannot establish
`h_root=h_replay`: the displayed exact bounds exclude that equality.
Consequently the same-field instantiations in Q6--Q10 require repair before
promotion. Their general conditional residual, localization, energy and
comparison implications are not thereby refuted.

## 2. The frozen comparison stage cannot simply be reused

Write `delta=1/2-3*kappa`, with the replay's `kappa=1/100000`. In the
complete-potential diagonal ledger, at residual order m,

\[
 L=7+2h+2m,\quad H_{loss}=8+2h+2m,
\]

so each nonlinear cross has the bound exponent

\[
 E_{cross}(h,J,m)=\frac{h\delta(J+1)}2-15-4h-4m.
\]

Q6 freezes `J=12817` for `m=N=0`. At its diagnostic h=1/200 this gives
`20030773/20000000>1`, exactly the published result. The slope in h is
positive at this J. But at the root parameter,

\[
 E_{cross}(h_{root},12817,0)
 <E_{cross}(1/513,12817,0)<0.
\]

Thus this ledger cannot certify even the requested exponent one at the
unchanged frozen stage after substituting the actual h. A negative upper-
bound exponent does NOT prove that the actual residual grows: cancellations
or stronger estimates could exist. It refutes the proposed certification
using this bound and this stage.

The exact necessary and sufficient stage condition for this cross bound is

\[
 J\ge\left\lceil\frac{2(N+16+4m+4h)}{h(1/2-3\kappa)}-1\right\rceil.
\]

For any fixed positive h and positive delta it is finite. All other residual
groups, source requests, log reserves and cutoff-selection conditions must
also be checked. Increasing J alone is not a proof of those inputs.

## 3. Endpoint geometry does not repair parameter identity

The current_endpoint_jets observation cylinder already restricts active bands
above the admitted threshold and uses the exact heat exterior near positive-
radius equatorial points. Thus absence of a historical first-band certificate
alone need not obstruct that particular cylinder. Nevertheless the needed
closed-time physical jets still concern one fixed field and its full pressure.
The endpoint cover cannot identify two different h values, nor prove uniform
closed-time bounds merely from interior smoothness or residual flatness.

Next: propagate the exact root parameter symbolically through the Q5 source
requests, Q6 stage/log accounts and Q7 chart, distinguishing numerical fixtures
from selected-field constants. Audit all scale and cutoff choices for that
same parameter before reevaluating the punctured-endpoint bound. Keep the
physical characteristic shift, covariance/axial terms, full pressure memory,
unit viscosity and every cutoff/interaction term. Do not reselect the roots,
width C^-20, gauges, profiles or actual field to fit an arithmetic fixture.

## 4. Evidence and validation

The audit inspects the actual Q6 call and Q7 literal via AST, checks the Q1
and root protocols, reproduces the frozen cross exponent and rehashes the
Q10 transitive evidence/source closure. Tests compare the independent loss
formula with the existing provider's pure arithmetic functions and verify
sharp stage thresholds for several rational parameters, including a tiny
positive fixture. These are exact arithmetic/interface/provenance checks,
not numerical simulations or continuum existence proofs. Historical scopes
and evidence are preserved, even where superseded by this correction.

At parent inspection on 29 October 2026 both push workflows succeeded:
[focused](https://github.com/johnchan713/maths_physics_showcase/actions/runs/40742854912)
and [broad](https://github.com/johnchan713/maths_physics_showcase/actions/runs/40742854759).
The previous empty connector response was insufficient to infer absent CI:
that connector filters to pull-request-triggered runs.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q11_parameter_obstruction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q11_parameter_obstruction/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q11_parameter_obstruction/evidence.json --output q11_fresh.json
```
