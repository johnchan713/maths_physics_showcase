# Q172: Q2 source audit and Q3 uniform-frame correction

Parent: `e4bfe868abf9d8c1276cb3780a35d0282e3580ea` (Q171).
Status: `Q172-Q2-source-audited-and-Q3-uniform-frame-corrected`.

Q171 correctly required inspection of the continuum sources beneath the Q2
and Q3 joins. This checkpoint performs that inspection without changing the
selected core, width `C^-20`, exact repairs, pressure gauge, tail, carriers,
roots or any historical frozen record.

## 1. Q2 survives the source audit

The ten pieces in Q2 trace to analytical continuum estimates, rather than to
the rational fixtures in its executable audit:

| Piece | Analytical source |
|---|---|
| activation | `compact_pressure_bridge` |
| reduced continuation | `reduced_continuation_margins` |
| shape and ideal prejoin | `shape_join_pressure` |
| five-moment join | `shape_join_pressure` and its exact annular inverse |
| ideal after join and first outer ramp | `regional_cone_join` |
| axial transition | `axial_stress_audit` |
| intermediate ramp and power | `intermediate_decay_audit` |
| compensated exterior | `heat_exterior_audit`, `outer_edge_stability` |

The first-ramp estimate uses the corrected negative-source kernel and yields
the required floor above `.586 XR`. The flat stress edges use the factored
weight `zeta`; no positive absolute edge infimum is invoked.

There is also no missing independent orientation inequality in the compact
box. With `t=-b/a`, `g=F(-a,b)` and
`T=F(p1-a,p2+b)`, exact algebra gives

\[
                 -T\cdot g=aF^2(P_c-v).                 \tag{1}
\]

Thus the proved positive `a,F,Pc-v`, together with the finite coordinate
box, gives a positive normalized directional margin. The compact modulation
and exact five-row repair retain the complete pressure memory before the
exterior provider is applied. Q2 is therefore source-audited relative to
these displayed provider derivations. This does not reprove every provider
from the manuscript's initial equations.

## 2. Q3's compactness sentence is not valid on an infinite tail

Q3 defines `CE=sup L||E||` and `Cd=sup L|d-dref|` and says they are finite
by compactness. Each accepted band has finitely many compact charts, but the
accepted tail contains infinitely many bands. Per-band compactness alone
does not give a uniform cross-band supremum: singleton compact sets carrying
the values `1,2,3,...` are an exact counterexample.

This corrects Q3 section 1. Historical Q3 bytes remain frozen. The error does
not invalidate the field or the constrained ODE; it invalidates only that
compactness-only justification of a band-independent upper propagator.

## 3. Q23 supplies the missing uniform argument on the same field

The later Q23 theorem separates a finite head from an analytic tail. For the
tail it introduces estimate constants derived from the fixed background
bounds and uses

\[
 S^2(\epsilon+\epsilon^2+k^{-1})
 \le {360(4/h)^5\over\ell}.                              \tag{2}
\]

An exact threshold

\[
 \ell_0\ge\max\{16,360(4/h)^5,L_{req}/(2r_0)\}          \tag{3}
\]

makes the error and length rows valid on every later band. The earlier
accepted set is finite, so compactness is legitimate there; the nonzero
rounded tangential carrier gives its normal floor. Taking the maximum of the
finite-head and analytic-tail constants yields a uniform upper propagator.
This uses the unchanged Q2 field and accepted labels.

The same-band provider supplies the positive phase, normal, spectral and
root floors, while the full-path provider supplies all fixed jets. Its
constrained pressure remains

\[
 \pi={i\{n\cdot Kt-n'\cdot t+n\cdot f\}\over km|n|^2}, \tag{4}
\]

including both moving-normal and unprojected source-normal terms. Potential
and full pressure are cut only after the whole-path solve, so both cutoff
commutators remain. With Q23 inserted, the Q3 upper-frame/path conclusion is
valid relative to these named providers; Q3's original compactness rationale
is superseded.

## 4. Limits and next obligation

This checkpoint establishes a justified correction and closes the particular
Q171 Q2/Q3 source-audit gate. It does not prove the fixed compensator and
projector supports/jets, compact or slow weighted inverse domains, generated
source membership, complete covariance/mean identities, or the six-group
physical residual. In particular it proves no unconditional pulse family,
residual improvement, force, exclusion or blow-up theorem.

Next audit the fixed compensator/projector family on this corrected Q2/Q3
tail, followed by the compact and slow inverse weighted domains. Preserve the
full axial velocity and `Pi1/Pi2` pressure memory in every later mean and
source calculation.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q172_q2_q3_source_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q172_q2_q3_source_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q172_q2_q3_source_audit/evidence.json \
  --output fresh.json
```
