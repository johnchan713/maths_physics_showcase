# Q14: restore the primitive expansion omitted by Q12

Parent: `f71c542e46948caebbf44d62981451ed615d8c75` (Q13).
Status: `compressed-interface-to-primitive-correction`.

Reviewing Q13's actual representative premise exposed a regression in Q12.
The quadratic formula there is correct for the historical recurrence, but
its endpoint is a first-real **interface**, not an original primitive.
The already-existing `core_certified_finite_stage_domain_join` explicitly
corrected this distinction. Q12 inadvertently bypassed that composition.
Its historical README, field names and evidence remain frozen; this new
checkpoint supersedes their interpretation as fully expanded primitive
requests. Q13's cutoff differentiation and Q12's exponent arithmetic are
unaffected. No continuum representative estimate follows from these counts.

## 1. Correct normalized and enlarged physical requests

At the retained kappa=1/100000, Q12 proved the interface identity

    I_j(n) = n+500038+99994*j^2+400035*j.

The actual first-real provider `high_jet_join.normalized_join(L)` requests
its complete radial graph through L+400009. Its axial graph through L+200018
is lower. The radial source is requested through L+399995; the downstream
nested graph adds fourteen derivatives. These are translation-invariant
finite call graphs: every order is L plus a fixed offset, so the identity
holds for every nonnegative integer L, not just tested examples.

Consequently the fully expanded primitive ceiling is

    P_j(n) = I_j(n)+400009
           = n+900047+99994*j^2+400035*j.

Q12's enlarged physical interface I_phys already carries the two extra
A0/A1 integration derivatives. Its primitive request is I_phys+400009.
The joint request is

    max(I_j(m), I_phys)+400009.

One must add the offset after the maximum (equivalently to both branches),
not add it twice. This correction is exact for the retained recurrence and
is conservative as a primitive request; it does not assert minimality.
For j=0,n=0 the count is 900047, and for j=1,n=0 it is 1400076,
matching the certified provider. The physical checks compare the enlarged
Q12 request with the certified physical request plus two derivatives.

All formulas also apply to Q12's symbolic root-parameter stage
J=ceil(K/h)+1, with h=exp[-8(exp(64)+10)] and
K=2(N+17+4m)/(1/2-3*kappa). No evaluation or rounding of h or J is needed.
The primitive ceiling remains finite for each fixed request. No constant
uniform in that enormous order is claimed.

## 2. Propagate the correction into pole domains

For original pole p and retained pole q_old, the provider rule is

    q_new=max(q_old,6*p,6*P+9,6*P+12*p+6),

where P is the fully expanded primitive order, not the smaller interface I.
If the order-dependent branch controls before and after correction, the
common pole increases by exactly 6*400009=2400054. If q_old or 6p controls,
the increase can be smaller or zero. In all cases it lies in [0,2400054].
This follows because the two order-dependent arguments increase by the
same amount while the other maximum arguments stay fixed.

At each resulting distinct q, retain the original fast/axial operators at
6q, pressure and moment operators at (6q,12q), and oriented operator at 6q.
The proof argument contains the sum of the archived constants. A larger
pole number does not imply a larger operator norm; replacing the archive
by only the largest-pole constant is unjustified. The implementation keeps
all distinct requested domains and supports a compressed symbolic sequence
of requests without claiming to enumerate the root-stage archive.

No actual root, field, profile, support, gauge or cutoff sequence is changed.
This expands the proof-domain requests. The full covariance, axial velocity,
moving pressure moments and both pressure memories stay in the providers.

## 3. Effect on the representative frontier

The independent bound still needed by Q13 is, for every fixed j,k, finite
closed-support, all-physical-band constants in

    |D^k(A_j,b_j e_theta,pi_j)|
      <= C_(j,k)(1+|log q|)^P_(j,k) q^[h*delta*j-(2A+2k)].

The first-real primitive contracts must be invoked at the expanded ceiling
and all corresponding archived poles. Finite graph composition at the
smaller interface order cannot certify those inputs. The conditional
closed-chart coefficient join already names this provider route; Q14 restores
its correct counts to the new symbolic calculation rather than assuming the
representative inequality from test success or source hashes.

Next: bind the actual complete representative formulas and their normalized
all-band seminorms to these corrected domains, keeping the exact root
parameter and full pressure. The remaining issue is the substantive
representative bound and support identification, not a numerical evaluation
of giant jets or cutoff scales. Endpoint closure, unconditional residual
improvement and blow-up remain unproved.

## 4. Validation

Checks compare compressed requests against the live certified state and
physical providers, including both extra integrations; exercise huge rational
fixture stages without a stage-by-stage traversal; verify the pole shift and
archive; and regenerate Q13 evidence and the inherited hash graph. Rational
fixtures are not substitutes for the transcendental root. Historical evidence
is unchanged. Q13's focused CI had succeeded on 29 October 2026; its broad
workflow was still running at inspection.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q14_primitive_ceiling/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q14_primitive_ceiling/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q14_primitive_ceiling/evidence.json --output q14_fresh.json
```
