# Q12: symbolic root-parameter stage, source and cutoff audit

Parent: `efb5d0777f22aa31cd78b6fad0dbd54a33636856` (Q11).
Status: `root-parameter-arithmetic-repaired-cutoff-identification-open`.

Q11 found that Q6--Q10 promoted a diagnostic parameter to the Q1 field.
This checkpoint repairs the parameterized arithmetic interface without
editing those frozen records. It proves the bounds below for every real
`0<h<1/100`, at the retained `kappa=1/100000`, and hence for the exact root

\[
h=\exp[-8(\exp(64)+10)]>0.
\]

The root is never evaluated, rounded to zero, or replaced by a rational
fixture. These are analytical arithmetic and request-graph results. The
actual cutoff constants and closed-time physical jets are not certified.

## 1. A fixed symbolic stage that closes all six exponent accounts

Fix nonnegative integer residual jet m and power N. Put

\[
\delta=1/2-3\kappa,\quad B_0=3/2-2\kappa,\quad
K=\frac{2(N+17+4m)}\delta,\quad
J_0=\lceil K/h\rceil,\quad J=J_0+1.
\]

These are finite integers, independent of q. There is no claim that their
decimal expansions or the actual stage fields have been computed. Set
`r=h J0-K`, so `0<=r<h`. The Q6 providers' regular residual after spending
its own half-stage log reserve has exponent

\[
E_{reg}=\delta K-5/2-2m+\delta r+h(B_0+\delta/2-2).
\]

The cutoff tail has exponent `E_tail=h delta (J+1)/2-7-2h-2m`.
Time, unit viscosity and pressure each use E_tail; both nonlinear crosses
use E_cross = E_tail - (8+2h+2m); the correction square uses 2 E_tail.
In particular

\[
E_{cross}=N+2+\delta r/2+h(\delta-4)>N+49/25>N+1,
\]

and

\[
E_{reg}>2N+63/2+6m-50007/20000000>N+1.
\]

The other groups are larger because `8+2h+2m>0` and E_tail>0.
Also J0>=m+2. Thus the *exponent ledger* closes at the root parameter,
unlike the frozen stage 12817. This is a conservative sufficient stage,
not the minimal stage and not a proof of its continuum source hypotheses.
The complete exponent is still the minimum of the fixed-state and six
tail exponents, never a sum of unrelated reserves.

## 2. Compact source requests without a list of J preceding stages

Write `mu=1-kappa`, `B_J=B0+delta J`, and

\[
T_{m,N}=(N+3/2+2m)/h+2,\qquad
p_e=\max\{1,\lfloor(T_{m,N}-b_e)/\kappa\rfloor+1\},
\quad (b_0,b_1)=(\mu,B_J).
\]

Then `h(b_e+p_e kappa-T)>0`. At the selected sufficient stage,
`h(B_J-T)>=N+65/2+6m-5 kappa h>0`, so p1=1. The A0
floor is positive and gives `0<h(mu+p0 kappa-T)<=h kappa`.
An additional integration for each map changes raw source orders from
`m+2p_e+3` to `m+2p_e+5`. A0 and A1 separately spend `h kappa/2`
on their own finite log degrees and retain `h kappa/2`. Neither pays the
regular reserve `h delta/2` or the cutoff half-gain reserve.

For any finite log degree P and positive reserve sigma, with `s=|log q|`,

\[
(1+s)^P e^{-\sigma s}\le \exp(P^2/(4\sigma)),\quad s\ge0.
\]

Indeed `log(1+s)<=sqrt(s)` and completing the square gives the bound.
Its constant may grow arbitrarily with the fixed stage, derivative order
and root parameter. Uniformity as h tends to zero is neither proved nor
needed for this fixed parameter. The degree P must be fixed after its
finite-stage source formula, before taking a scale supremum.

The existing recurrence computes `state_ceiling` using a list of all j
layers. Its exact closed form at this fixed kappa is

\[
S(j,n)=n+500038+99994j^2+400035j.
\]

To verify it, the sourced integration counts in layer i are
`(100000+49997i,4)` and the generic counts are
`(149997+49997i,50001)`. The actual nested call graphs give backward
offsets `200018+99994i` and `300011+99994i`. Their sum over
`0<=i<j` is `99994j^2+400035j`; the seed is n+500038.
Thus this compression is an identity for the inherited recurrence,
not a truncated induction or an assumption about actual high jets.

The enlarged physical request must also reach the first-real primitive.
For J>=1 set

\[
n_* = m+\max(2p_0+19,2p_1+9),\qquad
L_* = S(J-1,n_*+200018+99994(J-1)).
\]

The generic call graph gives n_* after both extra integrations. The
preceding sourced layer contributes the displayed offset. Request
`max(S(J,m),L_*)` from the original primitive. At J=0 the physical
formula is `n_*+300021`. Tests compare these expressions with the actual
recursive providers, including the two-derivative enlargement. No source
term is dropped: full pressure/moment memory, axial velocity, all six
covariance entries, physical shifts and fixed gauges remain the inherited
analytical obligations. Counting requests does not bound their values.

## 3. The remaining same-field cutoff inequality

Let c_j>=max(0,log C_j) and P_j>=0 bound *every* cutoff-uniform closed-chart
constant and log degree requested at stage j (jets through j). For the
**already selected** cutoff sequence, put s_j=log a_j. The sufficient
test from the diagonal provider, now at the root h, is

\[
s_j>\max\left\{s_{j-1}+1,\log(q_*^{-1})+1,1,
\frac{4(c_j+j)}{h\delta j}+\frac{4P_j^2}{h^2\delta^2j^2}\right\}
\qquad(j\ge1).
\]

Here a finite list of (c,P) pairs means taking the maximum over the list.
With gamma=h delta j/2, for s>=s_j the log of the target majorant is at
most `c+P^2/(2 gamma)-gamma s/2<-j`. Hence
`C(1+|log q|)^P q^gamma<exp(-j)<2^-j` on `0<q<=1/a_j`.
The formula is sufficient, not necessary; failing this test does not prove
that an actual tail is large.

The algebra establishes existence of admissible thresholds if the constants
are finite. It does **not** identify them with the selected field's a_j.
Q12 neither reselects a_j nor substitutes a newly constructed field. Before
restoring the same-field join, exhibit the original selection rule at the
root h, verify its actual closed-chart constants and supports, and verify
the displayed inequalities for its sequence. A diagnostic choice made at
1/200 does not certify these inequalities: the two terms scale as h^-1
and h^-2. If the original sequence was symbolically chosen at root h with
these full constants, document that identity rather than choosing it again.

The quantifiers must remain: one h and one field/cutoff sequence first;
then, for each (m,N), a finite comparison stage J and a sufficiently small
q neighborhood (including `q<=1/(2a_J)`). Different comparison stages do
not authorize different physical fields.

## 4. Chart propagation and remaining limits

Use the same h in `D=1/2-h`, `q=tau+z^2 q^(2h)` and
`L=1-2h eta^2`. The whole closed chart has L>=1-2h>49/50.
The inherited transfer `q<=4(tau+|z|^(1/D))` is valid in this regime;
thus no h=1/200 is needed in the geometry. This only transfers a bound
once its same-field closed-chart hypotheses are established. It does not
prove the punctured terminal all-order velocity/potential/pressure bound.

Next obligation: audit the original diagonal selection and its complete
cutoff-uniform physical constants/supports at root h against Section 3.
Then revisit the punctured endpoint inequality. Q11's same-field promotion
block remains open; no unconditional pulse, complete residual, force,
comparison, exclusion or blow-up theorem is claimed.

## 5. Reproduction and provenance

The code stores the root symbolically and proves rational coefficient
bounds for the whole h interval. Separate rational fixtures include very
small positive values; they are regression checks, never approximations
of the root. Q10's transitive graph plus Q11's record are rehashed, and Q11
is independently regenerated. Historical evidence is unchanged.

At parent inspection on 29 October 2026 Q11's [focused workflow](https://github.com/johnchan713/maths_physics_showcase/actions/runs/40943649504)
had succeeded; its [broad workflow](https://github.com/johnchan713/maths_physics_showcase/actions/runs/40943649713)
was still running. CI is a regression check, not a continuum proof.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q12_symbolic_parameter/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q12_symbolic_parameter/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q12_symbolic_parameter/evidence.json --output q12_fresh.json
```
