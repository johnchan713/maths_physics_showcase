# Q13: cutoff-uniform derivatives and the actual representative frontier

Parent: `18d8c52172856032ef031ad4302db0b791acbae7` (Q12).
Status: `all-order-cutoff-factor-proved-representative-bound-open`.

The original complete-potential diagonal, Section 3, and its certified
successor, Section 2, already select scales symbolically for arbitrary finite
lists of actual constants. Their selection rule agrees with Q12 after
substituting gamma=h delta j/2. No numerical list of a_j is required for
that conditional existence implication. Q12's missing numerical evaluation
is therefore not an additional mathematical obstruction. Identification still
requires the same h, complete representatives, fixed gauges and constants;
the diagnostic h=1/200 never becomes the root parameter. This checkpoint
does not select a replacement field or cutoff sequence.

## 1. Exact closed-chart differential operators

Keep the root h=exp[-8(exp(64)+10)], D=1/2-h, and work on q>0,
tau=1-t>=0, with q=tau+z^2 q^(2h). Define eta=z/q^D and s=aq.
Then |eta|<=1 and L=1-2h eta^2>=49/50. Implicit differentiation gives

\[
q_z=2\eta q^{1-D}/L,\quad q_t=-1/L,\quad
\eta_z=q^{-D}(1-\eta^2)/L,\quad
\eta_t=D\eta/(qL).
\]

Consequently, for any smooth F(eta,s),

\[
\partial_z(q^\beta F)=q^{\beta-D}
 [2\beta\eta F+(1-\eta^2)F_\eta+2\eta sF_s]/L,
\]
\[
\partial_t(q^\beta F)=q^{\beta-1}
 [-\beta F+D\eta F_\eta-sF_s]/L.
\]

Start at beta=0,F=chi(s). After r z derivatives and l time derivatives,
with n=r+l, induction gives

\[
\partial_z^r\partial_t^l\chi(aq)
 =q^{-Dr-l}\sum_{k=0}^n c_k(h,\eta)s^k\chi^{(k)}(s).
\]

Every c_k is a polynomial in h,eta divided by L^(2n); redundant powers
of L are allowed. The recurrence is explicit: a z step replaces c_k by
`[2(beta+k)eta c_k+(1-eta^2)c'_k+2eta c_(k-1)]/L`, and a t step
replaces it by `[-(beta+k)c_k+D eta c'_k-c_(k-1)]/L`.
Missing coefficients are zero. These formulas retain the derivative of
the denominator; dropping it would invalidate mixed-order estimates.

For a fixed smooth chi equal to one on s<=1/2 and zero on s>=1,
all positive derivatives are supported in 1/2<=s<=1. For n>0 the
coefficient c_0 vanishes. Thus the formula holds globally on q>0,
including both seams by smooth extension, and implies

\[
|\partial_z^r\partial_t^l\chi(aq)|
 \le C_{r,l,\chi}\,q^{-Dr-l},\qquad a>0.
\]

The constant is independent of a and can even be chosen uniformly for
0<=h<=1/100: the coefficient functions are continuous on that compact
(h,eta) rectangle with L>=49/50. This is an all-order induction, not an
inference from sampled derivatives. No claim concerns q=0 itself.
Cartesian x,y derivatives of this cutoff vanish. The h=0 endpoint used
in diagnostics is an auxiliary limit, not a replacement for positive h.

## 2. Joining the cutoff to complete representatives

Let Z_j=(A_j,b_j e_theta,pi_j) in the inherited fixed gauges. Suppose its
Cartesian mixed jets of total order k<=m+1 satisfy, on the complete closed
physical support,

\[
|D^k Z_j|\le C_{j,k}(1+|\log q|)^{P_{j,k}}
 q^{h\delta j-(2A+2k)},\qquad A=1/2+h,
\]

with finite constants uniform over physical bands. For 0<q<=1 each
cutoff derivative of order d costs at most q^-d since D<=1. Leibniz
then bounds a derivative of order m of curl(chi(aq) A_j) using
k+d=m+1. Its loss is `2A+2k+d<=2A+2(m+1)`.
The direct swirl and pressure use k+d=m and obey the same conservative
bound. Summing finitely many terms therefore proves

\[
|(V_{j,a},P_{j,a})|_m\le\widetilde C_{j,m}
 (1+|\log q|)^{\widetilde P_{j,m}}
 q^{h\delta j-L_m},\quad L_m=2A+2(m+1),
\]

uniformly in a, where V=curl(chi A_j)+chi b_j e_theta and P=chi pi_j.
The commutator grad(chi) cross A_j is included. Subsequent differentiation
also retains pi_j grad(chi), all mean and axial components, and all
pressure/covariance memory in the representatives. No pressure gauge is
changed. Multiplication and differentiation do not enlarge a smooth
representative's closed support, conditional on its inherited zero extension.

## 3. What remains to prove

This discharges the cutoff-factor part of the uniformity required by Q12.
The remaining substantive inequality is the representative hypothesis in
Section 2, for every fixed j,k, on the actual common support, with full
pressure and fixed potentials. Interior smoothness alone does not prove
its supremum finite over all physical bands or at closed angular endpoints.
Nor does a bound on curl A alone bound the fixed-gauge potential A needed
by the cutoff commutator.

Once those constants and supports are proved for the root-parameter field,
the existing symbolic diagonal recipe absorbs them without numerical scale
evaluation or stage-uniform constants. The chosen sequence must be held
fixed across all requested residual orders and powers. Q12 supplies the
comparison stages, not the missing representative theorem.

Next: trace the actual complete potential/direct-swirl/pressure bounds to
the Q1--Q5 providers and establish their closed-support, all-physical-band
suprema (or exhibit the precise missing one). The punctured endpoint bound,
unconditional residual improvement and blow-up remain unproved.

## 4. Checks and provenance

Exact rational polynomial diagnostics verify the recurrence's mixed-partial
commutation, first derivatives, h=0 limit and coefficient bounds. The proof
for all orders and the exact transcendental root is the induction above.
Q12 evidence is regenerated; the full inherited evidence graph is rehashed.
Historical frozen records are unchanged.

On 29 October 2026 Q12's focused CI succeeded; its broad CI was still running
at inspection. Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q13_cutoff_uniformity/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q13_cutoff_uniformity/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q13_cutoff_uniformity/evidence.json --output q13_fresh.json
```
