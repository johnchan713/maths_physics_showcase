# Q193: uncut Duhamel and phase-partition recurrence

Parent: Q192, `c1c00c61fdb34d01808b3c17d732b6427658a3d6`.
The live branch, Q191--Q192, the physical-pulse and historical actual-band
providers, their frozen evidence, and the active Q189--Q192 CI runs were read
first. The Q189 selected profile, complete pressure formula, supports,
frequency and unchanged `C^-20` activation width remain fixed.

Q192 closed the auxiliary rational frame transfer and restored the physical
`epsilon^-kappa` graph loss. This checkpoint carries that corrected ledger
through the fixed-coordinate differentiated Duhamel recurrence and proves the
full Faà di Bruno phase-exponent comparison. It stops before applying the
cutoffs and taking the complete curl.

## 1. Fixed-coordinate Duhamel recurrence

Divide the homogeneous nonzero-harmonic frame solution by its positive
Gaussian envelope `P`. The normalized propagator is bounded on the whole
closed pulse. If an order-`i` coefficient derivative costs `S^p_i`, one
Duhamel integral costs at most one additional `S`, so

\[
 b_0=0,\qquad
 b_j=\max_{1\le i\le j}\{p_i+1+b_{j-i}\}.            \tag{1}
\]

The historical provider uses the sharper input `p_i=i`; (1) then gives
`b_j=2j`. Q192 independently established the conservative frame degree

\[
 p_i=2i+4\quad(i\ge1).                              \tag{2}
\]

Using (2) in (1), induction gives

\[
 b_j=7j.                                            \tag{3}
\]

Indeed the `i=1` term gives `7+7(j-1)=7j`, while every other `i` gives
`2i+5+7(j-i)=7j-5(i-1)<=7j`. Thus every fixed auxiliary derivative has a
finite polynomial loss

\[
 |D^jt_m|\le C_jS^{7j}P.                            \tag{4}
\]

The exponent `7j` is deliberately nonsharp. It is independently sufficient
because the primary theorem permits a finite derivative-dependent power of
`S`. Q193 does not silently import `2j` as an independently proved sharp
coefficient ledger.

For a zero-data source with degrees `a_j`, multiplication by the projected
frame gives

\[
 g_j=\max_{0\le i\le j}(p_i+a_{j-i}),
\]
\[
 b_0=g_0+1,qquad
 b_j=\max\left\{g_j+1,
 \max_{1\le i\le j}(p_i+1+b_{j-i})\right\}.         \tag{5}
\]

Equation (5) is recorded only as a conditional mapping rule. Q193 does not
assert that every generated nonlinear source has the required full harmonic
bounds or smooth endpoint extension.

## 2. Every Faà di Bruno phase partition

For `a` spatial and `b` time derivatives, a term with `j` nonempty phase
factors has `1<=j<=a+b`. Before the carrier is inserted, the differentiated
phase factors cost

\[
 Q^{-a/2-b(1+h)}.
\]

The carrier satisfies `k<=2Q^{-h/2}`, so the term's total loss is

\[
 L_{a,b,j}={a\over2}+b(1+h)+{jh\over2}.             \tag{6}
\]

The selected target exponent is

\[
 T_{a,b}=a\left({1\over2}+{h\over2}\right)
       +b\left(1+{3h\over2}\right).                \tag{7}
\]

Subtracting (6) from (7) gives the exact nonnegative margin

\[
 T_{a,b}-L_{a,b,j}={h\over2}(a+b-j)\ge0.            \tag{8}
\]

The worst partition `j=a+b` saturates (8), so there is no undisclosed spare
power in the phase estimate. Fixed nonzero harmonics change constants, not
the domain or the `Q` exponent.

## 3. Conditional uncut primary conclusion

Combining (4) and (8), for each fixed request the homogeneous uncut frame
amplitude has

\[
 |D^{(a,b)}t_m|\le C_{a,b}S^{B_{a,b}}P,             \tag{9}
\]

and its oscillatory factor fits the unchanged exponents

\[
 s_x={1\over2}+{h\over2},\qquad
 s_t=1+{3h\over2}.                                  \tag{10}
\]

Fixed-coordinate Duhamel differentiation adds no new `Q` loss. Pulse-coordinate
derivatives follow from the ODE and retain a finite polynomial loss. The
primary pressure keeps

\[
 \pi={i\over km|n|^2}
 \left(n\!\cdot K_{bg}t-n'\!\cdot t\right),         \tag{11}
\]

while the `n.f` term must be restored for every sourced inverse. The positive
lower bound remains a statement about the uncut radial amplitude relative to
`P`; it is not a pointwise lower bound for the oscillating cutoff velocity.

## 4. Remaining boundary

Q193 has not assembled the physical fields. The next obligation is to put the
actual transverse, pulse and shell cutoffs on vector potentials and pressures,
take the complete cylindrical/Cartesian curl, and prove the velocity,
pressure and potential suprema and smooth closed-support extension on the same
domain. Radial graph coefficients, differentiated cutoffs, cylindrical
connections, axial velocity, the pressure denominator, and the one-sided
`T=0` endpoint must all remain present.

The named provider supplies that primary assembly relative to its inputs, but
Q193 does not independently reproduce it. No unconditional actual-pulse
hierarchy, complete residual improvement, admissible force, comparison result
or Navier--Stokes blow-up is claimed. Tests verify exact recurrences, exponent
arithmetic and frozen provenance; they are not a continuum proof.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q193_duhamel_phase_recurrence/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q193_duhamel_phase_recurrence/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q193_duhamel_phase_recurrence/evidence.json \
  --output q193-fresh.json
```
