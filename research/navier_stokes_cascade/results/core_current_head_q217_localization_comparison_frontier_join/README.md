# Q217: localization, compact force, energy and comparison frontier join

Parent: Q216, `21566beab39c60b2100f69805db2cb09488d0ce5`.
The live branch, Q216 and its active CI, Q203, Q206 and the frozen
Q168--Q170/Q182 localization and comparison records were read first.
Historical evidence, the selected leading/background field, roots, complete
pressure gauge, supports, Borel and diagonal scales, accepted tail and
unchanged `C^-20` width remain fixed.

Q216 source-audits the all-fixed-order and all-finite-pole quantifiers used by
the fixed-stage complete residual recurrence. This checkpoint binds that exact
object to the already proved localization, force, energy, pressure-comparison
and protected-core implications. It also isolates the first new hypothesis in
the downstream comparison chain: finite energy of a competitor on each fixed
preterminal interval. Local smoothness, zero data and compact smooth forcing
do not imply that hypothesis by themselves.

The conclusion remains relative to Q216's cited continuum provider theorems.
Executable ledgers check formulas, exact arithmetic and provenance; they are
not formal proofs of those analytical inputs.

## 1. Same source-audited physical object

Q216 consumes Q215 and the reconstructed primitive providers. Q215 joins the
same Q203/Q205 physical family, fixed pressure gauge, complete-potential
diagonal and residual recurrence. Q203 in turn consumes the frozen endpoint
and localization chain replayed in Q168--Q170 and Q181--Q182. The evidence
closure therefore fixes one field, one pressure, one support hole, one
localization cutoff and one force construction. No canonical-family transfer,
new root, changed cutoff scale or parallel pressure normalization is inserted.

For each fixed finite jet and pole, Q216 supplies the finite constants required
before the diagonal scale is chosen. This is sufficient for local finiteness
and fixed-jet residual flatness. It does not supply a constant uniform in
order or stage.

## 2. Exact localization and terminal force trace

With the Cartesian-smooth axisymmetric cutoff `c`, write

\[
 u=\nabla\times(cA)+cB=cu_0+d,
 \qquad d=\nabla c\times A,
 \qquad p=cp_0.                                      \tag{1}
\]

The field is divergence free, has zero initial data and equals the source-
audited field on the protected late core. Its unit-viscosity residual is

\[
\begin{aligned}
\mathcal R_1(u,p)={}&c\mathcal R_1(u_0,p_0)+c_tu_0+d_t
 -(\Delta c)u_0-2\sum_i c_i\partial_i u_0-\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
 +c(u_0\cdot\nabla c)u_0+c(u_0\cdot\nabla)d\\
&+c(d\cdot\nabla)u_0+(d\cdot\nabla)d .             \tag{2}
\end{aligned}
\]

These are exactly twelve nonzero groups. The apparent extra term vanishes
because

\[
 (\nabla c\times A)\cdot\nabla c=0.                 \tag{3}
\]

Equation (2) retains the complete pressure cutoff, all viscous cutoff terms,
both nonlinear cross orientations and the cutoff-curl square. For force jet
order `m`, direct differentiation requests

\[
 A:m+2,\qquad B:m+1,\qquad p:m,                     \tag{4}
\]

and the endpoint Cauchy step requests one more derivative. The frozen endpoint
provider supplies `A:m+4`, `B:m+3`, `p:m+2`, leaving one derivative unused in
each channel after the Cauchy step.

At the Cartesian corner `c=1`; the punctured endpoint provider controls the
other compact pieces. The arbitrary-power residual bound and the
one-higher endpoint jets consequently give every terminal force jet. The
indexwise Borel extension, whose scale is chosen after the finitely many
actual norms at that index are known, produces

\[
 f\in C_c^\infty(\mathbb R^3\times(0,\infty)).       \tag{5}
\]

The individual pressure equation still contains the force divergence,

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f. \tag{6}
\]

## 3. Smooth force does not regularize the terminal field

Q206 applies to this same protected core. With

\[
 A={101\over200},\quad X_{\rm in}={1\over101},
 \quad X_a={4\over101},
\]

the support-hole margin is `3/101`. Along

\[
 x_\tau=(\sqrt{2X_{\rm in}\tau},0,0),\qquad t=1-\tau,
\]

all annular diagonal corrections vanish exactly, while the positive-order
background is small rather than zero. For some `e0>0` and sufficiently small
positive `tau`,

\[
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-101/200}.                    \tag{7}
\]

The two-path pressure difference has the corresponding gauge-independent
lower bound of order `tau^-101/100`. Thus velocity and pressure value jets at
the terminal origin do not exist, even though (5) holds. The force trace is a
cancellation statement for the full residual; it is not a future velocity or
pressure construction.

## 4. Zero-safe energy control

For every `t<1`, the localized field is smooth and compactly supported. If

\[
 E=\|u\|_2^2,\quad D=\|\nabla u\|_2^2,
 \quad F(t)=\int_0^t\|f(s)\|_2\,ds,
\]

then regularizing by `sqrt(E+epsilon^2)` before division gives

\[
 E(t)+2\int_0^tD(s)\,ds\le F(t)^2.                 \tag{8}
\]

No terminal value of `u`, `p` or `E` enters (8). The energy supremum and the
integrated-dissipation estimate are separate consequences; their extrema
cannot be added under a single copy of `F(1)^2`.

## 5. Fixed-gauge comparison and the first new conditional premise

Fix `T<1`. Let `(v,P)` be a locally smooth solution with the same force and
zero data, and put `w=v-u`, `pi=P-p`. The two individual pressure equations
both retain `-div f`; subtracting cancels the common force exactly once:

\[
 w_t+\operatorname{div}g=\Delta w-\nabla\pi,
 \qquad -\Delta\pi=\partial_i\partial_jg_{ij},       \tag{9}
\]

where `g=w tensor w+w tensor u+u tensor w`.

The pressure recovery and expanding-ball proof require

\[
 \boxed{\sup_{0\le t\le T}\|v(t)\|_2<\infty}
 \qquad\text{for every fixed }T<1.                  \tag{10}
\]

Under (10), `g` is uniformly `L1`, its Riesz pressure lies in `H^-2`, and the
momentum equation puts the harmonic pressure-gradient difference in `H^-3`.
A harmonic `H^-3` distribution has zero gradient. The eighth-power cutoff
then gives

\[
 {1\over2}E_R'+{1\over2}D_R^2
 \le\|\nabla u(t)\|_\infty E_R+C_TR^{-2}.            \tag{11}
\]

At fixed `T`, sending `R` to infinity proves `v=u` on `[0,T]`.

Condition (10) is the first new, non-derived premise after accepting Q216's
source-audited construction. It is essential. Indeed

\[
 v(x,t)=t e_1,\qquad P(x,t)=-x_1,qquad f=0          \tag{12}
\]

is a smooth divergence-free Navier--Stokes solution with zero initial data
and compact smooth force `f=0`: `v_t+nabla P=0`, while viscosity and convection
vanish. But `v(t)` is not in `L2(R^3)` for any `t>0`. Thus local smoothness,
zero data and compact forcing alone cannot justify (10), and pressure Poisson
recovery alone cannot discard the harmonic gradient.

## 6. Protected-core limit and exact scope

For each fixed protected `tau`, take

\[
 t_\tau=1-\tau,\qquad T_\tau=1-{\tau\over2}.
\]

The legal order is: fix `T_tau<1`, recover pressure, send the comparison radius
to infinity, and only then send `tau` to zero. With (10), comparison transfers
(7) to the competitor. A globally smooth competitor that also satisfies the
fixed-interval finite-energy hypothesis would be bounded near `(0,1)`, giving
the existing conditional exclusion implication.

This checkpoint establishes no exclusion of all locally smooth competitors.
It does not discharge Q216's cited regional, Volterra, weighted-kernel,
Fourier and finite-composition continuum theorems, produce uniform actual-
pulse constants, regularize the terminal origin, construct a field after
time one, prove an unconditional residual or force theorem, or prove
Navier--Stokes blow-up.

The next justified obligation is to decide the intended competitor class.
Either derive (10) from explicit global spatial decay/finite-energy hypotheses,
or retain it in the theorem statement as essential. Independently, every
cited upstream continuum provider must be discharged before any unconditional
construction claim.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q217_localization_comparison_frontier_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q217_localization_comparison_frontier_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q217_localization_comparison_frontier_join/evidence.json \
  --output fresh.json
```
