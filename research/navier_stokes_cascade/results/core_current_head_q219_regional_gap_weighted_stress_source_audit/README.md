# Q219: regional-gap and weighted-relative-stress source audit

Parent: Q218, `413e9d6c2bfb0dff19d6a49f437bd94da9b4002b`.
The live branch, Q218 and its active CI, Q216, the finalized leading-profile
certificate, Q189, the regional-cone join, the weighted-edge audit and both
same-band gap joins were read first. Historical evidence, selected parameters,
roots, complete pressure gauge, supports, Borel scales and unchanged `C^-20`
activation width remain fixed.

Q218 sends the audit back to the first continuum leaf that Q216 still cites:
the leading regional gap infima and weighted relative-stress estimate. The
source chain does contain continuum arguments for these rows. They use a
finite analytic regional cover, exact flat-edge factorizations, a finite-head
plus geometric-tail Borel argument and one later value-level selector. They do
not infer continuum bounds from sampled fields or rational fixtures.

This is a source audit at the theorem-dependency level. It does not rederive
every primitive regional ODE lemma from the initial equations, numerically
resolve the selected giant profile, or promote passing tests to a continuum
proof.

## 1. Regionwise leading gaps on one full-pressure profile

The regional source covers activation, continuation, shape, exact repair,
axial, intermediate, heat-compensated and outer-edge pieces with exact shared
endpoints. It retains the same completed pressure primitive across every
join. On the activation and continuation input, the recorded powers are

\[
 a>C^{-22},\qquad P_c-v>C^{-401},
 \qquad 2(P_c-v)^2-(v-2)J_c^2>C^{-802}               \tag{1}
\]

where the last two gates are required only on the subset `v>=2`. The common
regional floor is `C^-1024`, while the compact envelope is `A=C^4096`; hence

\[
 22,401,802,1024<4096.                               \tag{2}
\]

The `v<1` continuation and shape pieces are not incorrectly required to pass
a `v>=2` gate.

The older first-outer-ramp proof applied the positive-data kernel lower bound
to a negative source term. The regional audit corrected that sign. Q189 gives
the exact elementary lower bound

\[
 {11\over10}{68\over125}-{2\over98}
 ={17701\over30625}>{1\over2}.                       \tag{3}
\]

Thus the weaker `p1>XR/2` floor actually consumed downstream survives, while
the superseded stronger expression is not reused.

The retained axial and intermediate controls are also strictly positive:

\[
 2(35)^2{47\over50}{417\over1000}>1,
 \qquad 2(12)^2{627\over2000}>1.                     \tag{4}
\]

The finite-frequency compact modulation and exact five-row repair preserve
their own `kappa/2` gaps. The heat exterior uses the compensated profile,
including nonzero pressure memory between separated supports, and the exact
global moment identities identify its backward stress with the same forward
pressure gauge.

Equations (1)--(4) are analytical region estimates and exact comparisons.
Interval arithmetic is used only for outward rounding of fixed scalar
expressions. Endpoint registries test coverage; they are not sampled
supremum proofs.

## 2. The absolute stress floor is false

The leading stress vanishes flatly at both ends of the active annulus, so

\[
 \boxed{\inf_{X_a<X<X_b}|T_0|=0}.                    \tag{5}
\]

Any same-band proof requiring a positive global absolute stress floor is
therefore invalid. The source proof instead uses the exact factorizations

\[
 T_0=e^{-t_1^2/y_a^2}A_a,qquad |A_a|\ge c_a>0       \tag{6}
\]

at the inner edge and

\[
 T_0=e^{-4/y_b^2}y_b^{-3}(b_0,y_b^6b_z),
 \qquad b_0>0                                       \tag{7}
\]

at the outer edge. With

\[
 \zeta=e^{-t_1^2/y_a^2-4/y_b^2},
 \qquad \delta=\min(1,y_a,y_b),                     \tag{8}
\]

overlapping edge collars and the compact middle give

\[
 |T_0|\ge c\zeta,
 \qquad |\partial^I T_0|+|\partial^I T_1|
 \le C_I\zeta\delta^{-3m-3}.                        \tag{9}
\]

The normalized leading direction consequently extends smoothly through both
zero-stress edges. The two exponential constants in (8) cannot be exchanged:
their quotient would have an uncontrolled exponential at one edge.

## 3. Why the relative background stress is valid

Define `W_n=T_n/|T_0|` in the active interior. The first correction vanishes
on an inner collar. At the outer edge, the matching factor gives

\[
 W_1=left({b_1\over\sqrt{b_0^2+y_b^{12}b_z^2}},0\right), \tag{10}
\]

which is smooth even at zero angle. Every higher stress has one common compact
support strictly inside the annulus. Hence every `W_n` extends smoothly over
the closed annulus.

This exact construction is essential. A generic absolute estimate would not
work: an angular perturbation

\[
 -q^{2h}e^{-4/\delta^2}\delta^{-4}
\]

eventually dominates and reverses the leading
`e^{-4/delta^2} delta^-3` stress as `delta` tends to zero. No fixed positive
`q` prevents that. Exact support and matching flat factors avoid the
obstruction.

At Borel index `n`, add the finitely many actual jets of `W_n` and
`T_n/zeta` through order `n` before choosing the free scale. For a later fixed
derivative order `m`, split into the finite head and the tail `n>=max(m,2)`.
The selected tail is geometric. This proves, on one domain for every fixed
order,

\[
 D^I\left[{\widehat T-T_0\over|T_0|}\right]
 =O_I(q^{2h}).                                       \tag{11}
\]

The constants may depend on the fixed order; no order-uniform seminorm or
convergent raw hierarchy is asserted. Coefficient one must remain in the
finite head, since treating it as diagonal tail would give only `q^h`.

## 4. One explicit perturbation radius and selector

On the leading closed pieces, the analytic regional arguments give positive
finite constants

\[
 f\le F_0\le\bar F,\quad -g_{0\theta}\ge t,
 \quad |g_0|\le\bar G,\quad D_0\ge d,               \tag{12}
\]

together with orientation and same-tilt margins `omega,gamma`. With the
finite Lipschitz constants recorded in the gap reconstruction, the valid
perturbation radius is

\[
 \rho_{gap}=\min\left\{{1\over2},{f\over2},{t\over2},
 {d\over2L_D},{\omega\over2L_H},{m_\phi\over2L_\phi}\right\}>0. \tag{13}
\]

It yields

\[
 F\ge f/2,\quad-g_\theta\ge t/2,\quad D\ge d/2,
\]

\[
 \lambda_{shear}^2\ge{ftd\over4(\bar G+1)^2},
 \qquad c_{shear}^2\ge{d\over4(\bar F+1)(\bar G+1)}. \tag{14}
\]

Only finitely many value, direction, chart, length and slow rows enter the
later integer selector. They are derivative-order independent, and every
constant is known before the selector is chosen. The selection order is
therefore acyclic. The full phase-normal floor remains on the same accepted
tail through the reference-relative estimate; axial transport and the
unprojected pressure source are retained.

## 5. Result, limits and next obligation

The same-band regional gap and weighted relative-stress leaf consumed by Q216
is supported by the cited continuum source chain. No diagnostic-only row or
circular scale choice was found. The critical correction is (5): the valid
input is the factored weighted floor and relative quotient, not a positive
absolute stress infimum.

This audit remains relative to the named regional ODE, exact-root and support
lemmas and does not reprove each one from first principles. It establishes no
uniform actual-pulse size or derivative bound, complete covariance correction,
global residual improvement, unconditional compact force, competitor
exclusion or Navier--Stokes blow-up.

The next justified obligation is Q216's next cited source theorem: audit the
shifted two-edge weighted-kernel estimate for every finite pole. That audit
must retain the physical translation in both directed integrals, prove that
the output keeps the same edge weight and pole, and decide whether its
symbolic exponential constant is genuinely finite without numerical
materialization.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q219_regional_gap_weighted_stress_source_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q219_regional_gap_weighted_stress_source_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q219_regional_gap_weighted_stress_source_audit/evidence.json \
  --output fresh.json
```
