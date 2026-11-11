# Q181: physical diagonal, endpoint and localized-force replay

Parent: Q180, `2e4249047efadff763bdaa03b72a2f2674f66b9d`.
The latest remote head, Q180, Q162 and the later Q166--Q168 source audits,
their frozen evidence and current CI were read first. Historical evidence,
the inherited physical field, gauges, support data, cutoff sequence, accepted
scales and `C^-20` width remain unchanged.

Q180 closes the coupled target/source/tail induction relative to its named
continuum maps and leaves the complete-potential diagonal, endpoint and force
replay next. Q162 contains that same-object replay. Q166--Q168 later sharpen
three points which are essential here: residual flatness and endpoint field
bounds are separate obligations; the localized residual has exactly twelve
nonzero groups; and the zero-datum energy estimate must not divide by
`sqrt(E)` before regularization.

The result below is conditional on the named continuum estimates. Exact
arithmetic, tests and evidence hashes certify interfaces and bookkeeping;
they do not prove those estimates.

## 1. Same physical object

Q180 and Q162 contain the byte-identical Q161 evidence digest. Q167 consumes
the frozen Q162 record, and Q168 consumes the frozen Q167 record. Thus the
field, potential and pressure gauges, support holes and diagonal cutoff
sequence are the ones fed by Q180's simultaneous physical induction. No
canonical-family transfer or parallel reconstruction is used.

At finite stage `j` retain

\[
 \Delta u_j=\nabla\times A_j+b_j e_\theta,
 \qquad \pi_j=\hbox{fixed-gauge pressure}.             \tag{1}
\]

The diagonal terms are

\[
 V_j=\nabla\times(\chi(a_jq)A_j)+\chi(a_jq)b_je_\theta,
 \qquad P_j=\chi(a_jq)\pi_j.                          \tag{2}
\]

Each velocity term is divergence free. The scale `a_j` is selected once,
after the finitely many constants requested at that index are known. One
diagonal handles the countable list of fixed-jet requests. Local finiteness,
not convergence of the raw uncut series, defines the resulting field.

## 2. Complete six-group residual and strict margins

Fix a Cartesian jet `m` and target power `N` before choosing a finite stage.
With

\[
 h={1\over200},\quad \kappa={1\over100000},\quad
 \delta={1\over2}-3\kappa,
\]

the comparison-stage condition is

\[
 {h\delta(J_0+1)\over2}-L_{m+2}-H_{m+2}\ge N+1.      \tag{3}
\]

The replay takes `J=J0+1`. The limiting regular and nonlinear-cross reserve is

\[
 {h\delta\over2}={49997\over40000000}>0.             \tag{4}
\]

The compact `A0` and `A1` logarithms each use their own reserve

\[
 {h\kappa\over2}={1\over40000000}>0.                 \tag{5}
\]

These three accounts are disjoint and are not added to manufacture a stronger
common exponent. The exact comparison residual contains all six groups

\[
 e_t-\Delta e+\nabla\pi
 +(U\cdot\nabla)e+(e\cdot\nabla)U+(e\cdot\nabla)e.   \tag{6}
\]

Thus time, unit viscosity, complete pressure, both nonlinear orientations and
the square beat the requested power for every fixed jet. All six covariance
entries, actual axial velocity, `A0/A1`, Gaussian and exterior tails,
`Pi1/Pi2`, the whole `partial_Z(J_z+c_rho P)` and every cutoff curl/pressure
derivative remain. The outcome is conditional flatness of every fixed local
Cartesian residual jet. It is not a stage-uniform estimate and does not use
an uncontrolled `J(q)` substitution.

## 3. Corner transfer and punctured endpoint are separate

The chart Jacobian has lower bound `99/100`. For arbitrary `N`, the exact
corner relation `q<=4 rho` transfers the complete residual estimate to
Cartesian order `rho^(N+1)`. This proves fixed-jet residual flatness at the
closed corner without a positive support-distance or edge-weight floor.

Residual flatness does not bound the potential, direct swirl or pressure at a
punctured endpoint. Q154/Q167 supply the separate estimate

\[
 \sup_{E\times[1-\delta_E,1)}
 |D_x^\alpha\partial_t^j(A,B,p)|<\infty              \tag{7}
\]

for every punctured terminal compact `E` and fixed orders. The full enlarged
band window and fixed pressure gauge remain. The pressure time jets retain

\[
 S_j=\sum_{k=0}^j {j\choose k}A_kA_{j-k},
 \qquad S_2-2A_2=2A_1^2>0.                           \tag{8}
\]

The successive exterior Taylor-coefficient ratio grows linearly, so the time
Taylor series has zero radius at fixed nonzero radius. Only smooth one-sided
jet gluing is licensed. The terminal origin is not included or regularized.

## 4. Exact localization and compact force

Let

\[
 u=\nabla\times(cA)+cB=cu_0+d,\qquad
 d=\nabla c\times A,\qquad p=cp_0.                   \tag{9}
\]

This is divergence free and has zero initial velocity. Its unit-viscosity
residual has exactly twelve nonzero groups:

1. `c R_1(u0,p0)`;
2. `c_t u0`;
3. `d_t`;
4. `-(Delta c)u0`;
5. `-2 sum_i c_i partial_i u0`;
6. `-Delta d`;
7. `p0 grad c`;
8. `c(c-1)(u0 dot grad)u0`;
9. `c(u0 dot grad c)u0`;
10. `c(u0 dot grad)d`;
11. `c(d dot grad)u0`;
12. `(d dot grad)d`.

The apparent thirteenth group vanishes exactly because

\[
 d\cdot\nabla c=(\nabla c\times A)\cdot\nabla c=0.   \tag{10}
\]

The force divergence remains in the pressure equation. For a force jet of
order `m`, the sharp direct requests are `A:m+2`, `B:m+1`, `p:m`; Cauchy
convergence uses one more derivative. Q167 supplies one further unused field
derivative in every channel.

The terminal force jets are extended by one Borel schedule whose dyadic scale
at index `j` is chosen after the finitely many actual norms at that index are
known. Every fixed differentiated tail is geometric. This constructs one
compact smooth force after time one. It does not construct a future velocity
or pressure and does not assume convergence of the raw Taylor series.

## 5. Zero-safe energy estimate and scope

For `t<1`,

\[
 {1\over2}E'(t)+D(t)=\langle f(t),u(t)\rangle,
 \quad E=\|u\|_2^2,\quad D=\|\nabla u\|_2^2.         \tag{11}
\]

Since `E(0)=0`, direct division by `sqrt(E)` is invalid. With
`Y_epsilon=sqrt(E+epsilon^2)`, one obtains `Y_epsilon'<=||f||_2`, then sends
`epsilon` to zero. If `F(t)=int_0^t ||f(s)||_2 ds`, this yields

\[
 E(t)+2\int_0^tD(s)\,ds\le F(t)^2.                  \tag{12}
\]

The energy supremum and integrated-dissipation bounds are separate
consequences; their extrema are not added under one copy of `F(1)^2`.

This closes P6 on Q180's physical object relative to the named continuum
providers. It establishes conditional local residual flatness, punctured
endpoint jets, exact localization, a compact smooth force and preterminal
energy bounds. It does not regularize the terminal origin, prove a global
residual improvement, construct future velocity/pressure, or prove comparison,
competitor exclusion or Navier--Stokes blow-up.

Next: replay fixed-gauge pressure recovery, expanding-ball comparison and the
protected-core lower bound on this exact localized field. The required limit
order is fixed `T<1`, then ball radius to infinity, then `tau` to zero.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q181_diagonal_endpoint_force_replay/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q181_diagonal_endpoint_force_replay/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q181_diagonal_endpoint_force_replay/evidence.json \
  --output q181_replay.json
```
