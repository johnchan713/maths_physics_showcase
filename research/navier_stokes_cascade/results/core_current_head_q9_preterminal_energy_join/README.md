# Current-head Q9 preterminal energy join

Parent: `ea485170aed8d1fd9550539cddbeb512f4e9317c`.
Status: `current-head-Q9-preterminal-energy-and-dissipation-joined-conditionally`.

Q8 constructs, conditional on its explicit punctured-endpoint premise, the
exact divergence-free compactly supported field `u` for every `t<1` and
one compact smooth force `f` satisfying

\[
 \partial_tu-\Delta u+(u\cdot\nabla)u+\nabla p=f,
 \qquad \nabla\cdot u=0,\qquad u(0)=0.             \tag{1}
\]

The Q8 representatives retain unit viscosity, all twelve localization terms,
the complete fixed-gauge pressure, actual axial velocity and
`-div(f)` in the pressure equation. This checkpoint proves the exact
preterminal energy consequence for those representatives. It does not remove
Q8's endpoint premise.

Executable tests freeze rational arithmetic, scope and provenance. They are
not a numerical construction or a continuum proof of the inherited endpoint
inequality.

## 1. Exact fixed-time identity

For each fixed `t<1`, `u` and `p` are smooth and compactly supported in
space. Incompressibility and integration by parts give

\[
 \int u\cdot\nabla p=0,\qquad
 \int u\cdot(u\cdot\nabla)u=0,\qquad
 -\int u\cdot\Delta u=\|\nabla u\|_2^2.          \tag{2}
\]

Define

\[
 E(t)=\|u(t)\|_2^2,\qquad D(t)=\|\nabla u(t)\|_2^2.
\]

Taking the `L2` inner product of (1) with `u` yields

\[
 {1\over2}E'(t)+D(t)=\langle f(t),u(t)\rangle.       \tag{3}
\]

The fixed pressure gauge remains essential in Q8's force and pressure
equation, although its contribution to (3) vanishes exactly.

## 2. Zero-safe square-root estimate

Division by `sqrt(E)` is invalid at the zero initial datum. For
`delta>0`, set `Y_delta=sqrt(E+delta^2)`. From (3) and `D>=0`,

\[
 Y_\delta'(t)
 ={E'(t)\over2\sqrt{E(t)+\delta^2}}
 \le \|f(t)\|_2
 {\sqrt{E(t)}\over\sqrt{E(t)+\delta^2}}
 \le \|f(t)\|_2.                                   \tag{4}
\]

With

\[
 F(t)=\int_0^t\|f(s)\|_2\,ds,
\]

integration gives
`sqrt(E(t)+delta^2)<=delta+F(t)`. Letting `delta` decrease to zero
proves

\[
 \sqrt{E(t)}\le F(t).                                \tag{5}
\]

No terminal trace of velocity, pressure or energy is assumed.

## 3. Factor-two energy and dissipation bound

Integrating (3) and using (5), for every `t<1`,

\[
\begin{aligned}
 E(t)+2\int_0^tD(s)\,ds
 &=2\int_0^t\langle f(s),u(s)\rangle\,ds\\
 &\le2\int_0^t\|f(s)\|_2F(s)\,ds
 =F(t)^2.                                             \tag{6}
\end{aligned}
\]

The factor two is exact. Conditional on Q8, `f` is compact and smooth, so
`F(1)<infinity`. Equation (6) yields the separate terminal estimates

\[
 \sup_{t<1}E(t)\le F(1)^2,\qquad
 \int_0^1D(t)\,dt\le {F(1)^2\over2}.               \tag{7}
\]

The dissipation limit follows by monotone convergence of the nonnegative
integrals as `t` increases to one. No value `u(1)`, `p(1)` or `E(1)`
is used.

## 4. Correction retained: endpoint extrema cannot be added

The fixed-time estimate (6) does not imply

\[
 \sup_{t<1}E(t)+2\int_0^1D(t)\,dt\le F(1)^2.        \tag{8}
\]

The energy supremum and total dissipation are different endpoint operations.
The exact scalar negative control remains:

\[
 y(t)=t/\epsilon\quad(0\le t\le\epsilon),\qquad
 y(t)=e^{-(t-\epsilon)}\quad(\epsilon\le t\le1),
 \qquad \epsilon={1\over100}.
\]

For `g=y'+y`, `g` is nonnegative, vanishes after `epsilon`, and

\[
 \int_0^1g={201\over200}.
\]

Using `e^{-s}<=1/(1+s)` gives

\[
 \sup y^2+2\int_0^1y^2
 \ge {37349\over22350}
 >{40401\over40000}
 =\left(\int_0^1g\right)^2.                         \tag{9}
\]

This scalar path refutes (8); it is not a Navier--Stokes solution. It
preserves the fixed-time structure and both separate conclusions in (7).

## 5. Remaining conditional premise and scope

The energy argument does not repair Q8's exact obstruction:

\[
 \boxed{
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p)|<\infty
 }                                                     \tag{10}
\]

for every fixed order on each punctured terminal compact. Thus (6)--(7) are
conditional consequences for the exact Q8 field and force.

Q9 proves neither a terminal field trace, future velocity or pressure,
fixed-gauge pressure recovery for competitors, expanding-ball comparison,
competitor exclusion nor blow-up.

The next justified obligation is fixed-gauge pressure recovery for this exact
finite-energy field and compact force, retaining `-div(f)`, followed by
expanding-ball comparison with every admissible smooth same-force
finite-energy competitor.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q9_preterminal_energy_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q9_preterminal_energy_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q9_preterminal_energy_join/evidence.json \
  --output fresh.json
```
