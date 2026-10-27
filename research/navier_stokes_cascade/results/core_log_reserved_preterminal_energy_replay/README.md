# Log-reserved same-field preterminal energy replay

Parent: `3bb6c5d1beba0f637d66695e87b75e618a2ee120`.

The parent constructs, relative to the frozen continuum contracts, the
same divergence-free compactly supported field `u` for every `t<1` and one
smooth compact force `f` satisfying

\[
 \partial_tu-\Delta u+(u\cdot\nabla)u+\nabla p=f,
 \qquad \nabla\cdot u=0,
 \qquad u(0)=0.                                         \tag{1}
\]

It retains the complete fixed-gauge pressure, all twelve localization terms
and unit viscosity.  This checkpoint replays the energy consequence for
those exact representatives.  Historical evidence remains unchanged.

The implication is conditional on the parent chain's continuum estimates.
The executable controls verify rational constants, dependency hashes and
scope; passing them is not an independent proof that the infinite field
exists with those estimates.

## 1. Exact identity before the terminal time

For each fixed `t<1`, the localized `u` and `p` are smooth and compactly
supported in space.  Hence no boundary-at-infinity limit is required:

\[
 \int u\cdot\nabla p=0,
 \qquad
 \int u\cdot(u\cdot\nabla)u=0,
 \qquad
 -\int u\cdot\Delta u=\|\nabla u\|_2^2.                 \tag{2}
\]

The first two equalities use `div u=0`; the last uses unit viscosity and
spatial integration by parts.  Define

\[
 E(t)=\|u(t)\|_2^2,
 \qquad D(t)=\|\nabla u(t)\|_2^2.
\]

Taking the `L2` inner product of (1) with `u` gives the exact identity

\[
 {1\over2}E'(t)+D(t)=\langle f(t),u(t)\rangle.           \tag{3}
\]

## 2. The zero-data square-root step

The formal division by `sqrt(E)` is singular at the zero initial datum.
For `delta>0`, set `Y_delta=sqrt(E+delta^2)`.  From (3) and `D>=0`,

\[
 Y_\delta'(t)
 ={E'(t)\over2\sqrt{E(t)+\delta^2}}
 \le \|f(t)\|_2{\sqrt{E(t)}\over\sqrt{E(t)+\delta^2}}
 \le \|f(t)\|_2.                                      \tag{4}
\]

With

\[
 F(t)=\int_0^t\|f(s)\|_2\,ds,
\]

integration of (4), followed by `delta` decreasing to zero, yields

\[
 \sqrt{E(t)}\le F(t).                                  \tag{5}
\]

This step uses only the actual zero datum and the preterminal smooth field;
it does not assume a terminal trace of `u` or `p`.

## 3. Fixed-time bound and endpoint consequences

Integrating (3) and using (5) gives, for every `t<1`,

\[
\begin{aligned}
 E(t)+2\int_0^tD(s)\,ds
 &=2\int_0^t\langle f(s),u(s)\rangle\,ds\\
 &\le2\int_0^t\|f(s)\|_2F(s)\,ds
 =F(t)^2.                                               \tag{6}
\end{aligned}
\]

The factor two on dissipation is essential.  Because the force is smooth
and compactly supported, `F(1)<infinity`.  Equation (6) implies the two
separate bounds

\[
 \sup_{t<1}E(t)\le F(1)^2,
 \qquad
 \int_0^1D(t)\,dt\le{F(1)^2\over2}.                    \tag{7}
\]

The second follows by monotone convergence of the nonnegative integrals
`integral_0^t D` as `t` increases to one.  No terminal velocity, pressure or
energy equality is used.

## 4. Correction: the endpoint bounds cannot be added

Equation (6) does not imply

\[
 \sup_{t<1}E(t)+2\int_0^1D(t)\,dt\le F(1)^2.            \tag{8}
\]

The two left-side extrema need not occur at the same time.  An exact scalar
diagnostic makes this failure quantitative.  Let `epsilon=1/100` and

\[
 y(t)=t/\epsilon\quad(0\le t\le\epsilon),
 \qquad
 y(t)=e^{-(t-\epsilon)}\quad(\epsilon\le t\le1).
\]

For `g=y'+y`, the forcing is nonnegative, vanishes after `epsilon`, and

\[
 \int_0^1g={201\over200}.
\]

Using `e^{-s}<=1/(1+s)` gives

\[
 \sup y^2+2\int_0^1y^2
 \ge {37349\over22350}
 >{40401\over40000}
 =\left(\int_0^1g\right)^2.                             \tag{9}
\]

This scalar path is a negative control for the inference in (8), not a
Navier--Stokes solution.  It preserves (6) and both separate conclusions in
(7).

## 5. Scope and next obligation

This checkpoint establishes the **conditional same-field preterminal energy
and total-dissipation bounds** (6)--(7).  It does not regularize the terminal
origin, construct future velocity or pressure, recover pressure for a
competitor, or prove comparison, competitor exclusion, or blow-up.

The next obligation is fixed-gauge pressure recovery with
`-div f` retained in the pressure Poisson equation, followed by the
expanding-ball comparison against every admissible smooth finite-energy
solution with the same zero datum and the same compact force.

Reproduce with

```sh
python3 test_audit.py -v
python3 audit.py --verify-record evidence.json --output fresh.json
```
