# Q218: energy-class entry to fixed-gauge comparison

Parent: Q217, `fe03f338f4c59c3ece385a4f7cd359b050977d30`.
The live branch, Q217 and its active CI, Q168--Q170, Q206 and their frozen
evidence were read first. Historical evidence, the source-audited physical
field, compact smooth force, complete pressure gauge, localization, support
hole, diagonal scales and unchanged `C^-20` width remain fixed.

Q217 identifies the competitor's uniform `L2` bound on each fixed
preterminal interval as the first new premise in the downstream comparison
chain. This checkpoint corrects the competitor statement. The premise is
automatic for a smooth solution in the standard energy class and has a direct
quantitative bound from the energy inequality. It is not automatic for a bare
locally smooth solution; the affine harmonic-pressure mode remains the exact
obstruction.

All conclusions remain conditional on Q216's upstream continuum provider
theorems. Tests and evidence hashes verify formulas, scope and provenance;
they are not proofs of those providers.

## 1. The correct comparison class

Fix `T<1`. Let `(v,P)` solve Navier--Stokes on `R^3 x [0,T]` with the Q217
force `f` and zero initial velocity. The energy comparison class consists of
solutions which are locally smooth and satisfy

\[
 v\in L^\infty(0,T;L^2(\mathbb R^3)),
 \qquad \nabla v\in L^2(0,T;L^2(\mathbb R^3)),       \tag{1}
\]

together with the energy inequality. No spatial decay condition is imposed on
`P`, and no terminal smoothness at `t=1` is assumed.

This is a real restriction on behavior at spatial infinity. Calling a
solution globally or locally smooth does not by itself imply (1). The
comparison theorem must therefore say “smooth energy-class competitor,” not
merely “smooth competitor.”

## 2. Quantitative entry from the energy inequality

Set

\[
 E_v(t)=\|v(t)\|_2^2,
 \qquad D_v(t)=\|\nabla v(t)\|_2^2.
\]

For an energy solution with zero datum,

\[
 E_v(t)+2\int_0^tD_v(s)\,ds
 \le 2\int_0^t\langle f(s),v(s)\rangle\,ds.          \tag{2}
\]

Young's inequality gives

\[
 2\langle f,v\rangle\le\|f\|_2^2+E_v,              \tag{3}
\]

so Gronwall yields the fixed-horizon estimate

\[
 \boxed{
 \sup_{0\le t\le T}\|v(t)\|_2^2
 \le e^T\int_0^T\|f(s)\|_2^2\,ds <\infty .}        \tag{4}
\]

The last quantity is finite because Q217's force is smooth and compactly
supported in space-time. Equation (4) is more than the Q169 comparison needs:
it supplies its uniform `L2` premise on every fixed `T<1`. No estimate uniform
as `T` approaches one is required.

Equation (4) does not derive energy-class membership from local smoothness.
The energy inequality and global spatial integrability in (1) are the stronger
boundary conditions that exclude harmonic acceleration modes.

## 3. Pressure recovery now enters without an extra hypothesis

Let `u,p` be the Q217 localized field and pressure and put

\[
 w=v-u,\qquad \pi=P-p,
 \qquad g=w\otimes w+w\otimes u+u\otimes w.          \tag{5}
\]

The individual pressure equations retain `-div f`; subtraction cancels the
common force once. The Q217 energy estimate for `u` and (4) for `v` imply

\[
 w\in L^\infty(0,T;L^2),
 \qquad g\in L^\infty(0,T;L^1).                     \tag{6}
\]

Thus the Riesz pressure of `g` lies in `H^-2`. The momentum equation places
the difference between the actual and Riesz pressure gradients in `H^-3`.
A harmonic `H^-3` distribution has zero gradient, so no competitor-pressure
decay assumption is needed. The complete force divergence remains in the two
individual pressure equations and cancels exactly once in their difference.

Q169's eighth-power expanding-ball cutoff then gives, at fixed `T`,

\[
 {1\over2}E_R'+{1\over2}D_R^2
 \le\|\nabla u(t)\|_\infty E_R+C_TR^{-2}.            \tag{7}
\]

Sending `R` to infinity yields `v=u` on `[0,T]`.

## 4. The obstruction outside the energy class remains exact

The field

\[
 v(x,t)=t e_1,\qquad P(x,t)=-x_1,
 \qquad f=0                                          \tag{8}
\]

is smooth and divergence free with zero initial velocity. Its time derivative
is `e1`, its pressure gradient is `-e1`, and its viscosity and convection are
zero, so the complete residual vanishes. Nevertheless `v(t)` is not in
`L2(R^3)` for any positive time.

Therefore neither local smoothness nor zero data plus compact smooth forcing
implies (1). Poisson's equation alone cannot remove the nonzero harmonic
pressure gradient in (8). Q217's finite-energy premise was not a technical
artifact; it was an omitted class boundary.

## 5. Protected limit and corrected conditional theorem

On the protected path, Q217 retains

\[
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-101/200}                    \tag{9}
\]

for sufficiently small positive `tau`. For each fixed `tau`, choose
`T_tau=1-tau/2`, apply (4), complete pressure recovery and send the comparison
radius to infinity. Only then send `tau` to zero. A competitor which is smooth
through the terminal point and belongs to the energy class on fixed
preterminal intervals would equal `u` along the path and contradict (9).

Hence, relative to the named upstream continuum providers, there is no
smooth-through-terminal same-force zero-data competitor in the standard
energy class. This is the exact corrected comparison statement. It does not
exclude every locally smooth solution, construct velocity or pressure for
`t>=1`, discharge the upstream providers, produce uniform actual-pulse
constants, or prove an unconditional force, residual or Navier--Stokes
blow-up theorem.

The next justified obligation returns to the first upstream continuum leaf
still cited by Q216: audit the same-band regional gap infima and weighted
relative-stress estimate on the finalized profile. That audit must prove the
continuum regional inequalities and flat-edge factorizations, rather than
promoting interval fixtures or sampled fields to a proof.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q218_energy_class_comparison_entry/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q218_energy_class_comparison_entry/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q218_energy_class_comparison_entry/evidence.json \
  --output fresh.json
```
