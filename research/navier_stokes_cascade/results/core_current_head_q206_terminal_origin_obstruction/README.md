# Q206: terminal-origin velocity and pressure obstruction

Parent: Q205, `b3a47a3e059b24dcbad8e9f43135d673f8ce4777`.
The live branch, Q205, its active CI, the Q203 endpoint/localization replay,
the Q182 protected-core comparison, Q170's same-field core path, the leading
pressure identification and the Borel background sum were read first.
Historical evidence, the selected profile, physical family, fixed pressure
gauge, supports, diagonal cutoffs, localization and unchanged `C^-20` width
are preserved.

Q205 asks whether the locally finite field has one-sided velocity and pressure
jets at the terminal origin. The answer is a rigorous obstruction relative to
the same named continuum providers. The value jet already fails: the protected
azimuthal velocity diverges like `tau^-A`. The pressure has a stronger
gauge-independent obstruction: its difference along two protected paths
diverges like `tau^-2A`, so no time-dependent spatial pressure constant can
repair continuity at the origin. Punctured one-sided jets and the smooth force
trace remain valid and are not promoted to field regularity at the corner.

The lower bounds below are analytical consequences of the frozen provider
estimates. Exact rational checks and evidence hashes certify the composition;
they do not prove those continuum providers or an unconditional blow-up
theorem.

## 1. Same protected core and velocity obstruction

Keep

\[
 h={1\over200},\qquad A={1\over2}+h={101\over200},
 \qquad X_a={4\over101},\qquad X_1={1\over101}.       \tag{1}
\]

For `tau=1-t>0`, put

\[
 x_{1,\tau}=(\sqrt{2X_1\tau},0,0).                    \tag{2}
\]

Then `q=tau`, `eta=0`, `X=X1`. The leading normalized swirl is

\[
 E_0(X_1,0)=e_0>0.                                   \tag{3}
\]

Every sourced/signed diagonal potential, direct swirl, pressure difference
and synchronized annular initialization has closed support beginning at
`Xa`, so it vanishes in a neighborhood of (2). The positive-order Borel
background is not zero, but on the same profile it satisfies

\[
 \tau^A u_\theta(x_{1,\tau},1-\tau)
 =e_0+R_u(\tau),\qquad |R_u(\tau)|\le C_u\tau^{2h}.   \tag{4}
\]

After the already selected localization plateau restrictions, (4) gives

\[
 |u_\theta(x_{1,\tau},1-\tau)|
 \ge {e_0\over2}\tau^{-101/200}.                     \tag{5}
\]

The points `(x1_tau,1-tau)` tend to `(0,1)`. A continuous one-sided terminal
velocity would be bounded on a compact corner neighborhood, contradicting
(5). Hence even the velocity value jet at the origin does not exist. No
higher-jet question can repair this failure.

## 2. Leading pressure gives a two-path obstruction

The retained fixed gauge is

\[
 \Pi_0(X,\eta)=-\int_X^\infty {E_0(x,\eta)^2\over2x}\,dx. \tag{6}
\]

Continuity and (3) give an `X2` with

\[
 X_1<X_2<X_a,\qquad |E_0(X,0)|\ge e_0/2
 \quad(X_1\le X\le X_2).                             \tag{7}
\]

Therefore

\[
 d_0:=\Pi_0(X_2,0)-\Pi_0(X_1,0)
 ={1\over2}\int_{X_1}^{X_2}{E_0(X,0)^2\over X}\,dX
 \ge {e_0^2\over8}\log{X_2\over X_1}>0.             \tag{8}
\]

Put `x2_tau=(sqrt(2 X2 tau),0,0)`. Both paths remain strictly inside the
same annular support hole. The rule-defined background pressure, including
all retained positive orders and cutoff derivatives, has the fixed-profile
normalized remainder

\[
 \tau^{2A}\{p(x_{2,\tau},1-\tau)-p(x_{1,\tau},1-\tau)\}
 =d_0+R_p(\tau),\qquad |R_p(\tau)|\le C_p\tau^{2h}.  \tag{9}
\]

Every generated annular pressure is exactly zero on these paths, and the
localization is one there. Choosing the same coefficient-after-constant
threshold as in the Borel/core argument yields

\[
 |p(x_{2,\tau},1-\tau)-p(x_{1,\tau},1-\tau)|
 \ge {d_0\over2}\tau^{-101/100}.                     \tag{10}
\]

Equation (10) is invariant under `p -> p+c(t)`. If any such gauge had a
continuous one-sided terminal pressure at the origin, both path values would
approach the same finite number and their difference would tend to zero,
contradicting (10). Thus the pressure value jet also fails, independently of
the fixed representative's absolute divergence.

The normalized remainder in (9) is stagewise and conditional: each fixed
positive background order begins at `tau^(2h)`, the Borel selector absorbs
its finite constants, and Q205's complete-potential cutoff removes every
generated annular pressure in the protected hole. No uniform raw hierarchy
constant is used.

## 3. Endpoint and force statements remain correctly separated

On every punctured terminal compact `|x|>=d>0`, the closed atlas sees only
finitely many bands and stages and supplies compatible one-sided potential,
direct-swirl, velocity and full-pressure jets. The smooth localized force
trace follows from the arbitrary-power residual bound at the Cartesian
corner plus those punctured field jets. Neither fact controls `u` or `p` at
`x=0`.

This distinction is essential. A smooth residual can result from cancellation
among singular time, viscosity, pressure and nonlinear terms. Smooth force
gluing does not imply a smooth preterminal field. The protected lower bounds
(5) and (10) therefore correct any reading of Q203's endpoint/force replay as
a regular terminal field or a future velocity/pressure construction.

## 4. Scope and next obligation

Q206 establishes, relative to the frozen same-object providers, that regular
one-sided velocity and pressure jets at the terminal origin are impossible.
It does not turn the conditional construction into an unconditional pulse,
residual, admissible-force, competitor-exclusion or Navier--Stokes blow-up
theorem. It constructs no field for `t>=1`.

The next obligation returns upstream to the exact promotion frontier rather
than trying to regularize an intentionally singular corner: reconstruct the
fixed-stage complete residual constant `C_(j,m,N)` directly from the selected
profile, roots and compact/slow inverse maps, retaining pressure, viscosity,
all cutoffs, both cross orientations and the square. A finite symbolic request
or passing test is not that continuum estimate.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q206_terminal_origin_obstruction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q206_terminal_origin_obstruction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q206_terminal_origin_obstruction/evidence.json \
  --output fresh.json
```
