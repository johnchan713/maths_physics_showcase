# Q168: exact localization, compact force and energy frontier audit

Parent: `e273c3b389bd708b8b61287488b130d0e0d5c386` (Q167).
Status: `Q168-current-head-localization-force-energy-audited-conditionally`.

This checkpoint audits Q8's localization and force construction on the exact
Q164--Q167 object.  It reads Q8's original algebra, Q9's energy argument,
Q155's physical force/energy join and Q162's same-object binding.  The point
is not to repeat Q8: Q8 left a punctured endpoint estimate open, whereas Q154
later proves that estimate and Q162 freezes the same physical field, cutoff
sequence, supports and pressure gauge.  Historical evidence is unchanged.

The conclusion is conditional on the named continuum providers.  Executable
ledgers check formulas, finite derivative requests and provenance; they are
not proofs of those continuum providers.

## 1. Same-object chain and the corrected premise

The frozen evidence gives the exact chain

\[
 \mathrm{Q167}\to\mathrm{Q162}\to\mathrm{Q155}\to\mathrm{Q8},
 \qquad \mathrm{Q162}\to\mathrm{Q154},
\]

while Q167 and Q8 contain the same Q7 evidence digest.  Q9 also contains that
same Q8 digest.  Q162 explicitly binds the physical field, diagonal cutoff
sequence, supports and fixed pressure gauge across Q154 and Q155.  Therefore
Q154 supplies, for the object audited at Q167, the estimate that Q8 displayed
as its remaining premise:

\[
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p_0)|<\infty
 \tag{1}
\]

for every fixed order on each punctured terminal compact \(E\).  This closes
the premise conditionally.  It does not include the terminal origin and does
not make the exterior time Taylor series convergent.

## 2. Divergence-free localization and zero datum

Write

\[
 u_0=\nabla\times A+B,\qquad B=b(r,z,t)e_\theta,
 \qquad p_0=p_{\rm complete}.
\]

For the Cartesian-smooth axisymmetric cutoff
\(c=\chi_r(r^2)\chi_z(z)\chi_t(t)\), set

\[
 u=\nabla\times(cA)+cB=cu_0+d,
 \qquad d=\nabla c\times A,
 \qquad p=cp_0.                                      \tag{2}
\]

The curl term is divergence free, and
\(\operatorname{div}(cB)=r^{-1}\partial_\theta(cb)=0\).  The time cutoff is
identically zero near the initial time, so both \(c\) and \(\nabla c\) vanish
there and \(u(0)=0\).  On the protected late core \(c=1\), hence \(d=0\) and
the Q167 field is unchanged.

## 3. Exact twelve-group residual

For

\[
 \mathcal R_1(v,p)=v_t-\Delta v+(v\cdot\nabla)v+\nabla p,
\]

substitution of (2) gives

\[
\begin{aligned}
\mathcal R_1(u,p)={}&c\mathcal R_1(u_0,p_0)+c_tu_0+d_t
 -(\Delta c)u_0-2\sum_i c_i\partial_i u_0-\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
 +c(u_0\cdot\nabla c)u_0+c(u_0\cdot\nabla)d\\
&+c(d\cdot\nabla)u_0+(d\cdot\nabla)d.              \tag{3}
\end{aligned}
\]

These are twelve nonzero groups.  To see that the list is exact, expand

\[
(cu_0+d)\cdot\nabla(cu_0+d).
\]

It produces one additional-looking term
\((d\cdot\nabla c)u_0\), but

\[
 d\cdot\nabla c=(\nabla c\times A)\cdot\nabla c=0. \tag{4}
\]

Equation (3) retains unit viscosity, all three viscous cutoff groups, the
complete fixed-gauge pressure term \(p_0\nabla c\), both nonlinear cross
orientations and \((d\cdot\nabla)d\).  Taking divergence of the forced
equation retains the force divergence:

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f_- . \tag{5}
\]

No divergence-free-force assumption or support-preserving Leray projection
is inserted.

## 4. Derivative requests

For a mixed derivative of total order \(m\), (3) directly requests

\[
 A:m+2,\qquad B:m+1,\qquad p_0:m.                    \tag{6}
\]

The largest \(A\) request comes from \(\Delta d\) or the differentiated
gradient of \(u_0\); the largest \(B\) request comes from the gradient of
\(u_0\); the pressure cutoff contains only derivatives placed on \(p_0\).
Uniform Cauchy convergence uses one further time derivative:

\[
 A:m+3,\qquad B:m+2,\qquad p_0:m+1,
 \qquad \partial^{m+1}\mathcal R_1.                 \tag{7}
\]

Q167's endpoint ledger supplies \(A:m+4\), \(B:m+3\), \(p:m+2\), leaving
one unused field derivative in every channel.  For the frozen sample
\(m=5\), the direct, Cauchy and supplied requests are respectively
\((7,6,5)\), \((8,7,6)\) and \((9,8,7)\), with residual jet six.

## 5. Terminal force jets and one compact extension

Near the closed Cartesian corner \(c=1\) and \(d=0\), so the preterminal
force \(f_-=\mathcal R_1(u,p)\) is the complete Q167 residual and all of its
fixed spatial jets are flat there.  On punctured terminal compacts, (1) and
the one-higher requests in (7) give uniform Cauchy limits.  A finite compact
cover therefore yields, conditionally,

\[
 F_j=\lim_{t\uparrow1}\partial_t^j f_-
 \in C_c^\infty(K),
 \qquad D_x^\alpha F_j(0)=0.                         \tag{8}
\]

For \(\sigma=t-1\ge0\), Q8's extension is

\[
 f_+(x,1+\sigma)=\sum_{j\ge0}
 \chi_0(b_j\sigma){\sigma^j\over j!}F_j(x).         \tag{9}
\]

At index \(j\), the increasing dyadic \(b_j\) is selected once after the
finitely many actual terminal norms with \(|a|+m\le\lfloor j/2\rfloor\)
are known.  On the full support of a differentiated summand,

\[
 \|D_x^a\partial_\sigma^m T_j\|_\infty
 \le C_{j,a,m}b_j^{m-j}.                             \tag{10}
\]

The admitted diagonal has \(m<j\); choosing each bound at most \(2^{-j}\)
gives a geometric tail for every fixed derivative.  Thus one obtains
conditionally

\[
 f\in C_c^\infty(\mathbb R^3\times(0,\infty)),
 \qquad \operatorname{supp}f\subset K\times[1-\tau_0,3/2]. \tag{11}
\]

This uses actual finite norms index by index.  It assumes neither raw Taylor
convergence nor a uniform-in-\(j\) norm.  It extends only the force; no future
velocity or pressure is constructed.

## 6. Zero-safe energy estimate

For each \(t<1\), the localized fields are smooth and spatially compact.
Incompressibility and integration by parts give

\[
 {1\over2}E'(t)+D(t)=\langle f(t),u(t)\rangle,
 \quad E=\|u\|_2^2,\quad D=\|\nabla u\|_2^2.         \tag{12}
\]

Because \(E(0)=0\), dividing directly by \(\sqrt E\) would be invalid.  For
\(\varepsilon>0\), put \(Y_\varepsilon=\sqrt{E+\varepsilon^2}\).  Then

\[
 Y_\varepsilon'
 \le \|f\|_2{\sqrt E\over\sqrt{E+\varepsilon^2}}
 \le\|f\|_2.
\]

With \(F(t)=\int_0^t\|f(s)\|_2\,ds\), integration and
\(\varepsilon\downarrow0\) give \(\sqrt{E(t)}\le F(t)\).  Returning to
(12) yields the exact factor-two bound

\[
 E(t)+2\int_0^tD(s)\,ds
 \le 2\int_0^t\|f(s)\|_2F(s)\,ds=F(t)^2.            \tag{13}
\]

Consequently

\[
 \sup_{t<1}E(t)\le F(1)^2,
 \qquad \int_0^1D(t)\,dt\le {F(1)^2\over2}.         \tag{14}
\]

The dissipation limit uses monotone convergence.  The two endpoint extrema
in (14) are separate consequences; their left sides cannot be added while
keeping the same right side.  No terminal value of \(u\), \(p\), or \(E\)
is used.

## 7. Scope and next obligation

Q168 closes the current localization/force/energy frontier conditionally.
It does not regularize the singular terminal origin, construct a solution
after time one, or prove comparison, competitor exclusion, residual
improvement or blow-up.

The next justified obligation is to audit Q156's fixed-gauge pressure
recovery and expanding-ball comparison on this Q164--Q168 object.  That
audit must retain \(-\nabla\cdot f\) in both individual pressure equations,
cancel the common force exactly once in the difference equation, and take
the limits in the order: fix \(T<1\), send the ball radius to infinity, then
approach the terminal time.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q168_localization_force_frontier_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q168_localization_force_frontier_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q168_localization_force_frontier_audit/evidence.json \
  --output q168-fresh.json
```
