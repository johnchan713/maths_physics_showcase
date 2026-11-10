# Q155: physical localization, compact force and preterminal energy join

Parent: Q154, `6a3df11220e95520151ea6e09e7be2066e3d6bd8`.
The live branch, Q154 evidence and CI status, Q8's localization ledger, the
certified complete-potential force theorem, force foundation, current
localized-force record and scale-compatible corner replay were read first.
Historical evidence, physical field, roots, fixed pressure gauge, supports,
accepted bands, diagonal cutoffs and `C^-20` width are unchanged.

Q8's exact localization and force theorem was conditional on its punctured
endpoint bound. Q154 supplies that bound for the later Q149--Q153 physical
field, and Q153 supplies the all-fixed-jet Cartesian-flat complete residual
at the terminal corner. This checkpoint binds those inputs to the same
representatives. The conclusion is a conditional compact smooth force and
preterminal energy estimate, not a future velocity or comparison theorem.

## 1. Divergence-free localization of the same field

Keep Q154's fixed-gauge decomposition

\[
 u_0=\nabla\times A+B,qquad B=b(r,z,t)e_\theta,
 \qquad p_0=p_{\rm complete}.                           \tag{1}
\]

`A` contains the background, all complete-potential diagonal wave and mean
representatives and every cutoff-curl term. `B` is the separately
divergence-safe axisymmetric direct swirl. `p_0` contains both exceptional
background tails, every covariance pressure, axial velocity and moving
pressure-moment contribution in the unchanged gauge.

Choose one Cartesian-smooth axisymmetric cutoff

\[
 c=\chi_r(r^2)\chi_z(z)\chi_t(t)                       \tag{2}
\]

whose support and a neighborhood lie inside Q154's endpoint cylinder. It is
one on the protected late core and zero on an initial-time neighborhood.
Define

\[
 u=\nabla\times(cA)+cB=cu_0+d,qquad
 d=\nabla c\times A,qquad p=cp_0.                     \tag{3}
\]

Both velocity summands are exactly divergence free: the first is a curl and
`div(cB)=r^-1 partial_theta(cb)=0`. Moreover

\[
 d\cdot\nabla c=(\nabla c\times A)\cdot\nabla c=0.     \tag{4}
\]

The time cutoff gives `u(0)=0`. Every derivative of `c-1` vanishes on the
protected core, so Q153's conditional inner growth and the original lower
bound are unchanged there. Directly cutting `u_0`, using a theta-dependent
cutoff on `B`, or reconstructing a potential from velocity would break this
same-field identity.

## 2. Exact twelve-term unit-viscosity residual

For

\[
 \mathcal R_1(v,p)=v_t-\Delta v+(v\cdot\nabla)v+\nabla p,
\]

direct expansion of (3) gives

\[
\begin{aligned}
\mathcal R_1(u,p)={}&c\mathcal R_1(u_0,p_0)+c_tu_0+d_t
 -(\Delta c)u_0-2\sum_i c_i\partial_i u_0-\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
 +c(u_0\cdot\nabla c)u_0+c(u_0\cdot\nabla)d\\
&+c(d\cdot\nabla)u_0+(d\cdot\nabla)d .                \tag{5}
\end{aligned}
\]

These are twelve distinct nonzero groups: base residual; two time terms;
three viscous cutoff terms; full pressure cutoff; quadratic rescaling;
velocity-cutoff transport; both old/curl cross orientations; and the curl
square. The apparent thirteenth group `(d dot grad c)u0` is exactly zero by
(4). The viscosity is one throughout. The pressure term is the complete
`p_0 grad c`; no pressure component or background tail is removed.

For `f_-=R_1(u,p)`, the forced pressure identity is

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f_-. \tag{6}
\]

The force need not be divergence free. Dropping its divergence, using an
unforced pressure formula or applying a fictitious compact-support-preserving
Leray projection would change the equation.

## 3. Correct finite derivative request

After applying a mixed derivative of total order `m` to (5), the sharp
direct field request is

\[
                         A:m+2,\qquad B:m+1,\qquad p_0:m. \tag{7}
\]

`Delta d` and differentiated gradients of `u_0` give the first two maxima;
`p_0 grad c` gives the pressure maximum. Uniform terminal Cauchy convergence
reserves one further time derivative:

\[
 A:m+3,qquad B:m+2,qquad p_0:m+1,qquad
 \mathcal R_1:m+1.                                      \tag{8}
\]

Q154 supplies the larger safe ledger `A:m+4`, `B:m+3`, `p:m+2`, leaving one
unused field derivative in every channel. This corrects the older wording
that called (8) the direct force request without altering any frozen field.

## 4. Complete terminal force jets

Near the terminal corner, `c=1` and `d=0`, so

\[
                         f_-=\mathcal R_1(u_0,p_0).      \tag{9}
\]

Q153 gives, for every fixed force jet and arbitrary `N`, a complete
compact-`X` residual estimate `C_N q^(N+1)`. With

\[
 \rho=\tau+|z|^{1/D},\qquad q\le4\rho,
\]

this becomes

\[
 |D_x^\alpha\partial_t^j f_-|
 \le4^{N+1}C_N\rho^{N+1}.                              \tag{10}
\]

Because `N` is arbitrary, every fixed force jet is Cartesian flat at the
origin. On every punctured terminal compact, Q154 supplies the uniform
potential, direct-swirl and pressure bounds at (8), including the next time
derivative. The exact exterior contributes zero localized residual outside
the cutoff transition. A finite compact cover and uniform fundamental-
theorem-of-calculus limits therefore give

\[
 F_j=\lim_{t\uparrow1}\partial_t^j f_-
 \in C_c^\infty(K),qquad D_x^\alpha F_j(0)=0.          \tag{11}
\]

This regularizes the force trace. It does not make `u_0` or `p_0` regular at
the terminal origin.

## 5. One compact smooth future force

For `sigma=t-1>=0`, choose `chi_0=1` on `[0,1/2]` and zero on `[1,infinity)`.
With increasing dyadic `b_j`, define

\[
 f_+(x,1+\sigma)=\sum_{j\ge0}
 \chi_0(b_j\sigma){\sigma^j\over j!}F_j(x).             \tag{12}
\]

On the whole support of the `j`th term, a spatial order `a` and time order
`m` obey the finite bound

\[
 \|D_x^a\partial_\sigma^mT_j\|_\infty
 \le C_{j,a,m}b_j^{m-j}.                               \tag{13}
\]

At index `j`, impose `|a|+m<=floor(j/2)`. Then `m<j`, so the power of `b_j`
is negative. Choose `b_j` after this finite set of actual terminal norms is
known so each admitted term is at most `2^-j`. Every fixed mixed derivative
eventually lies in this single diagonal and has a geometric tail. Hence

\[
 f\in C_c^\infty(\mathbb R^3\times(0,\infty)),
 \qquad \operatorname{supp}f\subset K\times[1-\tau_0,3/2], \tag{14}
\]

and `f=f_-` before `t=1`. No convergence of the exterior velocity or pressure
Taylor series is used. Equation (12) constructs no future velocity or
pressure.

## 6. Preterminal energy with the correct endpoint claim

For every `t<1`, compact support and `div u=0` give

\[
 {1\over2}E'(t)+D(t)=\langle f(t),u(t)\rangle,qquad
 E=\|u\|_2^2,quad D=\|\nabla u\|_2^2.                 \tag{15}
\]

Let `F(t)=int_0^t ||f(s)||_2 ds`. Regularizing `sqrt(E)` at the zero initial
datum gives `sqrt(E(t))<=F(t)`. Substitution in (15) yields the fixed-time
estimate

\[
 E(t)+2\int_0^tD(s)\,ds\le F(t)^2.                    \tag{16}
\]

Consequently

\[
 \sup_{t<1}E(t)\le F(1)^2,
 \qquad \int_0^1D(t)\,dt\le\tfrac12F(1)^2.             \tag{17}
\]

The two separate extrema in (17) cannot be added while retaining the single
budget `F(1)^2`; their maxima need not occur at the same time. This preserves
the corrected energy statement from the certified force theorem.

This checkpoint conditionally establishes same-field localization, zero
initial velocity, one compact smooth global force and preterminal energy and
dissipation bounds. It does not construct a future solution, recover a
competitor's pressure, prove comparison, regularize the terminal origin or
establish blow-up.

The next obligation is fixed-gauge pressure recovery and an expanding-ball
comparison with every admissible smooth finite-energy solution having the
same zero datum and the same force. That proof must not assume competitor
pressure, derivative or boundary-flux bounds.

Passing tests checks exact localization algebra, derivative requests,
corner exponents, Borel scheduling, energy arithmetic and frozen provenance.
It is not an independent continuum proof.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q155_physical_force_energy_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q155_physical_force_energy_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q155_physical_force_energy_join/evidence.json \
  --output q155_replay.json
```
