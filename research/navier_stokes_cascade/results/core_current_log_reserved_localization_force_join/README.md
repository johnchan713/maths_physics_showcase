# Current log-reserved localization and compact-force join

Parent: `1cd71054cd5f222901d2e1480d41ccd23702be23`.
Status: `current-log-reserved-localization-and-compact-force-joined-conditionally`.

The parent joins Cartesian flatness of the complete current six-group
residual to the punctured one-sided potential, direct-swirl and fixed-gauge
pressure jets of the same field. This checkpoint uses exactly that package
to replay axisymmetric localization and to glue the resulting terminal force
jets to one smooth compact future force.

The conclusion remains conditional on the frozen core/root, regional/frame,
all-fixed-order physical and endpoint bounds. Exact tests freeze the
algebraic identity, derivative ledger, Borel diagonal and provenance; they
do not numerically construct the infinite field or prove those continuum
inputs.

## 1. Divergence-free localization on the current field

Keep

\[
 u_0=\nabla\times A+B,\qquad B=b(r,z,t)e_\theta,
 \qquad p_0=p_{\rm complete},                           \tag{1}
\]

where `p0` remains in its fixed gauge with every background, wave,
covariance, axial and moving-moment pressure contribution. With the inherited
Cartesian-smooth axisymmetric cutoff
`c=chi_r(r^2) chi_z(z) chi_t(t)`, define

\[
 u=\nabla\times(cA)+cB=cu_0+d,\qquad
 d=\nabla c\times A,\qquad p=cp_0.                      \tag{2}
\]

The curl term is divergence free. Since `c` and `b` are axisymmetric,
`div(cB)=r^{-1}partial_theta(cb)=0`. Also

\[
 d\cdot\nabla c=(\nabla c\times A)\cdot\nabla c=0.     \tag{3}
\]

The time cutoff gives zero initial velocity, while `c=1` on the protected
late core, so the current field there is unchanged.

## 2. All twelve residual groups

For unit viscosity, direct expansion gives

\[
\begin{aligned}
\mathcal R_1(u,p)={}&c\mathcal R_1(u_0,p_0)+c_tu_0+d_t
 -(\Delta c)u_0-2\sum_i c_i\partial_i u_0-\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
 +c(u_0\cdot\nabla c)u_0+c(u_0\cdot\nabla)d\\
&+c(d\cdot\nabla)u_0+(d\cdot\nabla)d .                \tag{4}
\end{aligned}
\]

These are exactly twelve groups: the base residual; cutoff and curl time;
two cutoff-viscosity and one curl-viscosity terms; the complete pressure
cutoff; quadratic rescaling; velocity-cutoff transport; both old/curl cross
orientations; and the curl square. The apparent thirteenth term
`(d dot grad c)u0` vanishes by (3).

The pressure term is the whole `p0 grad c`; no component is dropped and the
gauge is not reset. Taking divergence of the exact forced equation retains

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f_- .\tag{5}
\]

The force need not be divergence free, and no compact-support-preserving
Leray projection is assumed. Residual flatness at viscosity one is not
silently transferred to another viscosity.

## 3. Compact terminal force jets

For a mixed force derivative of total order `m`, formula (4) requests

\[
 A:m+3,\qquad B:m+2,\qquad p_0:m+1.                    \tag{6}
\]

Uniform terminal Cauchy convergence uses one further time derivative:

\[
 A:m+4,\qquad B:m+3,\qquad p_0:m+2,                    \tag{7}
\]

plus residual jet `m+1`. These are exactly the finite orders supplied by the
current parent; changing `m` changes only the finite comparison request.

Let `f_-=R_1(u,p)` for `t<1`. Near the corner, `c=1,d=0`, so `f_-` is the
complete current residual and every fixed jet is Cartesian flat. Away from
the origin, the punctured endpoint bounds and their next time derivatives
give uniform Cauchy limits. The exact heat exterior closes the equatorial
part of the compact cover. Therefore

\[
 F_j=\lim_{t\uparrow1}\partial_t^j f_-
 \in C_c^\infty(K),\qquad D_x^\alpha F_j(0)=0.         \tag{8}
\]

This proves smooth force traces. It does not regularize `u0` or `p0` at the
terminal origin. The punctured input remains conditional on

\[
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p_0)|<\infty              \tag{9}
\]

at every fixed order.

## 4. One smooth compact future force

For `sigma=t-1>=0`, let `chi0` equal one on `[0,1/2]` and vanish on
`[1,infinity)`. Define

\[
 f_+(x,1+\sigma)=\sum_{j\ge0}
 \chi_0(b_j\sigma){\sigma^j\over j!}F_j(x),            \tag{10}
\]

with increasing dyadic `b_j`. On the support of the `j`th term,

\[
 \|D_x^a\partial_\sigma^m T_j\|_\infty
 \le C_{j,a,m}b_j^{m-j}.                               \tag{11}
\]

At each `j`, impose the finite constraints `|a|+m<=floor(j/2)` and choose
`b_j` from the actual finite norms so every admitted differentiated term is
at most `2^-j`. Since `m<j`, the exponent in (11) is negative. Every fixed
mixed derivative eventually enters this single diagonal, and its tail is
geometric. Thus (10) matches all `F_j` and yields one force

\[
 f\in C_c^\infty(\mathbb R^3\times(0,\infty)),\qquad
 \operatorname{supp}f\subset K\times[1-\tau_0,3/2].   \tag{12}
\]

No raw Taylor convergence, uniform terminal-order constant, numerical
infinite scale list, future velocity or future pressure is used.

## 5. Scope and next obligation

This establishes a conditional localization and compact-force implication on
the newly reconciled current log-reserved chain. The older current-force
checkpoint also recorded energy, but that downstream result is not silently
promoted here. This checkpoint proves neither preterminal energy on the new
chain nor comparison, competitor exclusion, or blow-up.

The next obligation is to replay the fixed-time preterminal energy identity
for this exact localized field and force, retaining the dissipation factor
two and the separate terminal energy and dissipation limits. Fixed-gauge
pressure recovery and expanding-ball comparison follow only after that.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_log_reserved_localization_force_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_log_reserved_localization_force_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_log_reserved_localization_force_join/evidence.json \
  --output core_current_log_reserved_localization_force_join_fresh.json
```
