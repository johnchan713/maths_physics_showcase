# Current-head Q8 localization and compact-force join

Parent: `1b1014d62f3eb33fd776caa3b850ed4ac8f37994`.
Status: `current-head-Q8-twelve-term-localization-and-compact-force-joined-conditionally`.

This checkpoint applies the exact axisymmetric localization identity to the
Q7 field with selected parameters `h=1/200`, `kappa=1/100000`. It retains
the unchanged supports and cutoff sequence, unit viscosity, the complete
fixed-gauge pressure and covariance memory, actual axial velocity, both
nonlinear cross orientations and the curl-commutator square. It then joins
the conditional terminal force jets by a single compact Borel extension.

The result remains conditional on the punctured endpoint inequality in
Section 4. The executable checks freeze algebra, finite derivative requests
and provenance. They do not prove that continuum inequality or construct the
infinite field numerically.

## 1. Divergence-free localization of the Q7 field

Write

\[
 u_0=\nabla\times A+B,\qquad B=b(r,z,t)e_\theta,
 \qquad p_0=p_{\rm complete}.
\]

For the inherited Cartesian-smooth axisymmetric cutoff
`c=chi_r(r^2) chi_z(z) chi_t(t)`, define

\[
 u=\nabla\times(cA)+cB=cu_0+d,
 \qquad d=\nabla c\times A,\qquad p=cp_0.             \tag{1}
\]

Both velocity summands are divergence free: the first is a curl and
`div(cB)=r^{-1} partial_theta(cb)=0`. Moreover

\[
 d\cdot\nabla c=(\nabla c\times A)\cdot\nabla c=0. \tag{2}
\]

The time cutoff gives zero initial velocity. On the protected late core,
`c=1`, so the selected Q7 field is unchanged.

## 2. Exact twelve-term residual

For `R_1(v,p)=v_t-Delta v+(v dot grad)v+grad p`, direct expansion gives

\[
\begin{aligned}
\mathcal R_1(u,p)={}&c\mathcal R_1(u_0,p_0)+c_tu_0+d_t
 -(\Delta c)u_0-2\sum_i c_i\partial_i u_0-\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
 +c(u_0\cdot\nabla c)u_0+c(u_0\cdot\nabla)d\\
&+c(d\cdot\nabla)u_0+(d\cdot\nabla)d.              \tag{3}
\end{aligned}
\]

These are exactly twelve nonzero groups: the base residual; two time terms;
three viscous cutoff terms; the whole pressure cutoff; quadratic rescaling;
velocity-cutoff transport; both old/curl cross orientations; and the curl
square. The apparent thirteenth term `(d dot grad c)u0` vanishes by (2).
The viscosity is one throughout.

The pressure contribution is the complete `p0 grad c`, without discarding a
component or resetting the gauge. Divergence of the exact forced equation
retains

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f_- . \tag{4}
\]

No divergence-free-force premise or compact-support-preserving Leray
projection is used.

## 3. Derivative-ledger correction

The older localization checkpoints called

\[
 A:m+3,\qquad B:m+2,\qquad p_0:m+1                \tag{5}
\]

the direct force request. That wording is conservative, not sharp. Applying
a mixed derivative of total order `m` to (3), the largest field requests are

\[
 \boxed{A:m+2,\qquad B:m+1,\qquad p_0:m.}          \tag{6}
\]

The `A` maximum comes from `Delta d` or the differentiated gradient of
`u0`; the `B` maximum comes from the gradient of `u0`; and `p0 grad c`
contains no field derivative beyond those placed on `p0`. Uniform Cauchy
convergence reserves one further time derivative, giving

\[
 A:m+3,\qquad B:m+2,\qquad p_0:m+1,                 \tag{7}
\]

together with residual jet `m+1`. Q7 supplies the larger orders
`A:m+4`, `B:m+3`, `p:m+2`, so one unused field derivative remains in every
channel. Historical frozen evidence is preserved; Q8 records the corrected
interpretation and uses the Q7 reserve.

## 4. Conditional terminal force jets

Let `f_-=R_1(u,p)` for `t<1`. Near the Cartesian corner `c=1,d=0`, hence
`f_-` is the complete Q7 residual. For every fixed mixed jet and arbitrary
power, Q7 gives its Cartesian flatness. Away from the origin, the remaining
premise is still

\[
 \boxed{
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p_0)|<\infty
 }                                                        \tag{8}
\]

for every fixed order on every punctured terminal compact `E`. Applying (8)
at the orders in (7), with the residual `m+1` jet, gives uniform Cauchy
limits. Together with the exact exterior and a finite compact cover, this
yields conditionally

\[
 F_j=\lim_{t\uparrow1}\partial_t^j f_-
 \in C_c^\infty(K),\qquad D_x^\alpha F_j(0)=0.      \tag{9}
\]

This proves force-trace regularity under (8). It does not regularize `u0` or
`p0` at the terminal origin.

## 5. One compact smooth future force

For `sigma=t-1>=0`, take `chi0=1` on `[0,1/2]` and zero on
`[1,infinity)`. Set

\[
 f_+(x,1+\sigma)=\sum_{j\ge0}
 \chi_0(b_j\sigma){\sigma^j\over j!}F_j(x),          \tag{10}
\]

where the increasing dyadic `b_j` are chosen once from the actual finite
terminal norms. On the whole support of the `j`th summand,

\[
 \|D_x^a\partial_\sigma^m T_j\|_\infty
 \le C_{j,a,m}b_j^{m-j}.                              \tag{11}
\]

At index `j`, impose `|a|+m<=floor(j/2)` and choose `b_j` so every admitted
term is at most `2^-j`. Here `m<j`, so the exponent in (11) is negative.
Every fixed derivative eventually enters this one diagonal and has a
geometric tail. Thus (10) matches all `F_j` and gives conditionally

\[
 f\in C_c^\infty(\mathbb R^3\times(0,\infty)),
 \qquad \operatorname{supp}f\subset K\times[1-\tau_0,3/2]. \tag{12}
\]

No raw Taylor convergence, uniform-in-`j` norm, numerical infinite scale
list, future velocity or future pressure is assumed.

## 6. Scope and next obligation

Q8 conditionally closes the exact localization identity and compact-force
extension on the Q7 package, while correcting the derivative-ledger wording.
It does not prove (8), terminal velocity or pressure regularity, preterminal
energy, pressure recovery, comparison, exclusion or blow-up.

The next justified obligation is the fixed-time preterminal energy identity
for this exact localized field and force, retaining the factor-two
dissipation estimate and taking terminal energy and dissipation limits
separately. Fixed-gauge pressure recovery and expanding-ball comparison can
only be replayed after that.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q8_localization_force_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q8_localization_force_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q8_localization_force_join/evidence.json \
  --output fresh.json
```
