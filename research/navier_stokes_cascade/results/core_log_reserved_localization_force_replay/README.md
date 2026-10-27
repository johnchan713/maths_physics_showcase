# Log-reserved same-field localization and compact-force replay

Parent: `d6cbc11f6151c86ed7b85beaf6d03378a3e69198`.

The parent joins the corrected log-reserved Cartesian residual flatness to
all fixed punctured one-sided potential, direct-swirl and fixed-gauge
pressure jets on one already selected field.  This checkpoint uses exactly
those inputs to replay axisymmetric localization and the smooth future-force
extension.  Historical evidence is unchanged.

The conclusion is conditional on the frozen primitive, operator, endpoint
and residual contracts in the parent chain.  The executable checks reproduce
the finite derivative ledger and evidence hashes; they do not numerically
evaluate the infinite field or replace the continuum arguments below.

## 1. Divergence-free localization of the same field

Retain

\[
 u_0=\nabla\times A+B,\qquad B=b(r,z,t)e_\theta,
 \qquad p_0=p_{\rm complete},                              \tag{1}
\]

where `p0` remains in its fixed gauge and contains all background, wave,
covariance, axial and moving-moment pressure memory.  Choose the inherited
Cartesian-smooth axisymmetric cutoff
`c=chi_r(r^2) chi_z(z) chi_t(t)` inside the enlarged terminal cylinder and
define, before zero extension,

\[
 u=\nabla\times(cA)+cB=cu_0+d,qquad
 d=\nabla c\times A,qquad p=cp_0.                         \tag{2}
\]

Both summands in `u` are divergence free.  The curl is divergence free and
`div(cB)=r^{-1} partial_theta(cb)=0`.  Moreover

\[
 d\cdot\nabla c=(\nabla c\times A)\cdot\nabla c=0.         \tag{3}
\]

The time cutoff gives zero initial data.  On the protected late core `c=1`,
so localization leaves the selected field unchanged there.

## 2. Exact twelve-term residual

For unit viscosity, direct expansion of
`R_1(v,p)=v_t-Delta v+(v dot grad)v+grad p` gives

\[
\begin{aligned}
R_1(u,p)={}&cR_1(u_0,p_0)+c_tu_0+d_t
 -(\Delta c)u_0-2\sum_i c_i\partial_i u_0-\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
 +c(u_0\cdot\nabla c)u_0+c(u_0\cdot\nabla)d\\
&+c(d\cdot\nabla)u_0+(d\cdot\nabla)d.                   \tag{4}
\end{aligned}
\]

These are twelve nonzero groups: the base residual; cutoff and curl time
terms; two cutoff-viscosity terms and the curl-viscosity term; the complete
pressure cutoff; quadratic rescaling; velocity-cutoff transport; both
orientations of the old/curl cross; and the curl square.  The apparent
thirteenth group `(d dot grad c)u0` is exactly zero by (3).

The viscosity is fixed at one.  In particular, the residual-flatness result
cannot be silently reused at a different viscosity.  The pressure term is
the full `p0 grad c`; no pressure component is discarded or regauged.
Taking divergence retains the force source:

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f_-.
                                                               \tag{5}
\]

There is no divergence-free-force premise and no compact-support-preserving
Leray projection.

## 3. Closed-time force jets

For a mixed force derivative of total order `m`, the inherited sufficient
localization ledger requests

\[
 A:m+3,\qquad B:m+2,\qquad p_0:m+1.                       \tag{6}
\]

Uniform Cauchy convergence reserves one more time derivative:

\[
 A:m+4,\qquad B:m+3,\qquad p_0:m+2.                       \tag{7}
\]

The parent supplies (7) on every punctured terminal compact for each fixed
`m`.  This is a request-by-request all-orders statement about the same field;
the comparison stage changes, but the field and cutoff sequence do not.

Set `f_-=R_1(u,p)` for `t<1`.  Near the terminal corner `c=1` and `d=0`, so
`f_-=R_1(u0,p0)`.  The corrected log-reserved corner replay gives, for every
fixed spatial multiindex `alpha`, time order `j` and power `N`,

\[
 |D_x^\alpha\partial_t^j f_-|
 \le C_{\alpha,j,N}(\tau+|z|^{1/D})^N.                   \tag{8}
\]

This is uniform on the Cartesian corner, rather than only on one path.  It
joins the punctured Cauchy bounds and exact exterior heat residual through a
finite compact cover.  Therefore

\[
 F_j=\lim_{t\uparrow1}\partial_t^j f_-
 \quad\hbox{exists in }C_c^\infty(K).                    \tag{9}
\]

Spatial fundamental-theorem-of-calculus identities and the next time jet
pass to these uniform limits.  Equation (8), with arbitrary `N`, makes every
spatial derivative of every `F_j` vanish at the origin.  This proves force
trace regularity only; it does not regularize `u0` or `p0` at `(0,1)`.

## 4. One smooth compact future force

Put `sigma=t-1>=0`.  Let `chi0` be one on `[0,1/2]` and zero on
`[1,infinity)`, and set

\[
 f_+(x,1+\sigma)=\sum_{j\ge0}
 \chi_0(b_j\sigma){\sigma^j\over j!}F_j(x),              \tag{10}
\]

where `b0=2` and `b_j` is an increasing dyadic sequence.  The exact scalar
derivative is

\[
 \sum_{\substack{0\le k\le m\\m-k\le j}}
 {m\choose k}b_j^k\chi_0^{(k)}(b_j\sigma)
 {\sigma^{j-m+k}\over(j-m+k)!}.                          \tag{11}
\]

On its whole support `0<=sigma<=1/b_j`, each summand in (11) is bounded by
a fixed constant times `b_j^(m-j)`.  At index `j>=1`, impose the finite set
`|alpha|+m<=floor(j/2)` and choose `b_j` large enough that every associated
mixed derivative of the `j`th summand is at most `2^-j`.  The exponent is
negative because every admitted `m` satisfies `m<j`.

For every fixed mixed derivative, all indices

\[
 j\ge\max\{1,2(|\alpha|+m),m+1\}                         \tag{12}
\]

enter this single diagonal, and the differentiated tail from `J` is bounded
by `2^(1-J)`.  Thus (10) is smooth, its terminal time jets are exactly the
`F_j`, and gluing it to `f_-` gives one force satisfying

\[
 f\in C_c^\infty(\mathbb R^3\times(0,\infty)),\qquad
 \operatorname{supp}f\subset K\times[1-\tau_0,3/2].     \tag{13}
\]

No convergence of the raw terminal Taylor series is assumed.  The scale
sequence depends on the actual finite terminal norms but is chosen once;
no uniform bound in `j` is required.  This construction supplies no future
velocity or pressure.

## 5. Scope and next obligation

This establishes a **conditional same-field localization and compact-force
implication** with unit viscosity, all cutoff derivatives, both nonlinear
crosses, the curl square, complete pressure memory and force divergence.
It does not claim unconditional pulse estimates or residual improvement.
It does not prove terminal velocity or pressure regularity, construct a
future solution, or establish competitor exclusion or blow-up.

The next obligation is to replay the fixed-time preterminal energy identity
for this exact localized field and force, keeping the dissipation factor two
and the separate endpoint bounds.  After that come fixed-gauge pressure
recovery and expanding-ball comparison against an arbitrary admissible
finite-energy competitor with the same datum and force.

Reproduce with

```sh
python3 test_audit.py -v
python3 audit.py --verify-record evidence.json --output fresh.json
```
