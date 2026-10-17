# Same-field compact smooth force and preterminal energy foundation

Parent: `4899fa23a49f654c53bd3332c4c269bf9072e744` (`endpoint_foundation`).
This checkpoint joins spatial/time localization, zero initialization, complete
force jets, smooth future **force** extension and the preterminal energy bound
to the **same current realization**, relative to its retained exact-core and
finite defining/existence/inverse/gap lemmas. Those lemmas are dependencies,
not silently discharged here. No unconditional force, complete residual,
global comparison, breakdown or blow-up theorem is certified.

The earlier `current_localized_force` and `localized_force_extension` records
remain unchanged: they give conditional interfaces and diagnostic constructions.
The advance here is to supply those interfaces from the current endpoint,
realization, physical-pulse and complete-pressure foundations, without changing
the already chosen infinite sequences. The current enlarged chart margin and
the whole-support differentiated Borel estimate are made explicit below.

## 1. Current field, domain and localization

Write the actual preterminal field as

\[
 u_0=\nabla\times A+B,\qquad B=b(r,z,t)e_\theta,\qquad p_0=p_{\rm complete}.
\]

Here `A` includes every background, seed, wave, mean and diagonal potential;
`B` is the theta-independent direct swirl. `p0` has the fixed gauge, both
background tails `Pi1/Pi2`, and every four-operation wave/mean pressure
difference, axial contribution and moving-moment repair. It is not an
unforced Riesz pressure substituted for the constructed pressure.
All parameters, width `C_setup^-20`, labels, physical/mean supports, exact
Gaussian ledger and diagonal choices remain those of the parent.

Put `tau=1-t`, `Aexp=1/2+h`, `D=1/2-h` (the potential is denoted `A`, not
`Aexp`), and use the already fixed observation constants

\[
 q_e=\min(q_\star/16,2^{-\ell_0-2}),\quad
 \tau_0=q_e/16,\quad z_0=(q_e/16)^D.
\]

The parent proves `q<=3q_e/8` on that late cylinder. The **current** closed
band window is `1/4<=q/Q<=4`, so

\[
 4Q\le16q\le6q_e\le3q_\star/8<q_\star,
 \qquad Q\le(3/8)2^{-\ell_0}.
\]

The old narrow-window statement `4Q<=3q_e` cannot simply be reused; the
unchanged observation margin nevertheless suffices. Choose a Cartesian-smooth
axisymmetric cutoff `c=chi_r(r^2) chi_z(z) chi_t(t)` whose spatial support `K`
and its neighborhood lie inside this cylinder, equal to one for
`r<=r0/2`, `|z|<=z0/2`, `tau<=tau0/2`, and zero for `tau>=tau0`.
One can choose `r0` smaller than the parent radial radius before fixing `K`.
Define the fields on a neighborhood of the support closure **before** zero
extension; multiplying an undefined exterior field by zero is not a definition.

\[
 u=\nabla\times(cA)+cB=cu_0+d,\quad
 d=\nabla c\times A,\qquad p=cp_0.
\]

`div(cB)=0` because `c` and `b` are theta-independent, hence `div u=0`.
Also `d dot grad c=0` exactly. A naive velocity cutoff `cu0`, or a
theta-dependent cutoff on the direct swirl, need not be divergence free.
These fields are smooth and compactly supported for every `t<1`, and
`u(0)=0` because the time cutoff vanishes on an initial neighborhood.

## 2. Complete localized residual, including pressure memory

For `R_nu(v,p)=partial_t v-nu Delta v+(v dot grad)v+grad p`, direct expansion
gives the following **twelve** terms (the actual viscosity is `nu=1`):

\[
\begin{aligned}
R_\nu(u,p)={}&cR_\nu(u_0,p_0)+c_tu_0+d_t
-\nu(\Delta c)u_0-2\nu\sum_i c_i\partial_i u_0-\nu\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
+c(u_0\cdot\nabla c)u_0\\
&+c(u_0\cdot\nabla)d+c(d\cdot\nabla)u_0+(d\cdot\nabla)d.
\end{aligned}
\]

The apparent additional term `(d dot grad c)u0` is exactly zero. In
particular, **all of** `p0 grad c` is retained, not just a convenient pressure
piece; a gauge change before cutting pressure would change this term.
A force jet of total order `m` needs potential jets through `m+3`, direct
swirl through `m+2`, and pressure through `m+1`. Uniform terminal Cauchy
bounds reserve respectively `m+4,m+3,m+2`, supplied by the parent on
punctured compact sets. These are arbitrary fixed finite orders, not one
uniform constant for every derivative order.

Set `f_minus=R_1(u,p)` for `t<1`. Taking divergence gives

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f_{-}.
\]

There is no assumption `div f=0`, no unforced pressure substitution, and no
Leray projection asserted to preserve compact support. The official
[Fefferman formulation](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf)
(equations (4), (5), (7), alternative (C); checked 2026-10-17) requires smooth
rapidly decaying data and all force space/time derivatives, but does not
require divergence-free forcing. Satisfying the force criteria is not by
itself a breakdown result.

## 3. Global uniform force jets at the terminal time

In a fixed neighborhood of `(x,t)=(0,1)`, `c=1` and all its positive-order
derivatives vanish. Thus `f_minus` equals the **complete current** residual
there. The realization foundation supplies, for all fixed multiindices,
time orders and `N`, uniform `O(q^N)` residual jets on compact `X<=Xext`;
outside the exact support boundary the field is the exact heat background
and its residual is zero. These regions exhaust the corner, including the
closed support boundary. Using the parent's inequality

\[
 q\le2\tau+4|z|^{1/D}
\]

gives a uniform bound `C_(alpha,j,N)(tau+|z|^(1/D))^N` throughout that fixed
corner box, not just on a chosen approach path. These Cartesian jets and
their next time derivatives are bounded as `t->1-`. Away from the origin,
the parent bounds the full field, potential, pressure and their next time
derivatives uniformly on punctured terminal compacts. A finite compact
cover of the rest of `K`, together with zero outside `K`, gives global
uniform Cauchy limits

\[
 F_j=\lim_{t\uparrow1}\partial_t^j f_-\in C_c^\infty(K).
\]

Spatial fundamental-theorem-of-calculus identities pass to the uniform
limits, including on segments through the origin, and identify every
spatial derivative. The next time-jet estimate identifies consecutive
terminal jets. Every spatial jet of every `F_j` is zero at the origin:
set `tau=0` in the corner bound and then take `N` arbitrarily large.
This argument regularizes the **force**, not the singular velocity or
pressure at the origin. It never forward-inverts or newly sums the Gaussian
ledger, nor changes any current pulse or mean-correction sequence.

## 4. One smooth future force extension; no future velocity assertion

Choose `chi0` smooth on `[0,infinity)`, one on `[0,1/2]` and zero on `[1,infinity)`.
For `sigma=t-1>=0` set

\[
 f_+(x,1+\sigma)=\sum_{j\ge0}T_j,
 \qquad T_j=\chi_0(b_j\sigma)\frac{\sigma^j}{j!}F_j(x).
\]

Take `b0=2` and `bj=2^n_j` increasing. The exact differentiated term is

\[
 \partial_\sigma^m T_j=F_j\sum_{\substack{0\le k\le m\\m-k\le j}}
 {m\choose k}b_j^k\chi_0^{(k)}(b_j\sigma)
 \frac{\sigma^{j-m+k}}{(j-m+k)!}.
\]

On the **whole support** `0<=sigma<=1/bj`, not merely the transition annulus,
each summand is bounded by a constant times `bj^(m-j)`. Consequently

\[
 \|D_x^\alpha\partial_\sigma^m T_j\|_\infty
 \le C_{j,\alpha,m}b_j^{m-j}.
\]

For each `j>=1` impose the finite set `|alpha|+m<=floor(j/2)` and make these
bounds at most `2^-j`. For example choose the integer `n_j>n_(j-1)` with
`n_j log 2 >= (log C_(j,alpha,m)+j log 2+2)/(j-m)` for every nonzero
constraint. Every denominator is positive. No upper bound uniform in `j`
on the terminal norms is needed. For fixed `alpha,m`, all
`j>=max(1,2*(|alpha|+m))` are controlled, and the tail starting at `J` is
at most `2^(1-J)`. The finitely many earlier terms are smooth.

Thus **one** increasing sequence yields uniform convergence of every fixed
mixed differentiated series on the full half-line, including zero.
The plateau and factorial give `partial_sigma^m f_plus(x,1)=F_m(x)`.
Gluing to `f_minus` is therefore smooth in all variables, has support
`K x [1-tau0,3/2]`, and retains origin-flatness. All polynomially weighted
derivative bounds follow from compact smooth support; for spatial radius
`R`, a weight of degree `L` costs at most `(1+R+3/2)^L` times the relevant
finite sup norm. This is an existence argument for scales depending on
actual terminal norms, **not** their numerical evaluation.

No raw Taylor series convergence, analytic heat continuation, zero extension
of nonzero terminal pressure, or future velocity/pressure construction is
used. A zero-radius analytic heat Taylor series does not obstruct smooth
Borel gluing of force jets.

## 5. Preterminal energy without circularity or lost dissipation

Only **after** constructing the compact smooth force, put
`E(t)=||u(t)||_2^2`, `D(t)=||grad u(t)||_2^2`, and
`F(t)=integral_0^t ||f(s)||_2 ds`. For every `t<1`, actual compact support
justifies integration by parts; pressure and transport fluxes vanish:

\[
 \tfrac12E'(t)+D(t)=\langle f(t),u(t)\rangle.
\]

Regularize the initially zero square root: from this identity and `D>=0`,
`d_t sqrt(E+delta^2)<=||f||_2`. Letting `delta->0` yields `sqrt(E(t))<=F(t)`.
Integrating the exact identity then gives

\[
 E(t)+2\int_0^tD(s)\,ds
 \le2\int_0^t\|f(s)\|_2F(s)\,ds=F(t)^2.
\]

In particular `sup_(t<1) E(t)<=F(1)^2` and, by monotonicity of the
dissipation integral, `integral_0^1 D<=F(1)^2/2`. **Do not** replace the
left side of the preceding fixed-time inequality by `sup E+2 integral_0^1 D`
while keeping the same right side. The two extrema need not occur together.

An independent negative control uses the scalar equation `y'+y=f`:
let `y=t/epsilon` up to `epsilon=1/100`, then let
`y=exp(-(t-epsilon))` up to `t=1`. Its nonnegative force has
`integral f=1+epsilon/2=201/200`. But `sup y^2=1` and, using
`exp(-s)<=1/(1+s)`,

\[
 \sup y^2+2\int_0^1y^2\,dt
 \ge2-\frac1{1+2(1-\epsilon)}+\frac{2\epsilon}{3}
 =\frac{37349}{22350}>\frac{40401}{40000}=(\int f)^2.
\]

This piecewise absolutely continuous scalar diagnostic is not a Navier--Stokes
solution; its strict gap persists under smoothing. The historical fixed-time
energy estimate is preserved, not retracted.

On the inherited protected core path, for sufficiently small
`tau<=min(tau0/2,r0^2/(8 Xin))`, localization is identically one. The current
conditional lower bound `u_theta >= (e0/2) tau^-Aexp` is unchanged. Its
positive `e0` still belongs to the retained exact-core lemma; the energy
estimate does not prove that lemma or independently establish blow-up.

## 6. Evidence, meaning and next obligation

`localized_checks.py` independently expands the twelve terms using exact
rational multivariate polynomials at viscosities `3/5,1,17/10` (nine vector
components). Each of the twelve omitted terms is nonzero in the fixture.
It detects broken divergence from naive and nonaxisymmetric cutoffs and
checks the full forced pressure Poisson identity. Separate polynomial
transition controls check 98 exact Leibniz identities and 49 terminal
factorial identities. Exact cube integrals check zero pressure/transport
energy work, positive viscous dissipation, both energy coefficients and
zero initial velocity. The zero-extended cube fixture has finite regularity:
it is **not** the actual compact smooth flow. Synthetic extreme terminal
norms check the finite Borel scale constraints; they are not actual norms.
None of these tests supplies a missing continuum supremum or retained lemma.

The new evidence binds eight own sources, eight dependency records and two
imported helpers (18 direct bindings); 997 inherited provenance bindings are
preserved and replayed. The validator rejects scope promotion, altered gates,
changed source bytes, substituted history keys and evidence overwrite.

Reproduction from this directory:

```sh
python -m unittest -v test_audit.py
python audit.py --output /fresh/force-foundation.json --verify-record evidence.json
```

Next: recover pressure and close the expanding-ball comparison for **every
admissible global smooth finite-energy competitor with this same zero datum
and this same force**, without assuming unproved competitor pressure,
derivative or boundary-flux bounds. Separately discharge the retained exact-core
and finite existence/inverse/gap dependencies. Until both tasks are proved,
the admissible force and protected growth here remain a relative theorem,
not an unconditional exclusion of global solutions.
