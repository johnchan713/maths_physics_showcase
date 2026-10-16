# Current-field localization, smooth compact force and whole-space energy

Status: `conditional-current-field-localization-smooth-compact-force-and-energy`.

This checkpoint joins the [current complete field and terminal jets](../current_endpoint_jets/README.md)
to ONE smooth compactly supported force on all nonnegative times, with
zero initial velocity and a spatially compact smooth solution for t<1.
Its exact equation includes all wave/mean pressure, curl cutoffs, viscosity
and nonlinear interactions. Whole-space energy and integrated dissipation
are bounded up to t=1, and the conditional protected core growth survives.

Every field conclusion is CONDITIONAL on the retained leading-profile,
summed-background, all-order inverse and uniform physical-jet estimates.
This is not an independent proof of those premises or a global solution
exclusion theorem. Actual infinite fields, terminal norms and extension
scales are not numerically evaluated. Tests are manufactured diagnostics,
not continuum certificates. Historical frozen evidence remains unchanged.

## 1. Use the same current field, pressure and admitted domain

Keep viscosity one and the SAME h, N, positive lambda, width C^-20,
frozen primary roots and labels. The parent gives

\[
u_0=\operatorname{curl}\mathcal A+\mathcal B,
\quad\mathcal B=b e_\theta,\quad\partial_\theta b=0,
\quad p_0=p_{\rm current}.\tag{1}
\]

A contains the summed background potential, every seed and later wave
potential, all azimuthal mean potentials and every diagonal cutoff.
The direct swirl is retained separately; no compact full vector potential
for it is asserted. p0 contains the SAME full p_B and both exceptional
base-pressure tails, plus ALL changing wave and six-covariance mean
pressures and axial memory. No pressure gauge is changed.

Let Aexp=1/2+h, D=1/2-h, tau=1-t and
`q=tau+z^2 q^(2h)`, `X=r^2/(2q)`. Use the parent's strict observation
domain, not the ambiguous shortcut q<qstar/2:

\[
q_e=\min(qstar/16,2^{-\ell_0-2}),\quad
\tau_0=q_e/16,\quad z_0=(q_e/16)^D.\tag{2}
\]

For tau<=tau0 and |z|<=z0,
`q<=2tau+4|z|^(1/D)<=3q_e/8`. All active Q satisfy
`4Q<=3q_e<=3qstar/16<qstar` and
`Q<=3*2^-ell0/16`. Thus every enlarged chart and fixed first-band
threshold remains admitted. Restricting the observation/localization
region does not alter the original construction or its parameters.

The parent supplies compatible one-sided jets of A,B,p0 on every
punctured terminal compact, with all fixed physical derivatives uniform.
Near the axis at q>=c>0 the annular fields vanish uniformly; the inherited
Cartesian background factors supply regularity. The fractional auxiliary
map need not be smooth at the axis. Complete residual derivatives are
O(q^N) for every fixed N, UNIFORMLY on X<=Xext; for X>=Xext the
uncut residual is exactly zero. This compact-X uniformity, not pointwise
flatness at a fixed X, is used below. It is part of the retained current
physical-residual estimates, not proved by numerical samples here.

## 2. Compact support, exact divergence and zero initialization

Choose r0 strictly inside the parent's local radial domain. Take smooth
axisymmetric cutoffs c=chi_x(r^2,z)chi_t(t) such that

\[
\chi_x=1\quad(r\le r_0/2,\ |z|\le z_0/2),\qquad
\operatorname{supp}\chi_x\Subset\{r<r_0,|z|<z_0\},
\]
\[
\chi_t=0\quad(\tau\ge\tau_0),\qquad
\chi_t=1\quad(\tau\le\tau_0/2).
\]

With smooth zero extension outside this supported domain, set

\[
u=\operatorname{curl}(c\mathcal A)+c\mathcal B,
\quad p=cp_0,\quad d=\nabla c\times\mathcal A,
\quad u=cu_0+d.\tag{3}
\]

The dependence on r^2 is Cartesian-smooth at the axis. Curl is exactly
divergence-free; `div(cB)=r^-1 partial_theta(cb)=0`. A theta-dependent
cutoff for the direct swirl would fail this argument. The fixed domain
margin ensures that all fields are smooth on a neighborhood of each
nonterminal cutoff boundary, so flat cutoffs genuinely give smooth zero
extensions; multiplying an undefined field by zero is not used.

For every T<1, u,p are smooth on R^3 x [0,T], with one fixed compact
spatial support K=supp chi_x. They vanish in a neighborhood of t=0,
since tau0<1/256; in particular u(0)=0. Near (0,1) the cutoff is
IDENTICALLY one, including all derivatives: the complete field, potential
and pressure are unchanged there.

## 3. Exact complete localization residual and forced pressure

For algebraic clarity let `R_nu(v,p)=v_t+(v.grad)v-nu Delta v+grad p`.
The actual construction uses nu=1; other positive diagnostic viscosities
do not inherit its flat-residual conclusion. Direct expansion of (3) gives

\[
\begin{aligned}
R_\nu(u,p)={}&cR_\nu(u_0,p_0)+c_tu_0+d_t
-\nu(\Delta c)u_0-2\nu\sum_i(\partial_i c)\partial_i u_0
-\nu\Delta d+p_0\nabla c\\
&+c(c-1)(u_0\cdot\nabla)u_0+c(u_0\cdot\nabla c)u_0\\
&+c(u_0\cdot\nabla)d+c(d\cdot\nabla)u_0+(d\cdot\nabla)d.
\end{aligned}\tag{4}
\]

These are all twelve retained terms. The apparent additional term
`(d.grad c)u0` is exactly zero, since `d=grad c cross A`.
d_t includes mixed time/spatial cutoff derivatives; Delta d includes
third spatial cutoff derivatives. Pressure contributes p0 grad c,
with the complete p0, not merely an angular or leading approximation.
All accumulated covariances, flat velocity pieces and exact designated
Gaussian tails remain inside the actual fields. No new inverse or
independent summation of Gaussian residual tails is performed.

An m-th force jet requires A through m+3, B through m+2, p0 through
m+1; a uniform terminal Cauchy bound requires one additional time
derivative. Every fixed such requirement is supplied by the current
endpoint theorem, conditionally on its physical estimates. Define
`f_-=R_1(u,p)` for t<1. Taking divergence of this EXACT equation gives

\[
-\Delta p=\sum_{i,j}\partial_i\partial_j(u_i u_j)-\operatorname{div}f_-.
\tag{5}
\]

The force need not be divergence-free. Dropping its divergence or replacing
the pressure with an unforced Riesz formula would change the equation.
No unproved compact-support-preserving Leray projection is applied.

## 4. Terminal force jets, including the singular corner

Away from x=0, the parent's compatible terminal jets and (4) give the
terminal jets of f_-. Near the corner c=1,d=0, so f_- is EXACTLY the
complete current residual. Combining its compact-X uniform flatness with
the exact exterior zero gives, for every fixed alpha,j,N,

\[
|D_x^\alpha\partial_t^j f_-(x,1-\tau)|
\le C_{\alpha,j,N}(\tau+|z|^{1/D})^N.\tag{6}
\]

This holds throughout a fixed neighborhood where all cutoff derivatives
vanish, not only along a special core path. To obtain UNIFORM convergence
of a given force derivative on R^3, choose a small ball around x=0 so
(6) is small at late times. On K outside that ball use the punctured
endpoint theorem and a finite cover; outside K the derivative is zero.
Thus all spatial/time derivative limits exist uniformly on R^3.

Pass to the limit in spatial coordinate integrals, including segments
crossing x=0. Uniform convergence of every spatial derivative identifies

\[
F_j=\lim_{t\uparrow1}\partial_t^j f_-(\cdot,t)\in C_c^\infty(\mathbb R^3),\quad\operatorname{supp}F_j\subset K,
\qquad D_x^\alpha F_j(0)=0.\tag{7}
\]

Compatibility of consecutive time jets follows by integrating the next
time derivative and passing to its uniform limit. This attaches a smooth
ONE-SIDED force through the entire terminal slice, not a smooth velocity
or pressure at the singular origin. Origin flatness alone would not
control limits elsewhere; the punctured endpoint input is essential.
F_j generally need not vanish away from the origin. Neither the full
force nor its jets may simply be set to zero immediately after t=1.

## 5. ONE compact smooth future force, without Taylor convergence

For sigma=t-1>=0 let chi0 be one on [0,1/2] and zero on [1,infinity).
Choose b0=2 and increasing dyadic b_j for j>=1, then define

\[
f_+(x,1+\sigma)=\sum_{j\ge0}
\chi_0(b_j\sigma)\frac{\sigma^j}{j!}F_j(x).\tag{8}
\]

For m time derivatives, the EXACT Leibniz formula for a scalar factor is

\[
\sum_{\substack{0\le k\le m\\m-k\le j}}
\binom mk b_j^k\chi_0^{(k)}(b_j\sigma)
\frac{\sigma^{j-m+k}}{(j-m+k)!}.
\]

On its support 0<=sigma<=1/b_j, each term is bounded by a constant
times b_j^(m-j). For any fixed spatial order a, therefore,

\[
\|D_x^a\partial_\sigma^m T_j\|_\infty\le C_{j,a,m}b_j^{m-j}.
\tag{9}
\]

For each j>=1 impose this bound <=2^-j for the FINITE list
`a+m<=floor(j/2)`. Its exponents m-j are negative. With
`b_j=2^n_j`, it suffices to choose an integer n_j>n_(j-1) exceeding
`(log C+j log 2+2)/((j-m)log 2)` for each listed padded C>=1.
This works for arbitrarily large finite norms; the synthetic scale
diagnostic is not an evaluation of the actual F_j or constants.

For each fixed a,m all terms j>=max(1,2(a+m)) are controlled. Their
differentiated tail from j=J is at most `sum_(j>=J)2^-j=2^(1-J)`,
uniformly on the WHOLE future half-domain, including sigma=0 and every
cutoff transition. Thus the sum is C-infinity with uniform differentiated
convergence. At sigma=0 the plateau makes its m-th time derivative
exactly F_m. Equation (7) then gives smooth matching to f_-.
All terms retain support in K and all spatial jets vanish at the origin.
Since b_j>=2, the sum vanishes for sigma>=1/2, smoothly at that boundary.

Consequently ONE force f, equal to f_- for t<1 and f_+ for t>=1, satisfies

\[
f\in C_c^\infty(\mathbb R^3\times(0,\infty);\mathbb R^3),\qquad
\operatorname{supp}f\subset K\times[1-\tau_0,3/2].\tag{10}
\]

Its zero initial neighborhood also gives smoothness on nonnegative times.
No convergence of the raw terminal Taylor series, uniform constants over
jet order, or analytic continuation is assumed. The parent's zero-radius
heat/pressure Taylor series do not obstruct this smooth-force construction;
they do not by themselves prove nonanalyticity of the resulting FORCE.
No future velocity or pressure is supplied by (8).

## 6. Whole-space energy, dissipation and protected growth

No energy estimate was used to establish force smoothness above.
On every [0,T], T<1, multiply the exact preterminal equation by u and
integrate in space. Smooth compact support and div u=0 remove pressure
and transport fluxes. At viscosity one,

\[
\frac12\frac{d}{dt}\|u(t)\|_2^2+\|\nabla u(t)\|_2^2=\langle f(t),u(t)\rangle.
\tag{11}
\]

Put `F(t)=integral_0^t ||f(s)||_2 ds`. Do not divide by ||u||2,
which is initially zero. Instead differentiate
`sqrt(||u||2^2+delta^2)`, discard its nonpositive dissipation term,
and apply Cauchy--Schwarz. Its derivative is <=||f||2. The zero datum
and delta->0 give ||u(t)||2<=F(t). Inserting that bound into (11) gives

\[
\|u(t)\|_2^2+2\int_0^t\|\nabla u(s)\|_2^2\,ds
\le2\int_0^tF'(s)F(s)\,ds=F(t)^2.\tag{12}
\]

The compact smooth force has F(1)<infinity. This bounds kinetic energy
uniformly and total dissipation on [0,1); only monotone convergence of
the nonnegative dissipation integral is used. No regularity of u at t=1
or energy statement about a future velocity is inferred.

The current inner path `x_tau=(sqrt(2Xin*tau),0,0)` is in the unchanged
plateau once `tau<=min(tau0/2,r0^2/(8Xin))`. Thus the parent's conditional
positive core coefficient still gives
`|u_theta(x_tau,1-tau)|>=(e0/2)tau^(-Aexp)` for sufficiently small tau.
This includes all current waves, means and cutoffs. It does not independently
verify e0, and bounded energy does not itself exclude concentration.

## 7. Admissibility and the remaining comparison obligation

The official [Fefferman/CMI formulation](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf),
PDF pages 1-2, equations (4),(5),(7) and alternative (C), was checked on
2026-10-16. It requires all derivatives of the data and force to have
arbitrary polynomial decay; it does not impose div f=0. Here u(0)=0
satisfies the data requirement. If K is contained in a radius-R ball,
(10) gives for every alpha,m,L the explicit finite bound

\[
\sup_{x\in\mathbb R^3,\ t\ge0}(1+|x|+t)^L|D_x^\alpha\partial_t^m f(x,t)|
\le(1+R+3/2)^L\|D_x^\alpha\partial_t^m f\|_\infty.
\]

These are the smooth-data/force requirements, CONDITIONALLY on the current
construction. They are not all requirements of a breakdown theorem.
The next bounded obligation is to join pressure recovery and the
expanding-ball comparison lemma to THIS same u,p,f and zero datum, for
an arbitrary putative global smooth bounded-energy competitor. No decay
of its pressure or global derivative/dissipation bound may be assumed
without proof; all spatial-infinity fluxes must be controlled. Independent
review of the inherited analytic foundation remains mandatory before
unconditional admissible forcing, exclusion or blow-up claims.

## 8. Reproduce

```sh
python3 research/navier_stokes_cascade/results/current_localized_force/test_audit.py -v
python3 research/navier_stokes_cascade/results/current_localized_force/audit.py \
  --verify-record research/navier_stokes_cascade/results/current_localized_force/evidence.json \
  --output /tmp/current_localized_force_recheck.json
```

Fresh manufactured checks use the current actual-q endpoint fixture for
all twelve localization terms, plus an independent exact polynomial curl
for forced pressure. Old terminal mixed-jet and compact-energy diagnostics
are replayed without substituting their old field for the current one.
The historical failed 64-node energy quadrature and successful 96-node
check retain the SAME 2e-9 gate. No tolerance is weakened and no numerical
diagnostic is promoted to a continuum proof.
