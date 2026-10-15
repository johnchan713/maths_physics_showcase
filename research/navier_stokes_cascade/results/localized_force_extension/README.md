# Conditional physical localization, smooth compact force and energy

This checkpoint continues the
[full-potential realization](../full_potential_realization/README.md).
Under its inherited uniform physical derivative, profile and support
hypotheses, it proves the following transfer to physical space:

* The COMPLETE local potentials have smooth one-sided terminal limits
  away from (0,1). The existing physical jet bounds imply these limits;
  a new choice of auxiliary labels on the terminal slice is unnecessary.
* Spatial and temporal localization of those potentials gives a real,
  exactly divergence-free whole-space field with zero initial velocity.
  Its exact equation error is the restriction of a smooth compactly
  supported force on R^3 x (0,infinity), including through t=1.
* This same field has uniformly bounded whole-space kinetic energy and
  finite integrated dissipation before t=1. The inherited inner growth
  survives localization and the complete correction sequence.

These are CONDITIONAL analytic conclusions. The actual infinite field,
terminal jets and extension scales are not numerically evaluated. The
global auxiliary atlas and earlier profile inequalities are not independently
verified here. Comparison with an arbitrary global smooth bounded-energy
solution remains a separate argument. No unconditional blow-up theorem is
claimed, and manufactured checks do not establish the inherited premises.

Parent: `cef7a529973f9ff43d9107cd132663f6bfee7b6a`. The proposed
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(9.15)--(9.21), (10.1)--(10.14), is the source under audit, not a premise
for its own conclusion. Its frozen SHA-256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Earlier source/evidence and GOAL.md are preserved.

## 1. The exact input interface

Write the previously realized local field at viscosity one as

\[
u_0=\operatorname{curl}\mathcal A+\mathcal B,\qquad
\mathcal B=b e_\theta,\quad\partial_\theta b=0,\qquad p_0=p_{\rm loc}.
\tag{1}
\]

It is smooth for tau=1-t>0 and 0<q<qstar<1. Keep the same fixed
0<h<1/100, A=1/2+h, D=1/2-h and

\[
q=\tau+z^2q^{2h},\quad\eta=z/q^D,\quad X=r^2/(2q).
\tag{2}
\]

The inherited inputs are the actual CONDITIONAL conclusions already
recorded in the preceding notes:

1. All Cartesian derivatives of the assembled representatives have uniform
   physical bounds on compact X ranges, with a fixed loss at each derivative
   order. For a fixed stage these are (5), (11) of the parent; the realized
   background has the corresponding Cartesian potential bounds. All constants
   are uniform in physical point and dyadic label, not merely pointwise in
   eta. The representatives agree on chart overlaps and have the stated
   smooth zero extensions. Bounded overlap is part of that inherited geometry.
2. Only finitely many correction-stage cutoffs are nonzero where q>=c>0.
   The base summation is locally finite there too. All full-field residual
   derivatives are O(q^N) for every fixed N as q->0, and the residual is
   exactly zero for X>=Xext.
3. In that heat exterior the chosen potential is A=0 and
   \[
   \mathcal B=K(r,\tau)e_\theta,\quad
   K=C r^{-1-2h}\mathcal H(4\tau/r^2),\quad
   p_0=-\int_r^\infty K(\rho,\tau)^2\,d\rho/\rho,
   \]
   \[
   \mathcal H(Z)=\frac1{\Gamma(1+h)}\int_0^\infty
            e^{-v}v^h(1+Zv)^{-h}\,dv.\tag{3}
   \]
   C is the fixed inherited exterior constant. The zero axial moments make
   the background streamfunction vanish in this exterior; all positive-order
   and annular potentials vanish there as well. At the axis the base potential
   has Cartesian form a(r^2,z,t)(-x2,x1,0); the direct base swirl has the same
   regular factor. Annular representatives vanish near the axis for each t<1.
4. At a fixed Xin<Xa and z=0 the complete local swirl obeys
   tau^-A(e0+O(tau^(2h))) with e0>0, as proved conditionally in the parent.

Input 1 is the important quantifier. Interior smoothness alone would not
give the endpoint result: a function supported away from the origin times
sin(1/(1-t)) is smooth for t<1 and zero near the origin, yet has no terminal
limit. It violates the UNIFORM derivative bounds required above. If the
inherited global geometry fails to supply those bounds after assembly, the
precise blocking inequality is

\[
\sup_{x\in E,\,1-\delta<t<1}
 |\partial_x^\alpha\partial_t^j(\mathcal A,\mathcal B,p_0)|<\infty
\quad\text{for every fixed }\alpha,j\text{ and }E\Subset\{x\ne0\}.
\tag{4}
\]

Section 2 derives (4) from the stated inherited uniform physical bounds.
It does not independently certify those earlier bounds or the global atlas.

## 2. Endpoint limits follow from the existing physical estimates

First work on a fixed compact neighborhood where c<=q<=c'<qstar. All
stage cutoffs with a_j c>=1 vanish. For each remaining stage, its uniform
physical bound is finite throughout this interval: powers of q and
polynomials in |log q| are bounded there. The same holds for the finitely
many nonzero background terms. The inherited bounded overlap and uniform
label estimates give a bound for the ASSEMBLED physical derivative.
There is no need to evaluate or choose a label at tau=0.

For any G in (A,B,p0), the fundamental theorem of calculus gives

\[
\sup_E|\partial_x^\alpha\partial_t^jG(t_2)
             -\partial_x^\alpha\partial_t^jG(t_1)|
\le |t_2-t_1|
 \sup_{E\times[t_1,t_2]}|\partial_x^\alpha\partial_t^{j+1}G|.\tag{5}
\]

The right-hand bound is uniform as t1,t2->1. Thus every fixed derivative
is uniformly Cauchy and has a one-sided limit. Apply the same theorem on
small spatial coordinate segments and pass to uniform limits to identify
the spatial derivatives of these limits. Its time version identifies the
consecutive terminal time jets. This proves compatibility to all orders.

For a terminal point with z!=0, q>=|z|^(1/D)>0, so that argument applies.
For z=0 but r>0, use the exterior. Differentiation under the positive
integral in (3) gives

\[
\mathcal H^{(j)}(Z)=\frac{(-1)^j(h)_j}{\Gamma(1+h)}
 \int_0^\infty e^{-v}v^{h+j}(1+Zv)^{-h-j}\,dv,
\quad |\mathcal H^{(j)}(Z)|\le(h)_j(1+h)_j\quad(Z\ge0).\tag{6}
\]

Hence |partial_tau^j K|<=C_j r^(-1-2h-2j), and

\[
|\partial_\tau^j(K(\rho,\tau)^2/\rho)|
       \le C_j\rho^{-3-4h-2j}.\tag{7}
\]

This pressure majorant is integrable for rho>=r_min>0. Dominated
convergence supplies its terminal time derivatives, while radial
derivatives follow by the lower endpoint and the chain rule. Thus K and
p0 have compatible smooth one-sided limits on positive-radius compact
sets. A=0 there. These two regions cover every compact terminal set
avoiding x=0 and prove (4) and (5).

This resolves the conditional endpoint step separated in the previous
checkpoint. It uses its already uniform physical estimates. Independent
review of the common charts, bounded overlap and global compatibility
remains a foundation requirement rather than a new boundary assumption.

## 3. Choose cutoffs inside the SAME local domain

Either q<=2tau, or (2) implies q<=(2z^2)^(1/(1-2h)). Since
1/(1-2h)<2 for every allowed h,

\[
q\le2\tau+4|z|^{1/D}.\tag{8}
\]

Choose tau0=qstar/16 and z0=(qstar/16)^D. On tau<=tau0,
|z|<=z0, equation (8) gives q<=3qstar/8<qstar/2, leaving a fixed
margin from the outer local boundary. Choose any sufficiently small r0>0
and smooth cutoffs satisfying

\[
\chi_x=1\ (r\le r_0/2,\ |z|\le z_0/2),\qquad
\operatorname{supp}\chi_x\Subset\{r<r_0,|z|<z_0\},
\]
\[
\chi_t=0\ (\tau\ge\tau_0),\qquad
\chi_t=1\ (\tau\le\tau_0/2),\qquad c=\chi_x(r^2,z)\chi_t(t).
\tag{9}
\]

The dependence on r^2 is smooth in Cartesian coordinates. Define on
t<1, with smooth zero extension outside the supported local domain,

\[
u=\operatorname{curl}(c\mathcal A)+c\mathcal B,
\qquad p=cp_0,\qquad d=\nabla c\times\mathcal A,
\qquad u=cu_0+d.\tag{10}
\]

Curl has zero divergence, and div(cB)=r^-1 partial_theta(cb)=0 because
c is axisymmetric. Therefore div u=0 exactly. The Cartesian base factors
and smooth annular zero extensions justify this also at r=0. A general
theta-dependent cutoff would spoil the direct swirl divergence and is
not used. The fixed domain margin justifies smooth extension elsewhere.

The fields are compactly supported in K=supp chi_x for every t<1 and
vanish near t=0, so their initial velocity is identically zero. They equal
the COMPLETE local fields near (0,1), including their potentials and pressure.

## 4. Full localization error, with no omitted terms

For clarity display arbitrary positive nu in the algebra; the inherited
construction and flatness implication use nu=1. Let
R_nu(v,p)=v_t+(v.grad)v-nu Delta v+grad p. Direct expansion of (10) gives

\[
\begin{aligned}
R_\nu(u,p)={}&cR_\nu(u_0,p_0)+c_tu_0+d_t
 -\nu(\Delta c)u_0-2\nu\sum_i(\partial_i c)\partial_i u_0-\nu\Delta d
 +p_0\nabla c\\
&+c(c-1)(u_0\cdot\nabla)u_0
 +c(u_0\cdot\nabla c)u_0\\
&+c(u_0\cdot\nabla)d+c(d\cdot\nabla)u_0+(d\cdot\nabla)d.
\end{aligned}\tag{11}
\]

The only omitted-looking transport term is (d.grad c)u0, which is
IDENTICALLY zero because d=grad c cross A. Every displayed term is
retained. In particular d_t includes the mixed time/spatial derivatives
of c, Delta d includes its third spatial derivatives, and the pressure
cutoff contributes p0 grad c. All covariances and mean memory are already
inside u0 and p0; (11) acts on the same realized field.

For an m-th force jet this expression requires A through m+3, B through
m+2 and p0 through m+1. Its next time derivative requires one more.
Thus the all-order interface bounds in Section 2 suffice on every compact
set avoiding x=0. In regions where c is zero, every term vanishes; near
(0,1), c=1 and d=0, so f=R_1(u,p) is exactly the inherited flat local
residual. No cutoff estimate has to overpower the singular core: the
cutoff derivatives vanish there.

The force need not be divergence-free. Taking divergence of its exact
definition gives

\[
-\Delta p=\sum_{i,j}\partial_i\partial_j(u_i u_j)-\operatorname{div}f.
\tag{12}
\]

This includes the force-divergence term; there is no separate unjustified
unforced pressure projection or replacement of the pressure in (10).

## 5. Uniform terminal force jets and a smooth compact continuation

Define f=R_1(u,p) for t<1. Away from x=0, Sections 2 and 4 give compatible
one-sided limits of every mixed derivative. Near (0,1), parent flatness
and (8) give, for every alpha,j,N,

\[
|\partial_x^\alpha\partial_t^j f(x,t)|
       \le C_{\alpha,j,N}(\tau+|z|^{1/D})^N.\tag{13}
\]

To get UNIFORM convergence on R^3, first choose a small ball where
(13) is small at late times. On the remainder of K use a finite compact
cover and the Cauchy bounds (5). The derivatives vanish off K. Their
limits are therefore uniform everywhere. Passing to the limit in spatial
coordinate integrals identifies smooth compact functions F_j with

\[
\lim_{t\uparrow1}\partial_x^\alpha\partial_t^j f(x,t)
       =\partial_x^\alpha F_j(x),\quad
\operatorname{supp}F_j\subset K,\quad
\partial_x^\alpha F_j(0)=0.\tag{14}
\]

Time compatibility follows from
partial_t^j f(t)=F_j-integral_t^1 partial_t^(j+1) f(s) ds, also after
every spatial derivative. In general F_j is NOT zero away from the
origin, so setting f=0 immediately after t=1 would be wrong.

For sigma=t-1>=0 take chi0=1 near zero, supported in (-infinity,1),
and define

\[
f(x,1+\sigma)=\sum_{j\ge0}
 \chi_0(b_j\sigma)\frac{\sigma^j}{j!}F_j(x),\qquad b_0=2.
\tag{15}
\]

Only nonnegative sigma is used; chi0 can be the same smooth plateau
cutoff as in the preceding sum. For |alpha|=a and m time derivatives,
Leibniz' rule gives the EXACT finite sum

\[
\partial_\sigma^m[\chi_0(b_j\sigma)\sigma^j/j!]
=\sum_{\substack{0\le k\le m\\m-k\le j}}
 \binom{m}{k} b_j^k\chi_0^{(k)}(b_j\sigma)
 \frac{\sigma^{j-m+k}}{(j-m+k)!}.
\]

On the support 0<=sigma<=1/b_j, every summand is bounded by a
constant times b_j^(m-j). Thus

\[
\|\partial_x^\alpha\partial_\sigma^m T_j\|_\infty
 \le C_{j,m}b_j^{m-j}\|F_j\|_{C^a}.\tag{16}
\]

For j>=1 choose increasing dyadic integers b_j=2^n_j so that the right
side is <=2^-j for every a+m<=floor(j/2). This is a finite list with
m-j<0, hence possible for arbitrarily large terminal jet constants.
`borel_log2_scale` makes this choice in logarithmic coordinates. Its
synthetic tests do not evaluate the actual F_j or their norms.

Every fixed mixed differentiated tail now converges uniformly, since its
requirements hold for all j>=max(1,2(a+m)). At sigma=0 the cutoff is
constant in every fixed summand, and its m-th derivative contributes only
when j=m, giving F_m. Uniform differentiated convergence and (14) prove
smooth matching through t=1. All summands keep spatial support in K and
vanish for sigma>=1/2, since b_j>=2. They are flat at x=0. The uniform
convergence also proves smoothness at the support boundaries.

Consequently the EXACT preterminal residual extends to

\[
f\in C_c^\infty(\mathbb R^3\times(0,\infty);\mathbb R^3),\quad
\operatorname{supp}f\subset K\times[1-\tau_0,3/2].\tag{17}
\]

Every fixed derivative is bounded; compact support then gives decay of
every polynomial order in space and time. This is an existence proof
conditional on the preceding continuum inputs. It supplies neither an
effective evaluation of the enormous constants nor a smooth velocity
after t=1. It does not assume convergence of the raw terminal Taylor series.

## 6. Whole-space energy, dissipation and surviving lower bound

On each [0,T], T<1, u and p are smooth with fixed compact spatial support.
Multiply their exact equation by u and integrate. Divergence and the
support make pressure and transport fluxes vanish, giving

\[
\tfrac12\frac{d}{dt}\|u(t)\|_2^2+\nu\|\nabla u(t)\|_2^2
       =\langle f(t),u(t)\rangle.\tag{18}
\]

With nu=1 set F(t)=integral_0^t ||f(s)||_2 ds. To avoid dividing by a
possibly zero norm, differentiate sqrt(||u||_2^2+delta^2), discard
dissipation and use Cauchy--Schwarz. Its derivative is <=||f||_2.
The zero datum gives ||u(t)||_2<=F(t) after delta->0. Inserting this
in the integrated identity yields

\[
\|u(t)\|_2^2+2\int_0^t\|\nabla u(s)\|_2^2ds
 \le2\int_0^t F'(s)F(s)ds=F(t)^2.\tag{19}
\]

The compact smooth force has F(1)<infinity. Thus the kinetic energy is
uniformly bounded and monotone convergence in the dissipation integral
gives finite total dissipation on [0,1). No regularity or continuation of
u itself at t=1 is used. This improves the preceding local energy bound
to a conditional WHOLE-SPACE energy and dissipation bound.

Finally x_tau=(sqrt(2 Xin tau),0,0) lies inside the unchanged cutoff
plateau whenever tau<=min(tau0/2,r0^2/(8 Xin)). Hence for sufficiently
small tau,

\[
|u_\theta(x_\tau,1-\tau)|
       \ge(e_0/2)\tau^{-A}.\tag{20}
\]

This lower bound includes all realized waves, means and physical cutoffs,
conditional on the inherited positive core coefficient. It is compatible
with bounded energy: its high velocity occupies a shrinking region.

## 7. What is checked and the next unclosed argument

The independent diagnostics cover all twelve rows of (11) at three
positive viscosities, cylindrical versus Cartesian residuals, nonzero
axis values, exact support/initial plateaus, the gamma-integral endpoint,
the mixed jets in (15), extreme synthetic scale constraints, and the
forced whole-space energy identity for a smooth compact manufactured curl.
Controls detect missing curl terms, theta-dependent localization of the
direct swirl, reversed tau time sign, an invalid zero force continuation,
and flatness at the origin without endpoint control elsewhere.

The energy quadrature at 64 nodes failed its integration-by-parts check
(about 3.67e-7). It is retained in evidence. Refinement to 96 nodes gives
about 1.12e-10 and passes the SAME 2e-9 tolerance; no tolerance was relaxed.
These are algebra/integration checks, not interval certification of the
actual global constructed flow.

```sh
python3 research/navier_stokes_cascade/results/localized_force_extension/test_audit.py -v
python3 research/navier_stokes_cascade/results/localized_force_extension/audit.py \
  --verify-record research/navier_stokes_cascade/results/localized_force_extension/evidence.json \
  --output /tmp/localized_force_extension_recheck.json
```

Next: audit the comparison with an arbitrary putative global smooth
bounded-energy solution having this force and zero datum. Its pressure
gradient cannot simply be assigned a decay condition; it must be recovered
from the equation, and every flux at spatial infinity must be controlled.
Only then can the lower bound contradict that solution's smoothness.
Independent review of the profile, background, uniform all-order physical
bounds and global support geometry remains a parallel prerequisite.
The current result is a conditional localization/forcing transfer theorem,
not an independently verified resolution of the Navier--Stokes problem.
