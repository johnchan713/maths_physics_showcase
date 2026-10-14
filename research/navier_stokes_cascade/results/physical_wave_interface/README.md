# Physical wave, covariance and full residual interface

This checkpoint derives exact physical curl and covariance identities,
the full cylindrical linearized momentum operator, the remaining mean
equations, and a conditional pulse propagator lemma. It checks them with
compact physical waves and a separate moving-normal reference pulse.
It does not construct the actual band waves against the summed background.
Uniform actual pulse bounds, mean corrections and full residual improvement
remain open. No forcing, complete-field energy or blow-up result follows.

Parent: `25f59628da5bc7c2efb66816075ddd98d7e340c1`. The preceding
[background summation](../background_summation/README.md) is conditional on
the inherited leading construction. Its augmented error is flat, but its
bare momentum residual retains the signed annular stress divergence.
That distinction is maintained here.

The primary source is the supplied
[proposed construction](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
equations (7.2)–(7.16), (7.24)–(7.41), (8.3), and (9.1)–(9.2).
Reviewed PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The arguments below are a scoped review, not independent verification of
the manuscript's complete construction.

## 1. Exact real curl and all six covariance entries

Work away from the axis. Let phi be real, n=grad(phi) nonzero, and a be a
real vector field with n.a=0. Define

\[
c=\frac{n\times a}{|n|^2},\qquad B=\nabla\times c,
\qquad A_\kappa=-\frac{c}{\kappa}\sin(\kappa\phi).
\]

The triple-product identity n cross (n cross a)=-|n|^2 a gives exactly

\[
w_\kappa=\nabla\times A_\kappa
 =a\cos(\kappa\phi)-\kappa^{-1}B\sin(\kappa\phi),
\qquad \nabla\cdot w_\kappa=0.                         \tag{1}
\]

No leading-order approximation is used. The cylindrical components of B
include

\[
B_r=r^{-1}\partial_\theta c_z-\partial_zc_\theta,
\quad B_\theta=\partial_zc_r-\partial_rc_z,
\quad B_z=(\partial_r+r^{-1})c_\theta-r^{-1}\partial_\theta c_r.
\]

The r^-1 c_theta term cannot be dropped. Smooth compact vector potentials
on an annulus extend by zero smoothly across the axis and outside their
support. No coordinate singularity is introduced by that extension.

Suppose a and c have theta-independent cylindrical components and
phi=p theta+psi(r,z,t), with kappa p a nonzero integer. Angular averaging
of cylindrical tensor components then gives the entire symmetric tensor

\[
Q_\kappa=\langle w_\kappa\otimes w_\kappa\rangle_\theta
 =\tfrac12a\otimes a+\frac{1}{2\kappa^2}B\otimes B.       \tag{2}
\]

The cosine/sine mixed average is zero. Thus the primary real curl changes
the mean covariance **quadratically**, with positive semidefinite error;
there is no term linear in kappa^-1. For a fixed derivative order j and
fixed theta-independent amplitudes,

\[
\|D^j(Q_\kappa-\tfrac12a\otimes a)\|
 \le C_j\kappa^{-2}\sum_{i\le j}\|D^iB\|\|D^{j-i}B\|.   \tag{3}
\]

This is a useful exact simplification for real primary pulses. If the
amplitude or phase depends on the band or frequency, the B bounds on the
right must still be proved uniformly. Equation (2) does not supply those
bounds. Nor does it remove the pointwise O(kappa^-1) velocity correction,
its derivatives, or interactions with earlier complex harmonics.

For multiple carriers, sum (2) only when cross products vanish for a
specified reason: disjoint supports, or exact angular-frequency
orthogonality. Distinct phases with the same angular frequency need not
have zero mean product. The actual manuscript uses auxiliary support
separation; the manufactured two-frequency fixture below does not verify
that separation or replace it. Nonzero-frequency nonlinear interactions
must still be estimated even when their mean vanishes.

## 2. Full momentum and the mean terms left after stress matching

Write R(u,p)=u_t+(u.grad)u-nu Delta u+grad p, with nu>0. For any
divergence-free w and an axisymmetric background v,

\[
R(v+w,p+\pi)-R(v,p)=L_vw+\nabla\pi+\nabla\cdot(w\otimes w),
\quad L_vw=w_t+(v\cdot\nabla)w+(w\cdot\nabla)v-\nu\Delta w. \tag{4}
\]

Here w is the COMPLETE curl velocity from (1). Let v=(b,V,G) in the
(r,theta,z) frame, F=V/r, and
Delta_s=partial_rr+r^-1 partial_r+r^-2 partial_theta_theta+partial_zz.
Set T=partial_t+b partial_r+F partial_theta+G partial_z. The three
components of L_v w+grad pi are exactly

\[
\begin{aligned}
&Tw_r+b_rw_r+b_zw_z-2Fw_\theta
 -\nu(\Delta_sw_r-r^{-2}w_r-2r^{-2}\partial_\theta w_\theta)+\pi_r,\\
&Tw_\theta+(2F+rF_r)w_r+V_zw_z+(b/r)w_\theta
 -\nu(\Delta_sw_\theta-r^{-2}w_\theta+2r^{-2}\partial_\theta w_r)
 +r^{-1}\pi_\theta,\\
&Tw_z+G_rw_r+G_zw_z-\nu\Delta_sw_z+\pi_z.                 \tag{5}
\end{aligned}
\]

These formulas retain the radial transport, cylindrical connections,
pressure, slow derivatives and viscosity that are not all in the
principal pulse ODE. The principal matrix's two 2F terms are visible in
(5). A principal ODE solution must still control the rest of (5), every
derivative of its curl correction, all cutoffs and the quadratic term in
(4). A covariance identity does not imply that control.

For a symmetric theta-independent cylindrical tensor Q,

\[
\nabla\cdot Q=
\begin{pmatrix}
\partial_rQ_{rr}+(Q_{rr}-Q_{\theta\theta})/r+\partial_zQ_{rz}\\
(\partial_r+2/r)Q_{r\theta}+\partial_zQ_{\theta z}\\
(\partial_r+1/r)Q_{rz}+\partial_zQ_{zz}
\end{pmatrix}.                                         \tag{6}
\]

Suppose the background residual equals E-D_tan S, where
D_tan S=(0,(partial_r+2/r)S_theta,(partial_r+1/r)S_z),
and the principal covariance Q0 matches Q0_rtheta=S_theta,
Q0_rz=S_z. The remaining mean source before mean velocity corrections is

\[
D=\begin{pmatrix}
\partial_rQ^0_{rr}+(Q^0_{rr}-Q^0_{\theta\theta})/r+\partial_zS_z\\
\partial_zQ^0_{\theta z}\\
\partial_zQ^0_{zz}
\end{pmatrix},                                         \tag{7}
\]

in addition to E and div(Q-Q0). For zero-mean oscillations, the angular
mean of the linear terms in (4) is zero. Angular averaging is averaging
of cylindrical components; it is distinct from the auxiliary average
used in the manuscript's mean constraints.

Choosing pi(r,z,t)=integral_r^infinity D_r(s,z,t) ds cancels the radial
entry but leaves D_z+pi_z and D_theta. It can also leave a nonzero
pressure constant inside the wave annulus. Imposing both an outer
pressure datum and an unchanged inner datum requires an integral
compatibility condition, or a mean correction that enforces it. Radial
balance by itself does not enforce the axial equation or these constraints.

The compact fixture gives, at (r,z,t)=(1.3,.2,.47),

| Quantity | Value (rounded) |
|---|---:|
| D before a radial pressure correction | (-.173840, .144236, -.926738) |
| pi_z | .376415 |
| D after that correction | (0, .144236, -.550323) |
| pressure constant inside wave support | -.764082 |

These numbers are a manufactured example of the missing mean terms, not
the actual mean defect of the manuscript's field. Its physical mean
equations and integral corrections must be constructed for that field.

## 3. Positive primary covariance and signed increments

Fix orthonormal tangential K,N with K=(-N_z,N_theta), u>0 and A_c>0.
For sigma=+1,-1 put

\[
d_\sigma=e_r-\sigma uK-A_cN,\qquad
H_\sigma=\tfrac12f^2(-\sigma uK-A_cN).
\]

For a_sigma=sqrt(1-sigma tau), |tau|<1, principal real carriers with
amplitudes f a_sigma d_sigma have cross covariance
S=f^2(-A_c N+u tau K). The two columns H are independent, and a target
with S.N<0 and |S.K|<(u/A_c)(-S.N) has strictly positive squared
amplitudes. This finite-dimensional identity does not prove uniform edge
estimates for the actual corrected stresses or actual pulse columns.

For any signed target increment Sigma, solve H b=Sigma and set
delta a_sigma=b_sigma/(2a_sigma), holding the positive base a fixed.
Writing C(w)=<w_r(w_theta,w_z)> and
B(w,l)=C(w+l)-C(w)-C(l), the leading change is B(W,L)=Sigma.
There is no sign restriction on this linear equation. Its constants
depend on H^-1 and the base amplitudes and require bounds in applications.

Keep the complete old velocity U=W+E and new curl velocity V=L+R.
The exact actual-velocity increment is

\[
C(U+V)-C(U)=\Sigma+B(E,L)+B(U,R)+C(V).                    \tag{8}
\]

In particular the signed linear inverse does not eliminate C(V).
The fixture checks (8) with one negative inverse coefficient, all curl
remainders retained, and a signed target with positive N component.
The omitted quadratic term has norm .533 at one tested point.

## 4. Moving-normal pulse: exact equations and conditional estimate

For a smooth real nonzero normal n(v), put

\[
K_0=\begin{pmatrix}0&-2F&0\\2F+g_\theta&0&0\\g_z&0&0\end{pmatrix},
\quad A_n=-K_0+\frac{n(n^TK_0-(n')^T)}{|n|^2},
\quad d=\nu\epsilon k^2|n|^2.
\]

For nonzero integer m, the constrained principal equation is

\[
t'+K_0t+m^2dt+ikmn\pi=-f,\qquad n\cdot t=0.
\]

Differentiate the constraint and take its normal component. This gives
exactly

\[
t'=A_nt-m^2dt-\operatorname{proj}_{n^\perp}f,
\qquad \pi=\frac{i}{km|n|^2}(n\cdot K_0t-n'\cdot t+n\cdot f). \tag{9}
\]

The moving-normal term and factor i in the pressure are essential. The
equation preserves n.t=0 for compatible initial data. In the homogeneous
case its real energy identity is

\[
\tfrac12(|t|^2)'=-g\cdot[t_r(t_\theta,t_z)]-m^2d|t|^2.    \tag{10}
\]

The 2F terms cancel in this energy identity; viscosity does not disappear.

Here is a conditional propagator argument, with its input isolated. Suppose
a smooth frame B maps R^2 onto n-perp, B and its left inverse are uniformly
bounded, and in its coordinates the undamped operator is
diag(lambda,-lambda)+E, with lambda>=0,
norm(E)+|d-d_ref|<=C/L. Let L be comparable to the pulse length and
P'/P=lambda-d_ref. Then for m=1 the Euclidean norm inequality is
partial_v |z| <= (lambda-d_ref+C'/L)|z|. Integration proves

\[
\|V_1(v,w)\|\le e^{C'}P(v)/P(w),\quad 0\le w\le v\le L.
\]

Scalar damping commutes with the time-dependent undamped matrix, so

\[
V_m(v,w)=\exp\!\left(-(m^2-1)\int_w^v d(a)\,da\right)V_1(v,w). \tag{11}
\]

For d>=0 this proves the same forward bound for every nonzero integer m.
A zero-initial-data source bounded by C P(w) has response bounded by
C' L P(v) by Duhamel. At each fixed derivative order, polynomial bounds
on the coefficient/frame derivatives yield polynomial losses in L by
induction on the differentiated ODE. The slow/edge weights must be fixed
along the path and the coordinate changes must have corresponding bounds.
This last list is a hypothesis, not a consequence of (11).
Zero source on an entire parameter rectangle preserves support by
uniqueness; periodic/descent properties require the coefficients and
source to descend to the same torus. Smooth dependence alone does not
establish the actual support geometry.

For the reference normal n=B_s(s,K), s=sigma u(.5+v/L),
lambda=lambda0/sqrt(1+s^2), choose

\[
B_s^2=\frac{\lambda_0}{\nu\epsilon k^2(1+u^2)^{3/2}}.      \tag{12}
\]

This is B_s SQUARED. Let rho=u(.5+v/L). The envelope normalized to one
at v=L/2 has the closed form

\[
\log P=\frac{\lambda_0L}{u}\left[
\operatorname{asinh}\rho-\operatorname{asinh}u
-\frac{\rho+\rho^3/3-u-u^3/3}{(1+u^2)^{3/2}}\right].     \tag{13}
\]

Since
-(lambda-d_ref)'=(u lambda0/L)rho[(1+rho^2)^(-3/2)+2(1+u^2)^(-3/2)],
positive constants c,C follow by taking rho between u/2 and 3u/2.
Twice integrating from the midpoint gives the continuum bound

\[
e^{-C(v-L/2)^2/(2L)}\le P(v)\le e^{-c(v-L/2)^2/(2L)}\le1. \tag{14}
\]

On fixed relative endpoint windows this is exponential in -L. After a
cutoff psi is inserted before the curl, the homogeneous principal
coefficient error is psi' t; the curl and full residual still have their
additional derivatives. Exponential endpoint decay absorbs fixed
polynomial losses if L is comparable to the square of the band index.
The actual comparison constants and derivative bounds remain to be
proved before applying that observation to the summed background.

### Reference ODE diagnostics, distinct from the physical fixture

Use F=1, g=3N, N=(-1,1)/sqrt(2), u=1.5, nu epsilon k^2=1, k=16.
Then lambda0^2=-2F N_theta(2F N_theta+3)>0 and
c0=lambda0/(2F N_theta)<0. Integrate **y=t/P**:
y'=(A_n-lambda I)y, with y(0)=e_r-s(0)K+c0 sqrt(1+s(0)^2)N.
This retains both the moving normal and viscosity in t=P y, while
avoiding exponentially small initial data being accepted as zero by a
solver's absolute tolerance. It does not assume y follows the instantaneous
growing eigenvector exactly.

For the compact pulse cutoff, compute
h_sigma=(2L)^-1 integral psi^2 P^2 y_r^2 and
H_sigma=(2L)^-1 integral psi^2 P^2 y_r y_tan.

| L | sqrt(L) h (either sign) | max norm(H/h+A_c N+sigma uK) | sampled max norm(psi' t) |
|---|---:|---:|---:|
| 40 | .327583 | .0212253 | .0249530 |
| 100 | .325715 | .00854221 | .000128000 |
| 240 | .324935 | .00355762 | 3.924e-9 |

Here A_c=-c0 sqrt(1+u^2), different numerically from the physical fixture's
A_c. Positive primary weights and a signed linear inverse are checked
using these INTEGRATED columns, rather than substituting reference limiting
directions into the inverse. A second, tighter ODE solve tests numerical
stability. Samples and solver tolerances are diagnostics, not interval
enclosures or proofs of actual band-uniform estimates.

## 5. Compact physical fixture and a useful failure of kinematic matching

For the independent physical-coordinate checks use p=1, beta=-1, u=1.5,
c0=-.75, A_c=-c0 sqrt(1+u^2), and set

\[
a_t=(r^{-2}+1)^{1/2},\quad K=(1/(ra_t),-1/a_t),
\quad N=(K_z,-K_\theta),
\]
\[
H(r)=\sqrt{1+r^2}+\log\frac r{1+\sqrt{1+r^2}},\quad H'=a_t,
\quad\phi_\sigma=\theta-z+\sigma uH(r)-.37t.
\]

Thus n_sigma=(sigma u a_t,1/r,-1) is an EXACT phase gradient and
n_sigma.d_sigma=0. The amplitudes are f sqrt(1-sigma tau)d_sigma,
tau=.25+.1 sin z+.05 cos t, with the standard C-infinity bump b,

\[
f=e^3 b((r-1.4)/.55)b(z/1.2)b((t-.5)/.45),\quad
b(x)=\begin{cases}e^{-1/(1-x^2)},&|x|<1,\\0,&|x|\ge1.\end{cases}
\]

Use angular frequencies k and 2k for the two signs and take the exact
curls. Support is in .85<r<1.95, |z|<1.2, .05<t<.95. The background
fixture is the divergence-free field
(-.17r,r Omega,.34z), Omega=.43+.31r^2/2+.19z+.11t,
with a nonconstant manufactured pressure. A supported oscillatory pressure
is also included, so pressure terms are exercised in the residual.

Third-order Taylor arithmetic differentiates each vector potential.
Its curl is available through second derivatives, sufficient for physical
viscosity. Independent Cartesian curls and Cartesian momentum are compared
with the cylindrical formulas, through the required derivative orders.
No finite differences or numerical PDE evolution enter this check.
Angular trapezoidal averaging with 8k+1 nodes is exact in exact arithmetic
for the finite Fourier products being sampled.

The principal cross covariance is the same S for all k, while the full
residual increment at (r,z,t)=(1.3,.2,.47) behaves as follows:

| k | angular RMS full residual increment | RMS divided by k^2 |
|---|---:|---:|
| 8 | 2344.09 | 36.6265 |
| 16 | 9233.14 | 36.0669 |
| 32 | 36820.2 | 35.9572 |
| 64 | 147173 | 35.9310 |

There is also an analytic reason for this failure, rather than just four
samples. With frequency-independent amplitudes and phases, -nu Delta w
has leading term nu sum_sigma kappa_sigma^2 |n_sigma|^2 a_sigma
cos(kappa_sigma phi_sigma). Every other term in (4), including the
Laplacian of the curl remainder, is at most O(k) at a fixed interior
point. The two distinct angular frequencies are orthogonal. Therefore
the RMS increment divided by k^2 tends to the strictly positive value

\[
\nu\left[\tfrac12\sum_{\sigma=+,-}
 j_\sigma^4|n_\sigma|^4|a_\sigma|^2\right]^{1/2},
\qquad (j_+,j_-)=(1,2).                                 \tag{15}
\]

This proves that kinematic matching and exact incompressibility ALONE
do not improve the full residual. It is not an obstruction to a correctly
evolved shear–viscosity pulse. The manuscript's actual pulse must cancel
the leading dynamics and then bound the rest of (5). The manufactured
physical field here does not satisfy that pulse equation, while the
separate reference ODE is not the actual global physical wave field.

## 6. Reproduction, scope and next acceptance condition

Run from the repository root:

```sh
python3 research/navier_stokes_cascade/results/physical_wave_interface/test_audit.py -v
python3 research/navier_stokes_cascade/results/physical_wave_interface/audit.py \
  --verify-record research/navier_stokes_cascade/results/physical_wave_interface/evidence.json \
  --output /tmp/physical_wave_interface_recheck.json
```

NumPy and SciPy are required. The output path must be fresh. Evidence
freezes the protocol, source hashes, numerical observations and scope.
All preceding audit source bytes and GOAL.md are checked and preserved.
Controls deliberately omit the curl correction, vector-Laplacian and
advection connections, pressure, moving-normal term, pressure's factor i,
viscous damping, support/frequency justification for mean cross terms,
or the signed inverse's quadratic remainder. Those omissions must fail.

| Kind of result | What is supplied | What is not supplied |
|---|---|---|
| Exact identities | (1)–(10), scalar damping identity (11), reference envelope and its Gaussian bounds | Uniform constants for the actual band fields |
| Conditional estimate | Forward propagator and differentiated Duhamel bounds given frame, phase, support and coefficient hypotheses | Verification of those hypotheses on every actual band and edge |
| Physical diagnostics | Compact real curls, complete tensor, pressure and full NS residual in two coordinate systems | The actual summed-background wave solution |
| Reference ODE diagnostics | Moving normal, pressure, viscosity, integrated covariance and signed inverse | Uniform actual pulse/edge estimates or auxiliary-torus descent |
| Actual construction | The inherited background argument remains conditional | Actual mean corrections, full residual improvement, smooth forcing, energy or blow-up |

The next decisive step is to verify the phase/frame comparison and
weighted pulse derivative bounds for the ACTUAL summed background,
with support/descent and fixed discrete frequencies retained. Then use
the resulting complete Q in the physical mean system, enforcing the
integral pressure constraints rather than canceling only two entries.
Acceptance requires a positive gain for the full residual in fixed
physical derivative norms, retaining every term in (4)–(8). Until that
estimate is proved, a successful covariance inverse or reference pulse
is progress on the interface, not a completed physical wave correction.
