# Compact mean foundation and the sharp axial conditioning loss

Parent: `372574a3e9f08bba9d31c351ce04298cd709a0ab`.
This independently audits the compact mean maps still retained as analytical
dependencies by the [primary-pulse foundation](../physical_pulse_foundation/README.md).
It supplies their all-fixed-order, common-domain estimates and quantifies
which part of the five-map deteriorates as lambda tends to zero. It does
not repeat or promote the older conditional finite residual iteration.
Historical sources, evidence, fixed positive lambda, h, N, `C^-20`, supports,
carriers and the background's complete pressure memory are unchanged.

The [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`,
(8.3)--(8.27), is the object being reviewed, not a theorem used to assume
these estimates. The exact-core and remaining finite background/pulse
existence and gap dependencies remain. The new assertions below are
analytical estimates; exact finite fixtures check identities, not continuum
suprema or the actual enormous parameter values.

## 1. Uniform inverse estimates, rather than high-norm smallness

Use the actual inputs now supplied by the parent. On every accepted band
Q, epsilon=Q^h, S=ell^2, the enlarged physical chart has 1/4<=s=q/Q<=4.
All fixed jets of s,X,eta are bounded on the same domain, one-sided at
T=0. Write L=vr.partial_y, Nt=vt.partial_y, vr=(1,-beta),
vt=(beta,1), beta=sqrt(2)-1. The radial operator is

\[
D_r=\partial_R+M d_rR^{d_r-1}L,\quad D_e=D_r+e/R,
\quad D_z=\epsilon\partial_Z,\quad t_*=-\epsilon\partial_T+c_iN_t.
\tag{1}
\]

All shells stay in a fixed positive compact radial interval. From the
EXACT definitions, with rho_g=log(Lambda_g)/log(Tg),

\[
M=c_i^{\rho_g}\epsilon^{-\kappa},\quad
(TgS)^{-1}<c_i\le S^{-1},\quad
M^{-1}<(TgS)^{\rho_g}\epsilon^\kappa. \tag{2}
\]

Both sides of the flooring estimate matter: an upper bound for M alone
would not justify flatness of a nonstationary remainder.

For e=0,1,2, define T_e and A_e by the SAME shifted primitive and compact
cutoff (4) of [compact means](../compact_mean_corrections/README.md).
Put U=R^dr and F(U)=R^e f(R)/(dr R^(dr-1)), extended by zero. Then

\[
R^eT_ef=(1-\chi_m)\int_{-\infty}^0F(U+a,y+Ma v_r)\,da
             -\chi_m\int_0^\infty F(U+a,y+Ma v_r)\,da. \tag{3}
\]

The shift is independent of U,Z,T at fixed a. On the fixed radial hull,
Leibniz differentiation therefore bounds each order-m coefficient seminorm
by finitely many order-at-most-m input and cutoff seminorms, with no M
power. Integrals have bounded effective length. Coordinate conversion
between R and U has bounded fixed jets. Near each shell edge only the
primitive toward that edge occurs; the fixed flat edge weight times any
fixed inverse-distance power is increasing sufficiently near that edge.
Its integral is bounded by the same type of edge weight, with constants
depending on the fixed order. The remaining interior has a positive edge
weight lower bound. This proves T_e:M_alpha -> M_alpha for every fixed
order. It preserves slow projection/support and smooth zero extension,
but may fill the auxiliary torus. Direct differentiation gives
D_e T_e f=f-A_e f; the unshifted primitive would not satisfy this identity.

For nonzero k in Z^2, each divisor is A+B sqrt(2), with conjugate of
size at most (1+sqrt(2))|k|_1. Its norm A^2-2B^2 is a nonzero integer.
Thus |vr.k| and |vt.k| are at least [(1+sqrt(2))|k|_1]^-1.
If int R^e bar(f)=0, the zero Fourier mode of A_e is EXACTLY zero.
Integrating a nonzero Fourier mode p times by parts on the whole U line
gives, for output order m,

\[
\|A_ef\|_m\le C_{m,p}M^{-p}\|f\|_{m+2p+3}
 \le C_{m,p}\epsilon^{\alpha+p\kappa}S^{B_{m,p}}.
\tag{4}
\]

The displayed derivative count is sufficient, not optimal: p radial
derivatives, p divisor losses and three extra torus derivatives make the
two-dimensional Fourier tail summable. Edge weights are harmless because
A_e is supported strictly in the interior. For any fixed physical loss
epsilon^-elloss and desired epsilon^K choose p with
alpha+p*kappa-elloss>K. The strict margin absorbs every polynomial in S
on the SAME starting-band domain; it changes constants, not the domain.
There is no high-norm smallness or derivative-dependent first band.
A nonzero radial zero-mode moment is a genuine obstruction, not flat.

The zero-Haar temporal multiplier (2 pi i c_i vt.k)^-1 similarly obeys
||u||_m<=C_m S ||f||_(m+4), preserves reality/slow support, and is defined
only when bar(f)=0. It may fill the torus. For axial desired gamma_d of
zero weighted moment, Psi=T_1 gamma_d, beta=-epsilon partial_Z Psi and
gamma=D_1 Psi give exact divergence zero and gamma=gamma_d-A_1 gamma_d.
The actual A_1 velocity remains in all later products. For slow inputs
with zero radial moment A_1 is identically zero, not merely flat.

These are single physical maps, not unrelated chart choices: the radial
shift is Y+[(r')^dr-r^dr]vr before rescaling; chi_m is a function of
r^2/(2q). For the integer covering Pi_i f=f(Jg^i Y),
Nt,abs Pi_i=Tg^i Pi_i Nt,i. Haar is preserved; hence its zero-mean inverse
intertwines with Tg^-i, and physical normalization gives exactly
Q^(-1-h)Tg^-i=c_i^-1. No averaging of distinct preimages is substituted.

## 2. Five fixed profiles: only the axial block loses 1/lambda

On the full convex hull of the reserved I4 profiles the actual summed base
has G_q=0 and V_q=a(eta)x^(-1-2lambda), x=r/sqrt(q), a>=a0>0.
Its radial flow is NOT zero. Choose the same narrow nonnegative bump and
three disjoint geometric dilates with g=exp(dpatch)>1, as in the existing
five-map. Put mu_p=int x^p eta0(x)>0. After dividing rows by their
nonzero factors, the angular and axial blocks are the Vandermonde matrices

\[
V_\theta=\big[1,z_i,z_i^2\big]_{i=1}^3,
\quad z=(g^2,g^{-2}q_\lambda,q_\lambda),\qquad
V_z=\begin{pmatrix}1&g\\1&gq_\lambda\end{pmatrix},
\quad q_\lambda=e^{-2d_{patch}\lambda}. \tag{5}
\]

Row factors are (mu_2,2a mu_(-2-2lambda),-a mu_(-2lambda)) and
(mu_1,a mu_(1-2lambda)). All moments have positive lower bounds and
bounded lambda derivatives on any fixed compact lambda interval, since
the original bump has fixed positive compact support. For
0<=lambda<=lambda_bar=1/4, the THREE angular nodes remain separated,
including at zero. Its inverse is uniformly bounded after row factors
are controlled. There is no angular 1/lambda singularity.

The axial zero-flux row forces s0=-g*s1, while its target row gives

\[
s_1={J_{\theta,q}\over a\mu_{1-2\lambda}g(1-q_\lambda)},
\quad s_0=-g s_1,\qquad
2d_{patch}\lambda e^{-2d_{patch}\lambda_{bar}}
 \le1-q_\lambda\le2d_{patch}\lambda. \tag{6}
\]

Thus every fixed S_alpha target has slow increment bounds

\[
\|\Delta v\|_m\le C_m\|(P_q,J_{z,q})\|_m,\qquad
\|\gamma_d\|_m\le {C_m\over\lambda}\|J_{\theta,q}\|_m. \tag{7}
\]

Here target seminorms include bounded slow derivatives of a^-1 and the
normalization factors. Lambda is fixed, NOT differentiated or retuned.
Constants may depend on the actual fixed patch and background. Formula
(7)'s uniform-in-lambda interpretation is for fixed geometric profiles
and a uniformly controlled positive a; it does not assume the entire
huge construction is uniform in lambda.

The loss is unavoidable, not a poor choice of finite matrix. If the zero
flux constraint int x gamma_d=0 holds on [l,u], choose any x0 there. Then

\[
|J_{\theta,q}|=|a\int x(x^{-2\lambda}-x_0^{-2\lambda})\gamma_d dx|
 \le 2|a|\lambda C_{l,u,\lambda_{bar}}\int x|\gamma_d|dx. \tag{8}
\]

For example take C=max(1,l^(-2lambda_bar),u^(-2lambda_bar))*log(u/l):
apply the mean-value bound to exp(-2lambda*y) in y=log(x). At
lambda=0 a nonzero J_theta cannot be corrected while preserving axial
flux. For the actual fixed lambda=e^(-4T)>0 the bound is finite, though
enormous. No small-target condition is needed for this signed LINEAR map.
It does not make the recomputed nonlinear defects identically zero.

## 3. Physical normalization, support and complete pressure moments

Use the exact q-normalized targets of the older five-map. Relative to
Q-normalized targets their factors are s^(2A), s^(2A-3/2), s^(2A-1)
for P,J_theta,J_z respectively, A=1/2+h. The normalized swirl base in
the R chart is a(eta)s^(lambda-h)R^(-1-2lambda), NOT simply
a(eta)R^(-1-2lambda). The physical increment is q^-A times the fixed
x-profiles; all those s powers and every fixed derivative are bounded.
Multiplying each row and its target by the same factor proves the exact
Q-chart five identities without a Q loss. The slow axial primitive is
exact. Its potential and induced radial velocity can fill the intervals
between bumps; the entire convex hull is inside I4. Therefore
(Delta v,Delta gamma) in M_alpha and Delta beta in M_(alpha+1) on the
same domain, with finite fixed-lambda constants. Actual b=O(epsilon)
and both pressure tails persist on this hull.

For ANY current mean velocities and the entire complete-wave tensor W,
use the full conservative equations of compact means, not a swirl-only
model. Reconstruct g_r, P=int bar(g_r), p_m=T_0(g_r-rho P) at EACH state.
The physical unit-integral bump becomes
rho=s^-1/2 rhohat(R/sqrt(s)); c_rho=(1/2)int R^2 rho. Then exactly

\[
\int R\bar p_m=-\tfrac12\int R^2\bar g_r+c_\rho P,\quad
\int R^2\bar E_\theta=\epsilon\partial_ZJ_\theta,\quad
\int R\bar E_z=\epsilon\partial_Z(J_z+c_\rho P). \tag{9}
\]

The last two use the two zero velocity moments; without them retain
-epsilon partial_T M-epsilon^3 partial_ZZ M. Radial fluxes and cylindrical
viscous terms vanish by compact radial integration by parts. Formula (9)
keeps W_theta z,W_zz,W_rz, axial velocity, radial flow, and the WHOLE moving
product c_rho P. The A0 remainder has zero Haar mean and is retained in
pointwise equations; Haar pressure is compact, giving the first identity.

For the initial zero mean, independent integration gives

\[
P=\int(\bar W_{\theta\theta}-\bar W_{rr})/R
 -\epsilon\partial_Z\int\bar W_{rz},\quad
J_\theta=\int R^2\bar W_{\theta z},
\]
\[
J_z=\int R[\bar W_{zz}-(\bar W_{rr}+\bar W_{\theta\theta})/2]
 +\tfrac\epsilon2\partial_Z\int R^2\bar W_{rz}. \tag{10}
\]

These are S_1 inputs from the actual primary bounds, without a new graph
loss after Haar averaging. The pressure subtraction uses the SAME full
p_B, including Pi1=K1-mu^2/(4X) and Pi2=K2; pressure linearity cancels
these against that exact base, not against leading p0. A radial constant
Pi2 generally has a physical axial gradient. No pressure tail is erased.

Equations (1)--(10) now supply the compact/temporal/five-map assumptions
of the older initial mean update relative to the parent inputs, including
all fixed jets and chart consistency. Its already derived S_1 -> S_2
initial compatibility gain and M_(2-3kappa) angular-mean residual remain
scoped to a fixed primary and recomputed means, with flat summands retained.
They do NOT close the nonzero-angular residual or all later source inputs.

## 4. Checks and next obligation

The exact rational controls exercise signed matrices with nonunit positive
row moments, the axial singular limit, integer quadratic norms, shifted
primitive jets, and complete pressure/covariance identities. Polynomial
bump fixtures have finite endpoint regularity only. Their amplitudes are
manufactured tensor data, not the actual pulse or a PDE solution; actual
complete-curl diagnostics are replayed through the historical chain.
Mutation tests prohibit replacing this scope by unconditional residual or
blow-up claims and bind the source bytes, protocol and inherited evidence.

Next independently verify every generated nonlinear harmonic source's
coefficient class, Gaussian envelope and FULL lifted endpoint extension,
then the complete complex covariance/mean iteration and residual including
viscosity, cutoffs and every old/new interaction. Finite tests cannot
discharge these continuum obligations. No unconditional all-family pulse
bound, full residual improvement, forcing result or blow-up is established.

Reproduce with fresh output (historical evidence must not be overwritten):

```sh
python3 research/navier_stokes_cascade/results/compact_mean_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/compact_mean_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/compact_mean_foundation/evidence.json \
  --output /tmp/compact_mean_foundation_replay.json
```
