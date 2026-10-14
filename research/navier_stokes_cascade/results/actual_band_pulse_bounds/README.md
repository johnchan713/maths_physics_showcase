# Conditional actual-band pulse and primary covariance bounds

This checkpoint instantiates the preceding wave-interface lemma using the
FORMULAS for the actual summed background, rather than assuming its phase
and frame comparisons. Conditional on the inherited leading-profile and
background bounds, it derives a common small-q domain, a positive pulse
lower bound, all fixed-order coefficient estimates, positive primary
covariance weights, a signed inverse with edge weights, and a gain for the
full LINEAR wave operator. It also applies the exact real-curl identity to
give a sharper primary covariance remainder.

The actual global leading profile, infinite coefficient sum, cutoff scales
and every auxiliary rectangle are not numerically materialized. Finite
diagnostics use explicitly manufactured smooth background data, including
the actual phase formula and huge integer carriers. They are not a numerical
construction of the global field or a proof of the continuum estimates.
Physical mean corrections, nonlinear interactions, a complete residual
improvement cycle, convergence, forcing, energy and blow-up remain open.

Parent: `9c490d8732da619500584458f9bf71f31079a8aa`.
Dependencies are the frozen [background-summation argument](../background_summation/README.md),
the [physical wave interface](../physical_wave_interface/README.md), and
their inherited leading-profile estimates. Those arguments have not
received an independent foundation review or a proof-assistant certificate.
No earlier evidence or GOAL.md is changed.

The supplied [proposed construction](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
is the primary source being reviewed: (4.26)–(4.27), (5.42), (6.1)–(6.32),
(7.2)–(7.40), and (9.1)–(9.2). Reviewed PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The source formulas are rederived below; its full conclusion is not assumed.

## 1. Inputs, quantifiers and the distinction from a reference pulse

Use viscosity one, as in the inherited normalized construction. Fix the
already selected leading profile, its finite radial frequency N,
and its construction exponent 0<h<1/100. Do not choose either again.
Put A=1/2+h, D=1/2-h. For a band ell set
Q=2^-ell, epsilon=Q^h, S=ell^2, k=ceil(epsilon^-1/2). A label fixes its
representative, sign, integer carrier, covering index and rectangle BEFORE
any differentiation.

The conditional inputs are:

1. On the closed active annulus and a fixed radial enlargement, the
   normalized background estimate from the preceding summation holds with
   every fixed composition of q partial_q, partial_X and partial_eta.
   Tangential errors are O(q^(2h)); the normalized radial velocity is O(q^h).
2. The leading profile has uniform positive shear/growth gaps on that
   compact set: F0>0, |g0|>0, N_theta<0, lambda0>0, and c0<0, bounded away
   from their relevant degeneracies. The normalized leading stress direction
   extends smoothly to the edges and has a strict cone margin. Its magnitude
   is at least c zeta, and every fixed derivative is bounded by
   C_I zeta delta^(-a_I). A bounded stress/zeta ratio is NOT assumed.
3. Prescribed slow boxes have diameter at most C S^-3. The auxiliary
   rectangle radius is fixed and positive; temporal speed ci obeys
   1/(Tg S)<ci<=1/S, so L=2r0/ci is comparable to S. Overlapping labels
   use the prescribed separated rectangles and bounded covering differences.

The source's separation construction is not numerically materialized here.
The path, common-lift and uniqueness arguments below are conditional on
that geometry. These hypotheses are not replaced by sampled success flags.

There is ONE qstar>0, fixed by value/gap comparisons and the primary cone,
on which every pulse inverse below is defined. Higher derivative orders
and finite harmonic sets change constants and polynomial degrees, not this
domain. Constants need not be bounded as the derivative order tends to
infinity. The background coefficients themselves can grow arbitrarily
fast with formal order; the inherited cutoff summation handles that issue.

The numerical choice h=1/200 later in this folder is a diagnostic parameter,
not a replacement for the actual construction's fixed h.

## 2. Transfer the actual background estimates into a fixed band chart

Write (R,Z,T)=(r/Q^(1/2),z/Q^D,(1-t)/Q), and s=q/Q. The implicit physical
coordinate equation becomes

\[
s-Z^2s^{2h}=T,\qquad X=R^2/(2s),\quad\eta=Zs^{-D}.          \tag{1}
\]

On a band 1/2<=s<=2 and tau>=0, |eta|<=1. The derivative of the left
side with respect to s is 1-2h eta^2>=1-2h>.98. Implicit differentiation
therefore bounds every fixed slow derivative of s,X,eta on the compact
radial enlargement, uniformly in the band. One-sided derivatives at tau=0
are retained. No power of Q is lost in this coordinate change.

Let b=Q^A ur_B, V=Q^A utheta_B, G=Q^A uz_B and F=V/R. If Fprof=E0/Rprof,
Rprof=R/sqrt(s), and gprof=(Rprof partial_Rprof Fprof,partial_Rprof U0),
the LEADING chart quantities satisfy

\[
F^{(0)}=s^{-(A+1/2)}F_{prof},\quad
G^{(0)}=s^{-A}U_0,\quad
g^{(0)}=(R F_R^{(0)},G_R^{(0)})=s^{-(A+1/2)}g_{prof}.       \tag{2}
\]

Here s is independent of R. Thus the shear direction and c0 are preserved,
and F, shear, growth and cone gaps stay uniform because s lies in a compact
positive interval. The inherited background estimates now give, for every
fixed slow derivative order j,

\[
|D^j b|\le C_j\epsilon,\quad
|D^j(F-F^{(0)},G-G^{(0)})|\le C_j\epsilon^2.               \tag{3}
\]

These are estimates for the actual conditional background, not for the
fixture below. If an asymptotic estimate at order j was initially stated
only below a j-dependent neighborhood, enlarge C_j on the remaining
closed q interval. There are only finitely many bands and locally finitely
many background summands there. Smoothness gives a finite enlarged
constant, so a common preliminary q domain can be used for all fixed j.

## 3. A uniform phase comparison, including rounding over the whole pulse

At the frozen leading representative define

\[
N=g_0/|g_0|,\quad K=(-N_z,N_\theta),\quad
\lambda_0^2=-2F_0N_\theta(2F_0N_\theta+|g_0|)>0,
\quad c_0=\lambda_0/(2F_0N_\theta)<0.
\]

The inherited leading cone gives rcone=sup |c0 (T.K)/(T.N)|<1, with
the limiting unit direction used at each edge. Choose one fixed
u>rcone/sqrt(1-rcone^2), which supplies a uniform strict margin. For sign sigma,

\[
B_s^2=\frac{\lambda_0}{\epsilon k^2(1+u^2)^{3/2}},\quad
s_\sigma(v)=\sigma u(1/2+v/L),\quad x_0=\sigma B_su/2,
\]
\[
(\widetilde p/R_0,p_z)=B_s\left(K-\frac{\sigma u g_0}{L|g_0|^2}\right).
\]

Take kp to be a nearest NONZERO integer to k ptilde, with a fixed tie rule.
The correct universal error is |p-ptilde|<=1/k. The often-used 1/(2k)
bound fails near zero: the nearest nonzero integer to .1 is 1, at distance
.9. No lower bound on |K_theta| is assumed. Every angular carrier kp is
nonzero, including representatives where K_theta vanishes.

Keep these numbers frozen. On a lifted rectangle, the actual phase and
its normalized spatial gradient are

\[
\Phi=p\theta+p_zZ/\epsilon+x_0R-v(pF+p_zG),
\quad n=\left(x_0-v(pF_R+p_zG_R),p/R,p_z-\epsilon v(pF_Z+p_zG_Z)\right). \tag{4}
\]

The auxiliary geometry gives Dr v=Dz v=0 and tstar v=1. Let Hphi=pF+pzG.
Since (ptilde/R0,pz).g0=-sigma Bs u/L, (3) and the box diameter give

\[
|n-B_s(s_\sigma,K)|
 \le C_\phi\{S^{-1}+S(S^{-3}+\epsilon+\epsilon^2+k^{-1})\},
\]
\[
|n'|\le C_\phi(S^{-1}+S^{-3}+\epsilon+\epsilon^2+k^{-1}).    \tag{5}
\]

The factor S multiplying rounding is essential: that error acts throughout
v<=L=O(S). The epsilon S term in the axial normal is also retained. Require

\[
S^2(\epsilon+\epsilon^2+k^{-1})\le1.                     \tag{6}
\]

Then (5) is O(1/S), uniformly over the whole pulse. Moreover

\[
E_{ik}=(t_*+bD_r+F\partial_\theta+GD_z)\Phi
 =\epsilon v[(H_\Phi)_T-G(H_\Phi)_Z]+bn_r.                \tag{7}
\]

Every fixed coefficient derivative retains the factor epsilon and has a
polynomial S bound. In fact the displayed expression is affine in v, so
the actual formulas give |D^I Eik|<=C_I epsilon S in ordinary slow/common
coordinates. These are coefficient derivatives, with the exponential
factored out and all discrete choices fixed.

### A sufficient threshold controlling every later band

Because 1/k<=sqrt(epsilon), the left side of (6) is at most
3 ell^4 exp(-a ell), a=h log(2)/2. For C>=1 and p>0, the inequality
log x<=x-1 applied to x=a ell/(2p) gives

\[
\log(C\ell^p e^{-a\ell})
\le \log C+p\log(2p/a)-p-a\ell/2.
\]

Thus the integer just above
max(1,2[log C+p log(2p/a)-p]/a) works for ALL later bands, not merely a
finite sampled set. Use p=4 and absorb the factor 3 and any desired finite
smallness constant into C. Increase the lower band further for the finite
normal, frame and cone gap conditions below and the preexisting q domain.
The actual C values have not been numerically evaluated; this is not a
numerical selection of the actual qstar.

## 4. Verify the actual moving-frame hypothesis

For n_tan=(n_theta,n_z), set Ka=n_tan/|n_tan|, Na=((Ka)_z,-(Ka)_theta),
sa=n_r/|n_tan|. Equation (5) and a sufficiently large lower band give a
uniform lower bound on |n_tan| and
Ka-K, Na-N, sa-s_sigma=O(1/S). Put

\[
B=[e_r-s_aK_a,N_a]
 \begin{pmatrix}1&1\\c_0\sqrt{1+s_a^2}&-c_0\sqrt{1+s_a^2}\end{pmatrix},
\]
\[
K_{bg}=\begin{pmatrix}0&-2F&0\\2F+RF_R&0&0\\G_R&0&0\end{pmatrix},
\quad A_n=-K_{bg}+\frac{n(n^TK_{bg}-(n')^T)}{|n|^2}.       \tag{8}
\]

B and its explicit left inverse are uniformly bounded. At the reference
normal, the algebraic projection in the (er-sK,N) coordinates is
[[0,2F0 Ntheta/(1+s^2)],[-(2F0 Ntheta+|g0|),0]], diagonalized by the
two columns above, with eigenvalues +/-lambda(v), where
lambda=lambda0/sqrt(1+s_sigma^2). The actual Kbg differs by
O(S^-3+epsilon^2); the normal difference, moving-normal term and B' are
O(1/S). Rational dependence on uniformly nondegenerate variables proves

\[
B^\ell(A_nB-B')=\operatorname{diag}(\lambda,-\lambda)+E,
\quad\max_{ab}|E_{ab}|\le C_E/S,\quad |d-d_{ref}|\le C_d/S, \tag{9}
\]

with d=epsilon k^2|n|^2 and
dref=lambda0(1+s_sigma^2)/(1+u^2)^(3/2). This verifies the FRAME INPUT of
the earlier interface lemma from the actual background estimates.

For each slow derivative order j>=1, normal derivatives cost at most S
per differentiated factor. Product/inverse chain rules therefore give
|D^j B|+|D^j Bleft|+|D^j E|+|D^j d|<=C_j S^j. Pulse derivatives are
bounded uniformly before ordinary auxiliary-coordinate conversion; n is
affine in v and |n'|=O(1/S). Constants may depend on j. The positive
spectral lower bound mu=min lambda is independent of band and label.

## 5. Upper bounds alone are insufficient: keep a positive radial pulse

Write the real homogeneous m=1 amplitude as t=B(z_plus,z_minus), with
z_plus(0)=P(0), z_minus(0)=0, P'/P=lambda-dref and P(L/2)=1.
As long as z_plus>0, rho=z_minus/z_plus satisfies

\[
\rho'=E_{21}+(-2\lambda+E_{22}-E_{11})\rho-E_{12}\rho^2.    \tag{10}
\]

Let eta=C_E/S and r=eta/mu. If eta<=mu/8, at either boundary of [-r,r]
the outward derivative is at most
-2 mu r+eta+2 eta r+eta r^2<0. Starting at zero, rho stays in that
interval. No numerical integration is needed for this continuum statement.
Meanwhile

\[
(\log(z_+/P))'=-(d-d_{ref})+E_{11}+E_{12}\rho,
\quad |\log(z_+/P)|\le (L/S)(C_d+C_E+C_Er).              \tag{11}
\]

This prevents z_plus from reaching zero and justifies continuation of the
ratio argument. The radial amplitude x=z_plus(1+rho) consequently obeys
c P<=x<=C P throughout the closed pulse interval. Its tangential ratio is

\[
t_{tan}/x=-s_aK_a+c_0\sqrt{1+s_a^2}\frac{1-\rho}{1+\rho}N_a
 =-s_\sigma K+c_0\sqrt{1+s_\sigma^2}N+O(S^{-1}).           \tag{12}
\]

This positive LOWER bound is what makes the primary covariance invertible.
An upper bound alone would also permit a pulse that disappears.

## 6. All fixed derivatives, sources and the common-torus path

In frame coordinates the homogeneous propagator obeys
norm V1(v,w)<=C P(v)/P(w), by the norm inequality from (9). Higher
harmonics have the exact scalar factor

\[
V_m(v,w)=e^{-(m^2-1)\int_w^v d(a)da}V_1(v,w),\quad m\ne0. \tag{13}
\]

Since d>=0, increasing a finite harmonic set does not worsen the value
bound or require a smaller qstar. To differentiate at FIXED v, divide
the frame solution by P. Its propagator is bounded and a coefficient
derivative of order i costs S^i. A differentiated Duhamel integral costs
one length factor S. Starting with the constant initial vector (1,0),
the recurrence
b0=0, bj=max_(1<=i<=j)(i+1+b_(j-i)) gives bj=2j. Pulse derivatives
follow directly from the ODE and add no further S power at fixed slow
derivative order. Multiplication by B gives the same bound for t.

For a source whose fixed-v derivatives have S powers a_j, multiplying by
Bleft and the normal projection gives forcing powers
g_j=max_(0<=i<=j)(i+a_(j-i)). The zero-data solution has
b0=g0+1 and bj=max(gj+1,max_(1<=i<=j)(i+1+b_(j-i))). Every fixed order
has a finite polynomial loss. The edge weight and epsilon exponent of the
source are retained because the slow point and X are fixed along the path.
This proves the weighted inverse mapping given the source bounds; it does
not assert bounds for every source of an unconstructed correction cycle.

On a common torus H with bounded covering difference Delta, the path is

\[
H_w=H+T_g^{-\Delta}c_i(w-v(H))v_t,\quad
D_HH_w=I-v_t\lambda_t,\quad D_Hv=T_g^\Delta\lambda_t/c_i.   \tag{14}
\]

Its higher derivatives vanish. Slow/transverse parameters stay fixed,
ordinary H derivatives cost at most S per pulse derivative, and no new
epsilon loss appears. Hence the homogeneous coefficient estimate is

\[
|D^I t|\le C_I S^{2|I|}P,
\quad |D^I\pi|\le C_I\epsilon^{1/2}S^{2|I|}P.              \tag{15}
\]

The pressure follows from
pi=i(n.Kbg t-n'.t+n.f)/(km|n|^2), retaining the moving normal and factor i.
Deck translations reindex copies and commute with (14). Uniqueness gives
descent to the common torus. It DOES NOT identify distinct preimages on
a finer band torus: an allowed source can differ on those preimages.
Vanishing source on an entire slow/transverse parameter rectangle implies
vanishing solution there; differentiated uniqueness gives its smooth zero
extension. No temporal zero extension is inferred before a cutoff.

## 7. Actual primary columns, quantitative positivity and edge weights

Let chi_g be the fixed nonzero transverse cutoff and psi the fixed relative
pulse cutoff, equal to one for |v-L/2|<=L/5. The auxiliary area measure is
|det(vr,vt)| ci dxi dv. Angular averaging contributes exactly 1/2, because
kp is a nonzero integer. Thus for a positive fixed constant cchi,

\[
h_\sigma=c_\chi c_i\int_0^L\psi^2x^2dv,\quad
H_\sigma=c_\chi c_i\int_0^L\psi^2x t_{tan}dv.              \tag{16}
\]

The reference damping gives a Gaussian bound
exp(-C(v-L/2)^2/(2L))<=P<=exp(-c(v-L/2)^2/(2L)). From (11),
the lower bound on x, and psi=1 on a central interval,
h_sigma is between c S^-1/2 and C S^-1/2. The weighted average of
|v-L/2|/L is O(S^-1/2): its unnormalized numerator is O(1), whereas
the denominator is at least c sqrt(L). Applying (12) yields

\[
H_\sigma=h_\sigma(-A_cN-\sigma uK+e_\sigma),\quad
A_c=-c_0\sqrt{1+u^2},\quad |e_\sigma|\le C S^{-1/2}.       \tag{17}
\]

The numerator statement uses integral (|v-L/2|/L)P^2 dv=O(1).
All constants are uniform over the actual labels under the stated inputs.

Here is an explicit robustness inequality for the positive inverse. In the
(N,K) frame let D0=[[-Ac,-Ac],[-u,u]] and each column error have Euclidean
norm at most e0. Factor H=(D0+Ecol) diag(h_plus,h_minus). Then

\[
\|D_0^{-1}E_{col}\|_\infty\le\beta
 :=e_0(A_c^{-1}+u^{-1}).
\]

For a target with A=-T_N/Ac>0 and |T_K|/u<=(1-eta_c)A, reference
weights lie between eta_c A/2 and A. If
beta<=eta_c/(4+eta_c), the Neumann-series error is at most eta_c A/4.
Both ACTUAL normalized weights therefore stay positive. Smooth extension
of the leading stress direction makes the frozen-representative cone
comparison uniform even at the flat edges. Increase the SAME lower band
to meet this finite value condition. Consequently

\[
\|H^{-1}\|\le C S^{1/2},\quad |\det H|\ge c/S,
\quad y=H^{-1}T_{0,*},\quad y_\sigma\ge c S^{1/2}\zeta.    \tag{18}
\]

There is no assumption that higher corrected stresses stay in this cone.
They are handled by a SIGNED inverse at this fixed positive primary wave.

For slow derivative order j, (15) integrated against P^2 gives
|D^j H|<=C_j S^(2j-1/2). Differentiating H H^-1=I gives
|D^j H^-1|<=C_j S^(2j+1/2). The stress estimates then imply

\[
|D^j y_\sigma|\le C_jS^{2j+1/2}\zeta\delta^{-a_j},\quad
|D^j\sqrt{y_\sigma}|\le C_jS^{2j+1/4}\sqrt\zeta\delta^{-b_j}. \tag{19}
\]

Each root derivative has y^(1/2-n) times n derivatives of y. Their
zeta factors cancel to sqrt(zeta); the lower bound in (18) justifies every
division. The actual edge weight is exponential in minus inverse squared
log-distance, so sqrt(zeta) times any fixed inverse power of delta is flat.
The roots extend smoothly by zero. All I=0 upper bounds may still have
inverse delta powers; bounded stress/zeta is not needed.

After multiplication by sqrt(epsilon), the pulse, transverse/time cutoffs
and the slow cutoff eta_beta, whose j derivatives cost S^(3j), the primary
harmonic amplitudes obey

\[
|D^I a_{primary}|\le C_I\epsilon^{1/2}S^{3|I|+1/4}
 \sqrt\zeta\delta^{-b_I}P.                               \tag{20}
\]

This establishes the conditional ACTUAL primary wave coefficient bounds.
For a real slow target increment Sigma in M_alpha, independent of theta
AND the auxiliary variable, set dSigma=H^-1(Sigma/epsilon) and
delta a_sigma=(dSigma)_sigma/(2 sqrt(y_sigma)). The same weighted
calculus proves the W_(alpha-1/2) bound. No sign or relative-size
restriction on Sigma is needed, but its weighted input estimates are
required. The exact covariance change retains the quadratic remainder
and interactions with all old/new curl corrections, as in the parent
checkpoint; the inverse is linear, not an exact nonlinear covariance solve.

## 8. Cutoff tails, exact curls and a primary covariance gain

P is positive at the finite pulse endpoints, however tiny. It cannot be
abruptly extended by zero there. Apply psi to potentials and pressures
BEFORE taking the curl. The principal cutoff error is psi' t for the
homogeneous pulse, and (1-psi)f+psi't for a sourced inverse. Both are
supported away from the central fifth, where P<=exp(-cL/50).

For the GLOBAL sourced linear residual, the original source harmonic must
already have smooth zero extension at its pulse endpoints. A general
inverse input satisfying only the coefficient envelope need not do so.
This additional source-support requirement is retained for any later cycle.

For every fixed derivative order, conversion to physical variables costs a
fixed exponential exp(B ell) and a polynomial ell^p. Gaussian endpoint
decay is exp(-a ell^2). To bound C ell^p exp(B ell-a ell^2) by q^M,
use log ell<=ell and take ell larger than
max(1,2[p+B+M log2]/a,sqrt(2 log C/a)). Thus ALL fixed physical
derivatives of these smoothly supported tails are flat. The tail threshold
may depend on derivative/decay order; the pulse inverse's domain need not.

In the normalized spatial operators, a radial coefficient derivative costs
at most epsilon^-kappa times a polynomial in S, where kappa=10^-5.
The exact real potential is -[n cross a/(k|n|^2)] sin(k Phi). From (20),
its curl correction r satisfies the sufficient bound

\[
|D^I r|\le C_I\epsilon^{1-\kappa}S^{3|I|+13/4}
 \sqrt\zeta\delta^{-b_I}P.                               \tag{21}
\]

The parent identity gives pointwise in the auxiliary variable
<w tensor w>_theta=(a tensor a+r tensor r)/2 for each real carrier:
the mixed cosine/sine terms vanish EXACTLY. Distinct labels have separated
supports. The auxiliary integral of P^2 contributes S^-1/2. Hence the
double-averaged primary cross-flux covariance in one slow box satisfies

\[
\mathcal C(w_{\beta})=\eta_\beta^2\epsilon T_{0,*}+\mathcal E_{C,\beta},
\quad |D^I\mathcal E_{C,\beta}|\le C_I\epsilon^{2-2\kappa}
 S^{3|I|+6}\zeta\delta^{-b_I}.                            \tag{22}
\]

The complete tensor correction is positive semidefinite before taking
particular cross entries. This sharper bound is for REAL PRIMARY carriers;
it is not asserted for arbitrary later complex corrections. Physical stress
assembly uses the squared slow partition: summing
Q^-2A eta_beta^2 epsilon (Q/q)^(A+1/2) T0 gives exactly
q^(-A-1/2) T0 in the principal covariance. The same bounded overlap and
chart comparisons preserve the remainder class. Physical stress
normalization Q^-2A gives covariance remainder power -1-2h kappa,
which is STILL NEGATIVE. Only the relative stress gain h(1-2kappa)>0
has been established. It is not absolute residual decay.

## 9. Full linear remainder, with the omitted terms kept

The principal amplitude equation handles fast auxiliary time, Kbg, leading
epsilon k^2|n|^2 damping and normal pressure. The complete linearized
physical operator from the parent checkpoint still has all the terms in
the following ledger, after factoring out a harmonic and normalizing:

| Term acting on a W_alpha principal amplitude | Gain above alpha |
|---|---:|
| slow time, b Dr and G Dz transport | 1-kappa |
| ik Eik a | 1/2 |
| background derivatives and radial frame term | 1 |
| pressure amplitude Dr pi, Dz pi | 1/2-kappa |
| epsilon amplitude Laplacian and vector connection | 1-2kappa |
| epsilon ik mixed phase/amplitude and phase-divergence terms | 1/2-kappa |
| epsilon ik angular vector-Laplacian connection | 1/2 |

For example epsilon k Dr a has exponent alpha+1/2-kappa, whereas
epsilon Dr^2 a has exponent alpha+1-2kappa. The slow pressure coefficient
has exponent alpha+1/2 before differentiation. The curl correction starts
in W_(alpha+1/2-kappa). Its principal operator preserves that exponent;
the remaining terms are no worse after their derivative losses. A common
conservative lower gain is therefore 1/2-3kappa>0. For the fixed primary
amplitude this is a W_(1-3kappa) LINEAR remainder plus additive flat tails.

`linear_remainder_review` independently differentiates a real harmonic with
the full fast radial operator and nonconstant pressure, then checks its
principal/remainder decomposition. It retains the derivatives of the radial
multiplier, mixed phase/amplitude viscosity and angular frame terms.
Earlier physical-coordinate tests compare the same full operator with
Cartesian momentum. The gain just derived uses the continuum pulse bounds,
not the numerical identity test. The quadratic velocity products and the
mean correction equations are not covered by this linear estimate.

## 10. What the diagnostics do, and the next unresolved estimate

The fixture uses s from (1), F0=R^-3 s^D(1+.1 eta+.03 eta^2),
G0=-1.5 R^-2 s^D(1+.1 eta+.03 eta^2), an O(epsilon) radial velocity
and specified O(epsilon^2) tangential perturbations. It is smooth data
for testing (2)–(9); it is not the actual global leading flow.

At diagnostic h=1/200 and band ell=51600, k=2^129 and L is about
7.69e8. Arbitrary precision retains the sub-mesh offsets and integer
rounding. The pulse is integrated through rho=zminus/zplus and
log(zplus/P), with a stiff solver; no exponentially tiny physical initial
value is supplied to a solver's absolute tolerance. Covariance is sampled
in the resolving coordinate (v-L/2)/sqrt(L), not on a uniform v grid.

The positive/negative-sign diagnostic gives sqrt(L) h about .403 and .318,
and covariance direction errors about 5.1e-9 and 5.8e-10. A tighter solve
checks stability. This quadrature is truncated to a central Gaussian
window and is recorded as a diagnostic; it is not an interval enclosure of
the full actual covariance. The continuum column estimates are (16)–(18).
The ratio diagnostic uses supplied entry/damping comparison constants
10000; sampled frame values do not certify those constants on a continuum.
The analytic proof instead establishes finite uniform constants from the
stated background inputs, without numerically selecting their values.

Other checks cover third-order profile/chart differentiation, phase
transport, analytic moving-frame derivatives, the all-later-band threshold,
192 perturbed-cone cases, square-root edge derivatives, the common-torus
path and exact rational epsilon budgets. Controls reject zero angular
frequency, a false universal half-integer rounding bound, missing fast time,
the axial phase derivative, moving basis, cover factor, mixed viscosity,
angular viscosity and slow pressure gradients. Inputs outside the certified
barrier/cone conditions are rejected, rather than being called proved.

```sh
python3 research/navier_stokes_cascade/results/actual_band_pulse_bounds/test_audit.py -v
python3 research/navier_stokes_cascade/results/actual_band_pulse_bounds/audit.py \
  --verify-record research/navier_stokes_cascade/results/actual_band_pulse_bounds/evidence.json \
  --output /tmp/actual_band_pulse_bounds_recheck.json
```

The new continuum step is a CONDITIONAL uniform actual-band pulse and
primary covariance estimate, with a linear wave gain. Its hypotheses are
the specific inherited background/leading bounds and support geometry;
they have not become unconditional by passing tests. Actual numerical
coefficient values or a usable numerical qstar are not supplied.

Next: use the COMPLETE angular-mean tensor of these exact curl velocities
in the physical mean system. Construct its pressure and velocity corrections,
enforce the integral constraints, and recompute all nonlinear and old/new
interactions for that same field. A successful next checkpoint must prove
a positive gain after a COMPLETE residual cycle. The present covariance
and linear gains do not establish that result or the blow-up conclusion.
