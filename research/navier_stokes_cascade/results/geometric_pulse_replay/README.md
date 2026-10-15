# Pulse and source estimates with the constructive global geometry

The global radius from the preceding support audit is compatible with the
conditional pulse argument. This batch replaces its hidden radius constants
by a pulse-length calculation and keeps the transverse area in the actual
Haar covariance. It supplies, conditional on the inherited background and
leading-profile estimates:

* One all-later-band threshold for phase/frame comparisons, a positive radial
  pulse, primary cone inversion and every finite-harmonic sourced inverse.
* Explicit radius costs: covariance weights scale as r0^(3/2) S^(-1/2),
  their inverse as r0^(-3/2) S^(1/2), and primary amplitude roots as
  r0^(-3/4) S^(1/4). Shrinking the supports does not come for free.
* Fixed-order derivative, source-support and cutoff-tail estimates on the
  common tori, with no additional epsilon-power loss from the geometry.

The finite diagnostics now use the actual globally admissible radius,
not the earlier r0=.04 fixture. Their background remains manufactured:
neither these computations nor the preceding geometry construct the actual
global leading profile. No unconditional complete residual improvement,
Navier--Stokes solution or blow-up theorem is asserted.

Parent: `1d7a71fd249cdb136f2f8658f5fa85afb740f5dc`.
Dependencies: [constructive geometry](../auxiliary_geometry_audit/README.md),
[conditional actual-band estimates](../actual_band_pulse_bounds/README.md),
and their stated profile/background inputs. The primary source under audit
is the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(6.10)--(6.22), (7.2)--(7.29), Proposition 7.2 and Lemmas 7.1, 7.4--7.5.
Its previously checked SHA-256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
All inherited evidence, source files and GOAL.md are frozen.

## 1. Quantifiers and the surviving analytic hypotheses

Fix the inherited h in (0,1/100), finite leading radial frequency N and
leading profile. Set kappa=10^-5. Use the explicit 1,250-color support
construction with a fixed

\[
0<r_0\le r_{\rm geo}=1/6057072.
\tag{1}
\]

The radius, centers, representatives and every integer carrier are frozen
before differentiation or correction. For a band ell put Q=2^-ell,
epsilon=Q^h, S=ell^2, k=ceil(epsilon^-1/2). The exact covering scales give

\[
\frac1{T_g S}<c_i\le\frac1S,\qquad
L=\frac{2r_0}{c_i},\qquad
2r_0 S\le L<2r_0T_g S<S.\tag{2}
\]

The last inequality follows from T_g<6 and 12*rgeo<1. Pure geometry is
valid from ell>=16, but those early bands need NOT support a long pulse.
We will increase the first analytic band, without changing h or N.

The remaining hypotheses must be retained explicitly:

1. On the closed active annulus and its fixed radial enlargement, the
   actual normalized background b,F,G has all fixed slow jets, uniformly
   in band and label, with b=O(epsilon) and
   (F,G)-(F^(0),G^(0))=O(epsilon^2). This includes the one-sided jets at
   tau=0 for q>0. The common domain is the one justified in the earlier
   conditional band transfer.
2. The frozen leading representatives have bounded coefficients, F0>0,
   nonzero shear, lambda0>=lambda_min>0, |c0| bounded above and away from
   zero, and a strict two-column stress cone margin eta_c>0. Their limiting
   directions extend smoothly to the radial edges. Fix the same tilt u>0
   chosen from that cone. Bounds may depend on this finite leading profile.
3. The inherited stress and source edge bounds hold. A source for the
   inverse has fixed closed slow, transverse and shell support with smooth
   zero extension, and all coefficient jets bounded by
   C_I epsilon^alpha S^(b_I) sqrt(zeta) delta^(-a_I) P(v).
   It need not be orthogonal to the phase normal, divisible by cutoffs,
   or descend to the finer band torus.

Small q cannot repair a failed leading cone or a missing edge bound.
This replay establishes an implication from these hypotheses, not a new
independent proof of them. The constants below are not numerically certified
for the actual inherited leading profile.

## 2. The phase/frame errors are naturally O(1/L)

Use the original phase, including nearest NONZERO integer rounding,

\[
\Phi=p\theta+p_zZ/\epsilon+x_0R-v(pF+p_zG),
\quad n=(x_0-v(pF_R+p_zG_R),p/R,p_z-\epsilon v(pF_Z+p_zG_Z)).\tag{3}
\]

At the representative the target normal is B_s(s_sigma,K), with
s_sigma=sigma*u*(1/2+v/L), and
(ptilde/R0,pz)=B_s(K-sigma*u*g0/(L*|g0|^2)). The bound for the rounded
angular parameter is |p-ptilde|<=1/k, including near ptilde=0.
Let d_box<=C*S^-3 and assume L>=1. The formulas and hypothesis 1 give

\[
|n-B_s(s_\sigma,K)|
\le C\{L^{-1}+(1+L)(S^{-3}+\epsilon+\epsilon^2+k^{-1})\}.\tag{4}
\]

The 1/L term is the deliberate tilt in the initial tangential normal;
it cannot be replaced by a radius-independent C/S. The L/k rounding and
epsilon L axial terms act over the WHOLE pulse. Choose the all-band condition

\[
S^2(\epsilon+\epsilon^2+k^{-1})\le1.\tag{5}
\]

Since 1<=L<=S, L^2*S^-3<=1 and
L^2*(epsilon+epsilon^2+k^-1)<=1. Thus (4) and its pulse derivative imply

\[
|n-B_s(s_\sigma,K)|+|n'|\le C/L.\tag{6}
\]

These value constants can be chosen independently of r0 in (1), provided
the length conditions hold. In particular eps*k^2 is between 1 and 4,
so B_s, |n_tan| and |n| stay uniformly away from zero on a long enough
pulse. All representative coefficients stay in the same compact gap set.

Use the actual moving frame B, including B', and projected operator
A_n=-Kbg+n(n^T Kbg-(n')^T)/|n|^2 from the earlier checkpoint. Smooth
dependence on the nondegenerate variables, (6), and the background comparison
give constants C_E,C_d independent of r0 such that

\[
B^\ell(A_nB-B')=\operatorname{diag}(\lambda,-\lambda)+E,
\quad |E_{ab}|\le C_E/L,\quad |d-d_{ref}|\le C_d/L.\tag{7}
\]

Here lambda=lambda0/sqrt(1+s_sigma^2), d=eps*k^2*|n|^2 and
dref=lambda0*(1+s_sigma^2)/(1+u^2)^(3/2). B and its left inverse are
uniformly bounded. Their jth fixed-v parameter jets, and those of the
ODE coefficients, are bounded by C_j S^j: n is affine in v<=L<=S and
the product/inverse chain rules have at most j differentiated factors.
Constants depend on j and a finite harmonic set, not on correction stage
or a new q threshold. No derivatives of rounded choices are taken.

## 3. Positive homogeneous pulse and the common inverse domain

Let mu=lambda_min/sqrt(1+(3u/2)^2). For L>=8*C_E/mu the Riccati ratio
rho=zminus/zplus stays in [-a,a], a=C_E/(mu L)<=1/8. At either endpoint
of this interval its outward derivative is bounded above by

\[
-2\mu a+(C_E/L)(1+2a+a^2)<0.\tag{8}
\]

Starting from zplus(0)=P(0), zminus(0)=0, integrate the logarithmic equation:

\[
|\log(z_+/P)|\le C_d+C_E(1+a),\qquad
c_xP\le x=z_++z_-\le C_xP.\tag{9}
\]

For example take cx=(7/8)*exp(-Cd-9CE/8),
Cx=(9/8)*exp(Cd+9CE/8). Both are independent of radius and band. The
reference envelope is positive even at the finite endpoints; zero extension
is postponed until after the smooth time cutoff.

For a real nonzero harmonic m, the frame propagator has the exact factor

\[
V_m(v,w)=e^{-(m^2-1)\int_w^v d}\,V_1(v,w),\qquad
\|V_m(v,w)\|\le e^{2C_E+C_d}P(v)/P(w).\tag{10}
\]

The exponential with m is at most one. A finite harmonic set changes
derivative constants, while the domain is unchanged. A normalized
source bound |f|<=F0*P gives the value estimate

\[
|t_m(v)|\le C L F_0 P(v),\qquad
|\pi_m(v)|\le C_M\epsilon^{1/2}(|t_m(v)|+|f_m(v)|).\tag{11}
\]

The pressure is
pi=i*(n.Kbg*t-n'.t+n.f)/(km*|n|^2). Its normal source component and
moving-normal term are essential. The projected constraint evolves as
(n.t)'=-m^2*d*(n.t), hence zero initial data retain n.t=0.

At fixed v, differentiated Duhamel contributes at most one S per integral.
If coefficient jets cost S^i and forcing jets cost S^g_j, the degrees satisfy
b0=g0+1 and
bj=max(gj+1,max_(1<=i<=j)(i+1+b_(j-i))). All are finite for fixed j;
the epsilon order alpha and the edge factors are unchanged. For homogeneous
data, division by the frozen P gives b_j=2j. Pulse derivatives follow
from the ODE; no derivative of its initial data is lost.

On a common lift, D_H H_w=I-vt*lambda_t, higher derivatives vanish,
and D_H v=Tg^Delta*lambda_t/ci=O(S). Thus returning to H introduces
only polynomial S losses, with Delta<=3. A mixed-band source is evaluated
on EACH lift. Averaging distinct preimages can erase a valid source and
is not a solution operator. Zero-data uniqueness and the fixed path prove
the slow/transverse/shell support and overlap identities, including all
boundary jets. Temporal support is not inferred before the time cutoff.

## 4. Actual Haar covariance: keep the small transverse area

Make the transverse cutoff explicit:
chi_g(xi)=beta(2xi/r0), where beta(s)=exp(1-1/(1-s^2)) for |s|<1 and
zero otherwise. Its closed support [-r0/2,r0/2] lies strictly inside
the rectangle's transverse interval, as required by the source convention.
Then

\[
\tfrac12r_0e^{-2/3}\le\int\chi_g^2d\xi\le r_0.\tag{12}
\]

Let psi(v/L) equal one for |v/L-1/2|<=1/5 and vanish for
|v/L-1/2|>=2/5. Angular averaging supplies 1/2. Normalized Haar
averaging on every common torus supplies exactly the same integral as
on the band torus, not a factor 14^Delta:

\[
h_\sigma=\frac{1+b_g^2}{2}\left(\int\chi_g^2d\xi\right)c_i
 \int_0^L\psi^2 x^2dv.\tag{13}
\]

The reference net growth has derivative between -CG/L and -cG/L.
Explicit valid constants are

\[
c_G=u\lambda_{min}(u/2)
 [(1+(3u/2)^2)^{-3/2}+2(1+u^2)^{-3/2}],
\]
\[
C_G=u\lambda_{max}(3u/2)
 [(1+(u/2)^2)^{-3/2}+2(1+u^2)^{-3/2}].\tag{14}
\]

Twice integrating around the midpoint yields
exp(-CG*(v-L/2)^2/(2L))<=P<=exp(-cG*(v-L/2)^2/(2L)). If L>=25,
the central interval |v-L/2|<=sqrt(L) lies where psi=1. Consequently

\[
2c_x^2e^{-C_G}\sqrt L
 \le\int\psi^2x^2dv\le C_x^2\sqrt{\pi L/c_G}.\tag{15}
\]

Combine (12)--(15) with ci=2r0/L and (2):

\[
c\,r_0^2 L^{-1/2}\le h_\sigma\le C\,r_0^2L^{-1/2},\qquad
h_\sigma\asymp r_0^{3/2}S^{-1/2}.\tag{16}
\]

The constants in these comparisons can be independent of r0. Omitting
the transverse integral would give the wrong radius exponent. There is
no radius-uniform lower bound after dividing only by S^-1/2.

The tangential ratio differs from -s_sigma*K+c0*sqrt(1+s_sigma^2)*N
by C/L. Also integral |v-L/2|*P^2/L dv<=1/cG. Dividing by the positive
lower mass in (15) proves

\[
H_\sigma=h_\sigma(-A_cN-\sigma uK+e_\sigma),\quad
A_c=-c_0\sqrt{1+u^2}>0,\quad |e_\sigma|\le C_{dir}/\sqrt L.\tag{17}
\]

Use the same quantitative cone inverse as before. It suffices that

\[
\frac{C_{dir}}{\sqrt L}(A_c^{-1}+u^{-1})
 \le\frac{\eta_c}{4+\eta_c}.\tag{18}
\]

The two actual primary weights stay positive. Equations (16)--(18) give
||H^-1||<=C*r0^-3/2*S^1/2, |det H|>=c*r0^3/S and
y_sigma>c*r0^-3/2*S^1/2*zeta for the inherited target. Upper and
derivative bounds retain their inverse-distance powers; stress/zeta
need not be bounded at a flat edge. Smooth root calculus gives

\[
|D^j\sqrt{y_\sigma}|\le
C_j r_0^{-3/4}S^{2j+1/4}\sqrt\zeta\,\delta^{-b_j}.\tag{19}
\]

The uniform covariance derivative bounds follow by integrating pulse jets
against P^2; labels, ci and L do not depend on the differentiated slow
point. No derivative of Haar coordinates or a discrete carrier is inserted.

## 5. One finite all-band choice and fixed-order coefficient bounds

Let Lmin be the maximum of 25, the finite normal/frame gap length,
8CE/mu, and the square of
Cdir*(Ac^-1+u^-1)*(4+eta_c)/eta_c. Take the first band above

\[
\max\{\ell_{background},16,\sqrt{L_{min}/(2r_0)},
\ell_{exponential}(h)\}.\tag{20}
\]

The earlier all-half-line inequality
3*ell^4*exp(-h*log(2)*ell/2)<=1 constructs a finite
ell_exponential by log x<=x-1. This proves (5) for EVERY later band.
Together (2) and (20) give all value and cone conditions. Derivative
orders, source orders and finite harmonic sets change bounds, not (20).
The actual Lmin and background threshold are not numerically evaluated.
Use the uniform lower bound of Ac and uniform constants over all leading
representatives in (20); a different threshold for each label would not
establish this conclusion.

The normalized fixed cutoff derivatives cost r0^-j in common auxiliary
coordinates. Since L>=25 and L<2r0*Tg*S, r0^-1<(2Tg/25)S<S;
these derivatives cost at most C_j S^j. Slow mesh derivatives cost S^(3j),
which dominate the pulse S^(2j) and cutoff S^j losses under Leibniz' rule.
Thus the primary amplitudes satisfy

\[
|D^I a_{primary}|\le
C_I r_0^{-3/4}\epsilon^{1/2}S^{3|I|+1/4}
 \sqrt\zeta\,\delta^{-b_I}P.\tag{21}
\]

In particular the derivative constants can be taken uniform in the
admissible radius family with (20); the explicit baseline r0^-3/4
remains. Smooth flat edge extension is the inherited exponential-weight
argument, and every fixed derivative is supported in the same closed
slow, shell and auxiliary supports.

The radial graph multiplier is still
Mi<=epsilon^-kappa*S^-rho; covering changes have bounded norm. Thus the
earlier full linear remainder gain 1/2-3kappa survives at each fixed
radius. In the real primary covariance remainder, the square of (21)
contributes r0^-3/2 and Haar pulse integration contributes r0^3/2.
These baseline radius powers cancel. This is a coefficient estimate and
conditional LINEAR/covariance gain, not an independently verified nonlinear
complete-residual improvement or absolute physical residual decay.

## 6. Sources and cutoff tails on the actual geometry

The source inverse preserves slow and transverse support because each
path fixes those parameters and zero data give zero solution outside them.
The all-band support construction ensures that distinct actual wave labels,
including derivatives of their smooth extensions, have zero products.
For one label, two wave factors with weights sqrt(zeta)*P have product
zeta*P^2<=sqrt(zeta)*P; a mean times a wave retains the wave support and
envelope. These facts justify using the same source classes in the
inherited finite-cycle algebra. They do not independently re-prove the
profile/background/mean inverse bounds assumed there.

For a sourced pulse, the original harmonic must additionally have temporal
smooth zero extension before it is treated as a GLOBAL residual. The
envelope alone is insufficient. After multiplication by psi the principal
uncancelled source is (1-psi)f+psi't, supported where

\[
P\le e^{-c_GL/50}\le e^{-(c_Gr_0/25)\ell^2}.\tag{22}
\]

Cutoffs are applied to potentials and pressures before the exact curl;
all their commutators are retained. At any fixed derivative order and
fixed epsilon order the physical conversion costs at most
C(r0)*ell^p*exp(B*ell). Every fixed target q^M is dominated by (22)
on all sufficiently late bands; `flat_tail_threshold` gives a sufficient
integer using log ell<=ell. These asymptotic tail thresholds may depend
on derivative, stage and target order; they do not redefine the common
inverse domain (20). Uniform flatness in the limit r0->0 is not claimed.

## 7. Independent diagnostics and the next actual obstruction

The numerical fixture uses the frozen actual phase/formula code, the
globally admissible rgeo and rgeo/2, and several huge bands. It keeps
arbitrary precision for integer carriers and normal rounding. The ODE is
solved after dividing by the exact P, avoiding underflow at the seed.
Covariance quadrature uses the resolving coordinate (v-L/2)/sqrt(L)
over the FULL pulse, includes both cutoffs and transverse area, and is
repeated more tightly. A separate three-component projected solve checks
the two-frame sourced solve, its normal constraint and normal pressure.
The covering-copy counterexample follows paths from the bottom of the
actual rectangle, at transverse coordinate zero. Multiplying by fixed
band-transverse and slow bumps equal to one there supplies the required
compact source support without changing any checked path values.
Finite data remain diagnostics; the all-band bounds are (2)--(22).

Deliberate failures include using the old test radius as global geometry,
accepting short geometry-only bands, dropping transverse area, claiming a
radius-uniform inverse, removing a normal source from pressure, averaging
distinct common-torus preimages, and inferring endpoint zero extension
from an envelope bound alone.

The explicit geometry introduces no new epsilon-exponent obstruction:
it imposes a finite later starting band and the displayed radius costs.
The next unverified analytic premise is now the actual leading profile
and summed-background estimate in hypothesis 1, together with the strict
gaps and all-order inverse inputs in hypotheses 2--3. A direct dependency
audit must justify those bounds on their common closed collar at the SAME
fixed N and h. In particular, for every fixed multi-index I it must justify

\[
\sup_{\ell\ge\ell_{bg},\gamma}\sup_{\mathcal D_\gamma}
\left\{\epsilon^{-1}|D^Ib|+
\epsilon^{-2}|D^I(F-F^{(0)},G-G^{(0)})|\right\}<\infty,
\tag{23}
\]

with the stated nonzero infima of the leading value/cone gaps and the
weighted edge bounds. D_gamma is the closed band domain with its fixed
radial enlargement and one-sided tau=0 boundary. Numerical phase fixtures
cannot supply (23). This is the remaining input to an UNCONDITIONAL pulse
estimate, rather than a failure of (20) for positive fixed radius.
The downstream conditional exclusion chain is unchanged.

```sh
python3 research/navier_stokes_cascade/results/geometric_pulse_replay/test_audit.py -v
python3 research/navier_stokes_cascade/results/geometric_pulse_replay/audit.py \
  --verify-record research/navier_stokes_cascade/results/geometric_pulse_replay/evidence.json \
  --output /tmp/geometric_pulse_replay_recheck.json
```
