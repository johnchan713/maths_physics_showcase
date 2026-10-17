# Current foundation join: punctured terminal jets

Parent: `b87f53726641b240050a71789c025401e359abaa`.
Status: `current-source-punctured-terminal-jets-relative-to-retained-core-and-finite-analytic-lemmas`.
On 2026-10-17 parent hosted CI was running without a failed step; its
preceding finite-induction CI had passed. Historical evidence is unchanged.

The SAME field from [the current realization](../realization_foundation/README.md)
has compatible smooth one-sided potential, direct-swirl, velocity and
pressure jets at t=1 on every punctured terminal compact in a fixed late
observation cylinder. This joins the actual reviewed background, all-scale
physical pulse, compact mean and finite-source estimates to the terminal
interface, relative to the retained exact-core and finite analytic lemmas.
It does not retain a new blanket endpoint-bound hypothesis in place of
those estimates. Those core and finite existence/inverse/gap dependencies
still remain; no unconditional residual or blow-up theorem follows.

The [older conditional endpoint argument](../current_endpoint_jets/README.md)
is not promoted wholesale. We explicitly include the CURRENT enlarged
1/4<=q/Q<=4 chart range, derive the seed's regular Cartesian potential
factors, and give all-fixed-order MIXED physical radial/time majorants
for both heat swirl and the complete exterior pressure. This closes the
endpoint obligation of the current single-field checkpoint. Localization,
initialization, smooth global forcing, energy and comparison for this
same field remain next.

Source: the retained manuscript, Sections 6, 9 and 10 and its heat integral;
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its final theorem is not an input. Md64, h,N, positive lambda, width
C_setup^-20, roots, labels, geometry, pressure gauge and both existing
diagonal schedules are unchanged. Restricting an observation cylinder
does not retune the field or discard any previously accepted label.

## 1. Complete seed representatives, including the Cartesian axis

Keep A=1/2+h, D=1/2-h, tau=1-t, q=tau+z^2*q^(2h),
eta=z/q^D, X=r^2/(2q), L=1-2h*eta^2>=1-2h>49/50.
For the leading axis fields let M=int_0^X U and m=M/X. The actual
azimuthal potential and direct swirl have the regular Cartesian factors

\[
\Psi e_\theta=\frac12q^{-A}m(X,\eta)(-y,x,0),\qquad
B=q^{-A-1/2}F(X,\eta)(-y,x,0),\qquad p=q^{-2A}\Pi(X,\eta). \tag{1}
\]

Indeed r*Psi=q^(1-A)M, and its COMPLETE curl gives

\[
u_z=q^{-A}U,\qquad
r u_r=\frac{2\eta XU-2D\eta M-(1-\eta^2)M_\eta}{L}. \tag{2}
\]

The exact-core and background foundations supply every fixed jet of F,
U,m,Pi on the closed core. The regular quantity is M/X, not 1/E or a
fractional radial auxiliary map at zero. The positive-order background
corrections and finite seed wave/mean pieces use their ACTUAL cut
potentials; all annular auxiliary-dependent pieces vanish near the axis.
For q>=c this is a uniform hole r^2<2Xa*c. Thus (1), the reviewed compact
background jets and the annular physical bounds supply full seed
potential/direct-swirl/pressure jets, not just a bound on curl velocity.

Write that actual seed as curl A_seed+B_seed. The same realized field is

\[
\mathcal A=A_{seed}+\sum_j\chi(a_jq)A_j,\quad
\mathcal B=B_{seed}+\sum_j\chi(a_jq)b_je_\theta,
\quad u=\operatorname{curl}\mathcal A+\mathcal B,
\quad p=p^{[0]}+\sum_j\chi(a_jq)(p^{[j]}-p^{[j-1]}). \tag{3}
\]

All curls and cutoff commutators stay. B is separately divergence-safe
axisymmetric swirl; no compact full vector potential for it is invented.
Full seed pressure, both Pi1/Pi2 tails and all axial/moving-moment memory
stay in the SAME gauge. Neither gauge reset nor time-zero extension is
allowed at t=1. The global backward pressure identities, not mere velocity
support, make the exceptional pressure corrections vanish after Xv; choose
the already fixed heat exterior Xext beyond this full support hull.

## 2. Closed enlarged chart coverage on the q-positive region

The physical foundation estimates hold on the enlarged closed range
1/4<=q/Q<=4. On c<=q<=qmax, all such dyadic charts Q=2^-ell lie in

\[
\max\{\ell_0,\lceil\log_2(1/(4q_{max}))\rceil\}
 \le\ell\le\lfloor\log_2(4/c)\rfloor. \tag{4}
\]

The older 1/2<=q/Q<=2 window may describe the narrower active supports;
it must not be copied as coverage of ALL current enlarged charts. For
c=2^-24,qmax=2^-20, (4) includes ell=18,...,26; the narrow formula only
gives 19,...,25. These are exact algebraic boundary controls, not actual
accepted pulse thresholds or a numerical failure of the historical field.

Use the same observation margin

\[
q_e=\min(q_\star/16,2^{-\ell_0-2}),\quad
\tau\le q_e/16,\quad |z|\le(q_e/16)^D. \tag{5}
\]

The implicit equation implies q<=2tau+4|z|^(1/D)<=3q_e/8. Therefore,
with the FULL enlarged Q<=4q,

\[
4Q\le16q\le6q_e\le3q_\star/8<q_\star,
\qquad Q\le\tfrac32q_e\le\tfrac38\,2^{-\ell_0}. \tag{6}
\]

This corrects an attempted narrow-window transfer without changing any
physical domain parameter. q<qstar/2 alone would not establish (6).

For each finite band, compact slow R,Z,T boxes with T=tau/Q>=0 have
finitely many mesh centers: ceil(lo/delta)-1 through floor(hi/delta)+1.
At T=0 retain the already selected inclusive indices -1,0,1. Normalized
physical time differentiation has slow coefficient -epsilon=-Q^h.
No carrier, color or label is reselected. Closed support separation,
absolute-Y pullbacks, ODE uniqueness and the SAME physical q-based mean
bumps identify chart formulas and their zero-extension jets. The fractional
Y map is used only away from the axis; the uniform axis hole and (1)
handle the axis. T=0 has one-sided jets, not a temporal zero extension.

Since a_j>2^j/qstar, an active stage on q>=c must satisfy
j<log2(qstar/c). The independently chosen background schedule is locally
finite there too. All finite-stage constants are finite; no uniformity
as the band, stage or derivative order tends to infinity is needed.
The actual primary bounds hold uniformly to |eta|=1 by their positive
implicit denominator and one-sided Gamma estimates. The compact mean
maps have the same one-sided slow-parameter bounds; the finite generated
source grammar and complete curls preserve them. It follows for every
fixed multi-index alpha and time order j that, on an enlarged q-positive
punctured compact,

\[
\sup_{E^+\times[1-\delta,1)}
|D_x^\alpha\partial_t^j(\mathcal A,\mathcal B,p)|<\infty. \tag{7}
\]

Equation (7), including the NEXT time derivative, is the precise required
inequality. Interior smoothness, bounded values alone, and flatness of
the residual at the origin do not imply it. Here it is supplied by the
reviewed all-fixed-order physical estimates and the finite closed atlas,
relative to the retained core/finite-map lemmas.

## 3. Punctured cover and all-order exterior mixed jets

Let a terminal compact satisfy |x|>=d>0 in the observation cylinder.
Put b=min(q_e/4,d^2/(16Xext),1/4),
zeta=min(d/4,(b/16)^D), and tau<=b/16. On |z|>=zeta/2,
q>=|z|^(1/D)>=(zeta/2)^(1/D)>0 and (7) applies. On |z|<zeta,
r>d/2, q<=3b/8 and
X>=d^2/(3b)>=16Xext/3>Xext. Thus a WHOLE neighborhood of the equatorial
plane is unchanged heat exterior. The overlapping regions cover the
punctured compact, including axis points with z!=0. The origin is excluded.

In the heat region the actual fixed-gauge representatives are A=0,
B=K e_theta and

\[
K=C r^{-1-2h}H(Z),\quad Z=4\tau/r^2\ge0,\quad
p=-\int_r^\infty K(\rho,\tau)^2\frac{d\rho}{\rho},
\quad H(Z)=\frac1{\Gamma(1+h)}\int_0^\infty e^{-v}v^h(1+Zv)^{-h}dv. \tag{8}
\]

Let F_j=(-1)^j H^(j), A_j=(h)_j(1+h)_j,
G_j=sum_(k=0)^j binomial(j,k)F_k F_(j-k), and
S_j=sum binomial(j,k)A_k A_(j-k). Positive Gamma domination gives
0<=F_j<=A_j, F_j'=-F_(j+1), 0<=G_j<=S_j and G_j'=-G_(j+1).
The last identity is the FULL Leibniz convolution; 2F0*Fj omits terms
already at j=2. Differentiation under the radial pressure integral is
dominated by C_j*rho^(-3-4h-2j), integrable at infinity for every fixed j.

For any exponent a define c_(0,0)=1 and

\[
c_{n+1,k}(a)=-(a+n+2k)c_{n,k}(a)+2c_{n,k-1}(a). \tag{9}
\]

Differentiating r^-a*sum c_nk Z^k F_(j+k) proves (9): BOTH the changed
radial homogeneity a+n and the derivative of Z must stay. For
a=1+2h+2j and any fixed Zmax>=Z,

\[
|\partial_r^n\partial_t^jK|
\le |C|4^j r^{-a-n}\sum_{k=0}^n|c_{n,k}(a)|Z_{max}^k A_{j+k}. \tag{10}
\]

Put b_j=2+4h+2j. The pressure majorants are

\[
|\partial_t^jp|\le C^2 4^j S_j r^{-b_j}/b_j,
\]
\[
|\partial_r^n\partial_t^jp|\le C^2 4^j r^{-b_j-n}
\sum_{k=0}^{n-1}|c_{n-1,k}(b_j+1)|Z_{max}^k S_{j+k},\qquad n\ge1. \tag{11}
\]

The n>=1 formula differentiates p_r=K^2/r, keeping every product. On
r>=rmin>0, 0<=tau<=tau0 choose Zmax=4tau0/rmin^2. Cartesian and frame
derivatives are finite sums of these radial derivatives with bounded
inverse-radius factors; axial derivatives vanish. Thus (10)--(11) give
(7) for every fixed Cartesian space-time jet on the heat region as well.
At tau=0 the exact mixed jets are

\[
\partial_r^n\partial_t^jK=C4^j A_j(-1)^n(a)_n r^{-a-n},\qquad
\partial_r^n\partial_t^jp=-C^2 4^j S_j(-1)^n(b_j)_n r^{-b_j-n}/b_j. \tag{12}
\]

These follow by dominated endpoint limits, not by analytic continuation.
For C!=0 the time Taylor coefficient ratio for K is
(4/r^2)(h+j)(1+h+j)/(j+1)>(4/r^2)j. The full pressure convolution
also has factorial lower growth. Their time Taylor radii are zero at
fixed h>0. This does not prevent one-sided smooth jets or later smooth
jet gluing; it prevents assuming a convergent Taylor future or using the
positive-Z integral at negative Z.

## 4. Uniform limits, finite Taylor remainders and what is still open

For each representative G and each fixed alpha,j, (7) at order j+1 gives

\[
\|D_x^\alpha\partial_t^jG(t_2)-D_x^\alpha\partial_t^jG(t_1)\|_E
\le M_{\alpha,j+1}|t_2-t_1|. \tag{13}
\]

Every jet has a uniform one-sided terminal limit. Uniform limits in
spatial coordinate integrals on enlarged punctured boxes identify all
spatial derivatives; time integration identifies successive time jets.
Write T_j^G for the terminal trace of partial_t^j G, distinct from the
heat functions F_j above. For every FIXED finite k, time integration gives

\[
\left\|D_x^\alpha\left[G(1-\tau)-\sum_{j=0}^k
\frac{(-\tau)^j}{j!}T_j^G\right]\right\|_E
\le M_{\alpha,k+1}\tau^{k+1}/(k+1)!. \tag{14}
\]

Curl commutes with these finite limits. A velocity m-jet needs potential
bounds through m+2 and direct swirl through m+1 for its Cauchy estimate;
pressure needs m+1. All are fixed finite requirements in (7). No uncut
infinite series, uniform-in-order constant, new diagonal sequence or
future velocity is invoked. The pressure gauge and both base tails stay
fixed throughout the transfer. Actual flat velocities remain in the field;
the exact Gaussian residual ledger is not newly summed or inverted.

The result is compatible one-sided C-infinity jets on the PUNCTURED
terminal cylinder for the SAME current field. It supplies neither regular
velocity/pressure at the singular origin nor a global admissible force.
Next join full axisymmetric spatial/time localization and zero
initialization, keeping all curls, pressure, viscosity and interactions;
then smooth force extension, whole-space energy/dissipation and the exact
global solution comparison. Retained exact-core and finite analytic lemmas
must be discharged before unconditional claims.

## 5. Reproduction and diagnostic scope

Exact rational controls cover the full enlarged band boundary, closed mesh
endpoints, 36 mixed terminal pressure identities, the seed's Cartesian
curl and 21 finite Taylor plus 21 Cauchy bounds. Independent positive Gamma
integrals at two diagnostic h values compare 24 nonterminal mixed heat
derivatives and 18 complete pressure products against (9)--(11).
An independent quadratic-kernel chain rule is used inside those integrals.
Omitted radial chain terms, pressure cross terms, radial homogeneity
updates and the wrong terminal time sign are detected. Polynomial fixtures
are manufactured q=1 controls, not the actual fluid. Finite sampled
integrals and passing tests are not a continuum proof; the all-order
result comes from the derivations above.

The parent frozen chain is reproduced with 979 inherited provenance checks
and 18 direct source bindings. Use fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/endpoint_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/endpoint_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/endpoint_foundation/evidence.json \
  --output /tmp/endpoint_foundation_fresh.json
```
