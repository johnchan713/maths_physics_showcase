# An exact analytic core attached to the reference exterior

This note replaces the finite Taylor pilot by a contraction for an **exact
infinite analytic series**. It bounds the continuation and incoming integrals,
then applies the existing five-moment annulus. The amplitude and transition widths
are specified by finite formulas; their actual numerical values are too large
or small to materialize. Computation uses logarithms.

Status: `analytic-core-and-reference-attachment-bounded`. This is an analytic
construction with outward numerical checks of its constants, not a formal
proof-assistant certificate or an independent peer review.

The result here is a matched **leading profile**, using the frozen exterior
construction and its exact pressure datum. The later Appendix C stress
modulation, full momentum corrections, smooth forcing, and velocity blow-up
are separate obligations. In particular, zero leading stress in the core does
not mean zero full Navier-Stokes residual.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pages 26-29 and 144-157, especially (4.9), (4.16), and (B.1)-(B.40).
PDF SHA-256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent: `7a258d209025806db0897843cb1abe92d2c4b2be`.
The previous audits and their evidence are frozen.

## 1. Fixed data and a complex neighborhood

Use the completed exterior's exact pressure datum. No sampled pressure fit is
substituted. Its representation, established in the preceding annulus note, is

\[
-\Pi_0(\eta)/P^2=\int_0^1(1+\eta^2)^{-2\theta}\,d\mu(\theta),
\qquad \mu\ge0,\quad \mu([0,1])<5,
\]

with an ideal contribution \((5/2)(1+\eta^2)^{-2}\). The argument is uniform
for \(P\ge16\), \(0<h\le.01\), and these pressure data. It therefore applies to
the selected Md64 exterior without evaluating its astronomically large P.
The Md4 diagnostic does not repair that family's historical outer-cone failure.

Set \(\epsilon=10^{-16}\), \(j=10^{-18}\), \(\sigma_*=10^{-20}\), and use
the axis data (B.1)-(B.3), including \(U_*=4\eta+j\). Put

\[
r_0=\sigma_*/256,\qquad \rho=r_0/8,\qquad M=200/\sigma_*^2,
\qquad Z_b=1024P^2/r_0.
\]

The source notation remains \(A=.5+h\), \(D=.5-h\),
\(d=1-\eta^2\), \(L=1-2h\eta^2\). In particular
\(H_*=D\eta+dU_*\), \(\chi=H_*^2/(H_*^2+\sigma_*^2)\),
\(\zeta_*=-LH_* /(H_*^2+\sigma_*^2)\), and
\(\phi_*=\exp(\Lambda\int_0^\eta\zeta_*(w)\,dw)\).

On the rectangle \(|\Re\eta|<1+r_0, |\Im\eta|<r_0\), \(|\eta|<2\),
\(|L|>.9\), and \(|H_*\pm i\sigma_*|\ge3\sigma_*/4\). Indeed the
vertical variation of H is at most \(53r_0<\sigma_*/4\).
Consequently \(|\zeta_*|\le M\) and
\(|\chi|=|1-\sigma_*^2/(H_*^2+\sigma_*^2)|<3\).
The positive-mixture pressure has a single analytic branch here:
\(\Re(1+\eta^2)>1-r_0^2\), so \(|\Pi_0|\le10P^2\).
On the half-size rectangle Cauchy's bound gives
\(|\Pi_{0\eta}|\le20P^2/r_0\) and \(|Z_*/L|\le Z_b\).

For a coefficient holomorphic on the half-size rectangle and bounded there
by b, its B.4 norm at radial degree zero is at most b. Cauchy's inequality
leaves the factor \((k+1)^2/4^k\le1\). This includes every angular derivative,
not just the derivatives computed in a finite Taylor jet.

## 2. An explicit contraction for the entire axis series

Use the space B_rho with weights

\[
a_{nk}=\frac{20^{-n}\rho^{-k}k!\binom{n+k}{k}}{(n+1)^2(k+1)^2}.
\]

Its algebra constant is at most \(\mathcal A=256\), since
\((8\sum_{i\ge1}i^{-2})^2<256\). Radial averaging has norm at most one.
The radial integrals J1, J2, I and multiplication by Y have norm at most 80.
Set \(\mathcal D=80/\rho\). A J-integrated monomial containing m unknown factors, one
fixed analytic coefficient of norm b, and at most one angular and one
logarithmic radial derivative on distinct unknown factors has bound
\(\mathcal D\mathcal A^m b\) times the product of the unknown norms. Omitting either derivative
only improves the estimate. An angular derivative on a radial average is
allowed. These statements follow from (B.9)-(B.10) and convolution; no
unintegrated angular differentiation is declared bounded on B_rho.

The order-one angular operator must be inverted before contraction. For
\(T=J_2\chi/2\), the degree-raising estimate gives

\[
\|T^k\|\le\frac{(40\cdot256\cdot3)^k}{k!(k+1)!},\qquad
\|(1+T)^{-1}\|\le V:=e^{30720}.
\]

The center is \(\Phi_0=f_0(Y\chi)\), \(u_0=-YZ_*/(2L)\). On its unit
ball both unknown norms are bounded by \(S=2+V+40Z_b\).
All estimates below hold for \(\Lambda\ge1\) and \(\|g\|_\rho\le1\).

Expanding the **original** source equations gives the following positive
majorants. The coefficient bounds already include \(L^{-1}\).

| Remainder | Linear coefficient sum | Quadratic coefficient sum |
|---|---:|---:|
| J R1 | 386 | 40 + 10 M |
| J R2, excluding p | 442 | 32 |
| J R2 terms containing p, p_eta, Dp | 34 times the norm of p | |

For example the two derivative products in R1 and R2 are
\((\partial_\eta A_Xu)D_X\Phi\) and \((\partial_\eta A_Xu)D_Xu\).
Their derivatives act on distinct factors, which is essential. The pressure
\(p=I(g^2\Phi^2)\) has norm at most \(80\mathcal A^3S^2\), and is retained.
Define

\[
M_1=\mathcal D[386\mathcal AS+(40+10M)\mathcal A^2S^2],
\]
\[
M_2=\mathcal D[442\mathcal AS+(32\mathcal A^2+2720\mathcal A^4)S^2].
\]

Let L1, L2 be the derivatives of these two positive polynomials in S, and set
\(F_b=V(M_1+M_2)\), \(L_b=V(L_1+L_2)\). Subtracting a product one factor
at a time proves these Lipschitz bounds in the maximum norm on the pair.
The actual fixed-point map has displacement at most F_b/Lambda and Lipschitz
constant at most L_b/Lambda. Factors of one half have only been discarded
in the direction of increasing the bounds.

For all mixed derivatives with radial order at most two and angular order
at most six on \(0\le Y\le4.1\), use \(E_b=10000\rho^{-6}\).
This follows by summing the coefficient weights on \(|Y|\le5\) and then
using Cauchy's inequality on a radial circle of radius .9. Choose

\[
K=\{10^6E_b(1+S+F_b+L_b+M+Z_b+P^2)\}^{16},
\qquad \boxed{\Lambda=10^{12}K^8\epsilon^{-2}}.
\]

The closed unit ball is invariant and the map contracts by less than 1/2.
Completeness therefore gives an exact analytic solution of (B.15) with
norm error at most F_b/Lambda. The usual contraction iteration specifies
the solution; a finite polynomial is not substituted for its limit.

The scalar alternating-series bound gives \(f_0(Y\chi)>.265\) on the
real rectangle. The derivative evaluation above, products, and reciprocals
with this positive denominator give, in particular,

\[
\Phi>.264,\qquad
\|U(4/\Lambda,\cdot)-U_*\|_{C^1}\le K/\Lambda<\epsilon/100.
\]

K also bounds the needed fixed-data derivatives, the comparison function,
its logarithm, u0, and the coefficients of their error bounds through the
stated finite orders. To see that the power 16 suffices, mixed derivatives
of reciprocals/logarithms through total order eight involve at most nine
factors of the evaluation bound and denominators at least .264; taking a
difference adds one error factor. The finite combinatorial coefficients
are absorbed by \((10^6)^{16}\), leaving room for these ten factors.

The normalized solution keeps the pressure term nonzero: C will be chosen
so that g=phi_star/C is analytic and bounded by exp(-Lambda) on the original
rectangle. Its derivatives are controlled by Cauchy's inequality, not by
setting g or its pressure contribution to zero.

## 3. Exit from the core and bounds for a reference family

The angular source needs an estimate that remains valid when chi=0.
The identity \(|H_*\chi'|\le2|H_*'|\chi\) controls the comparison profile.
Using \(|\zeta_*|\le2\sigma_*^{-1}\sqrt\chi\) on the real interval,
the remainder in the source is bounded by

\[
100K^2\chi+100K^2\sqrt\chi+100K^2/\Lambda.
\]

Absorb the first term in .02 Lambda L chi. Young's inequality absorbs the
second in .03 Lambda L chi with constant error below
\(10^6K^4/\Lambda<.01\). Since \(-W_*>2.8\) and the h term is at most .1,
this proves \(S_q\ge2.5+.95\Lambda L\chi\) throughout the core.
Thus \(p_1=-2Y\Phi_Y/\Phi>0\), and \(n_s=-2u_Y\).

At Y=4, the alternating polynomial in (B.19) gives p1>2.3 wherever
chi>.99. On its complement the previously proved separation gives
\(|Z_*|>\delta_*=jP^2/100\). The error bound gives
\(|n_s|\ge\delta_*/2\ge j\). No division by chi is used.

Put X0=4/Lambda, Xb=100, Xi=110. For a width t1, construct the reference
phi_r,U_r exactly by (B.22): agree with the analytic solution up to t1 in
\(y=\log(X/X_0)\), taper both native logarithmic slopes over [t1,2t1],
then freeze them. Only cutoffs with values in [0,1] are integrated.
In particular there is no inverse-width loss in angular norms.

For later use choose the deliberately large, amplitude-independent number

\[
\boxed{B=10^{20}(1+\Lambda+K)^{64}}.
\]

For widths with \(100K^2t_1<j/100\), and amplitudes large enough that
\(C^{-2}e^B<j/100\), the same source argument gives
\(S_{q,r}\ge2.4+.94\Lambda L\chi\) and \(l_r\le1\).
The axial source differs from Zstar by at most
\(100K^2/\Lambda+100K^2t_1+C^{-2}e^B\), preserving its sign on the
complementary region: the native term \(100K^2/\Lambda\) is also below
j/100. The forward equation for ns preserves its lower bound. This j-scale
test is stronger than the constant error test needed for Sq.

Here and below use the angular algebra norm
\(\|v\|_k=\sum_{i=0}^k\sup|\partial_\eta^iv|/i!\). At k=1 this is
the exact norm used by the earlier joining criterion.
The following family bounds through k=3 hold on [X0,Xi], uniformly in the
subsequent C and the smaller reference widths:

\[
\|p_{1,r}\|_3,\quad\|n_{s,r}\|_3,\quad
\|p_{1,r}^{-1}\|_3,\quad\|\log(CE_r)\|_3 < B/100,
\]
\[
\|CE_r\|_3,\quad\|(CE_r)^{-1}\|_3<e^B,\qquad p_{1,r}\ge X/10.
\]

The bounds can be obtained directly from the forward source integrals.
The factor phi_star cancels from p1's integral. The remaining normalized
reference Phi and its reciprocal lie between fixed positive constants and
have angular norms bounded by \(10^3K^3\). The source bound gives
\(\|p_{1,r}\|_3\le10^6X\Lambda K^8\) and
\(\|n_{s,r}\|_3\le10^6K^8\).
For the reciprocal, combine p1>=X/10 with the differentiated reciprocal
formula; its bound is at most \(10^{24}\Lambda^4K^{24}\).
The logarithm of CE is bounded by integrating Lambda*zeta and adding
log Phi and half log(2X), at most \(10^6\Lambda K^8\).
Each of these bounds is below B/100. Exponentiating the last bound and
differentiating at most three times proves the displayed exponential bounds.

Integrating \(D_Xp_{1,r}=XS_{q,r}/L-l_rp_{1,r}\) also gives

\[
p_{1,r}(y)\ge p_{1,r}(0)e^{-y}+1.88\chi(e^y-e^{-y}),
\qquad p_{1,r}(X)\ge1.2(X-X_0^2/X).
\]

Consequently p1,r>3 for X>=100. The first comparison preserves p1,r>2.3
on chi>.99. Elsewhere \(|n_{s,r}|\ge j\), and the choice of C below makes
the axial contribution large. Uniformly on [X0,Xi],

\[
2.3<v_r=p_{1,r}+p_{2,r}^2/p_{1,r}\le C^2e^{6B}<C^3/2,
\qquad |p_{2,r}|<Ce^{2B}.
\]

Explicitly, on the complementary region
\(p_{2,r}^2/p_{1,r}\ge16j^2C^2e^{-2B}/(\Lambda^2B)>2.3\).

## 4. Parameter order and an activation estimate without inverse kappa

Choose the transition length **before C**:

\[
T_{sh}=2000B,\qquad
\boxed{\log C=10000B+100\log(1/\epsilon)},
\qquad \delta=C^{-20}.
\]

Then \(\log C\ge\Lambda(2M+1)\), so the complex bound on g required
by the core contraction holds. The core, reference, and final continuation
may depend on C, but all the preceding bounds use the same K,Lambda,B.
Use \(t_1=\kappa_0=\delta\), and give each of the two final cutoffs
logarithmic width delta. These cutoffs fit between X=100 and X=110.

For 0<y<t1 use the actual fields defined by (B.26):

\[
e_a=(1-\kappa_0)\sigma(y/t_1),\quad\kappa=1-e_a,\quad
\partial_y\log\phi=-\kappa p_{1,r}/2,\quad
\partial_y U=-\kappa Xn_{s,r}/2.
\]

Afterward keep kappa=kappa0 until the final cutoffs. Forward integration
from the exact axis values defines actual fields, not manufactured moments.
Subtracting the reference slopes gives field differences bounded in C2 by
\(B^2ye_a\) during activation and by \(B^2(t_1+\kappa_0)\) up to Xb.
The logarithmic interval has length less than B; the axial integral is
bounded using X<=110. Derivatives of the cutoff itself are not used here.

Write s=ye_a or s=t1+kappa0 as appropriate. The original five-moment
formulas give, for \(r=E_r/E\),

\[
|p_1-p_{1,r}|\le e^{20B}s,\quad
|p_2-p_{2,r}|\le Ce^{30B}s,\quad |r-1|\le e^Bs.
\]

One way to track these generous constants is to scale I,J and E by C
before taking differences. Their cumulative differences cost at most
\(e^{4B}s\); M,S cost the same, and Cp costs \(C^{-2}e^{3B}s\).
On X>=X0, X inverse is at most B, and the scaled H inverse and its first
two angular derivatives cost at most exp(2B). The Q formula then costs
at most exp(9B)s. The N formula followed by division by E costs at most
C exp(10B)s. The displayed 20B and 30B bounds include these products
and both terms of the reciprocal-difference identity. This argument uses
one more angular derivative of the fields, but no radial derivative of
their differences. It retains the common pressure datum.

There is an exact cancellation that must occur **before** estimating:

\[
t_s=p_{2,r}r/p_{1,r},\quad
v_s=\kappa(p_{1,r}+p_{2,r}^2r^2/p_{1,r}).
\]

Thus no 1/kappa0 appears. With \(D_c=C^4e^{100B}\), exact algebra gives

\[
|P_c-v_r|,\quad |v_s-\kappa v_r|,\quad |J_c|\le D_cs,
\qquad D_c<C^6.
\]

For example
\(P_c-v_r=\Delta p_1+(p_{2,r}^2/p_{1,r})(r-1)
+(p_{2,r}/p_{1,r})r\Delta p_2\).
The polynomial oracle checks this identity and the corresponding Jc,vs
identities with exact rational coefficients.

On activation, \(P_c-v_s\ge2e_a\), \(P_c>2.2\), and

\[
\frac{(v_s-2)_+J_c^2}{2(P_c-v_s)^2}
\le C^3D_c^2t_1^2/8<1.
\]

This proves the admissible quadratic test when vs>2 and the strict relaxed
test when vs<=2. A first collar with vs>2 is explicit: take
\(t_c=t_1/(10\sqrt{\log C})\). On it
\(e_av_r\le e^4C^{-97}<.01\), using the flat step's endpoint bound.
The factorization at the inner edge also survives: near zero
\(e_a=e^{-t_1^2/y^2}g_a(y)\), with ga smooth and positive, and

\[
e_a(y)^{-1}\int_0^y e_a(u)b(u,\eta)\,du\in y^3C^\infty.
\]

This follows by the change of variable v=t1^2/u^2, or the integration
identity in Lemma A.9. The field differences and then their cumulative
moments therefore divide smoothly by ea. At the edge,
\(T_0/e_a=F p_{s,r}\ne0\), with the stated positive cone direction.
Analyticity in eta follows from the core and forward formulas; the cutoffs
depend only on y. A smaller common complex neighborhood excludes their
nonzero denominators. Later edits start strictly after this collar.

For the rest of [X0,Xb], the error is below 6 Dc delta <.01 and
\(v_s\le\kappa_0C^3+.01<1\). At Xb, first taper the axial prescription
to zero and then interpolate a to .8. During the first taper the reference
Pc is \(p_{1,r}+\beta p_{2,r}^2/p_{1,r}>3\). During the second, bs=0
and Pc=p1>2. The added field changes are bounded by the same width budget.

Afterward a=.8 and U is constant in X. Its parameter derivative stays
close to the native one, since the new -.4 logarithmic phi slope is
parameter independent. The source exceeds 1:
\(.6\cdot2.8-.1-.02>1\). At a hypothetical downward crossing p1=2,
\(D_Xp_1=X S_q/L-.6p_1>X-1.2>0\). Thus p1>2 up to Xi.
The actual endpoint satisfies

\[
\|G_i-U_{nat}(4/\Lambda)\|_{C^1}\le6B^2\delta<\epsilon/100,
\quad \|\log(CE(X_i))\|_3<B.
\]

## 5. Shape transition, core moments, and the earlier exact join

Use (B.34) with the fixed Tsh and the actual endpoint functions ell_i,Gi.
Since \(\|\sigma'\|_\infty\le64\),
\(64(B+1)/(2000B)<.05\). Therefore \(.55\le l\le.65\),
\(.7\le a\le.9\), and bs=0 throughout the shape transition.

The old angular gradient is
\(-H_c\ell_i'\ge .95\Lambda L\chi-o(1)\), where the negative bound
comes from the explicit source absorption and width budgets above.
The new gradient obeys
\(-H_c(\log f)'\ge-2j^2-2K/\Lambda-o(1)\), by completing the square
in \(\eta H_*=(D+4d)\eta^2+dj\eta\).
Their convex combination, \(-W>2.8-o(1)\), and l>=.55 give Sq>1.
The same barrier preserves p1>2 throughout this transition and the following
ideal-power interval. In these source comparisons, the absorption and
continuation budgets are each below one millionth. The combined negative
errors are below .02 by the explicit parameter inequalities;
o(1) is descriptive shorthand, not a new parameter choice.

Set \(X_{sep}=110e^{T_{sh}}\) and \(X_R=110(CP)^{10}\). Hence
\(x_{sep}=e^{T_{sh}}/(CP)^{10}<e^{-8}\), with XR>10000.
On [0,Xsep], \(\|U\|_1<10\); in the shape transition
\(\|CE\|_1\le e^{203B}\). Near the axis CE is sqrt(X) times a smooth
function bounded in this norm by exp(B). Integrating the five original
densities, rather than inferring moments from the small radius, gives

\[
\|\widehat M\|_1\le10x_{sep},\quad
\|\widehat S\|_1\le101x_{sep},
\]
\[
\|\widehat I\|_1+\|\widehat J\|_1\le e^{210B}C^{-1}x_{sep}^{3/2},
\qquad \|C_p\|_1\le e^{410B}C^{-2}.
\]

For Cp the axis factor makes the integral finite. On the rest of the core
the logarithmic interval contributes less than exp(3B); on the shape
transition its length and the CE bound give less than exp(410B).
Subtract the ideal primitives and apply exactly the normalized row
combinations of the earlier annulus. Their denominators involve fixed
xc=exp(-6), P>=16 and f>=1/2. All ideal powers vanish at zero. Uniformly,

\[
\max_i\|r_i(X_{sep})\|_{C^1}
\le e^{1000B}(x_{sep}^{1/5}+C^{-2})
\le 2e^{1400B}C^{-2}<\epsilon/4.
\]

Also \(\|G_i-4\eta\|_{C^1}\le j+\epsilon/100+\epsilon/100
=3\epsilon/100<\epsilon/4\). The positive integrals between Xsep and
the fixed entry y=-8 give the previous annulus's entry ball. Its existing
exact five-moment correction can then be applied. This uses actual incoming
integrals of the constructed fields; the numerical manufactured fixtures
of the earlier checkpoint are not used to establish this implication.

Specifically, the five entry bounds are no larger than

\[
\epsilon/4+e^{-2}g_b,\quad
\epsilon/4+(5/8)e^{-16/5}g_b,\quad
\epsilon/4,\quad
\epsilon/4+(12/256)e^{-4/5}g_b^2,\quad \epsilon/4,
\qquad g_b=3\epsilon/100.
\]

Every row is below \(2.55\,10^{-17}\). The cross terms in the fourth row
cancel because \((4\eta+g)^2-8\eta(4\eta+g)=-16\eta^2+g^2\).
This calculation is fed into the earlier annulus's unchanged acceptance
function, using uniform upper bounds rather than manufactured samples.

## 6. What the checks can and cannot establish

The Banach contraction is an analytic existence argument. The finite
coefficient-index checks exercise its weight algebra, but do not replace
the all-index proof. Interval arithmetic checks every displayed parameter
inequality at P16, Md4 and Md64. The pressure shape remains the exact
positive-mixture datum, so those checks do not certify a quadrature fit.
The two independently assembled polynomial identities check the inner
sources and the stress cancellation. They cannot prove a missing analytic
estimate. The continuation bounds above are the mathematical argument;
small numbers alone cannot establish them. They remain open to independent
review.

Large constants are intentional overestimates, not estimates of the
physical size of a realizable experiment or of proximity to a blow-up proof.
In particular even degree 24 has no useful tail certificate from the
generic exact-core norm here. This prevents re-labelling the earlier finite
Taylor run as the newly proved analytic fixed point.

No claim about the full corrected PDE, admissible smooth force, or actual
finite-time singularity follows from this checkpoint alone.

Two limitations of earlier checks remain explicit. A constant .01 error
bound alone cannot preserve the sign of a source of size j=1e-18; the
reference-source checks now use j/100. Also, the generic coefficient norm
does not make a degree-24 tail small. Neither control is represented as a
counterexample to the actual profile or to the manuscript.

Reproduce from the repository root:

```sh
python research/navier_stokes_cascade/results/axis_core_attachment/test_audit.py -v
python research/navier_stokes_cascade/results/axis_core_attachment/audit.py \
  --verify-record research/navier_stokes_cascade/results/axis_core_attachment/evidence.json \
  --output axis_core_attachment_reproduced.json
```

The next obligation is Appendix C: convert the strict relaxed stress on the
remaining region into the fully admissible stress needed by the wave
construction, preserving the protected inner collar and matching moments.
Only after that can the full momentum corrections and forcing argument be
completed and tested.
