# Background jets: derivative conventions, axis values and the coefficient inverse

Status: `background-jet-construction-audited-relative-to-leading-profile`.

This audit derives the normalized background-jet bound required by equation
(23) of the [pulse replay](../geometric_pulse_replay/README.md) from an explicit
leading-profile input, rather than assuming estimates for the summed field.
It checks the regular axis expressions and the all-order coefficient inverse,
then specifies the seminorms that the cutoff selection must control. In
particular, Euler and frozen scale derivatives have different constants.

The remaining input is the finalized LEADING profile: its smoothness and
analytic-parameter collar, its reserved power-law patch and exterior
identities. Its strict growth/cone gaps and weighted stress edges also remain
to be independently verified. Thus this is not an unconditional pulse bound,
complete residual improvement, or blow-up result. The finite diagnostics
below do not evaluate the actual enormous-frequency profile.

Parent: `94e1f3376d6813376d5369bba9c3b949a9802dd4`. This follows the
[external review](../../EXTERNAL_REVIEW_2026-10-15.md). The defining equations
are (5.2)--(5.8), (5.27), (5.42)--(5.46), and (6.10)--(6.11) in the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its conclusion and its formalization are not premises. All earlier checkpoint
files and GOAL.md remain unchanged.

## 1. The input that is still needed

Keep the SAME fixed h in (0,1/100), finite leading frequency N, and profile.
Write A=1/2+h, D=1/2-h, L=1-2h eta^2, d=1-eta^2 and Rprof=sqrt(2X).
Separate the following leading-profile certificate from its consequences:

1. F0=E0/Rprof, U0 and Pi0 are smooth on every required compact closed
   profile rectangle, including eta=+/-1. The radial primitive M0 and V0/X
   have their regular axis extensions. On one interval
   0<=X<=a^2=Xa exp(tc/4), all fixed radial derivatives extend holomorphically
   in eta to a neighborhood of [-1,1], uniformly in X. The neighborhood
   may depend on the derivative order. This interval precedes all moment
   patches and modulation.
2. The retained cutoff has the exact ordered radii Xkeep=Xa exp(tc/8),
   Xcut=Xa exp(3tc/16), a^2=Xa exp(tc/4). On I3, with X=rho3 exp(y),
   0<y<5, the actual leading fields satisfy U0=0 and
   E0=e(eta) exp(-(1/2+ell_patch)y), where e is smooth and positive on
   [-1,1] and 0<ell_patch<=.01. The five fixed nonnegative unit-mass bumps
   have centers 1,2 for U and .5,1.5,2.5 for E, each of width .1.
   The later, disjoint I4 is the reserved mean interval.
3. The leading identities U0=M0=V0=0 for X>=Xv hold, with smooth joins;
   the compact pressure-source property and finite integrals used below
   follow. Zero leading U on I3 alone does NOT imply zero M0 or V0 there.

These are properties of defining maps, not numerical tolerances. The
previous profile checkpoints propose their proofs. This batch does not
recheck all the nonlinear core contraction, joining, modulation, repair and
heat-exterior arguments that establish this certificate. A finite C3 bound
alone is insufficient to substitute for item 1 at arbitrary order.

The pulse's strict positive gaps and weighted leading/first-stress estimates
are additional inputs, not consequences of this regularity certificate.
In particular a small absolute perturbation need not preserve a cone at a
flat stress edge.

## 2. Axis quotients have prescribed limits

For positive order n put lambda_n=2nh. Reconstruct

\[
 M_n=X\int_0^1 U_n(tX,\eta)\,dt,\qquad
 v_n:=V_n/X=
 \frac{2\eta U_n-2(D+\lambda_n)\eta (M_n/X)-d\partial_\eta(M_n/X)}{L}.
 \tag{1}
\]

All derivatives of M_n/X are obtained under this fixed integral. At X=0
its value is U_n(0,eta), not zero by a convention about division. Hence

\[
 v_n(0,\eta)=
 \frac{2(A-\lambda_n)\eta U_n(0,\eta)-d\partial_\eta U_n(0,\eta)}{L}.
 \tag{2}
\]

It vanishes at positive order when the prescribed zero axis trace is used.
The leading trace generally does not vanish. Likewise, if
Z_k=[2k eta-2 eta X partial_X+d partial_eta]/L, then

\[
 Z_k(Xv)=X Z_{k-1}v,\qquad
 X^{-1}Z_{k-D}Z_k(Xv)=Z_{k-1-D}Z_{k-1}v.\tag{3}
\]

The regular formula for Omega_k/X contains

\[
 \frac{(1-\lambda_k)v_k+Xv_{k,X}+D\eta v_{k,\eta}}{L}
 -2\nu(2v_{k,X}+Xv_{k,XX})
 +\sum_{i+j=k}\{v_i(v_j/2+Xv_{j,X})+U_i Z_{\lambda_j-1}v_j\},
 \tag{4}
\]

and, for k>=1, the additional term
-nu Z_(lambda_(k-1)-1-D) Z_(lambda_(k-1)-1) v_(k-1).
These expressions give the axis values and all fixed derivatives. They
retain the radial axial-viscosity source starting in the order-two pressure.

A simple control is U0=eta, independent of X, with the divergence-compatible
v0 from (1). At X=eta=0, (4) gives Omega0/X=-1/2. Assigning zero to that
quotient would give the wrong first pressure slope: the zero-trace equation
Pi1_X=2F0 F1-Omega0/(2X) gives Pi1_X(0,0)=1/4. This is a smooth manufactured
axis control, not the actual leading solution. The repo already uses the
regular formula; no frozen numerical value is corrected by this observation.

## 3. Why every coefficient exists on the same radial collar

The exact source decomposition in
[background summation](../background_summation/README.md), sections 1--2,
gives W_n=(F_n,U_n,M_n/X-U_n,Pi_n,partial_xi F_n,partial_xi U_n), xi=sqrt(X):

\[
 \partial_\xi W_n+\xi^{-1}\operatorname{diag}(0,0,2,0,3,1)W_n
 =A_{0,n}W_n+A_1\partial_\eta W_n+f_n,\qquad W_n(0)=0.\tag{5}
\]

All lower orders must be globally finalized before forming f_n. By (1)--(4),
the sources have no singular axis division. The derivative matrix and every
eta derivative of it map the first four coordinates into the last two and
annihilate the last two. The diagonal Volterra kernels preserve these
subspaces. Thus A1(xi) D0 A1(s)=A1(xi) D0 partial_eta A1(s)=0 for every
intervening diagonal D0, even at different radii. This is stronger than
pointwise A1(xi)^2=0; the latter alone would not justify the argument.

For fixed n, intersect the finitely many complex neighborhoods of the lower
coefficients. The eta-independent cutoffs preserve analyticity there, and
the moment bumps lie to the right. On a smaller neighborhood all coefficients
and f_n have a finite bound C_n on the SAME radial interval [0,a]. Let
0<Delta<=1 be a further strip loss and z=2C_n a. In a Picard word of length
k at most p=ceil(k/2) derivative factors survive. Assign loss Delta/p to
each derivative; the ordered radial simplex has volume a^(k+1)/(k+1)!.
The terms are bounded by

\[
 \frac{z^{k+1}}{(k+1)!}\max(1,p/\Delta)^p,
 \qquad \sum_{k\ge0}\frac{z^{k+1}}{(k+1)!}\max(1,p/\Delta)^p
 \le (z+\beta)e^\beta,\quad\beta=e z^2/\Delta.\tag{6}
\]

For k=2j use j^j<=e^j j! and (2j+1)!>=(j!)^2; for k=2j+1 use
(2j+2)!>=((j+1)!)^2. The two sums are at most z exp(beta) and
exp(beta)-1. This proves convergence for ANY finite C_n; no bound uniform
in n or common complex width for all n is required. The same word estimate
with k radial integrations proves analytic-class uniqueness.

The nonsingular formula Gg=xi integral_0^1 t^c g(t xi)dt gives all fixed
xi derivatives on successively smaller parameter neighborhoods. The first
four components are even and the last two odd, by the parity of (5) and
uniqueness. Even smooth functions of xi give smooth functions of X>=0.
The source and coefficients are smooth in X on the whole collar, although
the flat activation need not be analytic in X. Consequently cutting F_n,U_n
between Xkeep and Xcut preserves the needed parameter analyticity on [0,a^2].

## 4. The same actual moment matrix remains invertible at every order

Forward reconstruction, not a cutoff of pressure, defines Pi_n. With lower
orders fixed, the five variations on I3 are exactly

\[
 (\delta m_1,\delta m_2,\delta m_3,\delta m_4,\delta m_5)=
 \left(\int\delta U,\int\sqrt{2X}\delta E,
 \int E_0\delta E/X,\int\sqrt{2X}E_0\delta U,-\int E_0\delta E\right).
 \tag{7}
\]

No current-order quadratic term occurs: its formal order is 2n, greater
than n. Lower-order products, radial memory, viscosity and the full outer
pressure-source integral belong to the incoming discrepancy.

Normalize the five moments by
(rho3,sqrt(2)rho3^(3/2),e,sqrt(2)rho3^(3/2)e,rho3 e), and reorder them as
(M,(M-J)/ell_patch,I,-S,Cp). The matrix is block diagonal. Its U block has
weights exp(y) and exp(y) w(y), where
w(y)=(1-exp(-ell_patch y))/ell_patch and w'(y)>0. Its determinant is an
integral of exp(y1+y2)[w(y2)-w(y1)] against two positive bumps with ordered
disjoint supports, hence is strictly positive.

Its E block has weights exp(3y/2), exp((1/2-ell_patch)y),
exp((-1/2-ell_patch)y). For three distinct real rates, a nonzero linear
combination of their exponentials has at most two distinct real zeros:
divide by the smallest exponential, differentiate, and apply Rolle's theorem
inductively. The evaluation determinant is therefore nonzero on ordered
triples and has constant sign there. Integrating against the three ordered
positive bumps preserves that nonzero sign. This proves invertibility for
the ACTUAL supports and every fixed positive ell_patch in the stated range,
without a sampled determinant.

There is also a direct quantitative bound for these exact supports. Write
x_j=exp(y_j), p=2+ell_patch. After factoring exp(-(1/2+ell_patch)y_j)
from the E columns, its determinant magnitude is the Vandermonde product
times the second divided difference of x^p. The latter is at least one:
f''(x)/2=p(p-1)x^(p-2)/2>=1 for x>=1. Unit bump mass therefore gives

\[
 |\det B_E|\ge d_E:=e^{-.51(4.65)}
 (e^{1.45}-e^{.55})(e^{2.45}-e^{.55})(e^{2.45}-e^{1.55})>0.\tag{7a}
\]

Every E entry is at most exp(1.5*2.55), so the cofactor formula gives
||B_E^-1||_infinity<=6 exp(3*2.55)/d_E<1000. Similarly
det B_U>=.9 exp(2.9-.01*2.05), while each U entry is at most
2.05 exp(2.05), giving ||B_U^-1||_infinity<2. The displayed numerical
inequalities are checked with outward interval arithmetic. This proves
the previously used normalized bound 1000 on the whole positive patch-
exponent range; physical row normalizations remain separate.

The normalized matrix is independent of eta and n. Its inverse is finite.
The positive e has a positive minimum on the closed eta interval. Thus
c_n(eta)=-B^{-1}z_n(eta) is smooth for every finite incoming discrepancy;
there is no small-data condition for positive orders. The normalization
still costs ell_patch^-1 for arbitrary M-J, and e^-1; neither is discarded.

The zero first moment makes M_n=V_n=0 after the velocity bumps. The zero
third moment makes Pi_n=0 after Xv, where its source is zero. Inductively,
all positive-order coefficients have a common compact support. Pressure
can be nonzero on I4; cutting it off there is invalid. Since the bumps are
right of a^2, they leave (5) and its parameter regularity intact. This
finishes the coefficient construction needed for VELOCITY jets, relative
to the leading certificate. Stress tails still require their separate
conservative moment and heat-edge arguments.

## 5. Convert the derivative convention before selecting cutoffs

Let E=q partial_q and Q_j=q^j partial_q^j, acting at fixed X,eta. The exact
operator identities are

\[
 E^m=\sum_{j=0}^m S(m,j)Q_j,\qquad
 Q_m=\sum_{j=0}^m s(m,j)E^j,\tag{8}
\]

where S and s are Stirling numbers of the second and signed first kinds.
This follows by E Q_j=j Q_j+Q_(j+1), or by applying both sides to monomials.
With all derivative orders through m included, the norm conversion costs
at most B_m=sum_j S(m,j) (a Bell number) in the first direction and m! in
the second. The mixed X,eta derivatives commute with both families.

For example E^2 f=Q_2 f+Q_1 f. A frozen-derivative bound 2^-n does NOT
transfer with that identical constant; through m it transfers with B_m 2^-n.
To insist on the same numerical diagonal budget, impose 2^-n/B_n on all
frozen jets through n. Alternatively select the cutoffs directly in Euler
seminorms, as follows. Both choices are finite at each step and cost no
extra q power. This makes explicit a constant which must not be omitted
when translating a formal or computational estimate.

Keep the particular Stokes potentials from (5.27); do not change their gauge.
For sigma=a_n q, chi=1 on sigma<=1/2 and chi=0 on sigma>=1, the exact
normalized increments are

\[
 q^{2nh}\chi(\sigma) E_n,\quad q^{2nh}\chi(\sigma) U_n,\quad
 q^{2nh}\left[\chi(\sigma)W_n^{r}
             -(E\chi)(\sigma)J_n^{r}\right],\tag{9}
\]

where W_n^r=V_n/Rprof and J_n^r=2 eta M_n/(L Rprof). These symbols are
scalars, distinct from the six-vector in (5). The last expression is the
increment of q^(1/2) ur. It retains the cutoff curl commutator.

On a fixed enlarged annulus Xlo<=X<=Xhi, Xlo>0, let M_(n,m) bound all
X,eta derivatives through m of E_n,U_n,W_n^r,J_n^r. Sections 2--4 give
FINITE such constants at every fixed n,m. Put
K_j=max(1,max_(i<=j) sup_sigma |(sigma partial_sigma)^i chi|). Then

\[
 C_{n,m}:=\max(1,2K_{m+1}(1+2nh)^m M_{n,m})
 \quad\Longrightarrow\quad
 |E^j\partial_X^a\partial_\eta^b(\text{increment})|
 \le C_{n,m}q^{2nh},\quad j+a+b\le m.\tag{10}
\]

The extra cutoff derivative m+1 in the radial entry is necessary. K_j is
independent of a_n. If cutoff bounds were supplied using frozen derivatives,
(8) gives K_j<=max_(i<=j) B_i max(1,||sigma^r chi^(r)||_(r<=i)).

An explicit recursive choice, for q0<=1, is

\[
 \log a_n\ge 1+\max\left\{\log a_{n-1}+\log2,
 \log(q_0^{-1})+\log2,
 \max_{m\le n}\frac{\log C_{n,m}+n\log2}{nh}\right\}.\tag{11}
\]

Combine this finite list with the inherited physical and weighted-stress
requirements by taking their maximum. On the ENTIRE interval q<=a_n^-1,
(10) is then at most 2^-n q^(nh) for m<=n. There is no endpoint-only
sampling and no bound uniform in n is assumed. The additional constants
can be enormous. This specifies admissible choices in the already free
cutoff schedule; it does not claim that actual numerical a_n were evaluated.

For fixed m choose J>=max(m,1). Below q<(2a_J)^-1 the initial J cutoffs
are one. Their finite sum is O_m(q^(2h)); the tail is at most
2^-J q^(h(J+1)), also O(q^(2h)). All fixed normalized velocity jets
therefore have the claimed q^(2h) error. On the remaining compact q
interval up to q0, local finiteness and closed-rectangle smoothness enlarge
the constant. This gives ONE common q domain for every fixed jet, with
constants depending on the jet order. It does not require one uniform
constant over all orders. The radial leading term gives q^A ur=O(q^h).

## 6. Pull back to the actual band, including its one-sided boundary

For Q=2^-ell set s=q/Q, X=R^2/(2s), eta=Zs^-D. The implicit relation is
s-Z^2 s^(2h)=T; on 1/2<=s<=2, |eta|<=1, L>=.98. On a normalized scalar
f(q,X,eta), the exact chart derivatives are

\[
 \partial_R=\frac R s\partial_X,\qquad
 \partial_Z=\frac{s^{-D}}L(2\eta E-2\eta X\partial_X+d\partial_\eta),\qquad
 \partial_T=\frac1{sL}(E-X\partial_X-D\eta\partial_\eta).\tag{12}
\]

Repeated derivatives also differentiate these coefficients. All their fixed
jets are bounded on the same compact chart domain, by the implicit function
theorem and L>=.98. No Q power occurs in (12). At eta=+/-1, T=0, the
coefficient d vanishes but L remains positive; all one-sided jets remain
finite. The profile extension, not a zero-extension in T, defines those jets.

For b=Q^A ur, F=Q^A utheta/R and G=Q^A uz, (9)--(12) give

\[
 b=Q^h s^{-1/2}(V_0/Rprof+O(q^{2h})),\quad
 F-F^{(0)}=s^{-A}O(q^{2h})/R,\quad
 G-G^{(0)}=s^{-A}O(q^{2h}).\tag{13}
\]

Here R is bounded away from zero on the fixed enlarged annulus. Every
fixed slow derivative obeys the same Q orders, since q^(2h)=Q^(2h)s^(2h).
Thus the left side of pulse equation (23) is finite, RELATIVE TO the leading
regularity/support certificate in section 1 and the cutoff selection (11).
Auxiliary-label derivatives of this axisymmetric background are zero;
label representatives do not introduce new constants.

## 7. What this changes, and the precise remaining gate

The normalized background estimate no longer needs to be taken as a separate
black-box assumption once the stated leading certificate is proved. Its
all-order inverse, axis limits, derivative conversion, cutoff commutator and
common-domain passage are supplied above. No change to the old formulas was
needed. The conversion of frozen to Euler jets requires the displayed Bell
factor; declaring the operators equal would be a genuine error.

The next audit must verify the actual leading defining maps on this collar
and the pulse's STRICT value and edge bounds at the same N,h. In particular,
for the finalized leading representatives, put N=g0/|g0|,
K=(-N_z,N_theta), lambda0 equal to the positive growth root, and
c0=lambda0/(2F0 N_theta). With That0 the smooth limiting unit stress
direction, the required inequalities include

\[
 \inf F_0>0,\quad \inf|g_0|>0,\quad
 \inf\{-2F_0 N_\theta(2F_0N_\theta+|g_0|)\}>0,\\
 \inf(-\widehat T_0\cdot N)>0,\qquad
 \sup\left|c_0\frac{\widehat T_0\cdot K}{\widehat T_0\cdot N}\right|<1.\tag{14}
\]

The orientation inequality selects the positive covariance cone; an absolute
ratio alone would also accept its negative. The required weighted T0,T1
bounds are separate. Neither taking q smaller nor the algebra above proves
(14) for an unverified leading field. The nonlinear core/continuation and
final cone certificate remain the next foundation work, followed by the
reserved physical mean-pressure bump and localization assembly checks.

The diagnostics test (8) against independently differentiated fractional
monomials and cutoff polynomials, (12) against implicit-coordinate
differentiation, the nonzero axis limit, actual bump ordering and matrix
refinement, and failure controls for dropped terms. They do not turn finite
data into continuum inequalities. Existing background tests are also rerun
to check the source and reconstruction algebra on which sections 2--4 rely.

```sh
python3 research/navier_stokes_cascade/results/background_jet_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/background_jet_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/background_jet_foundation/evidence.json \
  --output background_jet_foundation_recheck.json
```
