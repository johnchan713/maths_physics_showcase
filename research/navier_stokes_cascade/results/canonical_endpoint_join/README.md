# Canonical punctured terminal jets and sharp exterior Gevrey order

Parent: `4b9bdd40944f851fc0aebc10fb47bf22f30b4e37`.
Status: `conditional-canonical-punctured-endpoint-join-and-sharp-exterior-Gevrey-two`.
Parent hosted CI was still running without a failed step when inspected
on 18 October 2026. Historical fields and frozen evidence are unchanged.

The SAME [canonical local realization](../canonical_local_realization/README.md)
has compatible one-sided smooth potential, direct-swirl, velocity and
pressure jets at t=1 on every punctured terminal compact in a fixed late
observation cylinder. This is relative to the retained exact-core/background
and finite analytical defining/existence/inverse/gap premises. It is not an
unconditional regularity, global forcing or blow-up theorem. No enormous
actual coefficient, cutoff schedule or terminal field is numerically evaluated.

We match the canonical family, not the older realized fields, to the
reviewed [endpoint mechanism](../endpoint_foundation/README.md). A new sharp
estimate identifies Gevrey order TWO for the unchanged nonzero exterior
heat and FULL pressure jets. This does not assert that order for the full
canonical interior. A separate exact-residual counterexample explains why
flat residual and bounded velocity alone cannot supply terminal traces.

The retained manuscript Sections 6,9,10 and positive heat integral are
objects of review, not final theorems used as inputs. PDF SHA256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.

## 1. SAME canonical representatives up to the closed chart boundary

Keep Md64, viscosity one, original h,N, fixed positive lambda, radius,
width C_setup^-20, all frozen carriers/labels/supports, common qstar and
BOTH existing diagonal schedules. Set A=1/2+h, D=1/2-h, tau=1-t,

\[
q=\tau+z^2q^{2h},\quad\eta=z/q^D,\quad X=r^2/(2q),\quad
L=1-2h\eta^2\ge1-2h>49/50,\qquad0<h<1/100.                     \tag{1}
\]

At tau=0,z!=0, q=|z|^(1/D), eta=sign(z), q_t=-1/(1-2h).
There is no inverse 1-eta^2. With the actual canonical seed and increments,

\[
\mathcal A_c=A_{c,0}+\sum_{j\ge1}\chi(a_jq)A_{c,j},\quad
\mathcal B_c=B_{c,0}+\sum_{j\ge1}\chi(a_jq)b_{c,j}e_\theta,
\quad u_c=\operatorname{curl}\mathcal A_c+\mathcal B_c,
\quad p_c=p_c^{[0]}+\sum_{j\ge1}\chi(a_jq)(p_c^{[j]}-p_c^{[j-1]}).
                                                                    \tag{2}
\]

The pressure gauge and operator-defined potentials are EXACTLY the ones
fixed before summation in the parent. No endpoint gauge reset is permitted.
Every curl commutator, full six-entry mean covariance, axial velocity and
moving pressure moment is present. Both exceptional base tails stay in
p_c^[0]; only SAME-base changing parts telescope. Direct swirl remains
separately axisymmetric, not an invented compact full vector potential.

For the leading core let M=int_0^X U, m=M/X. Its regular Cartesian factors
are

\[
A_\theta e_\theta=\tfrac12q^{-A}m(X,\eta)(-y,x,0),\quad
B=q^{-A-1/2}F(X,\eta)(-y,x,0),\quad p=q^{-2A}\Pi(X,\eta).          \tag{3}
\]

The complete curl gives u_z=q^-A U and
r*u_r=[2eta XU-2Deta M-(1-eta^2)M_eta]/L. The reviewed core/background
estimates control every fixed jet of the regular m, not a bare 1/E or
fractional auxiliary map at r=0. On q>=c>0, all annular corrections have
a uniform axis hole r^2<2Xa*c. Thus (3) and the summed background handle
the axis; auxiliary pullbacks are only used at positive radius.

The physical charts are the enlarged CLOSED 1/4<=q/Q<=4,T=tau/Q>=0
charts, with discrete choices frozen before taking limits. Br, its normal
projection and left inverse are smooth functions of the closed-chart
coefficients on the SAME positive value floors. Differentiating the
canonical finite-interval primary ODE, including its canonical entrance,
gives finite one-sided parameter jets by Duhamel induction. No division
by T or 1-eta^2 occurs. Its own Hc integrals, inverse and positive ac roots
have the same fixed-jet bounds; weighted flat targets control root edges.
This is not the old primary or its old roots renamed.

The zero-entry sourced ODE, fixed-root signed map, shifted radial primitive,
zero-Haar temporal inverse and fixed slow five-map preserve these one-sided
jets, with finitely many input derivatives at each fixed stage. All source
and moment/pressure targets are recomputed in the canonical family. Full
moving-normal/source-normal pressure terms and actual compact A0/A1 pieces
remain. Thus the parent's canonical expression induction is also a closed
T>=0 parameter induction, under the retained analytical input estimates.
No arbitrary new blanket endpoint regularity premise is inserted.

## 2. Finite atlas and the precise bound needed for traces

For c<=q<=qmax, every accepted enlarged dyadic chart has

\[
\max\{\ell_0,\lceil\log_2(1/(4qmax))\rceil\}
 \le\ell\le\lfloor\log_2(4/c)\rfloor.                            \tag{4}
\]

Each band has finitely many mesh boxes on a compact slow range. Closed
mesh endpoints keep the preselected T=0 indices -1,0,1. Frequencies,
colors and common lifts are not reselected at the terminal time. The
same absolute-Y pullbacks, ODE uniqueness, physical radial primitives,
and q-based mean bumps identify the formulas and zero-extension jets.
For (4) at c=2^-24,qmax=2^-20 the enlarged window is 18,...,26;
copying the narrower historical window 19,...,25 would miss its edges.

Because the canonical scale recipe gives a_j>2^j/qstar, on q>=c an
active correction stage must satisfy j<log2(qstar/c). The independent
background schedule is locally finite there too. Finite band/mesh/stage
sets and all-fixed-order closed physical estimates therefore imply

\[
\sup_{E^+\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(\mathcal A_c,\mathcal B_c,p_c)|<\infty
 \quad\hbox{for every fixed }\alpha,j,\ \hbox{including }j+1.      \tag{5}
\]

This is the precise blocking inequality if any retained closed-chart
profile, root or inverse estimate fails. Interior smoothness, bounded
velocity and origin-flat residual do NOT replace (5). Constants may depend
arbitrarily on the fixed compact, band, stage and derivative order.

For a domain with explicit margins choose
q_e=min(qstar/16,2^(-ell0-2)), tau<=q_e/16, |z|<=(q_e/16)^D.
The implicit equation gives q<=2tau+4|z|^(1/D)<=3q_e/8. Since Q<=4q,
4Q<=16q<=6q_e<=3qstar/8 and Q<=3*2^-ell0/8. This restricts only the
observation cylinder, not parameters, labels or the actual field.

To cover a terminal compact |x|>=d>0, take
b=min(q_e/4,d^2/(16Xext),1/4), zeta=min(d/4,(b/16)^D), tau<=b/16.
On |z|>=zeta/2 there is a positive q lower bound. On |z|<zeta,
r>d/2, q<=3b/8 and X>=d^2/(3b)>=16Xext/3. That ENTIRE equatorial
neighborhood is the unchanged heat exterior. The overlapping sets cover
axis points with z!=0 and positive-radius equatorial points, not the origin.
The full global pressure moments and compact changing representatives,
not velocity support alone, justify the same exterior pressure.

## 3. Complete exterior mixed jets; a new sharp order estimate

Keep its fixed gauge and constant C:

\[
K=C r^{-1-2h}H(4\tau/r^2),\quad
H(Z)=\Gamma(1+h)^{-1}\int_0^\infty e^{-v}v^h(1+Zv)^{-h}dv,
\quad p=-\int_r^\infty K(\rho,\tau)^2\,d\rho/\rho.              \tag{6}
\]

For F_j=(-1)^j H^(j), A_j=(h)_j(1+h)_j, positive domination gives
0<=F_j<=A_j on Z>=0. Set
G_j=sum binomial(j,k)F_k F_(j-k), S_j=sum binomial(j,k)A_k A_(j-k).
Then G_j'=-G_(j+1); ALL Leibniz products stay. Pressure differentiation
is dominated by an integrable rho^(-3-4h-2j) density. The reviewed radial
recursion c_(n+1,k)(a)=-(a+n+2k)c_(n,k)(a)+2c_(n,k-1)(a) yields all fixed
mixed radial/time bounds on r>=rmin>0,0<=tau<=tau0, with
Zmax=4tau0/rmin^2. In particular,

\[
|\partial_t^jK|\le |C|4^j A_j r^{-1-2h-2j},\qquad
|\partial_t^jp|\le C^2 4^j S_j r^{-2-4h-2j}/(2+4h+2j),            \tag{7}
\]

and equality in absolute value holds at tau=0. The time sign is positive
for K and negative for p. No positive-Z integral is used at negative Z.

For EVERY j>=1, elementary product bounds give

\[
\frac h j(j!)^2\le A_j\le h(j+1)(j!)^2\le2^j(j!)^2,\qquad
2A_j\le S_j\le3\,2^j(j!)^2.                                    \tag{8}
\]

Indeed (h)_j lies between h(j-1)! and h*j!, while (1+h)_j lies between
j! and (j+1)!. In the pressure convolution use A_k<=2^k(k!)^2, also at
k=0, so each summand is at most 2^j(j!)^2/binomial(j,k). The reciprocal
binomial sum is at most 3: endpoints give 2, and the j-1 interior terms
are at most 1/j each. The j=0,1 cases are immediate.

Thus pure temporal derivatives of K and p obey C_E(8/rmin^2)^j(j!)^2
on every such exterior compact. For fixed radial n, the recurrence
coefficients have degree at most n in a=1+2h+2j. Also
(j+k)!/j! <= max(1,n)^n(j+1)^n for k<=n. The full mixed bound is at
most C_nE*(8/rmin^2)^j*(j+1)^(3n)*(j!)^2, hence at most
C_nE*[2^(3n+3)/rmin^2]^j*(j!)^2 using j+1<=2^j. Pressure uses degree
n-1 for n>=1 and the same coarser bound suffices. Cartesian/frame
derivatives add only fixed inverse-radius costs; axial derivatives vanish.

The order TWO is sharp at any fixed r>0,h>0 when C!=0. At tau=0,
the lower bounds in (8) give K jets bounded below by a positive multiple
of (4/r^2)^j(j!)^2/j, and pressure jets by a positive multiple of
(4/r^2)^j(j!)^2/[j(2+4h+2j)]. These cannot be bounded by M R^j(j!)^s
for ANY fixed M,R and s<2: even j!>=(j/2)^(floor(j/2)) defeats every
exponential after taking the power 2-s>0 and dividing by a polynomial.
This proves sharp Gevrey-2 terminal growth, not just zero analytic radius.
It concerns ONLY the unchanged exterior; arbitrary canonical interior jet
constants do not imply a Gevrey bound. Later force gluing must be smooth,
not analytic heat continuation or pasting nonzero pressure jets to zero.

## 4. Compatible terminal traces and an exact-residual obstruction

Equation (5), including the NEXT time derivative, gives the uniform Cauchy
estimate ||D_x^alpha partial_t^jG(t2)-D_x^alpha partial_t^jG(t1)||
<=M_(alpha,j+1)|t2-t1| for G=(Ac,Bc,pc). Uniform limits in enlarged
spatial coordinate integrals identify spatial derivatives; time integration
identifies successive time jets. For each fixed finite k, with terminal
trace T_j^G,

\[
\left\|D_x^\alpha\left[G(1-\tau)-\sum_{j=0}^k
 (-\tau)^jT_j^G/j!\right]\right\|
 \le M_{\alpha,k+1}\tau^{k+1}/(k+1)!.                            \tag{9}
\]

Curl commutes with these finite limits. A velocity m-jet Cauchy estimate
needs potential order m+2 and direct swirl m+1; pressure needs m+1.
No uncut infinite series or uniform-in-order constant is invoked. Actual
flat velocities and the finite Gaussian residual ledger remain untouched.

The necessity of the endpoint bounds is not merely formal. On t<1 take

\[
u=\sin(1/(1-t))e_1,\qquad
p=-x_1\cos(1/(1-t))/(1-t)^2.                                    \tag{10}
\]

This is smooth, divergence-free, uniformly bounded in velocity, and has
EXACTLY zero full Navier--Stokes residual: advection and viscosity vanish,
and grad p=-u_t. It even has the bounded-on-compacts potential
A=(0,0,x2*sin(1/(1-t))). Yet along tau_n=(pi/2+2pi*n)^-1 and
sigma_n=(3pi/2+2pi*n)^-1 the velocity equals +e1 and -e1, so there is
no terminal trace. It violates (5); it is neither the canonical field nor
a compact finite-energy counterexample. Exact differentiated quadrature
identities check the cancellation and chain-rule signs. The two sequences
prove nonconvergence analytically, not by time sampling.

## 5. What is closed and what remains

The canonical one-field punctured endpoint implication is closed relative
to the retained core/background and finite analytical premises. It asserts
no regular velocity or pressure at the singular origin and no future
velocity. New gauges, labels or diagonal schedules are not introduced.

Next join axisymmetric spatial/time localization and zero initialization
to THIS canonical field, retaining all potential curls, direct swirl,
pressure, viscosity and nonlinear interactions. Then construct the
admissible smooth compact force, whole-space energy/dissipation and exact
global comparison. Older global fields are not silently identified.
Independent discharge of the retained exact-core/background and finite
defining/existence/inverse/gap assumptions remains before unconditional
pulse, complete-residual, forcing or blow-up claims.

Exact rational checks cover the canonical matching, enlarged atlas,
factorial bounds, complete pressure convolution, finite trace requirements
and oscillatory cancellation. Reused independent positive Gamma integrals
test mixed heat/pressure identities at diagnostic parameters. Those finite
controls are not continuum proofs or evaluations of actual terminal fields;
the scoped proofs are the derivations above.

```sh
python3 research/navier_stokes_cascade/results/canonical_endpoint_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/canonical_endpoint_join/audit.py \
 --verify-record research/navier_stokes_cascade/results/canonical_endpoint_join/evidence.json \
 --output /tmp/canonical_endpoint_join_fresh.json
```
