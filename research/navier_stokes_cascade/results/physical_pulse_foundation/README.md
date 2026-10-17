# Actual primary pulses: one domain, all physical jets, complete curls and tails

Parent: `585342bc4225cd69cf07a0a9a24d65ef01729804`.
The [weighted-stress foundation](../weighted_stress_foundation/README.md)
now supplies the actual background inputs previously separated in
[physical pulse transfer](../physical_pulse_transfer/README.md). This
checkpoint independently reviews the phase/frame, radius, graph derivative
and support arguments and joins them to those inputs. It derives uniform
PRIMARY pulse bounds on every used physical scale, relative to the retained
exact pressure-family core and finite existence/inverse/gap lemmas. These
remain explicit dependencies, not numerical assertions or the manuscript's
final theorem. All historical sources and evidence stay unchanged.

The scope is important. The positive lower bound is for the uncut radial
ODE amplitude, NOT the complete oscillating velocity. The sourced inverse
is a mapping for inputs satisfying its stated coefficient/support bounds;
those bounds are NOT proved here for every source in the nonlinear hierarchy.
Complete-wave covariance mean corrections and the complete residual still
require independent review. No unconditional all-pulse-family, residual,
forcing, comparison or blow-up theorem is promoted by this checkpoint.

Sources: (6.1)--(6.23), (7.2)--(7.29), (7.37)--(7.40), (9.19) in the
retained [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The exact geometry and radius costs are recorded in
[constructive supports](../auxiliary_geometry_audit/README.md) and
[geometric pulse replay](../geometric_pulse_replay/README.md). We use unit
viscosity and keep Md64,h,N,t1=C^-20 and the existing supports/carriers.

## 1. One actual background/chart domain, with frozen discrete choices

Keep A=1/2+h, D=1/2-h, kappa=10^-5. For a band ell let
Q=2^-ell, epsilon=Q^h, S=ell^2, k=ceil(epsilon^-1/2). In its physical chart

\[
 R=r/Q^{1/2},\quad Z=z/Q^D,\quad T=(1-t)/Q,\quad s=q/Q,
 \qquad s-Z^2s^{2h}=T,\quad X=R^2/(2s),\quad eta=Zs^{-D}. \tag{1}
\]

The cancellation 2D+2h=1 removes Q from the implicit equation.
On the closed physical band 1/4<=s<=4, T>=0, its selected positive branch
has derivative L=1-2h eta^2>=.98. Uniform implicit differentiation gives
every fixed jet of s,X,eta on a fixed compact radial enlargement, without
a Q loss. At T=0 these are one-sided jets; the Gamma heat formula is never
evaluated at negative Zheat. No global analytic heat extension is assumed.

The parent gives every fixed Euler/profile jet on one q-domain. Consequently,
with b=Q^A ur,B, F=Q^A utheta,B/R, G=Q^A uz,B,

\[
 |D^I b|<=C_I epsilon,\quad
 |D^I(F-F^{(0)},G-G^{(0)})|<=C_I epsilon^2.             \tag{2}
\]

Here F^(0) and its radial shear both acquire s^(-A-1/2) from the leading
profile, so their direction, c0 and strict growth/cone comparisons transfer.
These are actual inputs now supplied relative to the retained core/finite
lemmas, not assumptions inferred from the manufactured wave fixture.
Representatives are chosen in each box's intersection with the CLOSED
active chart, not at arbitrary exterior mesh centers. Continuity of the
nondegenerate leading quantities gives a fixed small radial enlargement;
the limiting stress direction is used only on the closed active shell.
All doubled slow boxes lie in that enlargement for sufficiently large ell.

Fix the constructive radius r0=1/6057072, Tg=4+sqrt(2),
Lambda_g=4-sqrt(2), rho_g=log(Lambda_g)/log(Tg),
dr=2[(1+h)rho_g-h kappa]. Exact flooring yields

\[
 (Tg S)^{-1}<c_i<=S^{-1},\quad
 M_i<=epsilon^{-kappa}S^{-rho_g},\quad
 L_s=2r0/c_i,\quad 2r0 S<=L_s<2r0 Tg S<S.          \tag{3}
\]

Take one ell0 large enough that 4Q<qstar, the doubled boxes stay away from
R=0 and in the fixed chart, and L_s exceeds the finite acceptance lengths
below. Require S^2(epsilon+epsilon^2+k^-1)<=1 and any padded finite
constant version needed for the value/frame gaps. Since 1/k<=sqrt(epsilon),
the left side is at most 3 ell^4 exp[-(h log2)ell/2]; the all-half-line
log x<=x-1 estimate supplies one finite threshold for ALL later bands.
No high derivative norm imposes a new ell0. Work on q<=2^(-ell0-2),
so every active dyadic label is accepted. The actual constants and first
integer are not numerically materialized. Representatives, signs, mesh,
cover, rectangle copies and integer carriers are frozen before derivatives.

## 2. Recheck the frame and the nonvanishing ODE pulse

At each leading representative define N=g0/|g0|, K=(-Nz,Ntheta),
lambda0^2=-2F0 Ntheta(2F0 Ntheta+|g0|)>0 and c0=lambda0/(2F0 Ntheta)<0.
Choose the SAME positive tilt u from the global strict leading cone.
Let Bs^2=lambda0/[epsilon k^2(1+u^2)^(3/2)],
s_sigma(v)=sigma*u*(1/2+v/L_s), x0=sigma*Bs*u/2, and
(ptilde/R0,pz)=Bs[K-sigma*u*g0/(L_s|g0|^2)]. Nearest NONZERO integer
rounding gives kp in Z minus {0} and |p-ptilde|<=1/k, not 1/(2k).
The phase and normal are exactly

\[
 Phi=p theta+pz Z/epsilon+x0R-v(pF+pzG),
\quad n=(x0-v(pF_R+pzG_R),p/R,pz-epsilon*v(pF_Z+pzG_Z)). \tag{4}
\]

Keeping the full interval v<=L_s and the fixed radius, direct subtraction
from Bs(s_sigma,K) gives errors bounded by
C[L_s^-1+(1+L_s)(S^-3+epsilon+epsilon^2+k^-1)]<=C/L_s.
The normal derivative has the same bound. The L_s/k rounding error cannot
be omitted. Also Eik=epsilon*v[(pF+pzG)_T-G(pF+pzG)_Z]+b*n_r retains
epsilon in every fixed coefficient derivative. Because epsilon*k^2 is in
[1,4], Bs stays bounded above and below. The floor is for |n_tan| and |n|,
NOT for |p/R|, which may tend to zero as k grows. No inverse p is used.

Set Ka=n_tan/|n_tan|, Na=((Ka)z,-(Ka)theta), sa=n_r/|n_tan|, and
B=[er-saKa,Na] [[1,1],[c0 sqrt(1+sa^2),-c0 sqrt(1+sa^2)]]. Its left
inverse and every fixed coefficient jet are bounded by powers of S. For
Kbg=[[0,-2F,0],[2F+R F_R,0,0],[G_R,0,0]],

\[
 A_n=-Kbg+{n(n^T Kbg-(n')^T)\over|n|^2},\quad
 Bleft(A_n B-B')=diag(lambda,-lambda)+E,
 \quad lambda={lambda0\over\sqrt{1+s_sigma^2}},\quad
 |E|<=C_E/L_s,\quad |d-dref|<=C_d/L_s,              \tag{5}
\]

where d=epsilon k^2|n|^2 and
dref=lambda0(1+s_sigma^2)/(1+u^2)^(3/2).
The reference algebra is checked without discarding n' or B': its
undamped coordinates are [[0,2F0Ntheta/(1+s^2)],[-(2F0Ntheta+|g0|),0]],
diagonalized by the two displayed columns. Nondegenerate smooth rational
dependence and (2),(4) bound the remaining terms in (5). For every fixed
j, coefficient jets in (5) cost at most C_j S^j, with no new denominator
smallness condition. This proves the frame premise rather than assuming it.

Let mu=min lambda>0. Pad C_E>0 and choose L_s>=8C_E/mu. In the m=1
homogeneous frame solution take zplus(0)=P(0), zminus(0)=0,
P'/P=lambda-dref, P(L_s/2)=1. For r=C_E/(mu L_s)<=1/8 the Riccati
equation for zminus/zplus has outward derivative at either boundary <=

\[
 -2mu r+(C_E/L_s)(1+2r+r^2)
 <=-{47\over64}{C_E\over L_s}<0.                  \tag{6}
\]

The logarithm of zplus/P has total variation <=C_d+C_E(1+r). Its scalar
exponential cannot reach zero, closing the continuation argument. Therefore
the radial amplitude x=zplus+zminus satisfies cP<=x<=CP on the WHOLE
closed pulse. Differentiating the reference growth rate gives positive
fixed Gaussian constants cG,CG, with

\[
 e^{-CG(v-L_s/2)^2/(2L_s)}<=P(v)
 <=e^{-cG(v-L_s/2)^2/(2L_s)}.                       \tag{7}
\]

The fixed-v differentiated Duhamel recurrence b0=0,
bj=max_(1<=i<=j)(i+1+b_(j-i))=2j supplies all slow jets relative to P.
Pulse derivatives follow from the ODE; ordinary auxiliary derivatives add
only polynomial S costs through c_i^-1. The m!=0 propagator has the exact
factor exp[-(m^2-1)int d] times the m=1 propagator. Thus a fixed finite
harmonic set needs no smaller domain. A sourced zero-data inverse retains
the source's epsilon exponent and sqrt(zeta) edge weight with finite
polynomial losses, PROVIDED the input has the stipulated bounds and smooth
support along the full lifted path. This is not a theorem about every
source of an unreviewed nonlinear iteration.

Pressure remains part of the principal system:

\[
 pi_m={i\over km|n|^2}(n\cdot Kbg*t_m-n'\cdot t_m+n\cdot f_m). \tag{8}
\]

The moving-normal and source-normal terms are necessary for compatibility
of n.t=0 with the projected evolution. Neither term is removed because
the velocity is transverse. For primary f=0, (8) gains 1/k.

## 3. Positive covariance coefficients with the actual Haar area

The fixed transverse cutoff and dv-to-Haar Jacobian contribute a factor
comparable to r0*c_i. For L_s>=25, integrate the lower bound (7) over
|v-L_s/2|<=sqrt(L_s), where the temporal cutoff equals one; integrate its
upper bound over the entire line. Along with cP<=x<=CP this gives

\[
 h_sigma\asymp r0^2/\sqrt{L_s}\asymp r0^{3/2} S^{-1/2}. \tag{9}
\]

The Gaussian first centered moment is O(sqrt(L_s)) times its mass.
Hence the two columns are h_sigma[-Ac N-sigma*u K+O(L_s^-1/2)],
Ac=|c0|sqrt(1+u^2)>0. One finite Lmin makes this error smaller than the
strict target-cone gap, uniformly over the frozen representatives. The
normalized two-column inverse is then bounded; its physical Haar inverse
retains r0^-3/2 S^1/2. The target is the ORIGINAL leading stress
(Q/q)^(A+1/2)T0, not an undocumented replacement of it by That.
The same positive roots have lower bound c r0^-3/4 S^1/4 sqrt(zeta),
and every fixed jet of the primary coefficient satisfies

\[
 |D^I a|<=C_I r0^{-3/4}epsilon^{1/2}
 S^{3|I|+1/4}\sqrt{zeta}\,delta^{-B_I}P.            \tag{10}
\]

Here zeta=exp(-t1^2/ya^2-4/yb^2), delta=min(1,ya,yb). Matrix inverse
and square-root differentiation retain sqrt(zeta), with polynomial edge
losses; they do NOT replace it by zeta. The exact flat constants in
sqrt(zeta) are t1^2/2 and 2. All fixed derivatives tend to zero at both
shell edges, so these roots extend smoothly by zero there.

Put cutoffs on vector potentials and pressures BEFORE taking the curl.
The primary pressure and potential coefficients have epsilon exponent 1.
The complete curl amplitude is a+remainder with remainder exponent
1-kappa>1/2, retaining all graph, cutoff and cylindrical connection terms.
There is no positive lower bound claimed for the oscillating/cutoff curl.

## 4. All physical derivatives, including the graph and frame terms

On physical evaluation the exact operators on coefficients are

\[
 partial_r=Q^{-1/2}(partial_R+M_i dr R^{dr-1}L_i),\quad
 partial_z=Q^{-D}partial_Z,\quad
 partial_t=Q^{-1-h}(-epsilon partial_T+c_iN_i).      \tag{11}
\]

All repeated derivatives also differentiate R^(dr-1). Common-torus
changes have fixed covering difference <=3 and bounded eigenvector
factors, not new epsilon losses. Since M_i<=epsilon^-kappa S^-rho_g,
spatial COEFFICIENT derivatives cost at most Q^(-1/2-h kappa), and time
derivatives Q^(-1-h), times fixed polynomial S factors. Cartesian frame
derivatives cost Q^-1/2 since R stays away from zero. These do not cancel.

For the PHASE alone, L_i Phi=0 because L_i v=0. Exactly grad_x v=0 and
partial_t v=Q^(-1-h). The explicit axial phase term costs
Q^-D epsilon^-1=Q^-1/2. Repeated chain rules give
|grad_x^a partial_t^b Phi|<=C_ab Q^(-a/2-b(1+h)) S^B for a+b>=1.
For a Faadi-Bruno partition into j phase factors, j<=a+b. Retaining
k<=2Q^-h/2 gives the UNIFORM harmonic losses

\[
 s_x=1/2+h/2,\quad s_t=1+3h/2,
\quad |grad_x^a partial_t^b e^{ikmPhi}|
 <=C_abm Q^{-a s_x-b s_t}S^B.                     \tag{12}
\]

These dominate coefficient/frame losses because kappa<1/2. Spatial curls
have the common physical Q^-1/2 factor, with normalized Dz=epsilon*dZ.
Velocity, pressure and fixed potential rescale by Q^-A,Q^-2A and
Q^(1/2-A), respectively. Thus every fixed physical mixed derivative has
the following bounds, uniformly over ALL accepted bands, labels and copies:

| Primary physical field | Bound for a spatial and b time derivatives |
|---|---|
| Complete curl velocity | C Q^(-(a+1)sx-b st) S^B sqrt(zeta) delta^-B P |
| Wave pressure | C Q^(-(a+2)sx-b st) S^B sqrt(zeta) delta^-B P |
| Fixed vector potential | C Q^(-a sx-b st) S^B sqrt(zeta) delta^-B P |

Constants include fixed radius, width, h,N and derivative/harmonic order.
They are not uniform as these vary. The graph cancellation applies to
phase derivatives, not to amplitudes or full curls; omitting it in the
phase overestimates losses, while omitting the graph in coefficients is
an incorrect identity. The background's radial velocity O(q^h) and complete
pressure, including the first two I4 pressure tails, have not been erased.

## 5. Support and Gaussian tails on the SAME domain for every order

Retain the replay's frozen cutoff: transverse |xi|<=r0/2 and
|v/L_s-1/2|<=2/5, with cutoff one for |v/L_s-1/2|<=1/5. No narrower
manuscript cutoff is silently substituted. The potentials and pressures
are smooth with support strictly inside the fixed auxiliary rectangles.
Curls and their derivatives cannot enlarge these supports. The exact
1250-color separation of CLOSED enlarged preimages gives at most one
nonzero physical label; all cross-label derivative products are zero.
Local finiteness holds at q>0. Axis smoothness comes from zero support
near r=0 at each such point, not smoothness of r^dr at the axis. The
T=0 boundary has one-sided jets, not temporal zero extension.

P is positive at finite pulse endpoints and cannot supply support by
itself. For a sourced global residual the ORIGINAL full harmonic source
must also have smooth endpoint zero extension. This hypothesis remains.
On supports of psi' and (1-psi)f, (7) gives P<=exp(-a0 ell^2),
a0=cG*r0/25>0. After any fixed physical derivative and source order,
the remaining upper bound is C ell^B Q^-Lambda exp(-a0 ell^2).
For ANY fixed desired K>=0 and ell>=ell0>=1,

\[
 ell^B Q^{-(Lambda+K)}e^{-a0 ell^2}
 <=\exp(D ell-a0 ell^2)<=\exp(D^2/(4a0)),\quad
 D=B+max(0,Lambda+K),                              \tag{13}
\]

using log ell<=ell and log2<1. Thus every such tail is O_(I,K)(Q^K)
on the ENTIRE SAME band domain, not just after a derivative-dependent
starting band. Polynomial edge losses are absorbed by sqrt(zeta) with
the exact positive constants t1^2/2 and 2; the constants may be huge.
Equation (13) gives a direct all-half-line constant, rather than scanning
bands or retuning h,N,r0,ell0 when a new derivative is requested.

## 6. Checks and next obligation

Exact rational multivariate Taylor jets compare direct physical pullbacks
with repeated graph/Cartesian operators for three labelled wave fixtures,
including complete curls, pressure and frame derivatives through total
physical order three. Separate controls retain the radial graph,
its differentiated coefficient and the fast auxiliary-time term; phase
radial-graph cancellation is exact. Rational moving-frame and normal-pressure
identities, implicit band-coordinate jets, all physical exponent partitions,
Riccati barriers and common-domain Gaussian square completion are checked
independently. Fixtures are not the actual gigantic background or a numerical
PDE solution. The analytical statements are (1)--(13), not a continuum
supremum certificate obtained from tests.

Next: independently audit COMPLETE-wave covariance means with radial and
axial velocities and pressure, then the actual generated source classes
and full residual including viscosity, all cutoffs and interactions. The
mapping lemma for a prescribed supported source is not that final source
verification. Earlier conditional hierarchy/realization results remain
conditional; no unconditional blow-up or complete residual improvement.

```sh
python3 research/navier_stokes_cascade/results/physical_pulse_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/physical_pulse_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/physical_pulse_foundation/evidence.json \
  --output /tmp/physical_pulse_foundation_fresh.json
```
