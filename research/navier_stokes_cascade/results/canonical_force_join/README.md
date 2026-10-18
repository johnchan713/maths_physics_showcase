# Same-canonical-field compact force, protected-disk flatness and energy

Parent: `99b5b2a72c16324a0c9389defa5e1d7a9cdfe76f`.
Status: `conditional-canonical-compact-force-protected-disk-flatness-and-preterminal-energy`.
On 18 October 2026 the parent's canonical endpoint CI step passed; the
complete hosted workflow was still running without a failed step.
Historical sources, fields and frozen evidence are unchanged.

We localize the SAME field from the [canonical endpoint join](../canonical_endpoint_join/README.md),
construct one compact smooth global FORCE with zero initial velocity, and
derive preterminal energy/dissipation bounds, relative to retained exact-core,
background and finite analytical defining/existence/inverse/gap premises.
The reused [force mechanism](../force_foundation/README.md) is applied to
canonical representatives and targets, not its historical field or force.
Neither future velocity/pressure nor global competitor exclusion is proved.
There is no unconditional admissible-force, residual or blow-up claim.

A new consequence of the uniform corner estimate is retained explicitly:
ALL terminal force jets are flat along the entire protected equatorial
disk, not just at the origin. The support-preserving future force can keep
that flat disk for every future time. Exact rational upper bounds choose
its one dyadic Borel schedule without floating-point acceptance tests.
Actual infinite coefficients, force norms and extension scales are NOT
numerically evaluated. Finite controls are not continuum proofs.

Sources: retained manuscript Section 10, (10.1)--(10.9), PDF SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`,
and the analytical canonical/foundation checkpoints. The manuscript's
final theorem is never an input.

## 1. Canonical localization with a domain margin

Write the actual canonical preterminal field as
u0=curl Ac+Bc, Bc=bc(r,z,t)e_theta, partial_theta bc=0, p0=pc complete.
Ac includes the actual background, seed, wave, fast/slow mean and diagonal
potentials, with its fixed operator gauge. Direct swirl stays separate.
The full pc keeps both Pi1/Pi2 base tails, all six-entry covariance pressure,
axial velocity and moving pressure moment, and every changing intermediate
pressure. No unforced Riesz pressure or new gauge is substituted.

Keep Md64, viscosity one, original h,N, fixed positive lambda, radius,
width C_setup^-20, all labels/carriers/supports, qstar and both existing
diagonal schedules. For Aexp=1/2+h,D=1/2-h,tau=1-t, use the already
proved observation constants q_e=min(qstar/16,2^(-ell0-2)),
tau0=q_e/16,z0=(q_e/16)^D. The enlarged 1/4<=q/Q<=4 charts obey
q<=3q_e/8, 4Q<=16q<=6q_e<=3qstar/8 and Q<=3*2^-ell0/8.

Choose a Cartesian-smooth axisymmetric c=chi_r(r^2)chi_z(z)chi_t(t),
supported spatially in K inside that observation cylinder, and zero for
tau>=tau0. Take c=1 for r<=rL/2, |z|<=z0/2, tau<=tau0/2, where rL is a
new localization radius inside the observed domain, NOT a retuned pulse
radius. Define fields on a NEIGHBORHOOD of the support closure before
extending by zero; zero times an undefined field is not a definition. Set

\[
u=\operatorname{curl}(cA_c)+cB_c=cu_0+d,\quad
d=\nabla c\times A_c,\qquad p=cp_0.                              \tag{1}
\]

Then div u=0, because the first summand is a curl and the direct swirl
remains axisymmetric. Also d dot grad c=0 exactly. These are genuine smooth
compact fields for every t<1, zero on an initial time neighborhood. Neither
naive velocity cutting nor theta-dependent cutting of direct swirl is used.

## 2. The COMPLETE localized residual and pressure equation

For R_nu(v,p)=v_t-nu Delta v+(v.grad)v+grad p, exact expansion gives

\[
\begin{aligned}
R_\nu(u,p)={}&cR_\nu(u_0,p_0)+c_tu_0+d_t
-\nu(\Delta c)u_0-2\nu\sum_i c_i\partial_i u_0-\nu\Delta d\\
&+p_0\nabla c+c(c-1)(u_0\cdot\nabla)u_0
+c(u_0\cdot\nabla c)u_0\\
&+c(u_0\cdot\nabla)d+c(d\cdot\nabla)u_0+(d\cdot\nabla)d.
\end{aligned}                                                    \tag{2}
\]

All TWELVE terms remain at actual viscosity one. The apparent thirteenth
(d dot grad c)u0 is exactly zero, not estimated away. Every part of p0 grad c
remains, including exceptional pressure memory. A force m-jet requires
Ac through m+3, Bc through m+2 and pc through m+1; uniform Cauchy estimates
reserve one more time derivative. The canonical endpoint proof supplies
these fixed finite orders on punctured compacts, not one uniform constant
for all orders. No actual flat velocities or Gaussian ledger terms are dropped.

Define f_minus=R_1(u,p). Its exact divergence equation is

\[
-\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f_{minus}.       \tag{3}
\]

Force divergence need not vanish. A compact-support claim is not justified
by projecting the force with a nonlocal Leray map. The official
[Fefferman formulation](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf),
equations (4),(5),(7) and alternative (C), was checked on 18 October 2026:
it requires smooth rapidly decaying data/force and bounded-energy global
competitors, not divergence-free force. Meeting its force conditions alone
does not establish the no-global-solution conclusion.

## 3. Uniform force jets and a flat protected equatorial DISK

On a fixed plateau box r<rL/2, |z|<z0/2, tau<tau0/2, every positive
derivative of c vanishes, so f_minus is the COMPLETE canonical residual.
Its all-fixed-jet O(q^N) bound is uniform for bounded X; for X beyond the
whole support hull the exact heat residual is zero. These regions exhaust
the box, INCLUDING paths with r changing with tau. Since
q<=2tau+4|z|^(1/D), for every fixed alpha,j,N,

\[
|D_x^\alpha\partial_t^j f_{minus}(x,1-\tau)|
 \le C_{\alpha jN}(\tau+|z|^{1/D})^N.                            \tag{4}
\]

This is uniform on the plateau box, not merely along the protected core
path. Outside a small corner ball, the canonical punctured endpoint bounds
control every term of (2), including its next time jet. A finite cover of
the rest of K, together with exact zero outside K, makes every force jet
uniformly Cauchy on R3 as t increases to 1. Let

\[
F_j=\lim_{t\uparrow1}\partial_t^j f_{minus}\in C_c^\infty(K).
                                                                    \tag{5}
\]

Uniform limits in spatial coordinate integrals identify all derivatives,
including across the origin. Time integration identifies successive time
jets. The force is regularized, not terminal velocity or pressure.

Taking tau to zero in (4) gives |D_x^alpha F_j|<=C|z|^(N/D), uniformly
on smaller closed plateau boxes. Hence for ALL alpha,j,

\[
D_x^\alpha F_j(x_1,x_2,0)=0\qquad(r<rL/2).                        \tag{6}
\]

Thus the traces are infinitely flat along the protected equatorial disk.
No new force-free OPEN neighborhood is asserted. The pressure/heat endpoint
jets away from the plateau can be nonzero. The Gaussian residual ledger
is used only through the one fixed-stage comparison defining the same
canonical realization; it is not independently resummed or forward-inverted.

## 4. One compact future FORCE; exact dyadic scale acceptance

For sigma=t-1>=0 choose smooth chi0=1 on [0,1/2], zero on [1,infinity),
and define

\[
f_+(x,1+\sigma)=\sum_{j\ge0}T_j,\quad
T_j=\chi_0(b_j\sigma)\frac{\sigma^j}{j!}F_j(x),\quad b_0=2.         \tag{7}
\]

The exact differentiated term is

\[
\partial_\sigma^m T_j=F_j
 \sum_{\substack{0\le k\le m\\m-k\le j}}
 {m\choose k}b_j^k\chi_0^{(k)}(b_j\sigma)
 \frac{\sigma^{j-m+k}}{(j-m+k)!}.                                \tag{8}
\]

On the WHOLE support 0<=sigma<=1/bj, each summand has bound proportional
to b_j^(m-j), not only where chi0 has positive derivatives. For each j>=1,
collect the finitely many spatial multi-indices alpha and m with
|alpha|+m<=floor(j/2); equivalently pad C_jam>=1 uniformly for their
finite order classes a=|alpha|. Choose rational c_jam>=max(0,log2 C_jam).
With integer n0=1 and b_j=2^n_j, any integer n_j strictly above

\[
\max\{n_{j-1}+1,\ 2,\ \max_{a,m}(c_{jam}+j+2)/(j-m)\}            \tag{9}
\]

satisfies C_jam*b_j^(m-j)<=2^(-j-2). Denominators are positive since
m<=floor(j/2)<j. The code takes exact rational ceilings and adds one;
it does not materialize an enormous 2^n_j or evaluate actual terminal norms.
Finite rational upper bounds exist for any finite real norms. No derivative-
growth rate, Gevrey bound, convergence of a raw Taylor series or uniform
constant across j is assumed. The parent's sharp EXTERIOR Gevrey-2 result
does not give such an assumption for the complete force, and none is needed.

For fixed alpha,m, all j>=max(1,2(|alpha|+m)) are controlled. Their mixed
derivative tail is at most the coarser 2^(1-J) for a tail starting at J.
Every earlier term is smooth. Thus ONE increasing sequence yields uniform
convergence of every fixed mixed differentiated series, including sigma=0.
The plateau and factorial make its mth time jet exactly F_m. Gluing f_plus
to f_minus gives

\[
f\in C_c^\infty(\mathbb R^3\times(0,\infty)),\quad
\operatorname{supp}f\subset K\times[1-\tau0,3/2].                 \tag{10}
\]

All polynomially weighted derivative norms are finite by compact smooth
support. Zero initial velocity has the required data bounds. This is a
conditional existence construction of one force, NOT future velocity or
pressure, and not analytic or zero continuation of nonzero heat-pressure jets.

Since each F_j has (6), every spatial/time derivative of each term in (7)
is zero on that disk for every sigma. Uniform differentiated convergence
preserves this property for the full force. A finite normal Taylor formula
at z=0 therefore gives, for every fixed alpha,m,N,

\[
|D_x^\alpha\partial_\sigma^m f_+(x,1+\sigma)|
 \le C_{\alpha mN}|z|^N,\quad
r\le rL/4,\ |z|\le z0/4,\ 0\le\sigma\le1/2.                     \tag{11}
\]

The constant is a finite uniform higher normal-derivative norm on that
compact box. This proves quantitative future preservation of the flat disk
without a new weighted cutoff sequence or uniform-in-order bound.

## 5. Preterminal energy without circularity

Only AFTER constructing this compact smooth force, let
E(t)=||u(t)||_2^2, D(t)=||grad u(t)||_2^2, F(t)=int_0^t ||f(s)||_2 ds.
For every t<1 the actual compact support allows integration by parts:

\[
\tfrac12 E'(t)+D(t)=\langle f(t),u(t)\rangle.                     \tag{12}
\]

Pressure and transport fluxes vanish; force divergence need not. Starting
at E(0)=0, regularize sqrt(E+delta^2), whose derivative is at most ||f||2.
Let delta decrease to zero to obtain sqrt E<=F. Integrating (12) then yields

\[
E(t)+2\int_0^tD(s)ds\le F(t)^2,
\quad\sup_{t<1}E(t)\le F(1)^2,\quad
\int_0^1D(s)ds\le F(1)^2/2.                                    \tag{13}
\]

The last two bounds are SEPARATE: do not add their two extrema while
keeping the same right side. The historical scalar rise/decay counterexample
is retained and replayed. Nonnegative dissipation uses monotone convergence;
no global smooth velocity or a priori finite energy is assumed to prove (13).

For tau<=min(tau0/2,rL^2/(8Xin)), the protected core path
r=sqrt(2Xin*tau),z=0 is inside the unit plateau. The inherited conditional
u_theta>=(e0/2)*tau^-Aexp remains unchanged. This does not prove the retained
positive core lemma, nor exclude a different global solution.

## 6. Meaning, checks and next obligation

The same-canonical-field localization/compact-force/energy implication is
closed relative to retained analytical premises. The force is not assumed
solenoidal; no regularity of terminal velocity at the origin, future velocity,
global uniqueness or unconditional blow-up conclusion is inferred.

Next recover pressure and justify the expanding-ball comparison for EVERY
admissible smooth finite-energy global competitor with THIS same zero datum
and THIS compact force, without unproved competitor pressure, derivative
or boundary-flux hypotheses. Old global fields are not automatically the
canonical construction. Independent discharge of retained exact-core,
background and finite defining/existence/inverse/gap premises also remains.

Reused exact polynomial checks expand all twelve localization terms and
the full forced pressure equation, keep all Borel Leibniz/factorial terms,
and integrate the compact energy identities. New exact integer/rational
scale checks use extreme SYNTHETIC log2 norms. The new plane fixture has
288 zero mixed component controls but a NONZERO eighth normal jet: it is
only a finite-order polynomial control, not an infinite-flatness or compact
C-infinity fluid proof. All-order assertions follow from (4)--(11) and
their retained analytical hypotheses, not from passing finite tests.

```sh
python3 research/navier_stokes_cascade/results/canonical_force_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/canonical_force_join/audit.py \
 --verify-record research/navier_stokes_cascade/results/canonical_force_join/evidence.json \
 --output /tmp/canonical_force_join_fresh.json
```
