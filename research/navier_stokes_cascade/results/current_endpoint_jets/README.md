# Current complete field: closed-time atlas and punctured terminal jets

Status: `conditional-current-field-closed-time-atlas-and-punctured-terminal-jets`.

The [current realized field](../current_potential_realization/README.md),
not the earlier hypothetical hierarchy, has compatible smooth one-sided
potential, direct-swirl, velocity and pressure jets at t=1 AWAY FROM x=0,
conditionally on its retained uniform all-order physical estimates.
This checkpoint audits the closed-time geometry and proves the transfer
of those estimates. It also derives every exterior heat/pressure time
jet and proves their ordinary terminal Taylor series have ZERO radius
of convergence when the exterior constant is nonzero. Smooth jet gluing,
not analytic continuation or pasting zero, is needed in the later step.

This is not a current global-force, global-energy or blow-up theorem.
The leading profile, background and all-order inverse estimates remain
analytical premises requiring independent review. Manufactured diagnostics
are not actual terminal fields or continuum certificates. All historical
sources and frozen evidence remain unchanged.

## 1. Assemble the complete representatives, including the seed

Keep the SAME h, N, positive lambda, width C^-20, common qstar, frozen
geometry and primary roots. Put A=1/2+h, D=1/2-h, tau=1-t and

\[
q=\tau+z^2q^{2h},\quad X=r^2/(2q),\quad\eta=z/q^D,
\qquad 0<h<1/100.\tag{1}
\]

The positive branch has |eta|<=1 for tau>=0, except that q=0 at
z=tau=0 is excluded. Its implicit derivative denominator is
`L=1-2h eta^2>=1-2h>49/50`. In particular at z!=0, tau=0,
`q=|z|^(1/D)`, `q_t=-1/(1-2h)`; no inverse power of `1-eta^2` is needed.

The current finite seed has the exact representative
`u^[0]=curl A_seed+B_seed`. A_seed is the summed background azimuthal
potential plus ALL finite seed wave potentials and mean azimuthal
potentials; B_seed is the base direct swirl plus all finite mean swirls.
These are the actual representatives already used to construct the seed,
not potentials reconstructed from a velocity estimate. With the current
diagonal scales a_j and complete increments from the parent, set

\[
\mathcal A=\mathcal A_{seed}+\sum_{j\ge1}\chi(a_jq)\mathcal A_j,
\quad\mathcal B=\mathcal B_{seed}+\sum_{j\ge1}\chi(a_jq)b_j e_\theta,
\]
\[
p=p^{[0]}+\sum_{j\ge1}\chi(a_jq)(p^{[j]}-p^{[j-1]}),
\qquad u=\operatorname{curl}\mathcal A+\mathcal B.\tag{2}
\]

All curl cutoff terms are therefore present. B is axisymmetric; no full
compact vector potential for its direct swirl is assumed. The SAME p_B
and both exceptional base-pressure tails stay in p^[0], in its fixed
gauge. Every changing wave/mean pressure and complete axial memory is in
the pressure differences. No pressure gauge is reset on the terminal slice.
Actual flat velocity pieces and exact Gaussian cutoff tails stay in the
field and residual; no endpoint forward inversion is performed.

The inherited background potential/direct swirl have the stated Cartesian
regular factors at the axis. All annular corrections vanish for X<Xa>0.
In the heat exterior X>=Xext, the potentials vanish and only the unchanged
base direct swirl K and its full pressure remain. These support and
background statements are the explicit inherited inputs, not consequences
of numerical tests in this checkpoint.

## 2. Closed-time geometry is finite where q is bounded below

On `c<=q<=qmax<q_e/2` (the observation domain q_e is defined below), a dyadic label with Q=2^-ell and
`1/2<=q/Q<=2` must satisfy

\[
\left\lceil\log_2\frac1{2q_{max}}\right\rceil\le\ell
\le\left\lfloor\log_2\frac2c\right\rfloor,\quad\ell\ge\ell_0.\tag{3}
\]

This is a finite CLOSED band set, including the boundaries. The admitted
ell0 may be arbitrarily large but is fixed before the construction.
For each fixed band the slow active box has bounded R,Z,T, with
`T=tau/Q>=0`, including T=0. Thus `s-Z^2*s^(2h)=T` and normalized
physical time differentiation has slow coefficient `-Q^h=-epsilon`,
not -1. For a closed interval [lo,hi] and
mesh delta=ell^-6, all closed support boxes intersecting it have indices

\[
\left\lceil lo/\delta\right\rceil-1\le a
\le\left\lfloor hi/\delta\right\rfloor+1.\tag{4}
\]

In particular T=0 uses the already selected indices a_T=-1,0,1.
There is no new label selection, carrier rounding or color choice at t=1.
The frozen atlas explicitly selected labels on this closed active set.
Its squared partitions, 1,250 colors, closed support separation, at most
162 pointwise slow labels and 2,744 common lifts apply unchanged.
The exact geometry certificates are replayed; the small diagnostic band
window is NOT claimed to be the actual pulse-acceptance threshold.

The physical auxiliary map `Y=v_r r^d_r+v_t t mod Z^2` is smooth for
r>0. All common-chart formulas are pullbacks of the same absolute Y
field. Liftwise ODE uniqueness, physical radial primitives and the SAME
q-based mean/pressure bumps preserve compatibility. Independently chosen
fixed-R bumps would not. At support changes smooth zero extensions identify
all jets. The fractional map is NOT asserted smooth at r=0: for q>=c,
all annular auxiliary-dependent fields vanish uniformly on
`r<sqrt(2Xa c)`, and the inherited Cartesian background handles the axis.

Also `a_j>2^j/qstar` since `a_1>2/qstar` and scales double. An active
cutoff requires `a_j c<1`, hence `j<log2(qstar/c)`. Thus only finitely
many correction stages survive on this whole region. The background
summation is likewise locally finite. Finite-stage physical bounds from
the current parent and the fixed seed are uniform in physical point,
label and all fixed jets. On this compact q interval, their powers of
q and log(q) are bounded. Finite band/mesh/stage sets, or the inherited
bounded overlap, therefore give

\[
\sup_{E\times[1-\epsilon,1)}
|D_x^\alpha\partial_t^j(\mathcal A,\mathcal B,p)|<\infty
\quad\hbox{for every fixed }\alpha,j,\tag{5}
\]

where E lies in the q-positive part of a local punctured terminal compact.
This derives the endpoint bound from the retained physical estimates;
it does not independently prove the analytic profile/inverse estimates.
If those estimates fail to be uniform up to |eta|=1, (5) is the exact
blocking inequality. Interior smoothness or origin-flat residuals alone
cannot replace it.

## 3. Cover the punctured terminal set, including the equatorial plane

Restrict only the observation cylinder, not the actual field, parameters
or accepted labels. Set `q_e=min(qstar/16,2^(-ell0-2))`, then work on
`r<R0, |z|<z0`, with `tau0=q_e/16`, `z0=(q_e/16)^D`.
From (1), either q<=2tau or
`q<=(2z^2)^(1/(1-2h))`, so throughout this cylinder at late times

\[
q\le2\tau+4|z|^{1/D}\le3q_e/8<q_e/2.\tag{6}
\]

Every active Q satisfies Q<=2q, hence `4Q<=3q_e<=3qstar/16<qstar`.
Also `Q<=3q_e/4<=3*2^-ell0/16`, so all active labels satisfy the
fixed admitted first-band requirement. This explicit margin avoids a
domain ambiguity in a naive endpoint transfer: `q<qstar/2` ALONE does
not imply `4Q<qstar` (take q=qstar/4 and Q=qstar/2). A qstar already
defined as a working physical domain might contain its own margin;
we do not rely on that interpretation or alter historical frozen notes.

Let E be a compact terminal set in this cylinder with |x|>=d>0.
Use a slightly larger punctured compact when identifying spatial jets.
Choose

\[
b=\min(q_e/4,d^2/(16Xext),1/4),\quad
\zeta=\min(d/4,(b/16)^D),\quad \tau\le b/16.
\]

On |z|>=zeta/2, `q>=|z|^(1/D)>=(zeta/2)^(1/D)>0`, so (5) applies.
On |z|<zeta, |x|>=d implies r>d/2 and (6) gives q<=3b/8.
Consequently `X>=d^2/(3b)>=16Xext/3>Xext`: a whole neighborhood,
not just the plane z=0, is in the unchanged heat exterior. These two
regions overlap and cover E. The same argument covers positive-radius
equatorial points and axis points with z!=0; it does NOT cover x=0.

## 4. Exterior heat and COMPLETE pressure terminal jets

In that exterior, keep the fixed pressure gauge

\[
K=C r^{-1-2h}H(4\tau/r^2),\quad
H(Z)=\frac1{\Gamma(1+h)}\int_0^\infty e^{-v}v^h(1+Zv)^{-h}\,dv,
\quad p=-\int_r^\infty K(\rho,\tau)^2\frac{d\rho}{\rho}.\tag{7}
\]

Positive integration gives `|H^(j)(Z)|<=(h)_j(1+h)_j` for ALL Z>=0.
Put `A_j=(h)_j(1+h)_j` and
`S_j=sum_(k=0)^j binom(j,k) A_k A_(j-k)`. Differentiation under the
integrals and tau=1-t give the exact terminal jets

\[
\partial_t^j K(r,1^-)=C4^j A_j r^{-1-2h-2j},
\]
\[
\partial_t^j p(r,1^-)=-\frac{C^2 4^j S_j}{2+4h+2j}
r^{-2-4h-2j}.\tag{8}
\]

The pressure uses the FULL Leibniz convolution, not just 2K times its
jth derivative. Its density has integrable majorant
`C_j rho^(-3-4h-2j)` for rho>=rmin>0. Dominated convergence gives uniform
time jets; radial derivatives follow by the lower integration endpoint
and differentiation of K. Angular/Cartesian frame derivatives are bounded
on positive-radius compacts and z derivatives vanish. This supplies (5)
on the second region too. The exact identity `p_r=K^2/r` is retained at
every time-jet order. Reversing the tau/time sign or dropping pressure
cross terms fails the explicit checks.

There is a useful precise limitation. If C!=0, the ratio of successive
time Taylor coefficients of K, in absolute value, is

\[
\frac4{r^2}\frac{(h+j)(1+h+j)}{j+1}
=\frac4{r^2}\left[j+\frac{h(2j+1+h)}{j+1}\right]\longrightarrow\infty.
\tag{9}
\]

For j>=1, `A_j>=h (j-1)! j!` and `S_j>=2 A_j`. Equation (8) therefore
also forces zero Taylor radius for p. This holds for every fixed h>0,
including the unchanged very small research h; no small-h truncation is
permitted before j->infinity. It is NOT an obstruction to C-infinity
one-sided jets or smooth Borel-type gluing. It rules out using a convergent
terminal Taylor series as that gluing. The integral for H at Z>=0 does
not itself prescribe a real smooth future formula at Z<0.

## 5. Compatible terminal jets, with the correct quantifiers

For each G in `(A,B,p)` and each fixed alpha,j, the NEXT time-derivative
bound in (5) and the fundamental theorem of calculus give

\[
\|D_x^\alpha\partial_t^jG(t_2)-D_x^\alpha\partial_t^jG(t_1)\|_E
\le M_{\alpha,j+1}|t_2-t_1|.\tag{10}
\]

Every jet thus has a uniform terminal limit. Pass to uniform limits in
spatial coordinate integrals on enlarged compact boxes to identify all
spatial derivatives of F_j=lim partial_t^j G. Time integration identifies
consecutive time jets. For every fixed finite k this also yields

\[
\left\|D_x^\alpha\left[G(1-\tau)-
\sum_{j=0}^k\frac{(-\tau)^j}{j!}F_j\right]\right\|_E
\le\frac{M_{\alpha,k+1}}{(k+1)!}\tau^{k+1}.\tag{11}
\]

Curl and its derivatives commute with these limits: a velocity m-jet
needs potential derivatives through m+2 and direct-swirl derivatives
through m+1 for the Cauchy bound. Pressure needs m+1. Every fixed finite
requirement is covered, with NO uniform bound over derivative order,
no interchange with an uncut infinite series and no new endpoint atlas.
Neither (10) nor (11) asserts convergence of the infinite Taylor series.
Nonzero terminal pressure jets cannot be zero-extended in time. These
are one-sided compatibility statements, not a future velocity construction.

Next: join axisymmetric spatial/time localization and zero initialization
to this SAME current field, retaining every curl, pressure, viscous and
nonlinear term; then glue the resulting force smoothly, establish global
energy/dissipation and audit the comparison argument. The origin remains
singular conditionally on the inherited lower bound; no regular velocity
or pressure there, global admissible force, or unconditional theorem is
inferred from the present punctured endpoint result.

## 6. Reproduce

```sh
python3 research/navier_stokes_cascade/results/current_endpoint_jets/test_audit.py -v
python3 research/navier_stokes_cascade/results/current_endpoint_jets/audit.py \
  --verify-record research/navier_stokes_cascade/results/current_endpoint_jets/evidence.json \
  --output /tmp/current_endpoint_jets_recheck.json
```

Exact rational certificates cover closed band/mesh endpoints, all-order
jet requirements and heat/pressure identities. Manufactured third-order
Cartesian jets use the ACTUAL implicit q at both terminal poles and two
charts; wrong fixed-R bumps, time sign and zero extensions are controls.
Independent positive Gamma integrals check the heat moments and complete
pressure convolution. The current realization and historical endpoint
record are replayed, without importing the latter's force/energy conclusions
for the current field. Raw roundoff scores are tested but not re-frozen
under a different ancestor comparator. Passing checks are not a continuum proof.
