# First sourced wave: regular-source gain and a flat-tail obstruction

Status: `first-sourced-wave-core-gain-with-exact-retained-cutoff-tails`.
Parent: `6eeb414542ecf473826a6bc9174fc050fd73160b`.

This checkpoint joins the CURRENT nonzero-angular residual to one forward
pulse solve, cuts its potential AND pressure, and recomputes the complete
wave covariance and mean pressure. With the retained continuum lemmas the
regular residual improves from order B=1-3kappa to D=3/2-6kappa. Designated
Gaussian cutoff tails remain exact additive residuals, not forward sources.
An exact homogeneous-pulse example proves why a physically flat source
cannot simply be inverted as though its response were flat.

This is ONE finite step using the inherited background, leading-profile,
common-torus geometry and inverse estimates. It is not an independent
foundation certificate, an actual giant-band numerical solve, an all-cycle
construction, an unconditional complete residual estimate, or blow-up.
The frozen historical evidence and GOAL.md are unchanged.

The proposed manuscript under audit is the [supplied source](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(6.19)--(6.32), (7.39)--(7.42), (8.12)--(8.27), (9.1)--(9.10).
PDF SHA-256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
We use the rederived identities in the frozen checkpoints, not the source's
ultimate conclusion. In particular the [conditional cycle](../residual_cycle_closure/README.md)
already excludes Gaussian tails from forward sources; the obstruction here
explains and enforces that restriction, rather than revoking that record.

## 1. Current state and precise source split

Keep the SAME Md64 data, N, h in (0,1/100), positive lambda, width C^-20,
base, supports, labels and qstar. Write epsilon=Q^h, S=ell^2,
k=ceil(epsilon^-1/2), kappa=1/100000. Use W_alpha, M_alpha and S_alpha
with ALL fixed coefficient derivatives, Gaussian P, sqrt(zeta) wave
weight, zeta mean weight and fixed-order inverse edge-distance losses.
Constants/polynomial S losses are uniform in band, copy and label; they
may depend on this finite stage, finite harmonic set and derivative order.

The [parent](../post_signed_mean_relaxation/README.md) has the original
complete real primary plus one complete phase-locked real signed wave,
w in W_1/2 and w-W0^as in W_(1-kappa). Its current tangential means have
order mu=1-kappa, radial mean order mu+1, full tangential mean residual
M_(3-6kappa), and defects S_(3-3kappa), plus separately retained flat errors.
The old slow Haar means are NOT zero. The full background p_B (including
both exceptional tails), wave pressures and current mean p_m are retained.

Define the EXACT current nonzero-angular residual as

    G_total = G_core + T_old.

Here T_old is the sum of the full homogeneous potential/pressure cutoff
commutators defined in section 2, with their actual signs in the residual.
For the remaining regular expression, the [complete-curl interaction](../residual_cycle_closure/README.md)
and [actual-band estimates](../actual_band_pulse_bounds/README.md) give:

| Regular source contribution | Wave order |
|---|---:|
| Cut original full linear remainder | 1-3kappa |
| Cut first real signed full linear remainder | 2-6kappa |
| Complete old wave-wave terms | 1-kappa |
| Current actual mean-wave terms | mu+1/2-1/2 = 1-kappa |

The parent-generated mean-wave terms of orders 2-4kappa and 2-3kappa are
included, not substituted for the old residual. Thus G_core is in W_B,
B=1-3kappa. The COMPLETE cut curl velocities, including flat cutoff curl
pieces, remain in every nonlinear product. We do not expand a product and
discard its flat factors. Pure base/mean/pressure fields have angular
harmonic zero, but the actual mean coefficients depend on FULL pressure.

Each nonzero term contains a wave, hence retains its P envelope and slow,
transverse and temporal closed support, even though a torus mean can fill
the torus. Same-label quadratic harmonics add; distinct-label derivative
products vanish by the inherited separated CLOSED enlarged supports.
The harmonic set is {-2,-1,1,2}; every coefficient extends smoothly by zero
at both pulse endpoints. Product weights zeta P^2 are stronger than the
required sqrt(zeta) P. Exact incompressibility of each complete curl
amplitude gives n.a=i div_*a/(km), so same-label quadratic transport loses
only kappa, not 1/2. Mean transport really can lose 1/2. These are the
specific source hypotheses, not a numerical inference from tests.

## 2. Full cutoff commutator and the forward-inverse obstruction

Let L_b u = partial_t u+(b.grad)u+(u.grad)b-Delta u, viscosity one,
u=curl A, u_psi=curl(psi A)=psi u+d_psi,
d_psi=grad psi cross A. Cut pressure as psi pi. The exact commutator is

\[
\begin{split}
T_\psi(A,\pi)&=L_bu_\psi+\nabla(\psi\pi)
                  -\psi(L_bu+\nabla\pi)\\
&=(\psi_t+b\cdot\nabla\psi)u
 -2(\nabla u)\nabla\psi-(\Delta\psi)u
 +\pi\nabla\psi+L_bd_\psi .
\end{split} \tag{1}
\]

This Cartesian identity automatically includes cylindrical frame terms
after coordinate transfer. It keeps viscosity, pressure, transport AND the
cutoff curl. The actual relative temporal cutoff psi is one centrally and
zero on endpoint collars. Its derivatives are confined away from the
Gaussian maximum; P<=exp(-c S) there. All fixed physical differentiation
costs are powers of Q^-1 and S, so exp(-c ell^2) beats them: (1) is physically
flat. Nevertheless a primary principal tail has epsilon order 1/2, not B,
if the SAME P is retained: its P cancels when dividing by the envelope.
Trading that P for an epsilon gain loses the input required by the inverse.

There is a sharper obstruction than a missing estimate. For ANY actual
nonzero homogeneous constrained principal pulse (t_h,pi_h), choose a
genuine C-infinity slow step chi, zero initially, one before the core, with
chi' supported in [L/10,L/5]. In the principal equation with right side -f,

\[
f=\chi' t_h,\qquad (t,\pi)=(-\chi t_h,-\chi\pi_h) \tag{2}
\]

is its EXACT zero-entry solution. On that early interval the actual Gaussian
is exp(-cL), L comparable to S, so f is flat in every fixed physical jet.
It is zero on both endpoint collars. At the core chi=1, however, and the
actual homogeneous lower bound |t_h|>=cP makes the response nonflat.
Thus smooth endpoint zero extension and physical flatness of a source do
NOT imply a flat response or zero exit. Normal constraint and pressure in
(2) are preserved exactly, including moving n. No fixture is needed for
this analytical counterexample.

For a transparent finite diagnostic take P(v)=exp(-(v-L/2)^2/(2L)).
On [L/10,L/5], P<=exp(-9L/200). A quintic step has derivative
30x^2(1-x)^2<=15/8, integral one; scaling its width to L/10 gives
|P chi'|<=75/(4L) exp(-9L/200), yet the scalar inverse core value is -1.
The zero/one quintic extension is C2, NOT C-infinity; it illustrates (2),
and does not replace the genuine smooth step used in the proof.
T_old must therefore stay additive: it is neither set to zero nor fed to
this forward inverse. The flat base residual is also retained separately;
it need not have a wave envelope or annular support.

## 3. First sourced solve, pressure and exact endpoint repair

For each current nonzero multiplier m take f=(G_core)_m along the fixed
common path. The constrained equation and its projected form are

\[
\begin{split}
t'+Kt+m^2dt+ikm n\pi&=-f,\qquad n\cdot t=0,\\
t'&=A_n t-m^2dt-\operatorname{Proj}_{n^\perp}f,\\
\pi&=\frac{i(n\cdot Kt-n'\cdot t+n\cdot f)}{km|n|^2}.
\end{split} \tag{3}
\]

The source need NOT be normal-orthogonal. Both n'.t and n.f in pressure
are essential; i/m also matters for negative harmonics and reality.
Use zero entrance data. The sourced propagator estimate from the
actual-band checkpoint maps W_B to W_B; fixed coefficient derivatives
follow its differentiated ODE recurrence. The audited geometry and
background supply those inputs on the SAME qstar, with no derivative-
dependent shrinking domain. Along a path with zero source, uniqueness
gives zero response, so slow/transverse supports persist. The pre-cutoff
response generally does not vanish at exit, where P is small but nonzero.

The uncut source potential A and pressure pi have order C=B+1/2.
Set delta w=curl(psi A), delta pi=psi pi BEFORE global zero extension.
Then delta w is a real divergence-free W_B field, its curl remainder has
order C-kappa, and every fixed jet extends smoothly at both endpoints.
Conjugate +/-m solutions ensure reality; no phase-locked real shortcut is
assumed. The uncut full linear remainder R_lin has order
D=B+1/2-3kappa. From (1) and (3), exactly

\[
G_{core}+L_b\delta w+\nabla\delta\pi
 =\psi R_{lin}+(1-\psi)G_{core}+T_\psi(A,\pi). \tag{4}
\]

The last two terms define T_new. They are flat, separately retained, and
are not forwarded to later inverses. This is a definition by exact
residual identity, not an estimate that erases errors. With full old/new
velocities in all nonlinear terms the updated REGULAR wave orders are

| Contribution | Order |
|---|---:|
| Full linear remainder | D=3/2-6kappa |
| Complete old/new wave interaction | B+1/2-kappa |
| New wave square | 2B-kappa |
| Actual total mean/new wave interaction | B+mu-1/2 |

Their minimum is D. Angular-zero parts go to the mean equations below;
nonzero parts retain wave supports. New waves have |m|<=2, and the new
quadratic residual has |m|<=4. No distinct-phase cancellation is presumed.

## 4. Recompute ALL covariances and mean pressure

Let W_ab=<w_a w_b>_theta and recompute it with w+delta w. Keep the full
six-entry symmetric tensor, cross terms and new square. Its increment is
M_C: old W_1/2 times new W_B has order C, and 2B>=C. Generic complex
quadratures remove the initial real cosine/sine improvement. If a,b are
old cosine/sine coefficients and c,d are new ones, exactly

    delta W = (a c^T+c a^T+b d^T+d b^T+c c^T+d d^T)/2.

Matching angular frequencies are combined before averaging; different
frequencies average to zero. The actual separated labels have no cross
products. The finite algebra fixture uses this FULL expression, not the
initial phase-locked signed identity.

With bars denoting the auxiliary Haar mean, reconstruct the compact mean
pressure from the changed radial source, including its moment bump:

\[
\begin{split}
\Delta\bar g_r&=-(\partial_R+1/R)\Delta\bar W_{rr}
 -\epsilon\partial_Z\Delta\bar W_{rz}+\Delta\bar W_{\theta\theta}/R,\\
\Delta P&=\int\Delta\bar g_r\,dR,\qquad
\partial_R\Delta\bar p_m=\Delta\bar g_r-\rho\Delta P,\\
\int R\Delta\bar p_m\,dR&=-\tfrac12\int R^2\Delta\bar g_r\,dR+c_\rho\Delta P,\quad
c_\rho=\tfrac12\int R^2\rho\,dR,\\
\Delta J_\theta&=\int R^2\Delta\bar W_{\theta z}\,dR,\qquad
\Delta J_z=\int R\Delta\bar W_{zz}\,dR-\tfrac12\int R^2\Delta\bar g_r\,dR.
\end{split} \tag{5}
\]

Differentiate the WHOLE c_rho Delta P product in the axial integrated
equation. The unaveraged pressure reconstruction uses the same full radial
source. Its pointwise change and full tangential mean residual changes are
M_(C-kappa); averaging removes graph derivatives and gives M_C. Thus the
new full E_theta,E_z are M_(C-kappa) and all three recomputed defects are
S_C. Their previous orders 3-6kappa and 3-3kappa are NOT preserved through
this wave step. The radial residual is the recomputed moment obstruction
-rho P plus the actual flat A0 contribution. The compact-pressure sign
convention in (5) is the same as the frozen mean equations.

Wave corrections have angular mean zero, so both mean-velocity moments
are unchanged, as are the nonzero old slow Haar mean and current velocity
orders mu/mu+1. Full p_m remains M_mu because C-kappa>mu. Both background
pressure tails are still present; they cancel only inside the SAME full
background equation, never by replacing p_B or the axial pressure moment.
Actual A0/A1 flat velocities stay in EVERY product. Flat ledgers are updated
by exact subtraction when a mean reconstruction moves a term, not counted
twice. The parent mean bounds plus (5) now join the ONE finite regular
normalized residual at min(D,C-kappa,C)=D, gain 1/2-3kappa>0 over B.
This uses all the stated analytical lemmas and excludes the retained exact
flat summands; it is not the remaining all-stage residual theorem.

## 5. Physical bounds and remaining obligation

For all fixed nonnegative spatial/time derivative counts a,b, the inherited
band transfer applies to this actual first sourced correction:

\[
\begin{array}{c|c}
\text{quantity}&\text{power of }Q\text{ (up to fixed polynomial }S\text{)}\\
\delta u&-A+hB-a(1/2+h/2)-b(1+3h/2)\\
\delta p&-2A+h(B+1/2)-a(1/2+h/2)-b(1+3h/2)\\
\delta\mathcal A&1/2-A+h(B+1/2)-a(1/2+h/2)-b(1+3h/2)
\end{array} \tag{6}
\]

Support is the original closed enlarged shell/slow box/rectangle with
smooth temporal extension. Bounds hold on EVERY sufficiently late scale,
including the inherited enlarged-chart margin 4Q<qstar. The h=1/200
diagnostic zero-jet powers are -0.50000015, -1.00250015, +0.00249985.
The residual power -2A-1/2+hD=-1.5025003 is still negative: normalized
order improvement is NOT physical absolute decay or a fixed-scale norm
contraction. Derivative/harmonic constants need not be uniform in iteration.

Crucially w_new-W0^as is now only W_B, not W_(1-kappa): the complex
sourced wave dominates the old curl remainder. The NEXT obligation is to
use that rho=B in the GENERIC complete-wave signed covariance inverse,
then reapply current fast/slow mean maps, recomputing pressure and defects
after EACH operation and preserving the tail exclusion. The initial-real
shortcut cannot be reused. Later generated-source bounds, every finite
cycle, infinite potentials, forcing/energy and blow-up remain unverified.

## 6. Reproduction and scope of diagnostics

`sources.py` checks the full Cartesian cutoff commutator with manufactured
compact jets, exact rational forced moving-normal equations for +/-1,+/-2,
the Gaussian early-source obstruction, general two-quadrature covariance
and moving axial pressure moments, and an exact rational order/physical
ledger. These checks detect omitted curl, pressure, viscosity and transport
terms. They do NOT numerically evaluate the actual first band or prove the
continuum inequalities. `audit.py` recursively reproduces the parent,
checks 640 inherited hashes and binds the restricted claims to source hashes.

```sh
python research/navier_stokes_cascade/results/first_sourced_wave_update/test_audit.py -v
python research/navier_stokes_cascade/results/first_sourced_wave_update/audit.py \
  --verify-record research/navier_stokes_cascade/results/first_sourced_wave_update/evidence.json \
  --output first_sourced_wave_replay.json
```

Use a fresh output path. Frozen records are never overwritten.
