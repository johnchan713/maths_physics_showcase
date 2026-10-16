# Current complex-wave state: fresh fast/slow mean relaxation

Status: `current-complex-state-full-angular-mean-relaxation-audited`.
Parent: `71965da8e3759b0db0d02b05c298858121c61d10`.
Its [complete CI run](https://github.com/johnchan713/maths_physics_showcase/actions/runs/38936441478)
passed, checked 2026-10-16. This includes the scoped replay adapter and
the previously failing first-sourced-wave step. No historical checkpoint
source or frozen evidence is changed here.

Using the retained analytical mean/inverse lemmas, one CURRENT fast map
followed by a CURRENT slow five-map gives full angular-mean tangential
residuals M_(2-6kappa), before auxiliary Haar averaging, and recomputed
defects S_(5/2-3kappa), plus retained flat errors. The current nonzero
regular wave residual stays W_(3/2-6kappa). That residual still limits the
complete finite-state bound; no further overall order improvement or
all-stage closure is claimed. The sharp initial zero-Haar compatibility
argument does not apply to the accumulated old slow mean.

The [supplied proposed manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
is the source under audit, (8.12)--(8.27), (9.3)--(9.10), not an assumed
proof. SHA256: `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The exact identities and analytical inputs below are the frozen
[compact mean maps](../compact_mean_corrections/README.md),
[previous mean derivation](../post_signed_mean_relaxation/README.md),
[sourced update](../first_sourced_wave_update/README.md) and
[current generic signed state](../complex_signed_covariance_update/README.md).
Manufactured rational checks are not a numerical giant-band construction,
an independent proof of the continuum premises, or a proof-assistant certificate.

## 1. Current inputs and accumulated old Haar means

Keep the same Md64, N, fixed h in (0,1/100), positive lambda, setup width
C^-20, leading construction, exact summed base, supports and discrete
labels. Use epsilon=Q^h, S=ell^2, kappa=1/100000, the common qstar and
enlarged-chart requirement 4Q<qstar on every accepted later scale. W/M/S
orders include all fixed coefficient derivatives, edge weights and fixed
inverse-distance powers; polynomial S losses and finite-stage constants
are retained. A larger exponent is a stronger estimate for epsilon<1.

Put

\[
\alpha=\tfrac32-3\kappa,\quad H=\alpha-\kappa,
\quad\nu=2-6\kappa,\quad\mu=1-\kappa,\quad D=\tfrac32-6\kappa. \tag{1}
\]

The CURRENT complete complex wave w is W_1/2; its remainder about the
unchanged original primary is W_(1-3kappa). The full tangential mean rows
are M_H, their auxiliary Haar averages M_nu, and current defects
(P,J_theta,J_z) are S_alpha. All six entries of W_ab=<w_a w_b>_theta
are held FIXED through both mean operations here. Angular and auxiliary
averaging are distinct. Full tangential mean velocities have order mu,
radial mean mu+1; their weighted velocity moments are zero.

Every preceding fast velocity has exact zero Haar mean, including its
actual compact curl remainder. Every preceding slow update is
auxiliary-independent. Thus the current old Haar velocity is the SUM of
the original slow mean and the preceding slow increment, not the original
mean alone. The latter increment has order 2-3kappa>1, so the accumulated
orders remain

\[
\bar v,\bar\gamma\in M_1,\qquad \bar\beta\in M_2. \tag{2}
\]

They are NOT zero. Both intervening wave updates left mean velocities
unchanged. These exact facts justify using (2) now, rather than importing
an initial mean state. All coefficients of the actual background have
their audited fixed-jet bounds; on the reserved I4 potential-support
convex hull G=0 and V is the same nonzero power-law velocity. The original
compact/temporal inverses and fixed five-map remain analytical premises;
no new cone, radius, roots or frequency is chosen.

## 2. Fast map at H and postfast pressure/defects

Apply the zero-Haar temporal inverse to the CURRENT full tangential rows:

\[
u_v=-c_i^{-1}N_t^{-1}E_\theta^\circ,\quad
\gamma_d=-c_i^{-1}N_t^{-1}E_z^\circ,\quad
\Psi=T_1\gamma_d,\quad
u_\beta=-\epsilon\partial_Z\Psi,\quad
u_\gamma=D_1\Psi=\gamma_d-A_1\gamma_d. \tag{3}
\]

Haar invariance of T1 and zero-mean inversion give zero Haar mean for all
three ACTUAL increments, including A1. They have orders H,H,H+1 for
tangential/tangential/radial components, preserve divergence and both
velocity moments, and leave (2) exactly unchanged. Keep A1 in every
nonlinear product. Its remaining fast-time error is flat, not zero.

Recompute radial source, compact pressure and all three defects from m+u.
Writing m=(beta,v,gamma), the exact changes are

\[
\begin{split}
\Delta\bar g_r={}&-(\partial_R+1/R)
 (2\overline{\beta u_\beta}+\overline{u_\beta^2})\\
&-\epsilon\partial_Z\overline{\beta u_\gamma+\gamma u_\beta+u_\beta u_\gamma}
 +(2\overline{v u_v}+\overline{u_v^2})/R,\\
\Delta J_\theta={}&\int R^2\overline{\gamma u_v+v u_\gamma+u_\gamma u_v},\\
\Delta J_z={}&\int R\overline{2\gamma u_\gamma+u_\gamma^2}
                   -\tfrac12\int R^2\Delta\bar g_r .
\end{split} \tag{4}
\]

Linear base, time and viscosity increments average to zero; OLD/NEW FAST
correlations do not. Their leading defect-change order is mu+H=5/2-5kappa,
so current defects remain S_alpha but CHANGE. Unaveraged radial source and
pressure changes are M_H, with axial pressure effect H+1. The sufficient
post-cancellation full tangential remainder orders are

| Term | Order |
|---|---:|
| Slow time | H+1 |
| Base radial transport and induced radial velocity | H+1-kappa |
| Radial viscosity | H+1-2kappa = 5/2-6kappa |
| Old mean transport | mu+H+1-kappa |
| Changed axial pressure | H+1 |
| New mean self-products | 2H+1-kappa |

Every row is STRICTLY above nu. Thus the old averaged covariance residual
M_nu, not fast viscosity, limits the postfast FULL tangential rows to M_nu,
modulo exact flat summands. This differs from the earlier real signed mean
checkpoint, where viscosity and the averaged target had the same order.

## 3. Slow map at the RECOMPUTED postfast defects

Use the same fixed five-map at order alpha, with two zero velocity-moment
rows and the three current linear cancellation rows

\[
\int (2V/R)d_v=-P,\quad
\int R^2(Gd_v+Vd_\gamma)=-J_\theta,\quad
\int (2RGd_\gamma-RVd_v)=-J_z. \tag{5}
\]

The slow d_v,d_gamma are M_alpha and d_beta is M_(alpha+1), with exact
curl realization and smooth potential support in the same I4 convex hull.
Recompute pressure and defects AGAIN. Let B=bar beta, v0=bar v, g0=bar gamma
at the postfast state. If R_g=g_r,new-g_r-(2V/R)d_v, exactly

\[
\begin{split}
\bar R_g={}&\epsilon\partial_Td_\beta
 -(\partial_R+1/R)(2bd_\beta+2Bd_\beta+d_\beta^2)\\
&-\epsilon\partial_Z(bd_\gamma+Gd_\beta+Bd_\gamma+g_0d_\beta+d_\beta d_\gamma)
 +(2v_0d_v+d_v^2)/R\\
&+\epsilon(\partial_{RR}+R^{-1}\partial_R+\epsilon^2\partial_{ZZ}-R^{-2})d_\beta.
\end{split} \tag{6}
\]

Only zero-Haar parts disappear against slow d after averaging. The
accumulated v0 term survives at alpha+1; the new square has order 2alpha,
time/base/viscosity at least alpha+2 and old radial mean terms alpha+3.
Consequently the exact recomputed defects, after the current rows (5), are

\[
P_{new}=\int\bar R_g,\quad
J_{\theta,new}=\int R^2(g_0d_v+v_0d_\gamma+d_\gamma d_v),
\]
\[
J_{z,new}=\int R(2g_0d_\gamma+d_\gamma^2)-\tfrac12\int R^2\bar R_g,
\qquad (P_{new},J_{\theta,new},J_{z,new})\in S_{\alpha+1}. \tag{7}
\]

Here alpha+1=5/2-3kappa<2alpha: the nonlinear defects need not vanish.
Using prefast targets instead would leave (4), whose available order
mu+H=5/2-5kappa is BELOW alpha+1 by 2kappa. This precise inequality blocks
the claimed final defect estimate with stale targets. Smooth smallness or
passing a fixture is not a replacement for this exact recomputation.

Full slow tangential changes have orders alpha+1 (slow time, base flux,
viscosity, changed axial pressure), mu+alpha+1-kappa (old means), and
2alpha+1 (new square), all above nu. The final FULL angular-mean tangential
rows therefore remain M_nu. The radial mean residual is the new moment
obstruction -rho P_new plus its actual flat pressure-cutoff term, with
stronger scalar order alpha+1. The old accumulated Haar velocity remains
M_1/M_2, and total means remain M_mu/M_(mu+1).

## 4. Full pressure memory and the remaining nonzero residual

At all three states use the complete FIXED complex-wave tensor and the
same full background pressure. The exact compact identities are

\[
\int R\bar p_m=-\tfrac12\int R^2\bar g_r+c_\rho P,\quad
\int R^2\bar E_\theta=\epsilon\partial_ZJ_\theta,\quad
\int R\bar E_z=\epsilon\partial_Z(J_z+c_\rho P). \tag{8}
\]

Differentiate the WHOLE moving product. Neither axial velocity nor W_zz,
W_theta z or p_m is deleted by leading polarization. Both exceptional p_B
tails cancel only in subtraction of the SAME full base equation. Actual
A0/A1 velocities remain in every product. Flat base residual is separate
and need not have shell support or shell edge weight. Flat residual
ledgers moved by reconstruction are changed by exact subtraction, not
counted twice. No designated Gaussian tail is forwarded to a pulse inverse.

With w fixed, the only nonzero-angular changes are
(w.grad)d+(d.grad)w, with ACTUAL increments, including flat velocities.
The mean-wave half-frequency loss gives W_t for a mean increment of order
t since w is W_1/2 and radial mean increments gain one. Hence the added
sources are W_H and W_alpha, both above D. A nonzero product keeps a wave
factor, its Gaussian P envelope, half-edge weight, closed slow/transverse/
temporal supports and smooth endpoint extension even if its mean fills
the torus. No new harmonic is introduced beyond the current finite set.
Thus the old W_D regular nonzero residual remains W_D; it is not cancelled.

For every fixed physical spatial/time jet (a,b), mean increments of order
t=H or alpha have sufficient powers -A+ht-a sx-b st for velocity,
-2A+ht-a sx-b st for pressure, and 1/2-A+ht-a sx-b st for the radial/axial
azimuthal curl potential. Here A=1/2+h, sx=1/2+h/2, st=1+3h/2; edge weights
and polynomial S factors are retained across ALL accepted later bands.
This partial potential is not a full direct-swirl or infinite realization.
At h=1/200 the fast zero-jet powers are -0.4975002, -1.0025002, +0.0024998;
the slow powers are -0.49750015, -1.00250015, +0.00249985. The complete
regular residual power remains -1.5025003, so normalized mean gain is not
absolute physical decay or a fixed-scale contraction.

## 5. Diagnostics and next obligation

`relaxation.py` starts with the manufactured means AFTER the previous slow
step, and the complete complex tensor AFTER the new generic signed step.
It then makes fresh fast and slow compact curl increments. Exact rational
identities check current correlations, accumulated Haar means, all three
pressure reconstructions, both integrated fluxes, six velocity moments and
the nonzero old slow cross terms. The fixture uses a one-circle graph and
finite endpoint regularity; it is not the actual two-torus inverse, an
actual numerical five-map solve or an evaluation of the giant profile.
Actual compact/five-map estimates remain explicitly inherited analytical
lemmas. All 671 inherited source digests and the parent's scoped roundoff
adapter are replayed without modifying historical evidence.

Next use the CURRENT nonzero REGULAR source of order D=3/2-6kappa for the
next sourced solve, retaining the exact excluded tail ledger. Join its new
source, covariance and pressure estimates to a full finite-stage induction;
the fixed-primary old-wave remainder is still only W_(1-3kappa). Later
stages, all-stage residual closure, infinite potentials, forcing/energy and
blow-up remain unverified.

```sh
python research/navier_stokes_cascade/results/complex_mean_relaxation/test_audit.py -v
python research/navier_stokes_cascade/results/complex_mean_relaxation/audit.py \
  --verify-record research/navier_stokes_cascade/results/complex_mean_relaxation/evidence.json \
  --output complex_mean_replay.json
```

Use a fresh output path; frozen records are never overwritten.
