# Post-signed full angular-mean relaxation

Parent: `a68e4e5d204e51dcc9ec6d561e10fecaf4985445` (first real signed
covariance update). Its [CI run](https://github.com/johnchan713/maths_physics_showcase/actions/runs/38912643624)
passed, checked 2026-10-16. All historical source and evidence bytes remain
unchanged. Keep the same fixed Md64, h, N, positive lambda, `C^-20` width,
leading construction, common geometry and accepted later-band domain.

Using the retained analytical mean/pulse lemmas, the next fast/slow update
brings BOTH full angular-mean tangential residuals, before auxiliary Haar
averaging, to `M_(3-6kappa)`. The recomputed compatibility defects belong
to `S_(3-3kappa)`. The radial mean residual has that stronger order plus
its retained flat cutoff error. This closes this post-signed MEAN step,
not the nonzero-angular residual, a complete later cycle or blow-up.

## 1. Current state and the nonzero OLD slow mean

Use the fixed Q-chart, `epsilon=Q^h`, `S=ell^2`, `kappa=1/100000`,
mean operators `D_r`, `D_z=epsilon partial_Z`,
`t_*=-epsilon partial_T+c_i N_t` and `D_e=D_r+e/R` from
[compact means](../compact_mean_corrections/README.md). Ordinary fixed
coefficient jets preserve M orders, radial graph operators cost kappa
each, and axial/slow-time operators gain one. Put

\[
\alpha=2-3\kappa,\quad H=\alpha-\kappa=2-4\kappa,
\quad\nu=3-6\kappa,\quad\mu=1-\kappa.
\tag{1}
\]

The [parent](../real_signed_covariance_update/README.md) supplies the
actual updated complete curl wave `w in W_(1/2)`, full tangential mean
residuals `M_H`, their auxiliary averages `M_nu`, and CURRENT defects
`(P,J_theta,J_z) in S_alpha`. The total tangential means are `M_mu`,
radial mean `M_(mu+1)`, with zero weighted velocity moments. The full wave
tensor `W_ab=<w_a w_b>_theta` keeps all six entries; HOLD w fixed during
both mean operations below. Angular mean and auxiliary Haar mean differ.

There is additional exact structure inherited from initialization:
`m=m_s+m_f`, where `m_f` has zero Haar mean, while `m_s` is the original
auxiliary-independent slow five-map velocity. The first signed WAVE
update did not change either mean velocity. Consequently

\[
\bar v,\bar\gamma\in M_1,\qquad \bar\beta\in M_2.
\tag{2}
\]

These averages are NOT zero. The zero-Haar argument that gave the initial
`S_1 -> S_2` compatibility gain cannot set them to zero in this update.
All estimates mean every fixed coefficient derivative order, with fixed
construction constants and polynomial S factors, not a numerical solve
of the giant actual profile or a proof-assistant certificate.

The actual summed base has `b=O(epsilon)`, `V,G=O(1)` on the active shell.
On the convex hull of the reserved I4 slow-potential supports its exact
values remain `G=0`, `V_q=a(eta)x^(-1-2lambda)`, `a>=a0>0`.
The support and weighted-jet checkpoints supply these actual inputs;
the fixed compact inverse, two-torus inverse and five-map estimates are
explicitly retained analytical lemmas. No new cone, roots or parameters
are chosen. Inverses may have huge finite constants at the fixed lambda;
no uniformity as lambda tends to zero is asserted.

## 2. Fast map at the CURRENT order H, then recompute

Use the current complete tangential residuals in

\[
u_v=-c_i^{-1}N_t^{-1}E_\theta^\circ,\quad
\gamma_d=-c_i^{-1}N_t^{-1}E_z^\circ,\quad
\Psi=T_1\gamma_d,\quad
u_\beta=-\epsilon\partial_Z\Psi,\quad
u_\gamma=D_1\Psi=\gamma_d-A_1\gamma_d.
\tag{3}
\]

The temporal inverse gives zero Haar mean pointwise in (R,Z,T).
Haar invariance of the compact radial inverse gives zero Haar mean for
Psi and the ACTUAL three increments, including A1. Thus
`u_v,u_gamma in M_H`, `u_beta in M_(H+1)`, divergence and both velocity
moments are preserved, and (2) is unchanged EXACTLY. Keep the actual
A1 velocity in all products. Its uncancelled fast-time derivative is a
separate flat residual, not zero. Mean support can fill the auxiliary
torus, but remains in the same shell/slow projection with smooth extensions.

Recompute `g_r`, `p_m=T_0(g_r-rho P)` and ALL three defects from m+u.
With `u=(u_beta,u_v,u_gamma)`, exact Haar averaging gives

\[
\begin{aligned}
\Delta\bar g_r={}&-(\partial_R+1/R)
 [2\overline{\beta u_\beta}+\overline{u_\beta^2}]\\
 &-\epsilon\partial_Z\overline{
 \beta u_\gamma+\gamma u_\beta+u_\beta u_\gamma}
 +[2\overline{v u_v}+\overline{u_v^2}]/R.
\end{aligned}\tag{4}
\]

Linear base, time and viscosity changes average to zero, but correlations
between the OLD fast mean and NEW fast increment need not vanish.
Their leading order is `mu+H=3-5kappa`, better than alpha. The flux changes
are exactly

\[
\Delta J_\theta=\int R^2\overline{\gamma u_v+v u_\gamma+u_\gamma u_v},
\quad
\Delta J_z=\int R\overline{2\gamma u_\gamma+u_\gamma^2}
             -\tfrac12\int R^2\Delta\bar g_r.
\tag{5}
\]

All three defect changes are `S_(mu+H)`; current defects therefore remain
`S_alpha`. They are changed, NOT frozen at their previous values.
Before Haar averaging the radial source/pressure changes are `M_H`:
the term `2V u_v/R` is the lowest, while u_beta has an extra power.
The fixed compact pressure map preserves this order; its axial effect
has order H+1. After principal fast-time cancellation the full tangential
error terms have the following sufficient orders:

| Term | Order |
|---|---:|
| slow time | H+1 |
| base radial transport, induced radial velocity | H+1-kappa |
| radial viscosity | H+1-2kappa = nu |
| old mean transport | mu+H+1-kappa |
| changed axial pressure | H+1 |
| new mean self-products | 2H+1-kappa |

Every row is at least nu. Combining with the old averaged residual
`M_nu` gives FULL tangential residuals `M_nu` after (3), modulo the
retained flat summands. Flat velocities remain in (4)--(5), even though
their values are arbitrarily small in asymptotic bounds.

## 3. Slow five-map at CURRENT defects, keeping the old Haar mean

Use the same fixed I4 five-map on the RECOMPUTED `(P,J_theta,J_z)` at
order alpha. Its two moment rows are zero and its last three rows are

\[
\int (2V/R)d_v=-P,\quad
\int R^2(Gd_v+Vd_\gamma)=-J_\theta,\quad
\int(2RGd_\gamma-RVd_v)=-J_z.
\tag{6}
\]

It gives slow `d_v,d_gamma in M_alpha`, `d_beta in M_(alpha+1)`, with
the EXACT curl realization of d_gamma and zero velocity moment changes.
Its potential support includes the intervals between bump supports;
that convex hull stays inside I4. Reconstruct pressure and defects AGAIN.

Let B=bar(beta), v0=bar(v), g0=bar(gamma) at the postfast state. By (2)
their orders are (2,1,1). For `R_g=g_r,new-g_r-(2V/R)d_v`, the exact
averaged remainder is

\[
\begin{aligned}
\bar R_g={}&\epsilon\partial_Td_\beta
 -(\partial_R+1/R)(2bd_\beta+2Bd_\beta+d_\beta^2)\\
 &-\epsilon\partial_Z(bd_\gamma+Gd_\beta+Bd_\gamma+g_0d_\beta+d_\beta d_\gamma)
 +(2v_0d_v+d_v^2)/R\\
 &+\epsilon(\partial_{RR}+R^{-1}\partial_R+\epsilon^2\partial_{ZZ}-R^{-2})d_\beta.
\end{aligned}\tag{7}
\]

The fast parts drop from these SLOW cross products only after Haar
averaging. Crucially `2v0 d_v/R` survives at order alpha+1. Time, base
flux and viscosity have order at least alpha+2; old radial mean products
have order at least alpha+3; the new square has order 2alpha. Because
alpha>1, every term is at least alpha+1, not generally 2alpha.

After cancellation of the CURRENT linear rows (6), the EXACT defects are

\[
P_{new}=\int\bar R_g,\quad
J_{\theta,new}=\int R^2(g_0d_v+v_0d_\gamma+d_\gamma d_v),
\]
\[
J_{z,new}=\int R(2g_0d_\gamma+d_\gamma^2)
                 -\tfrac12\int R^2\bar R_g.
\tag{8}
\]

Thus ALL THREE lie in `S_(alpha+1)=S_(3-3kappa)`; they need not be
identically zero. Dropping the old slow mean would falsely promote this
to order 2alpha. Using prefast targets instead of CURRENT ones would
leave the nonzero changes (4)--(5) uncancelled.

The full tangential change under this slow map has orders at least
alpha+1: slow time, base fluxes, viscosity and changed axial pressure.
Products with the total old means have order at least
`mu+alpha+1-kappa`, and the new square gives `2alpha+1`.
All exceed nu. Therefore the final FULL angular-mean tangential rows
remain `M_nu`; the radial residual is `-rho P_new` plus its actual flat
pressure-cutoff term, with scalar order `alpha+1>nu`.

## 4. Complete pressure memory, supports and the uncancelled wave residual

At every state keep the entire six-entry covariance of the updated wave
and use the SAME fixed physical pressure bump. The exact identities are

\[
\int R\bar p_m=-\tfrac12\int R^2\bar g_r+c_\rho P,\quad
\int R^2\bar E_\theta=\epsilon\partial_ZJ_\theta,\quad
\int R\bar E_z=\epsilon\partial_Z(J_z+c_\rho P).
\tag{9}
\]

The whole moving product remains inside partial_Z. Neither W_zz nor
W_theta z nor the axial mean velocity vanishes by leading polarization.
Both background pressure tails remain in p_B and cancel only against
that SAME full base pressure, not a truncated leading pressure. The flat
base residual stays a separate summand; it need not have shell support
or the shell edge weight. The nonzero A0/A1 errors are never erased
inside later products. Both means represent the same physical maps on
chart overlaps, by the retained compact/temporal inverse lemmas.

During these mean changes the nonzero-angular residual changes ONLY by
`(w.grad)d+(d.grad)w`. The mean--wave half-frequency loss is retained:
`w in W_(1/2)` and a tangential mean increment of order t, with radial
order t+1, give a source in W_t. Hence these two steps add wave sources
`W_H` and `W_alpha`, conserving each wave's closed support, smooth endpoint
extension, Gaussian P envelope and half-edge weight. Means may fill the
auxiliary torus, but a nonzero source retains its wave factor. This does
not cancel the OLD nonzero-angular residual. No complete residual-cycle
gain or general later complex-amplitude covariance estimate follows.

For every fixed physical space/time jet (a,b), the same harmonic losses
`sx=1/2+h/2`, `st=1+3h/2` dominate coefficient graph/frame losses.
An increment of mean order t=H or alpha has sufficient physical powers
`-A+ht-a sx-b st` for velocity and `-2A+ht-a sx-b st` for its pressure
change. The radial/axial curl potential Psi has power
`1/2-A+ht-a sx-b st`. This is ONLY that azimuthal potential; a complete
potential realization of the direct swirl and an infinite sum are not
asserted. All Q, covers and discrete choices stay fixed in differentiation;
the estimates retain their weights and polynomial S factors on every
accepted later scale. Physical residual exponents can still be negative.

## 5. Reproduction and remaining scope

The rational diagnostics use radially compact finite-regularity polynomial
means with an explicit ONE-circle graph derivative, a full six-entry
updated curl covariance and a moving pressure bump. They independently
check temporal inversion, graph curl divergence, nonzero old/new fast
correlations, exact slow defect differences, nonzero old slow cross terms,
velocity moments and pressure/flux identities at all three states.
They do not evaluate the actual two-torus inverse, actual five-map targets
or giant continuum profile. The inherited two-torus/cutoff and fixed
five-map diagnostics are replayed; actual all-orders estimates are the
analytical consequences above, not numerical proof. All 625 inherited
source digests are preserved.

Next: audit the CURRENT nonzero-angular residual and its sourced pulse
inverse, then the later general complex-wave covariance/mean cycles,
without using the special first-real-update shortcut. Full residual
improvement, infinite realization, forcing, energy and blow-up stay open.

```sh
python3 research/navier_stokes_cascade/results/post_signed_mean_relaxation/test_audit.py -v
python3 research/navier_stokes_cascade/results/post_signed_mean_relaxation/audit.py \
  --verify-record research/navier_stokes_cascade/results/post_signed_mean_relaxation/evidence.json \
  --output /tmp/post_signed_mean_relaxation_replay.json
```
