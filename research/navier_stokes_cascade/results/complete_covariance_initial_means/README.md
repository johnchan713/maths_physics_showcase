# Actual-primary complete-covariance INITIAL mean update

Parent: `0cb55230b7b80332266cab01daed83d3e336c79e` (physical primary-pulse
transfer). Its [CI run](https://github.com/johnchan713/maths_physics_showcase/actions/runs/38898489926)
completed successfully, checked 2026-10-16. All historical records are
unchanged. This checkpoint sharpens and supplies the actual inputs for the
INITIAL mean step in the older conditional cycle ledger; it does not
certify that ledger's entire nonlinear iteration.

Result, using the retained analytical lemmas: hold the actual complete
primary curl wave fixed. Starting with zero mean velocity, apply the
compact pressure map, one zero-Haar fast update and one slow five-map to
the RECOMPUTED defects. The final three defects belong to `S_2`. The
angular-mean annular residual belongs to `M_(2-3kappa)` plus the explicitly
retained flat errors. No new covariance wave has yet been added. This is
not a bound for the complete nonzero-angular residual or later cycles.

## 1. Actual inputs, notation and scope

Use exactly the fixed `Md=64`, positive lambda, h, N, width `C^-20`,
common geometry, covers and sufficiently late physical bands from the
[parent](../physical_pulse_transfer/README.md). On each fixed Q-chart,
`epsilon=Q^h`, `S=ell^2`, `kappa=1/100000`. All derivative statements mean
EVERY fixed coefficient derivative order, with polynomial S losses and
constants depending on the fixed construction; neither the first band nor
the actual enormous profile coefficients are numerically evaluated.

`M_alpha` and `S_alpha` are the inherited weighted shell/normalized-moment
classes with factor epsilon^alpha; `W_alpha` is the wave class. Radial
graph differentiation costs kappa, axial differentiation adds one power:

\[
D_r=\partial_R+M d_rR^{d_r-1}L,\quad D_e=D_r+e/R,
\quad D_z=\epsilon\partial_Z,\quad t_*=-\epsilon\partial_T+c_iN_t.
\]

An overbar is auxiliary Haar average, NOT the physical angular average.
For the actual complete primary `w`, set `W_ab=<w_a w_b>_theta`. Keep all
six entries, including the complete curl remainder. Separated closed wave
supports and the all-order primary bounds give `W_ab in M_1`. The real
primary sine/cosine cancellation gives
`bar W_(r theta,rz)-Sigma_0 in M_(2-2kappa)`; the actual summed background
stress satisfies `Sigma-Sigma_0 in M_3` by the relative `O(q^(2h))`
comparison. These use the actual-background/geometry-to-pulse derivations,
not a finite scan. The primary covariance statement is special to real
primaries and is not silently used for subsequent complex waves.

Use the exact summed slow base `(b,V,G)`, not just its leading profile.
On the convex radial hull of the reserved I4 mean supports, the
[support audit](../coefficient_support_audit/README.md) gives
`G=0`, `V_q=a(eta)x^(-1-2lambda)`, `a>=a0>0`, and the local bound
`b=O(epsilon)` with all fixed jets. In particular b is NOT zero.
All radial moment operations preserve slow support, have smooth shell
extensions and incur only polynomial S costs; means may fill the
auxiliary torus. The compact inverse, temporal inverse, moving-support
and fixed five-map estimates are the explicitly retained analytical
lemmas from [compact means](../compact_mean_corrections/README.md).

### Both background pressure tails remain

Write total pressure as `p_B+p_m+pi_w`. The primary wave pressure has zero
angular mean, not zero pointwise value. The background on I4 still has

\[
q^{2A}p_B=\Pi_0+\chi_1q^{2h}[K_1-\mu^2/(4X)]
                  +\chi_2q^{4h}K_2,\qquad A=1/2+h.
\tag{1}
\]

Here `D=1/2-h`, `d=1-eta^2`, `L=1-2h eta^2`. Both tails are retained.
Subtracting the SAME exact base residual with
the SAME `p_B` removes its gradient by pressure linearity; it does not
assert either tail vanishes. Replacing `p_B` by leading `p_0` would leave
their actual gradients. For example with `p_2=-2A+4h=-1+2h`, `sigma=c_2q`,

\[
\partial_z[q^{p_2}\chi_2(\sigma)K_2(\eta)]
 ={q^{p_2-D}\over L}
 [2\eta(p_2\chi_2+\sigma\chi_2')K_2+d\chi_2K_2'].
\tag{2}
\]

Thus `Pi_2,X=0` does not imply zero physical axial pressure gradient.
Conversely the particular first-tail centrifugal term on a cutoff
plateau is the z-independent physical function `-mu^2/(2r^2)`; its axial
gradient DOES cancel. A cutoff transition need not cancel. We do not
overclaim nonzero gradients for every individual pressure piece.

## 2. Exact compact pressure and complete axial moment

For `m=(beta,v,gamma)`, use the full conservative mean equations (2) in
the compact-means checkpoint, including cylindrical connections, all
viscosity, and all `W_ab`. In particular

\[
g_r=-t_*\beta-D_1(2b\beta+\beta^2+W_{rr})
 -D_z(b\gamma+G\beta+\beta\gamma+W_{rz})
 +(2Vv+v^2+W_{\theta\theta})/R
 +\epsilon(\Delta_0-R^{-2})\beta,
\quad \Delta_0=D_r^2+R^{-1}D_r+D_z^2.
\tag{3}
\]

At EACH updated velocity reconstruct

\[
P=\int\bar g_r\,dR,\quad p_m=T_0(g_r-\rho P),\quad
D_rp_m-g_r=-\rho P-A_0(g_r-\rho P).
\tag{4}
\]

Here the fixed PHYSICAL bump has chart expression
`rho=s^(-1/2)rhohat(R/sqrt(s))`, `s=q/Q`, and unit radial integral.
The A0 term has zero Haar mean and is flat to any fixed requested
order; it remains in the actual pressure and residual. Define

\[
J_\theta=\int R^2\overline{Gv+V\gamma+\gamma v+W_{\theta z}}\,dR,
\quad J_z=\int R\overline{2G\gamma+\gamma^2+W_{zz}}\,dR
              -\tfrac12\int R^2\bar g_r\,dR.
\tag{5}
\]

With the exact two velocity moments zero, the identities are

\[
\int R^2\bar E_\theta\,dR=\epsilon\partial_ZJ_\theta,\quad
\int R\bar E_z\,dR=\epsilon\partial_Z(J_z+c_\rho P),
\quad c_\rho=\tfrac12\int R^2\rho\,dR.
\tag{6}
\]

Indeed `int R bar p_m=-1/2 int R^2 bar g_r+c_rho P` by radial
integration by parts. The WHOLE moving product remains inside partial_Z;
`(partial_Z c_rho)P` is generally nonzero. These formulas retain axial
velocity, `W_theta z`, `W_zz`, `W_rz`, the full pressure moment and both
background tails via the exact base subtraction. Dropping any one is not
justified by leading transverse-wave polarization. PDF pp.94 and99,
equations (8.15)--(8.17),(8.26)--(8.27), were rendered and inspected to
check the signs (source SHA256 in protocol).

## 3. Sharper actual initial raw defects

Initially `m=0`. Pointwise `g_r^0 in M_(1-kappa)`, but Haar translation
invariance gives EXACTLY `bar D_r f=partial_R bar f` and
`bar t_*f=-epsilon partial_T bar f`. Hence `bar g_r^0,bar p_m^0 in M_1`.
Compact integration by parts gives

\[
\begin{aligned}
P^0&=\int{\bar W_{\theta\theta}-\bar W_{rr}\over R}\,dR
       -\epsilon\partial_Z\int\bar W_{rz}\,dR,\\
J_\theta^0&=\int R^2\bar W_{\theta z}\,dR,\\
J_z^0&=\int R[\bar W_{zz}-\tfrac12(\bar W_{rr}+\bar W_{\theta\theta})]\,dR
       +\tfrac\epsilon2\partial_Z\int R^2\bar W_{rz}\,dR.
\end{aligned}\tag{7}
\]

All three are in `S_1`, without a graph loss. The old conservative
`S_(1-kappa)` remains valid, but is unnecessary for this initial step.
The two Haar-averaged tangential residuals are

\[
\bar E_\theta^0=(\partial_R+2/R)(\bar W_{r\theta}-\Sigma_\theta)
                  +\epsilon\partial_Z\bar W_{\theta z},\quad
\bar E_z^0=(\partial_R+1/R)(\bar W_{rz}-\Sigma_z)
                  +\epsilon\partial_Z(\bar W_{zz}+\bar p_m^0).
\tag{8}
\]

Thus both are in `M_(2-2kappa)`; the pointwise zero-Haar residuals are
only in `M_(1-kappa)`. Averaging removes the radial graph operator in
(8), not its earlier losses in the actual curl covariance.

## 4. Initial fast update: exact zero Haar means

Set `v_f=-c_i^-1 N_t^-1(E_theta^0)^circ`,
`gamma_d=-c_i^-1 N_t^-1(E_z^0)^circ`,
`Psi_f=T_1 gamma_d`, `beta_f=-epsilon partial_Z Psi_f`,
`gamma_f=D_1 Psi_f=gamma_d-A_1 gamma_d`.
The inverse acts on zero-mean functions. Haar invariance of the radial
primitive and zero weighted integral give

\[
\bar\beta_f=\bar v_f=\bar\gamma_f=0
\tag{9}
\]

exactly, including the actual flat A1 correction. Divergence and the two
velocity moments vanish exactly. The tangential classes are
`M_(1-kappa)` and the radial class is `M_(2-kappa)`. Keep the A1 velocity
in every subsequent product; its fast-time residual is retained flat.

Recompute (3)--(5) with this SAME full primary wave. All linear changes
have zero Haar mean; in particular

\[
\bar g_r^f-\bar g_r^0
 =-(\partial_R+1/R)\overline{\beta_f^2}
  -\epsilon\partial_Z\overline{\beta_f\gamma_f}
  +\overline{v_f^2}/R\ \in M_{2-2\kappa}.
\tag{10}
\]

The other flux changes are quadratic too, so current defects still belong
to `S_1`. Fast-time cancellation leaves FULL tangential mean residuals
in `M_(2-3kappa)` modulo retained flat terms: the limiting term is
`epsilon D_r^2 v_f` of order `1+(1-kappa)-2kappa`.
Slow time, b-fluxes, beta_f V/G, and the changed axial pressure are no
worse. This retains viscosity rather than promoting the bound to order2.

## 5. One slow five-map, then recompute AGAIN

Use the fixed I4 map with current targets `(P^f,J_theta^f,J_z^f) in S_1`:

\[
\int R^2\delta v=\int R\delta\gamma=0,\quad
\int(2V/R)\delta v=-P^f,\quad
\int R^2(G\delta v+V\delta\gamma)=-J_\theta^f,\quad
\int(2RG\delta\gamma-RV\delta v)=-J_z^f.
\tag{11}
\]

The angular powers `(2,-2-2lambda,-2lambda)` and axial powers
`(1,1-2lambda)` have invertible geometric-bump Vandermonde maps at fixed
positive lambda. The inverse has no small-target/sign restriction, but
its constants need not be uniform as lambda tends to zero.
It gives slow `delta v,delta gamma in M_1`, `delta beta in M_2`.
The compact curl realizes delta gamma EXACTLY, preserving divergence,
moments and support on the convex radial hull inside I4.

Keep w fixed, update the ACTUAL velocity and reconstruct pressure again.
Write `R_g=g_r^new-g_r^f-(2V/R)delta v`. Because of (9), its exact
Haar average simplifies to

\[
\begin{aligned}
\bar R_g={}&\epsilon\partial_T\delta\beta
 -(\partial_R+1/R)(2b\delta\beta+(\delta\beta)^2)\\
 &-\epsilon\partial_Z(b\delta\gamma+G\delta\beta+\delta\beta\delta\gamma)
 +(\delta v)^2/R
 +\epsilon(\partial_{RR}+R^{-1}\partial_R+\epsilon^2\partial_{ZZ}-R^{-2})\delta\beta.
\end{aligned}\tag{12}
\]

Old fast times slow cross terms vanish only AFTER Haar averaging. They
are not zero pointwise. Slow time, b radial flux, axial flux and viscosity
in (12) have orders at least3; the slow square has order2. Every fixed
coefficient derivative obeys the same exponents by the product and
inverse lemmas. On the fixed interior patch the weight has a positive
lower bound; conversion to the normalized scalar moments costs no epsilon.
After (11) the recomputed exact defects are

\[
P^{new}=\int\bar R_g,\qquad
J_\theta^{new}=\int R^2\delta\gamma\delta v,\qquad
J_z^{new}=\int R(\delta\gamma)^2-\tfrac12\int R^2\bar R_g.
\tag{13}
\]

Therefore ALL THREE are in `S_2`, not identically zero. This initial
one-power compatibility gain is sharper than the inherited generic
`alpha+0.9-2kappa` estimate. It requires the exact initial zero-Haar fast
state (9); later means with nonzero averages do not enjoy this cancellation.
Historical evidence and its conservative valid exponents remain untouched.

Changes of tangential residual under this slow update have order at
least2: slow time/transport/pressure cost one epsilon and products with
the fast means are higher. Thus its complete angular-mean increment
residual is bounded by `M_(2-3kappa)`, with radial part in `M_2`, plus
the retained flat pressure/axial reconstruction remainders and the exact
flat base residual. The latter is a SEPARATE flat summand; it need not
have the shell edge weight or shell support. No flat velocity is set to
zero inside nonlinear products. No assertion about the nonzero-angular
wave residual follows from averaging.

## 6. Checks, retained limitations and next obligation

`means.py` supplies exact rational manufactured complete-curl/covariance,
radial pressure moment, recomputed defect and zero-Haar product identities,
with controls that detect missing axial velocity/pressure, a dropped
moving-bump derivative, and invalid extension to a nonzero-old-mean state.
Its polynomial fixtures have finite endpoint regularity, not the actual
C-infinity profiles. The exponent ledger is a rational check of the
analytical derivation above, not a standalone continuum proof. Parent
actual-background/physical bounds and the older compact-map diagnostics
are replayed, preserving 595 inherited source digests.

The next obligation is to realize signed compact covariance targets

\[
H_\theta=-T_2[\bar E_\theta-\sigma_\theta\epsilon\partial_ZJ_\theta],\quad
H_z=-T_1[\bar E_z-\sigma_z\epsilon\partial_Z(J_z+c_\rho P)]
\]

by the fixed primary-wave inverse and bound the COMPLETE updated curl
covariance, axial pressure and all generated nonzero-harmonic sources
through every finite cycle. A signed target is not yet an actual added
wave covariance. Every sourced inverse still needs its own input-class,
Gaussian envelope and smooth endpoint-extension verification; the
mean--wave half-frequency loss and all old/new interactions must remain.
Only then can one claim full residual improvement including cutoff
tails, viscosity and all scales. Forcing, energy, infinite summation and
blow-up remain unverified. No parameters or claims are retuned to obtain
passing tests.

Reproduce (fresh output paths; never overwrite frozen evidence):

```sh
python3 research/navier_stokes_cascade/results/complete_covariance_initial_means/test_audit.py -v
python3 research/navier_stokes_cascade/results/complete_covariance_initial_means/audit.py \
  --verify-record research/navier_stokes_cascade/results/complete_covariance_initial_means/evidence.json \
  --output /tmp/complete_covariance_initial_means_replay.json
```
