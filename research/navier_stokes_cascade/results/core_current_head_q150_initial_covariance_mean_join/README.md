# Q150: initial complete covariance/mean join and later feedback gap

Parent: Q149, `dd7aa09ebc33f588e165b3b3b61bc62a87afe3a2`.
The live branch, Q149 evidence and current CI, the complete initial covariance
and mean audit, covariance and dynamic-mean rechecks, and the quantitative
multiharmonic residual budget were read first. Historical evidence is
unchanged. The finalized profile, physical primary family, pressure gauge,
supports, carriers, accepted bands and `C^-20` width are unchanged.

Q149 correctly showed that the primary covariance remainder
`Q^(-1-2h*kappa)` cannot be bounded by absolute decay. Its stated frontier
was nevertheless one stage too early: the repository already contains an
audited exact cancellation for the **initial** complete covariance/mean
update on the same physical primary family. This checkpoint binds that
result to Q149 and identifies the first open inequality in the later
signed/dynamic cycle.

## 1. Same-primary identity and complete pressure memory

Q149 obtains its complete real primary from `physical_pulse_transfer` on
every sufficiently late used dyadic band. The initial covariance audit
replays that same provider, with the same `Md=64`, positive `lambda,h,N`,
global geometry and width. Thus its fixed primary `w` is the physical family
joined in Q149; no canonical-family identification is used.

Write

\[
 W_{ab}=\langle w_aw_b\rangle_\theta,
 \qquad a,b\in\{r,\theta,z\}.
 \tag{1}
\]

All six symmetric entries are retained, including `W_rz`, `W_theta z`,
`W_zz`, radial velocity and the complete curl remainder. The real-primary
cancellation gives the two matched averaged stresses at order
`2-2*kappa`; it is not assigned to later complex waves.

The full background pressure remains

\[
 q^{2A}p_B=\Pi_0+\chi_1q^{2h}
 \left[K_1-{\mu^2\over4X}\right]+\chi_2q^{4h}K_2.
 \tag{2}
\]

Both tails cancel only when subtracting the same exact base equation.
The wave pressure is zero only in angular mean. At every updated mean
velocity the compact pressure is reconstructed from

\[
 P=\int\overline{g_r}\,dR,
 \qquad p_m=T_0(g_r-\rho P),
 \tag{3}
\]

and the axial moment contains the whole moving product

\[
 \int R\overline{p_m}\,dR
 =-\frac12\int R^2\overline{g_r}\,dR+c_\rho P.
 \tag{4}
\]

Consequently `partial_Z(c_rho P)` includes the generally nonzero
`(partial_Z c_rho)P`. Axial velocity, `W_theta z`, `W_zz`, `W_rz`, the
compact pressure and both exceptional background tails are therefore part
of the joined identities. No leading-polarization truncation is made.

## 2. Exact initial cancellation and recomputation

The initial old mean is zero. Haar translation invariance removes the graph
loss from the averaged radial derivative, so the raw scalar defects
`P^0,Jtheta^0,Jz^0` lie in `S_1`. The zero-Haar temporal inverse and compact
curl give a fast update

\[
 \overline{\beta_f}=\overline{v_f}=\overline{\gamma_f}=0
 \tag{5}
\]

exactly, including the flat compact-curl correction. Fast linear changes
therefore vanish only after Haar averaging; the quadratic changes remain.
Reconstructing (3)--(4), applying the fixed slow five-map, and reconstructing
again gives

\[
 P^{new},\ J_\theta^{new},\ J_z^{new}\in S_2.
 \tag{6}
\]

These three quantities are generally nonzero. The complete angular-mean
annular residual obeys

\[
 \overline E_{\theta,z}^{new}\in M_{2-3\kappa}
 \quad\hbox{modulo the separately retained flat errors},
 \tag{7}
\]

while the radial part is in `M_2`. The limiting loss in (7) is the retained
fast viscosity term `epsilon D_r^2 v_f`; it prevents promotion to order two.
The exact base flat residual is additive and need not carry the shell weight.

With `epsilon=Q^h`, (7) has the positive shell exponent

\[
 Q^{h(2-3\kappa)}.
 \tag{8}
\]

Thus Q149's negative absolute covariance power is canceled at this initial
stage by structure, not dominated term by term. For every fixed coefficient
derivative order the inherited product and inverse lemmas give the same
orders uniformly on Q149's accepted later tail. This completes the initial
complete covariance/mean join. It does not estimate the full nonzero-angular
residual or add a signed covariance wave.

## 3. Why the initial gain cannot be reused later

Set

\[
 B=1-\kappa,\qquad C=B+\frac12,\qquad
 D=C-\kappa=\frac32-2\kappa.
 \tag{9}
\]

The initial residual order is stronger than `D`:

\[
 (2-3\kappa)-D=\frac12-\kappa>0.
 \tag{10}
\]

Equation (10) is an initial-state fact. Later signed waves change all six
covariance entries and meet nonzero old Haar means. The exact reduction used
in (5)--(6) then acquires old-mean/new-increment crosses. The dynamic recheck
shows these crosses have their natural correlation rate; same-phase complete
waves likewise retain the old/new covariance cross. Neither support,
divergence nor the two velocity moments remove those terms.

The quantitative current residual ledger keeps all signed first and second
input harmonics, full uncut source-normal pressure, complete cut velocities,
both orientations of fresh mean/wave feedback, the zero-mode pressure, and
all three scalar moments. Every group reaches order `D` except the fresh
mean/current-wave feedback. From the proved coarse mean size it has only
order `B`, missing the desired order by

\[
 D-B=\frac12-\kappa.
 \tag{11}
\]

The first genuinely open later-cycle estimate is therefore, for every fixed
`n` on the original accepted band,

\[
 \left\|\mathbb P_\perp
 \left[N(d,w_{cur})+N(w_{cur},d)\right]\right\|_
 {n,\,\mathrm{half\ weight},\,9p}
 \le \epsilon^D Q_n(S).
 \tag{12}
\]

One sufficient route is to prove the actual post-wave zero-Haar fast row
inputs at order `D` and the fresh postfast `P,Jtheta,Jz` at order `C`, with
the extra `Jtheta` derivative needed by the radial slow potential. This is
an unproved sufficient route, not an added hypothesis. A correlated proof
of (12) could also suffice.

The inherited compact class witness shows that coarse size, support,
divergence and moment conditions alone cannot prove (12): a mean of order
`B` advecting a half-frequency wave produces exactly order `B`. That witness
is not an output of the prescribed actual fast/slow maps, so it does not
refute those maps or the conditional finite-cycle construction.

## 4. Scope and reproduction

Established analytically by composition of frozen providers: Q149's uniform
physical primary family enters the complete six-entry initial covariance and
compact mean identities; both pressure tails and the moving pressure moment
remain; the recomputed initial scalar defects satisfy (6), and the
angular-mean residual satisfies (7). The rational tests check identities,
exponents, provider contracts and immutable provenance. They are not an
evaluation of the enormous profile or a continuum proof by themselves.

Still open: (12) for the actual later increments, or its sharper sufficient
row/scalar estimates; uniform later dynamic norm comparisons; every sourced
inverse envelope and endpoint support; the full nonzero-angular residual;
viscosity, cutoff and all interaction terms through every finite cycle; the
infinite field, forcing, energy and blow-up obligations.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q150_initial_covariance_mean_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q150_initial_covariance_mean_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q150_initial_covariance_mean_join/evidence.json \
  --output /tmp/q150_replay.json
```
