# Regular harmonic-source foundation: a sharper linear bound and exact tail split

Parent: `10b5bdff358cd6c54474be533dec1618f829d940`.
The live [mean-map foundation](../compact_mean_foundation/README.md) and
[primary-pulse foundation](../physical_pulse_foundation/README.md) now supply
the inputs used below, relative to their retained exact-core and finite
background/pulse existence and gap lemmas. This checkpoint independently
audits the regular source operators, not the manuscript's final theorem.
All historical source/evidence bytes, h, N, positive lambda, `C^-20`,
supports, discrete labels and primary roots stay unchanged.

Result: the COMPLETE regular linear remainder of a pulse of order alpha
has the sufficient stronger order alpha+1/2-kappa. The previous bound
alpha+1/2-3kappa remains valid; its extra two losses are unnecessary under
the supplied ALL-fixed-coefficient-jet estimates. Complete same-label curl
transport has order alpha+alpha'-kappa, while mean--wave transport keeps
the genuine half-power loss. These results supply the regular nonzero
source of the actual primary-plus-initial-mean state at order 1-kappa,
including full endpoint extension. They do NOT independently close all
later complex covariance/mean states or the complete residual iteration.

Reference under audit: supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(7.38)--(7.40), (9.1)--(9.4), SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The estimates below are analytical derivations. Exact rational manufactured
jets check operator identities, not continuum bounds or actual giant fields.

## 1. Fixed chart, classes and the complete linear operator

Write epsilon=Q^h, S=ell^2, k=ceil(epsilon^-1/2), kappa=10^-5. Fixed
nonzero harmonic multipliers m may change constants but not epsilon powers.
Use the frozen normalized operators
Dr=partial_R+M dr R^(dr-1)L, Dz=epsilon partial_Z,
tstar=-epsilon partial_T+c_i Nt. All ordinary coefficient derivatives
preserve a class W_alpha, with polynomial S costs, weight sqrt(zeta)P and
fixed inverse edge-distance powers. A displayed Dr costs at most kappa;
Dz and slow time gain one. Requesting another COEFFICIENT jet does not
apply an additional physical Dr. Background b is O(epsilon), V=RF and
G are O(1), and their fixed coefficient jets obey the same orders.

Let omega=km, n=(Dr Phi,p/R,Dz Phi), n.t=0, with n and all fixed
coefficient jets polynomially bounded and |n| bounded below. The reviewed
phase defect Eik=(tstar+b Dr+F partial_theta+G Dz)Phi is epsilon times
polynomial S. The projected pulse solve retains both moving-normal and
source-normal pressure terms; t in W_alpha gives pi in W_(alpha+1/2).

For ANY theta-independent amplitude a and pressure pi, the full cylindrical
linearization about B=(b,V,G), after removing exp(i omega Phi), equals

\[
\mathcal L(a,\pi)=c_iN_ta+Ka+\epsilon\omega^2|n|^2a
                          +i\omega n\pi+\mathcal E(a,\pi),
\quad K=\begin{pmatrix}0&-2F&0\\2F+RF_R&0&0\\G_R&0&0\end{pmatrix}, \tag{1}
\]
\[
\begin{split}
\mathcal E(a,\pi)={}&(-\epsilon\partial_T+bD_r+GD_z)a+i\omega Eik\,a
 +a_r b_R e_r+a_zD_zB+(b/R)a_\theta e_\theta+(D_r\pi,0,D_z\pi)\\
&-\epsilon\{(D_r^2+R^{-1}D_r+D_z^2)a-R^{-2}(a_r,a_\theta,0)\\
&\hspace{13mm}+2i\omega(n_rD_r+n_zD_z)a
 +i\omega(D_rn_r+n_r/R+D_zn_z)a
 +(2i\omega n_\theta/R)(-a_\theta,a_r,0)\}. \tag{2}
\end{split}
\]

This follows directly by differentiating each scalar harmonic and using
the cylindrical vector Laplacian and transport connections. There is no
discarded radial advection, axial pressure or angular viscous connection.
Repeated Dr differentiates its radial multiplier. The source's pressure
normal term includes n.f even when the source is not transverse.

The exact potential C=i(n cross t)/(omega |n|^2) yields COMPLETE amplitude
a=t+r, r=curl_* C in W_(alpha+1/2-kappa). The coefficient curl includes
(Dr+1/R)C_theta, not just Dr C_theta. Its all-fixed-coefficient-jet bound
comes from one displayed spatial derivative and 1/k; it is NOT a statement
that the complete amplitude remains transverse.

## 2. Sharper sufficient linear gain, with each loss charged once

For t in W_alpha, the gains in (2) are:

| Retained group | Gain above alpha |
|---|---:|
| slow time, b Dr, G Dz | 1-kappa |
| i omega Eik | 1/2 |
| base derivatives and b/R | 1 |
| pressure-amplitude gradient | 1/2-kappa |
| amplitude vector viscosity | 1-2kappa |
| mixed viscosity and phase divergence | 1/2-kappa |
| angular vector-Laplacian connection | 1/2 |

For the principal velocity operator acting on r, c_i Nt, K and
epsilon omega^2 |n|^2 ALL preserve its order alpha+1/2-kappa. In
particular Nt is a coefficient derivative: its fixed jets have already
been bounded by polynomial S in the differentiated pulse solve, including
source derivatives. It is not a second spatial graph derivative. Since
c_i<=1/S, it introduces no negative epsilon power. The remaining terms
acting on r (with pressure zero) have order at least alpha+1-2kappa,
which is stronger when kappa<1/2. Thus

\[
\mathcal L(t+r,\pi)+f\in W_{\alpha+1/2-\kappa}
\quad\hbox{before temporal cutoff},\qquad
\delta_{lin}=1/2-\kappa>0. \tag{3}
\]

The case f=0 includes homogeneous and fixed signed pulses. The case f!=0
uses the supplied source's ALL fixed coefficient derivatives along the
whole lifted path, not just its value. No parameter or frequency is retuned.
The exact term ledger checks this sufficient bound; the two limiting
groups are pressure gradient and mixed phase/amplitude viscosity, alongside
the principal operator on the complete curl remainder. We do not claim
that every pulse saturates the bound, or that this improves the entire
finite-cycle theorem without checking its other substeps.

## 3. Nonlinear sources: complete divergence, not principal transversality

For a COMPLETE curl harmonic a exp(i km Phi), exactly

\[
div_*a+i km(n\cdot a)=0,\qquad
div_*a=D_ra_r+a_r/R+D_za_z. \tag{4}
\]

Hence n.a=i div_*a/(km) belongs to W_(alpha+1/2-kappa); it generally
is NOT zero. For another same-label harmonic a' exp(i km' Phi), the
factored transport is

\[
a_rD_ra'+a_zD_za'+R^{-1}(-a_\theta a'_\theta,a_\theta a'_r,0)
                       +i km'(a\cdot n)a'. \tag{5}
\]

The apparent half-power loss in the last term is cancelled by (4),
leaving W_(alpha+alpha'-kappa) for a nonzero output and the corresponding
M class for m+m'=0. Negative harmonics cause no exception: m'/m is a
fixed finite multiplier and division by the OUTPUT m+m' is never used.
Both old/new orders and the new square are retained. Product weights are
zeta P^2; since P<=1 and zeta<=1 these dominate the required wave weight.
All fixed jets follow by Leibniz' rule with finite polynomial/edge losses.

For a mean (beta,v,gamma) with orders (mu+1,mu,mu), the symmetric
mean--wave interaction is in W_(alpha+mu-1/2). Tangential mean advection
of the phase really costs k; no wave incompressibility identity cancels
it. Other groups have sufficient orders alpha+mu-kappa,
alpha+mu+1-kappa and alpha+mu+1/2. For example a constant axial mean
advecting cos(kz)e_x has nonzero derivative of size k, so a stronger
generic mean--wave bound cannot be inferred from amplitude derivatives.

For DIFFERENT overlapping phases, (4) does not bound a.n'. The exact
plane-wave pair cos(kx)e_y, cos(ky)e_x has symmetric transport with mean
squared norm k^2/2. Actual different-label derivative products instead
vanish by the parent's separated CLOSED enlarged supports. A distinct
phase is never silently treated as a same-label harmonic.

## 4. The supported-source grammar and what cannot enter the inverse

At each finite expression level:

* Linear regular expressions are the cutoff psi times (2) and the
  principal operator on the uncut curl remainder. The potential AND wave
  pressure are cut before taking the actual curl and global zero extension.
* Quadratic and mean--wave expressions always use COMPLETE cut velocities,
  including all cutoff curl pieces and actual compact A0/A1 velocities.
  Every nonzero term retains a wave factor. Even if a mean inverse fills
  the torus, that product retains the wave's closed slow/transverse/temporal
  support, Gaussian P and smooth zero extension at the full pulse endpoints.
* Same-label harmonics add, zero outputs enter the full mean system, and
  separated labels have zero derivative products. No division by a cutoff,
  partition or updated amplitude is used to assert support. A signed
  target uses the ORIGINAL positive primary roots and its retained half
  edge-weight calculus. That map's later target bounds still need their
  own covariance/mean induction, not this source lemma alone.
* A sourced zero-entry inverse fixes slow/transverse parameters along its
  common-torus path. Uniqueness preserves support in those parameters.
  Its UNCUT exit need not be zero. Cutting its potential and pressure
  restores all endpoint jets. A Gaussian envelope by itself is not support.

The full physical cutoff commutator, with u=curl A and
d_psi=grad psi cross A, is exactly

\[
T_\psi=(\psi_t+B\cdot\nabla\psi)u
 -2(\nabla u)\nabla\psi-(\Delta\psi)u
 +\pi\nabla\psi+\mathcal L_Bd_\psi. \tag{6}
\]

After the prescribed principal source cancellation the exact split is
psi times the uncut linear remainder PLUS (1-psi)f+T_psi. The last
two are designated additive tails, not forwarded sources. Their physical
jets are flat on the SAME starting domain by the parent's Gaussian square
completion; constants depend on derivative and requested power. Their
actual cutoff velocities still enter every later nonlinear product.
The base flat residual is a separate additive term and need not have a
wave envelope, shell support or shell edge weight.

Why this split is essential, rather than optional bookkeeping: for any
nonzero homogeneous pulse (t_h,pi_h), take a genuine smooth step chi
turning on before the core and constant on endpoint collars. The source
f=chi' t_h is physically flat and zero on both endpoint collars, yet
its exact zero-entry response is (-chi t_h,-chi pi_h), nonflat in the
core. Normal pressure and the moving constraint remain satisfied.
Physical flatness of an input therefore cannot certify a flat inverse
response. The older tail obstruction remains in force.

For a finite harmonic set H, sourced corrections retain H, signed
corrections add +/-1, and quadratic residuals have |m|<=2 max(max|H|,1).
All constants may grow with the finite expression/stage. This is a closed
source-operator grammar GIVEN its coefficient inputs, not a proof that
every actual later covariance/mean target has those inputs.

## 5. Join the actual initial state, and stop at the next obligation

Keep the complete actual primary W_(1/2) and initial fast/slow mean state
of [initial means](../complete_covariance_initial_means/README.md), whose
compact-map inputs are now supplied by the parent. Its tangential means
are M_(1-kappa), radial mean M_(2-kappa), with exact velocity moments.
The slow means are not set to zero. Its full pressure uses all six wave
covariances, the moving c_rho P moment and the SAME background pressure,
including both Pi1/Pi2 tails. Pure mean/base pressure is angular zero;
this does not remove it from the mean equations or truncate the base.

The regular nonzero residual has exactly these sufficient orders:

\[
\begin{array}{c|c}
\hbox{complete primary linear remainder}&1-\kappa\\
\hbox{complete primary quadratic terms}&1-\kappa\\
\hbox{actual mean--primary terms}&1-\kappa
\end{array} \tag{7}
\]

Thus the first REGULAR supported inverse input is W_(1-kappa), rather
than the older sufficient W_(1-3kappa), with all fixed coefficient jets,
the original half-edge/Gaussian weight, all nonzero multipliers +/-1,+/-2,
and full smooth endpoint extension on every accepted scale. Retained flat
tails are excluded from that source but remain in the exact residual.
This closes a specific initial source foundation, not a complete residual
gain, the construction of all subsequent sourced waves, or blow-up.

Next join the actual first sourced correction, full complex covariance
change and recomputed mean pressure/defects to this source estimate; then
independently close the later nonlinear induction and ALL residual groups.
Do not substitute the sharper exponent into historical frozen records or
claim a stronger full-cycle rate without checking its other limiting rows.
No unconditional all-pulse-family or full residual improvement is asserted.

## 6. Exact checks and reproduction

`source_checks.py` compares full differentiated cylindrical harmonics with
(1)--(2) using rational jets, nonconstant radial graph coefficients,
nonzero axial/background/pressure terms and positive/negative multipliers.
Separate complete-curl complex coefficients check (4)--(5), including
zero and negative outputs. Omission controls detect the angular viscous
connection, mixed viscosity, pressure gradient and differentiated radial
multiplier. Parent physical Cartesian/cylindrical checks are replayed.
Finite polynomial endpoint fixtures are only C7 after zero extension;
the analytical source proof requires the actual C-infinity cutoffs.
Neither those fixtures nor exact exponent arithmetic are continuum proofs.

```sh
python3 research/navier_stokes_cascade/results/regular_source_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/regular_source_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/regular_source_foundation/evidence.json \
  --output /tmp/regular_source_foundation_replay.json
```

Use fresh output paths; preserve historical evidence.
