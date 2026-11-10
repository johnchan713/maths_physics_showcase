# Q152: inherited physical stress identity and feedback join

Parent: Q151, `387edc730980a28fe72d61a5d33d2d50402c56a6`.
The live branch, Q151 evidence and CI, the inherited actual-band pulse,
physical-pulse foundation, complete initial covariance, same-weight primary
covariance, Haar-stress criterion, scalar-gain and fast-row records were read
first. Historical evidence, the finalized profile, roots, supports, pressure
gauge, accepted scales and `C^-20` width are unchanged.

Q151 correctly rejected any transfer of the canonical stress identity to the
inherited physical primary by `O(1/L)` frame or root closeness. It overlooked
a different retained route: the inherited physical family already has its own
same-family stress columns, positive roots, physical square normalization and
squared slow partition. This checkpoint records that correction, closes the
physical principal-stress gate, and joins the already derived scalar and
fast-row estimates to that same family. It does not identify the physical and
canonical entrances.

## 1. The physical family has its own exact principal assembly

For each active physical slow box `beta`, let `H_beta` be the two selected
Haar columns computed from that family's cut homogeneous pulses. The
transverse and temporal cutoffs, angular factor `1/2`, Haar Jacobian and both
signs are already inside `H_beta`; the slow partition `eta_beta` is outside.
The original positive roots satisfy

\[
 H_\beta y_\beta=T_\beta,
 \qquad T_\beta=(Q_\beta/q)^{A+1/2}T_0.                 \tag{1}
\]

This is the inherited actual-band root equation, rather than a root copied
from the separately rebuilt canonical entrance. With

\[
 \epsilon_\beta=Q_\beta^h,\qquad
 a_\beta=Q_\beta^{-2A}\epsilon_\beta,
 \qquad \Sigma_0=q^{-A-1/2}T_0,
 \quad A=\tfrac12+h,
\]

the physical velocity-square normalization gives

\[
 a_\beta T_\beta
 =Q_\beta^{-2A+h+A+1/2}q^{-A-1/2}T_0
 =\Sigma_0.                                             \tag{2}
\]

The exponent of `Q_beta` is `h-A+1/2=0`. The prescribed physical cover has

\[
 \sum_\beta\eta_\beta^2=1.                              \tag{3}
\]

Therefore the selected principal covariance of the complete localized
**physical** primary obeys the tensor identity

\[
 \overline{(W_{0,\mathrm{principal}})}_{r e}
 =\sum_\beta\eta_\beta^2 a_\beta H_\beta y_\beta
 =\Sigma_{0,e},\qquad e\in\{\theta,z\}.                \tag{4}
\]

Equation (4) differentiates only after the full assembly. Derivatives of the
partition, target, roots and columns cancel through the exact identity; no
partition is divided out and neither sign is counted twice. Recutting a
column after solving (1), omitting a box, or using the velocity scale only
once would leave a nonzero root, cover or normalization defect.

## 2. The complete curl gives an independent quantitative margin

The actual real primary has exact angular covariance

\[
 W_0={\epsilon\over2}
   \left(A\otimes A+k^{-2}B_{\rm curl}\otimes B_{\rm curl}\right). \tag{5}
\]

Real sine/cosine quadrature removes the mixed term, while every entry of the
curl square remains. The inherited all-band estimate gives, for every fixed
coefficient jet,

\[
 \|W_0-W_{0,\mathrm{principal}}\|_n
 \le \epsilon^{2-2\kappa}Q_n(S)\zeta\delta^{-p_n}.      \tag{6}
\]

Thus Q151's principal discrepancy `L_e` is exactly zero by (4). Even if the
complete selected covariance is used instead, division by `epsilon` leaves
order `1-2kappa`, whereas the required gate has order

\[
 D-1=\tfrac12-2\kappa.
\]

The complete-curl estimate has a strict margin of `1/2`. Hence

\[
 \|\epsilon^{-1}L_e\|_n
 \le \epsilon^{D-1}Q_n(S)\zeta\delta^{-p}               \tag{7}
\]

holds on the inherited physical family without any polynomial comparison to
the canonical family. Q151's asymptotic objection to `O(1/L)` transfer
remains correct and is simply unnecessary here.

## 3. Same-family scalar, row and feedback consequences

Q149 supplies the same inherited physical primary on every used dyadic band
with uniform fixed physical derivatives, complete curl, pressure and closed
supports. Q150 supplies its initial six-entry covariance/mean join. Equations
(1)--(7) now discharge the principal selected-stress input isolated by the
Haar-row criterion for this family.

Keep

\[
 B=1-\kappa,\qquad C=B+\tfrac12,
 \qquad D=C-\kappa=\tfrac32-2\kappa.                    \tag{8}
\]

The existing exact scalar subtraction gives postfast
`P,J_theta,J_z in S_C`, retaining axial velocity, `W_theta z`, `W_zz`,
`W_rz`, both background pressure tails and the derivative of the whole
moving pressure moment. The exact fast-row subtraction retains every
old/new mean cross, new square, cylindrical term, viscosity term, compact
pressure difference and the nonzero `A1` remainder. The complete sourced
covariance contributes at `C-kappa=D`. With (7), the physical postwave
zero-Haar rows and fresh fast increment have order `D`; the fresh slow
increment has order `C`. Both fresh mean/current-wave orientations therefore
have order

\[
 \min(D,C)=D.                                           \tag{9}
\]

Consequently the next regular nonzero source has the conditional order `D`
and pole `9p` on the same physical family. The exact additive tails remain
outside that source and every inverse. At `kappa=10^-5` the exposed compact
request still reaches very large but finite jets; this checkpoint does not
replace their raw profile/source bounds by a generic operator norm.

## 4. The next precise obligation

The stress-family gate is no longer the frontier. The physical sourced
inverse still requires every generated regular harmonic to belong to its
full-path Gaussian class and to extend smoothly by zero at both pulse
endpoints. For every fixed output jet and nonzero generated harmonic, the
next join must establish a bound of the form

\[
 \sup_v P(v)^{-1}|D^I f_m(v)|
 \le \epsilon^D Q_{I,m}(S)\sqrt\zeta\,\delta^{-p_I},
 \qquad f_m\text{ has smooth endpoint zero extension}, \tag{10}
\]

at all finite high jets requested by `A0/A1`. Separately, every excluded
Gaussian, spatial, scale and compact-pressure tail must satisfy, after the
actual physical derivative losses,

\[
 |D_{\rm phys}^I T_{\rm add}|\le C_{I,N}q^N
 \quad\hbox{for every fixed }I,N.                       \tag{11}
\]

Only after (10)--(11) are bound to this same physical family can the existing
finite-stage/complete-potential residual diagonal be invoked. The current
checkpoint does not prove the full angular-zero residual, all generated
source supports, a complete residual theorem, admissible forcing, or blow-up.
Passing tests checks exact identities, exponents and frozen provenance; it
is not an independent continuum proof.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q152_physical_stress_feedback_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q152_physical_stress_feedback_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q152_physical_stress_feedback_join/evidence.json \
  --output q152_replay.json
```
