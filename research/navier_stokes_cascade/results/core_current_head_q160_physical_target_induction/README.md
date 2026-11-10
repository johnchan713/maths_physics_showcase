# Q160: physical all-finite-cycle target and moment induction

Parent: Q159, `d7337e7714aa866e9d361e1ba69c9afbc4d701dc`.
The live branch, Q159 evidence and CI state, Q149/Q150/Q152, the current
moment audit, current-target cycle bridge, source-edge audit, finite-stage
induction and Q153 were read first. The hosted commit has no workflow run or
combined status. Historical evidence, the finalized profile, roots, physical
primary family, pressure gauge, supports, later tail, diagonal and `C^-20`
width are unchanged.

Q159 correctly isolated the universal weighted target statement, but it
classified the existing current-target recurrence too narrowly. The recurrence
is an analytical all-half-line theorem relative to a current state and its
four actual update maps. Q149, Q150 and Q152 provide its seed and every map
premise on the inherited physical family. This checkpoint performs that
same-object binding and proves the P3 statement for every finite cycle.

The result has constants depending on the finite cycle and fixed jet. It is
not uniform in the cycle number and does not justify an infinite uncut sum.
Passing symbolic tests verifies identities, exact arithmetic and provenance;
it is not the continuum induction below.

## 1. The physical seed is exactly the recurrence seed

Set `kappa=10^-5`. Q152 gives the first physical regular-source order

\[
 B_0={3\over2}-2\kappa .                              \tag{1}
\]

The corresponding current signed target order is

\[
 C_0=B_0+{1\over2}=2-2\kappa .                       \tag{2}
\]

This is not copied from the separate canonical entrance. Q150 gives the
initial physical six-entry covariance/mean target at order `2-3 kappa`, with
the complete axial velocity, compact pressure, moving pressure moment and both
exceptional base pressure tails. Q152 supplies on that same physical family
the exact real principal identity

\[
 H_\beta y_\beta=T_\beta,\qquad
 a_\beta T_\beta=\Sigma_0,\qquad
 \sum_\beta\eta_\beta^2=1,                            \tag{3}
\]

and the complete-curl covariance remainder at order `2-2 kappa`. Applying the
exact real prefix and then the first generic signed current-state difference
therefore gives the physical target `Xi_0 in M_(2-2 kappa)`. The half-order
complete-curl margin is positive. No physical/canonical frame comparison,
partition division or transplanted target is used.

In coefficient form, for every fixed mixed jet `I`,

\[
 |D^I\Xi_0|\le C_{0,I}\epsilon^{C_0}\zeta
        y_a^{-p_{0,I}}y_b^{-q_{0,I}} .                 \tag{4}
\]

The constants include the fixed physical radius, width, profiles, pressure
gauge and finite edge losses, and are uniform over accepted physical bands
and labels for this fixed jet.

## 2. Exact target difference, with complete pressure memory

For `e=2` in the angular row and `e=1` in the axial row, define

\[
 I_ef=R^{-e}\int_{R_-}^Rr^ef(r)\,dr,\qquad
 Q_ef=f-\sigma_e\int R^ef,\qquad
 X_e(E_e)=-I_eQ_e\overline{E_e}.                       \tag{5}
\]

The unit moments of `sigma_e` and the two velocity moments make (5) the
actual current signed target. Since this map is linear in the equation,
every update obeys the identity

\[
 X_e^{new}-X_e^{old}=-I_eQ_e\overline{\Delta E_e}.      \tag{6}
\]

The cancellation must be taken in (6) before terms are estimated. For a
complete six-entry covariance change `D`, radial pressure is reconstructed
from

\[
 \delta g=-D_1^0D_{rr}-D_zD_{rz}+D_{\theta\theta}/R,
 \quad \delta p=I_0(\delta g-\rho\!\int\delta g),       \tag{7}
\]

and the actual target changes are

\[
 \delta X_\theta=-D_{r\theta}-I_2Q_2D_zD_{\theta z},
 \qquad
 \delta X_z=-D_{rz}-I_1Q_1D_z(D_{zz}+\delta p).        \tag{8}
\]

Thus `D_rr,D_rtheta,D_rz,D_thetatheta,D_thetaz,D_zz` all remain. The axial
velocity, `W_thetaz`, `W_zz`, compact pressure, and derivative of the whole
moving product remain. In particular

\[
 Q_eD_zH=D_z(Q_eH)+\epsilon(\partial_Z\sigma_e)
              \int R^eH                                  \tag{9}
\]

is not silently commuted. `Pi1/Pi2` cancel only against the same complete
base equation; neither is removed from current nonlinear products.

Write the signed complete curl as `V=L+R` about the fixed physical principal
wave `W0`. The exact covariance identity is

\[
 \Delta W=\operatorname{Bil}(W_0,L)
 +\operatorname{Bil}(U-W_0,L)+\operatorname{Bil}(U,R)
 +\operatorname{Cov}(V).                              \tag{10}
\]

The selected entries of the first term equal the current `Xi`; the other
four entries and all error terms remain. Substituting (10) into (8) cancels
the old target and gives, for `C>=C0`, the following complete orders:

| Current contribution | Order |
|---|---:|
| old physical-wave error crossed with `L` | `C+1/2-kappa` |
| old wave crossed with full curl remainder | `C+1/2-kappa` |
| signed square, including all curl terms | `2C-1` |
| axial covariance and reconstructed pressure | `C+1` |
| fast old/new means | `C+2-2kappa` |
| fast square | `2C+1-2kappa` |
| slow linear/pressure/viscosity | `C+1` |
| slow old/new means and square | `C+2-kappa`, `2C+1` |

Subtracting `C+1/2-kappa`, every row is affine in `C` with nonnegative
slope and nonnegative value at `C0`; the two signed cross rows tie. Hence

\[
 \Xi\in M_C\quad\Longrightarrow\quad
 \Xi^{new}\in M_{C+1/2-\kappa}                       \tag{11}
\]

on the whole half-line. This is the exact state-difference theorem, rather
than an exponent assigned from compact support or a manufactured fixture.

## 3. The edge weight is propagated, not inferred from compactness

Every physical wave factor has `sqrt(zeta)`, so each complete covariance
product has `zeta`. Mean fields have full `zeta`; since `0<zeta<=1`, their
products retain it. Differentiation introduces only finite inverse powers of
the two edge variables at each fixed jet. The compensation bump is supported
on the fixed interior hull where `zeta` has a positive floor.

The input of each primitive in (5), (7) has exact zero weighted total. The
two-edge oriented primitive estimate therefore retains the full `zeta` at
both ends, with finite powers `y_a^-p y_b^-q`. Without that zero total, a
positive input such as `zeta` has a primitive with nonzero terminal value and
the conclusion is false. This is why ordinary compactness cannot replace the
moment argument.

Consequently (11) holds in the weighted fixed-jet class:

\[
 |D^I\Xi^{new}|\le \widehat C_I
 \epsilon^{C+1/2-\kappa}\zeta
 y_a^{-\widehat p_I}y_b^{-\widehat q_I}.              \tag{12}
\]

## 4. Both complete moment identities persist at every operation

For the actual current means, let `gbar` be the full radial source,
`Pmean=int gbar`, `Qg=int R^2 gbar`, and

\[
 \bar p=I_0(\bar g-\rho P_{mean}),\qquad
 \int R\bar p=-{1\over2}Q_g+c_\rho P_{mean}.           \tag{13}
\]

Direct integration of the two conservative tangential equations, using
compact radial extension and the two velocity moments, gives

\[
 \int R^2\overline{E_{\theta,j}}
 =\epsilon\partial_ZJ_{\theta,j},                     \tag{14}
\]
\[
 \int R\overline{E_{z,j}}
 =\epsilon\partial_Z(J_{z,j}+c_\rho P_{mean,j}).      \tag{15}
\]

Here `Jtheta` contains the full angular-axial flux, while `Jz` contains axial
velocity squares, `W_zz` and `-Qg/2`. Equation (15) differentiates the whole
moving pressure product.

Sourced and signed waves have nonzero angular frequency and do not alter the
two Haar velocity moments. The fast temporal-inverse curl has exact zero Haar
mean in all three components and preserves both moments. The slow five-map
imposes its two zero-moment rows and exact divergence. Therefore (13)--(15)
apply after the sourced, signed, fast and slow updates, and inductively after
every finite physical cycle. Zero-Haar correlations that survive products
remain in the recomputed covariance and pressure; only the linear zero-Haar
response vanishes.

## 5. Close the finite-cycle quantifier

Let

\[
 \delta={1\over2}-3\kappa,qquad
 B_j=B_0+j\delta,\qquad C_j=B_j+{1\over2}.             \tag{16}
\]

Assume (4), (14), (15) at cycle `j`. Sections 2--4 prove them after the whole
current physical cycle at order `C_j+1/2-kappa`. Since

\[
 (C_j+1/2-\kappa)-C_{j+1}=2\kappa>0,                  \tag{17}
\]

and `0<epsilon<=1`, this stronger bound implies the required order
`C_(j+1)`. The seed in section 1 and induction prove, for every finite `j`
and fixed mixed jet `I`,

\[
 |D^I\Xi_j|\le C_{j,I}\epsilon^{C_j}\zeta
 y_a^{-p_{j,I}}y_b^{-q_{j,I}},                         \tag{18}
\]

together with (14)--(15). Constants and edge powers may grow with `j` and
`I`; no stage-uniform bound is needed or claimed.

This discharges P3 on the inherited physical family and supersedes the open
flag in the frozen source-edge record. That historical flag remains unchanged.
Q153 already contains a family-local generated-source/endpoint provider, but
it invoked the finite-stage chain while P3 was still conditional. The next
checkpoint must replay that provider from (18), verify exact residual and tail
allocation, and only then reconsider downstream diagonal, endpoint, force,
comparison and core-exclusion claims. Q160 does not establish an unconditional
complete residual, admissible force, competitor exclusion or blow-up.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q160_physical_target_induction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q160_physical_target_induction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q160_physical_target_induction/evidence.json \
  --output q160_replay.json
```
