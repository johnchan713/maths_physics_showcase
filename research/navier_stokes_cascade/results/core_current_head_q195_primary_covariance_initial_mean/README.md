# Q195: exact primary covariance and complete initial mean join

Parent: Q194, `477a0d22eb94e4332017e98418d0f6874507035c`.
The live branch, Q194, the physical-pulse transfer, the frozen complete
initial covariance/mean evidence, Q183's generated-source frontier and active
CI were read first. Historical evidence, the selected profile, pressure gauge,
roots, supports, accepted scales, frequencies and the `C^-20` activation width
are unchanged. Q190 CI has completed successfully. Q191--Q194 CI had no failed
step and was still in progress when this checkpoint began.

Q194 independently assembled the complete real primary velocity, pressure and
potential on one support domain. This checkpoint performs the next calculation:
it takes the angular covariance of that exact assembled field, retains every
cylindrical component, and joins the resulting tensor to the complete initial
compact-pressure and mean identities. The result is an initial-primary
statement. Later complex waves and generated sources are not inferred from it.

## 1. Exact complete covariance of the cut curl

The potential and pressure cutoffs in Q194 are independent of the physical
angular coordinate. For one nonzero integer angular carrier, applying the full
real cylindrical curl therefore leaves the exact form

\[
 u_p=u^c\cos(k\Phi)+u^s\sin(k\Phi),
 \qquad u^c,u^s\in\mathbb R^3.                         \tag{1}
\]

Both coefficient vectors contain radial, angular and axial velocity. In
particular, `u^s` includes every derivative of the carrier and every cutoff,
as well as the `1/r` cylindrical connections. No leading-polarization
truncation is made. Distinct physical labels and their derivatives have zero
products on the separated closed supports. Since the angular carrier is a
nonzero integer,

\[
 W_{ab}=\langle u_{p,a}u_{p,b}\rangle_\theta
 ={1\over2}(u^c_au^c_b+u^s_au^s_b),
 \qquad a,b\in\{r,\theta,z\}.                         \tag{2}
\]

The mixed sine/cosine mean vanishes exactly. Equation (2) supplies all six
symmetric entries

\[
 W_{rr},\ W_{r\theta},\ W_{rz},\ W_{\theta\theta},
 \ W_{\theta z},\ W_{zz}.                            \tag{3}
\]

It is important that (2) applies after the cut potential has been curled. It
therefore includes the complete cutoff curl, rather than adding a cutoff to a
previous principal covariance formula.

Let

\[
 s_x={1\over2}+{h\over2},\qquad
 s_t=1+{3h\over2}.
\]

Leibniz applied to Q194's complete velocity estimate gives, for every fixed
spatial/time order `(a,b)`,

\[
 |\nabla^a\partial_t^b W_{ij}|
 \le C_{a,b}Q^{-((a+2)s_x+bs_t)}S^{B_{a,b}}
       \zeta\,\delta^{-B_{a,b}}P^2.                  \tag{4}
\]

Thus every covariance entry has the same physical `Q` loss as the Q194 wave
pressure, while carrying the squared shell and Gaussian weights. The wave
pressure itself has angular mean zero but is not zero pointwise. It remains in
the nonzero-angular equations and in later generated sources.

## 2. A pointwise envelope is not a stress identity

Equation (4) is only an upper bound. It cannot imply either selected stress
match or its sign. For the same physical family used in Q194, the unchanged
positive roots and actual Haar columns supply the extra identity

\[
 H_\beta y_\beta=T_\beta,\qquad
 a_\beta T_\beta=\Sigma_0,\qquad
 \sum_\beta\eta_\beta^2=1.                           \tag{5}
\]

Consequently the two principal averaged entries match the selected stress
exactly, and the complete curl remainder obeys

\[
 \overline W_{r\theta,rz}-\Sigma_0
 \in M_{2-2\kappa}.                                  \tag{6}
\]

This is a same-family identity. No closeness argument transfers roots or mean
fields from the canonical family. The correction matters: Q194's field bounds
alone establish (4), not (5)--(6). The root/Haar identity is the structural
input that makes the initial covariance cancellation possible.

## 3. Complete initial compact pressure and scalar targets

Set the old mean velocity to zero and retain the exact slow background. The
full conservative radial source contains `W_rr`, `W_rz` and
`W_thetatheta`. Reconstruct at the current state

\[
 P=\int\overline{g_r}\,dR,\qquad
 p_m=T_0(g_r-\rho P).                                \tag{7}
\]

Radial integration by parts gives the exact initial scalar targets

\[
 P^0=\int{\overline W_{\theta\theta}-\overline W_{rr}\over R}\,dR
      -\epsilon\partial_Z\int\overline W_{rz}\,dR,    \tag{8}
\]

\[
 J_\theta^0=\int R^2\overline W_{\theta z}\,dR,       \tag{9}
\]

\[
 J_z^0=\int R\left[\overline W_{zz}
       -{\overline W_{rr}+\overline W_{\theta\theta}\over2}\right]dR
       +{\epsilon\over2}\partial_Z\int R^2\overline W_{rz}\,dR. \tag{10}
\]

The sixth component `W_rtheta` appears in the angular tangential stress row.
Hence (6) and (8)--(10) use all six entries. Axial velocity cannot be removed:
it contributes to `W_rz`, `W_thetaz` and `W_zz`.

The compact pressure moment is

\[
 \int R\overline{p_m}\,dR
 =-{1\over2}\int R^2\overline{g_r}\,dR+c_\rho P,
 \qquad c_\rho={1\over2}\int R^2\rho\,dR.            \tag{11}
\]

Therefore the axial identity differentiates the whole moving product:

\[
 \partial_Z(c_\rho P)=(\partial_Zc_\rho)P
                         +c_\rho\partial_ZP.          \tag{12}
\]

The first term in (12) is generally nonzero. Both exceptional background
pressure tails `Pi1/Pi2` cancel only through subtraction of the same exact
base equation. The primary wave pressure vanishes only under angular
averaging and is not deleted from the pointwise residual.

## 4. Initial complete mean gain

Haar translation invariance removes the graph part of the averaged radial
derivative. The raw targets (8)--(10) are in `S_1`, while (6) places the two
averaged tangential rows at order `2-2*kappa`.

Apply the actual zero-Haar fast inverse and compact curl. Its three velocity
means vanish exactly, including the flat axial `A1` correction. Recompute
(7)--(11); the surviving changes are quadratic. Apply the fixed slow five-map
to those current targets, then recompute once more. The exact moment identities
remain

\[
 \int R^2\overline E_\theta
 =\epsilon\partial_ZJ_\theta,\qquad
 \int R\overline E_z
 =\epsilon\partial_Z(J_z+c_\rho P).                  \tag{13}
\]

The final three scalar defects belong to `S_2`. The radial mean residual is
of order two, and the complete tangential angular-mean residual has order

\[
 2-3\kappa={199997\over100000}.                       \tag{14}
\]

The loss in (14) is the retained fast viscosity term
`epsilon D_r^2 v_f`. It is not promoted to order two. Flat `A0/A1` and exact
base summands remain separately retained. This closes Q194's stated
complete-primary covariance/initial-mean obligation on the same field and
same support domain.

## 5. Precise remaining boundary

This calculation does not rederive a later complex covariance cycle. Such a
wave can have a principal sine quadrature, so the real-primary square
cancellation in (2) cannot simply be reused. It also does not prove the
generated-source bound required by the finite-cycle induction.

The first fresh direct estimate remains Q183's fixed-cycle gate. For the first
actual nonzero generated harmonic `f_(0,m)` and every fixed physical jet `I`,
one must derive finite `r,s` and a constant on the same domain such that

\[
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-B_0}
 \sup_{\Omega_{0,\ell}}\zeta^{-1/2}y_a^ry_b^s
 \left|D^I\left({f_{0,m}\over P}\right)\right|<\infty, \tag{15}
\]

together with smooth mixed-jet zero extension at both physical pulse
endpoints. The assembled source in (15) must retain source-normal and
moving-normal pressure, `Pi1/Pi2`, vector and mixed viscosity, cutoff
derivatives, cylindrical connections, mean-wave terms, both cross
orientations and the new square.

The exact tests check (2)--(14), rational exponents, provider identities and
immutable provenance. They do not evaluate the enormous selected profile and
are not a continuum proof. No unconditional pulse hierarchy, complete
residual improvement, admissible force, comparison exclusion or
Navier--Stokes blow-up is claimed.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q195_primary_covariance_initial_mean/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q195_primary_covariance_initial_mean/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q195_primary_covariance_initial_mean/evidence.json \
  --output q195-fresh.json
```
