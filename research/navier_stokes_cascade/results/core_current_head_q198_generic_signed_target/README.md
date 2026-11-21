# Q198: corrected generic complex signed target and weighted-input obstruction

Parent: Q197, `17a7a424c30e6e208696b4d98960d835f15c1439`.
The live branch, Q197, its active CI, the frozen generic-signed covariance
provider, the two-edge target budget and Q184 logarithmic-loss correction
were read first. Historical evidence, the original positive roots, pressure
gauge, supports, accepted scales, carriers and unchanged `C^-20` activation
width are preserved. Q197 CI had completed 54 steps without failure and was
running step 55 when this checkpoint began.

Q197 leaves the complete current wave error at order `B=1-kappa`, its
recomputed scalar/pressure changes at order `C=3/2-kappa`, and its regular
nonzero source at order `D=C-kappa`. This checkpoint substitutes those
corrected rates into the existing generic complex signed identity, retains
the full pressure reconstruction, and isolates the first actual-input
inequality still missing. The exponent conclusion is analytical relative to
the named weighted root/operator inputs; diagnostics do not certify those
continuum inputs.

## 1. The current target uses the original roots

Let `Xi=(Xi_theta,Xi_z)` be the normalized current projected target formed
from the complete current Haar tangential rows and moments. For one localized
box keep the original positive amplitudes `a` and original two-column map
`H`. The signed linear update is

\[
 d=H^{-1}(\Xi/\mathrm{scale}),\qquad
 \delta a_j={d_j\over2a_j}.                         \tag{1}
\]

The denominator is the original `a_j`, not an updated amplitude. The two
selected principal cross-covariance entries equal `Xi`. Summing the original
squared slow partition gives the physical target without dividing by a
partition. Both cutoff derivatives and both signs stay inside the columns.
The signed square is not part of the cancellation, and the other four tensor
entries generally remain nonzero.

At the corrected Q197 rates put

\[
 B=1-\kappa,\qquad C=B+\tfrac12=\tfrac32-\kappa,
 \qquad D=C-\kappa=\tfrac32-2\kappa.                \tag{2}
\]

The signed principal field has order `B`; its complete curl remainder has
order `D`; its potential and wave pressure have order `C`. This is a signed
homogeneous correction around the original roots. It is not another sourced
zero-entry solve.

## 2. Complete generic complex covariance

Let `W0` be the original principal wave, `U` the current complete wave,
`E=U-W0`, and `V=L+R` the new complete signed curl. With `Cov` the angular
zero mode and `Bil(x,y)=Cov(x,y)+Cov(y,x)`, exactly

\[
 \Delta W=\operatorname{Cov}(U+V)-\operatorname{Cov}(U)
 =\Lambda+\mathcal E,                               \tag{3}
\]
\[
 \Lambda=\operatorname{Bil}(W_0,L),\qquad
 \mathcal E=\operatorname{Bil}(E,L)
             +\operatorname{Bil}(U,R)+\operatorname{Cov}(V). \tag{4}
\]

Equations (3)--(4) retain both quadratures, all equal harmonics and all six
symmetric entries. Unequal harmonics are not merged before averaging.
Separate closed labels have zero products by support separation.

The orders are

\[
 W_0,U:\tfrac12,\qquad E,L:B,\qquad R:D.
\]

The principal cross has order `C`. Every error term in (4) has order

\[
 2B=\tfrac12+D=2-2\kappa=:C_0.                     \tag{5}
\]

Thus the complete covariance change is still order `C`, while the part left
after the selected principal cancellation is order `C0`. The strict gain is

\[
 C_0-C=\tfrac12-\kappa>0.                           \tag{6}
\]

Equivalently the generic gain condition is

\[
 \min\{B-\tfrac12,\ \tfrac12-\kappa,\ C-1\}
 =\tfrac12-\kappa>0.                                \tag{7}
\]

This uses the Q197 correction that the old wave error is `B`, rather than the
improved regular residual order `D`.

## 3. Full pressure and projected target

Hold the current mean velocities fixed during this signed wave operation.
For the complete six-entry Haar covariance change write `D_ab`. Recompute

\[
 \delta g_r=-D_1D_{rr}-\epsilon\partial_ZD_{rz}+D_{\theta\theta}/R,
 \quad \delta P=\int\delta g_r,
 \quad \delta p_m=I_0(\delta g_r-\rho\delta P),       \tag{8}
\]
\[
 \delta J_\theta=\int R^2D_{\theta z},\qquad
 \delta J_z=\int RD_{zz}-{1\over2}\int R^2\delta g_r. \tag{9}
\]

The unchanged full pressure still contains its base/mean components and both
`Pi1/Pi2` tails. The axial identity differentiates the whole product
`c_rho delta P`, including `(partial_Z c_rho) delta P`.

For the moving moment projection `Q_e`, differentiation does not commute:

\[
 Q_e\epsilon\partial_ZF
 =\epsilon\partial_Z(Q_eF)
  +\epsilon(\partial_Z\sigma_e)\int R^eF.            \tag{10}
\]

Projection is applied after differentiation. Once the principal cross in
(3) is cancelled, the current projected targets are

\[
 X_{\theta,new}=-\mathcal E_{r\theta}
   -I_2Q_2\epsilon\partial_ZD_{\theta z},
\]
\[
 X_{z,new}=-\mathcal E_{rz}
   -I_1Q_1\epsilon\partial_Z(D_{zz}+\delta p_m).      \tag{11}
\]

The axial response in (11) has order `C+1>C0`; hence both targets have order
`C0`. This does not promote the Haar pressure/scalar changes, which remain
order `C`, or the full pointwise pressure change, which remains order `D`.
The radial Haar remainder can therefore still be only order `C`. Fresh fast
and slow mean corrections remain necessary.

## 4. The pressure proof domain must expand

Let `p` be the prior fixed pole, `n` the target output jet, and `d_(n+3)`
the original two-edge target allowance. A sufficient common wave pole is

\[
 q=\max\{3p,\ 6(n+2)+9,\ 6(n+2)+d_{n+3}+6\}.        \tag{12}
\]

The covariance and pressure source requires pole `2q`. For the diagnostic
example `p=1,n=0,d_3=12`, (12) gives `q=30` and source pole `60`; the old
largest declared pole `12` is insufficient. The existing exact shift and
oriented primitive recipes may be reevaluated at these finite domains, but
that changes constants rather than the physical field, width or band.

## 5. Precise remaining inequality

Q197 and the frozen provider establish the algebra and rate closure
conditionally. They do not supply an independent actual-field certificate
for every expanded request. For each mixed jet `I`, the missing same-family
bound is

\[
 |D^I\Xi|\le K_I\epsilon^C\zeta
 y_a^{-p_I}y_b^{-q_I},                               \tag{13}
\]

together with the separate original-principal, current-error, reciprocal-root,
common-chart, support and compact-operator bounds through the `2q` pressure
pole. Ordinary unweighted compact bounds do not imply (13): the bounded
function `1` has quotient `1/[x^2(1-x)^2]`, which diverges at both flat edges.
This is an obstruction to an inference, not a counterexample to the intended
canonical field.

The fixed polynomial `S` losses can again be absorbed once, after a finite
request is collected. The parallel pure ledger is

\[
 \widetilde B=1-2\kappa,\quad
 \widetilde C=\tfrac32-2\kappa,\quad
 \widetilde D=\tfrac32-3\kappa,\quad
 \widetilde C_0=2-3\kappa.                           \tag{14}
\]

The physical zero-jet regular residual remains
`-15025001/10000000<0`. Therefore (6) is neither absolute physical decay nor
a complete residual theorem.

## 6. Scope and next obligation

This checkpoint closes the corrected exponent algebra, exact complex tensor
identity, complete pressure routing, projection commutator and finite domain
expansion. It does not certify (13), construct an unconditional signed pulse,
perform the post-signed fast/slow relaxation, or improve the complete
nonzero residual beyond `D`.

Next validate (13) and its separate original-root/chart/support/operator
inputs at every finite request on the expanded domain. Then apply the actual
post-signed fast map, recompute pressure and all defects, apply the current
slow five-map, and recompute them again. Complete residual, force,
comparison, infinite diagonal and blow-up obligations remain open.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q198_generic_signed_target/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q198_generic_signed_target/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q198_generic_signed_target/evidence.json \
  --output q198-fresh.json
```
