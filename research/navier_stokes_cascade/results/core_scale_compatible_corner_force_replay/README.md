# Scale-compatible closed-corner and compact-force replay

Parent: `9a7ef975a3c5dab72cb01aaad64c01c10f157d47`.
Status: `scale-compatible-corner-flatness-and-force-replayed`.

The parent corrects the false requirement that support distances, root radii
and flat edge weights have one positive physical minimum.  This checkpoint
replays that corrected certificate through the closed-chart coefficient,
logarithmic tail, Cartesian corner, localization and force arguments.  The
result is a justified analytical implication relative to the certified
provider chain.  None of these steps uses a positive physical support floor
or a positive pointwise edge-weight floor.

The replay does **not** promote the local leading-core growth calculation to
a lower bound for the same complete diagonal field.  That distinct estimate
is isolated in section 5.  Historical evidence is unchanged.

## 1. Correct primitive inputs suffice for the fixed-stage constants

On the closed terminal chart

\[
 T+\eta^2=1,\qquad |\eta|\le1,\qquad
 L=1-2h\eta^2\ge49/50,                                  \tag{1}
\]

the corrected primitive certificate consists of uniform background gaps,
the normalized support separation `3/8012`, dimensionless annular and repair
contraction slack, fixed normalized support hulls, flat zero extensions and
finite weighted seminorms at each fixed order.  Finite composition therefore
gives, for every fixed stage `J`, increment `j`, jet `m` and compact normalized
radial interval `K`,

\[
 C^{\rm stage}_{J,m,K}
 =\sup q^{-(hB_J-K_m)}|\mathcal R(U^{[J]})|_m<\infty,     \tag{2}
\]

\[
 C^{\rm rep}_{j,m,K}
 =\sup {q^{-(g_j-L_m)}|Z_j|_m\over
               (1+|\log q|)^{P_{j,m}}}<\infty.           \tag{3}
\]

The suprema include `|eta|=1`.  Fixed normalized supports give compact-`X`
uniformity.  Physical support distances may shrink like
`3s/8012 -> 0`; only disjointness at each positive scale is used.  The flat
edge weight may vanish; its upper weighted seminorm is exactly the domain
needed by the compact inverse and zero-extension operations.

## 2. Half-gain and one fixed comparison stage

Put `g_j=h delta j>0`.  For `s=-log q>=1`, (3) and
`log(1+s)<=sqrt(s)` give

\[
 (1+s)^P e^{-g_js/2}\le \exp(P^2/g_j).                  \tag{4}
\]

The construction chooses the cutoff scale after this finite constant is
known, so on the whole selected interval

\[
 |(V_{j,a_j},P_{j,a_j})|_m
 \le2^{-j}q^{g_j/2-L_m}.                                \tag{5}
\]

Fix a derivative order `m` and target power `P`.  Choose one finite `J`, held
fixed while `q` tends to zero, satisfying

\[
 hB_J-K_{m+2}\ge P+1,\qquad
 {h\delta(J+1)\over2}-L_{m+2}\ge P+1.                  \tag{6}
\]

Both left sides have positive slope in `J`.  Equations (2), (5), (6), and
the exact residual difference give

\[
 \sup_{X\in K,\ |\eta|\le1}|\mathcal R(u^*,p^*)|_m
 \le C_{m,P,K}q^{P+1}.                                  \tag{7}
\]

Here the residual contains time differentiation, unit viscosity, the full
fixed-gauge pressure, both nonlinear cross terms and the correction square.
The complete potentials retain axial velocity, all six covariance entries,
moving pressure moments and both exceptional pressure memories.  The stage
is never replaced by an uncontrolled `J(q)`.

## 3. Cartesian corner and the extra force derivative

Let `rho=tau+|z|^(1/D)`.  The exact implicit-coordinate estimate is

\[
 q\le4\rho.                                               \tag{8}
\]

Thus (7) gives `|R|_m <= 4^(P+1) C rho^(P+1)`.  Since `P` is arbitrary,
every fixed residual jet is flat at the Cartesian terminal corner.  Asking
for the next time derivative is simply another finite request in (2)--(7),
so it supplies the Cauchy compatibility needed for the force jets.

For an order-`m` localized force jet, the exact twelve-term localization
formula requires the potential through order `m+3`, direct swirl through
`m+2`, pressure through `m+1`, and their next time derivatives.  The
punctured endpoint checkpoint supplies those fixed orders.  Corner flatness,
the punctured compact cover and the exact heat exterior therefore give one
smooth compact terminal force jet.  The already proved Borel cutoff diagonal
extends it for `t>=1`; the preterminal energy identity gives

\[
 \|u(t)\|_2^2+2\int_0^t\|\nabla u(s)\|_2^2\,ds
 \le\left(\int_0^t\|f(s)\|_2\,ds\right)^2,
 \qquad t<1.                                             \tag{9}
\]

This closes the earlier corner obstruction and the compact-force/energy
implication relative to the same certified fixed-gauge provider chain.  It
does not construct a future velocity or pressure.

## 4. What was corrected in the replay

Replacing the heterogeneous scalar margin by the scale-compatible
certificate changes no exponent in (2)--(7).  The proof only needs:

* positive background and cone gaps for reciprocals;
* exact normalized disjointness for products of distinct labels;
* contraction constants below one for the two root maps; and
* upper weighted edge seminorms for compact maps and smooth zero extension.

It never takes the reciprocal of a physical support distance, an absolute
root radius, or the flat weight.  Consequently the parent correction is
compatible with the complete corner and force chain; a fictitious scalar
minimum is neither required nor silently reintroduced.

## 5. Precise remaining lower-bound obligation

The local growth checkpoint proves positive leading-core shear/growth and
the first activation collar.  Its stated scope does not prove the following
same-field estimate for the complete-potential diagonal:

\[
 \exists X_{\rm in},e_0,\tau_*>0\quad
 \inf_{0<\tau<\tau_*}
 \tau^A\left|u^*_{\theta}
   (\sqrt{2X_{\rm in}\tau},0,1-\tau)\right|\ge e_0.      \tag{10}
\]

A sufficient proof must identify the physical leading coefficient on that
path and establish, on the same selected cutoff sequence, either

\[
 V_{j,a_j}(x_\tau,1-\tau)=0\quad(j\ge1)                 \tag{11}
\]

from the exact protected support hole, or the quantitative alternative

\[
 \sum_{j\ge1}|(V_{j,a_j})_\theta(x_\tau,1-\tau)|
 \le\tfrac12 e_0\tau^{-A}.                              \tag{12}
\]

The local stress/cone inequalities do not by themselves imply (10), and a
velocity upper bound cannot supply it.  Until (10) is proved on the actual
same field, the expanding-ball comparison remains a conditional implication
and no competitor exclusion or blow-up conclusion follows.

## 6. Reproduction

The exact controls check the two positive-slope stage gates, logarithmic
half-gain, Cartesian factor, shrinking physical support distances and force
jet order ledger.  They do not evaluate continuum suprema or prove (10).

```sh
python3 research/navier_stokes_cascade/results/core_scale_compatible_corner_force_replay/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_scale_compatible_corner_force_replay/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_scale_compatible_corner_force_replay/evidence.json \
  --output fresh.json
```
