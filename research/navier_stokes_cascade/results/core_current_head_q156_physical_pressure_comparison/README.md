# Q156: physical fixed-gauge pressure and expanding-ball comparison

Parent: Q155, `c4955037ccbc0cd039d8df071c93c26f85eaf6c7`.
The live branch, Q155 field/force/energy record, the frozen Q10 comparison,
comparison foundation and current CI status were read first. Historical
evidence, the Q155 physical field, fixed pressure gauge, compact force,
supports, accepted bands, diagonal cutoffs and `C^-20` width are unchanged.

Q155 supplies, conditionally on its named continuum providers, one
divergence-free field `u`, complete fixed-gauge pressure `p=cp0`, compact
smooth force `f`, zero datum and finite preterminal energy. This checkpoint
applies the refined pressure/comparison lemma to those exact representatives.
It also corrects the downstream scope: the comparison closes, but Q155's
frozen interface does not bind constants and a path for a quantitative
same-physical-field core lower bound. No terminal exclusion or blow-up is
therefore claimed here.

## 1. Difference equation with complete pressure memory

Fix `T<1`. Let `v,P` be smooth on compact subsets of
`R^3 x [0,T]`, solve Navier--Stokes with the same force and zero datum, and
satisfy

\[
             \sup_{0\le t\le T}\|v(t)\|_2<\infty .       \tag{1}
\]

No global derivative, pressure-decay, dissipation or boundary-flux bound is
assumed. Put

\[
 w=v-u,\qquad \pi=P-p,\qquad
 g=v\otimes v-u\otimes u
   =w\otimes w+w\otimes u+u\otimes w .                 \tag{2}
\]

The common force cancels exactly once:

\[
 w_t+\operatorname{div}g=\Delta w-\nabla\pi,
 \qquad \operatorname{div}w=0,
 \qquad-\Delta\pi=\partial_i\partial_jg_{ij}.          \tag{3}
\]

The individual equations retain the complete force divergence,

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f,
 \qquad
 -\Delta P=\partial_i\partial_j(v_i v_j)-\nabla\cdot f. \tag{4}
\]

Thus neither individual pressure is replaced by an unforced Riesz formula,
and `div f` is not subtracted a second time.

## 2. Finite energy removes the harmonic pressure gradient

On the fixed interval, `w` is uniformly `L2` and `g` uniformly `L1`. Define

\[
                  \pi_*=\sum_{i,j}R_iR_jg_{ij}.         \tag{5}
\]

The Fourier transform of an `L1` stress is bounded, so `pi_*` is uniformly
`H^-2`; the underlying embedding requires `s>3/2`, not the endpoint. For a
compact smooth time test `a`, (3) gives

\[
 \int a\nabla\pi\,dt
  =\Delta\int aw\,dt+\int a'w\,dt
   -\operatorname{div}\int ag\,dt\in H^{-3}.           \tag{6}
\]

Subtracting `grad pi_*` leaves a harmonic `H^-3` distribution. Its Fourier
transform is a weighted `L2` function supported at the single point zero,
and hence vanishes. Therefore

\[
                      \nabla\pi=\nabla\pi_*             \tag{7}
\]

as a space-time distribution. Spatial stress truncation and local
smoothness justify the compact pressure-flux tests. A time-dependent
spatial constant remains possible, but its flux against `w grad chi_R`
vanishes by incompressibility. This recovers the competitor pressure only
to the extent required by comparison and assumes no pressure temperedness.

The finite-energy condition is essential. Without it,
`v=t e1, P=-x1` and `u=p=0` have the same zero datum and force but a nonzero
harmonic pressure gradient.

## 3. Refined compact mixed-pressure split

Write

\[
 g=w\otimes w+h,\qquad h=w\otimes u+u\otimes w.       \tag{8}
\]

Q155 gives one fixed compact spatial support `K subset B_L` for `u`, hence
also for `h`. For `R>=max(1,2L)` and `R<=|x|<=2R`, kernel separation gives

\[
             |\pi_h(x,t)|\le C H_T R^{-3}.             \tag{9}
\]

Pairing with `R^-1 w` on the annulus and using only the global `L2` bound
gives an `O(R^-5/2)` flux. No isotropy of the compact mixed moment is used.
This is the refinement that removes the older unsplit
`R^-1 B_R^(1/2)` term and its `R^-4/3` Young remainder.

Choose `chi_R=phi(x/R)^8`, where `phi=1` on the unit ball and vanishes
outside the radius-two ball, and set

\[
 E_R=\int\chi_R|w|^2,\qquad
 D_R^2=\int\chi_R|\nabla w|^2,\qquad
 B_R=\|\phi_R^4w\|_6 .                                \tag{10}
\]

Local Sobolev gives `B_R<=C(D_R+R^-1 W_T)`. For the quadratic pressure,
split

\[
 \phi_R^4\pi_q
 =\sum R_iR_j(\phi_R^4w_iw_j)
  +\sum[\phi_R^4,R_iR_j](w_iw_j).                      \tag{11}
\]

The standard `L^(3/2)` Riesz bound controls the first term. The commutator
kernel is bounded by

\[
             C|z|^{-3}\min\{|z|/R,1\}.                 \tag{12}
\]

At exponent `4/3`, its exact radial integral is

\[
 R^{-4/3}\int_0^Rr^{-2/3}\,dr
  +\int_R^\infty r^{-2}\,dr={4\over R},                \tag{13}
\]

so its norm is `O(R^-3/4)`. Both the near cancellation and far cap are
required. Interpolation of the remaining `phi_R^3 w` factor yields

\[
 \left|\int\pi_qw\cdot\nabla\chi_R\right|
 \le C_T\left(R^{-1}B_R^{3/2}
              +R^{-7/4}B_R^{3/4}\right).              \tag{14}
\]

The eighth cutoff power is sufficient because
`k/2-1>=3k/8` exactly requires `k>=8`; the seventh power fails this pairing.

## 4. Complete expanding-ball estimate

The exact weighted difference identity is

\[
\begin{aligned}
 {1\over2}E_R'+D_R^2={}&
 -\int\chi_R(w\cdot\nabla)u\cdot w
 +{1\over2}\int|w|^2\Delta\chi_R\\
 &+{1\over2}\int|w|^2v\cdot\nabla\chi_R
 +\int\pi w\cdot\nabla\chi_R .                       \tag{15}
\end{aligned}
\]

On the cutoff annulus `u=0`, so `v=w` there. This support separation is
made before estimating transport. The strain term is at most
`||grad u||_infinity E_R`; the cutoff-Laplacian term is `O(R^-2)` from the
uniform `L2` bound. For `R^-b D_R^a`, `0<a<2`, Young absorption leaves
radius power `2b/(2-a)`. The two nonconstant rows in (14) and transport are

| term | `a` | `b` | Young decay | Sobolev offset |
|---|---:|---:|---:|---:|
| quadratic pressure/transport | `3/2` | `1` | `4` | `5/2` |
| quadratic commutator | `3/4` | `7/4` | `14/5` | `5/2` |

Together with the `R^-5/2` compact mixed flux, the coarse `R^-2`
cutoff-Laplacian term is limiting. After absorption,

\[
 \boxed{
 {1\over2}E_R'+{1\over2}D_R^2
 \le\|\nabla u(t)\|_\infty E_R+C_TR^{-2},
 \qquad E_R(0)=0 .}                                   \tag{16}
\]

This is an analytical consequence of the stated hypotheses. The executable
ledger checks the exponents; it is not a continuum proof. The `R^-2` rate
is the result of this coarse uniform-`L2` estimate, not a universal optimal
rate for an individual solution.

## 5. Radius limit and the precise remaining obstruction

Let `L_T=sup_[0,T] ||grad u||_infinity`, finite by Q155's compact
preterminal smoothness. Gronwall gives

\[
 E_R(t)\le2C_TT\exp(2TL_T)R^{-2}.                     \tag{17}
\]

Fix `T<1` first and send `R` to infinity. Every fixed ball eventually lies
where `chi_R=1`, so `w=0` on that ball. The cutoff family need not be
monotone and the constants need not stay bounded as `T` approaches one.
Since `T<1` is arbitrary,

\[
                  v=u\quad\hbox{on }\mathbb R^3\times[0,1). \tag{18}
\]

Equation (18) is the completed Q156 comparison. It does not by itself give
a terminal contradiction. The exact missing same-field interface is

\[
 \boxed{\exists e_0,A,\tau_0>0:\quad
 x_\tau\in\{c=1\},\qquad
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-A}
 \quad(0<\tau<\tau_0).}                               \tag{19}
\]

Q155's prose says the original lower bound is unchanged, but its frozen
record does not identify a quantitative path, constants and same-physical-
field provider for (19). Importing an older canonical/current-field bound
would conflict with Q149--Q151's explicit distinction between that family
and the inherited physical primary. Until (19) is bound to the exact
Q153--Q155 field, Q156 claims neither global smooth competitor exclusion
nor blow-up.

The next obligation is to derive or bind (19), including its path and
localization plateau, for this exact physical diagonal. Only then may (18)
be combined with local boundedness of a global smooth competitor near the
terminal point. No future velocity or pressure is constructed.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q156_physical_pressure_comparison/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q156_physical_pressure_comparison/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q156_physical_pressure_comparison/evidence.json \
  --output q156_replay.json
```
