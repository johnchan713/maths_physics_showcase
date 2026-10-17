# Current compact-force comparison with an R^-2 infinity remainder

Parent: `44fc3767b90e07f11df62b4726d2d48a3d8e85ec` (`force_foundation`).
This checkpoint joins the SAME current field, full pressure, zero datum and
compact smooth force to pressure recovery and comparison with every global
smooth finite-energy competitor. The special construction and its protected
core growth remain relative to the retained exact-core and finite
defining/existence/inverse/gap lemmas. Those dependencies are not discharged
by this comparison; no unconditional breakdown or unforced singularity is
established.

The general argument in `global_solution_comparison` and its older
conditional application in `current_solution_comparison` are preserved.
The new analytical estimate uses the compact support of the mixed stress
to improve the complete expanding-ball error from **R^-4/3 to R^-2**.
It also checks the hypotheses on the current force foundation, rather than
substituting the historical localized field or asserting that a later
construction inherits a comparison interface automatically.

## 1. Current hypotheses and the exact pressure difference

Fix `T<1`. Relative to the retained construction lemmas, the parent supplies

\[
 u=\operatorname{curl}(cA)+cB,\qquad p=cp_0,\qquad
 u_t+(u\cdot\nabla)u=\Delta u-\nabla p+f,\qquad\operatorname{div}u=0.
\]

The full potential, direct swirl, both background pressure tails, all wave
and mean pressure differences, axial terms and moving-moment repairs stay
in the fixed gauge. `u,p` have one fixed compact spatial support `K`, are
smooth on `R3 x [0,T]`, and `u(0)=0`. The force is the parent's ONE smooth
compact global force, obtained from the complete twelve-term localized
residual and smooth terminal force-jet gluing. The parent's energy bound
supplies `sup_(t<1) ||u||2<infinity`. Compact smoothness also gives finite
`sup_[0,T] ||u||6` and `L_T=sup_[0,T] ||grad u||infinity`.

Let `v,P` be any smooth solution with this same force and datum and with
`sup_[0,T] ||v||2<infinity`. These are ordinary explicit hypotheses for the
general comparison lemma. The special `u,p,f` satisfy them only relative
to the retained construction lemmas. No global derivative or dissipation
bound, pressure growth condition, or energy inequality for `v` is assumed.
Smoothness of `v,P` means smoothness on compact sets.

Set `w=v-u`, `pi=P-p`, `g_ij=v_i v_j-u_i u_j`, and write

\[
 g=w\otimes w+h,\qquad h=w\otimes u+u\otimes w.
\]

On this interval `||w||2<=W_T` and `sum ||g_ij||1<=G_T` uniformly for finite
constants. Since `u=0` outside `K`, **h is compactly supported in K**, even
though `w` need not decay pointwise or have compact support. Cauchy--Schwarz
gives `H_T=sup_t sum ||h_ij||1<infinity`. The exact difference equations are

\[
 w_t+\operatorname{div}g=\Delta w-\nabla\pi,
 \qquad-\Delta\pi=\partial_i\partial_jg_{ij}.
\]

Individually `-Delta p=DD(uu)-div f`, and likewise for `P`.
The common force divergence cancels **once** on subtraction. It need not
be zero. The individual constructed pressure is not replaced by an
unforced Riesz formula.

## 2. Recover the pressure gradient before estimating its flux

Let `T_ij=R_i R_j` have multiplier `-xi_i xi_j/|xi|^2`, and put
`pi_star=sum T_ij g_ij`. An L1 stress has bounded Fourier transform, hence
`pi_star` is uniformly in `H^-2`; the embedding needs `s>3/2`, not the
critical endpoint `s=3/2`. Poisson's equation says only that
`Delta(pi-pi_star)=0`, so it does not alone remove a harmonic gradient.

For every compact smooth time test `a` inside `(0,T)`, momentum gives

\[
 \int a\nabla\pi\,dt
 =\Delta\int aw\,dt+\int a'w\,dt-\operatorname{div}\int ag\,dt.
\]

The right side is in `H^-3`: its respective terms are in `H^-2`, `L2` and
`H^-3`. Also `int a grad pi_star` is in `H^-3`. Their difference is
harmonic, and its Fourier transform is a weighted L2 **function** supported
at `{0}`. It therefore vanishes. This proves
`grad pi=grad pi_star` as a space-time distribution without assuming that
the competitor pressure itself is tempered. Spatial truncations of the L1
stress and local smoothness identify the Riesz expression with the local
integrable pressure flux. Testing against `a(t) chi_R w` is then legitimate;
a spatially constant pressure difference contributes zero flux. At fixed
radius, compact smoothness permits integration down to `t=0`.

This is the pressure-recovery argument of the preserved general lemma;
its source and frozen controls are reproduced. The standard external input
for the next step is Riesz boundedness on `Lp`, `1<p<infinity`, used at
`p=3/2`. It was checked in the author-hosted
[Knapp, Advanced Real Analysis](https://www.math.stonybrook.edu/~aknapp/download/a2-1-realanal-clickable.pdf),
III.5 Theorem 3.26 and III.6's Riesz example, printed pp83,93, on 2026-10-17.
The support refinement below is derived here.

## 3. Estimate compact mixed pressure directly on the outer annulus

Choose `K subset B_L`, and `R>=max(1,2L)`. Put
`pi_h=sum T_ij h_ij`, `pi_q=sum T_ij(w_i w_j)`. Outside `K`, the local
identity distribution in `T_ij` is zero and its ordinary kernel is

\[
 k_{ij}(z)=\frac{3z_i z_j-\delta_{ij}|z|^2}{4\pi |z|^5}.
\]

Each component is at most `2/(4*pi*|z|^3)`. On
`A_R={R<=|x|<=2R}` and `y in K`, `|x-y|>=R-L>=R/2`; therefore

\[
 |\pi_h(x,t)|\le\frac{16}{4\pi}R^{-3}
       \sum_{i,j}\|h_{ij}(t)\|_1.
\]

For a fixed smooth ball cutoff `phi`, one on the unit ball, supported in
the radius-two ball, with `0<=phi<=1`, set `chi_R=phi(x/R)^8`.
Its gradient is supported in `A_R` and bounded by `C/R`. Thus

\[
 \left|\int\pi_h w\cdot\nabla\chi_R\right|
 \le C H_T R^{-4}|A_R|^{1/2}W_T
 \le C_T R^{-5/2}. \tag{1}
\]

This needs only the mixed tensor's L1 norm and support, with no pointwise
decay of `w`. In particular it avoids applying the old main-pressure bound
`C_T(B_R+1)` to the entire stress: that unnecessary constant term produced
the old `R^-1 B_R^1/2` contribution and its `R^-4/3` Young remainder.
The older estimate remains valid; (1) is a sharper use of its hypotheses.

## 4. Quadratic pressure, transport and all cutoff terms

Define

\[
 E_R=\int\chi_R|w|^2,\quad
 D_R=\left(\int\chi_R|\nabla w|^2\right)^{1/2},\quad
 B_R=\|\phi_R^4w\|_6,\quad\phi_R=\phi(x/R).
\]

The compact Sobolev estimate gives
`B_R<=C(D_R+R^-1 W_T)` using only derivatives on a fixed radius ball.
Split only the quadratic pressure:

\[
 \phi_R^4\pi_q=\sum T_{ij}(\phi_R^4w_iw_j)
                 +\sum[\phi_R^4,T_{ij}](w_iw_j).
\]

Riesz boundedness and Holder give main norm
`||sum T_ij(phi_R^4 wi wj)||_(3/2)<=C B_R W_T`.
The commutator kernel is bounded by
`C |z|^-3 min(|z|/R,1)`. Its `4/3` radial norm power is the sum

\[
 R^{-4/3}\int_0^R r^{-2/3}dr+\int_R^\infty r^{-2}dr=4/R,
\]

so the commutator norm in `L^(4/3)` is at most `C R^-3/4 W_T^2`.
Both the near cancellation and the far cap are required. This argument
uses neither an L1 pressure estimate from L1 stress nor a global gradient
bound for `v`. Interpolating the remaining factor `phi_R^3 w` between
`phi_R^4 w in L6` and `w in L2` gives

\[
 \left|\int\pi_q w\cdot\nabla\chi_R\right|
 \le C_T\left(R^{-1}B_R^{3/2}+R^{-7/4}B_R^{3/4}\right). \tag{2}
\]

The eighth cutoff power is retained: the L4 commutator pairing requires
`k/2-1>=3k/8`, or `k>=8`. This is a constraint of this factorization.

The COMPLETE weighted difference-energy identity is

\[
\begin{aligned}
 \tfrac12E_R'+D_R^2={}&-\int\chi_R(w\cdot\nabla)u\cdot w
 +\tfrac12\int|w|^2\Delta\chi_R\\
 &+\tfrac12\int|w|^2v\cdot\nabla\chi_R
 +\int\pi w\cdot\nabla\chi_R.
\end{aligned}
\]

On `supp grad chi_R` the current velocity is zero, so `v=w` there.
Only after this separation is the transport flux bounded by
`C_T R^-1 B_R^3/2`. The strain is at most `||grad u||infinity E_R`,
and the Laplacian cutoff flux is bounded by `C R^-2 W_T^2`.
The pressure flux equals that of `pi_q+pi_h` by gradient recovery.

For `0<a<2`, maximize in `D>=0` to obtain
`R^-b D^a <= delta D^2+C_(a,delta) R^(-2b/(2-a))`.
After substituting the Sobolev bound, the only nonconstant dissipation
powers in (2) and transport have `(a,b)=(3/2,1),(3/4,7/4)`, hence Young
decays `4,14/5`. The Sobolev offsets both decay as `R^-5/2`, as does (1).
The Laplacian remainder has power `2`. Absorb the finitely many quadratic
contributions with total coefficient at most `1/2`; every remaining term
is bounded by `C_T R^-2`. Thus, with constants independent of radius,

\[
 \boxed{\tfrac12E_R'+\tfrac12D_R^2
 \le\|\nabla u\|_\infty E_R+C_T R^{-2},\qquad E_R(0)=0.} \tag{3}
\]

The `R^-2` power is the limit of the **coarse uniform-L2 Laplacian bound
used here**, not a universal optimal comparison rate for each fixed
solution. Indeed the family `w_R=R^-3/2 W(x/R)` preserves L2 norm and
scales a cutoff Laplacian flux exactly as `R^-2`. A divergence-free compact
`W` can be supported in a patch of the annulus where `Delta(phi^8)` has
fixed nonzero sign. This varying family tests only what those coarse norm
data can imply, not the behavior of one fixed solution at all radii.

## 5. Limits and the conditional exclusion implication

For fixed `T<1`, compact smoothness makes `L_T` finite. Gronwall gives

\[
 E_R(t)\le2C_T R^{-2}\int_0^t
   \exp\left(2\int_s^t\|\nabla u(a)\|_\infty da\right)ds
 \le2C_T T e^{2TL_T}R^{-2}. \tag{4}
\]

The bound also covers `L_T=0` with no division. Each fixed ball lies where
`chi_R=1` for all sufficiently large `R`. Letting radius tend to infinity
in (4) forces `w=0` there for every `t<=T`; continuity identifies the fields.
No monotone family of cutoffs is needed. Since `T` was arbitrary,
the competitor agrees with the SAME current constructed velocity on `[0,1)`.

Only now let `t->1` along the parent's protected core path
`x_tau=(sqrt(2 Xin tau),0,0)`, with
`tau<=min(tau0/2,r0^2/(8 Xin))`. It remains on the current cutoff plateau,
and the inherited complete-field bound is
`|u_theta(x_tau,1-tau)|>=(e0/2) tau^-A`, with positive `e0,A`.
A global smooth competitor would be bounded on one compact neighborhood
of `(0,1)`, contradicting this growth. The implication is relative to the
retained construction and core-growth lemmas. Neither `C_T` nor `L_T` is
asserted bounded as `T->1`; no future velocity or pressure is constructed.

## 6. Checks, preserved evidence and the next foundation gate

Exact rational local polynomial controls check all nine stress components,
the full energy flux-divergence identity at three positive viscosities,
and the difference pressure equation with a manufactured residual source.
Every energy term is nonzero in the fixture. These local polynomials are
not same-force whole-space finite-energy trajectories; the comparison uses
the exact source cancellation only under its actual same-force hypothesis.
Sixty rational exterior point-tensor kernel controls check the separated
bound, Hessian sign and coefficient, with `1/(4*pi)` tracked separately.
Those point stresses are diagnostic measures, not actual L1 velocity
products. Exact exponent and support controls detect an insufficient
separation radius, critical kernel endpoints, a seventh cutoff power,
the old unsplit linear pressure term and lost Young powers.

The parent and general comparison records are reproduced with their own
comparators. The new record binds eight own sources, eight dependency
records and two imported helpers; 1,015 inherited provenance bindings stay
unchanged. Selected exact diagnostics and passing tests do not prove a
continuum supremum or discharge the retained analytic premises. The
continuum argument for (1)--(4) is given above.

Next audit the remaining construction premises independently, starting with
the exact core's infinite-space contraction and recovery of the original
equations, and its identification with the completed global backward
pressure. The older core results explicitly scoped that pressure-family
identification separately. Trace that dependency through the current
moment/background records before deciding whether it is closed; do not
label it either proved or obstructed merely from old scope text. Then
resolve the retained finite background/pulse defining, existence, inverse
and gap lemmas at the same parameters. The next substantive result must
prove the first outstanding identity/inequality or expose its precise
obstruction. Further downstream implication-only repetition is insufficient.

Reproduce using fresh paths:

```sh
python -m unittest -v test_audit.py
python audit.py --output /fresh/comparison-foundation.json --verify-record evidence.json
```
