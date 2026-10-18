# Same-canonical-field pressure recovery and global comparison

Parent: `634356dff6d7b6c0b350c714cd11403d52c216ae`.
Status: conditional canonical comparison, relative to retained analytical
construction and positive-core premises. Historical evidence is unchanged.
The parent's canonical-force CI step passed on 18 October 2026; its whole
workflow was still running. This is not an unconditional blow-up theorem.

We check the ordinary hypotheses of the general comparison argument on the
SAME canonical field and force from `canonical_force_join`. The historical
force and comparison fields are never identified with this field. We also
derive a compact mixed-pressure moment expansion and exhibit a smooth,
compact, divergence-free witness showing why its leading moment cannot be
discarded from support and incompressibility alone.

## 1. Check the canonical comparison interface

Fix T<1. Relative to retained exact-core/background and finite analytical
defining, existence, inverse and gap lemmas, the parent provides

\[
u=\operatorname{curl}(cA_c)+cB_c,\quad p=cp_c,\quad
u_t+(u\cdot\nabla)u=\Delta u-\nabla p+f,\quad\operatorname{div}u=0. \tag{1}
\]

These are its actual complete fixed-gauge representatives, own roots and
mean targets, not the earlier family's representatives. Both base pressure
tails, all six covariance entries, axial velocity, changing intermediate
pressure and every localization term are retained. On [0,T], u,p are smooth
with support in one fixed compact K; u(0)=0. The parent's ONE compact smooth
global force is the full residual for t<1. Only the force was extended past
one. Its energy bound gives uniform preterminal L2 control. Compact
smoothness gives finite sup_[0,T] ||u||6 and L_T=sup_[0,T] ||grad u||infinity.

Let v,P be any smooth whole-space solution with this SAME zero datum and
SAME force and with sup_[0,T] ||v||2 finite for every T<1. No global bound
on grad v, v_t, dissipation, or pressure growth is assumed. Smooth means
smooth on compact space-time sets, including t=0. A global bounded-energy
competitor in the stated problem meets these explicit hypotheses.

With w=v-u, pi=P-p and g=vv-uu, uniformly on the fixed interval,

\[
g=w\otimes w+h,\quad h=w\otimes u+u\otimes w,\quad
\|w\|_2\le W_T,\quad\sum\|g_{ij}\|_1\le G_T,\quad
\operatorname{supp}h\subset K,\quad\sum\|h_{ij}\|_1\le H_T.       \tag{2}
\]

These bounds follow from L2 products; compact support of w is not assumed.
The force divergence cancels exactly once on subtraction:

\[
w_t+\operatorname{div}g=\Delta w-\nabla\pi,\qquad
-\Delta\pi=\partial_i\partial_jg_{ij}.                          \tag{3}
\]

Individually -Delta p=DD(uu)-div f, and similarly for P. The individual
pressures are not replaced by unforced Riesz formulas; div f need not vanish.

## 2. Recover the gradient without a competitor pressure hypothesis

Let T_ij=R_iR_j, multiplier -xi_i xi_j/|xi|^2, and pi*=sum T_ij g_ij.
L1 stress has bounded Fourier transform, so pi* is uniformly in H^-2:
the Fourier weight is integrable for s>3/2, strictly. Poisson alone leaves
a harmonic gradient and is insufficient. For a compact smooth time test a,

\[
\int a\nabla\pi\,dt=\Delta\int aw\,dt+\int a'w\,dt
                       -\operatorname{div}\int ag\,dt.         \tag{4}
\]

The terms are in H^-2,L2,H^-3. Subtract int a grad pi*: the result is
harmonic and belongs to H^-3. Its Fourier transform is a weighted L2
FUNCTION supported at the single point zero and therefore vanishes.
Thus grad pi=grad pi* as a space-time distribution. This argument derives
temperedness of the averaged gradient; it does not assume temperedness of P.

For each ball, stress truncations converge in L1 and their transforms in
H^-2. The localized main Riesz term and the integrable commutator below
identify pi* as locally integrable in space-time. Consequently testing the
gradient equality against a(t) chi_R w is legitimate, and the two pressure
fluxes agree almost everywhere. Spatial constants in pressure have zero
flux. Compact smoothness permits integration of the energy identity from
t=0. These steps apply the ordinary lemma in `global_solution_comparison`,
not its historical special field.

The sole singular-integral input is Riesz Lp boundedness, at p=3/2.
It was rechecked on 18 October 2026 in the author-hosted
[Knapp, Advanced Real Analysis](https://www.math.stonybrook.edu/~aknapp/download/a2-1-realanal-clickable.pdf),
III.5 Theorem 3.26 and III.6 Example 1, printed pp83,93.
The retained manuscript Section 10.3, Lemma 10.5, is a comparison source;
its final theorem is never an input. Manuscript SHA256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.

## 3. New mixed-pressure moment expansion; no automatic cancellation

Choose K subset B_L and R>=max(1,2L). On R<=|x|<=2R put
M_ij(t)=int h_ij(y,t)dy. The off-support Riesz kernel is

\[
k_{ij}(x)=\frac{3x_ix_j-\delta_{ij}|x|^2}{4\pi|x|^5},\qquad
\pi_h=\sum T_{ij}h_{ij}=k(x):M+\mathcal E(x),\qquad
|\mathcal E(x)|\le\frac{384}{4\pi}L H_T R^{-4}.                 \tag{5}
\]

Here the local identity distribution is zero outside K. For any unit e,
direct differentiation gives

\[
4\pi D_e k_{ij}(x)=\frac{3}{|x|^4}
 [e_i n_j+n_i e_j+\delta_{ij}(e\cdot n)-5n_i n_j(e\cdot n)],
\quad n=x/|x|.                                                 \tag{6}
\]

Each term is bounded, so |D_e k_ij|<=24/(4*pi*|x|^4). The segment
x-sy stays at distance at least R-L>=R/2; the mean-value formula gives
24*16*|y|/(4*pi)*R^-4. Integrate componentwise against h to prove (5).
The bound is uniform in t on this fixed interval. No Taylor expansion of
an actual infinite fluid or numerical far-field fit is being substituted.

Trace k=0, so the leading term depends on M°=M-tr(M)I/3. However divergence
freedom does NOT imply M°=0. Choose a nonzero smooth radial compact bump
psi(x)=rho(|x|^2), V=(-x2,x1,0)psi, and, solely as a tensor witness, u=w=V.
Then div V=0 and parity/rotation give int ViVj=a diag(1,1,0), a>0. Hence

\[
M=2a\operatorname{diag}(1,1,0),\qquad
k(Re_3):M=-\frac{4a}{4\pi R^3}\ne0.                            \tag{7}
\]

This is an actual smooth compact divergence-free tensor construction, not
a pair of same-force Navier-Stokes trajectories. It disproves only the
proposed inference from compactness and divergence freedom to moment
isotropy. It is not an obstruction to the comparison theorem. A moment
isotropy condition, if independently available, would improve the mixed
flux from R^-5/2 to R^-7/2; it is neither supplied nor needed here. Even
then the coarse Laplacian cutoff error below remains R^-2.

## 4. Join the complete expanding-ball bound to the canonical field

Take phi smooth, one on B1, supported in B2, 0<=phi<=1, phi_R=phi(x/R),
chi_R=phi_R^8. Define E_R=int chi_R|w|², D_R²=int chi_R|grad w|² and
B_R=||phi_R^4 w||6. Compact Sobolev gives B_R<=C(D_R+R^-1 W_T), using
only derivatives on the fixed ball. Retain the quadratic pressure split
phi_R^4 pi_q=sum T_ij(phi_R^4 wi wj)+sum[phi_R^4,T_ij](wi wj).
Its main L^(3/2) norm is <=C B_R W_T. The commutator kernel has absolute
bound C |z|^-3 min(|z|/R,1). Its L^(4/3) radial norm power is

\[
R^{-4/3}\int_0^R r^{-2/3}dr+\int_R^\infty r^{-2}dr=4/R.          \tag{8}
\]

Thus the commutator norm is <=C R^-3/4 W_T². Pair with the remaining
R^-1 phi_R^3 w in L3 and L4, respectively, to obtain quadratic pressure
flux <=C_T(R^-1 B_R^(3/2)+R^-7/4 B_R^(3/4)). The eighth power is needed
by this factorization; seven fails. No L1-to-L1 pressure assertion is used.
The mixed pressure estimate |pi_h|<=16/(4*pi) H_T R^-3, including its
nonzero leading moment, gives mixed flux <=C_T R^-5/2.

The COMPLETE difference energy identity is

\[
\tfrac12E_R'+D_R^2=-\int\chi_R(w\cdot\nabla)u\cdot w
 +\tfrac12\int|w|^2\Delta\chi_R
 +\tfrac12\int|w|^2v\cdot\nabla\chi_R
 +\int\pi w\cdot\nabla\chi_R.                                 \tag{9}
\]

For R>=2L, u=0 and v=w on the cutoff annulus. Transport has bound
C_T R^-1 B_R^(3/2), strain at most ||grad u||infinity E_R, and the
Laplacian term at most C_T R^-2. Young absorption of powers 3/2 and 3/4
gives decays 4 and 14/5; the Sobolev offsets decay 5/2. Retaining all terms,

\[
\tfrac12E_R'+\tfrac12D_R^2\le L_T E_R+C_T R^{-2},\quad E_R(0)=0,
\qquad E_R(t)\le2C_T T e^{2L_T T}R^{-2}.                        \tag{10}
\]

The integrated local dissipation is controlled as well: integrate the first
inequality and insert the energy bound, discarding only -E_R(t)<=0, to get

\[
\int_0^tD_R(s)^2ds
 \le 2C_T t[1+2L_T T e^{2L_T T}]R^{-2},\qquad t\le T.           \tag{11}
\]

Fix T first, then let R grow. Every fixed ball eventually lies in chi_R=1,
so (10) gives v=u throughout [0,T]. There is no monotonic-cutoff requirement
and no T-uniform control of L_T or C_T is asserted. Since T<1 was arbitrary,
the velocities agree for all t<1. No competitor energy inequality, global
derivative bound, pressure decay, future velocity or future pressure was used.

## 5. Conditional exclusion and the actual next proof gate

Only after the radius limit use the parent's protected canonical core path.
With tau=1-t, x_tau=(sqrt(2 Xin tau),0,0), it lies in the SAME plateau when
tau<=min(tau0/2,rL²/(8 Xin)); rL is the canonical localization radius, not
an old pulse radius. Retained positive-core growth gives

\[
|u_\theta(x_\tau,1-\tau)|\ge(e0/2)\tau^{-Aexp},\quad e0,Aexp>0.  \tag{12}
\]

A global smooth competitor is bounded on a compact neighborhood of (0,1),
contradicting (12). This closes the downstream comparison/exclusion
IMPLICATION for this canonical family, relative to the retained construction
and positive-core lemmas. It does not independently establish those lemmas,
an unconditional admissible-force construction, or a blow-up result.

Next audit the exact-core/global-backward-pressure defining-map identity.
Trace `core_contraction_audit`, `core_all_order_foundation`, the same-root
moment/background records and their canonical use before deciding whether
that identity is closed. Local contraction for an input pressure family is
not, on its own, identification with the completed global pressure. Then
discharge retained background and finite defining/existence/inverse/gap
premises at the unchanged parameters. The next checkpoint must establish
the first outstanding identity/inequality or document its precise failure;
repeating downstream conditional implications is insufficient.

Exact rational controls below check (6), (7), point-tensor remainders,
support separation, cutoff and radius exponents, and the complete local
energy identity. Point measures and local polynomial fixtures are explicitly
not whole-space same-force flows. Tests and source bindings are reproducible
engineering evidence, not continuum proofs; (1)--(12) supply the analytical
argument under its explicit hypotheses. No historical evidence is rewritten.

```sh
python3 research/navier_stokes_cascade/results/canonical_comparison_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/canonical_comparison_join/audit.py \
 --verify-record research/navier_stokes_cascade/results/canonical_comparison_join/evidence.json \
 --output /tmp/canonical_comparison_fresh.json
```
