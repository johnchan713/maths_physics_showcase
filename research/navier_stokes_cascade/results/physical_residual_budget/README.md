# Full base residual and the first correction budget

Status: `physical-residual-derived-with-open-C3-and-correction-bounds`.

The full **unlocalized axisymmetric base** residual and its first positive-order
correction equations are derived below. They retain radial acceleration,
radial and axial viscosity, pressure, cylindrical curvature, the mixed axial
derivatives and all products in a finite coefficient sum. A conditional
constant bound is given for each base component. The actual post-modulation
derivative constant has not yet been bounded. In particular the preceding
C2 estimate does **not** bound the third angular derivative entering radial
axial diffusion. Milestone 3 of the blow-up plan remains open.

Parent: `753501b7c205e6c9073e1c939bbc340dfbf65dc1`.
The fixed profile frequency remains **N=1+floor(H^32)** from the
[conditional repair estimate](../repair_state_bounds/README.md). No later
physical scale is used to select it. All numerical fields below are
manufactured fixtures; the actual enormous-frequency field is not resolved.
Earlier analytic dependencies have not received independent mathematical
review. There is no verified PDE correction sum, smooth force or blow-up.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The notation is aligned with (4.2)-(4.7), (4.11)-(4.12), (5.1)-(5.6) and
(5.26). The following is a direct chain-rule calculation and derivative
budget, not a verification of the manuscript's positive-order existence or
summation claims. Its tangential stress operator in (5.41) must be
distinguished from the divergence of a full symmetric wave covariance.

## 1. Coordinates and complete derivative operators

Write Aphys=1/2+h, D=1/2-h, 0<h<1/2, R=sqrt(2X), d=1-eta^2 and
L=1-2h eta^2. Aphys is unrelated to the earlier large envelope A=C^4096.
The actual construction uses h<.01. For t<1,

\[
q-z^2q^{2h}=1-t,\quad r=q^{1/2}R,\quad z=q^D\eta,
\quad\epsilon=q^{2h}.
\]

At fixed physical coordinates,

\[
\partial_t(q^k f)=q^{k-1}T_kf,\qquad
\partial_z(q^k f)=q^{k-D}Z_kf,
\]
\[
T_kf=\frac{-kf+Xf_X+D\eta f_\eta}{L},\qquad
Z_kf=\frac{2k\eta f-2\eta Xf_X+df_\eta}{L}.
\]

Crucially, the second axial derivative is
q^(k-2D) Z_(k-D) Z_k f. The power changes after the first derivative.
For a=-2eta X/L, b=d/L, c=2k eta/L and c'=2(k-D)eta/L, its coefficient is

\[
\begin{split}
Z_{k-D}Z_kf={}&a^2f_{XX}+2abf_{X\eta}+b^2f_{\eta\eta}\\
&+(aa_X+ba_\eta+(c+c')a)f_X
+(bb_\eta+(c+c')b)f_\eta\\
&+(bc_\eta+c'c)f.
\end{split}
\]

Here a_X=-2eta/L, a_eta=-2X(1+2h eta^2)/L^2,
b_eta=-2(1-2h)eta/L^2 and c_eta=2k(1+2h eta^2)/L^2.
This retains the mixed derivative and derivatives of L. Radial scalar
diffusion has coefficient 2(Xf_XX+f_X); radial vector diffusion also contains
the cylindrical curvature term -f/(2X).

## 2. Exact residual of the base field

Use the regular swirl coefficient F=E/R and M=integral_0^X U(x,eta) dx.
The regular radial coefficient v0 and radial flux V are **different**:

\[
V=ru_r=Xv_0=\frac{2\eta XU-2D\eta M-dM_\eta}{L}.
\]

Thus u_r=q^-1/2 V/R, u_theta=q^-Aphys RF, u_z=q^-Aphys U,
p=q^-2Aphys Pi. Incompressibility is exactly V_X+Z_-Aphys U=0.
The base pressure relation is Pi_X=F^2. Let b=-Aphys-1/2, c=-Aphys and

\[
G_\theta=T_bF+V(F_X+F/X)+UZ_bF-2\nu(XF_{XX}+2F_X),
\]
\[
G_z=T_cU+VU_X+UZ_cU+Z_{-2Aphys}\Pi
-2\nu(XU_{XX}+U_X),
\]
\[
\Omega_0=T_0V+V(V_X-V/(2X))+UZ_0V-2\nu XV_{XX}.
\]

The **entire** physical residual R(u,p)=u_t+(u.grad)u-nu Delta u+grad p is

\[
\boxed{\begin{aligned}
R_r&=\frac{q^{-3/2}}{R}
 [\Omega_0-\nu\epsilon Z_{-D}Z_0V],\\
R_\theta&=Rq^{-Aphys-1}
 [G_\theta-\nu\epsilon Z_{b-D}Z_bF],\\
R_z&=q^{-Aphys-1}
 [G_z-\nu\epsilon Z_{c-D}Z_cU].
\end{aligned}}
\]

Before imposing Pi_X=F^2, the radial formula additionally contains
q^(-2Aphys-1/2)(2X Pi_X-2XF^2)/R. This is the exact pressure/centrifugal
cancellation, not an estimate. The other radial terms start one order later
in the normalized pressure balance:

\[
rq^{2Aphys}R_r=\epsilon\Omega_0
-\nu\epsilon^2 Z_{-D}Z_0V.
\]

The leading tangential stress pair T=(Ttheta,Tz) obeys

\[
RG_\theta=-(R\partial_X+2/R)T_\theta,
\qquad G_z=-(R\partial_X+1/R)T_z.
\]

Its physical units are q^(-Aphys-1/2)T. This cancels only the leading
tangential residual when used with the prescribed cylindrical stress
operator. On an inner interval where T=0, axial diffusion remains. Its
physical size q^(-Aphys-1+2h) diverges for the construction's small h.
Likewise q^-3/2 Omega0 can diverge. A relative epsilon gain alone is not
smoothness of the absolute residual.

An illustrative symmetric tensor having only r-theta and r-z cross entries
also has radial divergence
q^-3/2 Z_b Tz. Omitting this term would be wrong. That cross-only tensor is
not positive covariance and is not asserted to be the manuscript's stress
operator. The eventual wave tensor has diagonal and other entries, whose
contributions and interactions must be computed in the correction cycle.

## 3. A conditional bound on every base term

On a compact annulus x0<=X<=x1, |eta|<=1, let K>=1 dominate x1, x0^-1,
L^-1, nu and every ordinary (X,eta) derivative through total order two of
F,U,Pi,V. K refers to the **final modulated and repaired** profiles. It is
not automatically the previous pre-modulation envelope.

The displayed operators give

\[
|T_kf|\le(|k|+2)K^3,\quad |Z_kf|\le(2|k|+3)K^3,\quad
|Z_{k-D}Z_kf|\le32(1+|k|)^2K^5.
\]

For the last inequality the individual coefficient bound sums to
(22+18|k|+4k^2)K^5, below 32(1+|k|)^2K^5. Consequently
|Gtheta|,|Gz|<=20K^4, |Omega0|<=10K^4, and the three axial viscous
coefficients including nu are at most 200K^6, 128K^6, 32K^6 respectively.
For 0<q<=1 this yields component bounds

| Physical component | Sufficient bound | Cancellation already used |
|---|---|---|
| Rr | 42 K^7 q^-3/2 | Pi_X=F^2 |
| Rtheta | 440 K^7 q^(-Aphys-1) | None assumed for Gtheta |
| Rz | 148 K^6 q^(-Aphys-1) | None assumed for Gz |

These are proved conditional inequalities, not an evaluated actual K.
At the axis use the regular Cartesian representation in F,U,v0,Pi; the
annular K contains x0^-1 and is not an axis bound. The identities extend
regularly: V=Xv0 makes Omega0/X and (Z_-D Z_0 V)/X smooth when all required
derivatives exist.

## 4. The missing derivative and the fixed-frequency cost

V contains M_eta. In Z_-D Z_0 V, the coefficient of M_etaetaeta is

\[
-\frac{d^3}{L^3}M_{\eta\eta\eta}.
\]

This comes from (d/L)^2 V_etaeta. The other terms involve at most two
angular derivatives of M. Thus a full base residual needs a bound on the
third angular moment derivative, or the stronger sufficient bound on U's
third angular derivative. A C2 value bound supplies neither.

For an exact control take U_m=m^-2 sin(m eta), independent of X. Its value,
first derivative and second derivative are bounded by m^-2, m^-1, 1 on the
whole angular interval. For its incompressible V at eta=0,

\[
Z_{-D}Z_0V=X(m+6/m).
\]

This tends to infinity although the whole family has uniformly bounded C2
input. This detects a missing estimate; it does not disprove the existence
of the actual smooth profiles or the proposed construction.

With y=log X and fixed phase Ny, an additive primitive modulation obeys

\[
\partial_y^j[N^{-1}\mathcal A(y,\eta,Ny)]
=\sum_{\ell=0}^j {j\choose\ell}N^{\ell-1}
 (\partial_y^{j-\ell}\partial_s^\ell\mathcal A)(y,\eta,Ny).
\]

Already j=2 contains N A_ss. The example sin(Ny)/N has value size 1/N,
first derivative size 1 and second derivative size N. Small state errors
do not imply small viscous errors. Exponential swirl modulation also has
product terms from differentiating exp(A/N). In converting y derivatives,
f_X=f_y/X and f_XX=(f_yy-f_y)/X^2. All these losses must be absorbed in
the actual K or its higher-order counterpart with N fixed.

## 5. First positive-order correction equations

Set lambda=2h and introduce coefficients F1,U1,Pi1 with physical factors
q^lambda relative to the leading field. Their divergence-free radial flux is

\[
V_1=\frac{2\eta XU_1-2(D+\lambda)\eta M_1-d(M_1)_\eta}{L},
\quad M_1=\int_0^X U_1(x,\eta)\,dx.
\]

The D+lambda term is essential. On the retained **stress-free inner
interval**, the equations cancelling order epsilon are

\[
\begin{split}
T_{b+\lambda}F_1+V_0(F_{1,X}+F_1/X)+V_1(F_{0,X}+F_0/X)
+U_0Z_{b+\lambda}F_1+U_1Z_bF_0\\
-2\nu(XF_{1,XX}+2F_{1,X})=\nu Z_{b-D}Z_bF_0,
\end{split}
\]
\[
\begin{split}
T_{c+\lambda}U_1+V_0U_{1,X}+V_1U_{0,X}
+U_0Z_{c+\lambda}U_1+U_1Z_cU_0+Z_{-2Aphys+\lambda}\Pi_1\\
-2\nu(XU_{1,XX}+U_{1,X})=\nu Z_{c-D}Z_cU_0,
\end{split}
\]
\[
\boxed{(\Pi_1)_X=2F_0F_1-\Omega_0/(2X).}
\]

This is a coupled linear system in the current unknowns. Its source terms
use finalized leading profiles. The first pressure correction changes the
radial centrifugal balance; correcting only tangential diffusion misses it.
No solution for the actual profile is asserted here. Regular inner traces,
angular derivative losses, a common interval of solvability, extension into
the annulus and five compatibility moments remain to be proved. On the
annulus positive-order tangential stress is retained rather than setting
these equations to zero everywhere.

The finite-series code also retains the quadratic F1^2, V1 F1_X,
U1 Z_(b+lambda) F1, V1 U1_X and U1 Z_(c+lambda) U1 terms, and the axial
diffusion of order one. Nothing is dropped on the basis of being formally
higher order.

## 6. Physical derivatives, nonlinear corrections and summation

For a smooth Cartesian profile q^k g(x_perp/sqrt(q),eta), a transverse
derivative costs q^-1/2, an axial derivative q^-D and a time derivative
q^-1. For a transverse, k axial and l time derivatives, the loss is
a/2+Dk+l. Constants contain profile derivatives and polynomial factors in
the exponent; they depend on fixed N, but this q-exponent loss does not.

If coefficients through order n are actually cancelled, after separating
every retained annular stress contribution, a sufficient worst-component
remainder budget is

\[
C_{n,m,N}q^{2h(n+1)-(3/2+2h)-a/2-Dk-l}.
\]

This is conditional on coefficient solvability, support/moment restoration
and all required derivative bounds. It is not a residual-improvement proof.
For an illustrative h=1/200, positive decay in this bound requires n>=151
at derivative order zero, n>=251 for one time derivative and n>=351 for
two. These are sufficient bookkeeping thresholds, not necessary physical
barriers or estimates for the construction's much smaller actual h. No
finite n supplies every derivative of a smooth force through t=1.

For any divergence-free perturbation w and pressure correction pi, the
exact nonlinear identity is

\[
R(u+w,p+\pi)=R(u,p)+w_t-\nu\Delta w
+(u\cdot\nabla)w+(w\cdot\nabla)u+\nabla\pi+(w\cdot\nabla)w.
\]

All mixed wave/background products and wave self-products remain here.
For a fixed physical derivative order m, the linear terms need w through
order m+2 and pi through m+1; the quadratic term needs w through m+1.
The all-order coefficient bounds, oscillatory phases and support cutoffs
must control these terms, rather than only a mean covariance.

Incompressibility is also affected by summation cutoffs. A coefficient has
streamfunction Psi_n=q^(D+lambda_n) M_n with
u_r,n=-Psi_n,z/r and u_z,n=Psi_n,r/r. Multiplying the streamfunction by
chi(c_n q) gives

\[
u_{z,n}^{cut}=\chi u_{z,n},\qquad
u_{r,n}^{cut}=\chi u_{r,n}
-\frac{2\eta}{Lr}q^{\lambda_n}(c_nq)\chi'(c_nq)M_n.
\]

The additional radial term preserves divergence exactly. Multiplying
velocity alone does not. The swirl can be cut directly because an
axisymmetric swirl has zero divergence. On cutoff transition regions,
c_n factors are converted using c_n q in the fixed transition interval;
the resulting product and residual commutators still need all derivative
bounds. No convergence of the uncut infinite formal series is assumed, and
no compact summation cutoff sequence has been constructed here.

## 7. Verification and the next obligation

Sixteen focused tests and nine scoped audit gates pass. The 324 manufactured
physical samples cover h=.005 and .2, three positive
viscosities, three q values, three angles eta, three radii and two rotations.
The h=.2 cases check the coordinate algebra outside the small-h construction
regime. Similarity coefficient evaluation agrees with the separate full
Cartesian derivative calculation within 5.5e-15 of each component's term
scale for the base and a finite positive-order sum. The leading stress and
full illustrative symmetric stress identities pass. A scalar-only physical
stencil refines to about 1.1e-9 normalized error.

Six controls detect omitted axial viscosity, a missing axial power shift,
the lost lambda in incompressibility, missing symmetric r-z axial divergence,
swirl curvature and a reversed stress sign. Exact C2 and fixed-N derivative
examples retain the budget failures. A deliberately noncompact exponential
cutoff checks the product rule; it is not the summation construction.

Run from the repository root, with NumPy and SciPy installed:

```sh
python research/navier_stokes_cascade/results/physical_residual_budget/test_audit.py -v
python research/navier_stokes_cascade/results/physical_residual_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/physical_residual_budget/evidence.json \
  --output /tmp/navier_stokes_physical_residual_budget_recheck.json
```

Choose a fresh output path. The evidence pins source files and preceding
audit bytes. Passing it certifies the stated algebraic/numerical checks only.

**Next:** bound the actual post-modulation radial jets and third angular
moment derivative with N fixed; then solve and bound the displayed first
inner correction and its five-moment annular extension. Later orders need
an all-derivative construction, divergence-preserving summation and the full
wave/mean correction cycle. Localization, admissible smooth forcing, energy
and a lower bound for the complete corrected velocity remain separate gates.
