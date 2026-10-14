# Compact mean corrections with the complete wave covariance

Status: `conditional-compact-mean-inverses-and-compatibility-defect-gain`.

The saved actual-band pulse checkpoint was not lost: its 32 tests and frozen
evidence reproduce. This continuation supplies the next **conditional
continuum step**: fixed linear maps for compact mean pressure, a
divergence-free axial increment, auxiliary-time inversion, and a signed
five-equation velocity correction. It derives the exact three compatibility
defects from the **complete** curl covariance and proves a positive gain for
the defects after the slow five-equation update.

This is a mean-correction result, not a positive gain for a complete
Navier--Stokes residual cycle. The current means must satisfy the stated
class hypotheses before the gain applies. The actual global background,
auxiliary rectangles and resulting global mean fields are not numerically
evaluated. The preceding continuum foundation and support geometry remain
conditional. Smooth forcing, complete-field energy and blow-up remain open.

Parent: `6498412c39cf120bea94b90c1f25fc33bf2cea8a`.
Dependencies: [actual-band pulses](../actual_band_pulse_bounds/README.md),
[physical wave interface](../physical_wave_interface/README.md),
[summed background](../background_summation/README.md). Earlier source and
evidence bytes and GOAL.md stay unchanged. The supplied
[proposed manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(6.2), (6.6), (8.1)--(8.27), is the reference being audited; its conclusion
is not a premise. The downloaded PDF has SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`,
matching the preceding record.

## 1. Fixed chart and the full mean equation

Keep the already chosen leading profile, radial frequency N, fixed h and
positive patch exponent lambda. None is chosen again. On a band Q use
epsilon=Q^h, S=ell^2 and s=q/Q in [1/2,2]. Normalize velocities by Q^A,
pressure by Q^(2A), residual by Q^(2A+1/2), and potentials by Q^(A-1/2).
These normalizations include the radial integration length.

In one fixed common torus y, let L=vr.partial_y and Ntime=vt.partial_y.
With M=Lambda_g^i Q^(d_r/2) and c_i=Tg^i Q^(1+h), the evaluated operators are

\[
D_r=\partial_R+d_r M R^{d_r-1}L,\qquad
D_z=\epsilon\partial_Z,\qquad
t_*=-\epsilon\partial_T+c_iN_{time}.
\tag{1}
\]

M, c_i and the common index are constants during coefficient differentiation
in this chart. The inherited graph bounds are M<=C epsilon^(-kappa) times
a power of S, M^(-1)<=C epsilon^kappa times a power of S,
c_i^(-1)<=CS, with kappa=10^-5. The argument below also rederives the
required torus divisor bound.

Write the complete actual velocity as
u=(b+beta,V+v,G+gamma)+w, where w is the full curl wave with angular mean
zero. Put W_ab=<w_a w_b>_theta, retaining every primary/corrector cross
product and every curl remainder. Angular mean is componentwise in the
cylindrical frame and is distinct from auxiliary Haar mean, denoted by a
bar. Define D_e=D_r+e/R and Delta0=D_r^2+R^-1 D_r+D_z^2. With the base's
retained radial--tangential stress pair Sigma, the mean momentum increment
is (D_r p_m-g_r,E_theta,E_z), up to the separately retained flat base error:

\[
\begin{aligned}
g_r={}&-t_*\beta-D_1(2b\beta+\beta^2+W_{rr})
 -D_z(b\gamma+G\beta+\beta\gamma+W_{rz})\\
 &+(2Vv+v^2+W_{\theta\theta})/R
   +\epsilon(\Delta_0-R^{-2})\beta,\\
E_\theta={}&t_*v+D_2(bv+\beta V+\beta v+W_{r\theta})
 +D_z(Gv+V\gamma+\gamma v+W_{\theta z})\\
 &-\epsilon(\Delta_0-R^{-2})v-D_2\Sigma_\theta,\\
E_z={}&t_*\gamma+D_1(b\gamma+\beta G+\beta\gamma+W_{rz})
 +D_z(2G\gamma+\gamma^2+W_{zz}+p_m)
 -\epsilon\Delta_0\gamma-D_1\Sigma_z.
\end{aligned}\tag{2}
\]

To derive (2), write each momentum transport as a cylindrical tensor
divergence, average in theta, and subtract the pure base balance. Exactly
one factor of w averages to zero. Two factors give W. The radial
centrifugal term supplies -(u_theta)^2/R; moving it to g_r gives the
displayed positive term. Cylindrical vector viscosity leaves -R^-2 in the
radial and angular components, with no such term in the axial component.
The mean's incompressibility is D_1 beta+D_z gamma=0. No assumption b=0 is
made; the inherited radial memory remains present.

The two required velocity constraints are

\[
M_\theta=\int R^2\bar v\,dR=0,\qquad
M_z=\int R\bar\gamma\,dR=0.\tag{3}
\]

All radial integrals keep physical (z,t) fixed. Since q=q(z,t) is independent
of radius, no new dyadic bands are introduced by these integrals.

## 2. A compact radial inverse must follow the phase path

For g=R^e f and e=0,1,2 set

\[
Ig(R,y)=\int_0^R g(R',y+M[(R')^{d_r}-R^{d_r}]v_r)\,dR',
\]
\[
Jg(R,y)=\int_0^\infty g(R',y+M[(R')^{d_r}-R^{d_r}]v_r)\,dR',
\quad I_cg=Ig-\chi_m Jg,
\]
\[
T_e f=R^{-e}I_c(R^ef),\qquad
A_e f=R^{-e}(\partial_R\chi_m)J(R^ef).
\tag{4}
\]

The fixed physical cutoff chi_m is zero near the inner shell edge and one
near the outer edge. It becomes a smooth moving slow cutoff in the chart.
Differentiating the integral along (1), not just at fixed y, gives exactly

\[
D_e T_e f=f-A_e f.\tag{5}
\]

I_c is zero on both exterior collars: below the source I=0, above it I=J.
It preserves the projection of support to (Z,T), although it need not
preserve auxiliary support. Omitting the shift in (4) adds a large fast
radial term and invalidates (5).

Here is why fixed coefficient derivatives cost no hidden powers of M.
Use U=R^(d_r) and extend F(U)=g(R)/(d_r R^(d_r-1)) smoothly by zero. Then

\[
A_-(U,y)=\int_{-\infty}^0 F(U+a,y+Ma v_r)\,da,\quad
A_+(U,y)=\int_0^\infty F(U+a,y+Ma v_r)\,da,
\]
\[
I_cg=(1-\chi_m)A_- -\chi_m A_+.\tag{6}
\]

At fixed a the shift is independent of U,Z,T. Derivatives hit F and chi_m,
and the zero extension creates no moving-boundary terms. On the positive
compact radial shell U and R have bounded coordinate changes. For an input
of class M_alpha, integration over bounded length preserves its epsilon
power and fixed-order polynomial S costs. Near an edge, only the primitive
toward that edge occurs. The exponential flat edge weight times any fixed
inverse-distance power is increasing sufficiently near that edge, so its
integral has the same type of bound. The output remains in M_alpha.

For the actual directions, put beta_g=sqrt(2)-1,
vr=(1,-beta_g), vt=(beta_g,1). Multiplying either divisor by its quadratic
conjugate gives a nonzero integer for any nonzero integer mode k. Hence

\[
|v_r\cdot k|,|v_t\cdot k|
 \ge [(1+\sqrt2)|k|_1]^{-1}.\tag{7}
\]

For example vr.k=(k1+k2)-sqrt(2)k2; the product with its conjugate is
(k1+k2)^2-2k2^2, which cannot be zero for nonzero integer k. The other
factor is at most (1+sqrt(2))|k|_1. The temporal case is the same argument.
Thus L^-p costs at most p Fourier derivatives; an additional three
derivatives leave a summable |k|^-3 tail in two dimensions.

If int R^e bar(f) dR=0, Haar translation invariance gives bar(A_e f)=0.
Integrate the full-line Fourier integral by parts p times in (6). Compact
support removes all endpoint terms, and every nonzero mode gains M^-p.
For each output derivative I and integer p>=1 this proves

\[
|\partial^I A_e f|\le C_{I,p}\epsilon^{\alpha+p\kappa}
 S^{B_{I,p}}\zeta\delta^{-B_{I,p}}.\tag{8}
\]

Only finitely many extra input derivatives are needed: p radial derivatives
and p+3 torus derivatives, in addition to the requested fixed derivatives.
The cutoff remainder is supported away from the radial edges. Given any
fixed physical derivative and desired q power, choose p large enough to
absorb its fixed graph/physical losses, then absorb the logarithmic S
power with a strict exponent margin. This proves flatness for a smooth
all-orders source; a finite regularity hypothesis would prove only a finite
gain. If f is auxiliary-independent and has zero weighted integral, A_e f
is exactly zero. A nonzero zero-mode moment instead leaves a nonzero
cutoff remainder: it cannot be called flat merely because f is smooth.

## 3. Pressure, divergence-free velocity, and all three defects

Choose the same fixed physical interior bump
rho_phys=q^(-1/2) rhohat(r/sqrt(q)), with int rhohat=1. In a chart,
rho=s^(-1/2)rhohat(R/sqrt(s)). Define

\[
P=\int\bar g_r\,dR,\qquad p_m=T_0(g_r-\rho P).
\tag{9}
\]

The input now has zero radial auxiliary-mean integral. Equations (5),(8)
give

\[
D_rp_m-g_r=-\rho P-A_0(g_r-\rho P),\qquad
\partial_R\bar p_m=\bar g_r-\rho P.\tag{10}
\]

The remainder is flat and has zero auxiliary mean; the compact radial
pressure obstruction is the scalar P. Pressure changes preserve the class
of changes in g_r because (9) is a fixed linear map.

For a desired axial increment gamma_d with int R bar(gamma_d)dR=0, set

\[
\Psi=T_1\gamma_d,\qquad
\Delta\beta=-D_z\Psi=-\epsilon\partial_Z\Psi,\qquad
\Delta\gamma=D_1\Psi=\gamma_d-A_1\gamma_d.\tag{11}
\]

This is a curl of an azimuthal vector potential: D_1 Delta beta+D_z Delta
gamma=0 exactly, since D_r and D_z commute. Its support is compact,
bar(Delta gamma)=bar(gamma_d), and its classes are M_(alpha+1),M_alpha.
For auxiliary-independent gamma_d of zero weighted integral the axial
increment equals gamma_d exactly. No slow derivative is commuted through
the cutoff inverse to assert these facts.

Define the remaining flux defects and bump moment

\[
J_\theta=\int R^2\overline{Gv+V\gamma+\gamma v+W_{\theta z}}\,dR,
\]
\[
J_z=\int R\overline{2G\gamma+\gamma^2+W_{zz}}\,dR
 -\tfrac12\int R^2\bar g_r\,dR,\qquad
c_\rho=\tfrac12\int R^2\rho\,dR.\tag{12}
\]

Multiply the two tangential equations in (2) by R^2 and R, respectively,
and integrate. Every radial flux is an endpoint term. The radial viscous
integrals also vanish, by two integrations by parts. In full generality,
before imposing (3), the identities are

\[
\int R^2\bar E_\theta\,dR
 =-\epsilon\partial_T M_\theta-\epsilon^3\partial_{ZZ}M_\theta
  +\epsilon\partial_ZJ_\theta,
\]
\[
\int R\bar E_z\,dR
 =-\epsilon\partial_T M_z-\epsilon^3\partial_{ZZ}M_z
  +\epsilon\partial_Z(J_z+c_\rho P).\tag{13}
\]

In particular (3) removes the first two terms of each row. To derive the
pressure contribution, integrate (10) once more with radial weight:

\[
\int R\bar p_m\,dR
 =-\tfrac12\int R^2\bar g_r\,dR+c_\rho P.\tag{14}
\]

The WHOLE product c_rho P stays inside the axial derivative. Indeed
c_rho=s/2 int x^2 rhohat(x)dx changes with (Z,T). Keeping c_rho partial_Z P
while dropping (partial_Z c_rho)P would be wrong even for perfectly smooth
waves and pressure.

With normalized interior bumps int R^2 sigma_theta=int R sigma_z=1, compact
covariance targets are therefore

\[
H_\theta=-T_2[\bar E_\theta-\sigma_\theta\epsilon\partial_ZJ_\theta],
\quad H_z=-T_1[\bar E_z-\sigma_z\epsilon\partial_Z(J_z+c_\rho P)].
\tag{15}
\]

Their sources are auxiliary-independent and have zero weighted integrals,
so the inverse is exact. They are signed stress targets for the fixed
primary-wave inverse; (15) alone does not realize a new wave covariance or
remove its quadratic error.

## 4. Fast mean update and the fixed five-equation velocity inverse

For the complementary zero-auxiliary-mean tangential residuals use

\[
\Delta v=-c_i^{-1}N_{time}^{-1}E_\theta^\circ,\qquad
\gamma_d=-c_i^{-1}N_{time}^{-1}E_z^\circ.\tag{16}
\]

The Fourier multiplier is (2 pi i vt.k)^-1 on nonzero modes and zero at
mode zero. From (7), output order m requires at most m+4 torus derivatives.
Together with c_i^-1<=CS this preserves M_alpha and reality. It preserves
slow/radial support and its smooth zero extension, but can fill the
auxiliary torus. Realize gamma_d by (11); the axial fast-time equation is

\[
c_iN_{time}\Delta\gamma=-E_z^\circ+c_iN_{time}(\Delta\gamma-\gamma_d).
\tag{17}
\]

The last term is flat by (8), using finitely more derivatives. It must
remain in the residual and the actual velocity must remain in all later
products. Both velocity moments (3) are preserved.

These definitions agree on chart overlaps. In physical variables the radial
shift in (4) is Y+[(r')^(d_r)-r^(d_r)]vr, with cutoff chi_m(r^2/(2q)).
Rescaling r=sqrt(Q)R and pulling back by Jg^i produces exactly (4) and M.
For time, let Pi_i F(Y)=F(Jg^i Y). The integer covering preserves Haar
average, and the eigenvector relation gives
Nabs Pi_i=Tg^i Pi_i N_i, hence
Nabs^-1 Pi_i=Tg^(-i)Pi_i N_i^-1 on zero-mean functions.
Define the physical temporal update as -Nabs^-1 E_phys^circ. Converting
residual to velocity normalization gives Q^(-1-h)Tg^(-i)=c_i^-1, exactly
as in (16). Thus the chart formulas represent one physical map. A common
covering may be used without averaging distinct preimages on a finer band
torus.

The three current defects are corrected next by slow velocities on the
reserved patch. There the SAME summed base has G_q=0 and
V_q=a(eta)x^(-1-2lambda), |a|>=a0>0, x=r/sqrt(q). This does not set the
radial base b to zero. Choose a nonnegative nonzero smooth bump eta0 in
that patch and disjoint geometric copies eta_j(x)=a_j^-1 eta0(x/a_j),
a_j=exp(jd), j=0,1,2. The support widths and d are fixed once.

Their power moments are mu_p exp(jdp), mu_p=int x^p eta0>0. The three
angular rows have powers (2,-2-2lambda,-2lambda); the two axial rows have
powers (1,1-2lambda). The matrices are

\[
A_\theta=[\,c_i\mu_{p_i}e^{jdp_i}\,]_{i,j=0}^2,
\quad c=(1,2a,-a),\qquad
A_z=[\,d_i\mu_{q_i}e^{jdq_i}\,]_{i,j=0}^1,
\quad d=(1,a).\tag{18}
\]

The determinant of each row-normalized matrix is the Vandermonde product
in exp(dp_i) or exp(dq_i). The powers are distinct for the fixed lambda>0.
This gives an explicit inverse with no small-target restriction:

\[
u=A_\theta^{-1}(0,-P_q,-J_{z,q})^T,\qquad
s=A_z^{-1}(0,-J_{\theta,q})^T.\tag{19}
\]

Use u's three coefficients for Delta v and s's two for gamma_d, with
physical velocity factor q^-A. Here each target uses its OWN physical
normalization: P_q=q^(2A)P_phys, Jtheta_q=q^(2A-3/2)Jtheta_phys and
Jz_q=q^(2A-1)Jz_phys. Multiplying each physical row and its target by the
same factor gives the chart equations

\[
\int R^2\Delta v=0,\quad\int R\gamma_d=0,
\quad\int(2V/R)\Delta v=-P,
\]
\[
\int R^2(G\Delta v+V\gamma_d)=-J_\theta,\qquad
\int(2RG\gamma_d-RV\Delta v)=-J_z.\tag{20}
\]

Because gamma_d is slow and its weighted integral is zero, (11) realizes
Delta gamma=gamma_d exactly. The potential and induced radial velocity
stay in the convex radial hull of the chosen bump supports, within the
reserved patch. Targets in S_alpha give Delta v,Delta gamma in M_alpha and
Delta beta in M_(alpha+1), for every fixed derivative. The inverse depends
smoothly on a and its derivatives, with |a|>=a0, on a compact parameter
set. Adjugates divided by the displayed nonzero determinants supply finite
bounds. The bounds may be VERY large for a small fixed lambda; no uniform
limit lambda->0 is claimed. At lambda=0 the axial powers coincide and the
map fails. No extra q-smallness is required to define this fixed linear
map for a finite target, so it does not shrink the pulse domain by stage.

## 5. Recompute from the same velocity: exact remainder and a mean-only gain

Hold the complete w fixed during the slow update (20). Write increments
(d_beta,d_v,d_gamma) and let g_r,new be recomputed from (2). Then

\[
R_g=g_{r,new}-g_r-(2V/R)d_v
\]
\[
\begin{aligned}
={}&-t_*d_\beta-D_1(2bd_\beta+2\beta d_\beta+d_\beta^2)\\
 &-D_z(bd_\gamma+Gd_\beta+\beta d_\gamma+\gamma d_\beta+d_\beta d_\gamma)\\
 &+(2v d_v+d_v^2)/R+\epsilon(\Delta_0-R^{-2})d_\beta.
\end{aligned}\tag{21}
\]

W cancels here because the wave is held fixed, rather than discarded from
either radial source. Expand (12), use all five rows (20), and obtain

\[
P_{new}=\int\overline{R_g}\,dR,
\quad J_{\theta,new}=\int R^2\overline{\gamma d_v+v d_\gamma+d_\gamma d_v}\,dR,
\]
\[
J_{z,new}=\int R\overline{2\gamma d_\gamma+d_\gamma^2}\,dR
 -\tfrac12\int R^2\overline{R_g}\,dR.\tag{22}
\]

These formulas retain the new velocity square. Solving the linear five-map
does NOT set the actual nonlinear defects to zero.

Assume the current means have v,gamma in M_0.9 and beta in M_1.9, and the
three targets are in S_alpha, alpha>=0.9. The inherited actual base has
b=O(epsilon), V,G=O(1) with the same estimates for every fixed coefficient
derivative on the potential support. The slow increments have orders
alpha,alpha,alpha+1. A conservative term ledger for (21) is:

| Term | Lower epsilon order |
|---|---:|
| slow time of d_beta | alpha+2 |
| radial derivative of b d_beta | alpha+2-kappa |
| radial products with old beta | alpha+2.9-kappa |
| radial square of d_beta | 2alpha+2-kappa |
| axial base fluxes | alpha+2 |
| axial products with old means | alpha+2.9 |
| axial new product | 2alpha+2 |
| old swirl times d_v | alpha+0.9 |
| square of d_v | 2alpha |
| radial viscosity of d_beta | alpha+2-2kappa |
| axial viscosity of d_beta | alpha+4 |

The fast-time derivative of d_beta is zero because this update is slow.
Radial graph differentiation costs at most kappa each time. Products
retain at least the required flat edge weight; the moment-update supports
are inside the shell. All rows, and all products in (22), have order at
least alpha+0.9-2kappa. Thus, under precisely these hypotheses,

\[
R_g\in M_{\alpha+0.9-2\kappa},\qquad
(P_{new},J_{\theta,new},J_{z,new})\in S_{\alpha+0.9-2\kappa}.
\tag{23}
\]

This is a positive **compatibility-defect gain**, including radial memory,
pressure's axial moment and viscosity. It is not a gain for nonzero angular
harmonics or for the whole momentum residual. The temporal update (16)
likewise has conditional tangential linear remainder order
alpha+1-2kappa after pressure is recomputed, from slow time, radial/axial
transport and viscosity. Its wave interactions still require the full
cycle estimate below.

## 6. Reproducible checks and the next inequality

The checks compare the full physical angular mean of momentum with (2),
using the previous compact exact-curl fixture with an added compact mean.
A separate normalized-coordinate calculation tests epsilon factors and
viscosity. All covariance entries, the induced radial velocity, old radial
memory, and the exact increment remainder are used.

The integrated-pressure diagnostic uses the complete curl covariance as
manufactured smooth tensor data. It is not a normalized realization of
the actual band flow. It integrates in radial subintervals split at EVERY
support edge; a coarse unsplit rule fails on the narrow bump and steep
cutoff derivatives. The Z derivative of g_r uses a recorded five-point
stencil, since differentiating an order-three Taylor residual a fourth
time would be invalid. A tighter stencil and quadrature must reproduce
the identity. These are approximate diagnostics, not interval enclosures.

The radial Fourier fixture uses d_r=2 and a degree-16 compact polynomial
in U. It has C7 zero extension, sufficient for the finitely many tested
identities, and is explicitly NOT the all-orders source assumed in (8).
Its oscillatory integrals use an exact finite integration-by-parts formula,
avoiding an unresolved huge-frequency mesh. Actual torus eigendirections,
signed bump matrices and their lambda=0 failure are checked separately.
Every source from prior checkpoints remains frozen.

```sh
python3 research/navier_stokes_cascade/results/compact_mean_corrections/test_audit.py -v
python3 research/navier_stokes_cascade/results/compact_mean_corrections/audit.py \
  --verify-record research/navier_stokes_cascade/results/compact_mean_corrections/evidence.json \
  --output /tmp/compact_mean_corrections_recheck.json
```

The next gate is the SAME-field complete residual cycle. For example,
with wave/mean/defect orders B_j, C_j=B_j+1/2, prove

\[
G_{wave}^{j+1}\in W_{B_j+1/10},\quad
(E_\theta^{j+1},E_z^{j+1})\in M_{C_j+1/10},\quad
(P^{j+1},J_\theta^{j+1},J_z^{j+1})\in S_{C_j+1/10}.
\tag{24}
\]

This must follow after, in order: sourced wave curls; the signed stress
increment from (15); the auxiliary-time update (16); the five-map (20);
and pressure recomputation after EVERY velocity change. The current
checkpoint proves (23) for the last map's stated input, not all three
inequalities (24) or their initialization. In particular the source/support
and Gaussian-envelope bounds for nonzero harmonics must survive mean
operations that fill the auxiliary torus. All old/new wave products, the
signed increment's quadratic covariance, mean/wave transport, flat cutoff
terms multiplied by later fields, and viscous curl terms must be included.
Failure to prove (24) is the present missing estimate, not a demonstrated
impossibility of the manuscript's construction. Converting a positive
normalized gain into smooth forcing still requires iteration, physical
fixed-derivative bounds, summation and localization.
