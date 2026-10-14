# First global background coefficient, with radial memory retained

Status: `conditional-first-global-background-coefficient-and-relative-gain`.

This checkpoint conditionally extends the first coefficient across the actual
annulus, restores its five compatibility moments, and derives its compact
signed tangential stress. It corrects an overstrong premise in the previous
work plan: **on I3, U0=0 does not imply M0=V0=Omega0=0**. That premise is
unnecessary for the linear positive-order moment solve. Its omission would,
however, give the wrong pressure and fifth incoming moment.

The argument depends on the inherited leading profile, its analytic-parameter
collar, exact moments, heat exterior, and continuum five-moment inverse.
It is not an independent verification of those foundations. The coefficient
is specified by convergent integrals and a finite linear solve, not by a
numerical evaluation of the enormous fixed-frequency profile. Manufactured
tests check the identities and reject missing terms; they do not certify
continuum existence.

Parent: `4a5a098862090231505d1976378574d6034ee9cb`. Source: the supplied
[Navier-Stokes manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Theorem 4.6, (5.2)-(5.6), Lemma 5.2, (A.9), and Corollary B.6; SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The source's proposed blow-up theorem is not treated as an established result.
No earlier audit or GOAL.md is changed. The same leading profile and
`N=1+floor(H^32)` are retained.

## 1. Continue into the collar BEFORE cutting off

The preceding inner coefficient exists on `X<=3/Lambda=3Xa/4`. A cutoff
there would create an error inside the stress-free core, so it is not an
admissible annular extension. Instead use the inherited common analytic-eta
rectangle, preserved through `Xa exp(tc/4)`, before modulation begins at
`Xa exp(tc/2)`. Choose radii by their exact logarithmic offsets:

| Purpose | Radius |
|---|---|
| Guaranteed stress-free coefficient interval | `Xminus=Xa exp(tc/16)` |
| Cutoff remains one through | `Xkeep=Xa exp(tc/8)` |
| Cutoff becomes zero at | `Xcut=Xa exp(3tc/16)` |
| End of coefficient's analytic collar | `a^2=Xa exp(tc/4)` |

All these offsets are positive and strictly ordered. Do not add the
astronomically small offsets to a rounded log(Xa); retain the common base
and rational multiples of tc separately.

The first six-variable system from the
[inner checkpoint](../first_correction_construction/README.md) is linear for
any smooth leading profile; it does not require zero leading stress in the
collar. Its coefficients on this rectangle are smooth in the radial variable
and holomorphic in eta. The inherited forward collar formulas have analytic
eta inputs and eta-independent radial cutoffs; a smaller common complex
neighborhood excludes their denominators' zeros. Every fixed radial
derivative has finite bounds on a smaller neighborhood. This is precisely
the regularity used here; radial analyticity across the flat activation is
NOT asserted.

On such a neighborhood define the finite actual supremum

\[
C_{\rm col}=\max(1,\sup\|A_0\|,\sup\|A_1\|,\sup|f_1|).
\]

This is an existence bound for the fixed inherited construction, **not** the
earlier core-only bound C^16 and not a uniform bound as parameters vary.
The same sparse derivative block annihilates adjacent derivative factors.
For any smaller strip loss `0<Delta<=1`, the Picard sum is bounded by

\[
(z+\beta)e^\beta,\qquad z=2C_{\rm col}a,\qquad
\beta=e z^2/\Delta.
\]

Thus the same proof constructs W1 on the full collar, however large its
finite coefficient bound is. Analytic uniqueness identifies its restriction
with the already constructed core solution. Smooth radial regularity follows
from the Volterra equation and its differentiated form on smaller eta
neighborhoods; at the axis it agrees with the earlier regular solution.
This extends the actual conditional coefficient rather than choosing a new
core solution.

Let kappa be a smooth radial cutoff equal to one through Xkeep and zero
from Xcut onward. Set tentative `F1=kappa F1_in`, `U1=kappa U1_in`, extending
by zero, and later add the five bumps. Do NOT cut off V1 or Pi1 directly.

## 2. The actual I3 geometry and the missing memory

The intermediate power stage begins at `Xw=XR exp(Td+2)`. With local
`y_stage=log(X/Xw)`, (A.9) reserves

\[
I_3:\quad T_w-14<y_{\rm stage}<T_w-9,\qquad
I_4:\quad T_w-8<y_{\rm stage}<T_w-3.
\]

Consequently use `rho3=XR exp(Td+Tw-12)` and `X=rho3 exp(y)`, `0<y<5`,
for the positive-order patch. The same five bump supports, between .45 and
2.55, lie strictly in I3. The leading cone repair occupies I1, the heat
compensation occupies I2, and both restore the relevant accumulated moments.
I3 and the later mean patch I4 remain unedited. Here

\[
U_0=0,\qquad E_0=e(\eta)e^{-(1/2+\ell)y},\qquad
e(\eta)>0,\quad \ell=\lambda_{\rm patch}>0.
\]

The leading accumulated M is the matched reference M, not zero by a support
argument. In the reference stages preceding the axial pulse, U is eta times
a nonnegative radial function, equal to 4 eta on its inner ideal region.
Therefore, on I3,

\[
M_0=\mu\eta,\quad \mu>0,\qquad
V_0=\frac{-2D\eta M_0-d(M_0)_\eta}{L}=-\mu,
\qquad \Omega_0=-\frac{\mu^2}{2X}.
\]

We used `2D eta^2+d=L`; V0 is constant in X and eta on this patch.
The axial pulse occurs later and only then restores the leading M to zero.
Even at eta=0, where M0 itself is zero, its angular derivative is not zero.
Thus the radial pressure source on I3 is `-Omega0/(2X)=mu^2/(4X^2)`.
The matching inverse below needs only U0=0 and the stated swirl, not zero
M0 or zero Omega0. More generally it remains valid for any fixed smooth
lower-order source there.

## 3. Define the five actual incoming discrepancies and solve them linearly

For every tentative extension reconstruct, from ZERO axis values,

\[
M_1=\int_0^X U_1(x,\eta)\,dx,\qquad
V_1=\frac{2\eta XU_1-2(D+2h)\eta M_1-d(M_1)_\eta}{L},
\]
\[
\Pi_1=\int_0^X\left(2F_0F_1-\frac{\Omega_0}{2x}\right)dx,
\qquad E_1=\sqrt{2X}F_1.
\]

The integrand Omega0/x is smooth at the axis. It is zero beyond the leading
axial support edge Xv, because U0=V0=0 there. All the following ordinary
integrals therefore converge and are smooth functions of eta:

\[
\begin{aligned}
m_1&=\int U_1,&m_2&=\int\sqrt{2X}E_1,&m_3&=\int(\Pi_1)_X,\\
m_4&=\int\sqrt{2X}(U_0E_1+U_1E_0),&&
m_5&=\int[2U_0U_1-X(\Pi_1)_X].
\end{aligned}
\]

In particular the baseline third and fifth moments contain respectively
`-integral Omega0/(2X)` and `+integral Omega0/2`. Integrate these through
the full leading source, including the later axial pulse, not only through
I3 or the inner cutoff. Denote these tentative values by m^base(eta).

Add two U1 bumps and three E1 bumps on the same disjoint log supports used
by the inherited matrix B. Lower-order Omega0 is unchanged. On I3 their
moment variations are exactly

\[
\delta m_1=\int\delta U_1,\quad
\delta m_2=\int\sqrt{2X}\delta E_1,\quad
\delta m_3=\int E_0\delta E_1/X,
\]
\[
\delta m_4=\int\sqrt{2X}E_0\delta U_1,\qquad
\delta m_5=-\int E_0\delta E_1.
\]

There are no quadratic current-order terms. Normalize `(m1,m2,m3,m4,m5)` by
`(rho3,sqrt(2)*rho3^(3/2),e,sqrt(2)*rho3^(3/2)*e,rho3*e)`. Call these
`(M,I,Cp,J,S)` and form `z=(M,(M-J)/ell,I,-S,Cp)`. The matrix of the
variations is exactly the inherited B, with U weights
`exp(y)`, `exp(y)*(1-exp(-ell*y))/ell`, and E weights
`exp(3y/2)`, `exp((1/2-ell)y)`, `exp((-1/2-ell)y)`.

Define the actual coefficient functions by

\[
c(\eta)=-B^{-1}z^{\rm base}(\eta).
\]

The inherited continuum inverse has norm at most 1000 in these normalized
coordinates. This specifies the actual finite solve, even when discrepancies
are enormous. Smoothness follows from the finite integrals, the positive
smooth e, and the constant invertible B. The ell^-1 cost in arbitrary data
has not disappeared. Reusing the nonlinear LEADING repair would introduce
spurious terms from the next expansion order.

## 4. Support, pressure tail, and preserved inner regularity

All five moments now vanish pointwise in eta. Above the velocity bumps,
F1=U1=0 and m1=0 gives M1=0, hence V1=0. They therefore vanish on I4.
Pressure need not vanish there: outside the bumps it obeys

\[
(\Pi_1)_X=-\Omega_0/(2X),\qquad
\Pi_1(X)=\int_X^{X_v}\frac{\Omega_0(x)}{2x}\,dx.
\]

The third moment and the zero axis datum imply Pi1=0 beyond Xv. Thus the
ENTIRE first coefficient, including pressure, has support in `[0,Xv]`,
not necessarily in the smaller velocity support. Directly cutting Pi1 off
at the end of I3 violates the radial momentum equation.

Forward integration and the late bumps leave the coefficient unchanged
through Xkeep. The radial cutoff is eta-independent, so analytic-eta
regularity survives throughout the chosen inner rectangle, including its
cutoff transition. M1/X is smooth at zero by the identity
`M1/X=integral_0^1 U1(tX,eta) dt`; V1/X is then smooth as well. All fixed
global profile derivatives have finite constants for the fixed construction.
No effective universal bound on those constants is asserted.

## 5. Close the two signed stress tails

Let rtheta1 and rz1 be the first tangential residual coefficients, with the
physical power `q^(-A-1+2h)` removed. Define

\[
T_{1,\theta}=-R^{-2}\int_0^R\varrho^2r_{\theta1}(\varrho)d\varrho,
\qquad T_{1,z}=-R^{-1}\int_0^R\varrho r_{z1}(\varrho)d\varrho.
\]

The conservative momentum equations give, with the inherited zero lower
viscous moments,

\[
\int R^2r_{\theta1}dR
=\frac{-(1+h)m_2+D\eta(m_2)_\eta+\eta m_4+d(m_4)_\eta}{L}=0,
\]
\[
\int Rr_{z1}dR
=\frac{-A m_1+D\eta(m_1)_\eta+d(m_5)_\eta}{L}=0.
\]

The leading angular viscous term uses the convergent SUBTRACTED exterior
power moment supplied by the heat audit. It is invalid to separate this
into two divergent integrals and differentiate them independently. The
pressure integration by parts uses compact Pi1 to give
`integral Pi1=-integral X Pi1_X`, which is why the fifth moment has its
specific pressure term.

The zero totals give equivalent backward stress integrals. The inner
coefficient equations make T1=0 through Xkeep. Outside Xv only the leading
azimuthal axial diffusion remains. It vanishes beyond Xb because the
physical leading heat field is z-independent. Hence
`supp T1 subset [Xminus,Xb]`. On the outer collar its flat cutoff gives,
for every fixed derivative, `|partial^I T1|<=C_I zeta delta^(-N_I)` by
the inherited backward heat-edge integral estimate. Interior compactness
and the identically zero inner portion complete this bound. T1 is signed;
this two-component stress is NOT a full positive wave covariance tensor.

## 6. Exact first-truncation remainder: relative gain, still singular

For the actual extended first coefficient, remove the prescribed tangential
stress `T0+epsilon T1`, where `epsilon=q^(2h)`. The radial normalized
coefficients of orders zero and one vanish globally by the pressure
identities. The first two tangential orders are canceled by their retained
stress. The remaining normalized coefficients are exactly

\[
\begin{aligned}
B_r&=\Omega_1-2XF_1^2-\nu Z_{-D}Z_0V_0,\\
C_r&=\Omega_2-\nu Z_{2h-D}Z_{2h}V_1,\\
B_\theta&=V_1(F_{1,X}+F_1/X)+U_1Z_{b+2h}F_1
-\nu Z_{b+2h-D}Z_{b+2h}F_1,\\
B_z&=V_1U_{1,X}+U_1Z_{c+2h}U_1
-\nu Z_{c+2h-D}Z_{c+2h}U_1.
\end{aligned}
\]

Here Omega1 is the linear cross coefficient including the time derivative
and radial viscosity of V1, and Omega2 is its quadratic self-interaction.
The physical remainder of this augmented BACKGROUND is

\[
\mathcal E_r=q^{-3/2+2h}(B_r+\epsilon C_r)/R,\quad
\mathcal E_\theta=R q^{-3/2+3h}B_\theta,\quad
\mathcal E_z=q^{-3/2+3h}B_z.
\]

The axis expression is interpreted in regular Cartesian coordinates, not
by an annular inverse-radius estimate. On every fixed compact profile range
there are finite constants C_m such that, for `0<q<=1`, all Cartesian
space-time derivatives through m satisfy the sufficient bound

\[
|\mathcal E|_m\le C_m q^{-3/2+2h-m}.
\]

Each transverse derivative loses at most q^-1/2, each axial derivative
q^-D, and each time derivative q^-1; the displayed coarser loss q^-m
covers their sums. This is one extra relative factor q^(2h) over the
leading augmented bound. For h=.005 its undifferentiated exponent is
still -1.49. **It does not prove absolute residual decay, let alone a
smooth force.** It also does not bound an actual wave realization, its
full tensor divergence, cross interactions, or localization errors.

## Reproducibility and next work

The checks independently evaluate 12 manufactured global moment repairs,
15 radial-source comparisons, 18 conservative-moment identities and 18
finite-remainder identities. The finite fixture has C4 leading memory and
C5 inner cutoff fields; it is not used to infer the C-infinity continuum
argument. Removing Omega changes the third moment off the symmetry plane
and the fifth moment everywhere tested. At eta=0 its TOTAL third-moment
contribution cancels, so that one scalar test alone cannot detect the error;
the nonzero local pressure tail and fifth moment still do. No threshold
was weakened to hide that cancellation.

```sh
python research/navier_stokes_cascade/results/annular_first_correction/test_audit.py -v
python research/navier_stokes_cascade/results/annular_first_correction/audit.py \
  --verify-record research/navier_stokes_cascade/results/annular_first_correction/evidence.json \
  --output annular_first_correction_reproduced.json
```

Next: extend the background induction to higher orders while preserving
this same collar and all five moments, control the fixed-derivative constants
and divergence-preserving summation, and construct the physical wave and mean
corrections that actually generate the retained annular stress. The next
full-error milestone must include those waves' complete covariance and
interactions. Smooth forcing, complete-field energy and blow-up remain open.
