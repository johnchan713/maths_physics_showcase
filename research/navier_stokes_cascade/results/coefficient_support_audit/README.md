# Formal coefficient supports and the reserved mean interval

Status: `formal-coefficient-supports-and-exact-reserved-mean-velocity-audited`.

This checkpoint closes the common-support and untouched-mean-VELOCITY
premises of the formal coefficient construction, using the fixed leading
profile and the preceding analytic collar and linear inverse lemmas. It
does not establish the outstanding all-order weighted leading edge
certificate, uniform physical pulse bounds, complete mean corrections or
the full residual. In particular, pressure is NOT equal to leading pressure
on the mean patch: the first two pressure orders must be retained.

Parent: `3adeeaa46c64ee2dc6d47dcaf951baaf946cfc05`. The supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`,
(5.2)--(5.8), (5.27), (5.42)--(5.46), (A.9), supplies defining equations,
not an assumed blow-up theorem. PDF pages 46 and 52 were visually checked
for the shifted viscosity operator, pressure moment, and velocity support
claim. Historical source and evidence are unchanged.

## 1. The actual remaining patches are disjoint

Keep the selected Md64 family, T=exp(64)+10, lambda=exp(-4T), h=lambda^2,
and N=1+floor(Hhier^32). Write D=1/2-h, A=1/2+h,
d=1-eta^2, L=1-2h eta^2 and lambda_n=2nh.
The power stage begins at Xw=XR exp(T+2), with duration Tw=240T.
Its reserved intervals in y_stage=log(X/Xw) are:

| Interval | Offsets from Tw | Use |
|---|---|---|
| I1 | (-25,-20) | Leading modulation repair |
| I2 | (-20,-15) | Leading heat compensation |
| I3 | (-14,-9) | Positive-order moments |
| I4 | (-8,-3) | Reserved mean correction |

The [regional cone join](../regional_cone_join/README.md) and exact
[global ledger](../global_moment_stress/README.md) leave I3 and I4 unchanged.
With rho3=XR exp(241T-12), I3 is 0<log(X/rho3)<5 and I4 is
6<log(X/rho3)<11. The two U bumps have centers 1,2 and the three E
bumps have centers .5,1.5,2.5, all of width .1. Set

    Xvel = rho3 exp(3).

This is only a support envelope, NOT a modification of the construction.
All five bump supports end by 2.55<3<5<6. The axial pulse is later than
I4; Xv denotes its completed support edge, before the outer collar.
The inner stress-free radius is Xinner=Xa exp(tc/16); the coefficient
cutoff remains one through Xa exp(tc/8), vanishes from Xa exp(3tc/16),
and lies inside the analytic collar Xa exp(tc/4), before modulation
starts at Xa exp(tc/2). These are exact relative offsets, not rounded radii.

## 2. Induction constructs every finite formal coefficient

Use the regular six-variable Volterra system and sparse analytic-parameter
lemma in [background jet foundation](../background_jet_foundation/README.md),
sections 2--4. On the SAME physical radial collar, the coefficients and
source have finite bounds C_n for each finite n. The complex eta strip may
shrink with n; no uniform C_n is asserted. The derivative blocks map the
first four components into the last two and annihilate the last two, even
with intervening diagonal kernels. The factorial Picard estimate therefore
converges for any finite C_n. Leading modulation vanishes on this collar;
the leading core/collar analytic premises are the replayed defining-map
lemmas, not a conclusion from numerical grids or finite C3 bounds.

Finalize all lower orders globally before forming the next source. Cut
only F_n and U_n with the fixed inner cutoff; add the five I3 bumps, and
reconstruct FROM ZERO AXIS VALUES

\[
 M_n=\int_0^X U_n\,dx,\quad E_n=\sqrt{2X}F_n,\quad
 V_n=\frac{2\eta XU_n-2(D+\lambda_n)\eta M_n-d(M_n)_\eta}{L},
\]
\[
 (\Pi_n)_X=\sum_{i+j=n}F_iF_j-\frac{\Omega_{n-1}}{2X}.
 \tag{1}
\]

The regular axis quotients in the inherited lemma make (1) smooth.
The five total moments are

\[
 \int U_n,\quad\int\sqrt{2X}E_n,\quad\int(\Pi_n)_X,\quad
 \int\sqrt{2X}\sum_{i+j=n}U_iE_j,\quad
 \int\left[\sum_{i+j=n}U_iU_j-X(\Pi_n)_X\right].\tag{2}
\]

All incoming integrals include the FULL lower-order pressure source through
Xv, not only I3. At current order, products involving n have only pairs
(0,n),(n,0), while Omega_(n-1) is already fixed; squares have order 2n>n.
On I3, U0=0 and E0=e(eta) exp[-(.5+lambda)y], with smooth e>0.
Thus the exact variations of (2) are

\[
 \int\delta U_n,\quad\int\sqrt{2X}\delta E_n,\quad
 \int E_0\delta E_n/X,\quad\int\sqrt{2X}E_0\delta U_n,
 \quad-\int E_0\delta E_n.\tag{3}
\]

Normalize (2) by (rho3,sqrt(2)rho3^(3/2),e,
sqrt(2)rho3^(3/2)e,rho3 e), call the results (M,I,Cp,J,S), and use rows
(M,(M-J)/lambda,I,-S,Cp). The continuum matrix B is exactly the inherited
fixed bump matrix. Its normalized infinity inverse norm is <1000
(U block <2, E block <787.562). The interval calculation reproduces the
analytic determinant lower bounds, not a sampled matrix invertibility test.
The raw costs lambda^-1 and e(eta)^-1 remain; they can be enormous but
are finite, smooth and bounded for this fixed compact parameter interval.
There is no small-target requirement: c=-B^-1 z_base is an EXACT linear
solution of all five equations. Consequently M_n=V_n=0 after Xvel,
as well as F_n=U_n=0. This proves the induction support claim for EVERY
finite order, not just the orders exercised by diagnostics.

## 3. Pressure support has two exceptional orders

For k>=1 the radial source Omega_k includes

\[
 -\nu Z_{\lambda_{k-1}-D}Z_{\lambda_{k-1}}V_{k-1},\qquad
 Z_a=\frac{2a\eta-2\eta X\partial_X+d\partial_\eta}{L}.\tag{4}
\]

All its other terms are linear in the kth velocity or quadratic products
of total order k. On X>=Xvel every positive-order velocity and its jets
vanish. Hence Omega_1=-nu Z_-D Z_0 V0, whereas Omega_k=0 for k>=2.
Equation (1) gives, throughout that exterior,

\[
 (\Pi_1)_X=-\Omega_0/(2X),\qquad
 (\Pi_2)_X=\nu Z_{-D}Z_0V_0/(2X),\qquad
 (\Pi_n)_X=0\quad(n\ge3).\tag{5}
\]

These first two sources vanish after Xv, where leading U0=M0=V0=0.
The exact THIRD TOTAL moment and zero axis pressure imply backward
pressure primitives. Therefore Pi1,Pi2 vanish after Xv, and Pi_n for
n>=3 vanish already after Xvel. Neither zero velocity nor zero local
pressure derivative is enough by itself.

On the actual I4, matching of the leading accumulated moment gives
M0=mu eta, mu>0, U0=0, V0=-mu and Omega0=-mu^2/(2X).
Thus, on I4,

\[
 \boxed{\Pi_1=K_1(\eta)-\frac{\mu^2}{4X},\quad
 \Pi_2=K_2(\eta),\quad\Pi_n=0\ (n\ge3).}\tag{6}
\]

The order-one radial derivative is strictly positive. The order-two radial
derivative is zero, but K2 depends on the complete later axial source and must NOT
be set to zero. A labelled manufactured control makes the distinction:
take M0=eta B(X), with B=mu on a plateau and smoothly returning to zero
later. At eta=0, Z_-D Z_0 V0=6X B', so Pi2_X=3nu B' and its backward
plateau value is 3nu mu, nonzero. This value is NOT asserted for the actual
leading profile, whose later axial pulse has its own angular dependence.

The exact conservative total-source identities from the background
summation turn the last moment equations in (2) and preceding-order
viscosity moments into zero total stress primitives. The inner equations
give T_n=0 through Xinner. In the exterior, use backward primitives:
order one contains the leading axial-viscous heat stress, vanishing only
at Xb; order two may contain the pressure tail through Xv; orders n>=3
have all same-order sources and shifted positive velocity viscosity zero
after Xvel. Thus the common annular support and the sharper right edges are:

| Positive order | Velocity right edge | Pressure right edge | Stress right edge |
|---|---|---|---|
| 1 | Xvel | Xv | Xb |
| 2 | Xvel | Xv | Xv |
| n>=3 | Xvel | Xvel | Xvel |

The [outer edge stress](../outer_edge_stability/README.md) identification
now has its common-support premise. This does not upgrade it to an
all-order weighted certificate at the inner flat edge.

## 4. The summed velocity is exactly leading on I4; pressure is not

For S_n=q^(D+lambda_n)M_n, multiply the STREAMFUNCTION, not each velocity,
by chi(c_n q). Its curl product rule is

\[
 r u_{rn}^{cut}=q^{\lambda_n}
 \left[\chi V_n-\frac{2\eta}{L}(c_nq)\chi'(c_nq)M_n\right].\tag{7}
\]

The axial and angular coefficients simply acquire chi. Since M_n and
ALL its spatial jets vanish on I4, both terms in (7) vanish there.
Every locally finite potential-cutoff sum therefore agrees EXACTLY in
velocity with the leading field on I4, at every q>0 in its domain,
including cutoff transitions. This identity requires no uniform-in-n
estimate and does not establish the weighted cutoff-selection inputs.
In normalized variables x=r/sqrt(q), the mean patch has

    q^A uz=0, q^A utheta=a(eta) x^(-1-2lambda),
    q^A ur=-mu q^h/x, a(eta)>0.

In particular, leading radial velocity is not zero. Meanwhile pressure is

\[
 q^{2A}p_B=\Pi_0+\chi(c_1q)q^{2h}
 [K_1-\mu^2/(4X)]+\chi(c_2q)q^{4h}K_2.\tag{8}
\]

No higher pressure coefficient survives on I4. The velocity premise for
the reserved mean matrix is supplied, but a complete covariance mean
update must retain (8), its axial derivatives, leading radial velocity,
and the full-wave axial velocity terms. None is removed by support.

## 5. Checks, interpretation and next obligation

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output NEW.json` here.
Dependencies are in requirements.txt. The replay checks the parent's
record and 551 inherited source hashes. New checks include 90 exact rational
velocity/cutoff identities, all-order algebra with finite grading controls,
outward determinant bounds, and the independent full-series radial
viscosity assembly. Manufactured diagnostics use h=.005, not the actual
astronomically small h, and do not prove continuum existence or evaluate
the actual coefficients. The analytical induction above carries that claim
relative to the explicitly cited regularity and inverse lemmas.

Next: verify all-order weighted leading stress/first-stress edges needed by
the background jet certificate, then establish actual-pulse size,
derivative and support estimates at every scale. Complete covariance means
must use (8). The final residual still requires viscosity, all cutoffs and
all interactions. No unconditional pulse estimate, residual improvement
or blow-up conclusion is claimed.
