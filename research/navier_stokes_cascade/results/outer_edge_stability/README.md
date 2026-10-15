# Direct outer-edge certificate and the first viscous stress

Status: `outer-heat-edge-and-first-viscous-stress-stability-audited`.

This batch checks the actual HEAT-TAIL FORMULAS directly. It proves strict
growth and oriented cone bounds on the terminal outer interval, and obtains
a sharper fact than the previous weighted estimate: the first axial-viscous
stress correction is positive in its angular component and zero in its axial
component there. Adding it with a nonnegative cutoff improves the cone.
Its ratio to the leading angular stress stays bounded even at the flat edge.

This closes an outer-edge estimate, not the complete leading-profile audit.
Identifying these backward stresses with the globally assembled forward
primitives still uses the exact moment and support identities. The inner
analytic construction, inner edge and compact modulation/repair certificate
remain to be independently rechecked. No unconditional pulse estimate,
complete residual improvement or blow-up result is asserted.

Parent: `f31516d84c8d46137e41c0882f5006e5d938d30a`. The source under audit is
the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.2), (4.11), (5.3), (5.9), (A.47)--(A.56), with SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The first-order outer-stress identity below is derived here. Earlier source,
evidence and GOAL.md bytes remain frozen.

## 1. Exact tail, parameters and the scale normalization

Keep the selected finite N,h and the same tail. Write A=1/2+h, D=1/2-h,
L=1-2h eta^2, d=1-eta^2. Put y=log(X/Xtail), Xb=Xtail exp(3), r=sqrt(2X),
and 1/2<=y<3. In the formulas in sections 1--5 we set q=1, so r is the
profile radius. For general q the physical radius is q^(1/2)r.

The defining heat integral and the tail are

\[
 H(Z)=\frac1{\Gamma(1+h)}\int_0^\infty e^{-v}v^h(1+Zv)^{-h}\,dv,
 \quad Z=\frac{2d}{X},\quad K=C_{tail}X^{-A}H(Z),\quad
 E_0=K f(y),\quad U_0=V_0=0,\tag{1}
\]
\[
 f(y)=1-\rho\sigma((3-y)/2),\qquad \rho=h/1000,
 \quad \sigma(s)=\frac{e^{-1/s^2}}{e^{-1/s^2}+e^{-1/(1-s)^2}}\ (0<s<1).
 \tag{2}
\]

Extend sigma by zero and one. Thus f'>=0, f=1 beyond y=3, and f'=0
through y=1. The heat edit is already complete by y=1/2. Assume Xtail>=100,
as follows from the reference schedule with XR>=100. Ctail is positive
and finite, although its actual numerical value and Xtail need not be
materializable.

The pressure is the backward centrifugal-pressure primitive. The exact
global moment identities make the previously defined forward T0 coincide
with the following backward stresses. They are a stated identification
premise; a smooth tail alone cannot remove an unknown forward integration
constant. This distinction is retained throughout this audit.

At general scale the leading stress is q^(-A-1/2)T0. The order-one stress
is q^(-A-1/2+2h)T1. They are different formal orders. Viscosity is one.

## 2. Bounds obtained from the defining heat integral

Let tH=-ZH'/H. Positivity of the Gamma weight gives 0<H<=1 and
0<=tH<h. Also

\[
 H\ge1-h(1+h)Z,\qquad
 -H'\le h(1+h),\qquad
 |H^{(j)}|\le(h)_j(1+h)_j.\tag{3}
\]

The last bound is for every fixed j and follows by differentiating the
positive integral. With X>=100 and 0<h<.01, (3) gives tH<h/4. The fixed
flat step has sigma'<=64: on [1/4,3/4] use sigma(1-sigma)<=1/4 and
2/s^3+2/(1-s)^3<=256. On s<=1/4, use the upper bound
exp(16/9-1/s^2)(2/s^3+128/27), whose two terms decrease as 1/s>=4
increases; its endpoint value is less than 64. The other tail follows by
reflection. Hence
f'/f<=.032h/(1-.001h)<h/4. Direct differentiation now gives

\[
 a=2+2h-2t_H-2f'/f,\qquad 2+h<a\le2+2h,\qquad b_s=0.\tag{4}
\]

These are whole-interval bounds including eta=+/-1; no division by d is
used. Let Hmin=1-2h(1+h)/Xtail>0. Then, on the closed terminal rectangle,

\[
 F_0=E_0/r\ge F_{min}:=
 \frac{C_{tail}X_b^{-A}(1-\rho)H_{min}}{\sqrt{2X_b}}>0.\tag{5}
\]

For the pulse representative g=(r partial_r F0,partial_r U0)=(-aF0,0),
so N=(-1,0), Kperp=(0,-1), and

\[
 |g|\ge(2+h)F_{min},\quad
 \lambda_0^2=2F_0^2(a-2)\ge2hF_{min}^2>0,\quad
 c_0=-\sqrt{(a-2)/2},\quad\sqrt{h/2}<|c_0|\le\sqrt h.\tag{6}
\]

This checks the shear normalization as well as the growth sign; a lost
factor of two in g would change the growth criterion.

For the actual fixed reference schedule, the left heat-compensation-patch
amplitude is

\[
 e_* =\exp[-119.5T+10.3-\lambda_{patch}(240T-19.5)],\quad
 h=e^{-8T},\quad T=e^{64}+10.\tag{7}
\]

The expression follows by integrating the first ramp, axial stage,
second ramp and power stage: their log-amplitude changes are -.2,
-T/2, -.5-lambda_patch/2 and
-(.5+lambda_patch)(240T-20), in addition to log P=T+1.
Since T>=64, e*/h^14<=exp(-7.5*64+10.3)<1. After that patch the unedited
reference swirl decreases; all compact leading edits end before the tail.
Dividing out f costs at most (1-rho)^-1<2 and H<=1. Therefore the actual
tail satisfies K<=2h^14<=2h^4. This is a schedule bound, not a sampled
small amplitude or a new choice of h.

## 3. Positive leading stress and its orientation

Put Hrad(v)=r_v^3 K_v and B(v)=2+2h-2tH(v)>0. The exact backward formulas
at q=1 are

\[
 T_{0\theta}=\frac{2K}r f'
 +\frac1{r^2}\int_y^3 r_vK_v B(v) f'(v)\,dv
 +\frac1{2Lr^2}\int_y^3 H_{rad}(v)f'(v)\,dv,\tag{8}
\]
\[
 T_{0z}=\frac\eta{Lr}\int_y^3(r_v^2-r^2)K_v^2 f(v)f'(v)\,dv.\tag{9}
\]

Every term of (8) is nonnegative, and its last integral is positive for
every y<3. Since K decreases with radius, compare the integrands of (9)
with the last integral in (8) to obtain

\[
 T_{0\theta}>0,\quad |T_{0z}/T_{0\theta}|\le2K(y)\le4h^4,
 \quad 2-(a-2)(T_{0z}/T_{0\theta})^2\ge2-32h^9>1.99.\tag{10}
\]

The sign is part of the result. With n=T0/|T0|, -n dot N=n_theta>0.
Together with (6), the pulse cone ratio is at most 4h^(9/2)<1.
These bounds prove the oriented positive cone on the terminal interval,
subject only to the forward/backward identification stated in section 1.

## 4. The first axial-viscous stress has a positive representation

On this outer interval every positive-order velocity and pressure coefficient
vanishes. This is the common-support conclusion of the global coefficient
construction; it is not inferred from U0=0 alone. Thus the order-one
tangential source is r_theta,1=-partial_zz(Kphysical f), and r_z,1=0.
The physical heat velocity Kphysical(r,t) is independent of z. At q=1,

\[
 y_z=-2\eta/L,\quad y_z^2=\alpha:=4\eta^2/L^2,\quad
 y_{zz}=\beta:=-\frac{2(d-2D\eta^2)}{L^2}
                  -\frac{8h\eta^2d}{L^3}.\tag{11}
\]

For general q these two coefficients acquire factors q^(-D) and q^(-2D),
respectively. Consequently r_theta,1=-K(alpha f''+beta f') at q=1.
Both alpha and beta are constant DURING radial integration: q,eta are fixed
when r varies. Retaining the stress integration weight gives

\[
 T_{1\theta}=-\frac1{2r^2}\int_y^3 H_{rad}(v)
                         [\alpha f''(v)+\beta f'(v)]\,dv,\quad T_{1z}=0.
 \tag{12}
\]

The integrand in (12) can change sign; its sign does not decide the sign
of the integral. Integration by parts, using f'(3)=0 and
partial_y log Hrad=1-h+tH, gives

\[
 \boxed{T_{1\theta}=\frac{\alpha rK}{2}f'
       +\frac1{2r^2}\int_y^3 H_{rad}(v)C(v)f'(v)\,dv,}\tag{13}
\]
\[
 C=\alpha(1-h+t_H)-\beta
  =\frac2{L^2}+\frac{8h\eta^2d}{L^3}+\alpha t_H>0.\tag{14}
\]

The cancellation producing 2/L^2 uses D=1/2-h. In particular it retains
the beta f' term. At eta=0 alpha=0 and beta=-2, so dropping beta would
erase the ENTIRE nonzero correction. No division by eta occurs in (13).

Since L>=.98 and tH<h, (14) gives C<7. Compare the boundary term of (13)
with the first term of (8), and its integral with the last term of (8):

\[
 0<T_{1\theta}\le B_{tail} T_{0\theta},\qquad
 B_{tail}:=2X_b/(.98)^2+7.\tag{15}
\]

The large but finite Xb factor is necessary. Uniformity as Xb tends to
infinity is not asserted. Relative O(q^(2h)) smallness is valid uniformly
through this flat edge, with this fixed construction-dependent constant.

## 5. Remove the flat factor before taking an endpoint

Let delta=3-y and zeta_b=exp(-4/delta^2). The actual step gives
f'(3-u)=rho exp(-4/u^2)u^-3 Bflat(u), with Bflat smooth and positive,
Bflat(0)=8 exp(1). For any smooth coefficient J,

\[
 \int_0^\delta e^{-4/u^2}u^{-3}J(u)\,du
  =\frac{e^{-4/\delta^2}}8\int_0^\infty e^{-w}
            J\left(\frac\delta{\sqrt{1+\delta^2w/4}}\right)dw.\tag{16}
\]

Dominated differentiation is valid at every fixed order because the
remaining derivatives are bounded by exp(-w) times a polynomial in w.
Equations (8)--(9),(13) therefore give smooth coefficients on a fixed closed
collar, jointly in delta,eta:

\[
 T_{0\theta}=\zeta_b\delta^{-3}b_0,\qquad
 T_{0z}=\zeta_b\delta^3 b_z,\qquad
 T_{1\theta}=\zeta_b\delta^{-3}b_1,\tag{17}
\]
\[
 b_0(0,\eta)=16\rho e K_b/r_b>0,\quad
 b_z(0,\eta)=\eta\rho e r_bK_b^2/(8L),\quad
 b_1(0,\eta)=4\alpha\rho e r_bK_b.\tag{18}
\]

The axial factor uses r_v^2-r^2, which contributes another delta^3 after
(16). In particular the following ratio extends SMOOTHLY, including eta=0:

\[
 \left.\frac{T_{1\theta}}{T_{0\theta}}\right|_{\delta=0}
       =\frac{2\eta^2X_b}{L^2}.\tag{19}
\]

The value at eta=0 is zero, despite positive T1theta at every interior
radius. Neither taking the raw numerical quotient 0/0 nor dropping the
radial weight recovers this limit. Every fixed derivative of T0,T1 is
bounded by C_I zeta_b delta^(-m_I). The positive b0 gives a common leading
lower weight. This rechecks the outer weighted edge estimate directly.

## 6. Stability of the actual background stress on this interval

Choose the already free Borel cutoff with 0<=chi<=1; for example
chi(s)=1-sigma(2s-1). This is an allowed choice in the earlier cutoff
construction and changes none of its exponent budgets. Higher-order
stresses Tn, n>=2, vanish on the outer interval by their common support.
Thus, conditional on those global support and exact moment conclusions,
the summed normalized background stress there is exactly

\[
 \widehat T_B=T_0+c(q)T_1,\qquad
 c(q)=q^{2h}\chi(a_1q)\ge0.\tag{20}
\]

Its angular component is positive and its axial component equals T0z.
The ratio magnitude in (10) can only decrease. The velocity on this interval
is exactly the leading heat-tail velocity, so a and the growth bound (6)
are unchanged. Equations (17)--(18) give the same smooth limiting direction
(1,0). This proves outer cone preservation WITHOUT using an unproved
bound on delta^-m near the edge and without a smallness assumption on c.
The sign restriction on chi is essential to that last statement; a signed
cutoff would instead require a separate small-q estimate using (15).

This sharper fact is local to the prescribed outer heat interval. A generic
weighted perturbation zeta_b delta^-m need not preserve a cone: for example
DeltaTtheta=-q^(2h) zeta_b delta^-4 overwhelms
T0theta=zeta_b delta^-3 as delta tends to zero, for every fixed q>0.
The generic warning in the frozen background note remains correct.

## 7. Verification and remaining foundation work

The diagnostics evaluate the EXACT Gamma-integral heat function through
its confluent-hypergeometric representation, and crosscheck it by independent
positive integration. They compare the signed formula (12) with the positive
formula (13), retain physical axial differentiation at nonunit q, and check
the cylindrical stress-divergence sign. Edge diagnostics integrate (16)
after removing the flat factor and keep eta=0,+/-1 and tiny deltas. A tighter
quadrature is used for refinement. These are diagnostics; sections 1--6
give the continuum argument. The actual huge schedule is not numerically
evaluated, and no convergence of a raw infinite heat Taylor series is used.

The remaining leading-profile gates are now the INNER analytic collar and
its stress-direction factor, the compact modulation and final repair with
their exact pressure memory, and the global moment identification of the
tail primitives. In the pulse replay's notation one must still establish,
on the remaining profile domain at the same N,h,

\[
 \inf F_0>0,\quad\inf |g_0|>0,\quad\inf\lambda_0^2>0,\quad
 \inf(-\widehat T_0\cdot N)>0,\quad
 \sup\left|c_0\frac{\widehat T_0\cdot K}{\widehat T_0\cdot N}\right|<1,
 \tag{21}
\]

with the inner weighted edge bounds. An exterior calculation cannot certify
those compact or inner nonlinear constructions. The reserved physical mean
supports, full localization assembly and unconditional residual gate remain.

```sh
python3 research/navier_stokes_cascade/results/outer_edge_stability/test_audit.py -v
python3 research/navier_stokes_cascade/results/outer_edge_stability/audit.py \
  --verify-record research/navier_stokes_cascade/results/outer_edge_stability/evidence.json \
  --output outer_edge_stability_recheck.json
```
