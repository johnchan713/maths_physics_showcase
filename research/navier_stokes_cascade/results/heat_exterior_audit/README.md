# Compensated heat exterior for the reference flow

This checkpoint constructs the exact three-moment heat edit of the fixed
reference family. It gives explicit bounds preserving the stress cone from
the compensation patch through terminal coordinate 1/2, and a positive
backward-stress argument on the remaining open terminal interval. The stress
vanishes at the outer edge with a smooth limiting direction.

The construction is a **reference field for X>0**. A regular axis and its
annular attachment, the full PDE corrections and an admissible smooth force
remain unverified. It is not a blow-up solution. Numerical surrogates below
have explicit truncation budgets; small numerical moment residuals do not
establish exact closure.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.9), (4.11), (4.15)-(4.16), Lemma 4.5 and Appendix A.6-A.7, pages 138-144.
Frozen PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent: `63d325fb7d218c53bcb052c38c1b5cd06ea964eb`.
The quantitative bounds and relative-bump matrix below are an explicit
specialization to our fixed family, not a verification of the entire paper.

## 1. Define the edit and the three equations

Retain Md=64, T=exp(64)+10, lambda=exp(-4T), h=lambda^2, Tf=1000,
co=.001 and XR>=100. Let Xtail be the start of the terminal interval and
XK=exp(.2) Xtail. On that interval and its continuation, put

\[
\mathcal H(Z)=\Gamma(1+h)^{-1}\int_0^\infty
 e^{-v}v^h(1+Zv)^{-h}\,dv,\quad Z=2(1-\eta^2)/X,
\]
\[
\chi(y)=\sigma((y-.2)/.3),\qquad
E_{\rm new}=E_{\rm ref}\{1+\chi(y)[\mathcal H(Z)-1]\}.
\]

Here y=log(X/Xtail), and chi is zero before .2 and one after .5.
The exact positive integral defines the heat kernel; a finite Taylor
polynomial is used only for the numerical diagnostics.

Use the second reserved patch, whose original coordinate is
(Tw-20,Tw-15). Its start is

\[
X_*=X_R e^{241T-18},\qquad
e_*=\exp[-119.5T+10.3-\lambda(T_w-19.5)].
\]

With t=log(X/Xstar), f=(1+eta^2)^(-1), the original swirl is
`E=e_star f exp[-(.5+lambda)t]`. Add the relative correction

\[
E\mapsto E(1+\varepsilon),\qquad
\varepsilon=\sum_{j=0}^2 c_j(\eta)
 \sigma'((t-t_j)/.3+.5),\quad t_j=.5+2j.
\]

All supports lie strictly inside this five-unit patch. Both edited regions
have U=0, so M and J are unchanged at every radius. Restore the other three
total moments: pressure increment Cp, energy moment S, and renormalized I.

Normalize their compensation rows by `e_star^2 f^2`,
`Xstar e_star^2 f^2` and `sqrt(2) Xstar^(3/2) e_star f`, respectively.
Use -S for the second row. The exact equations are

\[
B c+q(c,c)=z,\qquad
z=\left(-\frac{\Delta C_p}{e_*^2 f^2},
\frac{\Delta S}{X_*e_*^2 f^2},
-\frac{\Delta I}{\sqrt2 X_*^{3/2}e_*f}\right).
\]

The matrix row rates are `-1-2lambda`, `-2lambda`, `1-lambda`:
`B[k,j]=integral exp(rate[k]*t) beta_j(t) dt`. The first two quadratic
rows are `integral exp(rate*t) (sum cj beta_j)^2/2 dt`; the last is zero.
Disjoint supports remove cross-products between different spatial bumps.
The normalized matrix and quadratic map are independent of eta.

## 2. Prove existence of the exact small correction

Use `||g||C1=sup|g|+sup|g_eta|`, a Banach algebra norm. The bump matrix is
a diagonally scaled Vandermonde matrix with nodes
`exp(-2-4lambda)`, `exp(-4lambda)`, `exp(2-2lambda)`. These remain separated
on 0<=lambda<=.01. Positive bump mass .3 and its width give an outward
inverse bound below 15; use the conservative bound **50**. Since
`sigma'^2<=64 sigma'`, the quadratic map has norm below 100.

Write alpha=chi(1-H). Positivity and differentiation of the Gamma integral
give, uniformly over eta,

\[
0\le\alpha\le3h/X,\quad |\alpha_\eta|\le5h/X,\quad
\|\alpha\|_{C^1}\le8h/X.
\]

Thus `||(1-alpha)^2-1||C1<=18h/X`. Put eK=Eref(XK).
The future reference satisfies `Eref(XK*x)<=2 eK x^(-.5-h)`.
Positive integration over x>=1 yields

\[
\|\Delta C_p\|_{C^1}\le18h e_K^2/X_K,\quad
\|\Delta S\|_{C^1}\le36h e_K^2,\quad
\|\Delta I\|_{C^1}\le16\sqrt2 e_K\sqrt{X_K}.
\]

The angular integral contains `h integral_1^infinity x^(-1-h) dx=1`.
Setting h=0 before integration would erase a finite contribution.
Along the reference, H=sqrt(2X)E decreases after Xstar, and its steep
hold multiplies H by h^4. Consequently
`rH=eK sqrt(XK)/(e_star sqrt(Xstar))<=h^4`, and
`rE=eK/e_star=rH/sqrt(XK/Xstar)`.

The C1 norms of 1/f and 1/f^2 are at most four and twelve. Substitution
in the three normalized rows therefore gives

\[
\|z\|_{C^1}\le\delta:=64h^4/X_*,\qquad
8\cdot50^2\cdot100\delta<1.
\]

The quadratic contraction supplies a unique small root, with

\[
\boxed{\max_j\|c_j\|_{C^1}\le C:=6400h^4/X_* .}
\]

Implicit differentiation gives a smooth root at every fixed derivative
order. No division by 1-eta^2 is used. In particular, the correction values
are zero at eta=+-1, while their eta derivatives generally are not zero.
The exact root restores all three totals. The axis pressure datum and all
profiles and cumulative data before the compensation patch are unchanged.
The original renormalized I total is zero: after the original terminal
interval its eta-independent Q=0 gives I=XHpow/(1-h).

## 3. Control the changed stress through terminal 1/2

All bounds here concern the exact edited profile. They include the interior
of the compensation patch, where the cumulative changes have not yet
cancelled. The inherited reference bounds on this interval are

\[
Q\ge h/2000,\quad 2+3h/2\le a\le4,\quad |b_s|<15,\quad |w|<1800,
\]
\[
a-b_sw>1.48,\qquad
2-2b_sw-b_s^2/a-(a-2)w^2>.82.
\]

For the late power-law part, the preceding source bound gives
`abs(w)<=60000 T e_star/lambda<1`; e_star<=h^14. The pulse audit supplies
the 1800 bound and its two strict margins. The post-pulse audit supplies
the remaining interval. Thus these are whole-interval, all-angle bounds.

Because both edits have U=0, bs is unchanged everywhere. Let delta denote
an edited-minus-reference quantity. The exact numerator changes are

\[
\Delta V=(1-h)\Delta I-D\eta\Delta I_\eta,\qquad
\Delta N=(4h\eta\Delta S-d\Delta S_\eta)/X
             +4A\eta\Delta C_p-d\Delta C_{p,\eta}.
\]

Here Q+W=V/(XH), w=N/(EQ), d=1-eta^2, D=.5-h and A=.5+h.
The following elementary bounds explain the constants used by the audit.

| Region | Bounds before taking stress ratios |
|---|---|
| Inside compensation | `abs(epsilon)<=64C`, `abs(epsilon_t)<=20000C/.3`; `||Delta I||C1<=400 sqrt(2) Xstar^(3/2) e_star C`; `||Delta S||C1<=10 Xstar e_star^2 C`; `||Delta Cp||C1<=10 e_star^2 C`. |
| After compensation, before heat | The cumulative changes equal the negatives of the complete heat discrepancies. `XH` is nondecreasing and `E>=eK`, so `abs(Delta Q)<=32h^4/Xstar` and `abs(Delta N/(EQ))<=2e5/Xstar`. |
| Within the heat edit | The cumulative changes are negatives of the remaining heat integrals. `abs(Delta Q)<=50/X`, and `abs(Delta N/(EQ))<=2e5/X`. |

For the first row, `E>=e_star/40`, `Q+W<=12`, and the exact reciprocal
formula gives `abs(Delta Q)<=5000C`. It also gives
`abs(Delta N/(EQ))<=4e6 C/h` and `abs(Delta w)<=1e12 C/h`.
For the last row, `abs(epsilon)<=3h/X`, Q+W<2 and the remaining-integral
bounds in section 2 apply from the current radius. The relative change in
Q is at most `1e5/(hX)`. These facts, and the unchanged-field middle row,
are all bounded by the following common budgets:

\[
\frac{|\Delta Q|}{Q}\le\epsilon_Q=\frac{10^{12}}{hX_*},\qquad
|\Delta w|\le\epsilon_w=\frac{10^{17}}{hX_*},\qquad
\frac{|\Delta a|}{h}\le\epsilon_Q.
\]

The compensation slope bound follows by differentiating log(1+epsilon).
For the heat edit, `abs(Delta l)<=700h/X` follows from
`chi'<=64/.3`, while its positive part is at most `4h/X`.
All reciprocal denominators above exceed .99.
Since `log(hXstar)=log(XR)+233T-18`, outward evaluation at T=64 gives
`0<epsilon_w<4.150e-6454`; both budgets are below 1e-6000.
The actual Md64 logarithms are stored separately.

It follows that the edited Q>=h/4000 and a-2>=h. Its first ratio margin
exceeds 1.47; its second exceeds .81. Indeed their losses are bounded by
`h epsilon_Q+15 epsilon_w` and `1e7(h epsilon_Q+epsilon_w)`.

For the finite-radius test, put p=ps1=XQ/L. Then
`p>=hXstar/4000>1e12`, `vs<120`, `c=Pc/p>.29`, `abs(c)<14000` and
`G=2c^2-(vs-2)(Jc/p)^2>.81`. Therefore

\[
P_c>v_s,\qquad
\frac{2(P_c-v_s)^2-(v_s-2)J_c^2}{p^2}
\ge .81-4\cdot14000\cdot120/10^{12}>.80.
\]

Lemma 4.5 establishes the strict admissible cone through terminal y=1/2.
The sufficient reference radius remains XR>=100.

## 4. Recover the remaining exterior stress

The restored moments give the tail identities
`I=XHpow/(1-h)-integral_X^infinity(H-Hpow) dx`,
`S=integral_X^infinity E^2 dx/2`, and the backward pressure integral.
These apply directly to our reference field on X>0. At large X,
`H/Hpow=1+O_h(1/X)`, hence `Q=O_h(1/X)` and `N=O_h(E^2)`.
The weighted stresses obey `r^2 Ttheta=O_h(X^-h)` and
`r Tz=O_h(X^-2h)`, so both integration constants at infinity vanish.
Integrating the exact radial stress identities from infinity therefore
gives the backward formulas below. This argument establishes reference
stress identities; Cartesian smoothness at the axis is still an open task.

Use physical scale q=1 to display the profile formulas. Set
`r=sqrt(2X)`, `K=Epow(X) H(2d/X)`, `f=fo(y)` and L=1-2h eta^2.
All q factors cancel from the directional ratio. For .5<=y<3,

\[
T_\theta=\frac{2K}{r}f'
 +\frac1{r^2}\int_y^3 r_vK_v B_v f'_v\,dv
 +\frac1{2Lr^2}\int_y^3 r_v^3K_v f'_v\,dv,
\]
\[
T_z=\frac{\eta}{Lr}\int_y^3(r_v^2-r^2)K_v^2f_vf'_v\,dv,
\quad B_v=2+2h+2Z_v\mathcal H'(Z_v)/\mathcal H(Z_v)>0.
\]

The angular stress is positive. K decreases with radius, f<=1, and the
last positive angular integral alone yields

\[
|T_z/T_\theta|\le2K(y)\le4h^4.
\]

Also `-ZH'/H<h/4` and `f'/f<h/4`, giving `a>2+h` and `a<=2+2h`.
Thus `Pc-a=Ttheta/F>0` and the directional cone margin is at least
`2-32h^9>1.99` on this entire open terminal interval. Beyond y=3 the
profile is the heat solution and its reference tangential stress is zero.

At eta=1 the heat edit's value is zero, but its tail angular derivative is
`I_eta/Hpow=-4(1+h)`. Consequently `ps1=2+2h=a`, as required for zero
heat stress. Dropping that derivative gives ps1=0 and a false nonzero
stress. The audit keeps this failure explicitly.

## 5. Factor the vanishing edge without division by zero

Put delta=3-y. For 0<u<2, write

\[
f'(3-u)=\rho e^{-4/u^2}u^{-3}B(u),\quad
B(u)=\frac{A(u)[8+u^3/(1-u/2)^3]}{[A(u)+e^{-4/u^2}]^2},
\quad A(u)=e^{-(1-u/2)^{-2}}.
\]

Here B is smooth and positive, B(0)=8e. Substituting
`u=delta/sqrt(1+delta^2 w/4)` in each backward integral removes the
common exponential and leaves the integrable positive weight exp(-w).
For the axial integral the potentially cancelling radius difference uses

\[
\frac{\delta-u}{\delta^3}
=\frac{w/4}{\sqrt{1+\delta^2w/4}[1+\sqrt{1+\delta^2w/4}]}.
\]

Dominated differentiation gives smooth coefficients btheta and bz with

\[
T_\theta=e^{-4/\delta^2}\delta^{-3}b_\theta,
\quad T_z=e^{-4/\delta^2}\delta^3 b_z,
\]
\[
b_\theta(0)=16\rho K_b e/r_b>0,\qquad
b_z(0)=\eta\rho r_bK_b^2e/(8L).
\]

The ratio extends as `delta^6 bz/btheta`, tending to zero. Its coefficient
limit is `eta Xb Kb/(64L)`. The size of this coefficient can depend strongly
on the fixed schedule; a uniform collar width is not inferred from samples.
At the edge the stress is zero and its extended direction is (1,0), with
directional margin two. No stress ratio 0/0 is evaluated.

## Numerical construction and error budgets

`construction.py` evaluates a finite Taylor expansion of the exact heat
kernel. If bn is the coefficient of Z^n in (H-1)/h, then
`b1=-(1+h)` and `b[n+1]=-b[n](h+n)(1+h+n)/(n+1)`.
The derivative integral bounds the degree-N remainder by
`h abs(b[N+1]) Z^(N+1)`. Its eta C1 bound has the additional factor 2N+3.
The audit propagates this positive error through all three moment integrals.

Each finite radial piece is integrated with positive quadrature; every
infinite power piece is integrated analytically. The angular n=1 term
keeps `h integral exp(-hy) dy` as a finite quantity. Quadrature refinement
is a diagnostic, not a certified quadrature error bar.

The actual Md4 member of the chosen family tests the enormous radial scales.
It retains its historical axial failure. A separate moderate-scale fixture
tests nonlinear root evaluation and independent source integration; it is
not a substitute for the Md64 reference.
The schedule and moment quadratures use eight panels; the earlier post-pulse
audit already demonstrated that four panels can miss a 1e-30 terminal
refinement tolerance. Edge samples are selected by their declared grid
positions, avoiding equality tests between decimals rounded at different
arithmetic precisions.

At the actual scales even hundreds of digits cannot resolve the quadratic
root correction against its linear part. `Compensation.graded` stores each
homogeneous correction separately. With a=max(abs(linear coefficients)) and
r=4*50*100*a<1/2, the tail after N grades is bounded by `a*r^N/(1-r)`.
The error in the exact-kernel target contributes at most 100 times its C1
error bound. Ordinary Newton iteration can return the linear coefficients
unchanged; the nonzero quadratic grade exposes that arithmetic failure.
These are truncation budgets. They do not include certified quadrature or
floating-point rounding error; the recorded refinements test those numerical
paths without replacing the analytic existence argument. The finite heat
polynomial also has a nonzero ODE residual, retained explicitly in the audit.

`exterior.py` evaluates both the direct backward stress and the factored
edge integrals. Tests compare these independent integration paths and keep
the positive flat factor after binary64 rounds it to zero. The exact
reference construction rests on sections 1-5, not these finite samples.

## Reproduction and next obligation

```bash
python research/navier_stokes_cascade/results/heat_exterior_audit/test_audit.py -v
python research/navier_stokes_cascade/results/heat_exterior_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/heat_exterior_audit/evidence.json \
  --output heat_exterior_reproduced.json
```

`bounds.py` checks the outward constants. `construction.py` defines the
heat approximation, tail moments, exact polynomial bump map and graded root.
`exterior.py` handles positive stress and its edge factors. `audit.py` pins
the protocol and inherited evidence; `test_audit.py` checks identities,
remainder budgets and deliberate failures.

The next construction is the regular axis and its annular attachment with
all five matching moments. The previously attempted direct heat attachment
failed and is not reused. Full momentum-residual cancellation and smooth
forcing remain separate requirements before any blow-up claim.
