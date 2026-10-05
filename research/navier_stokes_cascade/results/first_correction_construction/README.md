# First inner correction and the next derivative gap

This checkpoint gives a **conditional analytic construction of the first
background correction on the actual protected core**. It extends the inherited
angular transfer to order three before modulation, derives a regular six-variable
system, supplies an actual core operator envelope, and proves convergence of its
Picard series. The five positive-order patch moments have a linear inverse.
Finite manufactured profiles independently check the coefficient equations.

The conclusion depends on the frozen analytic leading-core construction. It is
not an independent review of that construction. The correction is constructed
on `0<=X<=3/Lambda`, strictly inside `Xa=4/Lambda`; its global annular extension
and interaction with physical waves are not verified here. The finalized
post-modulation third-derivative bound, all-order summation, admissible forcing,
energy, and blow-up of the complete field remain open.

Parent: `5a94f20dfc1de702a507cc9764bfc4a68ff01ab8`. Source manuscript:
[Navier–Stokes](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
especially (5.1)–(5.6), Lemma 5.1 and Appendix B. SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The source suggests the sparse Picard argument; the coefficients and estimates
below are derived from the already audited physical residual equations.
Earlier evidence files and GOAL.md remain unchanged.

## What is and is not established

| Result | Scope |
|---|---|
| Incoming five-moment transfer in factorial angular C3 | Actual unmodulated profile, conditional on inherited core/reference bounds |
| First coupled background coefficient | Convergent analytic construction on `X<=3/Lambda`, with zero axis traces |
| Operator envelope `C_op=C^16` | Actual protected core on an explicit complex angular neighborhood; no giant parameters are materialized |
| First-order five-moment inverse | Exact linear map on the reserved power-law patch, conditional on its leading support hypotheses |
| Residual refinement | Manufactured polynomial data; radial degrees 3, 5, 7 |
| Protected-core swirl lower bound | Leading background plus this one correction only |
| Post-modulation C3, global coefficient support, waves, infinite sum and full blow-up | Open |

## 1. The actual C3 transfer before modulation

Write

\[
\|v\|_3=\sum_{j=0}^3\frac{\sup_{[-1,1]}|\partial_\eta^jv|}{j!}.
\]

This is a Banach algebra. A norm bound R controls the raw third derivative by
`6R`. With the same S and rho as the exact axis-series certificate, its coefficient
weights and the generating function give, for `0<=Y<=5`,

\[
\|\Phi(Y,\cdot)\|_3,\|u(Y,\cdot)\|_3
\le S\rho^{-3}\sum_{j=0}^3\frac{(4/3)^{j+1}}{(j+1)^2}
=\frac{544}{243}S\rho^{-3}<D_3:=3S\rho^{-3}.
\]

No radial or angular tail was truncated in this inequality. Cauchy's bound
`|zeta_*''|<=8M/r0^2` on the inherited rectangle yields

\[
\|\log\phi_*\|_3\le
L_{*,3}:=\Lambda M(2+r_0^{-1}+4/(3r_0^2)).
\]

The exact lower bound `Phi>1/4` gives
`||log Phi||_3<=100 D3^3`. The earlier reference certificate already includes
angular C3 bounds for `p1,r` and `ns,r`; radial activation cutoffs are independent
of eta. Integrating those bounds just as in the C2 transfer gives

\[
\|\log(CE)\|_3\le
L_{*,3}+100D_3^3+55\,10^6\Lambda K^8+2\log(27.5\Lambda)+10<B.
\]

The axial change is at most `1e8 K^8 delta`, hence the terminal offset obeys
`||g||_3<=3e-18`. On disks of radius 1/4 about real eta in [-1,1],
`|log f|<2` for `f=(1+eta^2)^(-1)`; Cauchy's estimate gives
`||log f||_3<=2(1+4+16+64)=170`. The shape interpolation therefore still has
`||log(CE)||_3<=204B` and `||E||_3<=C^(-1) exp(204B)`.

The five moment integrals can now be differentiated three times. In the axis
pressure integral retain the factor `E^2/X`: the core estimate is
`||E||_3<=C^(-1) sqrt(X) exp(B)`, so its integrand has a finite bound at zero.
The continuation integral keeps the logarithmic factor `log(Xsep/Xa)`.
Normalizations cost

\[
\|f^{-1}\|_3=5,\qquad \|f^{-2}\|_3=4+8+8+4=24.
\]

The energy-row coefficient increases from 6500 to
`24(101+160+64)=7800`, still below the previously generous `1e6` cap. Each
normalized incoming row remains below `epsilon/4`, and the original core row
bound is below `.255 epsilon`. Restoring the terminal offset adds at most
`gb+gb^2`, where `gb=3e-18`. Consequently

\[
\|z\|_3<2.86\,10^{-17},\qquad
\|c\|_3<5.72\,10^{-14},\qquad
\sup|c_{\eta\eta\eta}|<3.432\,10^{-13}.
\]

The same small-root contraction works in this factorial norm; uniqueness
identifies this root with the previously constructed joining correction.
This changes neither the field nor N. Outward arithmetic at 80 and 110 digits
checks the enlarged constants.

For the raw derivative equation `Bc+Q(eta)[c,c]=z`, differentiation gives

\[
\begin{split}
(B+2Q[c,\cdot])c'''={}&z'''-6Q[c',c'']-6Q'[c,c'']\\
&-6Q'[c',c']-6Q''[c,c']-Q'''[c,c].
\end{split}
\]

The five terms are independently differentiated in the diagnostic. They must
not be omitted because a lower-order target is small.

**This is not the final modulated C3 bound.** The stress coordinates already
contain one angular derivative of the moments. Taking three derivatives of
the loop requires those coordinates through order three, hence the original
moment inputs through order four. The new incoming C3 transfer cannot supply
that C4 input. The corresponding moving-phase, repair and radial bounds remain
open, with `N=1+floor(H^32)` fixed.

## 2. A regular first-correction system

Use the physical-residual notation
`A=1/2+h`, `D=1/2-h`, `lambda=2h`, `d=1-eta^2`, `L=1-2h eta^2`,
`b=-A-1/2`, `c=-A` and viscosity `nu>0`. Write `xi=sqrt(X)` and

\[
K_1=\frac1X\int_0^X U_1(x,\eta)\,dx-U_1,
\quad P=\partial_\xi F_1,\quad Q=\partial_\xi U_1,
\quad W=(F_1,U_1,K_1,\Pi_1,P,Q).
\]

The first coefficient's exact radial flux is

\[
V_1=\frac{X}{L}
\left[2(A-\lambda)\eta U_1-2(D+\lambda)\eta K_1
-d\,\partial_\eta(U_1+K_1)\right].
\]

In particular `D+lambda` cannot be replaced by D. The coupled physical equations
become

\[
\boxed{\partial_\xi W+\xi^{-1}\operatorname{diag}(0,0,2,0,3,1)W
=A_0(\xi,\eta)W+A_1(\xi,\eta)\partial_\eta W+f_1(\xi,\eta).}
\]

`system.py` specifies every entry. The first four equations are

\[
F_{1,\xi}=P,\quad U_{1,\xi}=Q,\quad
K_{1,\xi}+2K_1/\xi=-Q,\quad
\Pi_{1,\xi}=4\xi F_0F_1-\xi\omega_0,
\]

where `v=V0/X` and `omega0=Omega0/X` are regular. Their explicit source is

\[
\omega_0=
\frac{v+Xv_X+D\eta v_\eta}{L}
+v(v/2+Xv_X)
+\frac{U_0[-2\eta(v+Xv_X)+dv_\eta]}{L}
-2\nu(2v_X+Xv_{XX}).
\]

Put `Hc=D eta+d U0`, `JF=F0+X F0_X`, `JU=X U0_X`,
`CF=Z_(b-D) Z_b F0`, `CU=Z_(c-D) Z_c U0`. The derivative matrix is

\[
A_1=\frac{2}{\nu L}
\begin{pmatrix}
0&0&0&0&0&0\\
0&0&0&0&0&0\\
0&0&0&0&0&0\\
0&0&0&0&0&0\\
H_c&-dJ_F&-dJ_F&0&0&0\\
0&H_c-dJ_U&-dJ_U&d&0&0
\end{pmatrix},
\quad
f_1=(0,0,0,-\xi\omega_0,-2C_F,
2\eta X\omega_0/(\nu L)-2C_U).
\]

Pressure coupling in the Q equation uses
`Pi1_X=2F0F1-omega0/2`; it is retained in A0 and f1. For an arbitrary
manufactured correction, let r1,t1,z1 denote the independently calculated
normalized radial, tangential and axial first residual coefficients. The
system defect is exactly

\[
(0,0,0,r_1/\xi,-2t_1/\nu,-2(z_1+\eta r_1/L)/\nu).
\]

This identity is checked at 54 combinations of radius, angle, h and viscosity.
Controls remove pressure coupling, lambda or the angular-derivative matrix;
all produce nonzero errors.

## 3. Actual core bound and convergent Picard construction

For the actual construction set `nu=1` and

\[
a^2=3/\Lambda=3X_a/4.
\]

This interval lies inside the exact analytic core. It precedes every activation,
modulation and repair patch. No tiny collar endpoints need to be added or
rounded. The core source series is holomorphic in eta, with all angular orders
controlled by its weighted norm. For complex `|Y|<=4` and angular distance
less than `rho/4` from [-1,1], the generating function gives

\[
|\Phi|,|u|\le
\frac{S}{1-4/20-1/4}<2S.
\]

The inherited complex axis factor has `|phi_*/C|<=1`. Thus F0 and U0 have
outer-neighborhood bounds `2S` and `6+2S`, respectively. Cauchy estimates on
`dist(eta,[-1,1])<rho/8` and a Y disk of radius one half bound the F0,U0 jets,
their radial averages and the regular v0 jets required above. A deliberately
generous common bound is

\[
K_{core}=10^{10}(1+S)(1+\Lambda)^2(1+\rho^{-1})^3\le C^4.
\]

The largest angular requirement is three, from differentiating the radial
average in v0 twice; the largest radial requirement is two. Multiplication by
eta,d and division by L cost fixed constants: here `|eta|<1.01` and
`|2h eta^2|<.03`, so `|L^(-1)|<2`. In the displayed matrix and regular source
formulas, each coefficient is a sum of fixed-constant multiples of at most two
such jets. Direct triangle bounds give, in maximum row-sum norm,

\[
\max(\|A_0\|,\|A_1\|,\|f_1\|)\le10^6K_{core}^2\le C_{op}:=C^{16}.
\]

`inner_operator_certificate` checks both parameter comparisons outwardly.
This is an explicit bound for the actual core, conditional on the inherited
analytic series; it is not an evaluation of its enormous parameter values.

Let `di=(0,0,2,0,3,1)` and define the regular Volterra inverse

\[
(Gg)_i(\xi)=\int_0^\xi(s/\xi)^{d_i}g_i(s)\,ds,
\qquad \mathcal K=G(A_0+A_1\partial_\eta).
\]

All exponents are nonnegative. The exact matrix structure maps the first four
coordinates into the last two and annihilates the last two. Consequently
`A1(xi) Dg A1(s)=0` for every diagonal matrix Dg. This also holds with any eta
derivative of the right A1. Thus two adjacent derivative factors in a Picard
word vanish, including the term where the derivative hits A1 itself. A length-k
word has at most `p=ceil(k/2)` surviving derivatives.

Fix a smaller angular neighborhood with strip loss `0<Delta<=1`. Cauchy loss
can be allocated as `Delta/p` to each derivative. Ordered radial integrals
have volume `a^(k+1)/(k+1)!`; the diagonal weights are at most one. Counting
words by `2^k`, a convenient upper bound is

\[
\|\mathcal K^kGf_1\|
\le\frac{(2C_{op}a)^{k+1}}{(k+1)!}
\max(1,p/\Delta)^p.
\]

Its kth root tends to zero: the surviving derivative loss grows like
`k^(k/2)`, whereas the radial factorial grows like `k^k`. More explicitly,
with `z=2C_op a` and `beta=e z^2/Delta`, even and odd lengths give

\[
\sum_{k\ge0}\|\mathcal K^kGf_1\|
\le (z+\beta)e^\beta.
\]

Indeed `j^j<=e^j j!` and `(2j+1)!>=(j!)^2`; the odd-length sum is bounded
by `exp(beta)-1<=beta exp(beta)`. This estimate does **not** require
`a C_op<1`. It proves uniform convergence on smaller angular neighborhoods of

\[
\boxed{W_1=\sum_{k=0}^{\infty}\mathcal K^kGf_1,\qquad W_1(0,\eta)=0.}
\]

The same estimate proves uniqueness by iterating the homogeneous integral
equation in the analytic class. The argument also holds for complex xi on
straight integration segments: the diagonal weights are t^di, `0<=t<=1`,
and `|Y|<=3`. Uniform convergence then gives radial analyticity. The spare
half-unit Y neighborhood permits a slightly larger radial disk, so the closed
interval endpoint is regular too. Angular derivatives converge on still
smaller neighborhoods. Their parity
is `(even,even,even,even,odd,odd)` in xi, so F1,U1,K1,Pi1 are smooth functions
of X, while P,Q are their regular xi derivatives. The radial coefficients also
follow the nonsingular recurrence in the next section. These regular traces
remove every apparent `1/xi` singularity.

This solves the **first inner background coefficient**. It does not include
the wave covariance or prove a global corrected solution.

## 4. Recurrence and manufactured computation

At the axis the exact zero-trace slopes are

\[
(F_1)_X(0)=-C_F(0)/4,\qquad
(U_1)_X(0)=-C_U(0)/2,\qquad
(\Pi_1)_X(0)=-\omega_0(0)/2.
\]

For ordinary radial coefficients `F1=sum f_n X^n`, the tangential diffusion
denominator is `2nu(n+1)(n+2)`; for U1 it is `2nu(n+1)^2`. Once lower rows
are known, the pressure recurrence is

\[
(n+1)\pi_{n+1}=2\sum_{i+j=n}(F_0)_i(F_1)_j
-\tfrac12(\Omega_0)_{n+1}.
\]

`series.py` implements these equations using angular Taylor rows, explicitly
retaining the derivative budget. The manufactured leading profiles are the
unchanged fixtures from the physical residual audit. Their higher radial rows
are genuinely zero because they are polynomials. They are not the actual core.
At three held-out points, degrees 3,5,7 reduce the first momentum coefficient
defects to roundoff; the largest final defect is below `2e-11`, and the ratio
of initial to final worst error exceeds `1e6`. The independent axis slopes
also agree. These diagnostics check recurrence algebra, not continuum existence.

## 5. The five positive-order patch moments

Put `E1=sqrt(2X) F1`. The total first coefficient moments are

\[
\begin{split}
m_1&=\int U_1\,dX,&m_2&=\int\sqrt{2X}E_1\,dX,\\
m_3&=\int(\Pi_1)_X\,dX,&m_4&=\int\sqrt{2X}(U_0E_1+U_1E_0)\,dX,\\
m_5&=\int[2U_0U_1-X(\Pi_1)_X]\,dX.
\end{split}
\]

Use the **positive-order patch I3**, not I1 already occupied by the leading
repair. Under its prescribed leading hypotheses, `U0=0`, accumulated leading
M0 is zero, and `E0=e(eta) exp[(-1/2-lambda_patch)y]` with `e>0` and
`lambda_patch>0`. Then V0 and Omega0 vanish there. Adding two U1 bumps and
three E1 bumps changes the five moments by the **linear** formulas

\[
\delta m_1=\int\delta U_1,\quad
\delta m_2=\int\sqrt{2X}\delta E_1,\quad
\delta m_3=\int E_0\delta E_1/X,
\]
\[
\delta m_4=\int\sqrt{2X}E_0\delta U_1,\qquad
\delta m_5=-\int E_0\delta E_1.
\]

After the earlier dimensional normalizations and the row transformation
`(m1,(m1-m4)/lambda_patch,m2,-m5,m3)`, this is precisely the already audited
matrix B. It has the inherited continuum inverse bound `||B^(-1)||<=1000`.
For any finite smooth normalized discrepancy z, choose `c=-B^(-1)z`.
There is no quadratic map and no small-root restriction for this coefficient.
The cost `lambda_patch^(-1)` in an arbitrary incoming discrepancy remains.

Independent integrand quadrature checks this linear map for four patch
parameters, including `lambda_patch=1e-50`. Reusing the nonlinear leading
repair produces a nonzero error. Quadrature is explicitly a diagnostic; the
continuum invertibility bound is inherited, not inferred from those samples.

This checkpoint has not verified the actual I3 support/memory hypotheses,
integrated the actual inner coefficient through an annular cutoff, or derived
the complete higher stress tails. Thus it establishes the compatible local
linear inverse, **not a finished global extension**.

## 6. A protected-core lower bound for the first partial background

Choose `X*=1/Lambda`, `eta=0`, and `Delta=rho/16` in the Picard bound. Write
`M1=(z+beta) exp(beta)` with the actual `C_op=C^16` and `a=sqrt(3/Lambda)`.
At this point `F0=Phi(1,0)/C>.264/C`. Consequently for

\[
0<q\le q_0:=\min\left(1,
\left[\frac{.132}{CM_1}\right]^{1/(2h)}\right),
\]

the first partial background has

\[
F_0+q^{2h}F_1\ge F_0/2>.132/C,
\qquad
u_\theta\ge \frac{.132\sqrt{2/\Lambda}}{C}\,q^{-A}.
\]

Along `z=0`, `r=sqrt(2q/Lambda)`, `t=1-q`, this lower bound diverges as
q decreases to zero. The threshold q0 is finite and positive but may be
enormously small; it is kept symbolic. This establishes survival of the leading
swirl against **one inner background correction only**. It proves neither
survival against the infinite correction sum nor energy or force admissibility.

## Next open work

1. Extend the actual original moment inputs to angular order four, carry three
   slow derivatives through the loop and moving phase map, and propagate the
   final repair. Retain fixed-N radial and cutoff costs to obtain a valid global
   residual constant.
2. Extend this actual coefficient beyond the core, verify I3's support and
   leading memory conditions, restore all five moments, and derive the retained
   annular stress and exterior tails.
3. Recompute the full residual for that same extended field, including physical
   waves, quadratic terms and axial diffusion of the coefficient. Repeat at all
   orders and prove divergence-preserving summation before claiming smooth force
   or complete-field blow-up.

Reproduce with the pinned dependencies:

```sh
python research/navier_stokes_cascade/results/first_correction_construction/test_audit.py -v
python research/navier_stokes_cascade/results/first_correction_construction/audit.py \
  --verify-record research/navier_stokes_cascade/results/first_correction_construction/evidence.json \
  --output first_correction_reproduced.json
```

Passing tests reject identifiable algebra mistakes and preserve provenance. They
do not constitute independent mathematical review or a proof-assistant certificate.
