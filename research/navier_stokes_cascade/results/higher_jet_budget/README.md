# Higher jets of the finalized compact background

This checkpoint supplies a **conditional higher-derivative budget for the actual
modulated and repaired background on the selected compact annulus**. It derives
the original incoming C4 transfer, a larger mixed input envelope, third-order
loop/inverse-phase bounds and a C3 bound for the same final moment repair.
The previously selected finite N is unchanged. A valid annular profile constant
is now `K=H^64`, conditional on the inherited analytic core and leading-profile
estimates.

This is a mathematical argument with finite diagnostic checks, not independent
review or a proof-assistant certificate. It does not solve the global first
correction, construct the physical wave tensor, sum all corrections, prove
smooth forcing, or prove complete-field blow-up. The explicit background
residual bounds below are still singular as q tends to zero.

Parent: `007e62724850fbe98707d5c50578e47b35768c4d`. Source: the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.9), (4.16), (B.22), (B.25)-(B.26), (C.5), (C.9)-(C.16).
The B.22 page was visually checked. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The source formulas motivate the construction; the explicit new estimates
are derived here. Earlier audit folders and GOAL.md remain frozen.

## 1. Domain, norms and unchanged parameter order

Keep all actual supports, axis pressure, profile data and parameters. Write

\[
A=C^{4096},\quad G=e^{A^{32}},\quad H=e^{G^4},\quad
\delta_L=e^{-A^{20}},\quad N=1+\lfloor H^{32}\rfloor.
\]

The radial activation width is the different quantity `delta=C^(-20)`.
Use factorial angular norms and mixed logarithmic-radial norms

\[
\|f\|_{C^k_\eta}=\sum_{s=0}^k\sup|\partial_\eta^sf|/s!,\qquad
\|f\|_{\mathcal J_k}=\sum_{r+s\le k}\sup|\partial_y^r\partial_\eta^sf|/(r!s!),
\quad y=\log X.
\]

They are product algebras. A derivative maps Jk to J(k-1) with norm at most k.
For periodic functions, joint jets additionally include phase derivatives,
with total order at most k. Angular endpoints are included one-sided.

The final jet budget applies on

\[
\mathcal R=[X_a e^{t_c/2},\rho_p e^5]\times[-1,1].
\]

This contains the entire modulation and first reserved repair patch. Earlier
geometry gives `A^(-1)<X<A`, radial length below A, and the original value
bounds and positive reciprocals below A. The original coordinates satisfy
`a>=1/A`, `2<v<=A`, and the variance construction has `mu,|mu*p2|,|t|<G`.
Those value bounds are unchanged.

We introduce only a larger **derivative envelope**

\[
A_{hi}=A^{256}.
\]

It does not replace A in G,H, the loop cap or N. This distinction prevents a
circular frequency selection. No enormous scale or actual high-frequency
field is materialized.

## 2. Original C4 transfer: derive the missing reference estimate

The exact core coefficient norm gives, for `0<=Y<=5`,

\[
\|\Phi\|_{C^4_\eta},\|u\|_{C^4_\eta}
\le\frac{14624}{6075}S\rho^{-4}<D_4:=3S\rho^{-4}.
\]

The corresponding generating sum through order five is also below three.
The inherited coefficient norm, radial Cauchy disks and the real lower bound
`Phi>.264` bound the native logarithm, reciprocal and required mixed jets
through order six by

\[
N_{nat}=10^{20}(1+S)^8(1+\rho^{-1})^{48}<K_{old}.
\]

Here `K_old` is the frozen core parameter, not the new physical K. To obtain
this deliberately generous bound, differentiate Phi and u by their coefficient
norm; use real logarithm/reciprocal formulas with denominator at least .264.
At most eight factors and factorials through order eight occur. The radial
disks have fixed radius, while the displayed rho power dominates all angular
costs. This extends the finite derivative bookkeeping, not the existence claim.
Separately Cauchy gives `||zeta_*||C4<K_old`; the positive-mixture axis pressure
has C5 norm below `1e5 P^2<K_old`.

For the reference continuation, put
`Rr=phi_r/phi_*`. Its defining slopes contain only radial cutoffs. Native
`Dy log Phi` is integrated for logarithmic length at most `2 delta`, so its
C5 integral norm is at most `200 K_old delta<1`. Thus the native exit values
and the same reciprocal formulas give

\[
\|R_r\|_{C^5_\eta},\|R_r^{-1}\|_{C^5_\eta}<3K_{old},\qquad
\|U_r-U_*\|_{C^5_\eta}<1000K_{old}/\Lambda.
\]

Radial averaging is a contraction in every angular norm. Use the original
source formulas, with `W=1-2D eta M/X-d M_eta/X`, `Hc=D eta+d U`,
`d=1-eta^2`, `L=1-2h eta^2`, and `lr=1+Dy log phi_r`:

\[
S_q=-Wl_r-h(1-2\eta U)-H_c(\log E)_\eta,
\]
\[
S_n=-WU_y-A_{phys}(1-2\eta U)U-H_cU_\eta
-d\Pi_\eta+4A_{phys}\eta\Pi+\eta E^2.
\]

Since `||L^(-1)||C4<2`, direct product bounds give
`||Sq,r||C4<=1e6 Lambda K_old^4` and
`||Sn,r||C4<=1e6 K_old^4`. Pressure uses the same axis datum and the regular
integrand `phi_r^2/C^2`; its C5 increment is below
`110 C^(-2) exp(2B)<1`. All fifth angular inputs used here have been supplied;
no C3 estimate is relabelled C4.

The equations for the integrated reference coefficients have exact solutions

\[
p_{1,r}(X)=\frac{X_aR_r(X_a)p_{1,r}(X_a)
+\int_{X_a}^X xR_r(x)S_{q,r}(x)/L\,dx}{XR_r(X)},
\]
\[
n_{s,r}(X)=\frac{X_an_{s,r}(X_a)+\int_{X_a}^X S_{n,r}(x)/L\,dx}{X}.
\]

The large factor phi_* cancels in the first formula. Bounding the positive
ratios and integrals, rather than exponentiating a large Gronwall coefficient,
yields safe bounds

\[
\boxed{\|p_{1,r}(X)\|_{C^4_\eta}\le10^{18}X\Lambda K_{old}^{12},\qquad
\|n_{s,r}(X)\|_{C^4_\eta}\le10^8K_{old}^8.}
\]

Both fit the existing `B/100` budget. For example the initial term is at most
`18 K_old^3`, and the source term at most `9e6 X Lambda K_old^6`; the
displayed bound dominates both since `X Lambda>=4`.

This cancellation does not remove the axis gradient from Sq:
`(log E)_eta=Lambda*zeta_*+(log Rr)_eta`. An independent manufactured check
compares the original cumulative-moment formula with this source integral
through four angular derivatives; dropping the gradient or the X factor fails.

The actual activation integrates these reference functions. Its C4 axial
change is at most `1e10 K_old^8 delta`, below `epsilon/100`. The axis factor
has

\[
\|\log\phi_*\|_{C^4_\eta}\le
\Lambda M(2+r_0^{-1}+4/(3r_0^2)+2/r_0^3).
\]

Together with `||log Phi||C4<=1000 D4^4`, radial slope integration gives
`||log(CE)||C4<B` on the continuation. The fifth-order reference version adds
`16 Lambda M/(5 r0^4)` and also fits B. The shape transition has
`||log f||C4<=682` for `f=(1+eta^2)^(-1)`, hence still
`||log(CE)||C4<=204B`.

Carry the five original densities through exactly the same angular algebra.
The pressure integral retains its integrable axis factor. Normalizations now
cost

\[
\|f^{-1}\|_{C^4}=5,\qquad \|f^{-2}\|_{C^4}=25,
\]

not the old C3 value 24. The energy-row coefficient is `25(101+160+64)=8125`,
still below the old generous `1e6` cap. The same root contraction in the C4
function algebra gives

\[
\|z\|_{C^4}<2.86\,10^{-17},\quad
\|c\|_{C^4}<5.72\,10^{-14},\quad
\sup|c^{(4)}|<1.3728\,10^{-12}.
\]

Uniqueness identifies this with the existing joining correction. Its exact
raw fourth differentiation contains eight lower terms with coefficients
`8,6,8,24,12,12,8,1`; the diagnostic checks each one independently.

## 3. Mixed inputs, including radial widths

The unchanged step satisfies `sup|sigma^(j)|<1e100`, `0<=j<=5`. Here is a
uniform bound, including its flat edges. On `[1/4,3/4]`, the derivatives of
`x^(-2)-(1-x)^(-2)` through order five are below `1e8`; logistic derivatives
are at most `2^j j!`. The fifth Bell number is 52, giving a bound below
`1e4*52*(1e8)^5`. Near zero expand the logistic denominator as a geometric
series. Its differentiated terms are bounded by

\[
52(1440n)^5x^{-15}e^{n[-x^{-2}+(1-x)^{-2}]}.
\]

The sum of `n^5 exp((n-1)a)` is below two because `a<=-14`. The remaining
bound is `104*1440^5*4^15*exp(16/9-16)`, also below `1e100`. Reflection
covers the other edge. Thus no uniform edge analyticity is assumed.

This loose fixed constant costs less than C12. Every radial derivative of a
width-delta cutoff still costs delta^(-1)=C20; a width-.1 bump of radial order
four includes its inverse fifth power. The actual angular C4 joining bound
from section 2 is used when differentiating the earlier joining bumps.

The resulting regional bookkeeping is:

| Region or input | Safe power of C |
|---|---:|
| Native nonzero log jets through order four | 4 |
| Reference field E mixed J2; U mixed J3 | 64; 64 |
| Reference source Sq,Sn mixed J2 | 110; 160 |
| Reference p1,ns mixed J3, from their radial equations | 128; 192 |
| Actual continuation nonzero log-E jets through order four | 256 |
| Actual continuation U mixed J4 | 512 |
| Shape, earlier join and outer prefix fields mixed J4 | 128 |
| All original fields mixed J4 | 2048 |
| All original cumulative moments mixed J4 | 8192 |

For the reference rows, two radial differentiations of Sq,Sn use at most
three field derivatives and the C5 angular bounds above; their cutoff costs
are at most C40, before the fixed cutoff constant. More precisely, W's J2
norm sees U_eta: two radial derivatives of U_eta differentiate its radial
cutoff only once. The bounds `W:C36`, `lr:C56`, `Hc:C36`,
`(log E)_eta:C36`, `U_y:C64`, `Pi_eta:C132` give source exponents at most
93 and 134, below the displayed 110 and 160. Pure angular bounds
`p1:C16`, `ns:C9`, `lr:C2` and the radial equations then give p1_y, p1_yy,
p1_yyy exponents 19,113,116 and ns_yyy exponent 163. This supplies the
128 and 192 padding without assigning every mixed derivative the worst
pure-radial width loss.

For the actual continuation, `F_y=1/2-kappa p1,r/2` and the final convex
interpolation need three cutoff derivatives, at most C75. The axial taper
can add a second such factor; `75+75+192+2<512`.

When passing from F=log E to E, retain the actual value bound E<1 on the
continuation. Its fourth Taylor jet is `E * sum_{j=0}^4 (delta F)^j/j!`,
not `exp(||F||J4)`. This gives a power below C2048. The same argument uses
E<2P on the outer prefix. Moments use their original integral formulas;
the largest product involves two field jets and the geometric X factors,
well below C8192. Pure angular pressure integration keeps `E^2/X` regular
at the axis; radial moment derivatives use their local density identities.

For a nonvanishing f with `|f|>=C^(-v)` and `||f||Jk<=C^s`, truncated
reciprocal Taylor expansion gives

\[
\|f^{-1}\|_{J_k}\le C^v\sum_{j=0}^k(C^{v+s})^j
<C^{v+k(v+s)+1}.
\]

Use the original floors `E>=C^(-4)` and `a>=C^(-22)`. The exact positive
exponent ledger propagates the fields, moments and reciprocals through the
original pressure/shear coordinates. Its largest listed exponent is below
`4096*256`. Consequently every required original J4 field/moment and J3
stress/loop input is bounded by A_hi. This does not change their sharper
original value bounds A.

## 4. Third-order loop and moving phase estimates

Write `q(z,theta)=exp(z sin(theta))/I0(z)` and `z=mu p`. Bounds on the
derivatives of I0 and its reciprocal give joint `(z,theta)` derivatives
through order four below `exp(4G)` whenever `|z|<G`. Define

\[
R(z)=\left\langle\left[\int_0^1q_z(tz,\theta)dt\right]^2\right\rangle.
\]

The positive Bessel numerator series gives `R>=exp(-2G)/2`; differentiating
the integral gives derivatives through order three below `exp(9G)`. Thus
`sqrt(R)` and its first three derivatives are below `exp(40G)`. Use the
regular signed root

\[
F(\mu,p)=d_0\mu\sqrt{R(\mu p)}.
\]

It has partial derivatives through order three below `exp(50G)`. The tilted
variance lower bound gives `V_mu>=d0^2 mu exp(-4G)`. Since
`F<=d0 mu exp(5G)`, it follows that `1/F_mu<exp(10G)`, including the exact
limit `F_mu(0,p)=d0/sqrt(2)`. Dividing by V_mu at zero would be invalid.

The cutoff target is `s=zeta(vs)*sqrt((2+delta_L/2-vs)/a)`. Its derivatives
through order three, and the needed input jets, are below `G^64<exp(G)`.
One template covers the algebraic and cutoff costs:

\[
10^{150}A^{10000}\delta_L^{-10}<G.
\]

The endpoint check and decreasing ratios prove this for every `A>=2^16`.
Differentiate `F(mu(x),p(x))=s(x)`, including all p derivatives. For example

\[
F_\mu\mu''=s''-F_pp''-F_{\mu\mu}(\mu')^2
-2F_{\mu p}\mu'p'-F_{pp}(p')^2.
\]

At order three the terms include `3 D2F[(mu',p'),(mu'',p'')]` and
`D3F[(mu',p')]^3`. The resulting safe bounds are `exp(70G)`, `exp(220G)`
and `exp(500G)` for the first, second and third derivatives of mu.

Use the removable quotient in t as an integral, not division by p. Joint
derivatives through order three of t at fixed theta are below `exp(2000G)`.
For the phase lift, `w=a(1+t^2)/(2pi v)` still satisfies
`G^(-3)<w<G^3`. Its third-order jets and those of the forward lift are below
`exp(5000G)` and `exp(5001G)` respectively. Differentiating the inverse
identity, including terms where derivatives hit w, gives inverse-phase
derivatives through order three below `exp(40000G)`. The lift value itself is
bounded only on one fundamental period; its positive-order derivatives and
all pulled-back periodic functions extend with uniform bounds to every phase.

Composition now gives t at fixed phase below `exp(150000G)` in joint J3.
The rational loop factors `1/(1+t^2)` and `t/(1+t^2)` have bounded scalar
derivatives through order three. Product and composition rules bound the
loop jets by `exp(500000G)`. Zero-mean periodic integration commutes with
slow differentiation; phase derivatives recover the loop sources. Therefore

\[
\boxed{\|\mathcal A\|_{J_3(y,\eta,\phi)},
\|\mathcal B\|_{J_3(y,\eta,\phi)}<e^{10^6G}<H.}
\]

All counts are finite, uniform in phase and include the zero-pressure and
zero-variance boundaries. They concern the same loop as before.

## 5. Same final repair, now in angular C3

Since phase Ny is independent of eta, the modulation has angular C3 field
error `d3<=4 A_hi H/N`. The five original density differences, with the
fixed geometric bounds, give moment/pressure error
`m3<=32 A^3 A_hi^2 H/N`. The normalization operator has C3 norm at most
`4 A_hi^2`, retaining the inverse-lambda cost. Hence

\[
\|z_N\|_{C^3}\le128A^3A_{hi}^4H/N<H^{16}/N.
\]

On the final leading repair patch, scaling both bumps by the original e(eta)
makes the matrix B and quadratic map Q independent of eta. The original
continuum bounds `||B^(-1)||<=1000`, `||Q||<=1e5` therefore apply in every
factorial function algebra. The same small-root contraction gives

\[
\|c_N\|_{C^3}\le2000H^{16}/N.
\]

Uniqueness identifies this with the already accepted C1 root. Its support,
field positivity, state cone and N are unchanged. Radial bump derivatives
retain their inverse widths. Propagating the differences through partial
moment integrals, including the incoming memory, yields

\[
\|\widetilde E-E\|_{C^3_\eta},\|\widetilde U-U\|_{C^3_\eta}<H^{18}/N,
\qquad
\|\Delta(M,I,J,S,C_p,\Pi)\|_{C^3_\eta}<H^{20}/N.
\]

Both are below one at the fixed N. Final angular C3 field and moment norms
are thus below `2 A_hi`. No claim of small radial derivative errors follows.

## 6. Fixed-frequency radial losses and an explicit annular K

For beta=script-A(y,eta,Ny)/N,

\[
\beta_y=\mathcal A_y/N+\mathcal A_\phi,\quad
\beta_{yy}=\mathcal A_{yy}/N+2\mathcal A_{y\phi}+N\mathcal A_{\phi\phi}.
\]

The last term grows with N. Retain `|beta|<=H/N` when bounding the exponential
and keep its squared first derivatives. Since `N>=H`, the modulated field J2
norm is below `32 A_hi N H<H^2 N`. The repairing bumps add less than one
in mixed J2; use `H^3 N` for the finalized fields.

Convert logarithmic to ordinary derivatives using
`f_X=f_y/X`, `f_XX=(f_yy-f_y)/X^2`. On the annulus X and its inverse are
below A. Set `F=E/sqrt(2X)` and use the exact radial flux

\[
V=\frac{2\eta XU-2D\eta M-dM_\eta}{L}.
\]

Its ordinary C2 jet uses the now-bounded third angular derivative of M.
Other differentiated primitives use `M_X=U`, so no extra unseen derivative
is assumed. Pressure uses `Pi_X=E^2/(2X)`; use the angular value bound on
one E factor to retain a linear N bound in Pi_XX. Direct product estimates
give ordinary mixed C2 jets of F,U,Pi,V below `H^8 N`. Since
`H^32<N<2H^32`, every required jet and geometric reciprocal is bounded by

\[
\boxed{K=H^{64}.}
\]

Apply the already derived physical background-residual inequalities, for
viscosity one and `0<q<=1`:

\[
|R_r|\le42H^{448}q^{-3/2},\quad
|R_\theta|\le440H^{448}q^{-A_{phys}-1},\quad
|R_z|\le148H^{384}q^{-A_{phys}-1}.
\]

These instantiate the formerly missing compact annular derivative constant.
They are **large singular upper bounds**, not a residual decay or improvement
result. The physical waves and their stresses are not included. Neither an
estimate outside this annulus nor a bound for the fully corrected field is
asserted here.

## 7. Checks and next obligation

The interval checks run at 80 and 110 digits. Finite manufactured diagnostics
independently check the regular variance formula, third implicit root jets,
third inverse-phase composition, mixed fixed-N chain rules, exponential
product terms, nonlinear C3 repair and all five fourth-angular moment
identities. The moment integrals use exact rational polynomial coefficients;
the older binary64 polynomial checker is separately compared to roundoff.
Counterexamples retain why C3 alone cannot bound C4 and why radial second
derivatives need not shrink with N. None samples the actual high-frequency
profile or replaces a continuum estimate with quadrature.

Next: extend the already constructed inner first coefficient through the
annulus, verify the actual I3 support/memory conditions, restore its five
linear compatibility moments, and derive its higher stress tails. Recompute
the residual for that same field including the physical wave tensor. All-order
solvability, summation, smooth compact forcing, energy and a complete-field
blow-up lower bound remain separate obligations.

```sh
python research/navier_stokes_cascade/results/higher_jet_budget/test_audit.py -v
python research/navier_stokes_cascade/results/higher_jet_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/higher_jet_budget/evidence.json \
  --output higher_jet_budget_reproduced.json
```
