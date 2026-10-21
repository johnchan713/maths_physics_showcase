# Independent all-index ledger for the nonlinear core fixed-point map

Parent: `6582f51d627b46fdcf244335b880223fcb6b59e3`.
Status: `core-fixed-point-all-index-ledger-independently-rederived`.

The preceding checkpoint independently identified the common pressure datum
and strengthened the pressure-to-growth partition.  This checkpoint addresses
the next part of the foundation frontier: the infinite coefficient-space
operator and every monomial in the normalized core fixed-point map.  It
rederives the product, integrated-derivative and order-one inverse bounds,
enumerates all 21 normalized source terms, retains all pressure memory, and
recovers the invariant-ball and Lipschitz majorants used by the selected
`Lambda`.

No historical evidence is edited.  The fixed pressure datum, `h,j,sigma`,
`rho`, `Lambda`, `C`, and coefficient weights are unchanged.  This closes the
operator/majorant part of the exact-core contraction provider; it does not
independently reprove the later scalar-growth comparison, flat activation
collar, compact roots or global construction.  The exact controls are
regressions for the displayed proof, not a formal-proof certificate.

## 1. Complete coefficient space and the product constant

For `F=sum_n F_n(eta)Y^n`, use

\[
 a_{nk}=\frac{20^{-n}\rho^{-k}k!\binom{n+k}{k}}
 {(n+1)^2(k+1)^2},\qquad
 \|F\|_\rho=\sup_{n,k,\eta}\frac{|\partial_\eta^kF_n|}{a_{nk}}. \tag{1}
\]

A norm-Cauchy sequence converges uniformly coefficient by coefficient and
derivative by derivative.  The fundamental theorem of calculus identifies
successive limits as derivatives of one function, and the supremum in (1)
then gives convergence in norm.  Thus the space is complete; no finite Taylor
cutoff is used in the fixed-point theorem.

For radial split `i,n-i` and derivative split `j,k-j`, Leibniz' binomial and
the factorial in (1) cancel.  The remaining binomial factor satisfies

\[
 \binom{i+j}{j}\binom{n-i+k-j}{k-j}\le\binom{n+k}{k},       \tag{2}
\]

because the left side counts one class of `k`-subsets from two disjoint
blocks of total size `n+k`.  Splitting the squared-denominator convolution at
its midpoint gives

\[
 \sum_{i=0}^n\frac{(n+1)^2}{(i+1)^2(n-i+1)^2}
 \le8\sum_{m\ge1}m^{-2}<16.                                \tag{3}
\]

The last strict inequality follows without evaluating zeta(2): for `m>=2`,
`m^-2<1/[m(m-1)]`, whose sum telescopes to one.  Applying (3) in both indices
gives

\[
                         \|FG\|_\rho\le A_0\|F\|_\rho\|G\|_\rho,
 \qquad A_0=256.                                           \tag{4}
\]

The code samples finite convolutions only as a failure control; (2)--(4) are
the all-index argument.

## 2. Fixed coefficients and every integrated derivative

If a degree-zero coefficient is bounded by `b` on the full complex rectangle,
Cauchy's estimate on the half rectangle gives the normalized multiplier

\[
                         (k+1)^2/4^k\le1.                  \tag{5}
\]

The inequality starts with equality at `k=0,1`; its ratio thereafter is
`((k+2)/(k+1))^2/4<=9/16`.  Hence the fixed coefficient has norm at most `b`.
This is why a complex bound, not only a real-axis supremum, is required.

Let `J_nu` solve

\[
 YG_{YY}+\nu G_Y=F,\quad G(0)=0,\qquad
 (J_\nu F)_{n+1}=\frac{F_n}{(n+1)(n+\nu)}.                 \tag{6}
\]

The exact weight shifts are

\[
 \frac{a_{nk}}{a_{n+1,k}}
 =\frac{20(n+2)^2}{(n+1)(n+k+1)}\le80,
\]
\[
 \frac{a_{i,k+1}}{a_{i+1,k}}
 =\frac{20(i+2)^2(k+1)^2}{\rho(i+1)(k+2)^2}
 \le\frac{80}{\rho}(i+1).                                 \tag{7}
\]

For the worst mixed word `J_nu[(partial_eta F)(D_Y G)]`, assign the extra
output radial degree to the differentiated factor.  What remains after (7)
is

\[
 \frac{(i+1)j}{(i+j+1)(i+j+\nu)}\le1.                     \tag{8}
\]

This includes `i=0` or `j=0`.  Omitting either derivative only improves the
bound, and averaging contributes a divisor.  Therefore every integrated
source word in the map is controlled by

\[
 D_0 A_0^m b\prod_r\|F_r\|_\rho,qquad D_0=80/\rho,        \tag{9}
\]

when it has `m` unknown factors and fixed coefficient norm `b`.  Bare
parameter differentiation is not asserted to be bounded.

## 3. The order-one angular inverse

Put `T=J_2(chi multiplication)/2`.  Multiplication by `chi`, whose fixed
coefficient norm is below three, preserves the lowest radial degree; `J_2`
raises it.  Starting at degree zero, repeated words obey

\[
 \|T^k\|\le\frac{30720^k}{k!(k+1)!}.                       \tag{10}
\]

Indeed the exact successive ratio is `30720/[(k+1)(k+2)]`.  The series is
absolutely summable and is bounded by `exp(30720)` after discarding the extra
`(k+1)!` denominator.  Its partial sums satisfy the two-sided geometric
identity, while the factorial remainder tends to zero.  Thus

\[
 (1+T)^{-1}=\sum_{k\ge0}(-T)^k,qquad
 \|(1+T)^{-1}\|\le V:=e^{30720}.                           \tag{11}
\]

No smallness assumption on `chi` is used.

## 4. Complete normalized source ledger

Use the same normalized variables

\[
 U=U_*+\Lambda^{-1}u,\quad
 W=W_*+\Lambda^{-1}B,\quad H_c=H_*+\Lambda^{-1}du,
\]
\[
 B=-2D\eta\,\bar u-d\,\partial_\eta\bar u,qquad
 p=I(g^2\Phi^2).                                           \tag{12}
\]

Direct substitution gives the retained sources

\[
 LR_1=[W+h(1-2\eta U)+du\zeta_*]\Phi
       +W D_Y\Phi+H_c\Phi_\eta,                            \tag{13}
\]
\[
\begin{aligned}
 LR_2={}&[A(1-4\eta U_*)+dU_{*\eta}]u-2A\eta\Lambda^{-1}u^2\\
 &+W D_Yu+H_*u_\eta+d\Lambda^{-1}uu_\eta\\
 &-4A\eta p+dp_\eta-2\eta D_Yp .
\end{aligned}                                              \tag{14}
\]

The independent term ledger contains 21 entries:

| Group | Terms after `1/L<=2` | Sum |
|---|---|---:|
| `R1` fixed linear | `W* Phi`, seed-`h` `Phi`, `W* D Phi`, `H* Phi_eta` | `19237/50 < 386` |
| `R1` ordinary quadratic | two `B Phi`, corrected-`U` seed, two `B D Phi`, `d u Phi_eta` | `952/25 < 40` |
| `R1` zeta quadratic | `d u zeta* Phi` | `10 M` |
| `R2` fixed linear | seed coefficient times `u`, `W* D u`, `H* u_eta` | `442` |
| `R2` ordinary quadratic | `-2A eta u^2`, two `B D u`, `d u u_eta` | `32` |
| `R2` pressure | `-4A eta p`, `d p_eta`, `-2 eta D p` | `16+10+8=34` |

The ledger retains the average in both `B` terms, every parameter derivative,
and all three pressure terms.  Applying (9) gives, for `||Phi||,||u||<=S`,

\[
 M_1=D_0[386A_0S+(40+10M)A_0^2S^2],                       \tag{15}
\]
\[
 M_2=D_0[442A_0S+(32A_0^2+2720A_0^4)S^2].                 \tag{16}
\]

## 5. Pressure count and Lipschitz count

The pressure primitive is the largest place an algebra power could be lost.
There are four factors in `g*g*Phi*Phi`, hence three product multiplications,
followed by radial integration:

\[
                         \|p\|\le80A_0^3S^2.              \tag{17}
\]

Each pressure source in (14) then costs one more fixed-coefficient product.
The exact sum is therefore

\[
                  34\cdot80 A_0^4=2720A_0^4,              \tag{18}
\]

not `2720 A0^3`.  For two states in the ball,

\[
 \Phi^2-\widetilde\Phi^2
 = (\Phi-\widetilde\Phi)(\Phi+\widetilde\Phi),
\]

so the pressure difference carries `2S`.  Its Lipschitz coefficient is

\[
                         5440A_0^4S.                       \tag{19}
\]

The other quadratic words have the same replace-one-factor rule.  Thus

\[
 M_1'=D_0[386A_0+2(40+10M)A_0^2S],                        \tag{20}
\]
\[
 M_2'=D_0[442A_0+2(32A_0^2+2720A_0^4)S].                  \tag{21}
\]

Tests that remove any one pressure row, lower the algebra power, or omit the
factor two fail the frozen exact ledger.

## 6. Invariant ball and strict contraction

The center is

\[
 \Phi_0=(1+T)^{-1}1,qquad u_0=-YZ_*/(2L).                 \tag{22}
\]

Because the coefficient of `Y` has reciprocal weight 80,
`||u_0||<=40 Z_b`; (11) gives `||Phi_0||<=V`.  Hence the closed unit ball
about the center lies inside the componentwise radius

\[
                         S=2+V+40Z_b.                      \tag{23}
\]

With `F_b=V(M_1+M_2)` and `L_b=V(M_1'+M_2')`, the exact map has displacement
at most `F_b/Lambda` and Lipschitz factor at most `L_b/Lambda`; discarding its
factors `1/2` is conservative.

Let `E=10000 rho^-6`, `epsilon=10^-16`, and write `H=K^(1/16)`.  The selected
formula is

\[
 H=10^6E(1+S+F_b+L_b+M+Z_b+P^2)>1,
\]
\[
 \Lambda=10^{12}H^{128}\epsilon^{-2}.                     \tag{24}
\]

Therefore

\[
 \frac{E\max(S,F_b,L_b)}{\Lambda}
 \le\frac{\epsilon^2}{10^{18}H^{127}}
 <10^{-50}.                                                \tag{25}
\]

Since `E>1`, (25) proves both ball invariance and strict contraction.  Banach's
theorem gives a unique infinite coefficient-space fixed point.  This is an
analytical conclusion from (1)--(25), not from the finite tests.

## 7. Result and remaining frontier

The all-index operator estimates, complete normalized monomial list, pressure
algebra power, Lipschitz factor, invariant ball and selected-parameter
contraction have now been independently rederived.  No omitted term or false
constant was found; the pressure count `2720 A0^4` is indispensable and is
retained.

This checkpoint does not promote the entire foundation chain.  Next one must
independently validate the scalar angular comparison and positive source,
then the complementary axial exit sector and flat activation collar with all
five pressure moments.  Continuation, matching roots, compact repair,
annular/exterior attachment, regional cone, global pressure, weighted edges,
coefficients, source cycles and the complete residual remain later
obligations.  No unconditional pulse, residual improvement, force theorem or
blow-up is claimed.

```sh
python3 research/navier_stokes_cascade/results/core_fixed_point_ledger_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_fixed_point_ledger_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_fixed_point_ledger_recheck/evidence.json \
  --output /tmp/core_fixed_point_ledger_recheck_fresh.json
```
