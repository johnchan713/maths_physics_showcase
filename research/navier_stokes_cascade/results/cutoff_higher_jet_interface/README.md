# Higher cutoff jets: exact repair propagation and the physical derivative boundary

Parent: `f84076cef9550e7bd146493ce9a5a4f7afdc5b3d`.
Status: `all-finite-repair-bump-jets-and-higher-cutoff-input-boundary`.

This checkpoint supplies **all finite derivative orders of the five actual
repair bumps**, an explicit interface from cutoff/collar jets to repaired
coefficient and physical radial jets, and a precise limit of the preceding
two-shape-norm cone/scale certificate. Even a strictly accepted low-order
certificate cannot uniformly control third velocity or second stress jets
over its admissible cutoff class. Exact moment repair does not remove the
local dependence on the missing higher cutoff derivatives.

This is NOT a claim that the fixed smooth cutoff has infinite derivatives,
that the original all-jet Borel prescription omitted them, or that an actual
Navier--Stokes pulse fails. The prescription already includes order-dependent
seminorms. The point is that its constants cannot all be replaced by the
recent value-gate Cstar. The actual higher cutoff jets still need bounds.
No field, cutoff, width, scale, root, tilt, support or pressure is changed.
Historical evidence is unchanged. Analytical implications and auxiliary
algebra controls below are distinguished from actual-profile verification.

## 1. Propagate every fixed radial order with the correct Jacobians

Keep the CORRECTED collar `tc=C^-20/(10 sqrt(log C))`, `Delta=tc/16`,
`s=log(X/Xkeep)/Delta`, and the unchanged smooth shape K. Let
`L_j >= sup |K^(j)|`, `L_0=1`, and let `J_j` bound the inner F1,U1 radial
jets through order m in any fixed product-algebra C_eta^a norm. These are
inputs for the SAME finalized leading/collar solution, not radial Cauchy
estimates across the nonanalytic activation.

Writing D=partial_y, y=log(X/Xkeep), the exact formula is

\[
 \partial_X^m[X^{-\alpha}b(y)]
 =X^{-m-\alpha}\prod_{i=0}^{m-1}(D-\alpha-i)b(y).
 \tag{1}
\]

For alpha=0 the coefficients are the signed first-kind Stirling numbers
`s(m,j)`; their recurrence is `s(m+1,j)=s(m,j-1)-m*s(m,j)`.
In particular the THIRD cutoff derivative is

\[
 \kappa_{XXX}=X^{-3}
  (\Delta^{-3}K'''-3\Delta^{-2}K''+2\Delta^{-1}K').       \tag{2}
\]

Valid bounds are `k_0=1` and

\[
 k_m=Xkeep^{-m}\sum_{j=1}^m|s(m,j)|\Delta^{-j}L_j,
 \qquad P_m=\sum_{j=0}^m {m\choose j}k_jJ_{m-j}.        \tag{3}
\]

For a retained logarithmic bump beta, let `B_j` enclose its y derivatives,
and let `p` be a positive lower X bound on all bump supports with
`0<rp<=sqrt(2p)`. Expand the polynomial in (1) with alpha=0 for a U bump
and alpha=1/2 for an F bump (`E/sqrt(2X)`). The corresponding budgets are

\[
 U_m=p^{-m}\sum_j|[D^j]\prod_{i<m}(D-i)|B_j,
 \quad F_m={p^{-m}\over rp}
       \sum_j|[D^j]\prod_{i<m}(D-1/2-i)|B_j.             \tag{4}
\]

If Cp bounds the FULL repair coefficients in C_eta^a, the finalized common
jet bound is

\[
 S_m=P_m+C_p\max(U_m,F_m).                              \tag{5}
\]

All five supports are disjoint, so there is no factor five. The angular
derivatives of the nonconstant moment normalizers belong to Cp, not to the
eta-independent bump. Equations (3)--(5) reproduce the earlier S,Sx,Sxx
exactly through order two with the same inputs. At higher order they require
the corresponding J_j, L_j and angular inverse bounds; these are not inferred
from the first-order constants. For positive coefficient orders after one,
the full lower-order source/moment ledger must likewise be retained.

## 2. Bound every derivative of the ACTUAL five bumps

The retained shape is
`sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2))`, extended flatly,
and `beta(y)=10 sigma'(10(y-c)+1/2)`. The earlier analytical bounds
`D_1=32,D_2=1216,D_3=220000` remain valid; none is overwritten.

Here is an explicit finite-order extension, without sampling a maximum.
On `0<s<=1/2` put `t=(1-s)^-2-s^-2`, `z=s^-2>=4`. Then
`sigma=logistic(t)<=exp(4-z)` and
`|t^(j)|<=2(j+1)! s^(-j-2)`. Let `P_k(v)` be the kth logistic derivative,
with `P_0(v)=v`, `P_(k+1)=v(1-v)P_k'`, and let `d_k` be the sum of
absolute coefficients of `P_k(v)/v`. Thus `|P_k(sigma)|<=d_k exp(4-z)`.

For every integer partition `(k_1,...,k_m)` with `sum j*k_j=m`, set
`k=sum k_j`, `p=ceil(m/2+k)`. Faà di Bruno and `z>=4` give the finite
integer bound

\[
 D_m^{gen}=\sum_{\sum j k_j=m}
 {m!\,d_k\over\prod_j k_j!}
 \prod_j[2(j+1)]^{k_j}\;\max(4,p)^p.                   \tag{6}
\]

Indeed `z^(m/2+k)<=z^p`; the maximum of `z^p exp(4-z)` on `z>=4` is
at 4 if p<=4, and at p otherwise, and is at most `max(4,p)^p`.
Symmetry covers the other half interval; exponential flatness gives zero
endpoint derivatives. Set `D_0=1`, take the smaller of (6) and a previously
proved D_m when available, and obtain, for EVERY finite m,

\[
 \boxed{B_m=10^{m+1}D_{m+1}.}                           \tag{7}
\]

This discharges bump-shape inputs in (4), not the unknown cutoff-shape or
inner solution jets. The formula can be very conservative. The constants
are actual-shape bounds with the unchanged width 1/10 and unit mass in dy.

## 3. Pull back physical radial jets without hiding scale losses

At fixed physical z,t, q and eta are independent of r and
`X=r^2/(2q)`, `R=r/sqrt(q)`. For a scalar component `H(q) f(X,eta)`,
where H is independent of r, the exact chain rule is

\[
 \partial_r^m[Hf]=Hq^{-m/2}
 \sum_{k=\lceil m/2\rceil}^{m}
 {m!R^{2k-m}\over2^{m-k}(m-k)!(2k-m)!}\partial_X^k f.
 \tag{8}
\]

For m=0 this is f. Hence (5), an upper normalized radius Rhi and (8)
give a finite explicit normalized radial bound with the factor
`|H| q^(-m/2)` RETAINED. There are no negative R powers in (8).
This statement is for scalar/cylindrical radial components at fixed z,t;
it does not purport to bound angular Cartesian derivatives, axial/time
chart coefficients, carrier derivatives, or all vector-potential curls.
Their existing product rules and extra cutoff-commutator derivative remain
separate. In particular it is not a complete physical-pulse theorem.
For the angular velocity, use f=E1=sqrt(2X)F1 in (8), with that extra
product differentiated; it is not legitimate to treat the radial factor
as H(q). On a fixed positive annulus this smooth nonzero multiplier does
not remove the third-jet principal term in the control below.

## 4. A low-order accepted gate still cannot supply L3

This continuum countercontrol is stronger than shrinking the transition,
which already made L1,L2 unbounded. Now they remain UNIFORMLY bounded.
Use an AUXILIARY `K0=1-sigma(s)` and

`psi(s)=sigma(8(s-1/4))*sigma(8(3/4-s))`.

The two transition supports are disjoint; psi is zero outside (1/4,3/4),
one on [3/8,5/8], and has derivative bounds P1=256,P2=77824. On its support,
`m0=1/(2*3^16)` bounds K0 and 1-K0 below: each exponential is at least
exp(-16)>3^-16 and the denominator is at most two. Also
`-K0'=sigma(1-sigma)t'>=16m0`, since t'>=32.
Set `a=m0/[4(1+P1+P2)]>0` and

\[
 K_n=K_0+{a\over n^2}\psi(s)\sin(n(s-1/2)),\quad n\ge1.\tag{9}
\]

All K_n are smooth, monotone, in [0,1], with the SAME flat endpoints and
allowed transition window. The following estimates hold for every n:

\[
 \|K_n-K_0\|_\infty\le a/n^2,\quad
 \|K_n'-K_0'\|_\infty\le a(n^{-1}+P_1n^{-2}),
\]
\[
 \|K_n''-K_0''\|_\infty
       \le a(1+2P_1/n+P_2/n^2).
 \tag{10}
\]

Thus `||K_n'||<33`, `||K_n''||<1217` uniformly. At s=1/2 the perturbation
has jets `(0,a/n,0,-an)`, whereas `K0'''(1/2)=3712`. Hence L3 is unbounded.
The continuum proof is (9)--(10), not the six finite algebra controls.

More generally ANY fixed smooth cutoff with STRICT norm-budget headroom
has a subinterval with values strictly between zero and one and slope
strictly negative (its total drop is one). A supported perturbation of the
form (9), with a sufficiently small fixed a and n sufficiently large,
preserves its M,N budgets and monotonicity on that subinterval, while L3
is unbounded. Outside that subinterval it leaves the base cutoff unchanged;
global monotonicity is preserved if the base cutoff was monotone. Thus a
strictly accepted fixed affine/mixed gate can coexist with this missing
higher input. A common manufactured affine budget 1218 accepts all the
explicit controls above; it is NOT the actual enormous-profile gate.

## 5. Exact moment repair and FULL pressure do not erase the local cost

At order one the leading source `w0=Omega0/X` is fixed. Each absolute
incoming moment still includes the original complete integrals through Xv,
in particular `W0/2,W1/2`; neither is truncated at Xcut. The difference
between two cutoff choices cancels that SAME source, rather than declaring
it zero. With `epsilon=a/n^2`, the five discrepancy differences are bounded
by (5) of `first_seminorm_bound` with seed `epsilon*J` and WITHOUT the
common fixed source term. The SAME linear inverse then gives
`||delta coefficients||_(C_eta^a)<=C*epsilon`. Exact ordinary moments and
the zero-axis forward pressure gauge continue to hold for each repaired
coefficient. Its forward pressure still has the complete original w0.
This comparison does not construct every higher-order background for K_n.
The difference envelope follows directly from the linear moment integrals;
it is not obtained by subtracting inequalities for two absolute upper bounds.

All repair bumps lie strictly to the RIGHT of the cutoff collar. Therefore
on that collar the finalized fields are EXACTLY
`F1=kappa F1in`, `U1=kappa U1in`; no repair coefficient cancels their local
jets. If, at a point s0 in such a perturbation plateau, the fixed inner
value `F1in(X0,eta0)` is nonzero, (2) and the product rule give

\[
 \delta F_{1,XXX}(X_0,\eta_0)
   =-{a n F_1^{in}(X_0,\eta_0)\over\Delta^3X_0^3}+O(1).
 \tag{11}
\]

For U1 the same statement uses its nonzero inner value. This is CONDITIONAL
nondegeneracy; no such actual-profile value is verified here. If both inner
coefficients vanish identically on the varying patch, this particular
control has no effect there and must not be claimed as an actual failure.

The exact first residual contains `-2nu(X F1XX+2F1X)` and
`-2nu(X U1XX+U1X)`. Its conservative stress primitives have viscous parts

\[
 T_{1\theta}^{vis}=2\nu X F_{1X},\qquad
 T_{1z}^{vis}={2\nu X\over\sqrt{2X}}U_{1X}.             \tag{12}
\]

Integration constants and the complete other residual terms remain; (12)
does NOT replace the full stress. The coefficients of F1XXX,U1XXX in its
second X derivative are `2nu X` and `nu sqrt(2X)`, respectively. All other
terms are uniformly bounded on the fixed collar for the C2-bounded family:
their residual integrands use at most first derivatives, their first
derivatives use at most second derivatives, and forward pressure/radial
primitives retain the same bounded full leading data. Thus exact repair
cannot cancel the indicated unbounded term when its inner value is nonzero.
Setting viscosity zero would remove this particular stress principal term,
but is forbidden for the positive-viscosity construction.

At any fixed q>0 with nonzero physical prefactor H, (8) gives
`partial_r^2(H T1)=H/q*(T1X+2X T1XX)`. A scale factor cannot uniformly
bound that family at a fixed scale. This is an input obstruction to a
derivative ESTIMATE, not a lower bound on the actual PDE residual, nor
proof that a complete summed field cannot have cancellations.

## 6. Scope, checks and next obligation

`jet_bounds.py` checks (1) and (8) independently through order eight,
compares the low-order interface to the old original calculator, derives
the actual bump bounds through order five, checks viscous primitive signs
and principal terms, and retains full-pressure difference controls.
Six rational oscillation controls keep a common accepted low-order bound
while their central third jets grow. They are auxiliary controls, not
samples of the unchanged leading/coefficient field or a continuum proof.

The precise additional obligation at radial order m is to supply the
unchanged L_1,...,L_m and J_0,...,J_m in (3)--(5), the complete corresponding
moment/inverse angular norms, and all other physical derivative inputs.
Second stress derivatives additionally require L3 unless a verified
same-profile structural cancellation or zero inner coefficient applies.
The original all-jet Borel schedule MUST retain these order-dependent
constants, as well as its extra curl-cutoff derivative. Actual low-order
cutoff/deficit/retained-scale acceptance, all physical pulses, complete
covariance including axial/pressure means, and full viscous interaction
residual remain unproved. No unconditional core/regional or blow-up claim.

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output /fresh/path.json`.
