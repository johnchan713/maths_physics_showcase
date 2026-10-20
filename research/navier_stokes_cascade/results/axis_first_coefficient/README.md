# First-coefficient axis jets from the analytic core and a collar join

Parent: `d93abfc5528b2b6b4826144973ba5aeed574eed9`.
Status: `complex-core-first-coefficient-axis-jets-and-conditional-whole-inner-collar-join`.

The preceding pressure-moment checkpoint leaves the actual second coefficient
open. Its local source needs first-coefficient jets on the ENTIRE inner
interval, not only at positive radius. This checkpoint supplies an explicit
axis-safe bound for the regular F1,U1 profiles through three X derivatives
and two angular derivatives, then joins it to the existing activated-collar
bound. The new ingredient is a quantitative complex-X extension of the
first-coefficient system on the analytic SEED core. It does not assume
radial analyticity across the flat activation or coefficient cutoff.

This is an analytical consequence of the retained infinite-core norm,
six-variable system and sparse Volterra/Picard hypotheses. They are not
independently established unconditional here. All gauges, roots, widths,
cutoffs, frequencies, schedules and supports remain unchanged. Historical
frozen evidence is not edited. No actual second-coefficient H2,3, complete
physical pulse, covariance/mean, residual or blow-up conclusion follows.

## 1. An explicit complex-core operator envelope

Keep S>=2, 0<rho<=1, Lambda>0 and delta=rho/(800S), Xa=4/Lambda.
The infinite norm majorant from `core_all_order_foundation` bounds the
analytic reference fields on |Y|<=5, dist(eta,[-1,1])<=rho/8 by

`|partial_Y^r H| <= ||H|| r! 20^-r (8/5)^(r+1)`.

Here Y=Lambda X, `||F0||<=256S` and `||U0||<=9+S/Lambda`.
Since 4delta<rho/8, the following common envelope holds for both reference
fields and their X derivatives through order two on |X|<=5/Lambda and
the outer 4delta angular strip:

\[
 B=2(256S+9+S/\Lambda)(1+\Lambda)^2.                    \tag{1}
\]

Angular Cauchy disks of radius delta give the mixed raw bound
`H=B(1+delta^-1+2delta^-2)` for X order<=2, eta order<=2 on 2delta.
The regular average Ubar=int_0^1 U0(sX,eta)ds is holomorphic on the SAME
complex disk, since the real segment sX stays in it. Its radial derivatives
gain factors s^r<=1. Consequently the regular

`v0=(2eta U0-2Deta Ubar-(1-eta^2)Ubar_eta)/L`

and all derivatives needed by the full first source have the same bounds
used in `collar_operator_envelope`: |v0|,|v0_X|,|v0_XX|,|v0_eta|<=12H.
Throughout 2delta, |eta|<1.01 and |L|>.97 for 0<h<=.01.

Set x=5/Lambda and b=1+x, so |xi|<=sqrt(x)<=b when X=xi^2.
The literal first-system formulas contain only F0,U0, their derivatives,
regular averages and L reciprocals; they do NOT require F0^-1 or Phi^-1.
Thus a larger complex disk does not need an unjustified zero-free Phi
extension. The earlier absolute triangle estimates apply equally to complex
X and xi. With viscosity ONE, define

\[
 W_0=512(1+x)H^2,\quad Z_0=128(1+x)^2H,               \tag{2}
\]
\[
 C=\max\{1,128(1+x+b)H,32(1+x)H,(b+3x)W_0+2Z_0\}.   \tag{3}
\]

These bound the full w0, both shifted axial operators, and respectively
the row-sum norms of A0,A1 and f in

\[
 W_\xi+D_0W/\xi=A_0W+A_1W_\eta+f,\quad
 D_0=\operatorname{diag}(0,0,2,0,3,1),\quad
 W=(F_1,U_1,K_1,\Pi_1,P,Q).                            \tag{4}
\]

The pressure source is STILL `f3=-xi w0`; `f4=-2CF` and
`f5=2eta X w0/L-2CU` retain both axial viscosities. A0's Pi1 coupling
is retained. Zero pressure is imposed only at its original axis datum,
not reset at a later radius. The selected leading core still carries its
full inherited pressure datum through F0,U0.

## 2. Complex rays, regular axis and reflection parity

The regular inverse on a complex ray is

\[
 (Gg)_i(\xi)=\xi\int_0^1t^{d_i}g_i(t\xi)\,dt.         \tag{5}
\]

It is holomorphic at xi=0, equals zero there, and has norm at most
`|xi| sup|g|/(d_i+1)`. There is no estimate dividing a nonzero value by
xi. Since every d_i>=0, the diagonal kernel weights have modulus<=1.
The same ordered radial integrals, factorial volumes and angular Cauchy
losses as in the retained REAL-ray Picard proof therefore work on complex
rays. All angular derivatives of A1 have only rows 5,6 and columns 1..4.
For every diagonal kernel K, `A1(xi) K partial_eta^j A1(s)=0`.
The sparse word count and convergence are unchanged, not inferred from
pointwise nilpotence alone.

With strip loss delta/2, put z=2Cb and

\[
 \beta_*={6z^2\over\delta},\qquad
 M_{core}=(z+\beta_*)\exp(\beta_*).                    \tag{6}
\]

Because e<3, this majorizes the inherited `(z+beta)exp(beta)`,
beta=e z^2/(delta/2). The Picard series converges normally on smaller
complex xi/eta domains and bounds W by Mcore on |xi|<sqrt(5/Lambda),
the 3delta/2 angular strip. The exponential is kept symbolic; it is NOT
evaluated at the selected enormous parameters. Uniqueness follows by
iterating the same vanishing Picard remainder on the homogeneous difference.

Let P0=diag(1,1,1,1,-1,-1). Directly from the exact matrix entries,

`A_j(-xi)=-P0 A_j(xi)P0`, `f(-xi)=-P0 f(xi)`, and D0 commutes with P0.

Uniqueness with the same zero axis data gives `W(-xi)=P0 W(xi)`.
Thus F1,U1,K1,Pi1 are EVEN holomorphic functions of xi and descend to
holomorphic functions of X=xi^2 on |X|<5/Lambda. P,Q are ODD and are
not themselves declared holomorphic functions of X. Parity here already
appears qualitatively in the original construction; the new step is its
quantitative complex-core use. Holomorphy in xi WITHOUT evenness would
not suffice: xi becomes sqrt(X), with a singular X derivative at the axis.

On the real 0<=X<=Xa, the analytic seed equals the actual leading core.
The two first-coefficient solutions satisfy the SAME zero-axis Volterra
problem there, so uniqueness identifies them. The analytic extension for
X>Xa is used ONLY for Cauchy estimates of this unchanged core restriction.
It is not identified with the actually activated leading or first profile
past Xa. No analyticity across a flat join is claimed.

## 3. Axis-safe jets and the whole inner collar

At each real 0<=X<=Xa use an X disk of radius 1/(2Lambda).
It satisfies |Lambda X'|<=4.5<5. For G=F1 or U1, Cauchy proves

\[
 |\partial_X^rG|\le r!(2\Lambda)^r M_{core},\qquad
 |\partial_X^r\partial_\eta^kG|
 \le r!k!(2\Lambda)^r(2/\delta)^k M_{core},            \tag{7}
\]

for r<=3,k<=2 on the real angular interval. The angular radius delta/2
is a conservative common choice. In factorial angular C2 norm the cost is
`c2=1+2/delta+4/delta^2`.

Independently let Ccol,Cxi bound the full ACTUALLY activated matrices/source
and their xi derivative on 2delta from `collar_third_coefficient`, and let
Mcollar bound its solution on 3delta/2. It is a DIFFERENT bound from Mcore.
Set `a=min(1,Xa)>0`, so a^2<=Xa. The earlier positive-radius formulas give
scalar jet bounds `p_r Mcollar+q_r` on a^2<=X<=Xa exp(tc/4), r<=3,
on the delta angular strip. Multiplying by c2 is valid by angular Cauchy.
Together with (7), the nonnegative SUM

\[
 \boxed{\sup_{0\le X\le X_a e^{t_c/4}}
 |\partial_X^r(F_1,U_1)|_{C_\eta^2}
 \le c_2\{r!(2\Lambda)^rM_{core}+p_rM_{collar}+q_r\},\quad r\le3} \tag{8}
\]

is an explicit common component bound. The intervals overlap; inverse
radius terms are used only on the positive side. `joined_jets` stores the
three rational coefficients of Mcore,Mcollar,1, not their exponentials.
Taking the componentwise maximum over these nonnegative triples supplies
one common FACTORIAL C_eta^2 norm bound for r<=3 of the PRE-CUTOFF inner
first profile. For a pointwise source envelope needing RAW mixed derivatives,
multiply that common triple by TWO: the second eta derivative carries 2!.
`common_raw_mixed_bound` makes this conversion explicit; the factorial-norm
constant must not be reused unchanged as a raw second-derivative bound.

The SAME constants also hold on the smaller complex eta strip
`dist(eta,[-1,1])<delta/2`, if the factorial angular seminorm is taken there.
Indeed its radius-delta/2 Cauchy disks remain inside the positive-side
delta strip and strictly inside the core's 3delta/2 strip. The bounds are
therefore usable as holomorphic lower-order source inputs on that smaller
strip, not merely as real-axis derivative data. This does not enlarge any
radial domain or assert analyticity of a flat radial cutoff.

The actual finalized first coefficient still multiplies the original
eta-independent cutoff kappa. Hence its collar jets obey

`|partial_X^r(kappa G)|Ceta2 <= sum_(j=0)^r binom(r,j) k_j J_(r-j)`,

with k_j the ACTUAL cutoff derivative bounds. `cut_product` implements
this exact positive propagation. In particular k3 J0 remains; no cutoff
shape or derivative bound is invented. The disjoint exact I3 repair
continues to use complete future pressure/moment memory. Its global
angular bounds are NOT furnished by this local join.

E1=sqrt(2X)F1 is not treated as an axis-regular scalar X profile, nor is
E1^-1 used there. These regular F1,U1 estimates feed source equations;
they do not establish whole-axis E1 X derivatives, full physical Cartesian,
axial/time or potential-curl pulse bounds.

## 4. Checks and the next obligation

Symbolic checks cover every structural matrix/source reflection entry,
variable diagonal kernels and derivatives of the right A1. An exact
12-order manufactured analytic six-system has 72 scalar recurrence checks;
its first four components have only even powers, while a parity-breaking
source creates odd powers. Omitting its pressure source changes the axial
coefficient. These are manufactured equation controls, not actual PDE
coefficient solutions or a proof by finite Taylor truncation. Original
real-system samples independently check the entire structural pattern;
those floating comparisons are diagnostics, not interval certification.

Next supply the unchanged cutoff jets and global repair/angular inputs,
then build and differentiate the second inner system with its complete
lower pressure/viscosity source. This checkpoint removes the missing
axis-side local FIRST-coefficient jet estimate relative to retained
foundations, not the second-coefficient H2,3 itself. Complete source moments,
reserved-mean Pi1/Pi2, all physical pulse supports and derivatives,
full-wave covariance with axial velocity and pressure mean corrections,
and the complete viscous/cutoff/interaction residual remain obligations.

```
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

The fresh output must not already exist. Historical evidence is byte-bound
and unchanged. Passing tests do not make the construction unconditional.
