# Quantitative third inner-coefficient jet on the actual cutoff collar

Parent: `2fcf54abdc2b3faf661d73e0feae79197f96c565`.
Status: `actual-activated-third-leading-jet-and-conditional-inner-coefficient-J3`.

The previous checkpoint exposes two separate higher-jet inputs: the
unchanged cutoff's shape derivatives and the inner coefficient's radial
jets. This checkpoint supplies an explicit third-order bound for the
SECOND input on the actual cutoff collar. It first bounds the actual
activated leading third jets, then differentiates the full six-variable
system with pressure and shifted viscosities retained. It uses angular
Cauchy estimates, not nonexistent radial analyticity at the flat join.

This is an analytical implication of the retained exact-core, sparse
Picard and regular-collar hypotheses. It does not independently certify
those foundations, evaluate the enormous selected solution, identify the
cutoff shape or accept the original first band. No field, scale, width,
cutoff, root, tilt, pressure gauge or historical frozen evidence is changed.

## 1. Bound the actual leading third radial jet

Keep `S>=2`, `0<rho<=1`, `Xa=4/Lambda`, `t=t1=C^-20`,
`ell=tc/4<t/20`, `tc=t/(10 sqrt(log C))`, and `delta=rho/(800S)`.
As in the retained collar estimate put `x=4.1/Lambda`, `b=1+x`.
The exact-core majorant on the outer 4delta eta strip supplies, for r<=3,

\[
 |\partial_X^r(F_r,U_r)|\le
 B_3:=2(256S+9+S/\Lambda)(1+\Lambda)^3.                  \tag{1}
\]

Indeed `r!20^-r(8/5)^(r+1)<2` for r=0,1,2,3. The same majorant bounds
the first three Y derivatives of Phi by S, and the unchanged zero-free
neighborhood gives `|1/Phi|<5`. For `a_r=-2Y Phi_Y/Phi`, Y<=4.1<5,
direct quotient differentiation therefore gives

`|a_r|<=50S`, `|a_r,y|<=2000S^2`, `|a_r,yy|<=100000S^3`.

For the last bound, `|(1/Phi)_Y|<=25S` and
`|(1/Phi)_YY|<=250S^2+25S`. Consequently
`|a_r,Y|<=280S^2`, `|a_r,YY|<=3000S^3`, and
`|a_r,yy|<=5*280S^2+25*3000S^3<=75700S^3`.

Use the SAME squared-flat activation `e_a=(1-t)sigma(y/t)`. Its earlier
actual-shape bounds give `|e_a,y|<=64/t`, `|e_a,yy|<=1216/t^2`.
For `J=(1/2)int_0^y e_a a_r`, define

\[
 j_1=25S,\quad j_2=1600S/t+1000S^2,\quad
 j_3=30400S/t^2+128000S^2/t+50000S^3.                    \tag{2}
\]

These bound J_y,J_yy,J_yyy on the SAME complex strip. Put
`E=exp(25ell S)` and

`e=(1,j1,j2+j1^2,j3+3j1*j2+j1^3)`,
`r=B3*(1,x,x+x^2,x+3x^2+x^3)`.

The exact formula `F_a=F_r exp J` yields logarithmic derivative bounds
`f_m=E*sum_(j=0)^m binom(m,j)r_(m-j)e_j` through m=3. The actual axial
identity `U_a,y=(1-e_a)U_r,y`, NOT substitution of the reference, yields

\[
 (u_0,u_1,u_2,u_3)=B_3\bigl(1+\ell x,x,
 x+x^2+64x/t,
 x+3x^2+x^3+128(x+x^2)/t+1216x/t^2\bigr).              \tag{3}
\]

Convert each by `partial_X^m=X^-m product_(i<m)(partial_y-i)` on the
activation side, using Xa in denominators, and sum the resulting positive
bounds with B3. Since E>=1, the code returns a factor K3 such that ALL
radial orders <=3 of both actual leading fields are bounded by K3*E.
The core side uses (1) instead of division by X. Flatness joins the two
sides, including the axis. Radial averages obey the same bounds by their
fixed-integral formula, with integration factors <=1.

Angular Cauchy disks of radius delta give

\[
 H_3=K_3E(1+\delta^{-1}+2\delta^{-2})                  \tag{4}
\]

for X order <=3, eta order <=2 on 2delta. Exponentials are separated:
`E<=2^ceil(50ell S)` since log(2)>1/2. Actual giant powers are never
materialized by the factor calculator, nor replaced by a polynomial in C.

## 2. Differentiate the COMPLETE matrix and source envelope

The retained system is

\[
 W_\xi+D_0W/\xi=A_0W+A_1W_\eta+f,\quad
 D_0=\operatorname{diag}(0,0,2,0,3,1),\quad
 W=(F_1,U_1,K_1,\Pi_1,P,Q),\quad(P,Q)=(F_{1\xi},U_{1\xi}). \tag{5}
\]

Here Q as a COMPONENT is distinct from the scalar envelope below.
`operator_pairs` implements direct positive-product majorants of every
nonzero entry in `first_correction_construction/system.py`, viscosity ONE.
It propagates a pair `(value bound,xi derivative bound)` by
`(v,d)*(w,e)=(vw,dw+ve)` and adds without cancellation. Its primitive inputs
are `X=(x,2b)`, `xi=(b,1)`, any required leading jet `(H3,2bH3)`, and each
of `v,v_X,v_XX,v_eta` bounded by `(12H3,24bH3)` using the regular radial
average. Thus the derivative of v_XX retains v_XXX, not a second-jet guess.
On 2delta use `|eta|<1.01`, `|1-eta^2|<2.021`, `|L|>.97`.

Every regular pressure-source term is included:

\[
 w_0=(v+Xv_X+D\eta v_\eta)/L+v(v/2+Xv_X)
 +U[-2\eta(v+Xv_X)+(1-\eta^2)v_\eta]/L
 -2(2v_X+Xv_{XX}).                                     \tag{6}
\]

For BOTH full shifted operators `Z_(k-D)Z_k` the calculator expands the
product, including coefficient derivatives. Writing `Z_k=a partial_X+
b_eta partial_eta+c_k` (b_eta here is a coefficient name, not its derivative),
the coefficients of FX,Feta,F in the second operator are respectively

`a a_X+b_eta a_eta+c_shift a+a c_k`,
`b_eta (b_eta)_eta+c_shift b_eta+b_eta c_k`,
`b_eta (c_k)_eta+c_shift c_k`.

The XX,Xeta,etaeta terms are `a^2,2a b_eta,b_eta^2`. The retained bounds
`|a|<=3X`, `|a_X|<=3`, `|a_eta|<=3X`, `|b_eta|,|(b_eta)_eta|<=3`,
`|c_k|,|(c_k)_eta|<=6` and their exact radial X dependence supply the pairs.
No coefficient is frozen when differentiated. This includes the required
leading mixed X/eta jets of order (3,2).

Take the maximum of row sums and source components to obtain

\[
 \|A_0\|,\|A_1\|,|f|\le C,\qquad
 \|A_{0\xi}\|,\|A_{1\xi}\|,|f_\xi|\le C_\xi.          \tag{7}
\]

The pressure source entry is STILL `-xi*w0`, and the other nonzero source
entries are STILL `-2CF` and `2eta X*w0/L-2CU`. Both matrix couplings to
Pi1 remain. Row/source pairs in the record are explicit rational expressions
in supplied inputs, not unspecified suprema. Every expression has degree at
most two in H3. If its separated exponential factor is E, evaluate the
positive expressions at H3/E and multiply BOTH final bounds by E^2;
lower-degree terms are then conservatively enclosed, since E>=1.

## 3. Obtain J3 by angular Cauchy and ordinary radial differentiation

Let C>=1,Cxi>=0 bound (7) on 2delta, and let M bound the SAME zero-axis
solution W on the 3delta/2 strip. The retained sparse Picard estimate gives

\[
 z=2Cb,\quad \beta=\mathrm e z^2/(\delta/2),\qquad
 M=(z+\beta)e^\beta.                                    \tag{8}
\]

The derivative-block sparsity already proved for (5) is used only here,
with strip loss delta/2. It is not inferred from pointwise nilpotence alone.
On the delta strip Cauchy gives `M1=2M/delta`, `M2=8M/delta^2` for Weta
and Wetaeta. Coefficient eta derivatives cost at most C/delta from 2delta.
Choose ANY `0<a<=sqrt(Xa)`; e.g. `a=min(1,Xa)` avoids evaluating a square
root. On the cutoff annulus xi>=a put

\[
 B=C+3/a,\quad B_\xi=C_\xi+3/a^2,\quad
 V_1=BM+CM_1+C,                                        \tag{9}
\]
\[
 V_{1\eta}=(C/\delta)M+(B+C/\delta)M_1+CM_2+C/\delta,
\]
\[
 V_2=B_\xi M+BV_1+C_\xi M_1+CV_{1\eta}+C_\xi.          \tag{10}
\]

These follow by differentiating `Wxi=(A0-D0/xi)W+A1 Weta+f` once in eta
and once in xi. In particular they retain the derivative of D0/xi, the
matrix derivatives, Wetaeta and fxi. No commutativity is assumed.
Consequently `|Wxi|<=V1`, `|Wxixi|<=V2` uniformly there.

For either scalar F1,U1, its component P or Q in (5) gives EXACTLY

\[
 \partial_X^3 F_1={P_{\xi\xi}\over8\xi^3}
                -{3P_\xi\over8\xi^4}+{3P\over8\xi^5}.
 \tag{11}
\]

The same formula holds for U1 with Q. Therefore the following C_eta^1
inner bounds are explicit, with `ceta=1+2/delta`:

\[
 J_0=c_\eta M,\quad J_1=c_\eta M/(2a),\quad
 J_2=c_\eta[V_1/(4a^2)+M/(4a^3)],
\]
\[
 \boxed{J_3=c_\eta[V_2/(8a^3)+3V_1/(8a^4)+3M/(8a^5)].}  \tag{12}
\]

The last angular loss uses a delta/2 disk centered on the real eta interval.
All formulas are affine in the unevaluated M; the calculator stores those
two rational coefficients instead of a gigantic exponential. These Jj
are required only on the cutoff collar, NOT estimates of singular inverse
xi expressions at the axis. Separate regular-axis formulas are unchanged.

## 4. Join, evidence and remaining boundary

Insert (12) into the previous higher-jet interface. In particular the cut
product needs

`P3=J3+3 k1 J2+3 k2 J1+k3 J0`,
`k3=Xkeep^-3*(L3/Delta^3+3L2/Delta^2+2L1/Delta)`.

Add the actual third repair-bump budget and the SAME complete pressure/
moment inverse norm. Thus J3 is no longer a free local-coefficient input
relative to the retained core/Picard hypotheses, but L3 of the unchanged
cutoff REMAINS necessary. Cp still carries all future w0 integrals through
Xv, not just this collar. A local source envelope neither truncates that
memory nor resets incoming pressure. No replacement cutoff is installed.

Symbolic checks cover the full third product, physical logarithmic factors,
and (11). Eighteen manufactured original-system comparisons check every
matrix row and source value/xi-gradient envelope; gradients use centered
finite differences and are explicitly diagnostics, not interval proofs.
Missing Cxi, smaller angular loss and omitted inverse-radius terms are
negative controls. A rational toy composition uses M<=2^4000 only for its
stated C=2,b=2,loss=1/8; it is NOT the selected actual solution.

Next supply the unchanged cutoff's required shape jets and all remaining
same-field global inverse/jet inputs, then check retained-scale acceptance.
Higher coefficient orders, all Cartesian/axial/time pulse derivatives,
complete covariance/axial-pressure means and full viscous-cutoff interaction
residual remain open. No unconditional pulse or blow-up follows.

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output /fresh/path.json`.
