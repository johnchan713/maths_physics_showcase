# Complete second inner source and a conditional positive-collar J3

Parent: `6bc32141ce0df0ebe51f6d1ea55c599c232166a5`.
Status: `full-second-inner-source-derivative-envelope-and-conditional-positive-collar-J3`.

The previous checkpoint supplies first-coefficient axis/collar jets. This
checkpoint builds the COMPLETE second-order local source and its xi derivative,
with all current/lower exponents distinguished, and obtains an explicit
conditional positive-collar third-jet estimate for the second inner solution.
It also identifies and bounds the EXTRA third-angular leading derivative in
the shifted radial axial viscosity. Reusing only the previous leading
eta<=2 raw envelope there would omit a needed term.

No previous checkpoint claimed this complete second-source estimate; this
is an advance, not an allegation that historical general-order identities
dropped viscosity. The construction, roots, fixed width C^-20, frequency,
cutoffs, gauges, schedules, supports and historical frozen evidence stay
unchanged. Actual cutoff norms, global repair, finalized H2,3, full physical
pulses, covariance/means and complete residual improvement remain open.
Analytical implications below retain the core/Picard/regional hypotheses;
manufactured diagnostics are not actual selected PDE solutions or proof.

## 1. Match the derivative input to the actual second source

Write A=1/2+h, D=1/2-h, L=1-2h eta^2, d=1-eta^2,
lambda_j=2hj, and Z_a=(2a eta-2eta X partial_X+d partial_eta)/L.
Use the original viscosity nu=1. The regular lower radial velocity is

`v_i=[2eta U_i-2(D+lambda_i)eta Ubar_i-d(Ubar_i)_eta]/L`,
`Ubar_i=int_0^1 U_i(sX,eta)ds`.

The complete w1=Omega1/X is

\[
\begin{split}
w_1={ (1-\lambda_1)v_1+X(v_1)_X+D\eta(v_1)_\eta\over L}
-2(2(v_1)_X+X(v_1)_{XX})\\
+\sum_{i+j=1}\left[v_i(v_j/2+X(v_j)_X)
+{U_i[2(\lambda_j-1)\eta v_j-2\eta X(v_j)_X+d(v_j)_\eta]\over L}\right]
-Z_{-1-D}Z_{-1}v_0.                                  \tag{1}
\end{split}
\]

It contains BOTH ordered pairs (0,1),(1,0). The last term is present now,
not postponed to pressure order three. Because v0 already contains
`-d Ubar0_eta/L`, its second eta derivative contains
`-d Ubar0_etaetaeta/L`. Differentiating (1) in xi also retains its X
derivative. A sufficient raw input on the whole inner radial interval is

\[
 |\partial_X^r\partial_\eta^k(F_0,U_0)|\le H_0
       \ (r\le3,k\le3),\qquad
 |\partial_X^r\partial_\eta^k(F_1,U_1)|\le H_1
       \ (r\le3,k\le2).                              \tag{2}
\]

The first profile here is the ACTUAL cut coefficient on that interval,
not just its uncut inner precursor. We may widen H0,H1 to be >=1.

There is no new free higher leading norm: the earlier leading radial
factor K3 bounds radial orders<=3 on the outer 4delta strip, with
E=exp(25 ell S). Angular Cauchy on 2delta gives

`H0=K3 E (1+delta^-1+2delta^-2+6delta^-3)`.

This raises ONLY the angular part of the earlier bound and preserves the
unchanged activation. For H1, apply the previous axis/collar factorial
C_eta^2 bounds to the SAME eta-independent cutoff by the full binomial
rule through radial order three. Multiply by TWO for raw eta-second
derivatives, as the previous checkpoint requires. `lower_cut_interface`
records the resulting coefficients of Mcore,Mcollar,1. The unchanged
cutoff's actual k1,k2,k3 remain inputs; no replacement cutoff is installed.
The bound holds on the smaller complex strip dist(eta,[-1,1])<delta/2.

## 2. Explicit regular-v costs, including the extra angular derivative

On the retained strip |eta|<=e=1.01, |d|<=2.021 and |L|>=l=.97.
Set c=D+lambda_i (c<=.51 for i=0,1), and let N be the numerator of vi.
For every needed X derivative the same radial-average estimate gives

`|N|/H <= 2e+2ce+2.021`,
`|N_eta|/H <= 2+2e+2c(1+e)+2e+2.021`,
`|N_etaeta|/H <= 4+2e+2c(2+e)+2+4e+2.021`.

The corresponding reciprocal bounds are

`R0=1/l`, `R1=4he/l^2`, `R2=4h/l^2+32h^2 e^2/l^3`.

Thus the raw costs are `R0 N0`, `R0 N1+R1 N0`, and
`R0 N2+2R1 N1+R2 N0`. Increasing h to .01 and c to .51 bounds EVERY
positive term over the entire retained interval. This proves the uniform
safe choices

`|partial_X^r partial_eta^k v0| <=20H0` for r<=3,k<=2,
`|partial_X^r partial_eta^k v1| <=12H1` for r<=3,k<=1.

These are regular even at X=0; no average is replaced by U at one point.
The upper endpoint majorant, not two sampled h values, proves the interval
claim. In the auxiliary profile U0=eta^3, F1=U1=0, at eta=0 the lower
pointwise U0 jets through eta order two vanish, but w1=6 from the retained
last term in (1). This is a derivative-demand control, not an actual PDE
profile family or a denial of the existing analytic Cauchy bounds.

## 3. Keep current and lower exponents separate

The second pressure and momentum sources in the exact hierarchy are

\[
 P_2=F_1^2-w_1/2,\quad
 S_{\theta2}=v_1(XF_{1X}+F_1)+U_1 Z_{-1+h}F_1
                                      -Z_{-1+h-D}Z_{-1+h}F_1,
\]
\[
 S_{z2}=Xv_1U_{1X}+U_1Z_{-1/2+h}U_1
                                      -Z_{-1/2+h-D}Z_{-1/2+h}U_1,
\]
\[
 f_2=(0,0,0,2\xi P_2,2S_{\theta2},
                       2[S_{z2}-2\eta XP_2/L]).       \tag{3}
\]

The current matrix is the same six-variable formula with the CURRENT
lambda2, whereas its leading coefficient fields still have order zero.

| Quantity | Exponent retained |
|---|---|
| Current matrix order | lambda2=4h |
| Lower first source and regular v1 | lambda1=2h |
| Leading radial axial viscosity | a=-1, a-D=-1-D |
| Lower theta axial operator | -1+h |
| Lower axial operator | -1/2+h |
| Current Pi2 matrix coupling | -1+2h |

Replacing lambda1 by lambda2 throughout the source is incorrect. So is
retaining the first-order matrix pressure coefficient -1 at order two.
The exact finite diagnostics detect both omitted previous viscosity and
blanket current-exponent substitution.

`operator_and_source` bounds EVERY matrix row and source entry with a pair
`(value,xi gradient)`. Products use `(v,d)(w,e)=(vw,dw+ve)` with no signed
cancellation. It uses X=(x,2b), xi=(b,1), leading jets=(H0,2bH0), first
jets=(H1,2bH1), v0=20 times the leading pair, v1=12 times the first pair.
Here X<=x, xi<=b and b^2>=x. The full shifted Z_a-D Z_a expansions retain
all coefficient derivatives, including the ones multiplying v0_etaeta.
Both pressure couplings to the current Pi2 equation remain. This yields

`|A0_2|,|A1|,|f2| <= Q2`,
`|partial_xi A0_2|,|partial_xi A1|,|partial_xi f2| <= Q2xi`.

The method is an explicit positive rational expression, not an unspecified
source supremum. Its value and gradient expressions have degree at most
two in H0,H1. Hence, for a SCALAR envelope T>=max(1,H0,H1), evaluate the
same pair map at H0=H1=1 to get q,qxi and conclude

\[
                         Q_2\le qT^2,\quad Q_{2\xi}\le q_\xi T^2. \tag{4}
\]

This keeps enormous actual constants symbolic. T is a fixed supremum bound,
not an eta-dependent field whose derivatives have been discarded.

## 4. A conditional positive-collar third jet for the second solution

Keep the actual delta=rho/(800S)<=1/1600 and set delta2=delta/4,
so the needed outer strip 2delta2=delta/2 fits the
first-coefficient input. The same sparse derivative-block/diagonal-kernel
Picard proof applies at order two: D0 and A1's structure are unchanged.
With the original zero axis data and a valid full source bound, put

`z2=2b q T^2`, `B2=6 z2^2/delta2`,
`M2=(z2+B2) exp(B2)`.

This bounds the second six-vector on 3delta2/2. In particular B2 grows as
T^4. M2 MUST NOT be evaluated at T=1; scaling only the subsequent radial
formula would lose its enormous source dependence.

Choose 0<a<=sqrt(Xa), e.g. a=min(1,Xa). The existing ordinary-xi
differentiation formula supplies the positive-radius bounds at
xi>=a. Let `(A_r,B_r)` be that formula's C_eta^1 affine coefficients
computed at q,qxi and delta2. Its coefficient degree in Q2,Q2xi is zero
for r=0,1, at most one for r=2, and at most two for r=3, with all terms
nonnegative. Therefore, with p=(0,0,2,4),

\[
 \boxed{|\partial_X^r(F_2,U_2)|_{C_\eta^1}
          \le T^{p_r}(A_rM_2+B_r),\qquad r\le3,\ \xi\ge a.} \tag{5}
\]

`second_collar` records the rational affine coefficients and the separated
z2/B2 factors, never the selected exponential. This is for the INNER
second solution, conditional on the actual cutoff-dependent T. It is NOT
an axis third-jet bound for the second coefficient, the globally repaired
H2,3, or a new accepted diagonal/physical pulse schedule.

## 5. Pressure memory, checks and remaining obligations

Exact I3 velocity repair is disjoint from the inner interval; its forward
velocity averages cannot feed backward into the inner source. This fact
does not erase later pressure memory or assign Pi1/Pi2 zero on I4.
The global second repair still needs the complete signed pressure-source
moments through Xv and its lower F1/U1 products, as in the preceding
`higher_pressure_source_moments` checkpoint. Neither (3) nor (5) bounds
those global discrepancies by truncating them at Xcut. The fixed axis
datum is never reset at an intermediate radius.

Independent symbolic full-source differentiation and nine manufactured
original-hierarchy comparisons check all matrix/source values and xi
gradient envelopes; finite-difference checks are diagnostics, not interval
certificates. Symbolic controls retain the regular axial exponent shift,
the third-angular leading term and the positive polynomial degree behind
(5). Manufactured fixtures are not the selected enormous solution.

Next quantify the unchanged cutoff jets and global repair/angular inputs,
complete the second coefficient's finalized third-jet bound H2,3 (including
its cutoff and exact repair), and check retained-scale acceptance. Second
axis/higher coefficient jets, every physical pulse derivative/support,
full-wave covariance with axial velocity and pressure means, full viscous,
cutoff and interaction residual, and independent inherited foundations
remain open. No unconditional pulse, residual improvement or blow-up claim
is made.

```
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

The fresh output must not exist. Historical evidence is unchanged.
