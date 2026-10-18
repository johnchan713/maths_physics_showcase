# Quantitative first-coefficient input for the compact background gate

Parent: `b29ca142b37c4200d5a90b824cb2824f70cc0688`.
This derives explicit finite upper bounds for the four first-coefficient
seminorms left as inputs by the [regional weighted gate](../regional_weighted_gate/README.md).
It propagates the actual collar Picard majorant through the unchanged radial
cutoff and exact five-moment linear repair, then through the full first
tangential residual, including pressure and shifted leading axial viscosity.
The resulting constant is a formula in primitive leading/collar data, not a
numerical evaluation of the selected enormous-parameter profile.

Everything below is an analytical implication of the retained exact leading
profile, regular analytic-eta collar, positive compact leading stress and
linear inverse lemmas. These hypotheses are not replaced by diagnostic
fixtures. Historical evidence, width C^-20, supports, roots and frequency are
unchanged. In particular the old core-only operator bound C^16 is NOT a bound
on the activated collar. No new heat-sign or canonical-seed result is claimed.

## 1. Collar values and radial jets without radial analyticity

Write Xa=a^2, Xc=c^2, 0<Xa<Xinner<Xkeep<Xcut<Xc, with the original radii.
The analytic collar ends at a_col^2=Xa exp(tc/4), after Xcut. On it the
actual six-variable system has maximum row-sum coefficient/source envelope
C_col>=1 on an eta strip of width sigma. Choose 0<Delta<=1 and 0<d<=1
with Delta+d<sigma. The retained sparse Picard argument gives

\[
 z=2C_{col}a_{col},\quad \beta=\mathrm e z^2/\Delta,\qquad
 M=(z+\beta)e^\beta.                                      \tag{1}
\]

Use |f|_1=sup|f|+sup|f_eta|, a product algebra norm on [-1,1].
Cauchy on the remaining strip bounds W and its eta derivative by

\[
 J=M(1+d^{-1}),\quad J_\eta=M(d^{-1}+2d^{-2}),\quad
 H=C_{col}(1+d^{-1})(J+J_\eta+1).                        \tag{2}
\]

Here W=(F1,U1,K1,Pi1,P,Q), P=F1_xi, Q=U1_xi, xi=sqrt(X).
The P,Q rows have diagonal coefficients 3,1, respectively. If HP,HQ denote
the corresponding rows of A0 W+A1 W_eta+f1, then exactly

\[
 F_{1,XX}=H_P/(4X)-P/X^{3/2},\qquad
 U_{1,XX}=H_Q/(4X)-Q/(2X^{3/2}).                         \tag{3}
\]

Consequently on [Xa,Xcut] both components have C_eta^1 bounds
J, Jx=J/(2a), Jxx=H/(4a^2)+J/a^3 for radial orders 0,1,2.
This uses the differential equations and angular Cauchy bounds, NOT radial
Cauchy estimates across the nonanalytic flat activation. Values are bounded
by J even on [0,Xa]. Constants are finite but may be much larger than every
fixed polynomial in C; (1) must not be silently replaced by such a polynomial.

## 2. Global moment costs retain the later pressure source

Tentatively cut only F1,U1 by the original eta-independent cutoff kappa,
0<=kappa<=1, with |kappa_X|<=k1, |kappa_XX|<=k2. Let f,u bound |F0|_1,|U0|_1
on [0,Xc]. Put w0=Omega0/X and retain its FULL support through Xv:

\[
 W_0=\int_0^{X_v}|w_0(x,\cdot)|_1dx,\qquad
 W_1=\int_0^{X_v}x|w_0(x,\cdot)|_1dx.                   \tag{4}
\]

The five ordinary incoming discrepancies obey, in their usual order,

\[
 (\mu_1,\mu_2,\mu_3,\mu_4,\mu_5)\le
 (X_{cut}J,\ X_{cut}^2J,\ 2fX_{cut}J+W_0/2,
 X_{cut}^2(u+f)J,\ 2uX_{cut}J+fX_{cut}^2J+W_1/2).       \tag{5}
\]

Indeed Pi1_X=2F0F1-w0/2, m2=int 2X F1,
m4=int 2X(U0F1+U1F0), and m5=int(2U0U1-2XF0F1+Xw0/2).
Only the tentative velocity terms stop at Xcut; the w0 terms do not.
The zero axis pressure gauge is preserved. On I3 and I4, U0=0 does NOT
imply M0=V0=w0=0. The historical core-only zero-memory wording was already
corrected by the later annular/support checkpoints and is not reinstated.

Let rho=rho3, e(eta)=e*/(1+eta^2), K=(sqrt(2)rho^(3/2))^-1,
0<lambda_patch<=.01. The actual linear solve cvec=-B^-1 z and inverse <800
give the C_eta^1 amplitude bound

\[
 C_p=800\max\{\mu_1/\rho,
 (\mu_1/\rho+4K\mu_4/e_*)/\lambda_{patch},
 K\mu_2,\ 4\mu_5/(\rho e_*),\ 4\mu_3/e_*\}.            \tag{6}
\]

The factor 4/e* is the full first reciprocal-amplitude jet norm, not 1/e.
The calculator conservatively replaces K by rho^(-3/2), with rho supplied
as the square of an exact positive radius. No small-target or quadratic
repair assumption occurs. Equations (4)--(6) expose the future-pressure
dependence of an earlier compact coefficient; they do not reset the pressure
tail at the end of the velocity bumps.

## 3. Cutoff and repair bounds for the actual velocity coefficients

Let p be a positive lower X bound for all five disjoint I3 bump supports,
0<rp<=sqrt(2p), and Bj=max sup|D_y^j beta_i| for j=0,1,2, y=log(X/rho).
These are the unchanged smooth unit-mass bumps, not new localizers.
The repaired F1,U1 have common C_eta^1 upper bounds S,Sx,Sxx on the required
compact region (S also holds on [0,Xc]):

\[
\begin{split}
S&=J+C_pB_0\max(1,r_p^{-1}),\\
S_x&=J_x+k_1J+C_p\max\{B_1/p,(B_1+B_0/2)/(pr_p)\},\\
S_{xx}&=J_{xx}+2k_1J_x+k_2J+
C_p\max\{(B_2+B_1)/p^2,(B_2+2B_1+3B_0/4)/(p^2r_p)\}.
\end{split}                                                   \tag{7}
\]

Disjoint bumps avoid a factor five. The E1 patches become F1 patches by
division by sqrt(2X); both derivatives of that factor are charged. The
logarithmic conversion D_X^2 beta=X^-2(D_y^2-D_y)beta is retained.

## 4. A compact stress bound with pressure and both axial viscosities

Set A=1/2+h, D=1/2-h, lambda=2h, b=-1-h, k=-A, Lmin=1-2h, X=Xc.
Let fx,ux bound the leading radial derivatives on [Xa,Xc], v bound |V0|,
w bound |w0|_1 on [0,Xc], and CF,CU bound the COMPLETE shifted leading axial
second derivatives Z_(b-D)Z_b F0, Z_(k-D)Z_k U0 on [Xa,Xc]. Define

\[
\begin{gathered}
 T(t)=((|t|+D)S+XS_x)/L_{min},\quad
 Z(t)=((2|t|+1)S+2XS_x)/L_{min},\\
 Z_F=((2|b|+1)f+2Xf_x)/L_{min},\quad
 Z_U=((2|k|+1)u+2Xu_x)/L_{min},\\
 V=X(3+2A)S/L_{min},\quad P_X=2fS+w/2,\quad P=XP_X,\quad Z_P=5XP_X/L_{min}.
\end{gathered}                                                 \tag{8}
\]

V bounds the reconstructed V1, using M1=int U1 and its eta derivative.
P bounds |Pi1|_1 from its FORWARD integral, with Pi1(0)=0. Since the first
pressure exponent -2A+lambda=-1, ZP bounds its full axial derivative;
neither its value nor eta derivative is dropped. Direct triangle estimates
on the exact order-one residual coefficients give

\[
\begin{split}
H_\theta={}&T(b+\lambda)+2\nu(XS_{xx}+2S_x)
+v(S_x+S/X_a)+V(f_x+f/X_a)\\
&+uZ(b+\lambda)+SZ_F+\nu C_F,\\
H_z={}&T(k+\lambda)+2\nu(XS_{xx}+S_x)+Z_P+vS_x+Vu_x
+uZ(k+\lambda)+SZ_U+\nu C_U.
\end{split}                                                   \tag{9}
\]

Here nu=1 for the construction. Both positive-order transverse diffusion
and shifted order-ZERO axial diffusion are included. The code retains
general positive nu for omission controls. These are coefficients of the
augmented background, not the complete wave-covariance residual.

The exact inner coefficient equations give T1=0 through Xinner>Xa.
Thus no minimum of |T0| at its flat zero Xa is taken. On the CLOSED interval
[Xinner,Xc], let tau>0 be a leading stress lower bound, uniform in eta.
With rtheta1=R*a_theta, rz1=a_z and R=sqrt(2X), the forward definitions are

\[
 T_{1\theta}=-X^{-1}\int_{X_{inner}}^X x a_\theta dx,
 \qquad T_{1z}=-(2X)^{-1/2}\int_{X_{inner}}^X a_z dx.      \tag{10}
\]

In particular |W1|=|T1|/|T0|<=(Xc Htheta/2+c Hz)/tau, using the sum bound
for the Euclidean norm and sqrt(Xc/2)<=c. Below Xinner the quotient is zero.
Compactness supplies a positive tau from the retained strict leading
positivity lemma, but its quantitative value has not been evaluated here.

## 5. Compose the previous gate and state the precise remaining inequality

Valid conservative first seminorm inputs are

\[
 A_E=2cS,\quad A_{EX}=2cS_x+S/a,\quad A_{UX}=S_x,\quad
 A_W=(X_cH_\theta/2+cH_z)/\tau.                         \tag{11}
\]

Here a<=sqrt(2Xa) and 2c>=sqrt(2Xc) are safe radial enclosures. Inserting
(11) into the preceding C_K formula, with Rlo=a, Rhi=2c, defines an explicit
finite C*. Its shared higher-order budget still contributes q^(2h)/2.
The whole 4Q chart therefore meets this compact value/cone interface if

\[
 \boxed{2h(\ell-2)\ \ge\ \log_2(C_*/r_{K_c}).}           \tag{12}
\]

The established exact outer sign/cone theorem completes this PARTICULAR
regional cover. All other value-domain, stress-magnitude, root, operator,
freezing and residual gates remain separate. Equation (12) is an explicit
sufficient inequality, not an actual failed-cone test and not an acceptance
certificate for the historical first band. Its primitive C_col, strip widths,
leading jet/source envelopes, bump/cutoff jets, tau and leading r_Kc inputs
still need justified quantitative evaluation on the selected construction.
Finiteness alone gives some later threshold; it does not prove a prescribed
old threshold, a single old field/force, uniform pulse bounds or blow-up.

## Checks and next obligation

`seminorms.py` has exact symbolic identities for (3), E-to-F derivatives and
both physical radial integration measures; exact rational pressure-memory
and triangle controls include eta=+-1. A separately assembled legacy full
first residual agrees at 54 manufactured interior points. These are algebra
checks/numerical diagnostics, not actual-profile suprema or a formal proof.
The composed rational example uses the conservative M<=2^800 for illustrative
C_col=2, a_col=2, Delta=1/2, d=1/4, sigma=1. Indeed e<3 gives beta<384 and
(8+beta)exp(beta)<2^777. These are consistent abstract envelope/radius inputs,
not the actual narrow-collar parameters or an evaluation of its solution.
The output binds the original source bytes and preserves all historical
evidence; fresh replays use their original comparators.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
Next replace the primitive envelopes and compact leading gap/floor inputs by
justified actual bounds and verify their combined accepted-tail inequalities.
Complete all-scale pulse/covariance/residual obligations remain open.
