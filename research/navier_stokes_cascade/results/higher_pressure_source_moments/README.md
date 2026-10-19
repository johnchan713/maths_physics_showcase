# Higher pressure-source moments with shifted axial viscosity

Parent: `39dca194ea7dbe3f8e505ecaaa42729ca5e64754`.
Status: `ordered-full-support-pressure-source-moments-and-conditional-second-repair-bound`.

The preceding finite-head checkpoint isolates the actual second coefficient's
third-jet norm H2,3. This checkpoint advances its GLOBAL pressure/moment
input: it generalizes the earlier **leading** signed-source identity to every
finite order, integrates the previous-order axial viscosity without discarding
it, and supplies explicit conditional second-coefficient repair discrepancies.
It does not supply the local second-coefficient solution or evaluate H2,3.
The earlier general conservative stress totals already retained shifted axial
viscosity; this is a new quantitative radial-source moment calculation, not
a correction alleging that those historical identities dropped it.
There is no change to the leading roots, construction, cutoff, width,
frequency, axis pressure gauge, support, diagonal schedule or frozen evidence.

## 1. Retain the entire ordered source

Set A=1/2+h, D=1/2-h, d=1-eta^2, L=1-2h eta^2,
lambda_j=2hj, M_j=int_0^X U_j and v_j=V_j/X, with its regular axis value.
The exact coefficient divergence relation is

\[
 LV_j=2\eta XU_j-2(D+\lambda_j)\eta M_j-d(M_j)_\eta,
 \quad L(V_j)_X=2(A-\lambda_j)\eta U_j+2\eta X(U_j)_X-d(U_j)_\eta. \tag{1}
\]

For k=n-1, the original hierarchy has

\[
\begin{split}
w_k={\Omega_k\over X}={ (1-\lambda_k)v_k+X(v_k)_X+D\eta(v_k)_\eta\over L}
-2\nu(2(v_k)_X+X(v_k)_{XX})\\
+\sum_{i+j=k}\left\{v_i(v_j/2+X(v_j)_X)
+{U_i[2(\lambda_j-1)\eta v_j-2\eta X(v_j)_X+d(v_j)_\eta]\over L}\right\}
-\nu Z_{a-D}Z_a v_{k-1},\qquad a=\lambda_{k-1}-1,                \tag{2}
\end{split}
\]

where Z_a=(2a eta-2eta X partial_X+d partial_eta)/L and the last term
is absent for k=0. It starts in pressure order n=2, not n=3.
All ordered pairs are present, including both (0,1) and (1,0) when k=1.
Define B_k=sum U_i v_j and C_k=sum v_i v_j over i+j=k. Product rules
and (1) give

\[
\boxed{w_k=(\mathcal J_k)_X+
{-\lambda_k v_k+D\eta(v_k)_\eta+d(B_k)_\eta
                   +2(\lambda_k-A)\eta B_k\over L}
+C_k/2-\nu Z_{a-D}Z_a v_{k-1},}                         \tag{3}
\]
\[
 \mathcal J_k={V_k-2\eta\sum_{i+j=k}U_iV_j\over L}
                     +X C_k-2\nu(V_k)_X.              \tag{4}
\]

For each ordered pair, the divergence substitution contributes
`(V_i)_X v_j-2(A-lambda_i)eta U_i v_j/L`; together with
`v_i(V_j)_X-v_i v_j/2` this gives
`partial_X(Xv_i v_j)+v_i v_j/2`. The exponent in the remaining mixed
term is lambda_i+lambda_j=lambda_k. This establishes (3) without
symmetry guesses, omission of a lower-order orientation or division by d.

## 2. Integrate axial viscosity, not just radial viscosity

Use the SAME full leading axial support Xv, which contains every finalized
lower coefficient. Exact functional M_j(infinity,eta)=0 makes V_j and its
fixed jets zero after their true supports. We require sufficient regularity
for the displayed derivatives and a common endpoint independent of eta.
Define a_{kj}=int_0^Xv X^j v_k, b_{kj}=int X^j B_k,
c_{kj}=int X^j C_k, j=0,1, and the scalar angular operator

\[
 T_b f={d f_\eta+2b\eta f\over L}.                     \tag{5}
\]

Integration by parts of `-2eta X f_X` proves, for compact regular f,

\[
 \int_0^{X_v}X^j Z_a f=T_{a+j+1}\int_0^{X_v}X^j f,
 \quad
 \int_0^{X_v}X^j Z_{a-D}Z_a f
       =T_{a+j+1-D}T_{a+j+1}\int_0^{X_v}X^j f.         \tag{6}
\]

There is no radial derivative of f in the final full moment, but the TWO
angular operators remain. Let
`H_kj=T_(lambda_(k-1)+j-D) T_(lambda_(k-1)+j) a_(k-1,j)` for k>=1,
and H_0j=0. Integrating (3), then X times (3), yields

\[
\boxed{ I_{k0}=2\nu v_k(0)+
{-\lambda_k a_{k0}+D\eta(a_{k0})_\eta+d(b_{k0})_\eta
                  +2(\lambda_k-A)\eta b_{k0}\over L}
                  +c_{k0}/2-\nu H_{k0},}             \tag{7}
\]
\[
\boxed{ I_{k1}=
{-(1+\lambda_k)a_{k1}+D\eta(a_{k1})_\eta+d(b_{k1})_\eta
                  +2(\lambda_k+D)\eta b_{k1}\over L}
                  -c_{k1}/2-\nu H_{k1}.}             \tag{8}
\]

Here I_kj=int_0^Xv X^j w_k. The radial viscous weighted total is zero,
whereas H_k1 is generally NONZERO. Smooth axis regularity alone allows
v_k(0) nonzero; the actual positive-order zero U_k axis trace and its eta
jets additionally imply v_k(0)=0 for k>=1. Equation (7) safely retains
the term without assuming those stronger data. For k=0 the formulas
reduce exactly to the historical signed-pressure-source theorem.

For a premature finite endpoint R, equations (7),(8) require respectively
`J_k(R)` and `R J_k(R)+2nu V_k(R)`, AND the axial corrections. For each
j put f=v_(k-1), beta=lambda_(k-1)+j,

\[
 E_1=-2\eta R^{j+1}f(R)/L,\qquad
 E_2=T_{\beta-D}E_1-2\eta R^{j+1}(Z_a f)(R)/L.        \tag{9}
\]

The moment correction is the stated radial boundary minus nu E2.
Angular differentiation of E1 is included. Thus (6) cannot be applied
at Xcut, I3 or I4 by discarding the future support. In particular local
zero positive-order velocity does not erase Pi2's inherited I4 constant.
These are signed total identities, not pointwise or absolute source bounds.

## 3. Explicit C_eta^1 bounds and the new angular input

Use |f|_r=sum_(s<=r) sup|partial_eta^s f|/s!, and 0<h<=.01. Write

`l=1-2h`, `R1=1/l+4h/l^2`,
`R2=1/l+6h/l^2+16h^2/l^3`.

They bound |1/L|1 and |1/L|2. If u_i>=sup_X|U_i|2,
v_i>=sup_X|v_i|2 on the ENTIRE [0,Xv], define
`B=sum u_i v_(k-i)`, `C=sum v_i v_(k-i)`.
For k>=1 also supply `p>=sup_X|v_(k-1)|3`.
Using |eta|1=2, |d|1=3, |d|2=4, and |f_eta|r<=(r+1)|f|_(r+1),

\[
 |T_{b-D}T_b f|_1
 \le R_1R_2(6+4|b-D|)(12+4|b|)|f|_3.                 \tag{10}
\]

Let V_j=Xv^(j+1)/(j+1), beta=lambda_(k-1)+j, and
`A_kj=nu R1 R2 (6+4|beta-D|)(12+4|beta|) V_j p`
for k>=1; set A_0j=0. Then

\[
 |I_{k0}|_1\le R_1X_v[(\lambda_k+4D)v_k
 +(6+4|\lambda_k-A|)B]+X_v C/2+2\nu v_k+A_{k0}=:Q_{k0},\tag{11}
\]
\[
 |I_{k1}|_1\le {R_1X_v^2\over2}[(1+\lambda_k+4D)v_k
 +(6+4|\lambda_k+D|)B]+X_v^2 C/4+A_{k1}=:Q_{k1}.       \tag{12}
\]

For k=1 this includes B=u0 v1+u1 v0, C=2v0 v1 and
beta=0,1 in the two axial terms. Radial viscosity cancellation does not
authorize setting A_11=0. These are analytic estimates conditional on
the actual full-support norms, not their numerical evaluation.

The regular average formula (1) eliminates axis divisions. For r=2,3,

\[
 \sup_X|v_i|_r\le R_r[4+4(D+\lambda_i)+4(r+1)]
                                  \sup_X|U_i|_{r+1}, \tag{13}
\]

where `R3=1/l+6h/l^2+32h^2/l^3+64h^3/l^4<2`.
The positivity of the displayed rational majorants in h verifies this
bound over the whole interval from its h=.01 endpoint. Averaging over
s in [0,1] does not enlarge the angular norm. No global radial derivative
is required FOR THESE TWO MOMENTS.

Consequently the explicit remaining second-pressure input is the actual
global U1 C_eta^3 norm, and the actual global U0 C_eta^4 norm for the
shifted viscosity. The preceding leading global C_eta^3 estimate alone
does not provide the latter; local positive-collar estimates alone do
not control the complete future support. Finite smoothness constants
from a qualitative construction are not an evaluated bound at the
retained parameter values. This checkpoint does not claim impossibility
of a sharper estimate using additional identities.

## 4. Feed the complete source to the unchanged second repair

The zero-axis-gauge pressure equation is
`Pi2_X=2F0 F2+F1^2-w1/2`. Suppose J bounds the tentatively cut inner
F2,U2 in C_eta^1 on [0,Xcut]. Let f0,u0 bound F0,U0 there, and f1,u1
bound the FINALIZED F1,U1 on their actual support [0,Xvel]. Ordinary
second-order discrepancy norms obey

\[
\begin{split}
\mu_1&\le X_{cut}J,\quad\mu_2\le X_{cut}^2J,\\
\mu_3&\le2f_0X_{cut}J+X_{vel}f_1^2+Q_{10}/2,\\
\mu_4&\le X_{cut}^2(u_0+f_0)J+X_{vel}^2u_1f_1,\\
\mu_5&\le2u_0X_{cut}J+f_0X_{cut}^2J+X_{vel}u_1^2
                         +X_{vel}^2f_1^2/2+Q_{11}/2. \tag{14}
\end{split}
\]

The lower-order F1 square, U1 square and U1F1 product are held fixed.
They are not a current-order quadratic repair. Substitution into the
existing inverse-800 bound gives

`|c2|1 <= 800 max(mu1/rho, (mu1/rho+K W1 mu4)/lambda_patch,
                  K mu2, W1 mu5/rho, W1 mu3)`,

with K=1/(sqrt(2)rho^(3/2)) and W1=4/e_star unchanged.
This controls the repair only AFTER a valid inner J and all global inputs
have been supplied. It does not construct that J, control H2,3 or license
replacing complete pressure by its local derivative on the mean patch.

## 5. Reproduction and limits

`moments.py` reuses the historical Fraction polynomial/angular-jet arithmetic,
not its k=0 flux formula. It directly expands (2) and separately checks
(3),(7)--(9) at k=0,1,3, both full and truncated support, eta=+-1 and
interior angles, with two manufactured C4 endpoint fixtures. One fixture
also has the actual positive-order zero axis traces. Nonzero shifted
viscosity and premature-support omission controls are retained. Independent
floating comparisons use the original `background_summation/hierarchy.py`
at the axis and positive radii. Those comparisons are diagnostics, not
interval certificates or selected PDE coefficient solutions.

Run from this directory:

```
python test_audit.py -v
python audit.py --verify-record evidence.json --output fresh.json
```

The fresh output must not exist. Historical evidence is byte-bound and never
overwritten. Numerical agreement and passing tests are not continuum proof.
Next: quantify the unchanged full-support angular inputs and the entire-axis
first-coefficient jets needed by the second inner system, derive its local
J3 and join exact repair to H2,3. Actual cutoff jets, retained scale/value
acceptance, every physical pulse derivative/support, full-wave covariance
with axial velocity and pressure means, complete viscous/cutoff/interaction
residual and independent checks of inherited analytical premises remain open.
