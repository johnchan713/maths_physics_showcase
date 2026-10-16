# Fixed-width reduced continuation and both final tapers

Status: `reduced-continuation-relaxed-margins-audited-for-reference-datum`.
Parent: `ed1bafb7a24150b8e9353ea25d7c8487e29bba10`.

For the same core, fixed reference axis datum, and delta=kappa0=t1=C^-20,
this checkpoint supplies the missing continuation inequality (8) in the
[compact pressure bridge](../compact_pressure_bridge/README.md). On
Jcont=[Xa*exp(delta),110], the unmodulated leading profile satisfies

\[
 a>C^{-22},\qquad v_s\le4/5<1,\qquad P_c>11/5.             \tag{1}
\]

Thus the growth>=2 subset of this interval is empty; the corresponding
third/fourth compact-input conditions are vacuous there. The loop is still
needed to produce positive pulse growth. This does not independently close
the subsequent shape/join/outer regional arguments or completed global
pressure normalization. Unconditional actual pulses, full residual
improvement, and blow-up remain unverified. Earlier evidence is unchanged.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.9), (4.16), (B.22)--(B.29), printed pages 150--153. Page 153 was
visually checked. SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The manuscript's theorem is not a premise. This is an analytical audit
with reproducible bookkeeping checks, not independent peer review or a
proof-assistant certificate.

## 1. Reference cutoff with no inverse-width derivative assumption

Use the exact core and constants proved in the
[core contraction](../core_contraction_audit/README.md) and
[core growth](../core_growth_activation/README.md) checkpoints.
Write Aphys=.5+h, D=.5-h, d=1-eta^2, L=1-2h*eta^2>=.98,
j=1e-18, sigma*=1e-20, Xa=4/Lambda, Xcut=100, Xi=110,
y=log(X/Xa), Ustar=4eta+j, Hstar=Deta+dUstar.
Set uerr=Ubar/Lambda and e=Eb*Fb/Lambda as in core growth.
These bound individual native u and Phi error jets at the required orders.

Let psi be one up to y=delta, between zero and one on [delta,2delta],
and zero afterward. The reference is defined by

\[
 (\log\phi_r)_y=\psi(\log\Phi_{nat})_y,\qquad
 (U_r)_y=\psi(U_{nat})_y,\qquad E_r=C^{-1}\sqrt{2X}\phi_r.   \tag{2}
\]

Here phi_r starts at phi_star*Phi(4). The native solution used in (2)
exists since 4exp(2delta)<4.1. The reference phi_r and U_r, not E_r,
are constant after 2delta. Differentiating the integral only in eta gives
|Ur-Ustar|, |Ur_eta-4|<=2uerr, |Ur_y|<=5uerr everywhere on the
reference interval. Radial averaging is a contraction, so

\[
 |W_r-Wstar|\le4uerr,\quad |H_{c,r}-Hstar|\le2uerr,
 \quad H_c=D\eta+dU.                                      \tag{3}
\]

W uses the radial average of U; Hc uses U itself. They must not be
identified. The native mixed log-Phi jet is at most 10K^2. Comparing
the integral in (2) with the native endpoint at min(y,2delta) gives

\[
 |Hstar(\log\phi_r-\log\phi_star)_\eta|
 \le500\chi+2000e+200\delta K^2,
\]
\[
 |(\log\phi_r-\log\phi_star)_\eta|
 \le500/\sigma_*+20\delta K^2.                             \tag{4}
\]

The 200 follows from |Hstar|<6, length<=2delta, and 10K^2;
20 bounds the unweighted integral. Also
ell_r=1+psi*Dy log Phi_nat<=1 and |ell_r|<6, because the exact
native shear is positive and |Dy log Phi|<=4chi+32e.
No radial derivative of psi is used in (2)--(4).

For subsequent integration use the factorial angular norms through order
three. The following coarse reference envelopes follow directly from
(2), the core jets, and the original forward moments, not from a sampled
profile. Put J=1+Lambda+K and keep B=1e20*J^64.

| Input or operation | Sufficient bound or reason |
|---|---|
| Native log Phi jets, U jets, fixed eta coefficients, Lambda*zeta jets through the needed finite orders | J^12 |
| log(CEr), away from X=0, through three eta derivatives | 10J^12; the cutoff integral has length <=2delta |
| CE_r and its reciprocal on Xa<=X<=110 | exp(B/20), including three eta derivatives |
| Axis Cp integral | use CEr=sqrt(2X)*phi_r before integrating E_r^2/(2X) |
| C^2 Cp through three eta derivatives | exp(B/5) |
| Sources and positive radial integrating factors | products of the preceding jets and exp(6log(110/Xa)) |
| p1,r, p1,r^-1, ns,r through two eta derivatives | exp(B) |

For the last row the inverse is taken only after the positive lower bound
below. Differentiating a reciprocal twice uses at most three inverse
factors. The integrating factor's eta derivatives require at most two
factors of J^12*log(110/Xa), not an exponential of that derivative norm.
These products fit inside exp(B); in particular 10J^12 and every
logarithmic integration length used here are below B/100.
For the axial source estimate below the sharper bounds
C^2*(|Cp|+|Cp_eta|), C^2*E_r^2<=exp(B) suffice. Pressure starts
from the unchanged Pi0, including its parameter derivatives.

## 2. Positive reference sources and growth on the whole reference interval

Use the original source, not the physical shear of the flattened reference:

\[
 S_{q,r}=-W_r\ell_r-h(1-2\eta U_r)
             -H_{c,r}(\Lambda\zeta_*+(\log\phi_r-\log\phi_*)_\eta).
\]

The native estimates give -Wstar>2.9, |Wr|<6, |Ur|<4.2, and
|zeta*|<=sqrt(chi)/sigma*. Equations (3)--(4) imply

\[
 S_{q,r}\ge\Lambda L\chi+2.8-524\chi
       -2\bar U\sqrt\chi/\sigma_* -2192e-4uerr
       -1000uerr/\sigma_*-200\delta K^2-40\delta K^2uerr.   \tag{5}
\]

Absorb 524chi in .02Lambda Lchi. Young's inequality bounds the square-root
loss by .03Lambda Lchi+35Ubar^2/(Lambda*sigma*^2). Hence it suffices that

\[
 E_{ref}=2192e+4uerr+1000uerr/\sigma_*
       +35\bar U^2/(\Lambda\sigma_*^2)+240\delta K^2<.1.    \tag{6}
\]

This holds analytically for the frozen parameters: if H=K^(1/16), then
Ubar,Eb*Fb,sigma*^-1<=H/1e6. The first four terms are at most
.003H^4/Lambda<3e-47. Since K<C and C>1e9, the last is at most
240C^-18. We obtain, with slack,

\[
                 S_{q,r}\ge2.4+.94\Lambda L\chi.           \tag{7}
\]

Retain every term in the axial source:

\[
 S_{n,r}=-W_rU_{r,y}-Aphys(1-2\eta U_r)U_r-H_{c,r}U_{r\eta}
         -d\Pi_{r\eta}+4Aphys\eta\Pi_r+\eta E_r^2.
\]

Subtract Zstar with the SAME Pi0. The three velocity differences cost
at most 30uerr, 22uerr and 22uerr. Pressure and E_r^2 cost at most
5C^-2 exp(B). Thus

\[
 |S_{n,r}-Zstar|\le100uerr+5C^{-2}e^B<j/100.               \tag{8}
\]

The final inequality follows from Lambda=1e44*K^8 and
log C=10000B+100log(1/epsilon), epsilon=1e-16.
It is checked outward at the frozen scales, not by evaluating C.

Now use p1,r_y=X*Sq,r/L-ell_r*p1,r and ns,r_y+ns,r=Sn,r/L.
Positive integrating factors and the initial exact core give

\[
 p_{1,r}\ge p_{1,r}(0)e^{-y}+1.88\chi(e^y-e^{-y}),
 \quad p_{1,r}\ge p_{1,r}(Xa)Xa/X+1.2(X-Xa^2/X)\ge X/10.  \tag{9}
\]

On chi>.99 replace the first initial value by 2.3. The lower comparison
increases for y>=0, since 2*1.88*.99>2.3. Thus p1,r>2.3 there.
On chi<=.99 the core separation is |Zstar|>jP^2/100>=2.56j.
The initial ns,r has the sign of Zstar and magnitude >=1.28j;
(8) and the positive forward ODE preserve that sign and |ns,r|>j.
With the reference envelopes, this gives

\[
 p_{2,r}^2/p_{1,r}\ge16j^2C^2e^{-3B}/\Lambda^2>2.3.
\]

Combining sectors, and bounding the same products from above,

\[
 2.3<v_r=p_{1,r}+p_{2,r}^2/p_{1,r}
       \le C^2e^{6B}<C^3/2,\qquad |p_{2,r}|\le Ce^{2B}.    \tag{10}
\]

At X>=100, (9) gives p1,r>119, using Xa<1. This wide margin will
protect both final tapers. None of these comparisons divides by chi.

## 3. Small comparison errors through the two narrow final tapers

After first activation, keep a=delta*p1,r and
Uy=-delta*X*ns,r/2 up to Xcut. On [100,100exp(delta)] multiply
Uy by beta in [0,1], decreasing to zero. On the next logarithmic interval
of width delta keep Uy=0 and interpolate a convexly from delta*p1,r
to .8. Denote the second endpoint by Xend=100exp(2delta)<110.

Until Xend compare actual fields and ALL forward moments to the reference
of section 1. In the factorial angular C2 norm, direct integration gives
|Delta log E|, |Delta U|<=delta*exp(3B). The first two native/cutoff
intervals have total length 2delta; later prescribed slopes are at most
delta*exp(B), integrated over length <B. The a interpolation adds at most
delta*(1+delta*exp(B)). No beta_y or inverse width is used for these
value/angular estimates. In particular the positive exponential change
has norm below two.

Here is a conservative complete error ledger, removing delta in each row:

| Quantity | Bound |
|---|---|
| Delta log E, Delta U | exp(3B) |
| C Delta E, C Delta H, Er/E-1 | exp(6B) |
| Delta M, C Delta I, C Delta J, Delta S, C^2 Delta Cp | exp(12B) |
| Delta W, C Delta NQ, Delta Ns | exp(20B) |
| Delta p1 | exp(30B) |
| Delta p2/C | exp(40B) |

For example J uses Ur*Delta H+Hr*Delta U+Delta U*Delta H;
S uses 2Ur*Delta U+(Delta U)^2-(2Er*Delta E+(Delta E)^2)/2.
Cp integrates the last swirl numerator in y with coefficient 1/2.
The supports begin at Xa>0, so the Cp difference has no divergent axis
integral. The pressure contributions to Delta Ns are
4Aphys*eta*Delta Cp-d*Delta Cp_eta, and -Delta(WU) is kept.
Positive denominator differences use the actual perturbed E,H.
The last row follows from
Delta p2=(X/L)*(Delta Ns/E-Ns,r*Delta E/(Er*E)).
The reference envelopes and integration length <B supply the displayed
exponents, including the additional eta derivative in (4.16).
The field and moment rows are angular C2 bounds; the W/NQ/Ns row is
through angular order one, and the p-coordinate rows are value bounds.
No unavailable third moment derivative is implicitly taken.

Put Dcmp=C^4*exp(100B)<C^6. Exact shear cancellation, valid also when
the axial factor beta is active, gives with r=Er/E

\[
 t_s=\beta p_{2,r}r/p_{1,r},\quad
 v_s=\delta[p_{1,r}+\beta^2p_{2,r}^2r^2/p_{1,r}],
\]
\[
 |P_c-(p_{1,r}+\beta p_{2,r}^2/p_{1,r})|\le Dcmp\delta,
\quad |v_s-\delta(p_{1,r}+\beta^2p_{2,r}^2/p_{1,r})|
          \le Dcmp\delta.                                \tag{11}
\]

For clarity the Pc error is exactly
Delta p1+beta*(p2,r^2/p1,r)*(r-1)
+beta*(p2,r/p1,r)*r*Delta p2. The beta and beta^2 factors are
different; the common delta cancels before estimating t_s.
Every product in (11) is below C^2*delta*exp(50B); Dcmp leaves slack
for sums. The bound 10Dcmp*delta<.01 follows from the frozen log C.

Before the axial taper, (10)--(11) give Pc>2.29 and vs<1/2.
During that taper the comparison Pc is >=p1,r>119, so Pc>3 and
vs<1/2. During the a taper, t_s=0, Pc=p1>3 and vs=a<=.8.
Also throughout these intervals a>=delta*X/10>C^-22, or a is a
convex interpolation with .8. The tiny positive floor is not rounded away.

## 4. The final constant-shear segment needs a source barrier

After Xend, a=.8 and Uy=0, so vs=.8, Pc=p1 and ell=.6.
It is invalid to extend the preceding O(delta) reference-field comparison
through this segment. The reference E has logarithmic slope .5, while
the actual E has slope .1. Exactly,

\[
 E/E_r=(E/E_r)(Xend)(X/Xend)^{-2/5}.                       \tag{12}
\]

At a fixed ratio X/Xend=1.1, the additional log error is
-(2/5)log(1.1), independent of delta. Shrinking the two tapers does
not shrink this error. This corrects that possible overextension of the
width-small comparison; no old frozen record is rewritten.

Instead compare the angular source directly with ell fixed at .6.
After the tapers, U is constant in X and log E_eta is constant in X.
Their differences from the frozen reference values remain at most
delta*exp(3B). Their radial averages obey the same supremum bound.
The source discrepancy from these changes is at most delta*exp(10B)<.1.
In the reference angular-gradient calculation (5), replace the constant
-Wstar by .6*(-Wstar)>1.74 and omit the native radial-slope loss.
The h loss is <.1 and the remaining absorption costs at most Eref<.1.
Consequently the actual Sq>1 on [Xend,110]. No axial or pressure
term is deleted from the integrated coordinates to obtain this estimate.

At a proposed first downward crossing p1=2.5, its ORIGINAL ODE gives

\[
 (p_1)_y=XS_q/L-.6p_1>100-1.5>0,                         \tag{13}
\]

contradicting such a crossing from p1(Xend)>3. Hence p1>2.5 up to
Xi. Together with section 3 this proves (1), and Pc-2>1/5 and
a>C^-22>C^-4096 give exactly the requested continuation input floor.

## 5. Scope, checks, and next obligation

This closes the post-activation continuation/taper gap in the parent's
equation (8), relative to the already audited exact core and its fixed
reference datum. It does not by itself verify every input on the longer
compact interval: the shape transition, earlier joining correction,
outer regional bounds, completed backward-pressure normalization, and
reserved supports must still be independently connected. The compact
modulation/repair theorem and pulse value bounds remain conditional on
those remaining inputs. No uniform actual-pulse derivative/support or
complete Navier-Stokes residual conclusion is promoted.

The checks comprise outward logarithmic budgets at 100/160 digits for
P16, Md4 and Md64, exact rational taper identities, a positive-ODE
comparison control, and original forward/source ODE checks on labelled
manufactured power profiles. Both pressure terms are individually tested.
The Md4 diagnostic does not repair its historical outer-cone failure.
Finite diagnostics are not continuum supremum certificates.
The source controls also reject replacing local U by its radial average
in Hc, and reject losing the half-offset in ell=.5+Dy log E.
The logarithmic sum routine replaces exp(z) for z<-200log(10) by
the outward interval [0,1e-200]. This avoids materializing a nested tiny
width while keeping its contribution; it is not a zero assignment.
For intervals with large absolute width it instead encloses the log-sum
between max(lower endpoints) and max(upper endpoints)+log(number of terms).
Exponentiating a wide difference of two nested logarithms is unnecessary
and can be computationally infeasible even at high relative precision.

Validation: 20 new and 202 inherited tests pass across ten relevant suites;
all ten evidence replays pass and 497 inherited source hashes match.
The suites cover the original attachment, core contraction/growth, inner
activation, compact envelope/bridge, loop, repair, outer edge, and this
continuation. Test success does not certify the remaining global hypotheses.

Next: audit the shape/join and global pressure/moment normalization at the
same parameters, then reserved coefficient/mean supports and physical
pulse/residual obligations. Keep (12) separate from width-small estimates.

```sh
python3 research/navier_stokes_cascade/results/reduced_continuation_margins/test_audit.py -v
python3 research/navier_stokes_cascade/results/reduced_continuation_margins/audit.py \
  --verify-record research/navier_stokes_cascade/results/reduced_continuation_margins/evidence.json \
  --output reduced_continuation_margins_recheck.json
```
