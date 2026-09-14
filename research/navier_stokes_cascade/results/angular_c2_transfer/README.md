# Actual second-angular-derivative transfer

Status: **`actual-C2-moment-transfer-and-matching-bounded`**.

This checkpoint closes the first estimate left open by the
[joined-stress review](../joined_stress_construction/README.md): transfer of
the actual core's higher angular derivatives through the shape transition,
all five incoming moments, and their joining correction. It is an analytic
extension of the existing core/continuation argument, with outward checks of
the constants. It is not independent peer review or a proof-assistant proof.

In plain language: the previous certificate controlled the profile's angular
slope. This one also controls how that slope changes. A small slope alone
does not guarantee small curvature; the previous counterexample remains valid.

For the fixed Md64 parameters and the same exact, implicitly defined profile:

| Quantity | Uniform bound |
|---|---:|
| Each normalized incoming moment at y=-8, factorial C2 norm | less than \(2.55\,10^{-17}\) |
| Joining target z after axial restoration, factorial C2 norm | at most \(2.86\,10^{-17}\) |
| Raw second derivative of z | at most \(5.72\,10^{-17}\) |
| Each exact joining coefficient, factorial C2 norm | at most \(5.72\,10^{-14}\) |
| Each exact joining coefficient, raw second derivative | less than \(5.721\,10^{-14}\) |

These are dimensionless matching quantities, not measures of distance to
blow-up. The full compact mixed-derivative envelope, modulation and repair
pressure-error constants, selected stress frequency, full PDE corrections,
smooth forcing, and blow-up remain unverified. No fields or parameters are
changed: the stronger bound applies to the same joining correction as before.

Source: the supplied [Navier-Stokes manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
especially (B.4), (B.12), (B.22), (B.26), (B.32), and (B.34)-(B.40).
The B.32/B.34 page was checked visually to retain its angular-only derivative
scope. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent: `86f5b77183ed9f93f8faa3fd4bee548596bf5b53`.
All preceding audit folders and evidence remain frozen.

## 1. Norm and exact core evaluation

Throughout this note, eta ranges over [-1,1], and

\[
\|v\|_2=\sup|v|+\sup|v_\eta|+\tfrac12\sup|v_{\eta\eta}|.
\]

This is a Banach algebra: the product rule, including its factor two, gives
\(\|uv\|_2\le\|u\|_2\|v\|_2\). Consequently
\(\|e^v\|_2\le e^{\|v\|_2}\). The factor 1/2 is essential: a bound R in
this norm bounds the **raw** second derivative by 2R, not R.

Use S, rho, K, Lambda, B, and C exactly as fixed in the
[analytic core attachment](../axis_core_attachment/README.md). In particular
\(\|\Phi\|_\rho,\|u\|_\rho\le S\), and
\(U=U_*+u/\Lambda\), \(U_*=4\eta+j\).
For \(0\le Y\le5\), the infinite coefficient norm gives

\[
\|F(Y,\cdot)\|_2\le\|F\|_\rho
\sum_{k=0}^2\frac{\rho^{-k}}{(k+1)^2}
\sum_{n\ge0}\binom{n+k}{k}(1/4)^n
\le\frac{496}{243}\rho^{-2}\|F\|_\rho.
\]

We discarded the positive denominator (n+1)^2 to increase the bound and
used the exact generating function (1-r)^(-k-1). Since rho<1, this is below
\(3\rho^{-2}\|F\|_\rho\). Thus define

\[
D_2=3S\rho^{-2},\qquad
\|\Phi\|_2,\|u\|_2\le D_2,\quad
\|U-U_*\|_2\le D_2/\Lambda<\epsilon/100,
\qquad\epsilon=10^{-16}.
\]

This is an all-index estimate on the exact analytic solution, not a finite
Taylor tail test. The last inequality is checked with the actual fixed
logarithmic parameters.

## 2. Carry angular derivatives through the fields

The source defines \(CE=\sqrt{2X}\phi_*\Phi\) in the core. The same complex
rectangle as before gives \(|\zeta_*|\le M\) and, by a Cauchy disk of radius
r0/2, \(|\zeta_*'|\le2M/r_0\). Therefore

\[
\|\log\phi_*\|_2\le L_*:=\Lambda M(2+r_0^{-1}).
\]

The known \(\Phi>.264>1/4\) and \(\|\Phi\|_2\le D_2\), with D2>=1, imply
\(\|\log\Phi\|_2\le20D_2^2\): the reciprocal derivative terms contribute
at most \(8D_2+8D_2^2\), and the value term at most \(\log D_2+2\).
It follows that

\[
\|CE/\sqrt X\|_2\le\sqrt2 e^{L_*}D_2\le e^B.
\]

For the continuation on [Xa,Xi], Xa=4/Lambda and Xi=110, the frozen reference
argument already supplies

\[
\|p_{1,r}(X)\|_2\le10^6X\Lambda K^8,\qquad
\|n_{s,r}(X)\|_2\le N_2:=10^6K^8.
\]

These are the C2 consequences of its explicit C3 bounds, not values sampled
from a manufactured reference. All cutoffs and kappa depend only on radius.
Before the last slope interpolation, \(D_X\log\phi=-\kappa p_{1,r}/2\);
in that interpolation the slope varies convexly to -.4. With
\(L_0=\log(110\Lambda/4)\), integration gives the uniform majorant

\[
\|\log(CE)\|_2\le
L_*+20D_2^2+55\,10^6\Lambda K^8+2L_0+10<B.
\]

Here the p1 bound was integrated as X d(log X)=dX. Replacing it by its
unweighted supremum over the entire logarithmic interval would introduce
an unnecessary large factor. The 2L0+10 term covers the parameter-independent
slope and the scalar factor log(sqrt(2X)). Thus \(\|CE\|_2\le e^B\) on
the continuation, and in particular \(\|\ell_i\|_2<B\).

For U, integrate \(D_XU=-\kappa Xn_{s,r}/2\), with the final cutoff factor
between zero and one. On the first interval of width delta=C^-20 use
kappa<=1 and e^delta-1<=2delta. Thereafter kappa<=delta. Since Xa<=4 and
Xi=110, the total change from U(Xa) is at most

\[
\tfrac12N_2(2X_a+X_i)\delta\le100N_2\delta
=10^8K^8\delta<\epsilon/100.
\]

Only angular derivatives have been taken, so no inverse cutoff width enters
this calculation. This statement would be false for unrestricted radial
derivatives. Hence the actual terminal axial offset satisfies

\[
\|g\|_2=\|G_i-4\eta\|_2
\le j+D_2/\Lambda+10^8K^8\delta\le3\epsilon/100=:g_b.
\]

In (B.34), \(0\le y_i\le T_{sh}=2000B\). Since
\(\|\log f\|_2<3\), where \(f=(1+\eta^2)^{-1}\),

\[
\|\log(CE)\|_2
\le200B+(1-\sigma)B+3\sigma\le204B,
\qquad \|CE\|_2\le e^{204B}.
\]

The exponential algebra estimate retains the (log E)_eta^2 term in E_etaeta.
The factor C^-1 stays outside this estimate; taking the norm of log E itself
and then exponentiating would lose the small amplitude.
Together these arguments give \(\|U\|_2<10\) and
\(\|E\|_2\le C^{-1}e^{204B}<1\) throughout [0,Xsep].

## 3. Integrate the five actual moments, including the axis pressure

Put \(q=x_{sep}=X_{sep}/X_R=e^{2000B}(CP)^{-10}<e^{-8}\) and x=X/XR.
Hats remove the radial factors M=XR Mhat, I=XR^(3/2) Ihat,
J=XR^(3/2) Jhat, and S=XR Shat. Let \(E_b=C^{-1}e^{204B}\).
The original densities give

\[
\|\widehat M\|_2\le10q,\qquad
\|\widehat I\|_2\le\frac{2\sqrt2}{3}E_bq^{3/2},\qquad
\|\widehat J\|_2\le10\frac{2\sqrt2}{3}E_bq^{3/2},\qquad
\|\widehat S\|_2\le101q.
\]

These follow by integrating U, sqrt(2x)E, U sqrt(2x)E, and U^2-E^2/2,
using the C2 algebra inequality. This also justifies differentiating these
integrals twice in eta.

For Cp, do not discard the axis factor. On [0,Xa],
\(\|E\|_2\le C^{-1}\sqrt X e^B\), so

\[
\left\|\int_0^{X_a}\frac{E^2}{2X}\,dX\right\|_2
\le2C^{-2}e^{2B}.
\]

On the rest, the logarithmic length is L=L0+2000B. Thus

\[
\|C_p\|_2\le C^{-2}e^{408B}(L/2+2)
\le C^{-2}e^{412B}.
\]

The last length inequality is checked outward. A uniform E bound without
sqrt(X) would instead give a divergent integral at zero and would not prove
this estimate. Angular derivatives retain the same integrable axis factor.

## 4. Normalize all five rows before claiming a small target

Use exactly the previous five normalized actual-minus-ideal rows, with
xc=e^-6. In this norm

\[
\|4\eta\|_2=8,\quad \|8\eta\|_2=16,\quad
\|16\eta^2\|_2=64,\quad
\|f^{-1}\|_2=5,\quad\boxed{\|f^{-2}\|_2=20.}
\]

The old C1 value for the last factor was 12. At eta=1 its value, first
derivative, and half-second-derivative are 4,8,8, respectively. Reusing 12
would be an incorrect C2 estimate. The derivative-aware row evaluator and
its independent differentiation check retain this normalization dependence.

At Xsep the rows satisfy the following bounds. Each expression bounds its
complete C2 norm, not only its value.

| Row | Bound |
|---|---|
| M | \(18q/x_c\) |
| J-4 eta I | \(60E_bq^{3/2}/(P x_c^{8/5})\) |
| I | \(10E_bq^{3/2}/(3P x_c^{8/5})+(5/8)(q/x_c)^{8/5}\) |
| -(S-8 eta M) | \(6500q/(P^2x_c^{6/5})+(5/12)(q/x_c)^{6/5}\) |
| Cp | \(20C^{-2}e^{412B}/(P^2x_c^{1/5})+(5/2)(q/x_c)^{1/5}\) |

For the second row the ideal J=4 eta I cancels exactly. For the fourth,
the norm of the remaining nonswirl numerator is bounded by
(101+16*10+64)q=325q, then multiplied by the reciprocal norm 20.
The ideal swirl contributions cancel their angular normalization exactly.
These operations include the derivatives of 4 eta I and 8 eta M.

Using P>=16, Eb<=1, and q<1, every row is bounded by

\[
10^6(q^{1/5}+C^{-2}e^{412B})
\le2\,10^6e^{412B}C^{-2}<\epsilon/4.
\]

All fixed coefficients in the preceding table are checked below 10^6;
the last inequality uses the actual log C=10000B+100 log(1/epsilon).
This is an integral estimate on actual fields, not an inference from the
small radius or the earlier C1 ledger.

From Xsep to y=-8 the swirl is exactly ideal and U=4 eta+g. The increments
of the five rows therefore have C2 bounds

\[
\left(e^{-2}g_b,\ \tfrac58e^{-16/5}g_b,\ 0,
\ \tfrac{20}{256}e^{-4/5}g_b^2,\ 0\right).
\]

In the energy row the exact identity
\((4\eta+g)^2-(4\eta)^2-8\eta g=g^2\) is used *before* estimating.
Adding the core budget epsilon/4 puts each row below .255 epsilon in C2.

## 5. Same exact correction, now bounded in C2

During axial restoration, the two linear rows gain at most gb each and
the energy row at most gb^2. The coefficient
\(\zeta=e^{6/5}f^{-2}/P^2\) has C2 norm at most
\(20e^{6/5}/256<1\). Thus the actual normalized target satisfies

\[
\|z\|_2\le .255\epsilon+g_b+g_b^2<.286\epsilon=:d_2.
\]

The original linear inverse remains bounded by beta=1000 in C2 because it
is independent of eta. For the quadratic map, use disjoint bump supports,
their exact mass .08, and beta_bump^2<=64 beta_bump. The same C2 algebra
argument, with the new zeta norm, bounds its bilinear norm by 100. No small
U-square term or normalization derivative is discarded.

The map \(c\mapsto B_{\rm mat}^{-1}(z-Q[c,c])\) is a contraction on the C2
ball of radius \(2\beta d_2=5.72\,10^{-14}\), since
\(8\beta^2(100)d_2=2.288\,10^{-8}<1\). That ball lies within the older
C1 uniqueness ball. Its root is therefore **the same exact correction**, not
a replacement fitted to new numerical data. Smoothness is inherited from
the original implicit-function argument.

In particular \(\|z''\|_\infty\le2d_2=5.72\,10^{-17}\). Substituting this
actual bound into the previous second-derivative equation gives

\[
\|c''\|_\infty\le
\frac{1000}{1-200000r_{\rm old}}
\left(5.72\,10^{-17}+700r_{\rm old}^2\right)
<5.721\,10^{-14},\qquad r_{\rm old}=6\,10^{-13}.
\]

This uses the exact equation
\(L_c c''=z''-2Q[c',c']-4\dot Q[c,c']-\ddot Q[c,c]\), retaining both
explicit derivatives of zeta. Here r_old is the older coefficient bound,
not the complex-strip radius r0 from section 2. The implementation reuses the
previous transfer helper with the newly supplied actual target bound.

## 6. Critical review and scope

The review checks four risks independently of the supremum arithmetic:

1. All five row jets agree with direct second differentiation, including
   at eta=+-1. Freezing angular normalizations fails a retained control.
2. The energy-row cancellation holds exactly in rational value/first/second
   jets. Dropping its mixed term fails.
3. Exponentiating the shape function retains both second derivative terms.
   Dropping the first-derivative square fails.
4. The axis sqrt(X) factor makes Cp integrable with its angular derivatives.
   Dropping it leaves a logarithmically divergent majorant.

These manufactured controls test bookkeeping; they do not establish actual
uniform bounds by sampling. Those bounds are the analytic argument above,
which depends on the frozen core existence and reference estimates. Earlier
files still report their historical open-C2 status; they have not been
rewritten to claim they already proved the new result.

Only `actual_target_C2_transfer` is closed here. The next obligation is the
**whole compact mixed-derivative and denominator envelope**, followed by the
actual modulation and repair pressure-error constants. The symbolic frequency
from the previous note remains a proposal. This result does not certify
the later stress repair, the full PDE solution, forcing, or blow-up.

Reproduce from the repository root with the pinned mpmath dependency:

```sh
python research/navier_stokes_cascade/results/angular_c2_transfer/test_audit.py -v
python research/navier_stokes_cascade/results/angular_c2_transfer/audit.py \
  --verify-record research/navier_stokes_cascade/results/angular_c2_transfer/evidence.json \
  --output angular_c2_transfer_reproduced.json
```
