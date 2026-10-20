# Sixth reference pressure input and fifth source transfer

Parent: `7a282097f781d13dc410e64d606109aff676f6c2`.
Status: `sixth-reference-pressure-and-fifth-source-transfer-with-full-memory`.

This checkpoint bounds the SAME reference continuation's fifth angular
sources and integrated coefficients. Its pressure source needs a sixth,
not fifth, angular pressure input; that input is derived here from the
positive-mixture axis datum and regular radial integration. A positive
averaging identity for the retained radial taper avoids requesting a
seventh mixed native derivative. These are conditional analytical estimates
relative to the retained core certificate, not a new existence proof.

This does NOT replace the actual prescribed continuation by the reference.
The actual continuation, shape, same annular joining root and outer prefix
still require fifth-order propagation before the full native compact
Volterra input is accepted. The parent's missing variance-cutoff budget is
also still missing; the explicit RADIAL step is not the variance cutoff.
No roots, gauges, supports, radial width C^-20, parameters or N are changed.
Historical evidence is untouched. Full pulse/covariance/residual closure
and blow-up are not claimed.

Source: the retained [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Lemma A.5, (B.22), and the reference source equations already transcribed
in `higher_jet_budget` and `prescribed_continuation_foundation`.
The local source bytes were rechecked; PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
In particular B.22 explicitly chooses 1-sigma, with the unchanged monotone
squared-flat radial step. Its shape is specified, unlike C.9's variance
cutoff. CI uses repository formulas, not a network PDF download.

## 1. Native packet and a positive sixth-order axis bound

Use factorial norms |v|k=sum_(j=0..k) sup|partial_eta^j v|/j!.
They are product algebras and |partial_eta v|5<=6|v|6. Write K=K_old,
Xa=4/Lambda, Xi=110, Y=Lambda X, F=E/sqrt(2X), and
F_nat=(phi_*/C)Phi, U_nat=U_*+u/Lambda, U_*=4eta+j.
All bounds below are real eta in [-1,1]; h is in (0,.01].

The retained coefficient norm gives, on 0<=Y<=5,

`|Phi|6, |u|6 <= 1+8S(1+rho^-1)^6 <= D0`,
`D0=9(1+S)(1+rho^-1)^6`.

Indeed discard the two squared denominators in the coefficient weight and
sum binomial(n+k,k)(Y/20)^n=(1-Y/20)^(-k-1). The factor at k<=6 is
at most (4/3)^7<8, and sum rho^-k <= (1+rho^-1)^6.
The same sum with an n factor bounds Y u_Y and Y Phi_Y in C5 by
12S(1+rho^-1)^5. Native positivity Phi>.264 is used only on the retained
0<=Y<=4.1 interval, which contains the entire reference taper.

Real reciprocal/logarithm Taylor expansions with 1/Phi<5 give

`|log Phi|6 <= 100000 D0^6`, `|1/Phi|6 <= 546875 D0^6`.

For the logarithm value use |log Phi|<=5D0; for its positive-order jets
use the finite sum of (5D0)^m/m, m<=6. The inverse uses
5 sum_(m=0..6)(5D0)^m. For Y partial_Y log Phi in C5, multiply the
fifth inverse bound 6*5^6 D0^5 by the preceding differentiated numerator
bound. These constants and powers all lie below the SAME

`Nnat=1e20(1+S)^8(1+rho^-1)^48 < K`.

Thus the required pure sixth and radial-slope fifth packets follow from
the coefficient norm, not a relabelled C4 estimate. No seventh derivative
is required for this packet.

For the WHOLE axis datum, let f=(1+eta^2)^-1 and
Pi_axis=-P^2 integral_[0,1] f^(2theta) dmu(theta), with positive mass <5.
The two poles +/-i have distance at least one from the real interval, so
|partial_eta^j log f|/j!<=2/j for j>=1. Each normalized Taylor jet of
f^(2theta), including its value <=1, is majorized coefficientwise by

`exp(sum_(j>=1) 4z^j/j)=(1-z)^-4`.

Consequently, for every finite k,

`|Pi_axis|k <= 5 binomial(k+4,4) P^2`, hence `|Pi_axis|6<=1050P^2`.

This is an analytical positive-measure estimate, not a quadrature fit.
Finite atom/point derivative checks below are diagnostics only.

## 2. The same reference taper is a positive average

Put y=log(X/Xa), psi(y)=1-sigma((y-t1)/t1) on its transition,
extended by 1 before and 0 after. Here t1=C^-20, with
4exp(2t1)<4.1. Set Rr=phi_r/phi_* and ell(y)=log Phi(4exp(y),eta).
B.22 gives partial_y log Rr=psi ell_y, with the SAME native entrance.
Integration by parts yields

`log Rr(y)=psi(y)ell(y)+integral_0^y [-psi'(t)]ell(t)dt`.

The nonnegative weights have total mass one. Once psi=0 the integral
stops changing. The same identity holds for Ur with Unat in place of ell.
The weights are eta-independent, so they commute with every angular jet;
by positivity and the triangle inequality they contract the factorial
angular norm. In particular

`|log Rr|6<=K`, `.264<Rr<=K`, `|Ur-U_*|6<=K/Lambda`, `|Ur|6<K`.

There is NO t1^-1 cost in this pure-angular estimate. Directly bounding
six angular derivatives of the slope instead would request a needless
mixed seventh input. This argument would not apply unchanged to an
eta-dependent or nonmonotone cutoff. Such a cutoff has NOT been substituted.

Value-sensitive finite exponential jets now give

`|Rr|6<=3K^7`, `|Rr^-1|6<=15K^6`.

Indeed sum_(m=0..6) K^m/m!<3K^6, while the value of Rr is <=K
and its inverse value <5. One must not use exp(|log Rr|6).

The retained Cauchy strip for zeta_* gives

`|log phi_*|6 <= Lambda M [2+sum_(j=2..6) 2^(j-1)/(j r0^(j-1))]`
`                  <=10 Lambda K^8`.

This includes the primitive value, since its real integration length is
at most one. The original value |phi_*/C|<=exp(-Lambda) and the previous
Rr value, together with the positive-order log jets, therefore imply

`|Fr|6 <= T := 1e8 exp(-Lambda) Lambda^6 K^49 < .01`,
`|Er|6 <= sqrt(220) T <1`.

The coefficient 1e8 exceeds 3*11^6. The same regular F bound holds on the
native interval from the axis to Xa. This retains its regular integrand,
instead of bounding E^2/X by a singular uniform compact X floor.

Thus the COMPLETE reference pressure is

`Pi_r=Pi_axis+integral_0^X Fr(x)^2 dx`,
`|Pi_r|6 <=1050P^2+110T^2 <=1200P^2<K`.

The original axis datum is not zeroed or changed after tapering.

## 3. Fifth source and integrated coefficient bounds

Write d=1-eta^2, D=.5-h, Aphys=.5+h, L=1-2h eta^2,
M_r=integral_0^X Ur dx, W=1-2D eta M_r/X-d M_r,eta/X,
Hc=D eta+d Ur and lr=1+partial_y log phi_r. The exact source equations are

`Sq=-W lr-h(1-2eta Ur)-Hc (log Er)_eta`,
`Sn=-W Ur,y-Aphys(1-2eta Ur)Ur-Hc Ur,eta`
`   -d Pi_r,eta+4Aphys eta Pi_r+eta Er^2`.

Radial averaging is a contraction also on [0,X]. In C5, the preceding
bounds give W<=30K, Hc<=5K, lr<=10K, Ur,y<=5K,
Ur,eta<=6K, Pi_r,eta<=6K, (log Er)_eta<=100 Lambda K^8.
For the last estimate use positive-order logarithm jets, not the large
constant -log C. Also |d|5=4, |eta|5=2, and |L^-1|5<2:
the geometric series is bounded by sum_(m>=0)(8h)^m<=25/23.
All terms, including BOTH pressure terms and the axial square, are charged.
Safe uniform bounds are

`|Sq,r|5 <= 1e4 Lambda K^10`, `|Sn,r|5 <= 1e3 K^2`.

For the reference native entrance, the exact core cancellation gives
p1r(Xa)=-2Y partial_Y log Phi and nsr(Xa)=-2Ur,y/Xa.
The same packets give |p1r(Xa)|5<=20K and |nsr(Xa)|5<=10K.
Use the exact variation formulas, with no new pressure gauge:

`p1r(X)=[Xa Rr(Xa) p1r(Xa)+integral_Xa^X x Rr Sq,r/L dx]/[X Rr(X)]`,
`nsr(X)=[Xa nsr(Xa)+integral_Xa^X Sn,r/L dx]/X`.

The first initial term is <=900K^14 and its integral term is
<=450000 X Lambda K^23. Since X Lambda>=4, these imply

`|p1r(X)|5 <=1e8 X Lambda K^24`, `|nsr(X)|5 <=1e5 K^3`.

Both are <B_old/100 on Xa<=X<=110 at the unchanged Md64 parameters.
The phi_* factor cancels in the integrating ratio, not in Sq's angular
gradient. No exponentially large Gronwall factor is introduced.

For completeness, at the genuine CORE entrance a=Xa the five moment C5
norms obey

`M<=aK, I<=a^2, J<=a^2K, S<=aK^2+a^2/2, Cp<=a`.

These follow directly from the five original densities, using |F|5<1
and |U|5<K on [0,a]. Each is <K at the retained parameters. They bound
the normalized fifth trace as well, but are NOT the traces at the later
actual compact entry; propagation to that entry is still required.

## 4. Reproduction, interpretation, and next obligation

Run `python test_audit.py -v`, then `python audit.py --verify-record
evidence.json --output FRESH.json` from this directory. The frozen record
binds the prior sources transitively and refuses overwrite. Interval
comparisons at 80 and 110 digits check the actual parameter margins using
logarithms; C,G,H,N and the actual flow are not materialized.

Independent exact rational Taylor fixtures use radial-polynomial regular
F and U, with nonzero fifth fields and sixth axis trace at both angular
endpoints and an interior point. Their five genuine axis-to-X moments
satisfy both source ODEs through all mixed total-order-five coefficients.
Removing either the axis datum or accumulated radial pressure changes the
fifth source. These algebraic diagnostics are not continuum verification.

Next propagate these reference C5 estimates through the ACTUAL prescribed
slopes, shape transition, same earlier annular root and outer prefix, and
carry the genuine fifth entrance memories to the parent's Volterra gate.
Separately supply the unchanged variance-cutoff shape budget (or a valid
actual-composition estimate). Fourth pressure/q5/moving-phase constants,
compact acceptance, FINAL global jets, all physical scale/support bounds,
complete wave covariance including axial/pressure means, and the full
viscous/cutoff/interaction residual remain unproved.
