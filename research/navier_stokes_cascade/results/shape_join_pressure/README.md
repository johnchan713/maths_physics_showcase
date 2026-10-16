# Shape barrier, actual C2 join entry, and the leading pressure ledger

Status: `shape-join-entry-and-exact-leading-pressure-ledger-audited`.
Parent: `a8371b14e525bbb03656909f90120221ae63009f`.
All parameters, historical evidence, and exact corrections are unchanged.

This checkpoint extends the independently reviewed inner construction
through its shape transition and into the actual five-moment joining ball.
It derives Sq>1 and p1>2.5 before restoration, and shows that the same
exact joining coefficients have factorial C2 norm <=5.72e-14. The older
conditional annulus bounds therefore apply to these actual incoming data.
It also derives the leading backward centrifugal pressure identity from
the signed, exact Cp ledger, for profiles completed by the stated repairs.
This last conclusion uses the explicitly defined
exact angular, heat, joining, and compact repair roots, not their numerical
quadrature approximations.

The exterior cone inputs, all-order signed stress/moment normalization,
reserved supports, uniform physical pulses, complete residual improvement,
and blow-up remain separate obligations. Finite tests are not continuum
proofs, independent peer review, or a proof-assistant certificate.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.15)--(4.16), (A.11), (B.32)--(B.39), (C.17). Printed page 156
was visually checked. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its final theorem is not a premise.

## 1. Endpoint data and the entire shape source

Use the exact core and the [reduced continuation](../reduced_continuation_margins/README.md).
Write Xi=110, epsilon=1e-16, delta=C^-20, Tsh=2000B,
ell_i=log(CE(Xi,eta)), Gi=U(Xi,eta), f=(1+eta^2)^-1.
The factorial angular norm is ||v||2=sup|v|+sup|v'|+sup|v''|/2.
It is a product algebra; the half is not optional.

The endpoint bounds follow from that checkpoint's reference integrals,
not from a small field difference extended past the final taper:

\[
 \|\ell_i\|_2<B,\qquad
 \|G_i-4\eta\|_2\le j+3uerr+\delta e^{3B}
                     <3\epsilon/100=:g_b.                 \tag{1}
\]

Here uerr=Ubar/Lambda bounds individual native error jets and j=1e-18.
The factor three covers the three factorial C2 components of the native
error. Before the final constant-shear segment the reference log(CE)
norm is at most 10J^12 and its difference is <=delta exp(3B),
J=1+Lambda+K. The final segment adds only the scalar
-(2/5)log(Xi/Xend), of magnitude <1. Thus the first bound follows from
10J^12+1+delta exp(3B)<B. Its angular derivatives do not acquire that
scalar term. Neither endpoint estimate discards the final segment.
For the second estimate, start at the exact core value at Xa and
integrate the actual prescribed Uy directly; its C2 integral costs at
most delta exp(3B). This avoids multiplying a value-only cutoff estimate
by an unjustified C2 norm constant.

For yi=log(X/Xi), prescribe the actual shape field

\[
 \log E=-\log C+yi/10+(1-s)\ell_i+s\log f,
 \quad U=G_i,\quad s=\sigma(yi/Tsh).                       \tag{2}
\]

Since |sigma'|<=64 and 64(B+3)/(2000B)<.05, its logarithmic coefficient
ell=.5+Dy log E is in [.55,.65]. Therefore a=2-2ell is in [.7,.9],
bs=0, vs=a, and Pc=p1. Radial derivatives of narrow cutoffs are not
asserted small.

For the source, keep Hc=Deta+dU, W=1-2Deta*(M/X)-d*M_eta/X,
Aphys=.5+h, D=.5-h, d=1-eta^2 and L>=.98. The preceding reference
calculation and (1) give, uniformly along the shape interval,

\[
 -H_c\ell_i'\ge .95\Lambda L\chi-E_{old},\qquad
 E_{old}\le E_{ref}+\delta e^{10B}<1/100.                  \tag{3}
\]

Eref is exactly the positive source-error ledger in section 2 of the
parent. The leading term is -Hstar*Lambda*zeta*=Lambda Lchi.
Its sqrt(chi) error is absorbed by Young's inequality before taking the
minimum; there is no division by chi. The discrepancy of the endpoint
angular derivative from the frozen reference is <=delta exp(3B).
Multiplication by Lambda*zeta*, bounded by exp(B), fits the additional
term in (3). U and its radial average remain close to Ustar, so -W>2.8.

The new gradient has a different estimate:

\[
 -H_c(\log f)'=2\eta H_c/(1+\eta^2)
       \ge-2j^2-4uerr-2\delta e^{3B}>-1/100.               \tag{4}
\]

Indeed eta*Hstar=(D+4d)eta^2+d*j*eta>=-j^2/(4D), and D>=.49.
This is a uniform square-completion bound, including eta=+-1 and the
tiny neighborhood of zero. The angular gradient in (2) is the convex
combination of (3) and (4). With |h(1-2eta U)|<.1, the ORIGINAL source
obeys

\[
 S_q=-W\ell-h(1-2\eta U)-H_c(\log E)_\eta
       >.55(2.8)-.1-.01-.01=1.42>1.                      \tag{5}
\]

At a hypothetical first downward crossing p1=2.5, its defining ODE gives
p1_y=X*Sq/L-ell*p1>110-.65*2.5>0. The parent supplies p1(Xi)>2.5,
so p1>2.5 throughout the shape transition. After it, E=P f (X/XR)^.1,
U=Gi and ell=.6; the same barrier continues up to axial restoration.
Thus the shape and following ideal-power portion have a>=.7, vs<=.9<1
and Pc>2.5. These are relaxed, not fully admissible, cone inputs.

## 2. Integrate the actual fields before normalizing the join target

Set XR=110(CP)^10, Xsep=110exp(Tsh), q=Xsep/XR, xc=exp(-6).
The same log C=10000B+100log(1/epsilon) gives q<exp(-8), XR>10000.
In the norm of section 1, ||log f||2<3 and the exponential product rule
give ||CE||2<=exp(204B), ||U||2<10 up to Xsep. In the core retain
CE=O(sqrt(X)): ||CE/sqrt(X)||2<=exp(B). Write Eb=C^-1 exp(204B)<1.

Removing XR factors from M,I,J,S by hats, direct positive integration gives

\[
 \|\widehat M\|_2\le10q,\quad
 \|\widehat I\|_2\le(2\sqrt2/3)E_bq^{3/2},\quad
 \|\widehat J\|_2\le10(2\sqrt2/3)E_bq^{3/2},\quad
 \|\widehat S\|_2\le101q.                                 \tag{6}
\]

For Cp the axis contribution is <=2C^-2 exp(2B), not an integral of a
constant E bound divided by X. The remaining logarithmic length is
L0+2000B, L0=log(110/Xa). Hence

\[
 \|C_p(Xsep)\|_2\le C^{-2}e^{408B}[(L0+2000B)/2+2]
                  \le C^{-2}e^{412B}.                     \tag{7}
\]

Pure eta differentiation commutes with these integrals by the displayed
integrable majorants. The ideal comparison primitives are those of the
existing [five-moment annulus](../axis_matching_audit/README.md). The five
normalized actual-minus-ideal rows are M, J-4eta I, I, -(S-8eta M), Cp,
with denominators xc, sqrt(2)Pf xc^(8/5), sqrt(2)Pf xc^(8/5),
P^2 f^2 xc^(6/5), P^2 f^2 xc^(1/5), respectively, in hatted units.
Retain ||f^-1||2=5 and ||f^-2||2=20, including their derivatives.

| Row at Xsep | Factorial C2 upper bound |
|---|---|
| M | 18q/xc |
| J-4eta I | 60Eb q^(3/2)/(P xc^(8/5)) |
| I | 10Eb q^(3/2)/(3P xc^(8/5))+(5/8)(q/xc)^(8/5) |
| -(S-8eta M) | 6500q/(P^2 xc^(6/5))+(5/12)(q/xc)^(6/5) |
| Cp | 20C^-2 exp(412B)/(P^2 xc^(1/5))+(5/2)(q/xc)^(1/5) |

For example the nonswirl numerator in the energy row costs
(101+16*10+64)q=325q; multiplication by 20 gives 6500. This is not
inferred from q being small alone. Each row is at most

\[
 10^6(q^{1/5}+C^{-2}e^{412B})
       \le2\,10^6e^{412B}C^{-2}<\epsilon/4.                \tag{8}
\]

The relation q^(1/5)=exp(400B) C^-2 P^-2 is used before comparison.
There is no evaluation of q, C, or a giant actual profile in floating point.

## 3. The existing exact join accepts these C2 data

Between Xsep and y=log(X/XR)=-8 the swirl is ideal and U=4eta+g.
The added normalized row bounds are

\[
 (e^{-2}g_b,\ (5/8)e^{-16/5}g_b,\ 0,
      (20/256)e^{-4/5}g_b^2,\ 0).                         \tag{9}
\]

The energy identity U^2-(4eta)^2-8eta*(U-4eta)=g^2 holds before
estimating any jets. At entry every row has norm <.255epsilon.
During restoration on -8<y<-7, the two linear rows gain <=gb and the
energy row <=gb^2. Thus the joining target has factorial C2 norm
<.286epsilon=:d2. Its exact quadratic map is Bmat*c+Q_eta[c,c]=z.

The positive bump-mass/Vandermonde estimate of the annulus supplies
||Bmat^-1||<=1000, independently of eta and XR. Here the additive U and
relative E bumps have width .08 and mutually disjoint supports. Their
quadratic rows retain both U^2 and E^2. Since
||exp(6/5) f^-2/P^2||2<=20exp(6/5)/256<1, the same product estimate
gives ||Q||<=100 in C2. Therefore

\[
 8(1000)^2(100)d2=2.288\,10^{-8}<1,\qquad
 \max_j\|c_j\|_2\le2(1000)d2=5.72\,10^{-14}.              \tag{10}
\]

This ball lies inside the original C1 uniqueness ball. The correction is
the SAME exact root, not new coefficients fitted to numerical data. In
particular a raw second derivative is bounded by twice its factorial norm.
All higher fixed derivatives exist by the invertible pointwise implicit map;
no simultaneous small bound at all derivative orders is asserted.

Because these data meet the annulus's stronger entry hypotheses, its
whole-patch inequalities apply: a>.7, vs<1, Q>1, E/P>.2 and Pc>3.32
on -8<=y<=-5 at XR>=10000. All five moment functions, including Cp,
are exactly ideal at the right endpoint. The original forward pressure and
stress coordinates then agree with the outer reference at every later
radius. Derivatives are restored as functions of eta, not at isolated angles.
The argument uses the established finite-dimensional inverse and stress
lemmas; their frozen bounds and equation tests are replayed, not replaced
by numerical Newton residuals.

## 4. Exact leading backward pressure, and the scope of that identity

Let E_id be the unedited exterior schedule and define
T_id(eta)=integral_0^infinity E_id^2/(2X) dX, Pi0=-T_id.
This is the defining positive-mixture datum used by the exact core.
The reference construction uses the exact roots of its finite correction
maps. Their Cp bookkeeping is as follows; the rows refer to TOTAL changes,
not vanishing pressure memory at every intermediate radius.

| Edit | Change of total Cp |
|---|---|
| Two angular corrections (A.11) | 0, the exact pressure row |
| Axial pulse corrections | 0, E is unchanged |
| Heat edit plus its separate compensation | Delta Cp_heat - Delta Cp_heat = 0 |
| Axis/shape replacement plus the join of section 3 | 0, equality of Cp at the joining endpoint |
| Compact modulation plus its exact five-row repair | 0, equality of Cp at the repair endpoint |

Exact existence at the selected outer parameters is supplied by the
[pulse-moment audit](../pulse_moment_audit/README.md), section 3 (inverse<5,
quadratic<32, target<=700lambda^30), and the
[heat-compensation audit](../heat_exterior_audit/README.md), section 2
(inverse<=50, quadratic<=100, target<=64h^4/Xstar). Their explicit smallness
inequalities are replayed as dependencies. No old Md4 numerical pressure
fit is substituted for these exact maps. The last row uses the already
audited compact repair, retaining the inverse-lambda target cost. That
construction still has compact-input hypotheses: (11) applies to any
completed leading profile obtained once those hypotheses hold. This
pressure-identification implication does not prove the remaining outer
cone inputs or a full corrected PDE profile.

At the regular axis, E^2/(2X) is smooth and bounded. In the heat tail,
E=O(X^(-1/2-h)), so this density is O(X^(-2-2h)). Its first two eta
derivatives obey integrable bounds too: the finite intervals are smooth,
and differentiation of the positive heat-kernel integral bounds all fixed
tail parameter jets on [-1,1]. Thus all quantities below and their first
two eta derivatives are well-defined. For the completed LEADING swirl,

\[
 \int_0^\infty E^2/(2X)\,dX=T_id,\qquad
 \Pi(X)=\Pi_0+\int_0^X E^2/(2x)\,dx
             =-\int_X^\infty E^2/(2x)\,dx.               \tag{11}
\]

This identifies the local reference datum with the leading backward
centrifugal pressure for the construction defined by those exact roots.
It does not assert equality for an approximate repair or the full corrected
Navier--Stokes pressure. In general the exact defect is

\[
 \Pi_{forward}(X)+\int_X^\infty E^2/(2x)\,dx
    =\Delta_\infty(\eta):=\Pi_0+\int_0^\infty E^2/(2x)\,dx. \tag{12}
\]

A nonzero fifth-row residual gives that nonzero functional defect; making
it tiny does not make (11) exact. Even a zero value at eta=0 does not
imply a zero derivative there. Between heat compensation and the heat edit
the partial Cp difference is generally nonzero, despite their cancelling
total. This is the pressure memory that must not be erased.

Equation (11) does NOT close every global moment identity. In particular
H=sqrt(2X)E is asymptotic to a positive multiple of X^-h, 0<h<1.
Consequently I(X)=integral_0^X H diverges as X^(1-h); I(infinity)=0
would be false. Its proper next normalization is renormalized I,
I-XHpow/(1-h), together with the separate M,J,S and signed stress
identities. Those all-order/background conditions and reserved supports
remain to be audited. Nor does (11) supply outer cone bounds or physical
pulse derivative/support estimates.

## 5. Checks and next obligation

The new evidence distinguishes the analytical bounds above from exact
rational jet controls, outward logarithmic arithmetic at 100/160 digits,
and labelled manufactured pressure/ledger examples. It tests normalization
derivatives, the energy cancellation, both endpoint angles, an omitted
quadratic pressure term, and nonzero intermediate memory. No giant actual
profile is sampled and no numerical diagnostic is promoted to a continuum
supremum. All inherited evidence remains frozen.

Local validation on 2026-10-16: 368 tests passed across 15 suites
(core growth/contraction, inner activation, core attachment, annulus,
angular C2, angular pulse moments, heat compensation, compact jets,
compact pressure bridge, modulation, repair, outer edge, reduced
continuation, and this checkpoint). All 15 corresponding frozen evidence
records replayed successfully. The new checkpoint contributes 23 tests
and verifies 508 inherited provenance hashes. These are local results;
the published commit's CI result is reported separately, not presumed.

Next: join the retained outer regional cone inputs to this inner/shape/join
certificate, audit renormalized I and global signed stress/moment identities,
then reserved coefficient/mean supports. Uniform actual pulses, complete
covariance mean corrections (including axial velocity and pressure), the
complete viscous residual bound, and blow-up remain unverified.

```sh
python3 research/navier_stokes_cascade/results/shape_join_pressure/test_audit.py -v
python3 research/navier_stokes_cascade/results/shape_join_pressure/audit.py \
  --verify-record research/navier_stokes_cascade/results/shape_join_pressure/evidence.json \
  --output shape_join_pressure_recheck.json
```
