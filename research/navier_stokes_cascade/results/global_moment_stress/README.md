# Renormalized global moments and signed exterior stress memory

Status: `exact-leading-global-moments-and-signed-stress-identification-audited`.
Parent: `cf7a30cb7dddd82278f9e47031ad97458e3eed45`.
The fixed parameters, support endpoints and historical evidence are unchanged.

This checkpoint derives all five leading global moment identities for a
profile completed by the specified EXACT correction roots. It proves that
the weighted leading tangential residual totals vanish and hence identifies
the forward stresses with the backward stresses used by the outer-edge audit.
It also gives an exact five-defect formula when restoration is imperfect.
These implications retain the compact construction's input hypotheses;
they do not independently certify its remaining regional cone inputs, the
complete leading profile, or any positive-order Navier--Stokes correction.

The [previous checkpoint](../shape_join_pressure/README.md) established
the actual shape/join entry and leading pressure ledger, but left the
other global moments and this signed-stress identification open. In
particular raw I does not have a finite total. All assertions below use
the subtracted power, not a numerical cutoff pretending I(infinity)=0.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.7), (4.11), (4.15)--(4.16), (A.11), (A.34)--(A.46).
Printed page 141 was rendered and visually checked. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its final theorem is not an assumption. The PDF review skill was used to
check the integration signs, radial weights and subtraction on that page.

## 1. Convergence, including angular derivatives and both ends

Write A=1/2+h, D=1/2-h, d=1-eta^2, L=1-2h eta^2,
0<h<.01. The final heat exterior has U=0 and

\[
 E=c_\infty X^{-A}{\cal H}(2d/X),\qquad
 H=\sqrt{2X}E,\qquad H_{pow}=kX^{-h},\quad
 k=\sqrt2c_\infty>0.                                      \tag{1}
\]

The positive constant c_infinity is independent of eta. This is essential
for both the subtraction and the physical scaling argument below.
The Gamma integral gives |Hcal^(j)|<=(h)_j(1+h)_j. Therefore, for X>=R>=100,

\[
 \left|\int_X^\infty(H-H_{pow})\,dx\right|
       \le2k(1+h)X^{-h},\qquad
 \left|\partial_\eta\int_X^\infty(H-H_{pow})\,dx\right|
       \le4k(1+h)X^{-h},                                  \tag{2}
\]
\[
 \left|\partial_\eta^2\int_X^\infty(H-H_{pow})\,dx\right|
 \le kX^{-h}\{4(1+h)+16h(1+h)(2+h)/X\}.                    \tag{3}
\]

For example, differentiate Z=2d/X: Z_eta=-4eta/X and Z_etaeta=-4/X.
The two terms in (3) are both required. The respective integrable density
majorants are O(kh X^(-1-h)) and O(kh X^(-2-h)). For every further fixed
eta or logarithmic-radial derivative the same chain rule gives a finite
sum with integrable bounds. No division by d is made at eta=+-1.

At the regular axis H=2XF with smooth F. The subtracted power kX^-h is
integrable there since h<1, also after angular differentiation. Thus

\[
 R(\eta):=\int_0^\infty(H-H_{pow})\,dx
 =\lim_{X\to\infty}\{I(X)-XH_{pow}(X)/(1-h)\}               \tag{4}
\]

is a well-defined smooth function. The limit and each fixed angular
derivative commute by these majorants. Compact U support gives convergent
M and J; S converges since E^2=O(X^(-1-2h)), while total Cp converges
since E^2/(2X)=O(X^(-2-2h)). These too are functional identities in eta,
not matches at finitely many angles.

One nonuniform limit must be retained. Taylor's formula with its signed
second-derivative remainder gives

\[
 D_I(X):=\int_X^\infty(H-H_{pow})\,dx
 =-2dk(1+h)X^{-h}+\mathcal E,
 \quad0\le\mathcal E\le2d^2kh(1+h)(2+h)X^{-h-1}.            \tag{5}
\]

For fixed X,k,d this tends to -2dk as h decreases to zero, although the
pointwise heat edit tends to zero. Setting h=0 before the infinite
integration loses that finite angular moment. Bounds (2)--(5) require
h>0; no uniform h=0 interchange is claimed.

## 2. Trace the exact five global equations through the edits

Use row order (M,R,J,S,Cp). Here Cp has its positive-density definition;
the pressure defect is Pi0+Cp(infinity), not Cp(infinity) itself.

The reference construction's [pulse moment argument](../pulse_moment_audit/README.md)
supplies the two exact linear axial equations M=J=0 and the exact scalar
quadratic strength equation S(infinity)=0, including both U-square terms.
Its section 3 supplies the two angular equations: the pressure total is
unchanged, and I=XH/(1-lambda) at the subsequent uniform-power boundary.
Its section 4 chooses the finite waiting length using the scalar Q equation.
With U=0 and eta-independent release, that terminal condition is Q=0
beyond the cutoff, so (4.16) gives I=XHpow/(1-h), namely R=0 for the
unedited power tail. The strength equation uses the WHOLE future energy,
including the angular corrections, rather than a truncated S integral.
These explicit reference correction lemmas remain dependencies; their
quadrature residuals are not substituted for their exact defining roots.

| Edit | Exact global action |
|---|---|
| Axial pulse and its two end corrections | Set M=J=0; change S, which the strength root then sets to zero; leave I and Cp unchanged |
| Two angular corrections and terminal waiting choice | Restore the stated I boundary and preserve total Cp; waiting then sets R=0; their S contribution is included in the strength equation |
| Heat edit plus separate heat compensation | M,J unchanged because U=0 on both supports; signed changes in R,S,Cp cancel exactly |
| Axis/shape replacement and the five-row join | Match all five cumulative moment FUNCTIONS at the joining endpoint, hence preserve all five global totals |
| Compact modulation and its exact five-row repair | Restore the full incoming moment memory at the repair endpoint, hence preserve all five global totals |

For the heat row the normalized target signs are (-Delta Cp,+Delta S,
-Delta I) with compensation row -S, as defined in the
[heat audit](../heat_exterior_audit/README.md). Delta R=Delta I for these
edits because c_infinity and the final power subtraction are unchanged.
The finite-dimensional inverse bounds and smallness conditions are
those replayed by the previous checkpoint: angular (5,32,700lambda^30),
heat (50,100,64h^4/Xstar), and the exact joining/compact repair roots.
The compact root still depends on its stated input hypotheses. Supports
of the LEADING edits are fixed compact sets; no positive-order or reserved
mean support claim is inferred here.

Consequently any completed leading profile satisfying these dependencies has

\[
 M_\infty=J_\infty=S_\infty=R=0,\qquad
 \Pi_0+C_{p,\infty}=0.                                   \tag{6}
\]

Together with (4), its exterior cumulative functions are exactly

\[
 M=J=0,\quad I=XH_{pow}/(1-h)-D_I(X),\quad
 S=\frac12\int_X^\infty E^2dx,\quad
 \Pi=-\int_X^\infty E^2/(2x)dx.                            \tag{7}
\]

The heat compensation and heat edit need not have zero partial memory
between their disjoint supports. All equalities in (6)--(7) concern the
exact roots, never finite-precision approximate coefficients.

## 3. Exact exterior defect map: the inequality small errors cannot solve

Consider two profiles with identical local fields U=0,E at all X>=Xb and
the same c_infinity. Let their cumulative differences there be
(m,R,j,s,p), with p the constant pressure difference. Primes below mean
eta derivatives. Define Bm=2Deta*m+d*m'. At profile scale q=1, direct
subtraction of (4.7), (4.11), (4.16) gives

\[
 \Delta V_0=-B_m/L,                                      \tag{8}
\]
\[
 \Delta(r^2T_\theta)
 =\{HB_m+(1-h)R-D\eta R'-dj'+2(h-D)\eta j\}/L,             \tag{9}
\]
\[
 \Delta(rT_z)
 =\{D(m-\eta m')+4h\eta s-ds'
                  +X(4A\eta p-dp')\}/L.                  \tag{10}
\]

Here r=sqrt(2X). No asymptotic approximation is used. The M term in
(9), the angular derivative of pressure in (10), and the pressure's X
factor must all remain. (8) shows why U=0 alone is insufficient to give
V0=0. When p=0, the limits in (9)--(10) are the constant weighted stress
memories (the H Bm term tends to zero). A pressure error can instead grow
linearly in r*Tz. This extends the earlier local memory calculation by
identifying its GLOBAL normalization and its physical consequences.

For any proposed approximate repair, exact exterior support would require
the right sides of (8)--(10) to vanish IDENTICALLY for every eta and X>=Xb.
A bound by 1e-100 is not that equality. At eta=0, a pressure defect p=epsilon*eta
has value zero but contributes -X*epsilon to r*Tz. At eta=1, a constant
R defect contributes (1-h)R/L to r^2*Ttheta. These are exact negative
controls, not roundoff diagnostics. Zero five defects suffice for vanishing;
no assertion that they are the only possible coupled cancellation is needed.

## 4. Signed weighted residuals and the backward representation

Fix q>0 and eta (or equivalently t<1,z) and integrate over physical radius.
The power swirl Kpow=c_infinity (r^2/2)^(-A) is independent of t,z.
The exact profile-to-physical integral conversions are

\[
 \int_0^\infty r^2(u_\theta-K_{pow})dr=q^{1-h}R,\quad
 \int_0^\infty r^2u_zu_\theta dr=q^{3/2-2A}J_\infty,\quad
 \int_0^\infty ru_zdr=q^{1-A}M_\infty.                    \tag{11}
\]

Tonelli's theorem on the nonnegative pressure density and the backward
pressure in (7) give the EXACT signed flux identity

\[
 \int_0^\infty rp\,dr=-\tfrac12\int_0^\infty ru_\theta^2dr,
 \qquad\int_0^\infty r(u_z^2+p)dr=q^{1-2A}S_\infty.        \tag{12}
\]

All these integrals converge for fixed h>0. Time/axial differentiation
commutes with the subtracted angular integral: at infinity the heat
estimate gives r^2(|u_theta-Kpow|+|partial_t u_theta|)=O_h(r^(-1-2h));
the corresponding fixed axial derivatives are also integrable. At the
axis the subtracted weight is r^(1-2h), which is integrable, and Kpow
has no time/axial derivative. Compact U support treats its integrals.
There is no assertion of bounds uniform as h goes to zero or q goes to zero.

For clarity the residual below is LEADING tangential momentum residual:
it includes radial viscosity and omits axial viscosity, just as (4.12).
Integrating its conservative angular/axial equations, the radial transport
boundaries vanish by axis regularity and U=M=0 in the exterior. Radial
viscosity contributes [r^2 partial_r u_theta-r u_theta] and [r partial_r u_z];
these vanish at both ends, the first by O_h(r^-2h) at infinity. Equations
(6), (11)--(12) then give the signed identities

\[
 \int_0^\infty r^2\mathcal R^{(0)}_\theta dr=0,\qquad
 \int_0^\infty r\mathcal R^{(0)}_z dr=0.                   \tag{13}
\]

As a sign/scaling cross-check, with p=0 but nonzero moment defects, the
negative of these totals is q^-h times the limiting expression (9), and
q^-A times (10). Indeed apply the exact coordinate differentiation rules
to q^(1-h)R, q^(3/2-2A)j, q^(1-A)m and q^(1-2A)s; their coefficients are
precisely (1-h), Deta, d, 2(h-D)eta, D and 4heta in (9)--(10).

The stress was defined with a MINUS integral from the axis. Hence (13)
implies the PLUS backward formulas

\[
 T_\theta(r)=r^{-2}\int_r^\infty s^2\mathcal R^{(0)}_\theta(s)ds,
 \qquad T_z(r)=r^{-1}\int_r^\infty s\mathcal R^{(0)}_z(s)ds. \tag{14}
\]

In the exact heat exterior the leading residual is zero, so (14) proves
T0=0 for X>=Xb. On the terminal collar this identifies the forward stress
with the backward formulas whose positivity, orientation and flat-edge
factorization were derived in the [outer-edge audit](../outer_edge_stability/README.md).
Its moment-identification premise is thus supplied for these completed
leading profiles. Its separate positive-order common-support premise is
NOT supplied: the first axial-viscous correction and later coefficients
still require that support audit. Nor does (14) independently verify all
regional cone bounds between the inner and outer collars.

## 5. Checks, scope and next obligation

`ledger.py` compares the defect map with independently assembled original
coordinates using exact rationals, each basis defect, combined defects,
both endpoint angles and several radii. It retains the growing pressure
memory, the M-induced exterior radial velocity and nonzero angular jets.
The signed ledger examples check algebra only; the existence proof is the
exact-root argument in section 2 and its explicitly retained dependencies.

An independent Gamma/Tricomi-U heat diagnostic uses h=1/8, NOT the selected
actual h, at 50/70 digits. Infinite I integrals subtract their first Taylor
term analytically before quadrature. A finite degree-12 Taylor formula at
small Z has an explicit integrated remainder budget; it is not assumed to
be a convergent series. The resulting forward angular stress of the pure
heat tail vanishes within the labelled numerical tolerance. At eta=1,
D_I=0 but I_eta/Hpow=-4(1+h), so deleting angular derivatives would fail.
A separate manufactured pressure fixture checks the Fubini sign in (12).
None of these numerical checks proves a continuum supremum or exact repair.
Historical evidence is never rewritten; the new record adds its own hashes.

The diagnostic's 24/36-node refinement differs by less than 3.04e-28 in
D_I; the 36-node zero-angular-stress discrepancy is below 9.54e-40. The
small-Z Taylor contribution has a separate integrated error bound below
1.65e-58. Quadrature refinement is not a rigorous enclosure of quadrature
error and is not included in the analytical moment proof.

Local validation on 2026-10-16: 213 tests passed across nine suites
(this checkpoint, shape/join, annulus, outer edge, pulse moments, heat
compensation, compact pressure bridge, modulation and repair); their nine
frozen evidence records replayed successfully. This checkpoint adds 31
tests and verifies 521 inherited provenance hashes. The parent commit's
CI succeeded; the new commit's CI status is reported separately.

Next: audit the remaining outer regional cone inputs jointly with the
inner/shape/join bounds. Concretely, on the compact input rectangle outside
the protected edge collars, with the unchanged Aenv=C^4096, the retained
test is a>=Aenv^-1, Pc-2>=Aenv^-1, and, wherever vs>=2,

\[
 P_c-v_s\ge A_{env}^{-1},\qquad
 2(P_c-v_s)^2-(v_s-2)J_c^2\ge A_{env}^{-1}.                 \tag{15}
\]

The regional bounds must apply at the SAME selected parameters and radius,
with the exact ledgers above; a successful Md4 pulse fixture does not repair
its failed axial stage. This checkpoint removes the global integration
constant ambiguity, not the need to verify (15) through the remaining
outer regions. Then audit the positive-order and reserved-mean supports.
Uniform actual-pulse values/derivatives/supports across every physical scale,
complete-wave-covariance mean corrections including axial velocity and
pressure, and the complete residual including viscosity, cutoffs and every
interaction remain open. No blow-up conclusion follows from this checkpoint.

```sh
python3 research/navier_stokes_cascade/results/global_moment_stress/test_audit.py -v
python3 research/navier_stokes_cascade/results/global_moment_stress/audit.py \
  --verify-record research/navier_stokes_cascade/results/global_moment_stress/evidence.json \
  --output global_moment_stress_recheck.json
```
