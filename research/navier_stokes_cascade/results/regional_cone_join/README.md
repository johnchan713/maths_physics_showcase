# Regional cone inputs, corrected outer comparison, and the leading value join

Status: `regional-input-gaps-and-leading-value-cone-joined`.
Parent: `1cbfe7e7ca29321ce510ac6e07e5e24be9b9015c`.
All parameters, widths, exact correction roots and historical evidence are
unchanged. Md=64, T=exp(64)+10, lambda=exp(-4T), h=lambda^2,
XR=110(CP)^10>10000, delta=C^-20 and Aenv=C^4096 are retained.

This checkpoint supplies the remaining REGIONAL VALUE margins for the
compact modulation input and joins its finite-frequency cone to both
protected collars and the heat-compensated exterior. It corrects a
sign-sensitive comparison in the older first outer-ramp argument. The
weaker floor actually needed there, p1>.5XR, survives the correction.

The result is an analytical leading-profile value/cone construction using
the preceding core, finite moment inverse, jet-envelope, loop and repair
lemmas. Their bounds are replayed, not replaced by sampled actual fields.
It is not independent peer review, a proof-assistant certificate, a complete
all-order leading-profile theorem, or an unconditional physical pulse bound.
Positive-order and reserved-mean supports, all-scale physical derivatives,
complete covariance mean corrections and the full residual remain open.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.9), (4.16), (4.20)--(4.23), Lemma 4.5, Appendix C.
The PDF review skill was used to render and inspect printed page 32,
checking the finite-radius cone and orientation formulas. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The manuscript's final theorem is not an input.

## 1. Exact input tests and the common domain

Put t=-b/a, v=a+b^2/a, Pc=p1+t*p2, Jc=p2-t*p1. On the unmodulated
input rectangle [Xminus,Xrep] x [-1,1], we need

\[
 a\ge A_{env}^{-1},\quad P_c-2\ge A_{env}^{-1},\qquad
 P_c-v\ge A_{env}^{-1},\quad
 2(P_c-v)^2-(v-2)J_c^2\ge A_{env}^{-1}\quad(v\ge2).        \tag{1}
\]

The last two conditions apply only on the indicated subset. Demanding
v>2 of the original profile on its v<1 portions would be the wrong gate:
producing that growth is the purpose of modulation.

Set Xa=4/Lambda, tc=delta/(10sqrt(log C)), Xminus=Xa exp(tc/2),
Xi=110, Xsep=110exp(2000B), rhop=XR exp(241T-23), Xrep=rhop exp(5).
The table has exact shared endpoints; offsets at Xa are not added to a
rounded log Xa.

| Original input region | Endpoints | Established value margins |
|---|---|---|
| Activation | Xminus to Xa exp(delta) | Pc>2.2; a>C^-22; third/fourth gaps >C^-401, >C^-802 wherever v>=2 |
| Reduced continuation and final tapers | Xa exp(delta) to Xi | Pc>2.2; a>C^-22; v<=.8 |
| Shape | Xi to Xsep | Pc>2.5; a>=.7; v<=.9 |
| Ideal prejoin with axial offset | Xsep to XR exp(-8) | Pc>2.5; a=v=.8 |
| Five-moment joining annulus | XR exp(-8) to XR exp(-5) | Pc>3.32; a>.7; v<1 |
| Ideal after join | XR exp(-5) to XR | Pc>1.1X; a=v=.8 |
| First outer ramp | XR to e XR | Pc>.5XR; a=v in [.8,2] |
| Axial transition | e XR to XR exp(T+1) | a=2; Pc-2 and third/fourth gaps >1 |
| Intermediate ramp | XR exp(T+1) to XR exp(T+2) | a>=2; Pc>12; third/fourth gaps >1 |
| Intermediate power | XR exp(T+2) to Xrep | Same intermediate margins |

The first five rows follow from the [activation](../compact_pressure_bridge/README.md),
[continuation](../reduced_continuation_margins/README.md), and
[actual shape/join](../shape_join_pressure/README.md) audits. In activation
a=k*ar with k>=C^-20 and ar>=2/(5Lambda), Lambda<C, hence a>C^-22.
Its flat cutoff at Xminus satisfies ea>=C^-401. These are uniform estimates,
including the very small angles, not interpolation from a grid.

The sections below derive the remaining rows with finite-radius terms
retained. Thus a>C^-22 and Pc-2>.1 throughout. The required third/fourth
gaps exceed C^-1024 wherever v>=2; on v<1 rows that condition is vacuous.
Since C^-1024>Aenv^-1, this proves (1) at the SAME selected parameters and
radius as the compact jet envelope.

## 2. Correct the first outer-ramp comparison

After the exact join the ideal fields are E=Pf(X/XR)^.1, U=4eta,
f=(1+eta^2)^-1 and all five ideal primitives are restored. The exact
formula has Q0>1.1, so p1=XQ0/L>1.1X. The ideal segment after the join
therefore has p1>1.1*10000 exp(-5)>74.

On the next unit ramp use y=log(X/XR), 0<=y<=1. U=4eta,
ell=Dy log(sqrt(2X)E) is in [0,.6], a=2-2ell and b=0. Exact M matching
gives W=-3+(2+8h)eta^2<0. Also Hc=(D+4d)eta and
(log E)_eta=-2eta/(1+eta^2), so the original source has Sq>=-h.
The p1 equation and initial condition are

\[
 p_{1y}+\ell p_1=XS_q/L,\qquad p_1(0)>1.1X_R.             \tag{2}
\]

The kernel K(y,z)=exp(-integral_z^y ell) satisfies
exp(-.6(y-z))<=K<=1. For the POSITIVE initial value use the lower kernel;
for the NEGATIVE source lower bound use the upper kernel. Consequently

\[
 p_1(y)\ge1.1X_Re^{-.6y}-\frac h{.98}X_R(e^y-1)
 \ge X_R\{1.1e^{-.6}-(.01/.98)(e-1)\}>.586X_R>.5X_R.       \tag{3}
\]

The older [joined-stress note](../joined_stress_construction/README.md),
section 2, instead displayed e^(-.6) multiplying BOTH terms. That stronger
expression does not follow from Sq>=-h alone. We preserve its frozen bytes
and supersede that comparison here; the floor .5XR needed downstream is
unchanged. No actual source equation is changed.

A negative control solves the comparison model ell=.6, Sq=-h, L=.98,
p1(0)/XR=1.1, h=.01 exactly. At y=1 its normalized value is
.589856892864..., BELOW the old displayed .594070217442..., but ABOVE
the corrected .586159311657.... This refutes the stated kernel inference,
not the actual ramp (which has additional nonnegative source terms), or
the selected tiny h. Independent quadrature checks the negative-source sign.

## 3. Finite-radius axial and intermediate gaps

The [Md64 axial](../axial_stress_audit/README.md) and
[intermediate](../intermediate_decay_audit/README.md) bounds apply at
XR>=100; increasing to the SAME XR>10000 preserves them. On the axial
stage a=2, p=p1>35, |b*w|<.525247, b^2<=2.048e-25, w=p2/p1.
Put D0=b^2/2, A0=1-(2+D0)/p. Exactly,

\[
 P_c=p(1-bw/2),\quad v=2+D_0,\quad
 2(P_c-v)^2-(v-2)J_c^2
   =2p^2(A_0+D_0/2)(A_0-bw-D_0/2).                        \tag{4}
\]

The first factor exceeds .94 and the second .417. In lower-bounding the
first, DROP its positive D0/2 term: inserting an upper bound there would
give the wrong direction. Pc-2, the third and fourth gaps are all >1,
including b=0. At b=0, v=2 uses the relaxed condition.

On both intermediate pieces b=0, a=2+2alpha<=2.2 and p>12. The
whole-angle estimate alpha*w^2<1e-42, including its narrow-angle memory,
gives

\[
 P_c-v=p-a>9.8,\qquad
 2(P_c-v)^2-(v-2)J_c^2
 =2p^2\{(1-a/p)^2-\alpha w^2\}>2(12)^2(.3135)>1.           \tag{5}
\]

The reduced ratio G alone is not the finite-p fourth gap. That gap divided
by p^2 is G-4cv/p+2v^2/p^2, where c=1-bw/a and
G=2c^2-(v-2)(w+b/a)^2. All finite terms remain in (4)--(5).
A passing Md4 pulse fixture does not repair its failed axial stage.

## 4. Protected boundary growth and selected finite frequency

On 0<=y=log(X/Xa)<=tc the native reference is still the exact core, whose
exit growth exceeds 2.3. Use R=exp(B/100): ar, ar^-1, qbar=qr/C and
their required y jets are bounded by R. Differentiating
vr=ar+C^2*qbar^2/ar gives

\[
 |v_{ry}|\le R+C^2(2R^3+R^5)<4C^2R^5<K_{act}=C^2R^{80}<C^3. \tag{6}
\]

Native variation costs at most delta*Kact<C^-17. Since
y/delta<=1/(10sqrt(log C)), the flat step has
ea<=exp(4)C^-100<C^-99, while vr<=C^3. The exact activation identity
vs=k*vr+ea*delta*Vc, |Vc|<=Kact, yields

\[
 v_s>2.3-2C^{-17}-C^{-96}>2.29.                            \tag{7}
\]

The scalar comparison holds already for C>=1e9. On the right protected
power collar v_s-2=2lambda>C^-1. The loop cutoff scale
delta_loop=exp(-Aenv^20) is less than Aenv^-1<lambda, so its variance
cutoff is zero in both protected collars. Its zero-mean primitives then
vanish there EXACTLY. The modulation is supported between Xminus and
Xplus=XR exp(T+3); its separation from rhop is 240T-26>0.

Apply the [compact jet envelope](../compact_jet_envelope/README.md), with
cone-input hypotheses now supplied by (1) and (7), and the
[loop](../loop_modulation_bounds/README.md) and
[exact repair](../repair_state_bounds/README.md) lemmas. Corrected (3)
supplies the same .5XR floor used by that envelope; the older stronger
expression is NOT needed for its derivative estimates. The constants stay

\[
 G=e^{A_{env}^{32}},\quad H_{hier}=e^{G^4},\quad
 N=1+\lfloor H_{hier}^{32}\rfloor,\quad
 \epsilon=G^{-64},\quad\kappa=G^{-4},\quad R_{box}=4G^3.    \tag{8}
\]

The loop's four coordinate gaps are >=kappa. The COMPLETE modulation plus
repair error in (a,b,p1,p2) is <epsilon/2, retaining inverse-lambda moment
normalization and both pressure terms. The cone gradient bound
256 Rbox^10 and G^30>=512*4^10 leave, throughout the compact region,

\[
 a\ge G^{-3}/2,\quad v_s-2\ge\kappa/2,\quad
 P_c-v_s\ge\kappa/2,\quad
 2(P_c-v_s)^2-(v_s-2)J_c^2\ge\kappa/2.                     \tag{9}
\]

Thus the SAME finite integer N is accepted analytically. It is not
numerically resolved, and its large radial derivatives are not discarded.
This discharges the regional value/cone inputs, not every later actual
physical derivative/support obligation.

## 5. Join to the compensated exterior with complete memory

The endpoint identity is exact:

\[
 X_{rep}=\rho_p e^5=X_Re^{241T-18}=X_*,                    \tag{10}
\]

the left scale of the HEAT-compensation patch. Its bumps are strictly
inside the following five-unit patch. Before Xstar that future edit does
not alter forward primitives; its total pressure change with the heat
edit is zero, so the same Pi0 is used. At Xrep the compact repair restores
ALL FIVE cumulative moment functions and the fields. The
[global moment/stress ledger](../global_moment_stress/README.md) propagates
exact equality of forward pressure/stress at every later radius. Field
value matching alone would not do this.

Apply the [heat-compensated exterior audit](../heat_exterior_audit/README.md)
from Xstar through terminal y=1/2. It includes the intermediate power,
pulse and angular corrections, with nonzero compensation memory between
separated supports. Its actual-profile bounds are Q>=h/4000, a-2>=h,
first/second reduced margins >1.47, >.81, Pc>vs, and normalized finite-p
fourth margin >.80. Do NOT apply an unedited pulse estimate alone there.
On terminal 1/2<=y<3 the [outer-edge formulas](../outer_edge_stability/README.md)
give positive angular stress and squared cone ratio <=16h^9<1, with
direction tending to (1,0). Exact global moments identify these backward
stresses with the assembled forward stress. The first axial-viscous
correction still has a SEPARATE positive-order common-support premise.

On Xa<X<Xminus, the unchanged inner flat-quotient construction and (7)
give the admissible cone, with smooth limiting direction
(ar,-br)/sqrt(ar^2+br^2), opposite the physical shear. Equation (9)
covers the middle; the compensated exterior
covers the rest. The leading stress is nonzero throughout Xa<X<Xb and
zero in the exact core and exterior. The known inner/outer direction
limits and compactness join their positive directional margins into some
fixed profile-dependent positive margin. This is a value/direction result,
not an audited all-order weighted-edge or physical wave certificate.

## 6. Evidence and next obligation

The new record separates exact rational coordinate identities, outward
scalar inequalities at 100/160 digits, and a labelled comparison model.
Its ten-region endpoint registry checks coverage, not sampled suprema.
The negative-source model rejects the old kernel inference; finite-p and
v<2 controls reject treating reduced ratios as the full admissible cone.
Historical files, including the stronger comparison, remain frozen; this
note is their correction and the new status record.

Local validation on 2026-10-16: 334 tests passed across 14 suites (this
checkpoint, global moments, compact envelope, loop, repair, activation,
continuation, shape/join, core growth, inner collar, axial, intermediate,
heat exterior and outer edge). All 14 corresponding frozen evidence
records replayed successfully. This checkpoint adds 31 tests and checks
533 inherited provenance hashes. The parent commit's CI succeeded; the
new commit's CI result is reported separately, not presumed.

Next: audit positive-order coefficient support and the genuinely untouched
reserved-mean intervals. The leading repair and heat compensation consume
two moment patches; those patches must not simply be declared available
again for later corrections. Then audit all-scale actual-pulse size/derivative/support
and complete-wave-covariance mean corrections including axial velocity and
pressure. The complete residual still needs viscosity, cutoffs, every
interaction and a justified improvement estimate. Neither passing tests
nor (9) proves those obligations or blow-up.

```sh
python3 research/navier_stokes_cascade/results/regional_cone_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/regional_cone_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/regional_cone_join/evidence.json \
  --output regional_cone_join_recheck.json
```
