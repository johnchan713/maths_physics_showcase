# Explicit leading stress floor on the retained compact interval

Parent: `5136b026a5626b21a283fba4f74a9ee424d66060`.
The [first-seminorm estimate](../first_seminorm_bound/README.md) left a
quantitative positive lower bound tau for |T0| on [Xinner,Xc] as an input.
This checkpoint supplies that input from the SAME leading core, modulation,
exact repair and compensated exterior estimates:

\[
 \boxed{\ |T_0(X,\eta)|\ge\tau_*:=\frac{f_*}{64G^{10}}>0,
 \quad X_{inner}\le X\le X_c,\quad -1\le\eta\le1.\ }       \tag{1}
\]

Here Xinner=Xa exp(tc/16), Xc=Xtail exp(1/2)=Xb exp(-5/2) is the
terminal heat-coordinate y_tail=log(X/Xtail)=1/2,
Aenv=C^4096, G=exp(Aenv^32), and f* is the explicit GLOBAL leading swirl
floor from [leading-frame foundation](../leading_frame_foundation/README.md).
This is an analytical consequence of the retained exact leading/core and
regional lemmas, not a numerical sample, a new independent core theorem,
or a proof-assistant certificate. No parameters, roots, frequency, width,
gauge or supports are changed; historical evidence is preserved.

The new work is the quantitative THREE-REGION amplitude join and its
substitution into the first-coefficient estimate. The local core growth,
compact input/cone, exact moment matching, heat positivity and canonical
seed results were established earlier and are not recounted as new proofs.
In particular (1) is NOT a floor at either flat stress edge and is NOT
a same-tilt cone margin or a bound for the full corrected/wave stress.

## 1. Use the existing global swirl floor, with no circular coefficient input

Let Yb=log(Xb/Xa)>0 and Ab=4Aenv. The retained leading estimate is

\[
 f_* = \frac{\exp(-2\Lambda M)}{4C}\exp(-A_bY_b/2)
 \le F(X,\eta)\qquad (X_a\le X\le X_b).                    \tag{2}
\]

It follows by integrating -2 partial_logX log F=a, with 0<a<=Ab and the
exact core anchor F(Xa)>=exp(-2Lambda M)/(4C). Positive exact repair
neighborhoods and the heat multiplier preserve this SAME leading F.
No positive-order coefficient or small-q acceptance bound occurs in (2),
so using it in the coefficient bound does not feed that bound back into
its own leading hypothesis. Finite does not mean representable in binary64.

## 2. Protected inner segment: keep the actual offset and flat exponent

The first segment is

\[
 [X_{inner},X_-]=[X_a e^{t_c/16},X_a e^{t_c/2}],\quad
 t_c=\frac{t_1}{10\sqrt{\log C}},\quad t_1=C^{-20}.          \tag{3}
\]

This lies inside the protected collar y=log(X/Xa)<=tc where modulation is
exactly absent. The retained [core/activation bound](../core_growth_activation/README.md),
section 5, is

\[
 |T_0|\ge\frac{9f_a}{16C^2}\exp(-1/s^2),\quad
 s=y/t_1,\quad f_a=\exp(-2\Lambda M)/(4C).                 \tag{4}
\]

On (3), s>=1/(160sqrt(log C)). Since f_a>=f*,

\[
 |T_0|\ge(9/16)f_*C^{-25602}.                             \tag{5}
\]

The number 25600=160^2 comes from the ORIGINAL tc/16 offset. It is not
an adjustable flat-exponential constant or an asymptotic O(1). The factor
C^-2 in (4) remains. Taking the infimum down to Xa would instead give zero.

For C>=2, log C<=C-1 and C^13>=4096C imply
log G=C^131072>4096log C, hence G>C^4096. Therefore

\[
 \frac{(9/16)f_*C^{-25602}}{f_*/(64G^{10})}
 =36G^{10}C^{-25602}>36C^{15358}>1.                       \tag{6}
\]

Thus the same floor (1) holds on the first segment, without moving its
left endpoint or asserting a positive unweighted floor at the flat edge.

## 3. Compact modulation plus exact repair: derive the amplitude from D

The second segment is [Xminus,Xrep]. Use the existing compact variables

\[
 t=-b/a,\quad v=a+b^2/a,\quad
 D=p_1+tp_2-v,\quad J=p_2-tp_1.
\]

The [regional value join](../regional_cone_join/README.md), section 4,
gives for the SAME accepted finite N and exact repair root

\[
 a\ge G^{-3}/2\ge R^{-1},\quad
 |(a,b,p_1,p_2)|_\infty\le R,\quad D\ge\kappa/2,
 \quad R=4G^3,\quad\kappa=G^{-4}.                        \tag{7}
\]

The complete modulation plus repair error, not just field closeness,
is already included in (7). With the construction's normalized viscosity
one, the exact normalized leading stress (physical scale q=1) is
T0=F(p1-a,p2+b). Rotating the two coordinates gives

\[
 |T_0|^2=F^2\frac{D^2+J^2}{1+t^2}.                       \tag{8}
\]

Since |t|<=R^2 and sqrt(1+t^2)<=2R^2 for R>=1, (2),(7),(8) imply

\[
 |T_0|\ge\frac{f_*\kappa}{4R^2}
             =\frac{f_*}{64G^{10}}.                     \tag{9}
\]

Only the positive D margin is needed for this amplitude inequality. The
other growth and fourth-cone gaps remain necessary for their own purposes;
(9) does not replace them. In particular a positive swirl and growing
shear alone would allow the local target (p1,p2)=(a,-b), giving zero stress.

## 4. Same-profile compensated exterior, including its finite-radius terms

The final segment is [Xrep,Xc]. Exact geometry identifies
Xrep=Xstar=XR exp(241T-18), the compensation patch's left endpoint.
The [compensated exterior](../heat_exterior_audit/README.md), sections 2--3,
retains on this WHOLE segment, including the interior of its patches,

\[
 Q\ge h/4000,\quad hX_*/4000>10^{12},\quad
 a_{ref}\le4,\quad |a-a_{ref}|\le h\epsilon_Q<1.          \tag{10}
\]

Thus a<5, and p1=XQ/L>=hXstar/4000>10^12, since 0<L=1-2h eta^2<=1.
Consequently the angular stress alone gives

\[
 |T_0|\ge T_{0\theta}=F(p_1-a)>F\ge f_*\ge\tau_* .       \tag{11}
\]

This is an application of the already proved reference/edit bounds, not a
new heat integral or sign proof. It applies to the assembled leading
stress ONLY because the [current same-root pressure/moment join](../global_pressure_moment_foundation/README.md)
restores all five FUNCTIONAL moments in the original gauge and identifies
the forward and backward stresses. No pressure constant, eta derivative,
or memory between separated compensation supports is reset. Local field
matching or approximate numerical moment cancellation would not justify
transferring (10)--(11) to that stress.

The intervals (3), [Xminus,Xrep], [Xrep,Xc] have exact common endpoints and
cover [Xinner,Xc]. Endpoints eta=+-1 use the retained smooth one-sided jets.
There is no negative-time heat extension and no uncovered transition.
Equations (5),(9),(11) prove (1), relative to the stated leading lemmas.

## 5. Substitute into the first-seminorm gate, without promoting its scope

In the previous notation Xc=c^2 and Htheta,Hz bound the complete first
tangential residual coefficients, including pressure and shifted viscosity.
Its conservative Euclidean sum bound now has NO free tau input:

\[
 A_W\le\frac{64G^{10}}{f_*}
          \left(\frac{X_cH_\theta}{2}+\sqrt{X_c}H_z\right). \tag{12}
\]

All global w0=Omega0/X moment integrals, reciprocal-amplitude angular jets,
cutoff and bump derivatives in that estimate remain. Pi1/Pi2 pressure
memory in the summed background and every full-wave/covariance/residual
obligation remain unchanged. The exact inner coefficient equations still
give T1=0 before Xinner, so (12) bounds the weighted quotient there by zero,
not by dividing by T0(Xa)=0.

Insert (12) into the preceding C* formula. The sufficient whole-4Q-chart
test remains

\[
                  2h(\ell-2)\ge\log_2(C_*/r_{K_c}).       \tag{13}
\]

What changed is that tau is now the explicit (1), not a compactness-only
constant. The actual collar/leading jet and pressure-source envelopes,
bump/cutoff bounds and the quantitative SAME EXISTING TILT margin in r_Kc
still need to be supplied/checked together with all other accepted-tail
gates. An amplitude floor does not choose a new tilt or prove those gaps.
No actual scale test or historical first-band acceptance is asserted here.

For stable bookkeeping keep

\[
 \log\tau_*=-2\Lambda M-\log(4C)-A_bY_b/2-\log64
                                  -10C^{131072}.         \tag{14}
\]

This is finite for the frozen parameters. Floating underflow of exp(logtau)
is not an analytical zero. The floor is profile-dependent, not uniform as
C,h,N, the width or the compact left endpoint change. It is a leading
stress floor, NOT an unconditional actual-pulse or full residual estimate.

## Reproduction and next obligation

The new checks cover the exact offset/exponent ledger, 81 rational rotated
stress identities and amplitude controls, 48 compensated-exterior algebra
controls and 64 exact substitutions into the first-coefficient bound.
These controls are not selected-profile samples. The continuum regional
inequalities and their dependencies are given above, not inferred from a
finite grid. Frozen dependency evidence is replayed under its own original
comparators; no historical evidence or tolerance is edited.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
Next quantify the actual activated-collar operator/jet envelopes and the
same-tilt leading cone radius, then verify (13) with the other value gates.
Unconditional pulses, complete residual closure and blow-up remain open.
