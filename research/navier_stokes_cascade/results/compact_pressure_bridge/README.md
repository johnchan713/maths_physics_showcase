# Whole activation and the compact pressure-memory bridge

Status: `whole-activation-and-conditional-compact-pulse-bridge-audited`.

Parent: `08c4c30d8a9f8de4ebd00b7484db6ab55530fe24`.
No parameters, historical evidence, or earlier derivations are changed.
This extends the **same** activation estimate from its first collar to its
whole first transition. It also derives quantitative pulse value bounds
from the conditional compact cone certificate and spells out exactly what
the five-moment repair must restore. It does not independently establish
the relaxed-cone margins in the subsequent reduced-shear continuation and
final tapers. That input prevents promotion to a whole-profile certificate.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
(4.16), (B.26)--(B.29), and (C.11)--(C.17), SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Printed page 161 was checked visually. The theorem in the manuscript is
not a premise. The calculations below are analytical derivations with
regression checks, not a proof-assistant certificate or independent peer review.

## 1. The activation estimate works on the whole first transition

Use all definitions of [core growth and activation](../core_growth_activation/README.md).
In this section t=t1=C^-20, y=t*s, 0<=s<=1,
ea=(1-kappa0)*sigma(s), k=1-ea, and kappa0=C^-20.
The reference remains the exact native core here: 4<=Y<=4*exp(t)<4.1.
The basic reference-jet envelope R=exp(B/100), and
9/4<=vr<=C^3, hold on this entire much shorter interval.

The previous cubic flat quotient is needed for smooth edge directions;
it need not be used away from the edge. Monotonicity gives the simpler
positive operator bound, including angular C1 norms,

\[
 \left\|\int_0^y e_a(z)b(z,\eta)\,dz\right\|_{C^1_\eta}
 \le t e_a(y)\sup_{0\le z\le y}\|b(z,\cdot)\|_{C^1_\eta}.       \tag{1}
\]

The cutoff has no eta dependence. Applying (1) twice bounds a second
integral by t^2*ea times the input norm. At s=0 interpret the resulting
quotients by the earlier smooth flat-quotient extension. No bound on
physical radial derivatives is inferred from (1).

Subtract the reference equations exactly:
Delta log E_y=ea*ar/2, Delta U_y=ea*Er*qr/2.
Consequently the complete table in section 4 of the parent remains valid
after replacing **t*ea*s^3 by t*ea**, and **t^2*ea*s^6 by t^2*ea**.
Here is why the extension does not introduce a new unbounded coefficient:
the old normalized integral operator cost at most 3/2; (1) costs at most
one. All other steps are products, exprel, and positive reciprocal
differences on the same reference collar, with t*R^3<.1. The factor
t*ea*s^3 occurring inside nonlinear products is replaced by t*ea<=t.
Thus each of the parent's R powers is still an upper bound. The moments
are integrated with dx=t*X*ds; the Cp row instead uses Cp_y=E^2/2.
Both 4*Aphys*eta*Delta Cp and -d*Delta Cp_eta remain in Ns.

With the parent's Kact=C^2*R^80<C^3, the exact algebra therefore gives

\[
 (P_c-v_s)/e_a=v_r+tD_c,\quad J_c/e_a=tJ_c^*,\quad
 (v_s-kv_r)/e_a=tV_c,\qquad |D_c|,|J_c^*|,|V_c|\le K_{act}.  \tag{2}
\]

These are bounds on values and the required angular jets. The first
collar retains the stronger cubic vanishing previously proved. Since
t*Kact<C^-17, for C>=2 equations (2) imply on every 0<s<=1

\[
 P_c\ge9/4-2C^{-17}>11/5,\qquad
 P_c-v_s\ge e_a(9/4-C^{-17})>(11/5)e_a,                    \tag{3}
\]
\[
 { (v_s-2)_+J_c^2\over2(P_c-v_s)^2}
 \le {25\over121}C^{-31}<C^{-30}<1/2.                     \tag{4}
\]

For (4), use (vs-2)+<=2C^3. Each comparison becomes weaker at C=2:
2C^-17<1/20 and 25/121<C. Thus these are universal scalar inequalities,
not conclusions from the finite C fixtures. Also a=k*ar>0.
The strict relaxed cone follows on the whole transition, and the full
admissible cone wherever vs>2. At s=1, vs<=2C^-17<1: requiring the full
cone of the unmodulated transition everywhere would be the wrong gate.

At Xminus=Xa*exp(tc/2), ea>=C^-401. From there through activation,
wherever vs>=2, the third and fourth gaps exceed C^-401 and C^-802.
Indeed the fourth is 2(Pc-vs)^2*(1-loss)>4ea^2.
This independently supplies this activation portion of the compact input
margin, beyond the smaller first collar of the parent checkpoint.

## 2. Exact pressure memory between separated supports

Let two profiles share the same axis pressure datum and have identical
fields on a radial interval after an edit. Write their moment differences
as mu=(m,i,j,s,c), functions of eta. They are constant in X there because
the five original moment densities agree. They need not be zero.
Let Aphys=.5+h, D=.5-h, d=1-eta^2, L=1-2h*eta^2, H=sqrt(2X)*E.
Subtracting the original forward formulas, with a prime meaning eta
derivative, gives the exact identities

\[
 \Delta\Pi=c,\qquad
 \Delta W=-(2D\eta m+dm')/X,
\]
\[
 B_\mu=(1-h)i-D\eta i'-dj'+2(h-D)\eta j,
 \quad\Delta Q_s=-\Delta W+B_\mu/(XH),
\]
\[
 \Delta N_s=-U\Delta W+
 [D(m-\eta m')+4h\eta s-ds']/X+4A_{phys}\eta c-dc',
\]
\[
 \Delta p_1=X\Delta Q_s/L,\qquad
 \Delta p_2=X\Delta N_s/(LE).                              \tag{5}
\]

The -U*Delta W term matters outside the U=0 repair patch. Equality of
fields alone does not erase memory. Even m=i=j=s=0 leaves
Delta p2=X*(4Aphys*eta*c-d*c')/(LE). For example c=epsilon*(1+eta)
gives Delta p2=-epsilon*X/E at eta=0. On the reserved power law
E=e(eta)*(X/rhop)^(-1/2-lambda), its magnitude grows like
X^(3/2+lambda). A finite small error is not exact restoration at all
later radii. This is an algebraic failure control, not a counterexample
to the manuscript, which explicitly repairs all five rows.

There is a sharper support obstruction. For pressure-only memory the
stress difference is exactly
`Delta T= sqrt(X/2)*(0,(4*Aphys*eta*c-d*c')/L)`.
At eta=0 with the same c, this is `(0,-epsilon*sqrt(X/2))`.
If the unedited stress vanishes at its finite outer edge, this error has
a nonzero limit there for every epsilon>0. An arbitrarily accurate
numerical repair is therefore not a proof of the claimed stress support
or flat edge. The exact contraction root, rather than a numerical Newton
stopping tolerance, is needed for that conclusion.

The exact five-row repair sets mu=0 **as functions of eta**, hence also
mu'=0. Equation (5) proves equality of all forward pressure/stress
coordinates at every subsequent radius with equal fields, including gaps,
remaining reserved patches, and the outer collar. Four rows or agreement
at isolated angles are insufficient. Conversely this propagation identity
does not itself prove the original profile's backward pressure normalization.

## 3. Quantitative pulse bounds from the conditional compact cone

This section is a conditional implication. Assume the original compact
input hypotheses used by [the modulation](../loop_modulation_bounds/README.md)
and [repair](../repair_state_bounds/README.md), including their relaxed-cone
margins. Keep A=C^4096, G=exp(A^32), Hhier=exp(G^4),
N=1+floor(Hhier^32), epsilon=G^-64, beta=1000, qstar=10^5,
Dtarget=Cstate=Hhier^16, Ccorr=A^128, r=A^-128.
The transformed target retains its lambda^-1 loss.

The exact contraction and state-error estimates give coefficient norm
<=2*beta*Dtarget/N<r/2 and combined coordinate error <epsilon/2.
The existing cone gradient bound 256*R^10 with R=4G^3 and
epsilon<=G^-4/(512R^10) leaves all four gaps at least kappa/2,
where kappa=G^-4. Also a>=(1/2)*G^-3>=1/R, |(a,b,p1,p2)|<=3A<R,
and 1/a<=R. Positivity of both field perturbations gives E>=C^-4/2.
Since X<C^12, F=E/sqrt(2X)>=C^-11=:fmin.
These constants are fixed-profile constants, not uniform in C or N.

Here are explicit new pulse-value consequences. Put t=-b/a,
v=a(1+t^2), Dcone=Pc-v, J=Jc, shear g=F*(-a,b),
T=F*((p1-a),(p2+b)). Then

\[
 |t|\le R^2,\quad v\le2R^3,\quad |J|\le2R^3,
 \quad 0<D_{cone}\le4R^3.
\]

The exact identities are

\[
 -\widehat T\cdot(g/|g|)={D_{cone}\over\sqrt{D_{cone}^2+J^2}},
 \quad \lambda_0^2=2aF^2{v-2\over v},
 \quad \mathcal Q^2={(v-2)J^2\over2D_{cone}^2}.             \tag{6}
\]

They yield the scale-independent, strictly positive bounds

\[
 F\ge f_{min},\quad |g|\ge f_{min}/R,\quad
 \lambda_0^2\ge\kappa f_{min}^2/(2R^4),
\]
\[
 |T|\ge f_{min}\kappa/(4R^2),\quad
 -\widehat T\cdot(g/|g|)\ge\kappa/(10R^3),\quad
 \mathcal Q^2\le1-\kappa/(64R^6)<1.                       \tag{7}
\]

For the last bound, use the fourth gap >=kappa/2, not just Dcone>0.
Store its deficit symbolically: numerical subtraction from one can erase
it. For |T|, use |T|=F*sqrt(Dcone^2+J^2)/sqrt(1+t^2) and
sqrt(1+t^2)<=2R^2. For orientation use sqrt(Dcone^2+J^2)<5R^3.
The factor v in the growth denominator in (6) is indispensable.
The mean-preserving shear loop and its moving inverse phase are the
inherited analytic construction, not numerically resolved actual fields.

## 4. Joining and the precise remaining inequality

The unchanged inner region and the compact region overlap on the known
first collar [Xa*exp(tc/2),Xa*exp(tc)]. Equation (3) now also covers the
rest of the first activation. Exact repair restores (5) to zero before
the outer collar; hence the existing outer-edge estimates apply there
without an additional pressure-memory error. Compactness together with
(7) would join their value/cone bounds on one fixed annulus.

This implication must **not** promote the old compact input certificate
to an independently audited foundation. Beyond the first activation,
on Jcont=[Xa*exp(t1),Xi] including both final tapers, the next unresolved
actual-profile test is the quantified relaxed-cone input

\[
 \inf_{J_{cont},\eta}(P_c-2)\ge A^{-1},\quad a\ge A^{-1},
\]
\[
 \inf_{J_{cont},\eta:\ v_s\ge2}(P_c-v_s)\ge A^{-1},\qquad
 \inf_{J_{cont},\eta:\ v_s\ge2}
 [2(P_c-v_s)^2-(v_s-2)J_c^2]\ge A^{-1}.                  \tag{8}
\]

An infimum over an empty v_s>=2 subset is vacuous. A sufficient route is
to establish the reference growth bound uniformly along the reduced-shear
continuation, transfer all forward moment errors at the fixed widths,
and separately control both final tapers. The local Y=4 exit gap alone
does not establish (8) on this longer interval. The older continuation,
shape/join and outer regional arguments remain dependencies to audit;
no numerical fixture here supplies their missing suprema.

Next: prove or refute (8) from the actual continuation equations, then
audit the completed global backward-pressure/moment normalization and
reserved supports. Uniform actual-pulse derivatives/supports, mean
corrections including axial velocity and pressure, complete viscous
residual improvement and blow-up remain unverified. No full q-uniform
pulse estimate follows from value bounds alone.

## 5. Reproduction

`bridge.py` checks exact rational memory and pulse-coordinate identities,
positive whole-interval integral operators, universal scalar endpoint
budgets, and deliberately dropped pressure/velocity terms. These are
manufactured controls; their seed fields are not the exact giant core.
`audit.py` checks frozen provenance, records scope, and replays the evidence.

Validation of this checkpoint: 19 new tests and 154 inherited tests pass
across the core, inner activation, compact envelope, loop, repair, and
outer-edge suites. All eight corresponding evidence replays pass; all
486 inherited source hashes match. These checks validate the recorded
identities, diagnostics, and provenance, not every continuum hypothesis.

```sh
python3 research/navier_stokes_cascade/results/compact_pressure_bridge/test_audit.py -v
python3 research/navier_stokes_cascade/results/compact_pressure_bridge/audit.py \
  --verify-record research/navier_stokes_cascade/results/compact_pressure_bridge/evidence.json \
  --output compact_pressure_bridge_recheck.json
```
