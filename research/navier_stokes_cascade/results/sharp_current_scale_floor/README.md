# Sharp complete cutoff cost and the mixed current-scale obstruction

Parent: `bbfe73f645de9184bc4cd4e52bd8cb8d64a099c2`.
Status: `sharp-full-cutoff-cost-floor-and-coupled-current-scale-certificate`.

This checkpoint solves the full four-row cutoff-cost minimization, including
competition between the stress and first-derivative rows. It then eliminates
cutoff shape from the MIXED current-frame scale test and derives its sharp
best-case required-deficit floor. Equality is sometimes attainable when a
constant row dominates; it is otherwise excluded by smooth endpoint flatness.
The earlier stress-row infimum explicitly left this full-max problem open.

These are exact results for the specified sufficient-bound program. They
do not evaluate the unchanged cutoff or the selected enormous profile, and
do not establish a lower bound on the actual PDE error. No cutoff, field,
tilt, scale, root, carrier, support, pressure gauge or memory is changed.
All historical frozen evidence is retained. The elementary optimization
uses the proved smooth-capacity theorem; applying it to the Navier-Stokes
ledger still requires every retained leading/core/regional premise.

## 1. Optimize the complete affine bound, holding its inputs fixed

Keep the corrected logarithmic collar width and the complete pressure,
repair and viscosity costs. The bound is

\[
 C_\star(M,N)=\max(C_0,D_0+D_1M,A_0+A_1M+A_2N),\quad
 C_0=\max(1,R_v),\ D_1>0,\ A_2>0,\ A_1\ge0.              \tag{1}
\]

All coefficients here are FIXED. The retained first-coefficient ledger
bounds repair moments using the cutoff value bound 0<=K<=1, so these
envelope coefficients do not depend on M,N. Reoptimizing shape-dependent
exact moments or changing another input would be a different program.
In particular full future pressure and shifted axial viscosities remain in
A0; the positive A2 is the viscous second-cutoff-derivative cost.

For the normalized smooth shape K, equal to 1 before s=0 and 0 after s=1,
the earlier continuum capacity theorem gives feasible norm enclosures iff
`A(M,N)=integral min(M,Ns,N(1-s)) ds>1`. Hence M>1 and

\[
 N> M^2/(M-1)\quad(1<M\le2),\qquad N>4\quad(M\ge2).       \tag{2}
\]

For an EXISTENTIAL minimization, M>2 can be replaced by M=2 without
increasing any row, still with N>4. This can describe a different auxiliary
shape and must not be used to assign new bounds to the retained K.
Set t=M-1 in (0,1]. In the closure N=2+t+1/t define

\[
 P=A_1+A_2,\quad Q=A_2,\quad S_0=A_0+A_1+2A_2,\quad G_0=D_0+D_1,
\]
\[
 f(t)=S_0+Pt+Q/t,\quad g(t)=G_0+D_1t,\quad
 t_s=\sqrt{Q/P}\le1.                                    \tag{3}
\]

For every smooth shape the stress row is STRICTLY above f(t).
After ts, both f and g are increasing, so the closed-variable minimizer
lies in (0,ts]. On that interval f decreases and g increases. This gives
the complete solution:

* If g(ts)<=f(ts), set t*=ts and V=S0+2sqrt(PQ).
* If g(ts)>f(ts), t* is the unique crossing in (0,ts), and
  V=G0+D1*t*=f(t*).

The exact branch comparison is the sign of

\[
 G_0-S_0+(D_1-2P)\sqrt{Q/P}.                             \tag{4}
\]

For the crossing put a=P-D1,b=G0-S0,Delta=b^2-4aQ. If a=0,
t*=Q/b (b>0 in this branch). Otherwise

\[
 t_*={b-\sqrt\Delta\over2a}
     ={2Q\over b+\sqrt\Delta}.                          \tag{5}
\]

Here Delta>0 and the denominator in the second formula is positive.
Indeed f-g is strictly decreasing on (0,ts], starts at +infinity and
ends negative; its first root has nonzero derivative. The first formula
is useful for exact quadratic-field arithmetic, not rounded subtraction.

Thus the full infimum is

\[
 \boxed{I=\inf C_\star=\max(C_0,V).}                     \tag{6}
\]

Sharpness follows from (2) by taking N just above its frontier and using
the previous interior shrink/mollify/normalize construction. This proves
the continuum statement; finite optimizer controls do not supply its proof.

## 2. Attainment matters at the boundary

The variable maximum is always strictly larger than V for a smooth shape.
In the stationary case f(t)>=V and smoothness makes the stress row larger.
In the crossing case t<t* gives f(t)>V, t>t* gives g(t)>V, and at t*
the actual stress row is larger than f(t*). The M>2 case follows from the
same comparison at t=1.

Consequently, if C0<=V, every smooth shape has Cstar>I and the infimum
is unattained. If C0>V, choose t=t* and N slightly above its frontier so
both variable rows remain below C0. The smooth-capacity construction then
gives Cstar=C0=I. These results concern existence of SOME smooth shape.
They do not identify the original cutoff.

A complete cutoff certificate Cstar<=B therefore exists iff

\[
 \boxed{B>I\quad\hbox{or}\quad[B=I\ \hbox{and}\ C_0>V].}   \tag{7}
\]

This distinguishes the full maximum from the earlier always-unattained
stress-row infimum. No historical statement is changed; that distinction
was explicitly reserved there. For example, Rd=3M,Rw=N has stress infimum
4 but full infimum 9/2, unattained. Conversely Rd=M,Rw=3M+N with constant
row 10 has an attained full minimum 10 although its stress infimum is 9.
If the constant is exactly 9, equality remains impossible.

The implementation represents a+b sqrt(q) by exact rational coefficients.
Opposite-sign terms are compared only after their signs are known, by
comparing a^2 and b^2 q; it never squares an unguarded inequality or uses
a floating optimizer. Decisions are independently checked against the
older capped-quadratic feasibility algorithm and its smooth witnesses.

## 3. Eliminate shape from the CURRENT mixed freezing gate

Use the preceding current-frame constants Cm=Cmesh,Cp=Cpert>0, with
normalization 8 and orientation omega/3 already charged. Suppose the
uniform error enclosure is E<=Cstar/S, with S>0. For the exact physical
power S=2^(2h(ell-2)); rational diagnostics may use the smaller
2^floor(2h(ell-2)), giving a conservative certificate.

For a valid same-u raw deficit lower bound d0, let
gammaB=omega^2*d0/(9kappa_u). The previous gates are E<=rb/8,
E<=omega^2*d0/(16Lpsi), and

\[
 C_m\ell^{-6}+8C_pE\le\gamma_B/8.
\]

Solve the LAST inequality with its value-error term retained. Define

\[
 H=\gamma_B/8-C_m\ell^{-6},\qquad
 B_{eff}=S\min\left\{r_b/8,{\omega^2d_0\over16L_\Psi},
                         {H\over8C_p}\right\}.           \tag{8}
\]

Because Cstar>0 and Cp>0, H<=0 excludes all shapes from this certificate,
even if the leading/mesh-only comparison would accept equality. For H>0,
combine (7) at Beff with the INDEPENDENT positive-chart test
ell^6>=Dbox/dchart. This is an exact existence test for the supplied
error-bound program. If a rounded S is used, it is sufficient, not an
equivalence to the unrounded physical test. Neither version is an iff
criterion for the actual PDE cone or the original first band.

When accepted, the existing capped-quadratic algorithm provides rational
M,N with capacity>1 and every original row<=Beff. The resulting full
current-frame bound is checked directly, retaining the complete pressure
target at the evaluation point and the compact/outer sign cover. No
global outer stress-magnitude condition is introduced.

An auxiliary control has Cm=0 and chooses S so Beff=9 for Rd=M,Rw=3M+N.
The background-only shape and independent mesh tests pass, but the mixed
certificate fails because its limiting cost 9 is unattained. Increasing
the budget to 10 passes. A dominating constant row at 10 also passes
equality, in agreement with the attainment distinction above.

## 4. A shape-independent required-deficit floor

Substituting the full cost infimum I into the current-frame formulas gives

\[
 \boxed{d_{inf}(\ell,S)=\max\left\{
 {16L_\Psi I\over\omega^2S},
 {72\kappa_u C_m\over\omega^2\ell^6}
       +{576\kappa_u C_p I\over\omega^2S}\right\}.}       \tag{9}
\]

Every smooth cutoff requires d_req>=d_inf in THIS sufficient certificate.
The inequality is strict if C0<=V, because both affine functions of Cstar
have positive slopes. If C0>V the floor is attainable by some shape,
subject separately to base-error and chart requirements. Similarly the
base profile error has infimum 8I/S with the same attainment distinction.

Thus a supplied d0 below d_inf cannot pass for any smooth shape; equality
also fails in the unattained case. This can be assessed before estimating
the old cutoff. It is not a physical error lower bound: the affine rows
are upper estimates and may be conservative. Failure excludes this
certificate, not a Navier-Stokes solution or cone. Passing (7)-(8) still
does not show that the unchanged cutoff's own norms lie in the feasible
region, or that the selected first scale meets every other obligation.

## 5. Evidence, provenance and next obligation

Ten exact optimizer regimes include stationary and both quadratic-crossing
signs, the linear crossing, irrational roots, constant dominance, and tied
boundaries. There are 414 exact budget comparisons with the independent
older capped-quadratic decision. A complete global-jet/pressure/viscosity
composition uses the retained full affine rows, generates only an auxiliary
smooth-shape budget witness, and checks it in the original mixed gate.
The full-pressure negative control changes the constant stress cost;
pressure or shifted viscosities are never erased in the accepted program.

All are algebra controls or manufactured parameter examples, not actual
enormous-profile suprema or a numerical continuum proof. Frozen source
bindings, original comparators and historical evidence are preserved.
No new smooth cutoff is installed and no distinct field/force is relabelled
as the historical one.

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output /fresh/path.json`.
Next obtain the unchanged cutoff's valid M,N, the actual raw-deficit lower
bound and all endpoint constants; apply the full mixed gate, not just its
optimal existential relaxation. Geometry, length/carrier/inverse gaps,
uniform actual pulses on every physical scale, full covariance/axial-pressure
mean corrections and the complete viscous cutoff/interaction residual
remain open. Unconditional pulse bounds, residual improvement and blow-up
are not established.
