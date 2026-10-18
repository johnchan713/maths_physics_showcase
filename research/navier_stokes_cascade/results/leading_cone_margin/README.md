# Quantitative leading direction margins and the unchanged-tilt gate

Parent: `1f031a3c84cfafc8461bdb531d05adde904581b9`.
The preceding [stress-floor join](../regional_stress_floor/README.md) is an
amplitude estimate, not a cone margin. This checkpoint derives explicit
orientation and RAW cone-deficit bounds from the same regional leading
profile, including the limiting unit directions at both flat edges:

\[
 H_0=-N\cdot n_0\ge\omega_*:=\frac1{640G^{13}},\qquad
 \rho^2:=c_0^2\frac{(K\cdot n_0)^2}{H_0^2}
       \le1-\varepsilon_*,\quad
 \varepsilon_*:=\frac1{2^{18}G^{22}}.                     \tag{1}
\]

Here n0 is the normalized leading stress extended smoothly to Xa,Xb,
G=exp(Aenv^32), Aenv=C^4096, and G>C^4096 for C>=2. The result is an
analytical implication of the retained exact core, modulation/repair,
same-gauge global moment and compensated heat lemmas. It is NOT an
independent unconditional proof of those lemmas or of a physical pulse.
No root, width C^-20, support, frequency, pressure gauge or tilt is changed.
Historical evidence and its original comparators remain untouched.

An important remaining distinction is made explicit below: (1) supplies a
quantitative same-tilt margin only if an inequality involving the EXISTING
u is verified. An existential admissible-u choice does not automatically
satisfy a stronger, conservative numerical bound. Failure of that sufficient
bound does not establish failure of the actual cone.

## 1. Direction identities with the full leading stress

Use the original normalized coordinates

\[
 t=-b/a,\quad v=a+b^2/a,\quad D=p_1+tp_2-v,\quad
 J=p_2-tp_1,\quad T_0=F(p_1-a,p_2+b).
\]

For N=F(-a,b)/|F(-a,b)| and K=(-Nz,Ntheta), direct dot products give

\[
 H_0=\frac{D}{\sqrt{D^2+J^2}},\quad
 \rho^2=\frac{(v-2)J^2}{2D^2},\quad
 1-\rho^2=\frac{2D^2-(v-2)J^2}{2D^2}.                    \tag{2}
\]

These are profile-scale q=1 identities; common positive physical/chart
factors cancel from direction ratios. They use the COMPLETE p2 and its
pressure terms. Field matching or dropping a pressure constant/eta derivative
would not identify this same T0. The current exact functional five-moment
join supplies that identification; the previously proved heat-sign and
canonical-seed results are consumed, not re-proved or counted as new.

## 2. Quantitative regional direction bounds

On the protected inner part Xa<X<=Xminus, the retained
[core/activation estimate](../core_growth_activation/README.md), section 5,
gives H0>=9/(32C^4) and rho^2<=1/2, including the limiting direction at Xa.
Here modulation is exactly absent, as verified by the regional join.

On [Xminus,Xrep], retain R=4G^3, kappa=G^-4, all four final coordinate
gaps >=kappa/2 and the R coordinate/inverse-a budget. Then

\[
 \kappa/2\le D\le4R^3,\quad |J|\le2R^3,\quad
 2D^2-(v-2)J^2\ge\kappa/2.
\]

Consequently (2) gives

\[
 H_0\ge\kappa/(10R^3)=\omega_*,\qquad
 1-\rho^2\ge\kappa/(64R^6)=\varepsilon_* .               \tag{3}
\]

The denominator 2D^2 is retained. A positive fourth gap by itself is not
the same numerical deficit after normalization.

On the pressure-matched compensated exterior [Xrep,Xc], let p=p1 and
w=p2/p1. The existing heat/edit bounds retain

\[
 p\ge10^{12},\quad 2<v<120,\quad
 P_c/p>.29,\quad |P_c/p|<14000,\quad |w|<1801,\quad
 |b|<15,\ a>2,\quad
 [2D^2-(v-2)J^2]/p^2>.80 .
\]

The 1801 includes the small correction to the reference 1800 bound.
Thus, with ALL finite-p terms retained,

\[
 D/p>.29-120/10^{12}>1/4,\quad D/p\le14000,\quad
 |J|/p\le|w|+|b|/a<2000.
\]

Since sqrt(14000^2+2000^2)<16000, (2) yields

\[
 H_0>1/64000,\qquad
 1-\rho^2\ge\frac{.80}{2(14000)^2}=\frac1{490000000}.     \tag{4}
\]

These bounds apply inside the compensation patch as well as between its
supports. They use the already compensated reference, not the unedited
axial pulse alone, and retain the same-root pressure/moment identity.

On the terminal interval Xc<=X<Xb, the existing
[outer-edge formulas](../outer_edge_stability/README.md) give N=(-1,0),
|T0z/T0theta|<=4h^4, c0^2<=h. Hence for 0<h<.01,

\[
 H_0\ge(1+16h^8)^{-1/2}>1/2,\qquad
 \rho^2\le16h^9<1/2.                                    \tag{5}
\]

All these are relative direction inequalities, not positive amplitude
floors at Xa or Xb. Their smooth one-sided limits are supplied by the
retained flat-edge construction. Taking limits in (3)--(5) is legitimate.

For G>=2, omega*<1/64000 and epsilon*<1/490000000. The inner orientation
comparison is 180G^13/C^4>180C^53244>1. The compact constants therefore
dominate the weaker regional requirements and prove (1) on the CLOSED
annulus. No new compactness-only minimum is left in (1).

## 3. The precise compatibility test for the SAME existing tilt

The original definition fixes one u>0 satisfying
u/sqrt(1+u^2)>rho_actual, the true leading supremum. This is stated in
the manuscript preceding (7.2) and in the retained
[actual-band estimate](../actual_band_pulse_bounds/README.md), section 3.
It is not a numerical assignment of u, and this checkpoint does not
replace it by a newly chosen tilt.

The desired same-tilt squared cone ratio is

\[
 \frac{(1+u^2)c_0^2(K\cdot n_0)^2}{u^2H_0^2}
 \le(1+u^{-2})(1-\varepsilon_*)=1-\delta_u,
\quad \delta_u=\varepsilon_*-(1-\varepsilon_*)/u^2.       \tag{6}
\]

Thus the explicit sufficient condition and a resulting margin are

\[
 u^2>(1-\varepsilon_*)/\varepsilon_*,\qquad
 \gamma_u=\delta_u/2>0.                                  \tag{7}
\]

Indeed sqrt(1-delta)<=1-delta/2. Equality in (7) supplies no strict margin
in this upper-bound class. A convenient padded check is

\[
 \boxed{\ u^2\ge2/\varepsilon_*=2^{19}G^{22}\ }
 \quad\Longrightarrow\quad \gamma_*=\varepsilon_*/4.     \tag{8}
\]

This is a check on the unchanged u, not a construction change. The original
existential admissibility condition does NOT imply (7) or (8), because
rho_actual can be much smaller than the conservative bound in (1).
For abstract data rho_actual^2=1/4 and existing u^2=1, the actual tilted
ratio is 1/2<1, while a coarse epsilon=2^-40 fails (7). This is an exact
logical control, not a counterexample to the selected profile or manuscript.
If (7) is not certified, the result is **uncertified by this bound**, not
actual cone failure; no old first-band or old-tilt claim is silently promoted.

## 4. Join to the original background radius, retaining band factors

The existing leading-frame data are f*,f^a, alpha,Ab,Bb,nu>0. On a full
1/4<=s=q/Q<=4 chart, a safe input tuple for the UNCHANGED
[background gap certificate](../background_gap_certificate/README.md) is

\[
 f=f_*/5,\quad \bar F=5f^a,\quad t_*=\alpha f_*/5,
 \quad\bar G=5f^a(A_b+B_b),\quad d_*=\alpha\nu f_*^2/25.
                                                               \tag{9}
\]

The factors use 1/5<4^(-1-h) and 4^(1+h)<5, as already proved; the
discriminant scales by the SQUARE of the common factor. Direction ratios,
omega* and epsilon* do not scale. If (7) is checked, insert omega*,gamma_u
and the SAME u into the original r_star formula; under (8), gamma* suffices.
This gives an explicit positive r_Kc from the retained spectral inputs.
No pulse frame/carrier is redefined and no new source inverse is asserted.

With the previous C* (including its explicit stress floor and all complete
pressure/jet costs), the compact acceptance inequality remains
2h(ell-2)>=log2(C*/r_Kc), together with every other accepted-tail gate.
The unverified primitive collar/jet/source bounds and quantitative relation
of the EXISTING u to (7) still prevent asserting this actual scale test.
An amplitude floor, direction certificate or passing test is not a full
physical-pulse, complete-covariance or residual-improvement theorem.
The positive-order Pi1/Pi2 memories, full radial-pressure source integrals
and moving mean/covariance pressure terms are not altered or erased here.

## Checks and next obligation

The new evidence contains 81 exact compact controls, 72 finite-p exterior
controls, terminal/hierarchy comparisons, exact fixed-u boundary and
missed-sufficient-gate controls, and a rational substitution into the original
scaled gap map. The rational inputs are not actual enormous-profile data.
Continuum estimates are (1)--(9), relative to the listed leading lemmas;
finite tests are algebra diagnostics, not a formal proof of their premises.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
Next certify the actual activated-collar/jet envelopes and either verify
(7) for the retained u or obtain a sharper same-u raw deficit. Keep all
pressure, mean-covariance, viscosity, cutoff and interaction obligations.
Historical first-band acceptance and unconditional blow-up remain open.
