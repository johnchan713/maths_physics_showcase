# Canonical manuscript primary: separate construction and root stability

Parent: `817edc46cb8ba64db62f3e7c0287c45522ef5eab`.
Status: `conditional-canonical-primary-rebuild-and-relative-root-stability`.
The parent hosted CI was in progress without a failed step when inspected
on 18 October 2026. Historical pulse/evidence files are unchanged.

The preceding checkpoint found that manuscript (7.8) uses reference s(v)
in its mixing matrix, while inherited replays use actual sa. This checkpoint
constructs the manuscript-canonical homogeneous primary and its positive
PRINCIPAL-covariance roots separately, on an explicit common value domain.
It proves quantitative comparison, NOT equality, with the inherited family.
No old physical pulse/root is silently overwritten or identified as new.

All claims below are conditional on the reviewed current phase/background
jets, spectral floors, stress cone and weighted edge hypotheses. There is
no numerical evaluation of the enormous selected background/first band.
Independent finite controls are diagnostics, not a continuum proof.

## 1. Canonical entrance and its own radial and column bounds

Use all notation and primitive bounds of `frame_column_foundation`.
Distinguish `Ba=U Ma`, `wa=sqrt(1+sa^2)` (inherited) from
`Br=U Mr`, `wr=sqrt(1+sref^2)` (manuscript), with
`M(w)=[[1,1],[c0*w,-c0*w]]` and the SAME actual U and normal n.
Write CEa=CE, CEr=CEref from that checkpoint. Both frames have a complete
projected equation, including nprime, Bprime and scalar viscosity:

\[
 z_j'=[diag(\lambda,-\lambda)+E_j-dI]z_j,\quad
 \|E_j\|_2\le CE_j/L,\quad |d-dref|\le Cd/L,
 \quad j=a,r.                                                   \tag{1}
\]

The manuscript Lemma 7.4 specifies `zr(0)=P(0)*(1,0)` in Br.
Construct it by the finite-interval linear ODE; do NOT transport the
inherited entrance and then call that the canonical manuscript entrance.
The inherited pulse independently has `za(0)=P(0)*(1,0)` in Ba.
The reference envelope P and all phase/carrier/cutoff choices are the same.

On `L>=8CEr/mu`, the Riccati barrier gives
`|qr|<=CEr/(mu L)<=1/8`. At either boundary its outward derivative is at
most `-47CEr/(64L)<0`. Thus, with `Dr=Cd+9CEr/8`,

\[
 x_{r,-}P\le x_r\le x_{r,+}P,\qquad
 x_{r,-}={7\over8}e^{-Dr},\quad x_{r,+}={9\over8}e^{Dr}.
                                                               \tag{2}
\]

These are positive RADIAL ODE bounds, not a positive floor for a complete
oscillating cutoff velocity. Both coordinates contribute to direction:

\[
 {t_{r,tan}\over x_r}=-sa\,a+c0\,wr\,{1-qr\over1+qr}\,b.
                                                               \tag{3}
\]

Set Wr=1+Sref. The parent's normalization and reference-frame CE bounds
therefore give the explicit constants

\[
 D0r=Ds+Sref\,Da+cmax\,Wr
          [Da+16CEr/(7\mu)],\qquad D1=u(1+cmax).                 \tag{4}
\]

In contrast to actual-wa mixing there is no additional cmax*Ds cost from
replacing wa by wr in (3): wr is already the reference factor.
The direction error from its midpoint is bounded by
`D0r/L+D1*|v-L/2|/L`. With the parent's Gaussian constants c,C,

\[
 Cdir_r=D0r+D1\,{x_{r,+}^2e^C\over2x_{r,-}^2c},\quad
 h_{r,\sigma}=s0\,j_{r,\sigma},\quad s0={r0^2\over\sqrt L},
\]
\[
 Jmin\,e^{-2/3}x_{r,-}^2e^{-C}\le j_{r,\sigma}
 \le Jmax\,x_{r,+}^2\sqrt{\pi/c}.                                \tag{5}
\]

The angular factor 1/2, actual transverse integral, ci=2r0/L and temporal
cutoff are retained. Require L>=25 and the common cone acceptance
`L>=[4*Cdir_r*max(Amin^-1,u^-1)/gamma]^2`, as well as the corresponding
inherited acceptance and original normal floors. These are finite VALUE
conditions, not a claim that an arbitrary numeric first band is accepted.

The existing cone proof now applies to Hr built from THIS canonical
primary, not to an inherited H relabeled as Hr. With the unchanged target
T, V=(-N.T)/A>0, W=(K.T)/u, |W|<=(1-gamma)V, A=|c0|sqrt(1+u^2),
the unique positive squared amplitudes are

\[
 y_r=H_r^{-1}T,\quad a_{r,\sigma}=\sqrt{y_{r,\sigma}},\quad
 {\gamma V\over8s0\,j_r^+}\le y_{r,\sigma}
 \le {3V\over s0\,j_r^-}.                                      \tag{6}
\]

This is a conditional canonical version of manuscript (7.24)--(7.29).
The same leading principal covariance target is matched exactly. It does
not assert that y_r=y_a or that the complete curl covariance is already
the principal covariance.

## 2. Same physical operator, different entrance: an explicit difference

Let W0=1+S0 and
`B0=2W0*max(1,cmax W0)`, which bounds both frame operator norms.
Let Da*=CEa+Cd and D*=max(CEa,CEr)+Cd. Since
`lambda(v)>0`, the logarithmic norm bound in either frame gives

\[
 |t_j(v)|\le T_*P(v),\qquad T_*=B0e^{D_*}.                       \tag{7}
\]

At v=0 the physical entrance difference is
`P(0)*c0*(wr(0)-wa(0))*b(0)`. In Ba coordinates its two components are
`P(0)*(wr-wa)/(2wa)*(1,-1)`. The parent's |wa-wr|<=Ds/L and wa>=1 give
its Euclidean norm at most Ds*P(0)/L. Propagating this difference by the
SAME physical homogeneous operator, in the Ba coordinates, gives

\[
 |t_r(v)-t_a(v)|\le {C_\Delta\over L}P(v),\qquad
 C_\Delta=B0\,Ds\,e^{Da_*}.                                    \tag{8}
\]

This uses the forward bound
`||Va(v,w)||<=exp(Da* (v-w)/L)*P(v)/P(w)`.
The factor P(0) cancels; no exponentially small entrance is inverted
without its matching propagator factor. No diagonal/commuting error is
assumed. These estimates hold on the whole closed pulse, not just a grid.

Pressure is not discarded. For each homogeneous pulse,
`pi=i(n.Kbg*t-nprime.t)/(k|n|^2)`. With beta=bmin/2, set

\[
 K_*=4(Fmax+DF)+2(gmax+Dg),\quad
 P_*={K_*\over beta}+{Dr'+Dt'\over beta^2}.
\]
\[
 |\pi_r-\pi_a|\le {P_*C_\Delta\over kL}P(v).                     \tag{9}
\]

Here |n|>=beta, ||Kbg||<=K* and |nprime|<=Dr'+Dt' on L>=1.
Thus both the normal derivative and the pressure denominator remain.
For later forced modes, the normal-source term must ALSO be retained;
their zero-data inverse is frame-invariant only with the complete forcing
and time-dependent coordinate transformation already proved by the parent.

## 3. Physical Haar-scale comparison and relative root stability

Write Hj for the two PRINCIPAL columns formed using the same angular
average and actual cutoffs. From (7)--(8),
`|xr*tr,tan-xa*ta,tan|<=2T* CDelta P^2/L`.
The Gaussian square has integral at most sqrt(pi L/c), and
`(J/2)*(integral chi^2)*ci<=Jmax*r0^2/L`. Hence

\[
 |(H_r-H_a)_\sigma|_2\le {s0\over L}C_H,\qquad
 C_H=2Jmax\,T_*C_\Delta\sqrt{\pi/c}.                            \tag{10}
\]

No radius-independent absolute covariance bound is substituted for s0.
The inherited inverse bound implies
`||Ha^-1||infinity<=6/(s0*j_a^-*m)`, m=min(Amin,u).
Since `||Hr-Ha||infinity<=2s0 CH/L`, define

\[
 \theta={12C_H\over j_a^-mL},\quad
 \kappa={24j_a^+\over\gamma j_a^-},\qquad
 L\ge {12C_H\over j_a^-m}(2\kappa+1).                            \tag{11}
\]

The physical s0 cancels in the RELATIVE perturbation only.
If (11) holds, theta<=1/(2kappa+1)<1. A Neumann inverse and
`yr-ya=-(I+Ha^-1 DeltaH)^-1 Ha^-1 DeltaH ya` give

\[
 \max_\sigma\left|{y_{r,\sigma}\over y_{a,\sigma}}-1\right|
 \le R:={\kappa\theta\over1-\theta}\le\tfrac12.                  \tag{12}
\]

Indeed ||ya||infinity/min ya<=kappa by the inherited component floors.
Taking positive square roots and reciprocal roots then gives

\[
 |a_r/a_a-1|\le R,\qquad |a_r^{-1}/a_a^{-1}-1|\le2R.             \tag{13}
\]

The latter follows from
`R/[sqrt(1-R)*(1+sqrt(1-R))]<=2R` for 0<=R<=1/2.
All constants in (11) are explicit finite functions of retained inputs.
They may be enormous; the implementation keeps their logarithms rather
than interpreting floating underflow as a singular inverse.

Reusing ya with the new pulse does NOT recover exact covariance:
`Hr*ya-T=(Hr-Ha)*ya`. Rebuilding yr recovers `Hr*yr=T` exactly.
For the canonical principal waves, disjoint supports therefore give
`C(W0r)=epsilon*T`, and defining the corresponding signed linearization
with Hr and ar gives `B(W0r,Lr Sigma)=Sigma`. These are principal-wave
identities only. They do not identify either entire construction with the
other and do not remove generated complete-curl covariance corrections.

## 4. Fixed jets, supports and what is not yet joined

On the common positive VALUE domain of (1)--(6),(11), frozen discrete
carriers are not differentiated. Differentiating the finite ODE, covariance
integrals, inverse equation Hr yr=T and positive square roots gives the
same finite polynomial S and inverse-edge losses for each fixed jet.
The initial P(0) and its differentiated exponential have only the existing
polynomial losses relative to P(0). The forward Gaussian estimate then
propagates them relative to P(v). No derivative-dependent first band is
introduced. The flat weighted target hypotheses remain necessary at the
shell edges: the root comparison alone does not prove smooth zero extension.
Together with those hypotheses, the existing weighted inverse/root argument
gives smooth edge extension and the same cutoff supports for the canonical
family. The actual oscillating field may still vanish pointwise.

This supplies a canonical PRIMARY lemma and a comparison route, rather
than an unsupported exact identification. Later mean/covariance and repair
maps must consistently use this canonical family and its rebuilt roots,
or consistently use the inherited family with all identities re-established.
No later map is asserted automatically unchanged. Full axial velocity,
Pi1/Pi2 pressure memory, moving pressure moments, all physical scale losses,
Gaussian cutoff commutators, viscosity and interactions remain in that join.
Core/background analytical premises and the complete residual still block
an unconditional pulse-family theorem, complete residual improvement or
blow-up conclusion.

`canonical_bounds.py` checks exact budgets, Riccati acceptance and rational
matrix/root perturbations at two physical radii. `canonical_checks.py`
independently solves the full three-component physical equation and the
canonical two-coordinate equation, keeps the complete pressure balance,
integrates actual cutoff/Haar columns, and detects the stale-root defect.
Its small manufactured phase fixtures are NOT evaluations of the selected
giant construction or certifications of its sufficient value thresholds.

```sh
python3 research/navier_stokes_cascade/results/canonical_primary_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/canonical_primary_foundation/audit.py \
 --verify-record research/navier_stokes_cascade/results/canonical_primary_foundation/evidence.json \
 --output /tmp/canonical_primary_fresh.json
```
