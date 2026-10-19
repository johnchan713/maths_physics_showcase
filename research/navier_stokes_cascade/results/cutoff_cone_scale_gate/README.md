# Correct retained collar width and the combined cutoff/cone/scale gate

Status: `corrected-collar-width-and-max-affine-cutoff-cone-scale-gate`.

Parent: `4de29921a9af3fc05e5ade7badb9a0c2361919ba`.
The preceding global axial/source and swirl checkpoints supply all primitive
leading-field inputs to the first-stress bound, relative to their retained
lemmas. We now assemble that bound into four explicit affine inequalities in
the **unchanged** coefficient cutoff's shape norms. This gives a conditional
admissible region and a finite later-band witness, not acceptance of the old
cutoff, tilt or scale. No profile, root, support, bump, frequency, pressure
gauge, mean correction or construction parameter is changed.

## 1. Correction: the coefficient collar is smaller than the native width

The frozen [cutoff derivative audit](../cutoff_repair_derivative_budget/README.md),
section 1 and protocol, mislabels `tc=C^-20`. The original
[core activation](../core_growth_activation/README.md), (13), and
[axis attachment](../axis_core_attachment/README.md), activation section, give

\[
 t_1=C^{-20},\qquad s_c=(10\sqrt{\log C})^{-1},\qquad
 t_c=t_1s_c,\qquad \Delta=t_c/16=
 \frac{C^{-20}}{160\sqrt{\log C}}.                     \tag{1}
\]

The annular coefficient uses the retained offsets
`Xinner=Xa exp(tc/16)`, `Xkeep=Xa exp(tc/8)`,
`Xcut=Xa exp(3tc/16)`, before the analytic-collar edge `Xa exp(tc/4)`
and leading modulation start `Xa exp(tc/2)`. Its logarithmic keep/cut width
is **Delta in (1)**. These are not new locations. In particular this
checkpoint does not shrink or retune the prescribed native width t1.

The historical derivative *formulas* remain correct when their tc argument
is the actual tc. Its frozen example was explicitly auxiliary rational
data, not the actual construction; passing it did not verify the erroneous
identification. Historical README, protocol and evidence are left byte-for-
byte unchanged. Their source integrity is not an endorsement of every
statement. This new checkpoint supersedes that width identification.

For `K(s)=kappa(Xkeep exp(Delta*s))`, let L1,L2 bound `sup|K'|,sup|K''|`.
The old qualitative conditions do not provide these upper bounds. Retain

\[
 k_1=\frac{L_1}{\Delta X_k},\qquad
 k_2=\frac{L_2/\Delta^2+L_1/\Delta}{X_k^2},\quad X_k=X_{keep}.
                                                               \tag{2}
\]

At the same Xk and shape norms, replacing the erroneous width t1/16 by
Delta multiplies k1 by `g=10 sqrt(log C)`. If the old k2 terms are B2,B1,
the corrected value is `g^2 B2+g B1`, not just g times their sum. This
comparison holds Xk fixed to isolate the width error; it does not assign
the incorrect radii to the actual field. Any actual acceptance calculation
using the wider identification must be recomputed. None had been certified.

## 2. The exact affine jet budget

Use the explicit collar norms `J,Jx,Jxx` from the first-seminorm checkpoint,
the same full-memory linear repair coefficient bound P, and `r=sqrt(rho3)`.
The five original **additive**, unit-mass, logarithmic width-.1 bumps have
`B0=320,B1=121600,B2=220000000`. With `p=r^2`, put

\[
\begin{split}
 s_0&=J+PB_0\max(1,r^{-1}),\\
 x_0&=J_x+P\max(B_1/p,(B_1+B_0/2)/(pr)),\\
 z_0&=J_{xx}+P\max((B_2+B_1)/p^2,
                   (B_2+2B_1+3B_0/4)/(p^2r)),\\
 x_1&=J/(\Delta X_k),\\
 z_1&=2J_x/(\Delta X_k)+J/(\Delta X_k^2),\qquad
 z_2=J/(\Delta^2X_k^2).
\end{split}                                                     \tag{3}
\]

Then the inherited factorial angular-C1 jet bounds are exactly

\[
 S=s_0,\qquad S_x=x_0+x_1L_1,\qquad
 S_{xx}=z_0+z_1L_1+z_2L_2.                                    \tag{4}
\]

These are identities between upper-bound expressions, not equalities for
the actual field. P is independent of shape *derivatives*: the tentative
moments use `0<=kappa<=1`, J and Xcut, while full future source totals Q0,Q1
are bounded by the signed-pressure identities through Xv, including the
axis viscosity term. No signed total is confused with an integral of an
absolute value. The separately established pointwise w0 bound enters the
stress estimate. Neither pressure nor velocity is cut off directly.

## 3. Four affine rows, with all source terms included

Set `a=sqrt(Xa), c=sqrt(Xc), X=c^2, L=1-2h`, and let tau be the retained
positive stress floor **only on [Xinner,Xc]**, not at Xa or Xb. The first
stress is zero through Xinner. Write the complete source bounds from
[first_seminorm_bound](../first_seminorm_bound/README.md) as

\[
 \Theta=T_0+T_1L_1+T_2L_2,\qquad
 Z=Z_0+Z_1L_1+Z_2L_2.                                      \tag{5}
\]

Here T0,Z0 are the original `source_bounds` evaluated at
`(S,Sx,Sxx)=(s0,x0,z0)` with ALL leading inputs
`(f,u,fx,ux,v,w,CF,CU)` retained. Explicitly

\[
\begin{split}
 t_x&=X(1+2u)/L+4\nu+v,& z_x&=X(1+2u)/L+2\nu+v,\\
 T_1&=t_x x_1+2\nu Xz_1,& Z_1&=z_x x_1+2\nu Xz_1,\\
 T_2&=2\nu Xz_2,&Z_2&=2\nu Xz_2.
\end{split}                                                     \tag{6}
\]

The pressure term, including w, and the shifted leading viscosities CF,CU
are constants in the shape variables and occur ONCE, in T0,Z0. The
implementation computes slopes by the same linear operator with only those
three inhomogeneous inputs set to zero. This is affine decomposition, not
their omission from the full estimate. Actual viscosity is nu=1.

The regional constant is the maximum of these four rows:

\[
\begin{split}
 R_u&=1,\\
 R_v&=5(2cs_0+1/2)/a,\\
 R_d&=5\{2c[(2c+1)(x_0+x_1L_1)+s_0/a+1]
                     +(2cs_0+1/2)/a\},\\
 R_w&=\tfrac12+\{(c^2/2)\Theta+cZ\}/\tau,\\
 C_\star(L_1,L_2)&=\max(R_u,R_v,R_d,R_w).
\end{split}                                                     \tag{7}
\]

This is exactly the retained compact-gate bound after (4), not an additional
loss replacing it by an unspecified constant. All affine coefficients are
nonnegative. In particular the stress row has the strictly positive slope

\[
 A_2=\partial_{L_2}R_w=
 \boxed{\frac{\nu c^3(c+2)J}{\tau\Delta^2X_k^2}>0}.            \tag{8}
\]

Positivity assumes the actual positive J and viscosity. The full X factors,
radial stress integration measures, inverse radii, the shared higher-order
Borel half-budgets and the chart factors are retained. The global F, U,
radial, pointwise source and shifted-operator inputs come from the two
preceding checkpoints; neither a fourth mixed derivative nor a new root is
needed here. All of this remains relative to their core/regional lemmas.

## 4. Precise remaining inequality and a conditional later-band witness

For the SAME existing tilt u_tilt, first check

\[
 \delta_u=\varepsilon_*-(1-\varepsilon_*)/u_{tilt}^2>0,
 \qquad \gamma_u=\delta_u/2.                                \tag{9}
\]

Use the same chart-floor tuple and omega*,gamma_u,u_tilt in the original
background-gap certificate, giving r=r_Kc>0. Here u in (6) is a norm of the
axial velocity, **not** u_tilt. Their distinction matters. Original
existential tilt admissibility does not imply the conservative test (9);
failure of (9) is not proof of actual cone failure.

At an existing scale ell, define `B=r*2^(2h(ell-2))`. The four explicit
inequalities

\[
 \boxed{R_i(L_1,L_2)\le B\quad(i=u,v,d,w),\quad L_1\ge1,L_2\ge4}
                                                               \tag{10}
\]

are the exact feasible region of **this sufficient bound**. They are NOT
necessary and sufficient conditions for the actual PDE cone. If
`Rd=D0+D1 L1`, `Rw=A0+A1 L1+A2 L2`, the derivative/stress constraints are

\[
 L_1\le(B-D_0)/D_1,\qquad
 L_2\le(B-A_0-A_1L_1)/A_2.                                  \tag{11}
\]

Additionally Ru,Rv<=B. If any row at the necessary lower corner (1,4)
exceeds B, no upper enclosure satisfying those lower bounds passes this
certificate. Conversely (1,4) passing does not construct a smooth K or
identify the actual shape. The bounds 1,4 are necessary, not sufficient.

For ANY supplied finite valid shape upper bounds and verified (9), set

\[
 m=\max(0,\lceil\log_2(C_\star/r)\rceil),\qquad
 \ell_{suf}=2+\left\lceil\frac{m}{2h}\right\rceil.            \tag{12}
\]

Then Cstar<=r*2^m and `2h(ell_suf-2)>=m`, so the whole `q<=4Q`,
`Q=2^-ell_suf` chart satisfies the compact perturbation gate, provided the
unchanged Borel tail is valid there. Every later band does too. The +2
cannot be dropped. This rational dyadic witness is conservative, not the
smallest admissible integer scale and not a new selection of the existing
scale. It does not discharge other geometry/freezing/inverse/pulse gates.
For the huge actual hierarchy, (3)-(12) are symbolic finite expressions;
we do not evaluate or round them to numerical values.

## 5. Why the old qualitative cutoff specification cannot close (10)

The same smooth monotone family `K_n=1-sigma(n(s-1/2)+1/2)` fits the allowed
window and endpoints, yet `L1>=8n,L2>=16n^2`. Thus (8) gives

\[
 R_w(K_n)\ge A_0+8A_1n+16A_2n^2\longrightarrow\infty.        \tag{13}
\]

For any fixed finite B, choosing `n=1+ceil(B/(16A2))` already gives
`16A2*n^2>B`. This is a rigorously quantified obstruction to a **uniform
upper-bound certificate across unspecified shapes at that fixed scale**.
It is not divergence of the actual stress, not proof that the retained
cutoff is bad, and not an obstruction to constructing a suitable shape or
taking a later scale after it is fixed. No compressed K_n replaces kappa.
Likewise L1=32,L2=1216 below are auxiliary enclosures, not a newly selected
shape for the unchanged field.

## Checks and next obligation

Twelve exact rational controls independently compare the affine jet,
complete source and max-row expressions against the original interfaces.
Further controls test the corrected width factors, positive viscous slope,
full inhomogeneous pressure/viscosity costs, exact shape-region boundary,
strict exterior rejection, fixed-budget compression and a same-tilt gap/
dyadic composition. Symbolic checks retain the different g and g^2 costs.
These are algebra diagnostics, not samples of the actual enormous profile
and not a continuum proof. Recursive historical source hashes are checked;
all old evidence and original comparison tolerances are preserved.

Next identify valid enclosures for the UNCHANGED cutoff shape and either
verify (9) for the unchanged tilt or prove a sharper same-tilt margin, then
test (10) at the retained scale with the corrected tc. Conditional existence
of a finite later witness is not verification of the old first band.
Uniform actual physical pulse size/derivative/support and inverse-gap
bounds, complete wave-covariance means including axial velocity/pressure,
and the complete viscous/cutoff/interaction residual remain open. No
unconditional pulse bound, complete residual improvement or blow-up follows.

Sources: linked repository derivations and the retained
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Lemma 5.2 Step 1 and the fixed-tilt setup before (7.2); SHA-256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its proposed theorem is not accepted as established.

```sh
python3 research/navier_stokes_cascade/results/cutoff_cone_scale_gate/test_audit.py -v
python3 research/navier_stokes_cascade/results/cutoff_cone_scale_gate/audit.py \
  --verify-record research/navier_stokes_cascade/results/cutoff_cone_scale_gate/evidence.json \
  --output /tmp/fresh_cutoff_cone_scale_gate.json
```
