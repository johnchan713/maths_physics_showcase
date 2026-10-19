# Exact repair derivatives and the missing cutoff-shape budget

Status: `exact-repair-bump-jets-and-conditional-cutoff-shape-budget`.
Parent: `ec57a0079d884edb2b1be30a587f98455628b1ef`.

This checkpoint bounds the actual five positive-order repair bumps through
two logarithmic derivatives and inserts those bounds into the retained first
coefficient estimate. It also identifies a missing quantitative input: the
inner coefficient cutoff has specified locations and smoothness, but no
specified transition shape or derivative budget. Its existence and finite
derivatives are not in doubt. A numerical upper bound for its derivatives
does not follow from support locations alone, even if monotonicity is added.
This is an obstruction to that attempted *uniform estimate*, not an
obstruction to the existence of a suitable cutoff or to the whole approach.

The [preceding global axial estimate](../global_axial_angular_bound/README.md)
supplies the two signed full-pressure source moment bounds conditionally on
the retained regional/core lemmas. The [first seminorm estimate](../first_seminorm_bound/README.md)
already exposes separate cutoff and bump inputs; none of its frozen evidence
is corrected or overwritten. We now discharge the bump inputs and make the
remaining cutoff dependence explicit. No existing root, cutoff, support,
frequency, width, pressure gauge or mean correction is changed.

## 1. What is specified, and what is not

The [annular extension](../annular_first_correction/README.md), section 1,
fixes `Xkeep=Xa exp(tc/8)` and `Xcut=Xa exp(3tc/16)`, where `tc=C^-20`.
It chooses a smooth eta-independent radial kappa, one through Xkeep and
zero from Xcut. Manuscript Lemma 5.2, Step 1, adds `0<=kappa<=1`, but also
does not select a formula. Source: the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`,
printed page 50. This source's proposed theorem is not treated as established.

Set

\[
\Delta=t_c/16,\qquad s=\frac{\log(X/X_{\rm keep})}{\Delta},
\qquad K(s)=\kappa(X_{\rm keep}e^{\Delta s}),\quad 0\le s\le1.
\]

Let `L1=sup|K'|`, `L2=sup|K''|`. These are finite for a fixed smooth cutoff,
but are not supplied parameter-only upper bounds. Exact chain rules give

\[
\kappa_X=\frac{K'}{\Delta X},\qquad
\kappa_{XX}=\frac{K''/\Delta^2-K'/\Delta}{X^2}.
\]

Consequently valid conditional radial budgets are

\[
k_1=\frac{L_1}{\Delta X_{\rm keep}},\qquad
k_2=\frac{L_2/\Delta^2+L_1/\Delta}{X_{\rm keep}^2}.
\tag{1}
\]

The second term in k2 is the logarithmic-coordinate Jacobian cost and cannot
be discarded. Since K drops by one and K' vanishes at both ends, necessarily
`L1>=1`, `L2>=4`: indeed `1<=int_0^1 |K'|` and
`|K'(s)|<=L2 min(s,1-s)`, whose integral is L2/4. These are **lower** bounds,
not upper bounds and not sufficient to validate any candidate budgets.

## 2. Why the missing shape cannot be inferred from support

Use the retained squared-flat step sigma, extended by zero/one outside [0,1],
and consider the auxiliary family

\[
K_n(s)=1-\sigma\bigl(n(s-1/2)+1/2\bigr),\qquad n=1,2,\ldots.
\]

Every member is smooth, monotone, between zero and one, with the same
allowed keep/cut window. Since `sigma'(1/2)=8`,
`|K_n'(1/2)|=8n`. Its derivative is zero at `s=1/2+1/(2n)`, so the mean value
theorem gives `sup|K_n''|>=16n^2`. Thus no finite bounds for L1 or L2 follow
from those qualitative constraints and the fixed window alone. This family
is an exact analytical countercontrol, **not** a replacement for the retained
cutoff, and not evidence that its actual derivatives are infinite.

One permissible *additional choice* would be `K=1-sigma`, for which the
bounds below imply L1<=32,L2<=1216. The repository has not identified its
previously fixed kappa with that formula. We therefore do not assign those
values to the unchanged field or claim its remaining scale gate is closed.
Specifying such a choice would be a new construction specialization, and
would require tracing affected coefficients, evidence and acceptance gates.

## 3. Continuum bounds for the actual five repair bumps

Unlike kappa, the repair bumps have an explicit retained shape; see
[matching roots](../matching_root_foundation/README.md), section 4, and
`stress_realization_audit/repair.py`. Write

\[
\sigma(s)=\frac{e^{-1/s^2}}{e^{-1/s^2}+e^{-1/(1-s)^2}},\quad
\beta_c(y)=\frac1w\sigma'\left(\frac{y-c}{w}+\frac12\right),\quad w=1/10.
\]

The exact unit mass is in **dy**, not dX or dR. U centers are 1,2 and E
centers are 1/2,3/2,5/2. Their supports are disjoint and lie between .45 and
2.55 in `y=log(X/rho3)`, strictly inside I3=(0,5), before I4=(6,11).
Positive-order corrections are additive; they are not multiplied by the
leading eta-dependent amplitude e.

Here is a direct global bound, without quadrature or giant-parameter
sampling. For `0<s<=1/2`, put `t=(1-s)^-2-s^-2`, `z=s^-2>=4`.
Then sigma is the logistic of t, `sigma<=exp(4-z)`, and

\[
|t'|\le4s^{-3},\quad |t''|\le12s^{-4},\quad |t'''|\le48s^{-5}.
\]

The logistic derivative factors have absolute bounds one: in particular
`1-6sigma+6sigma^2=6(sigma-1/2)^2-1/2` lies in [-1/2,1]. Therefore

\[
\begin{aligned}
|\sigma'|&\le e^{4-z}4s^{-3}\le32,\\
|\sigma''|&\le e^{4-z}(12s^{-4}+16s^{-6})\le1216,\\
|\sigma'''|&\le e^{4-z}(48s^{-5}+144s^{-7}+64s^{-9})<220000.
\end{aligned}
\tag{2}
\]

To verify the last numerical enclosure, `z^(p/2)e^-z` is decreasing on
z>=4 for p=5,7. For p=9 its maximum occurs at 9/2, and is less than
`5^5 e^-4`. Thus the third bound is at most
`48*32+144*128+64*5^5=219968<220000`. The first two bounds use the same
calculus for powers 3,4,6. Symmetry handles s>=1/2, and flatness handles the
endpoints and exterior. These are analytical real-variable bounds, not
sampled maxima; they conservatively sharpen previous upper bounds without
altering the bumps or earlier evidence.

The resulting actual logarithmic bump budgets are

\[
(B_0,B_1,B_2)=\left(\sup|\beta|,\sup|\beta_y|,\sup|\beta_{yy}|\right)
\ \le\ (320,121600,220000000).
\tag{3}
\]

## 4. Insert into the same first-coefficient repair

Use the inherited C_eta^1 collar bounds `J,Jx,Jxx` and linear repair
coefficient norm Cp. Cp is still the inverse-800 estimate with the full
`4/e_star` reciprocal-amplitude derivative cost and the `lambda_patch^-1`
conditioning. The moments in Cp must retain the complete future pressure
source; the newer signed bounds Q0,Q1 may bound the signed incoming moments,
but not the absolute or pointwise source. Nothing here erases that memory.

Take `rho=rho3` as a lower X radius on every bump; then `sqrt(rho)` is a
valid lower bound for sqrt(2X). Write r=sqrt(rho). Combining (1)-(3) with
the unchanged `first_seminorm_bound.extended_jets` gives

\[
\begin{aligned}
S&=J+320C_p\max(1,r^{-1}),\\
S_x&=J_x+k_1J+C_p\max(121600/\rho,121760/(\rho r)),\\
S_{xx}&=J_{xx}+2k_1J_x+k_2J\\
&\quad+C_p\max(220121600/\rho^2,220243440/(\rho^2r)).
\end{aligned}
\tag{4}
\]

Disjoint supports avoid a factor of five. The extra inverse-radius terms
come from `F=E/sqrt(2X)`, with

\[
F_{XX}=\frac{\beta_{yy}-2\beta_y+3\beta/4}{X^2\sqrt{2X}}
\]

for an E bump. The U-bump second derivative is
`(beta_yy-beta_y)/X^2`. Angular coefficient derivatives are already in Cp;
all cutoffs and bumps here are eta-independent. S bounds values from the
axis; Sx,Sxx are used only beyond Xa as in the retained estimate. Only F1,U1
are cut and repaired. M1,V1,Pi1 are reconstructed from their exact integrals,
with the original zero-axis pressure gauge; no pressure cutoff is introduced.

Equations (1),(4) are explicit **conditional** inputs to the complete first
stress/source estimate and its regional Cstar. The precise missing cutoff
condition is to supply valid enclosures for the actual L1,L2 and then check
the same-field scale inequality

\[
2h(\ell_{\rm scale}-2)\ge\log_2(C_\star(L_1,L_2)/r_{K_c}),
\]

together with the separate leading pointwise source, shifted axial,
fixed-tilt and all other acceptance conditions. Existence of finite shape
norms alone does not verify this inequality at a previously fixed scale.

## Evidence and limits

`audit.py` checks symbolic chain rules, the exact central step jets
(8,0,-3712), rational derivative enclosures, all five retained support
intervals, the compressed-cutoff countercontrols and composition with the
existing first-jet interface. Rational propagation fixtures use auxiliary
values and are not the enormous selected profile. Historical byte bindings
and strict scope guards are reproduced. Passing code is not a continuum
existence proof; the analytical arguments and inherited premises above
determine the scope.

Next: identify or bound the unchanged cutoff shape, and quantify the
remaining actual pointwise leading w0, shifted axial and radial inputs before
joining Cstar to the fixed tilt and scale gates. Uniform physical-scale pulse
bounds, complete covariance/mean corrections and the full viscous residual
remain unclosed. No unconditional pulse bound, residual improvement or
blow-up is claimed.

```sh
python3 research/navier_stokes_cascade/results/cutoff_repair_derivative_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/cutoff_repair_derivative_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/cutoff_repair_derivative_budget/evidence.json \
  --output /tmp/fresh_cutoff_repair_derivative_budget.json
```
