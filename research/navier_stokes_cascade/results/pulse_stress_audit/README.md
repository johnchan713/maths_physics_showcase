# Corrected pulse stress: uniform finite-parameter bounds

This checkpoint supplies an explicit analytic reduction with outward-rounded
constants for the **entire corrected reference pulse**, at every
`eta in [-1,1]`. The sufficient stress tests (A.24) have margins

\[
a-b_sw>1.48,\qquad
2-2b_sw-b_s^2/a-(a-2)w^2>.82.
\]

The same `XR>=100` as the preceding three covered stages is sufficient for
the **actual finite-radius cone on this pulse**. Both axial end corrections,
the angular dependence of the moment-closing strength, the incoming memory,
pressure, energy, and finite positive `h` are included. The error in the
reduced stress relation below is smaller than `1e-48`.

This is not a proof-assistant certificate or verification of the complete
manuscript. The later profile interpolation/exterior stages, heat replacement,
five-moment axis-annulus attachment, higher-order PDE corrections and smooth
forcing remain unverified. No complete Navier-Stokes blowup candidate is
promoted. All previous failed candidates and historical evidence are retained.

## Source, inputs and what is new

Source: the user-supplied [OpenAI manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pages 26-31 and 134-136, especially (4.9), (4.16), Lemma 4.5 and
(A.24), (A.27)-(A.30). PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent: `856f4fff414aeaa27cc2e04a2f09a55222e7d4e9`.

The [moment checkpoint](../pulse_moment_audit/README.md) supplies the exact
corrected reference schedule, with

\[
M_d=64,\ T=e^{64}+10,\ \log P_*=T+1,\ \lambda=e^{-4T},\ h=\lambda^2,
\quad T_w=240T,\ T_f=1000,\ c_o=.001.
\]

The smallness bounds use `T>=64`. Thus the pulse estimate by itself also
holds for the `Md=4` family member; **that does not repair its failed earlier
axial stage**. Only the chosen `Md=64` inherits the preceding axial result.

On `0<=y<=13/lambda`, set `xi=lambda*y` and write

\[
E=\frac{e_b}{1+\eta^2}e^{-(1/2+\lambda)y},\quad
R_b=R+R_e,\quad R=A(\eta)R_0(\xi),\quad U=ER_b.
\]

`A` is pulse strength, not the paper's exponent `Ap=1/2+h`. The inherited
inputs are

\[
1.0100502<A<1.0100504,\quad |A_\eta|\le3000T\lambda,
\quad \|c_i\|_{C^1}\le2\lambda^{20}.
\]

Here `||g||C1=sup|g|+sup|g_eta|`; `Re` is the sum of the two disjoint,
width-.3 bumps `ci*sigma'((y-center_i)/.3+.5)`. Their centers are
`13/lambda-3` and `13/lambda-1`. Set

\[
\bar E=e^{-119.5T+.3},\qquad \mu=24e^{-120.5T+.7}.
\]

The exact `eb` is smaller than `Ebar`. The incoming normalized moment
`m0=M/(Xp Ep)` has C1 norm at most `mu`; the incoming
`s0=S/(Xp Ep^2)` satisfies `|s0|<=500T`, `|s0_eta|<=1`.
The parent proves these inequalities and M/J/S closure, rather than merely
observing small discrete residuals. Its source, code and evidence are hashed.

## 1. Shape and end-patch bounds

The fixed shape is

\[
R_0(\xi)=\phi(\xi)[1-\sigma(\xi-10)],\quad
\phi(\xi)=\int_0^\xi\sigma(v/.02)\,dv.
\]

Extend it by zero to the left. It is flat at both support endpoints.
The inherited global derivative bounds are `sigma'<=64` and
`|sigma''|<=20000`, including the flat endpoints. Consequently

\[
0\le R\le13.2<14,\quad R_\xi\le A,
\quad |R_\xi|\le846,\quad |R_{\xi\xi}|<270000.
\]

Indeed `0<=phi<=11`, `0<=phi'<=1`, `|phi''|<=3200`. Applying the product
rule gives `|R0'|<=1+11*64` and
`|R0''|<=3200+128+11*20000`; multiplication by `A<1.2` gives the bounds.
The one-sided bound `R_xi<=A` follows from the negative cutoff derivative,
not from an absolute derivative estimate.

Disjoint end supports give, for both the function and its eta derivative,

\[
|R_e|, |(R_e)_\eta|\le128\lambda^{20},\qquad
|(R_e)_y|\le(40000/.3)\lambda^{20}.
\]

These are fixed-`y` bumps: we do **not** silently treat their `xi` derivatives
as order one. We handle their convolution separately, so no third derivative
of sigma is needed. The full profile satisfies
`|Rb|<14`, `|Rb_eta|<=33000Tlambda+128lambda^20<1`.

## 2. Expand the averaged velocity with a finite remainder

Put `beta=1/2-lambda` and `m=M/(XE)=AX(U)/E`. Exactly,

\[
m_y+\beta m=R_b.
\]

Taylor-expand only the main profile in its exponential convolution. Since
`integral_0^infinity exp(-beta*v)*v^2/2 dv=beta^-3`,

\[
m=R/\beta-\lambda R_\xi/\beta^2+e_m,\qquad
|e_m|\le2400000\lambda^2+\mu+300\lambda^{20}.
\]

The three terms respectively bound the Taylor remainder, the nonzero initial
memory `m0 exp(-beta*y)`, and the end-bump convolution. These bounds hold
at the start, throughout the cutoff and end patches, and at pulse end.
Since `beta>.49`, also

\[
|m|<30,\quad |m_\eta|<1,\quad
|m_\eta|\le(33000/\beta)T\lambda+\mu+300\lambda^{20}.
\]

The last inequality differentiates the exact convolution once in eta; it
does not assume the amplitude is even or independent of angle.

## 3. Positive angular stress, including every small angle

Use `D=1/2-h`, `d=1-eta^2`, `Jp=2eta/(1+eta^2)` and
`c=D eta Jp=(1-2h)eta^2/(1+eta^2)`. Define

\[
B=2D\eta m+d(m_\eta-J_pm),\quad W=1-EB,\quad
q_0=\frac{\lambda-h+c}{1-\lambda}.
\]

Here `|B|<=61`. Substitution into the exact angular source (4.9) gives

\[
Q_y+(1-\lambda)Q=\lambda-h+c
 +E[-\lambda B+(2h\eta+dJ_p)R_b].
\]

The error source is bounded by `32 Ebar`. We must also retain its initial
condition. The preceding axial expression gives `0<Qa<10`. On the following
unit ramp, `W<=1` and the positive Q has source at most `lambda+c<1`, so
`Qb<11` and `|Qb-q0|<12`. On the intervening power interval the additional
source is `-lambda L K0 exp(-y)`, with `0<=K0<=4`. Integrating exactly gives

\[
|Q(X_p)-q_0|\le16e^{-(1-\lambda)T_w}
\le32e^{-240T}<\bar E.
\]

This uses `lambda*Tw<log(2)`. A second variation-of-constants estimate now
gives `|Q-q0|<=100 Ebar` everywhere in the pulse. Furthermore

\[
q_0\ge(\lambda+\eta^2)/3,\qquad
100\bar E<\lambda/12,
\qquad\boxed{Q\ge(\lambda+\eta^2)/4>0.}
\]

Thus the estimate covers the narrow region `|eta|~sqrt(lambda)`; no angular
grid is used to justify division by Q. The finite `h=lambda^2` is retained
in the exact equations, not replaced by zero.

## 4. Reassemble the full axial stress

Let `s=S/(XE^2)`. On this pulse the exact accumulated integral is

\[
s=e^{2\lambda y}\left[s_0+\int_0^y e^{-2\lambda v}
(R_b(v)^2-1/2)\,dv\right].
\]

Because `2lambda*y<=26`, `Tlambda<1`, `|Rb|<14` and `|Rb_eta|<1`,

\[
|s|\le600e^{26}/\lambda,\qquad |s_\eta|\le15e^{26}/\lambda.
\]

For example, the first bound uses
`exp(26)*(500T+197/(2lambda))`; the second uses
`exp(26)*(1+14/lambda)`. This includes the incoming energy memory and both
end bumps, without discarding any term on the ground that it is tiny.

Before the later angular edits their pressure contribution can be omitted
**only because their complete pressure increment is exactly zero**. On the
remaining unedited continuation `l<=0`, `0<=theta<=1`. The same positive
tail estimates as the prior audit therefore give

\[
|\Pi|/E^2\le1/2,\qquad |\Pi_\eta|/E^2\le |J_p|\le1.
\]

Using `M_eta/(XE)=m_eta-Jp*m` and
`S_eta/(XE^2)=s_eta-2Jp*s`, formula (4.16) becomes exactly

\[
\frac{N}{E}=-R_b+(D+c)m-D\eta m_\eta+\mathcal E_N,
\]
\[
\mathcal E_N=E\left[BR_b+4h\eta s-d(s_\eta-2J_ps)
 +4A_p\eta\frac{\Pi}{E^2}-d\frac{\Pi_\eta}{E^2}\right],
\quad |\mathcal E_N|\le2000e^{26}\bar E/\lambda.
\]

In particular the pressure and energy terms are present. The leading
numerator has absolute value below 45. Replacing the denominator Q by q0
therefore costs at most `54000 Ebar/lambda^2`; dividing `EN` by Q costs
at most `8000 exp(26) Ebar/lambda^2`. Their sum is below
`10000 exp(26) Ebar/lambda^2`.

## 5. A stable, uniform formula for w

Use the exact cancellation `D+c-beta=lambda-h+c`. The main terms give

\[
w=\frac{N}{EQ}=2R-C_dR_\xi+e_w,\qquad
C_d=\frac{\lambda(D+c)(1-\lambda)}{\beta^2(\lambda-h+c)}.
\]

The coefficient of R before rounding is `(1-lambda)/beta`, not exactly two;
its difference contributes at most `42lambda`. The eta derivative uses

\[
\frac{|\eta|}{\lambda+\eta^2}\le\frac1{2\sqrt\lambda},
\]

which follows from `(abs(eta)-sqrt(lambda))^2>=0`. This is why the known
`A_eta=O(Tlambda)` produces `O(Tsqrt(lambda))`, not an uncontrolled O(T)
error. Collecting the nonzero bounds gives

\[
|e_w|\le\epsilon_w:=8000000\lambda+4\mu/\lambda
+1600\lambda^{19}+51000T\sqrt\lambda
+10000e^{26}\bar E/\lambda^2.
\]

Every summand decreases for `T>=64`: the only polynomial factors are
`T exp(-2T)` and `T exp(-4T)`. Outward evaluation at T=64 gives

\[
0<\epsilon_w<8.396\cdot10^{-50}<10^{-48}.
\]

The quotient defining Cd decreases with c because
`lambda-h-D=lambda-1/2<0`. Its maximum is at c=0, where the exact choice
`h=lambda^2` yields `Cd=D/beta^2<2.001`. No 0/0 evaluation is used.
Thus `|w|<28+2.001*846+epsilon_w<1800`.

The exact shear is

\[
a=2+2\lambda,\quad
b_s=-(1+2\lambda)R_b+2(R_b)_y=-R+e_b^s,
\quad |e_b^s|\le1800\lambda.
\]

The last estimate includes `2(Re)_y`; it is not an asymptotic statement
about an uncorrected profile.

## 6. Maximize two scalar quadratics, then check the finite radius

Put `Bstar=2.001*1.0100504`, `delta_b=1800lambda`, `Wmax=1800` and
`epsilon_cross=14 epsilon_w+delta_b Wmax`. Since R is nonnegative and
`R_xi<=A`, even at the large negative cutoff derivative,

\[
b_sw\le-2R^2+B_*R+\epsilon_{\rm cross}
\le B_*^2/8+\epsilon_{\rm cross}<.510612<.52.
\]

The end bumps need not be nonnegative: they are already in the controlled
remainders. Since `a>=2`,

\[
2b_sw+b_s^2/a+(a-2)w^2
\le-3.5R^2+2B_*R+\epsilon_2
\le2B_*^2/7+\epsilon_2<1.167112<1.18,
\]
\[
\epsilon_2=2\epsilon_{\rm cross}+14\delta_b+\delta_b^2/2
+2\lambda W_{\max}^2.
\]

These exact quadratic maxima prove the two announced margins without a
sampled extremum. The negative cutoff derivative helps these upper bounds,
but its full absolute size was retained in all error and compactness bounds.

Positive ratio margins alone are not a finite-radius cone certificate.
Let `p=ps1=XQ/L`, `c*=1-bsw/a`, `j*=w+bs/a`, and `v=a+bs^2/a`.
The exact identity in Lemma 4.5 gives

\[
G:=2(c^*)^2-(v-2)(j^*)^2
=(1+b_s^2/a^2)[2-2b_sw-b_s^2/a-(a-2)w^2]>.82.
\]

The compact ranges above imply `c*>=.74`, `v<101` and `c*v<1300000`.
For `p>=20000000`, `Pc=p*c*>v` and

\[
\frac{2(P_c-v)^2-(v-2)J_c^2}{p^2}
=G-4c^*v/p+2v^2/p^2>.82-.26=.56>0.
\]

Here `X>=Xp=XR exp(T+2+Tw)`, so

\[
p\ge X_p\lambda/4=X_Re^{237T+2}/4.
\]

For `XR>=100`, its logarithmic lower bound at T=64 already exceeds
`log(20000000)`. Finally `v-2>=2lambda>0` everywhere, including R=0.
The audit records `log(2lambda)` separately; subtracting rounded values of
v and two is not a legitimate positivity test. Lemma 4.5 now supplies the
strict admissible cone on this finite corrected pulse.

## Independent checks and failures retained

The numerical controls are **not globally admissible reference candidates**.
Their purpose is to make the finite terms large enough to catch coding errors.

- A manufactured finite-parameter pulse uses `lambda=.02`, `h=lambda^2`,
  a non-even strength, signed end corrections and arbitrary smooth incoming
  moments. Its pressure is the exact infinite constant-slope tail.
  DOP853 evolves eight normalized moments/eta derivatives alongside Q and
  N/E from the separate source equations (4.9). The latter are compared with
  (4.16) across start, ramp, cutoff and both end patches, with two step sizes.
  The shape and required radial velocity are shared inputs, so this is an
  independent equation-path check, not a wholly independent PDE solver.
- At separate manufactured states, high-precision automatic differentiation
  of the full moment formulas is compared with the direct source formulas.
  Omitting finite h, eta derivatives, or energy produces detectable errors.
  Omitting signed end bumps changes the shear and is rejected.
- An exactly solvable affine **particular solution** at the `Md=4` lambda
  compares direct cancellation with its algebraically stable expression.
  At eta=0, 80 digits give the false value w=0 instead of about -1.01005;
  180/260 digits agree with the stable result. Omitting the first convolution
  derivative changes w by about 2.02010. This particular solution excludes
  startup transients by definition; its numerical value is not presented as
  the full pulse. The continuum bound above retains those transients.
- The narrow angular region is also exercised at
  `eta/sqrt(lambda)=0, .5, 2, 10`. Binary64 loses `a-2` even when lambda
  itself remains representable. At the actual Md=64, binary64 also underflows
  lambda to zero. Positive logarithmic parameter bounds remain the evidence.
- A deliberate small-radius example passes both ratio inequalities but has
  `Pc<v`, demonstrating why the finite-radius check cannot be skipped.
- A first unit-test run falsely failed the angular inequality at its exact
  equality point using 120-digit arithmetic. The comparison is now reduced
  to exact fractions, `t/(1+t^2)<=1/2`, and the original positive rounding
  overshoot is retained as a control. No tolerance or mathematical bound was
  relaxed to make that test pass.
- The old `Md=4` axial failure, old later-lambda failure, and failed unforced
  `.14` resolution comparison are not repaired or overwritten by this result.

Quadrature and ODE residuals are diagnostic errors, not rigorous integration
error bars. The continuum claims come from the analytic inequalities and
outward constants, conditional on the explicitly inherited moment argument.

## Reproduce and read the code

With the packages in `requirements.txt` installed, from the repository root:

```bash
python research/navier_stokes_cascade/results/pulse_stress_audit/test_audit.py -v
python research/navier_stokes_cascade/results/pulse_stress_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/pulse_stress_audit/evidence.json \
  --output pulse_stress_reproduced.json
```

Each command keeps the old record unchanged. The first exercises individual
identities, assumptions and failed controls. The second rebuilds the evidence,
checks every required gate and its source hashes, compares it with the frozen
record, then writes the new reproduction to the requested path.

`bounds.py` follows sections 1-6 in order. `point` constructs an outward
interval; `lo` and `hi` extract its endpoints. `pulse_bound` first rejects
unsupported inputs, then constructs universal parameter envelopes at T=64.
`w_parts` keeps each positive error source separate. `Cd`, `linear`, `bsw`
and `second` evaluate the two quadratic envelopes. The `checks` dictionary
tests the numerical constants in the derivation; a false or indeterminate
entry prevents a passed status. The returned logarithms describe the actual
Md, while the smallness constants apply uniformly for every T>=64.

`diagnostics.py` uses the following equation paths:

| Function | What each block computes |
|---|---|
| `step_jet`, `shape_jet` | Flat endpoints first; logistic step derivatives inside; short primitive only where needed; product-rule pulse derivative. |
| `moment_rhs` | Eight derivatives for m, m_eta, I/(XH), its eta derivative, J/(XHE), its derivative, s and s_eta. The different decay rates come from their different denominators. |
| `moment_stress` | D, d, Jp and W; the complete Q moment formula; energy and pressure terms; the complete N/E formula. |
| `source_stress` | The transport factors W and Hc; the separate angular and axial source equations, including all pressure terms. |
| `SyntheticPulse.profile` | Main strength and eta derivative; two signed bumps and their radial/eta derivatives; the decaying swirl E. |
| `initial`, `rhs`, `samples`, `integrate` | Consistent initial moments and stresses; source evolution; prescribed observation points; segmented integration and comparisons. |
| `algebra_controls` | Differentiate the moment expressions at fixed states; compare with source equations; quantify omitted-term failures. |
| `affine_control` | Construct tiny finite lambda and h; compute both cancelling and stable formulas; retain the first derivative term and positive a-2 separately. |
| `finite_cone_controls` | Check the quadratic identity; preserve the counterexample with insufficient radius. |
| `angular_equality_control` | Check the scaled boundary with exact fractions; record the false floating-point comparison. |

`audit.py` hashes its own instructions and implementation plus inherited
evidence, computes independent controls, applies the fixed gates, and writes
`evidence.json`. Its validator also recomputes the scientific predicates from
recorded data: changing a failed quantity without changing a boolean is not
accepted. `test_audit.py` includes these evidence-tampering checks. Tiny
physical quantities use relative comparison without a unit floor; numerical
residuals have an explicit noise allowance.

## Next gate

Verify the post-pulse interpolation and exterior transition with the same
finite parameter schedule, including the actual angular correction patches.
Only after all required reference stages pass should the heat replacement
and five-moment axis-annulus attachment be promoted to the next construction
gate. Full PDE residual correction and smooth forcing remain separate tasks.
