# Q191: physical-pulse chart, pressure and scaling seed audit

Parent: Q190, `f8d3b963ee0b684e22ea43c58b58dd1a7d59b65f`.
The live branch, Q187, Q189--Q190, the leading-frame, weighted-stress and
physical-pulse providers, their frozen evidence, and the current Q189/Q190 CI
runs were read first. Historical evidence and the selected Q189 five-moment
full-pressure profile remain unchanged.

Q190 identified the primary physical-pulse all-band estimate as the first
fresh independent derivation after the same-profile compact object. This
checkpoint reconstructs its low-order and scaling seeds exactly. It does not
promote the provider's all-fixed-jet conclusion: the actual mixed frame/graph
jet estimate stated in Section 7 remains to be derived independently.

## 1. Dimensionless chart and Jacobian

Keep

\[
 h={1\over200},\quad A={1\over2}+h,\quad D={1\over2}-h,
 \qquad s-Z^2s^{2h}=T.
\]

The physical rescaling loses no power of `Q` because

\[
 2D+2h=1.
\]

Writing `eta=Zs^{-D}`, the implicit derivative denominator is

\[
 L=1-2h\eta^2\ge 1-2h={99\over100}>0
 \qquad (|\eta|\le1).                              \tag{1}
\]

Thus this low-order chart inversion has a uniform rational floor. Equation
(1) is an analytical identity, while the executable checks it at exact
rational samples; those samples do not establish every higher chart jet.

## 2. Carrier and pulse-length geometry

Put `x=epsilon^(-1/2)>=1` and `k=ceil(x)`. Since `x<=k<=2x`,

\[
 1\le \epsilon k^2\le4.                            \tag{2}
\]

For `r0=1/6057072` and `Tg=4+sqrt(2)<6`, the actual length choice satisfies

\[
 2r_0S\le L_s<2r_0T_gS<12r_0S<S.                  \tag{3}
\]

The last inequality follows already from `12r0<1`. Nearest *nonzero* integer
rounding only gives `|p-ptilde|<=1/k`; the stronger `1/(2k)` estimate is not
available when the nonzero constraint changes the nearest integer.

## 3. Complete constrained pressure

For the moving normal constraint `n.t=0`, start from the retained amplitude
equation

\[
 t'+Kt+m^2dt+ikm\,n\pi=-f.
\]

Differentiating the constraint and solving without dropping the moving-normal
or source-normal terms gives

\[
 \boxed{\pi={i\bigl(n\!\cdot Kt-n'\!\cdot t+n\!\cdot f\bigr)
                   \over km|n|^2}}.                \tag{4}
\]

The exact rational fixture has both `n'.t` and `n.f` nonzero and verifies
`n'.t+n.t'=0`. It is a sign and term-completeness diagnostic for (4), not a
numerical proof of the continuum frame estimates.

## 4. Riccati continuation

Let `e=C_E/L_s`, `r=e/mu<=1/8`. At either boundary of the invariant interval
for the stable-to-unstable ratio, the outward derivative is bounded by

\[
 -2\mu r+e(1+2r+r^2)
 =e(-1+2r+r^2)
 \le-{47\over64}e<0.                               \tag{5}
\]

This proves the stated low-order nonvanishing barrier once the coefficient
error `C_E/L_s` is supplied. It does not itself derive that error for all
mixed derivatives.

## 5. Physical powers

Set

\[
 s_x={1\over2}+{h\over2},\qquad
 s_t=1+{3h\over2}.
\]

The phase derivative losses dominate the coefficient losses because

\[
 s_x>{1\over2}+h\kappa,\qquad s_t>1+h
 \quad (0<\kappa<1/2).                              \tag{6}
\]

For fixed derivative orders `(a,b)`, the reconstructed losses are

\[
 \begin{array}{c|c}
 \text{quantity}&\text{power lost}\
 \hline
 u_p&(a+1)s_x+bs_t\\
 p_p&(a+2)s_x+bs_t\\
 {\cal A}_p&as_x+bs_t.
 \end{array}                                       \tag{7}
\]

The zero-order rescalings close exactly:
`-A+h/2=-sx`, `-2A+h=-2sx`, and `1/2-A+h=0`.

## 6. One later-band threshold for each fixed request

Every fixed polynomial and physical loss can be absorbed by the Gaussian
band gain. For `a>0` and finite `d`,

\[
 d\ell-a\ell^2={d^2\over4a}
                 -a\left(\ell-{d\over2a}\right)^2
 \le {d^2\over4a}.                                 \tag{8}
\]

Therefore each fixed derivative request admits a finite uniform constant and
one threshold covering every later band. The constant and threshold may
depend on that fixed request; no stage-uniform all-index constant is claimed.

## 7. Precise remaining inequality

The next analytical obligation is to prove, for every fixed mixed index
`alpha` on the actual Q189 representatives and on one common later-band
domain,

\[
 |D^\alpha E_{\rm frame}|+|D^\alpha(d-d_{\rm ref})|
 \le {C_\alpha S^{p_\alpha}\over L_s},              \tag{9}
\]

with finite `C_alpha,p_alpha`, retaining derivatives of `n`, `B`, `Bleft`,
the graph coefficients and the complete pressure denominator. Equation (9)
is the input required by differentiated Duhamel and Faà di Bruno recurrences
for every fixed physical curl, pressure and potential jet. The named physical
pulse provider supplies it relative to its inputs; Q191 has not independently
reconstructed it.

Accordingly this checkpoint proves no unconditional actual-pulse hierarchy,
complete residual improvement, admissible force, competitor exclusion or
Navier--Stokes blow-up. Passing tests verify exact identities, dependency
hashes and numerical diagnostics; they are not a continuum proof.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q191_physical_pulse_seed_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q191_physical_pulse_seed_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q191_physical_pulse_seed_audit/evidence.json \
  --output /tmp/q191_physical_pulse_seed.json
```
