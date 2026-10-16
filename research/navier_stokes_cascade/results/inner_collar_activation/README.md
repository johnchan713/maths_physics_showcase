# Inner activation: a width-rescaled flat quotient and a cubic direction error

Status: `inner-flat-activation-and-conditional-cone-audited`.

This batch independently derives the inner activation estimate from the
forward five-moment formulas, retaining pressure and parameter derivatives.
After rescaling the activation interval, field differences have the factor
`t1*ea*s^3`, cumulative differences have `t1^2*ea*s^6`, and the normalized
stress-direction error is `O(t1*s^3)`. The resulting squared cone loss is
`O(t1^2*s^6)`, sharper than the earlier `O(y^2)` sufficient comparison.
The constants in these rescaled statements are uniform as the width and
kappa0 decrease through the stated ranges. Physical radial derivatives do
NOT have width-independent constants.

This is a theorem **relative to an exact analytic-core seed** with the
growth and regularity hypotheses below. It does not independently reprove
the nonlinear core contraction in the frozen axis attachment. Its fixtures
are not that core, and their nonzero baseline defect is explicitly retained.
Thus neither the entire inner foundation, the unconditional actual pulses,
the complete residual gain, nor blow-up is certified here.

Parent: `da0b3fe4a297f0367b75d095fac029103aa0153c`.
Defining equations: (4.7)--(4.11), (4.15)--(4.16), (4.20)--(4.23),
(B.26)--(B.31), and (7.1) of the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
The checked source SHA-256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Its theorem and formalization are not premises. Earlier evidence and GOAL.md
remain byte-for-byte unchanged.

## 1. Fixed seed, the actual activation, and the remaining premise

Keep the same fixed finite N,h and leading construction. Write
`A=1/2+h`, `D=1/2-h`, `d=1-eta^2`, `L=1-2h*eta^2`, `X=Xa*exp(y)`.
On the retained seed collar assume:

* The reference fields Er,Ur and their forward primitives are smooth,
  with Er>0. They agree with the actual inner core to the left of Xa.
  The common axis pressure datum is fixed, not refitted after activation.
* With `ar=1-2*partial_y log Er`, `qr=-2*partial_y Ur/Er`, the **exact**
  core equations give `ps,r=(ar,qr)`, hence `Tr=0`. On the closed collar,
  `ar>=p*>0`, `Fr=Er/sqrt(2X)>=F*>0`, and
  `2+gamma<=vr=ar+qr^2/ar<=V`, with gamma>0.
* All needed parameter derivatives are bounded. For the analytic statement,
  these fields and primitives, with each fixed radial derivative, are
  holomorphic and bounded on a common complex eta neighborhood; nonzero
  denominators stay nonzero on a smaller neighborhood.

Choose `0<t=t1<=tbar` so `y=t*s`, `0<=s<=1/4`, stays inside that collar,
and `0<=kappa0<=1/2`. Define the actual activation by forward integration:

\[
 \sigma(s)=\frac{e^{-1/s^2}}{e^{-1/s^2}+e^{-1/(1-s)^2}},\quad
 e_a=(1-\kappa_0)\sigma(s),\quad k=1-e_a,
\]
\[
 (\log E)_y=(1-k a_r)/2,\qquad U_y=-kE_rq_r/2.                 \tag{1}
\]

This is precisely the shear prescription: `a=k*ar` and
`-bs=k*qr*Er/E`. The common k cancels from `ts=-bs/a`; there is no
division by kappa0. Only the diagnostic widths use tbar=.01; the theorem
uses any fixed sufficiently short seed collar. The chosen enormous actual
amplitude and core are not numerically evaluated in this batch.

## 2. A flat Volterra operator with an honest endpoint

Put `g(s)=1/(exp(-1/s^2)+exp(-1/(1-s)^2))`, with `g(0)=e`.
For a smooth coefficient b and any nonnegative integer m define

\[
 {1\over\sigma(s)}\int_0^s\sigma(u)u^m b(u)\,du
       =s^{m+3}\mathcal V_m[b](s).                              \tag{2}
\]

The change of variable `v=1/u^2-1/s^2` gives the exact expression

\[
 \mathcal V_m[b](s)=\frac1{2g(s)}\int_0^\infty e^{-v}
 (1+s^2v)^{-(m+3)/2}
 g\left(\frac{s}{\sqrt{1+s^2v}}\right)
 b\left(\frac{s}{\sqrt{1+s^2v}}\right)\,dv.                    \tag{3}
\]

In particular `V_m[b](0)=b(0)/2`, not zero and not a numerical 0/0.
Each fixed s derivative of the integrand is bounded by a fixed polynomial
in v times exp(-v). The arguments lie in [0,1/4]; g and its reciprocal
have bounded derivatives there. Products and the chain rule therefore give
an integrable majorant for every fixed mixed derivative. Differentiation
under the integral proves a smooth extension at s=0.

The same proof works with b=b(t,s,eta,kappa0), uniformly on any compact
seed-parameter set including t=0 and kappa0=0. There are no derivatives of
a shrinking cutoff in (3). If b is holomorphic in eta on a common strip,
the integral is holomorphic on that strip, with uniform bounds on smaller
strips. These statements follow from (3), not a sampled Taylor series.

For b independent of t and eta, the first endpoint jets at m=0 are

\[
 \mathcal V_0[b](0)=b(0)/2,\quad
 \mathcal V_0[b]'(0)=b'(0)/2,\quad
 \mathcal V_0[b]''(0)=(b''(0)-3b(0))/2.                       \tag{4}
\]

Here g'(0)=2e contributes to the second identity's differentiated quotient
only by cancellation; in the second derivative the remaining g terms also
cancel. The coefficient of s^2 in the integral of a constant is -3/2;
this radial Jacobian must not be discarded.

## 3. Field differences and every forward moment

Write Delta for actual minus reference. Subtracting (1) from the seed
slopes, then using (2) with m=0, gives exact smooth factorizations

\[
 \Delta\log E=t e_a s^3 b_E,\qquad
 \Delta U=t e_a s^3 b_U.                                      \tag{5}
\]

For example `bE=V_0[ar(t*.,eta)/2]` and
`bU=V_0[(Er*qr)(t*.,eta)/2]`. The coefficient functions and every fixed
mixed derivative in the rescaled variables have finite uniform bounds as
t,kappa0 decrease. The entire function
`exprel(z)=(exp(z)-1)/z`, with value one at zero, gives

\[
 {E-E_r\over e_a}=t s^3 E_r b_E\,\operatorname{exprel}(t e_a s^3 b_E),
 \qquad {E_r/E-1\over e_a}=t s^3 b_R.                         \tag{6}
\]

Every difference of a smooth algebraic expression in the fields, excluding
radial derivatives, has the same factor. Since `dx=t*X*ds`, the forward
moments M,I,J,S and pressure increment Cp satisfy

\[
 (\Delta M,\Delta I,\Delta J,\Delta S,\Delta C_p)
            =t^2 e_a s^6 (b_M,b_I,b_J,b_S,b_C).               \tag{7}
\]

Indeed their normalized integrands, after removing `t*ea*s^3`, are
`X*bU`, `X*bH`, `X*(Ur*bH+Hr*bU+t*ea*s^3*bU*bH)`,
`X*(2Ur*bU+t*ea*s^3*bU^2-bE2/2)`, and `bE2/2`.
Here bH and bE2 are the corresponding coefficients of H=sqrt(2X)E and
E^2. A second use of (2), now m=3, proves (7). The pressure row has no
extra X: `Cp_y=E^2/2`. All eta derivatives commute with these formulas.

For clarity, use the original forward formulas rather than an abstract
continuity assertion. With NQ the numerator in Qs,

\[
 W=1-(2D\eta M+dM_\eta)/X,\quad
 N_Q=(1-h)I-D\eta I_\eta-dJ_\eta+2(h-D)\eta J,
\]
\[
 Q_s=-W+N_Q/(XH),\quad
 N_s=-WU+[D(M-\eta M_\eta)+4h\eta S-dS_\eta]/X
                  +4A\eta\Pi-d\Pi_\eta,                    \tag{8}
\]
\[
 p_1=XQ_s/L,\qquad p_2=XN_s/(LE),\qquad
 \Delta\Pi=\Delta C_p.
\]

Insert (5)--(7) and use the exact reciprocal-difference identity. This
proves

\[
 \Delta p_s=t e_a s^3 b_p.                                   \tag{9}
\]

For example `Delta(NQ/H)=Delta NQ/H-NQ,r*Delta H/(Hr*H)`.
The local denominator variation is order `t*ea*s^3`, even though the
moment variations have the smaller factor (7). Thus the integrated stress
coefficients cannot in general be assigned the stronger factor (7).
Both pressure contributions in (8) are retained; omitting `-d*Delta Pi_eta`
or `4A*eta*Delta Pi` changes the coefficient of (9).

## 4. Normalized stress and a quantitative cone collar

The exact zero-stress seed premise now gives

\[
 T_0=e_a B,\quad B=F\{(a_r,q_r)+t s^3 b_B\},\quad
 B(0,\eta)=F_r(0,\eta)(a_r(0,\eta),q_r(0,\eta)).              \tag{10}
\]

In particular `B-Fr*(ar,qr)` is divisible by `t*s^3` in the smooth
class. On a fixed-width y collar its first two radial jets vanish at y=0.
The direction matches the moving reference shear direction through second
order; it is not asserted constant away from the edge.

Set `ts=qr*(Er/E)/ar`, `vs=k*ar*(1+ts^2)` and use (4.20). Equations
(6), (9), and exact cancellation of k give smooth bounded coefficients

\[
 (P_c-v_s)/e_a=v_r+t s^3 D_c,\quad
 J_c/e_a=t s^3 J_c^*,\quad
 (v_s-kv_r)/e_a=t s^3 V_c.                                  \tag{11}
\]

These identities are relative, so they remain meaningful when the stress
is smaller than machine range. Suppose an actual seed bound gives
`|Dc|,|Jc*|,|Vc|<=K`. This K is obtained from the seed jets, lower bounds
for X,Er,ar,L, and (3), (5)--(9); it is independent of t,kappa0 in the
fixed compact ranges. We have not numerically certified such a K for the
actual enormous seed. Let `v*=2+gamma`, and put

\[
 d_* =\min\{1,v_*/2,v_*/(2\sqrt{\max(1,V-1)})\},
\]
\[
 s_* =\min\left\{\frac14,\left(\frac{d_*}{K\bar t}\right)^{1/3},
 [\max(16,16/9+\log(2(V+1)/\gamma))]^{-1/2}\right\}.          \tag{12}
\]

For K=0 omit the middle entry. On `0<s<=s*`,
`sigma(s)<=exp(16/9-1/s^2)`. Hence `K*t*s^3<=d*` and
`ea*(V+1)<=gamma/2`. Equations (11) imply

\[
 v_s\ge2+\gamma/2,\qquad (P_c-v_s)/e_a\ge v_*/2,
\]
\[
 0\le { (v_s-2)J_c^2\over2(P_c-v_s)^2}
    \le {2(V-1)K^2t^2s^6\over v_*^2}\le\tfrac12.             \tag{13}
\]

This supplies a strict admissible cone margin and its continuous edge
limit. The pulse quotient in (7.1) is at most `1/sqrt(2)` on this collar.
After also imposing sigma<=1/2 (already implied by s<=1/4),
`a>=p*/2`, `F>=F*`, `|g0|>=F*p*/2`. The correct growth formula is
`lambda0^2=2*a*F^2*(1-2/vs)`, so

\[
 \lambda_0^2\ge p_*F_*^2\gamma/(4+\gamma)>0.                 \tag{14}
\]

It is NOT `2*a*F^2*(vs-2)`; the denominator vs matters when the axial
shear is nonzero. Positive orientation follows from
`B dot (1,ts)=F*(Pc-vs)/ea>0`. B has a positive lower norm and smooth
unit direction on the closed smaller collar, including eta=+/-1.

Since `2<=g(s)<=6` and `1-kappa0>=1/2`, (10) yields a lower bound
`|T0|>=c*exp(-1/s^2)`. Each fixed mixed derivative has a bound

\[
 |\partial_y^j\partial_\eta^k T_0|
       \le C_{jk}t^{-j} e^{-1/s^2}s^{-3j}.                   \tag{15}
\]

The constants are uniform in rescaled t,kappa0 ranges under the seed
hypotheses. Once t is frozen they are finite profile constants, independent
of physical q, dyadic band, and correction stage. Equation (15) is not a
claim of uniform physical radial derivatives as t tends to zero.

The holomorphic seed neighborhood is preserved by (1)--(9): the cutoffs
have no eta dependence, forward integration preserves holomorphy, and
the exponential defining E never vanishes. On a smaller common strip the
remaining denominators stay nonzero. Fixed radial differentiation changes
bounds, not that parameter neighborhood. This verifies propagation of
the common collar, not existence of the analytic seed itself.

## 5. The obstruction: an approximately solved core is insufficient

For a reference with nonzero Tr the exact identity is instead

\[
 T_0=(F/F_r)T_r+e_a B,\qquad
 T_0/e_a=(F/F_r)(T_r/e_a)+B.                                 \tag{16}
\]

A fixed tolerance `|Tr|<=epsilon` cannot bound the first quotient near
s=0: `epsilon/ea` is unbounded for every epsilon>0. An exact core identity
(or a separately justified flat relative bound with the right direction)
is essential. Making physical q smaller does not correct this leading
profile error. The diagnostic seed has nonzero Tr, so tests of B are
explicitly tests of the activated excess, NOT tests of that seed's actual
stress cone. A `1e-100` baseline error already overwhelms the activation
at s=.02. No fixture is silently promoted to an actual core solution.

## 6. Checks, limits and the next item

`collar.py` implements the positive quotient, eta jets, all five moments,
the original stress map and its cancellation-safe difference. `review.py`
compares independent finite-interval quadrature, primitive derivatives,
literal moment-map subtraction at resolvable points, eta differentiation,
quadrature refinement, endpoint limits and deliberate omitted terms.
Tests also reject zero-assignment of the quotient, omitted Jacobian,
discarded pressure memory, baseline-defect erasure and evidence promotion.
The continuum argument is sections 1--5, not those sampled diagnostics.

The next foundation item remains independent verification of the nonlinear
analytic seed: its coefficient-space estimates, inverse, invariant ball,
exact residual cancellation and actual positive growth constants at the
same N,h. Compact modulation/repair, global moment identification, reserved
mean supports and the complete residual improvement remain open too.
This batch closes a relative activation/edge estimate, not those gates.

```sh
python3 research/navier_stokes_cascade/results/inner_collar_activation/test_audit.py -v
python3 research/navier_stokes_cascade/results/inner_collar_activation/audit.py \
  --verify-record research/navier_stokes_cascade/results/inner_collar_activation/evidence.json \
  --output inner_collar_activation_recheck.json
```
