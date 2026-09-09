# Pulse moment closure: finite bounds and stable correction diagnostics

The reference axial pulse now has a finite-parameter **moment-closure
argument**: two smooth corrections set `M=J=0`, and a unique smooth pulse
strength sets `S(infinity)=0`. Outward interval constants place that strength
between **1.0100502 and 1.0100504 for every eta in [-1,1]**. The normalized
finite-parameter remainder is below **1e-100**.

This is an analytic reduction checked with interval arithmetic, not a
proof-assistant certificate or a verification of the complete manuscript.
The **corrected pulse's stress cone remains unverified**. Heat compensation,
axis-annulus attachment, subsequent PDE corrections and smooth forcing also
remain open. No Navier-Stokes blowup candidate is promoted. The earlier
unforced `.14` failure and the failed `Md=4` axial stress condition remain.

## Scope and fixed source

Source: the user-supplied [OpenAI manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pages 28, 127 and 129-136, especially (4.15), (A.1)-(A.3) and
(A.14)-(A.20). PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent checkpoint: `7748c2652212096b58af60b101c2b7b98ee09580`.

Retain the preceding parameter family and now select the interpolation
length and terminal cutoff coefficient explicitly:

$$
M_d=64,\quad T=e^{64}+10,\quad \log P_*=T+1,\quad
\lambda=e^{-4T},\quad h=\lambda^2,\quad
T_w=240T,\quad T_f=1000,\quad c_o=.001.
$$

The proof below uses `T>=64`; it therefore also covers the moment equations
for the `Md=4` member. Only `Md=64` inherits the preceding all-angle axial
stress result. These statements must not be interchanged. The earlier
`Tf=200` numerical pressure pilot is retained unchanged; `Tf=1000` is a new
explicit choice, compatible with the earlier pressure-envelope bounds.

At pulse start let `Ep=eb*f`, `f=(1+eta^2)^(-1)`, and use the local logarithmic
coordinate `y`. The prescribed pulse is

$$
U=ER_b,\quad E=E_p e^{-(1/2+\lambda)y},\quad
R_b=A(\eta)R_0(\lambda y)+c_1\beta_1(y)+c_2\beta_2(y),
\quad 0\le y\le13/\lambda.
$$

Here `A` denotes **pulse strength**, not the paper's exponent `1/2+h`.
The flat step is (A.5),

$$
\phi(\xi)=\int_0^\xi\sigma(v/.02)\,dv,\qquad
R_0(\xi)=\phi(\xi)[1-\sigma(\xi-10)].
$$

The bumps are `sigma'((y-center)/.3+1/2)`, supported on intervals of width
`.3` centered at `13/lambda-3` and `13/lambda-1`. This is a permitted
rescaling of the paper's fixed bump. Each bump integrates to `.3`; the
coefficients absorb that normalization. Their supports are disjoint from
each other and from the main pulse, which ends at `11/lambda`.

## 1. Normalize the incoming moments before estimating them

By (4.15),

$$
M=\int U\,dX,\qquad J=\int UH\,dX,\qquad
S=\int(U^2-E^2/2)\,dX.
$$

Put

$$
m_0=\frac{M(X_p)}{X_pE_p},\quad
j_0=\frac{J(X_p)}{X_pH_pE_p},\quad
s_0=\frac{S(X_p)}{X_pE_p^2}.
$$

The exact accumulated energy scale is

$$
\log e_b=T/2+.3-T_w/2-\lambda(T_w+.5)
=-119.5T+.3-\lambda(T_w+.5).
$$

The factor `.3` follows by integrating the first slope ramp and the later
unit ramp; the entire preceding axial decay is retained. For the scalar
coefficients at the end of the axial stage, `0<=Ka,Ja<=4` and `0<=K2a<=16`.
Then

$$
m_0=\eta(1+\eta^2)K_a e^{-1-T_w}/e_b,
\qquad
j_0=\eta(1+\eta^2)J_a e^{-1-T_w+\lambda(T_w+.5)}/e_b.
$$

Use the norm `||g||_C1 = sup|g| + sup|g_eta|`; for vectors take the maximum
component norm. This is a Banach algebra norm. Since the norm of
`eta(1+eta^2)` is six, and `lambda(Tw+.5)<=1`,

$$
\|m_0\|_{C^1},\|j_0\|_{C^1}
\le24e^{-120.5T+.7}<1.
$$

The normalized positive `U^2` contribution to `s0` is at most
`64 exp(-T+.4)<1`; its eta derivative is at most
`256 exp(-T+.4)<1`. The normalized `E^2` memory is independent of eta.
The earlier first-ramp bound `F1<1.5` gives

$$
F_p\le e^{2\lambda(T_w+.5)}(T+T_w+2.5)
\le2(241T+2.5).
$$

Consequently

$$
|s_0|\le500T,\qquad |\partial_\eta s_0|\le1.
$$

`Prefix.moments` evaluates these same normalizations without their coarse
upper bounds. A regression check compares all three moment identities
against the preceding reference at matched quadrature. That identity check
does not claim convergence: the audit separately rejects four/eight incoming
panels and requires sixteen for its reported fixtures. The original
implementation and historical records are left untouched. The energy scale
also has a separate algebraic crosscheck.

## 2. Solve both axial moment equations with a nonzero determinant

The two weights in `dy` are `exp(s1*y)` and `exp(s2*y)`, where

$$
s_1=1/2-\lambda,\qquad s_2=1/2-2\lambda.
$$

Let `Y=13/lambda-3` be the first center and define
`b(s)=integral exp(s*t) beta(t) dt` for its centered bump. Normalize each
moment at `Y`. The exact system is

$$
B c=d(A,\eta),\qquad
B=\begin{pmatrix}b(s_1)&e^{2s_1}b(s_1)\\
b(s_2)&e^{2s_2}b(s_2)\end{pmatrix},
$$

$$
d_i=-m_{i,0}e^{-s_iY}
-A\int_0^{11/\lambda}e^{s_i(y-Y)}R_0(\lambda y)\,dy,
\qquad(m_{1,0},m_{2,0})=(m_0,j_0).
$$

The determinant factor is exactly

$$
e^{2s_2}-e^{2s_1}=e^{2s_1}\operatorname{expm1}(-2\lambda)<0.
$$

For `lambda<=.01`, positive bump masses satisfy
`.3 exp(-.075)<=b(si)<=.3 exp(.075)`, and
`1-exp(-2lambda)>=lambda`. Interval arithmetic gives

$$
\|B^{-1}\|_\infty<20/\lambda.
$$

The full main pulse lies at least `2/lambda-3` before the first center.
Use `R0<=11` and `si>.4` to obtain, uniformly for `.9<=A<=1.2`,

$$
\|d\|_{C^1_\eta},\ |\partial_A d|
\le200e^{-.8/\lambda},\qquad
\|c\|_{C^1_\eta},\ |\partial_A c|
\le4000\lambda^{-1}e^{-.8/\lambda}\le\lambda^{20}.
$$

The last comparison is evaluated without forming `exp(-1/lambda)`:
it is equivalent to
`lambda*(log(4000)+21*log(1/lambda))<=.8`, whose maximum on the chosen
range occurs at `lambda=exp(-256)`. Thus `M=J=0` exactly for each trial
strength. The coefficients are affine in `A` and smooth in eta.

The displayed C1 bound holds at fixed trial strength. After substituting
the root below, the chain rule gives the safe bound `2 lambda^20` for the
actual coefficient functions. This distinction matters for a future cone
audit.

## 3. Include the later angular corrections in the S budget

The `E` corrections in (A.11) are independent of the axial strength.
We establish their existence at the new parameters instead of treating the
old pressure pilot as a certificate for them.

Coarse global derivative bounds for the fixed step are
`|sigma'|<=64`, `|sigma''|<=20000`. On `[1/4,3/4]` these follow from
`sigma(1-sigma)<=1/4`,
`2/t^3+2/(1-t)^3<=256` and `|(-6/t^4+6/(1-t)^4)|<=3072`.
For `0<t<=1/4`, use `sigma(1-sigma)<=exp(4-1/t^2)`. The resulting
majorants for the two derivatives are

$$
(2t^{-3}+16)e^{4-t^{-2}},\qquad
(4t^{-6}+6t^{-4}+64t^{-3}+352)e^{4-t^{-2}}.
$$

Each summand increases up to `t=1/4`; evaluating there bounds the first by
one and the second by one. Reflection covers the other endpoint.
These deliberately loose bounds suffice at the selected parameters.

Write `beta=1-lambda`, `a=log(2/(1+eta^2))`. At interpolation end, the
deviation of `I/(XH)` from `1/beta` is exactly

$$
\delta_I=e^{-\beta T_f+a}\delta_0+
\frac a\beta\int_0^1
e^{-\beta T_f(1-z)+a\sigma(1-z)}\sigma'(z)\,dz.
$$

Here `delta0` is independent of eta, `|delta0|<=1/beta`,
`0<=a<=log(2)` and `|a_eta|<=1`. Hence
`|deltaI|<=4/beta`, `|deltaI_eta|<=6/beta`. At the first angular bump,
another `Tu-3` units have passed, with `Tu=30 log(1/lambda)=120T`, so

$$
\|d_I\|_{C^1}\le(10/\beta)e^{-\beta(T_u-3)}
\le700\lambda^{30}.
$$

The normalized angular matrix has slopes `+beta` and `-(1+2lambda)`;
the pressure row has an extra factor two. Direct outward bounds give
inverse norm below five. The quadratic pressure term has norm below 32,
using the positive bump integral and `sigma'^2<=64 sigma'`.
Lemma A.2 therefore applies because

$$
8\cdot5^2\cdot32\cdot700\lambda^{30}<1.
$$

It supplies exact (A.11) corrections with C1 norm at most `7000 lambda^30`.
Their relative edit and eta derivative are below `.01`. Their radial
slope change is bounded by

$$
\frac{(20000/.3)7000\lambda^{30}}
{1-64\cdot7000\lambda^{30}}<\lambda/2.
$$

Thus `E>0` and `l<=-lambda/2` persist on the edited patches.
The choice `Tf=1000` gives `l>=-lambda-.1` on interpolation using
`64 log(2)/Tf<.1`. With `co=.001`, the terminal slope is at most
`-3h/4`, including its small cutoff edit.

## 4. Bound every post-pulse energy contribution

Let

$$
R(\eta)=\frac1{X_{\rm end}E_{\rm end}(\eta)^2}
\int_{X_{\rm end}}^\infty E^2\,dX.
$$

The interpolation contributes at most `Tf`. The following uniform interval
contributes at most `2Tu`, allowing the relative angular edits. For the
release, the first unit interval contributes at most one, the `l=-1` hold
at most one-half, and its remaining tail at most `(2/3)h^7`. The last bound
uses the exact attenuation `XE^2 -> h^8 XE^2` and integrates the subsequent
rate `l<=-3h/4`. Thus the entire release contributes less than two.

$$
0<R\le1002+240T\le256T,\qquad
|R_\eta|\le4\cdot256T.
$$

The derivative bound follows from the squared normalized angular shape:
its logarithmic derivative is at most two, plus
`2*.01/(1-.01)` on the edited patches. Four is a safe common bound.
Stage endpoints do not depend on eta.

The prescribed terminal waiting length is finite and positive: the `l=-1`
hold and following unit transition leave `Q>1`, whereas
`Qp<=co*h*exp(3)/(1-co*h)<1`. With the exact angular correction imposed,
the scalar `Q` identity then gives the terminal angular moment in (A.8).
This argument concerns moments; it does not verify the stress cone there.

## 5. The last moment is a scalar quadratic with a unique root

The disjoint supports eliminate cross terms between the main pulse and its
bumps, and between the two bumps. Therefore, exactly,

$$
F(A,\eta):=\frac{\lambda S(\infty)}{X_pE_p^2}
=A^2K_b-\frac{1-e^{-26}}4+\mathcal E(A,\eta),
$$

$$
\mathcal E=\lambda s_0-\frac\lambda2e^{-26}R
+\lambda\sum_{i=1}^2 c_i(A,\eta)^2
\int e^{-2\lambda y}\beta_i(y)^2\,dy.
$$

Each last integral is below `.3*64<20`. Since `ci` is affine in `A`, the
whole expression is a quadratic, with tiny but nonzero coefficient
increments. In particular,

$$
|\mathcal E|,|\partial_\eta\mathcal E|\le1000T\lambda,
\qquad |\partial_A\mathcal E|\le80\lambda^{41}.
$$

The function `T exp(-4T)` decreases for `T>=64`. Outward arithmetic yields

$$
1000\cdot64e^{-256}<4.234408\cdot10^{-107}<10^{-100}.
$$

To enclose `Kb`, `bounds.py` uses monotone left/right sums for the initial
step primitive, an exact middle integral on `[.02,10]`, and monotone bounds
for the cutoff on `[10,11]`. At 1024 boxes,

$$
.2450496182718<K_b<.2450496221363.
$$

The displayed decimals are rounded outward; comparisons use the full
interval endpoints. These bounds and the finite remainder give

$$
F(1.0100502,\eta)<-3.08\cdot10^{-8},\qquad
F(1.0100504,\eta)>6.42\cdot10^{-8}.
$$

On the entire paper bracket `[.9,1.2]`, `F_A>.441>.35`. The intermediate
value theorem and strict monotonicity give a unique root in the narrower
interval, at every eta. The implicit function theorem gives smoothness,
including the one-sided endpoint derivatives, and

$$
|A_\eta|\le3000T\lambda.
$$

Only C1 bounds are claimed quantitatively. All fixed higher derivatives
exist after the parameters are fixed; the smallness bounds needed by later
PDE constructions have not been established here.

## Numerical diagnostics and deliberate failures

The interval argument is the continuum result. Numerical refinement is a
separate implementation check; it is not substituted for interval bounds.

| Diagnostic | What it establishes |
|---|---|
| `Md=4`, lambda `1e-4` and `1e-5`, h `exp(-2T)` | Actual scaled axial corrections, the full post-pulse swirl, its angular corrections and the S quadratic. These fixtures are not globally admissible candidates. |
| 32/48 quadrature orders and 100/140 digits | Refinement of weighted pulse integrals, correction coefficients and numerical roots. |
| Independent DOP853 integration on the end patch | Both moment ODEs reach zero to the declared numerical tolerance. The incoming main-pulse moments are shared inputs. |
| New-family `Md=4`, 280/320 digits | Retains the positive finite amplitude increment; it is not the `Md=64` numerical profile. |
| Positive stage-by-stage tail records | The infinite tail and tiny signed angular energy edits are accounted for separately. |

The flat cutoff creates a narrow maximum in the weighted pulse integrals.
`log_main_moment` locates it and integrates in its natural width. The middle
part is integrated exactly. The unevaluated `[0,.02]` part has a separately
recorded positive relative upper bound below `1e-100` for the fixtures.
The quadrature error itself is controlled numerically by refinement, not
rigorously enclosed. Accordingly, fixture moment residuals are numerical
diagnostics; exact closure belongs to the analytic argument above. The very
small linear residuals measure the quadrature-assembled system, not certified
errors against the exact weighted integrals.

At lambda `1e-4`, eta `.5`, the numerical strength is approximately
`1.0745601181`, not the leading value `1.0100502663`. The difference is
mostly the incoming S-memory. At lambda `.001`, even the upper estimate of
`F(1.2,0)` is negative: the paper's amplitude bracket fails for that
deliberately oversized parameter.

For the new-family `Md=4` diagnostic, lambda is approximately `6.047e-113`.
The finite amplitude increment is approximately `9.510e-109`. We compute it
both by direct high-precision subtraction and by rationalization. Ordinary
binary64 subtraction returns zero. This calculation uses an **exact-integral
proxy** that omits the two axial bumps' positive energy; their effect would
lower the strength by at most `120 lambda^41`. That positive omission bound
does not include numerical quadrature error and must not be represented as
a complete error bar for the printed root. The continuum amplitude interval
above includes all those contributions.

Additional retained controls:

- Four and eight panels under-resolve the incoming axial integrals. Their
  32/48-order discrepancies are about `1.13e-13` and `1.89e-18` for Ka,
  exceeding the fixed `1e-20` gate. Sixteen panels reduce that discrepancy
  below `3e-23`; K2a also passes, below `8e-22`. The earlier continuum
  stress bounds use analytic bounds on these integrals and are unaffected.
- Removing the axial end corrections leaves the moment discrepancy.
- Correcting only M leaves a large relative J-error, above `.8` in the fixtures.
- Subtracting the two exponential matrix entries in binary64 gives zero at
  lambda `1e-20`; `expm1` retains their negative difference.
- The affine dependence on the incoming moments and its nonzero quadratic
  energy contributions are stored even when ordinary arithmetic erases them.
- The original `lambda*w^2>1` failure from the earlier checkpoint is retained.

The axial prefix is odd in eta, while the fixed swirl is even. Tiny affine
cross terms mean the fully corrected strength is **not assumed even**.

## Code guide and reproduction

The suite contains 32 regression tests and 34 audit gates. The private
energy-only schedule skips unused pressure integrals and caches the common
incoming angular moment. Post-interpolation energy integrals are independent
of eta before normalization and are reused across angles. A regression check
compares this optimized calculation with the original schedule and its
previously evaluated total; the private schedule rejects pressure requests.

| File or function | Calculation and reason |
|---|---|
| `bounds.shape_enclosure` | Integrates positive interval bounds for Kb. Refining the boxes narrows the enclosure. |
| `bounds.closure_bound` | Checks derivative, inverse, contraction, energy and root-sign inequalities before emitting a closure status. |
| `Prefix` | Accumulates the initial moments and the full logarithmic energy scale. |
| `Continuation.energy_mass` | Integrates all eight post-pulse stages and records each angular edit separately. |
| `stable_two_row` | Divides out bump masses, uses the nonzero `expm1` determinant, then solves the two rows. |
| `MomentSolver.affine` | Separates the coefficient multiplying strength from the inherited-moment coefficient. |
| `MomentSolver.root` | Forms the S quadratic, evaluates its positive root, and measures both linear and quadratic residuals. |
| `MomentSolver.source_check` | Integrates `m_i' + s_i*m_i = Rb` over the end patch using a separate bump/ODE implementation. |
| `family_increment` | Keeps the tiny strength shift and a positive bound for omitted axial-bump energy. |
| `audit.py` | Refines calculations, checks failures, pins source hashes, and verifies a saved evidence record. |
| `test_audit.py` | Tests equations, scales, nonzero small terms, failure detection and evidence integrity. |

Function comments explain each normalization and calculation. No decimal
size of a small residual is interpreted as a distance to a proof.

From the repository root, with the listed requirements installed:

```bash
python3 research/navier_stokes_cascade/results/pulse_moment_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/pulse_moment_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/pulse_moment_audit/evidence.json \
  --output pulse_moment_reproduced.json
```

The first command runs the regression checks. The second reconstructs the
evidence, checks its protocol and source hashes, compares the new values to
the saved record, and writes a separate output. CI runs these alongside the
preceding research checks.

The next required research gate is the **stress cone throughout the corrected
pulse**, using the now-bounded strength and its eta derivative. Moment closure
alone does not establish that gate or the manuscript's blowup theorem.
