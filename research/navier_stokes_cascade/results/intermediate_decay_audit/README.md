# Intermediate decay audit: finite parameter choices and stress bounds

The next two reference stages now have an explicit parameter choice and a
bound over **every angle and their entire lengths**. Their dangerous stress
ratio must be smaller than one; the bound obtained here is **less than
`1e-42`**. This is an analytic inequality reduction evaluated with outward
interval arithmetic, rather than an inference from sampled trajectories.

Together with the preceding axial result, these three stages satisfy the
appropriate stress cone for any **`XR>=100`**. This is a sufficient radius
floor for these pieces. Later matching requirements may demand a larger
radius. The pulse, heat compensation, annulus and complete residual remain
unbuilt or unverified. This checkpoint does not prove the manuscript's
Navier–Stokes theorem or promote a blowup candidate. The analytic reduction
is documented for review, not certified by a proof assistant.

## What is covered

Set, in the order required by (A.6),

$$
M_d=64,\quad T=e^{64}+10,\quad \log P_*=T+1,\quad
\lambda=e^{-4T},\quad h=\lambda^2=e^{-8T},\quad T_w=240T.
$$

All parameters are finite and positive. These choices quantify the
requirements for the stages below; they do not assert every later
unspecified smallness condition in the paper.

| Reference stage | Bound | Cone scope at `XR>=100` |
|---|---|---|
| Axial transition, inherited from the previous audit | `abs(bs*w)<.526`; second ratio margin `>.949` | Strict where `bs!=0`; relaxed where `bs=0` |
| Unit transition from `l=0` to `l=-lambda` | `(-l)*w^2 < 1e-42` | Strict for `x>0`; relaxed at `x=0` |
| Constant `l=-lambda`, length `Tw` | `lambda*w^2 < 1e-42` | Strict throughout, including `eta=0` |

The numerical diagnostics use the computable **`Md=4` member of the new
parameter family**, with 280/320 digits. This checks the equations and
scale handling. It is not the `Md=64` construction and does not turn the
failed `Md=4` axial stage into a passing one. The continuum bounds at
`Md=64` follow from the inequalities below and the previous axial audit.

Four controls are retained in [evidence.json](evidence.json):

1. The old `lambda=1e-5`, `h=exp(-2T)` case still has
   `lambda*w^2=5.627489052...e27` at `eta=.5`, violating the necessary `<1`
   condition for every positive radius.
2. An initial version of this diagnostic omitted `exp(-T)` from `E^2`.
   At a small-angle checkpoint it reported about `1.92e-24` instead of
   `1.69e-52`, contradicting the proposed `1e-42` bound. The integrated
   energy-scale check caught the error. The corrected code retains the
   omitted-factor control.
3. A nonzero pole increment obeys `Delta v/(h*y)` approximately `-2`,
   where `v=N/(eta E^2)`. Dropping `h` makes that increment zero. An
   independently assembled source solution reproduces the nonzero value.
4. Binary64 makes the selected `lambda` and `h` zero. Wide-exponent
   arithmetic preserves them, but even 110-digit arithmetic loses `2lambda`
   when it is added to two. The strict shear excess is represented
   separately by `log(a-2)=log(2)-4T`.

The historical unforced `.14` grid-resolution failure is unchanged.

## Source and assumptions

Source: the user-supplied [OpenAI manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pages 26, 28–31 and 129–135; especially (4.9), (4.16), Lemma 4.5
and (A.23)–(A.27). PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The authoritative full hash is also pinned in `protocol.json`.
Parent checkpoint: `ceb373b15012fd4fcbd5e74aacb0c58d68ea9c28`.

The prefix retains the paper's fixed smooth step and shapes. To the right,
the unedited swirl has `l<=0`, angular exponent `0<=theta<=1`, and positive
amplitude. The later angular bumps may be omitted before their support
only because their total pressure increment is zero. Consequently the
backward pressure masses satisfy

$$
0\le B_\theta\le B_0\le1,\qquad
\Pi/E^2=-B_0/2,\qquad \Pi_\eta/E^2=J'_0B_\theta.
$$

These exact envelopes are used in the continuum proof. They do not require
the future axial pulse amplitude, since pressure depends on `E`, not `U`.
Future edits must preserve the required moments and bounds; their effect
is not automatically covered by this audit.

## 1. Positive angular stress at the next stage

Let `z=eta^2`, `delta=exp(-T)`, and use `Q=Qs`. The preceding audit gives
at the end of the axial stage

$$
Q_a\ge\frac{z+\delta}{10},\qquad
K_a\le4e^{-11},\qquad W\ge1-4e^{-11}>\tfrac12.
$$

The `Q` floor follows from the positive-Q identity there:
`Q0>=z/2+45 delta/(128 e)` and `Q>=(1-217*2e-28)Q0`.
The final eleven units have `k=0`, which supplies the bound on `K_a`.
Here `W=1-LK_a exp(-x)` and `0<L<=1`.

Use `x` for distance from the start of the unit ramp, and
`alpha=-l=lambda*sigma(x)`. Equation (4.9) becomes

$$
Q'+(1-\alpha)Q=\alpha W-h+c_\eta,
\qquad c_\eta=(1-2h)\frac{z}{1+z}\ge z/3.
$$

For `0<=x<=1`, its integrating-factor kernel lies between
`exp(-(x-t))` and one. Apply the lower kernel to positive terms and the
upper kernel to the negative `-h` term. This avoids assuming positivity
of `Q` before proving it. Since `h<=delta/(20e)`,

$$
Q(x)\ge\frac z{10}+\frac{\delta e^{-x}}{20}
+\frac\lambda2\int_0^x e^{-(x-t)}\sigma(t)\,dt.
$$

At `x=1`, symmetry gives `integral_0^1 sigma=1/2`, hence

$$
Q_b\ge z/10+\delta/(20e)+\lambda/(4e).
$$

On the constant-slope interval use its local coordinate `y`, with
`beta=1-lambda`. Since `h<=lambda/4`, the source is at least
`lambda/4+z/3`. Variation of constants and `e<3` give

$$
\boxed{Q(y)\ge z/10+\lambda/20+\delta e^{-\beta y}/60>0.}
$$

The memory term is retained, including at `eta=0`. Its small size is not
a reason to replace it by zero.

## 2. Bound the axial stress without dropping finite terms

Let `n=N/E^2`, `Jp=2eta/(1+z)`. The previous moment identity gives

$$
|n_a|\le3|\eta|(T+2).
$$

For completeness, its pressure and energy terms give
`|n_a/eta|<=2T+6+2h(T+2.5)+180delta`. The last two terms are less than
`T`, uniformly for `T>=64`. The continuous value is used at `eta=0`.
This uses `F1<1.5`, `E_a^(-2)<=4delta` and the bound 45 on the geometric
coefficient; all are retained from the axial derivation.

Here `U=0`, so the geometric terms in the axial **source** vanish. The
backward pressure envelopes imply `abs(Sn/E^2)<=3abs(eta)`. Put
`S(x)=integral_0^x alpha`. Then

$$
n'=S_n/E^2+2\alpha n,\qquad
|n(x)|\le3|\eta|(T+x+2)e^{2S(x)}.
$$

The exact energy scale is

$$
E^2=P_*^2 e^{-.4}\delta(1+z)^{-2}e^{-x-2S(x)}.
$$

The factor `delta` is the accumulated axial decay. The diagnostic checks
this scale independently by integrating `2(l-1/2)` from `E=P*f` through
the first ramp and axial interval.

## 3. Remove angle sampling from the stress test

For every `a>0` and `z>=0`,

$$
\frac{z}{(z+a)^2}\le\frac1{4a},
$$

because `(z-a)^2>=0`. Combining this inequality with the positive-Q
floors and the energy bound yields, on the unit ramp,

$$
w^2\le900P_*^2(T+x+2)^2e^{2S(x)}
\le900P_*^2(T+3)^2e^\lambda.
$$

On the power interval, `x=1+y`, `S=lambda(y+1/2)`, so

$$
w^2\le8100P_*^2(T+y+3)^2e^{-1+\lambda+\lambda y}.
$$

Dropping the favorable `exp(-1)` factor and using `y<=Tw=240T` gives a
single upper bound for both stages:

$$
\lambda w^2
\le8100e^2(241T+3)^2e^{-2T+e^{-4T}(1+240T)}
\le H(T):=8100e^3(241T+3)^2e^{-2T}.
$$

Both `exp(-4T)(1+240T)` and `H(T)` decrease for `T>=64`.
Outward interval arithmetic checks the derivative signs and evaluates

$$
H(64)<9.959505\cdot10^{-43}<10^{-42}.
$$

This argument covers the entire angular interval, including the moving
small-angle region where the memory term and angular term are comparable.
At `eta=0`, `N=0` and the stress ratio is exactly zero.
Since `bs=0`, `a=2+2alpha`, both expressions in (A.24) are positive:
`a>=2` and `2-(a-2)w^2>2(1-10^-42)`.

## 4. A sufficient finite radius on the three covered stages

Write `p=ps1=XQ/L`. Since `0<L<=1`, the preceding Q floors imply

$$
p_{\rm axial}\ge\frac{45X_R}{128}(1-217\bar\delta),\quad
p_{\rm ramp}\ge\frac{eX_R}{20},\quad
p_{\rm power}\ge\frac{e^2X_R}{60},\qquad \bar\delta=2\cdot10^{-28}.
$$

For example, on the power interval
`X=XR exp(T+2+y)` cancels `delta exp(-(1-lambda)y)` in its memory floor.
Thus no enormous radius needs to be materialized to check these bounds.
All three lower bounds exceed five when `XR>=100`.

On the two intermediate stages, `Pc=p`, `Jc=pw`, `vs=a<=2.2`.
For `alpha>0`, Lemma 4.5 reduces the cone to

$$
p>a,\qquad \alpha w^2<(1-a/p)^2.
$$

The right side is at least `(1-2.2/5)^2=.3136`, comfortably above
`1e-42`. At `alpha=0`, `vs=2` and `Pc>2` gives the relaxed condition.

For the axial stage let `C` be the prior upper bound on `abs(bs*w)` and
`D` its upper bound on `bs^2/2`. With `a=2`, the exact identity

$$
(v_s-2)(w+b_s/2)^2=(b_sw+b_s^2/2)^2/2
$$

shows it is sufficient that

$$
1-C-D/2-\frac{2+D}{p_{\rm axial}}>0.
$$

The interval evaluation in the record verifies this positive gap. Zero
shear uses the relaxed condition; nonzero shear uses the strict condition.
These bounds persist as `XR` increases. They do not verify the radius
conditions needed for the future axis attachment or the remaining stages.

## 5. Numerical diagnostics and a bounded future pressure approximation

The diagnostic uses `Md=4`, giving approximately `T=64.59815`,
`lambda=6.04671637e-113`, `h=3.65627788e-225` and `Tw=15503.556`.
It computes cumulative angular and energy moments, then evaluates (4.16).
The fixed-eta profiles retain the energy factor that the first draft missed.

For a tractable pressure crosscheck, continue the `-lambda` power beyond
the checked interval. The actual schedule retains that same slope and
angular shape for the following pulse of length `13/lambda`. At every
checked point, the difference between the actual pressure masses and the
continuing-power approximation is bounded by

$$
|\Delta B_0|,|\Delta B_\theta|\le e^{-13/\lambda}\le\lambda^8.
$$

The last inequality follows from
`8lambda log(1/lambda)<=8/e<13`. No doubly enormous exponential is
evaluated. The resulting absolute uncertainty in `v=N/(eta E^2)` is at
most `4lambda^8(1+x)exp(2S(x))`, retained as a positive number in the record.
The continuum proof above uses the exact pressure envelopes and does not
depend on this numerical approximation.

The independent ramp calculation solves the pressure-deficit equation
backwards with DOP853, then solves the Q and N source equations forwards.
It evolves `Z=(1-B)/(2lambda)` so the pressure deficit is not lost when
subtracting nearly equal numbers. Q is scaled by `eta^2+delta+lambda`.
Two step limits, seven angles and three ramp locations are checked.
On the full power interval an independently assembled scalar source
solution uses `expm1`; it is compared to the moment formulas at seven
locations, fixed angles, and an angle that follows the memory scale.
These methods share the prescribed prefix and bounded pressure approximation.
They are independent equation assemblies, not independent full-paper proofs.

At the pole, a cancellation makes `Delta v` proportional to `h*y`.
It is checked relatively on that tiny scale, using 280/320 digits.
The omitted-h and ordinary-arithmetic controls erase it. An absolute test
with a unit floor would miss this error, so the record comparison does not
use one for tiny increments.

## Reading the code and reproducing the result

| File | Calculation order |
|---|---|
| `bounds.py` | Check constants; compute logarithmic parameter intervals; bound both stress ratios; check the finite-radius cone. |
| `reference.py` | Build the first-ramp moments; carry them through the axial stage; propagate the two U=0 pieces; evaluate moments and independent source solutions. |
| `audit.py` | Refine interval precision; compare the diagnostic methods; retain all four controls; verify source hashes and scientific scope. |
| `test_audit.py` | Check algebra, scales, boundaries, small-angle behavior, the finite-h increment and evidence integrity. |

The numbered derivation above explains each mathematical calculation used
by the code. Function docstrings identify the state variables and every
change of normalization. The JSON stores rounded displays of interval
endpoints; the strict comparisons are made before serialization. The
claimed `1e-42` upper bound is an exact decimal bound.

From the repository root, with `requirements.txt` installed:

```bash
python3 research/navier_stokes_cascade/results/intermediate_decay_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/intermediate_decay_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/intermediate_decay_audit/evidence.json \
  --output intermediate_decay_reproduced.json
```

The audit recomputes 21 gates. CI also runs every preceding research audit.
The next construction task is the axial pulse's `M/J/S` closure and its own
stress-cone bounds, followed by the later exterior and annular matching work.
