# Axial stress audit: a bounded repair and the next parameter obstruction

This checkpoint reproduces the failed `Md=4` axial stress test and obtains
positive bounds for both ratios in (A.24) with **`Md=64` over the entire axial
interval and every `eta` in `[-1,1]`**. The calculation combines an analytic
inequality reduction with outward-rounded `mpmath.iv` arithmetic. It does not
sample angles to establish the bound, set `h=0`, or underflow away small terms.

This is a local result about one stage of the manuscript's reference profile.
It does not verify the manuscript's theorem, produce a complete matched
profile, or promote a blowup candidate. The analytic reduction below remains
subject to mathematical review; this is not a proof-assistant certificate.
The original unforced `.14` resolution failure and all previous evidence
remain unchanged.

## Results and scope

The [record](evidence.json) retains:

- The `Md=4`, `eta=.5`, `y=exp(2)-1` failure: `Q=1.060674272...`,
  `N/E^2=-5.732467126...`, and `Pc/ps1=-1.925705214...`. Since `Q>0`,
  `Pc` is negative for every positive `XR`. A separate integration of (4.9)
  reproduces the moment evaluation. Increasing `XR` cannot repair this sign.
- Outward interval bounds for `Md=4,8,16,32,64,128,256`. `64` is the first
  passing choice in this list, not an optimized threshold. An inconclusive
  upper bound is not a rejection of the corresponding parameter.
- Nested partitions with 256, 512 and 1024 boxes in transition coordinate
  `t=log(1+y)/Md`, and 30/50-digit interval refinement. Even the coarse
  `Md=64` bound gives `|bs*w|<.555`, hence a second (A.24) margin above `.89`.
  Exact reported bounds and refinements are in the record.
- A second failure in the old protocol: at the start of its constant
  `l=-lambda` interval, with `eta=.5` and `lambda=1e-5`,
  `lambda*w^2=5.627489052...e27`, whereas the necessary condition is `<1`.
  This failure also cannot be fixed by increasing `XR`.

The new axial bound implies that a sufficiently large radius can realize
the axial cone inequalities via Lemma 4.5. **No finite `XR` is selected here**,
and no later outer stage is certified. At zero shear, including `eta=0` and
the constant endpoints, `vs=2`; these locations use the relaxed condition.
We do not claim the strict `vs>2` condition there.

The next task is to choose and check `lambda` *after* the new `Md` and `P*`,
and then choose `h` sufficiently small relative to that `lambda`.
The old numerical value cannot be justified by the axial repair. After that
come the pulse's `M/J/S` closure, heat compensation and the five-moment axis
annulus. The new enormous scale also requires a fresh inner-amplitude
precision strategy; the preceding pilot's 120/160-digit `log(g)` check is
not transferred to `Md=64`.

## Source and finite parameters

Source: the user-supplied [OpenAI manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pages 26–31 and 129–135; equations (4.9), (4.11), (4.16),
(4.20)–(4.21), Lemma 4.5 and (A.5)–(A.7), (A.23)–(A.26).
The PDF SHA-256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The parent checkpoint is `1b3fda937cf89ff6fb9321748e26c302a748256a`.
The preceding schedule, protocol and evidence are hashed by this audit.

For every candidate the bound covers the following finite parameter family:

\[
T_d=e^{M_d}+10,\qquad \log P_*=T_d+1,\qquad 0<h\le e^{-2T_d}.
\]

Assume `h<lambda<.1` and the preceding swirl shapes (`Tf=200`, `co=.001`).
Choose `h` after `lambda`, as prescribed in (A.6). The historical numerical
controls alone retain the equality `h=exp(-2Td)`; that equality is not imposed
on a future globally matched profile. All estimates below hold uniformly
for smaller positive `h`, so making that later choice does not lose the
axial bounds.
No unspecified later-stage smallness threshold is asserted. The unedited
swirl to the right of the axial stage has `l<=0`, angular exponent
`0<=theta<=1`, and positive amplitude. The terminal correction also has
negative `l`: `sigma'<=8` gives `fo'/fo<=4 co h/(1-co h)<h`.
The pressure-preserving angular bumps can be omitted from the backward
pressure integral before their support, by (A.11)/(A.23). These facts,
rather than a numerical pressure extrapolation, supply the tail bounds below.

For all `Md>=4`, outward arithmetic verifies

\[
\delta=e^{-T_d}<2\cdot10^{-28}=:\bar\delta,
\qquad 0<h\le\delta^2,
\qquad E^{-2}\le4e^{-T_d-1.6}<4\bar\delta
\]

throughout the axial stage. Here
`E=P* exp(-.2) (1+eta^2)^(-1) exp(-y/2)`.
These positive upper bounds are carried through the calculation. Neither
`P*`, `h` nor the inner amplitude is materialized for the larger choices.

## Reduction over every angle

Write `e=eta`, `z=e^2`, `d=1-z`, `Jp=2e/(1+z)`, and let

\[
K=4e^{-y}+\int_0^y e^{s-y}k(s)\,ds,\qquad
K_2=16e^{-y}+\int_0^y e^{s-y}k(s)^2\,ds.
\]

The monotone profile has `0<=k<=K<=4`, `0<=K2<=16`. Define the first-ramp
constants `r1=I/(XH)` and `F1=(integral E^2 dX)/(XE^2)` at its end.
In the ideal region they start at `5/8` and `5/6`. The first ramp satisfies
`r'=1-(1+l)r`, `0<=l<=.6`; consequently

\[
5/8\le r_1<1,\qquad 3/(8\mathrm e)\le c:=1-r_1\le3/8.
\]

Indeed `r>=5/8` is an invariant lower barrier, and
`(1-r)'=-(1-r)+lr>=-(1-r)`. Also

\[
F_1=e^{-.6}\left(5/6+\int_0^1 e^{1.2(t-\int_0^t\sigma)}dt\right)
\le1+(5/6)e^{-.6}<1.5.
\]

The last bound uses `0<=t-integral_0^t sigma<=.5`.
No quadrature approximation is needed for these inequalities.

Substitution in the exact moment expression (4.16) gives

\[
Q_0=\frac{z[1+2dK]+c(3-6z+8z^2)e^{-y}}{1+z},
\]
\[
Q-Q_0=h\left[-\frac{1+3z}{1+z}+2zK
-ce^{-y}\left(-\frac{1+3z}{1+z}+16z\right)\right].
\]

Here `Q0` denotes the expression at `h=0`, not an approximation substituted
for `Q`. Since `3-6z+8z^2>=15/8`,

\[
Q_0\ge\frac{45}{128\mathrm e}e^{-y},\quad
|Q-Q_0|\le28h,\quad
Q\ge(1-217\bar\delta)Q_0>0.
\]

This retains positivity at `eta=0`, where the memory term is essential.
It avoids catastrophic cancellation in the original moment expression.

For pressure, define dimensionless positive tail masses

\[
B_0=\int_y^\infty E(s,e)^2/E(y,e)^2\,ds,\qquad
B_\theta=\int_y^\infty\theta(s)E(s,e)^2/E(y,e)^2\,ds.
\]

Both are between zero and one. Thus `Pi/E^2=-B0/2` and
`Pi_eta/E^2=Jp Btheta`. The exact normalized axial moment is

\[
n:=N/E^2
=\frac{e[-Wk+(4hz-2d)K_2]}{E^2}
-(dJ_p+2he)(F_1+y)-(1+2h)eB_0-dJ_pB_\theta,
\quad W=1-(1-2hz)K.
\]

In particular the absolute geometric coefficient is at most
`44+64h<45`. This identity contains all geometric, pressure and finite-`h`
terms; `M-eta M_eta=0` for this linear-in-eta prefix.

Set `r=1/(1+y)=exp(-Md t)`, where `t=log(1+y)/Md`. For `z>0`, combine the
last identities, cancel the common factor `z` only after using a positive
lower bound for `Q`, and use `K>=k`. The remaining angular expression is a
ratio of affine functions of `d`:

\[
\frac{2d(1+F_1r)+(2-d)r}{1+2dk}.
\]

Its derivative has constant sign, so its maximum over `0<=d<=1` occurs
at an endpoint. Using `F1<1.5` gives

\[
\boxed{\quad
|b_sw|\le
\frac{8\sigma'(t)}{M_d(1-217\bar\delta)}
\left[\max\left(2r,\frac{2+4r}{1+2k}\right)+370\bar\delta\right].\quad}
\]

The finite remainder is bounded by `10h+360 delta<=370 bardelta`:
the first term comes from finite `h`, the second from the geometric term
and `E^(-2)<=4 delta`. At `z=0`, `N=bs=0` and the product is exactly zero;
the displayed bound remains valid by the separate endpoint calculation.
It therefore also covers every small nonzero angle and its narrow memory
transition, without sampling that transition.

Finally `|k'|<=32/Md` gives

\[
b_s^2\le\frac{16384}{M_d^2}\bar\delta\le1024\bar\delta.
\]

Since `a=2`, the two lower margins checked are
`2-|bs*w|` and `2-2|bs*w|-bs^2/2`. They bound both expressions in (A.24)
from below. The corresponding lower bound on `Pc/ps1` is `1-|bs*w|/2`.
For `exp(Md)-1<=y<=Td`, `k'=0`, so both ratio inequalities are immediate.

## Outward interval calculation in the remaining coordinate

The closed-form step satisfies `sigma'<=8`. For completeness, put
`v=|2t-1|`, `u=v^2`. Then

\[
\sigma'(t)=\frac{8(1+3v^2)}{(1-v^2)^3}
\operatorname{sech}^2\!\left(\frac{8v}{(1-v^2)^2}\right).
\]

The inequality `cosh^2 x>=1+x^2` reduces the proposed bound to the
nonnegative polynomial
`58u+9u^2-4u^3+u^4`; it is nonnegative on `[0,1]` because
`9u^2-4u^3>=5u^2`. Equality at the midpoint gives the maximum eight.

On each exact rational box `[a,b]`, `bounds.py` uses:

- `r<=exp(-Md a)` and `k>=4 sigma(1-b)`;
- the convexity of `2/t^3+2/(1-t)^3` to bound its maximum by an endpoint;
- the monotonicity of `sigma(t)(1-sigma(t))` toward `1/2` to bound its
  maximum at the closest point to `1/2`;
- the smaller of that derivative enclosure and eight.

For boxes contained in `[0,1/4]`, the endpoint-safe majorant
`(2/t^3+16) exp(4-1/t^2)` is increasing and bounds `sigma'`; the right
endpoint follows by reflection. Each arithmetic/transcendental expression
is evaluated with outward interval rounding, and the upper endpoints feed
the bound above. Thus there is no interpolation assumption between boxes.
The small differences between interval precision settings concern enclosure
tightness; positivity is checked independently in each run.
The JSON decimals display rounded interval endpoints; the runner checks
the margins using the interval endpoints before serialization.

## Independent numerical controls and remaining failure

`diagnostics.py` evaluates (4.16) with 80-digit cumulative moments for the
old `Md=4` protocol. A separate double-precision DOP853 implementation
integrates (4.9), evolving
`log(E/E0), Pi/E0^2, Pi_eta/E0^2, K, Q, N/E^2` through the first ramp and
the axial stage. It uses an independently coded step and the source terms
directly, with two maximum step sizes, five angles including zero and
`1e-12`, and four axial locations. The pressure datum and ideal initial
moments are shared. This is independence of equation assembly, not an
independent construction of the entire paper. Long double-precision
pressure cancellation is explicitly excluded (`y<=10`).

The old constant-slope failure uses backward pressure to avoid that
cancellation. At its start `U=bs=0`, so `a=vs=2+2 lambda`, `Pc=ps1`
and `Jc=ps1 w`. If the cone held, `ps1>a` and

\[
2\lambda w^2<2(1-a/p_{s,1})^2<2.
\]

The measured `lambda*w^2>5e27` violates this necessary condition for all
positive radii. A 32/48-order quadrature refinement checks this value and
the axial counterexample. This retained failure explains why the next
checkpoint must retune the later stage before attempting moment closure.

## Reproduction

From the repository root, with `requirements.txt` installed:

```bash
python3 research/navier_stokes_cascade/results/axial_stress_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/axial_stress_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/axial_stress_audit/evidence.json \
  --output axial_stress_reproduced.json
```

The runner recomputes 17 gates, verifies source hashes, retains both
counterexamples and rejects scientific scope promotion. CI runs the same
commands alongside the preceding research checks. Reproduction tolerances
allow platform-dependent ODE roundoff; the interval margin gates themselves
are always recomputed and must be positive.
