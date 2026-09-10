# Post-pulse reference stress: finite interval, including angular patches

For the fixed reference family, an analytic reduction with outward constants
now bounds the stress on the **entire interval from pulse end through terminal
coordinate y=1/2**, for every eta in [-1,1]. The bounds are

\[
Q\ge h/2000>0,\qquad |w|<10^{-440},\qquad
a>2,\qquad 2-(a-2)w^2>1.99.
\]

Here `bs=0`. An explicit finite-radius check shows `XR>=100` suffices on this
new interval. The actual angular correction patches are included; their
pressure is not incorrectly set to zero inside their supports.

This is not a proof-assistant certificate or a verification of the complete
manuscript. **The bound does not extend to the entire terminal collar or
infinite reference tail.** At terminal coordinate 3, Q=0 and Pc=0, so the
reference tail fails the required Pc>2 test. The pending heat replacement
addresses that different endpoint regime. Axis attachment, subsequent PDE
corrections and smooth forcing remain unverified. No blowup candidate is
promoted; old Md4 and unforced `.14` failures are unchanged.

## Inputs and exact scope

Source: the supplied [OpenAI manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pages 26-31, 130-132 and 136-138: (4.9), (4.16), Lemma 4.5,
Proposition A.4, (A.10)-(A.18) and (A.31). The explicit endpoint `1/2` is on
printed page 131. PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent commit: `6e1104dcbd88f708e593bb7a5f61b52bd8a1dab6`.

Retain the preceding moment and pulse-stress arguments and their exact
reference schedule:

\[
M_d=64,\ T=e^{64}+10,\ \log P_*=T+1,\ \lambda=e^{-4T},\ h=\lambda^2,
\quad T_w=240T,\ T_u=120T,\ T_f=1000,\ c_o=.001.
\]

The numerical smallness reduction uses T>=64. This result is **not uniform
under arbitrary smaller independent h**: it uses h=lambda^2. The manuscript
uses a sharper exterior estimate for its more general parameter ordering;
our fixed family permits the coarser, explicit estimate below.

The inherited moment closure gives U=M=J=0 after the pulse, S(infinity)=0,
and the exact angular correction conditions (A.11). Therefore W=1 and bs=0.
At pulse end the preceding audit gives Q>=(lambda+eta^2)/4. All stage
boundaries are independent of eta and XR under logarithmic normalization.

## 1. Bounds on the actual corrected swirl

For either relative angular correction coefficient,
`||ci||C1<=7000lambda^30`, where `||g||C1=sup|g|+sup|g_eta|`. The disjoint
width-.3 bumps are translates of sigma'. Set

\[
\epsilon=448000\lambda^{30},\quad
\epsilon_\eta=\frac{\epsilon}{1-\epsilon},\quad
\epsilon_l=\frac{(20000/.3)7000\lambda^{30}}{1-\epsilon}.
\]

The parent bounds sigma' by 64 and |sigma''| by 20000. Thus the relative
edit and its eta derivative are <=epsilon, its logarithmic eta derivative
is <=epsilon_eta<1, and its slope change is <=epsilon_l<lambda/2.
Consequently E stays positive on the patches.

Throughout interpolation, the uniform interval, both angular patches and
**all later radii**, the exact slope and angular derivative obey

\[
-1\le l\le-3h/4<0,\qquad |g|:=|\partial_\eta\log E|\le1.
\]

To check each piece: interpolation has
`l=-lambda-sigma'(y/Tf) log(2/(1+eta^2))/Tf` and
`g=-theta*2eta/(1+eta^2)`; both bounds follow directly. The angular patches
have `-3lambda/2<=l<=-lambda/2`, and `lambda/2>=3h/4`. Exterior ramps
interpolate between -lambda, -1 and -h; the hold is -1. On the terminal
cutoff, `fo'/fo<=co*h*64/[2(1-co*h)]<h/4`. Beyond it l=-h. All pieces
after the angular patches have g=0.

These derivative bounds describe the **edited** E. Unlike before-support
pressure cancellation, they remain valid inside an angular patch.

## 2. A positive Q floor on the precisely stated interval

Since U=M=J=0, equation (4.9) reads

\[
Q_y+(1+l)Q=-l-h-D\eta g,\qquad D=1/2-h.
\]

During interpolation the last term equals `theta*c_eta>=0`; its source is
at least lambda-h. On a correction patch the source is at least
`lambda/2-h-epsilon_eta/2>=lambda/3`. Since `1+l<=1`, the barrier Q=lambda/4
has positive inward derivative. It therefore persists through interpolation,
the uniform interval and both patches, starting from the inherited pulse-end
floor. No angular samples are used for this argument.

At exterior-transition start, (A.11) gives I=XH/(1-lambda), independent of
eta. Substitution into (4.16), not a reset of Q, gives exactly

\[
Q_{\rm rel}=\frac{\lambda-h}{1-\lambda}=\lambda.
\]

The subsequent stages have the following lower bounds:

| Piece | Why Q stays positive |
|---|---|
| Unit ramp to l=-1 | Source >=lambda-h; Q>=lambda/4 remains a barrier. |
| Hold l=-1 of length 4 log(1/h) | Q_y=1-h, so Q increases by exactly 4(1-h)log(1/h). |
| Unit ramp to l=-h | Source >=0 and 0<=1+l<=1, so Q>=exp(-1)*Q_hold_end>1 throughout. |
| Waiting interval l=-h | Q=Q_before exp(-(1-h)y), decreasing exactly to Qp. |
| Terminal 0<=y<=1/2 | fo is still constant; Q=Qp exp(-(1-h)y). |

For rho=co*h, (A.13) simplifies without cancellation to

\[
Q_p=\frac{c_oh}{1-c_oh}\int_0^1
e^{(1-h)(1+2z)}\sigma'(z)\,dz.
\]

Since the integral of sigma' is one,
`co*h<=Qp<=co*h*exp(3)/(1-co*h)<1`. Thus the waiting length
`log(Q_before/Qp)/(1-h)` is positive and finite. On terminal [0,1/2],
`Q>=co*h*exp(-1/2)>h/2000`; every preceding floor is larger. This proves
the announced uniform Q floor on the covered interval.

## 3. Bound remaining energy and pressure without deleting patch terms

For any current post-pulse y and any t>=0, the **global future** slope bound
gives

\[
E(y+t)^2\le E(y)^2e^{-(1+3h/2)t}.
\]

By S(infinity)=0 and U=0 after pulse end,

\[
\frac{S(X)}X=\frac1{2X}\int_X^\infty E(x)^2dx,
\qquad \Pi(y)=-\frac12\int_y^\infty E(v)^2dv.
\]

Differentiating these actual integrals uses |g|<=1 and eta-independent
stage boundaries. Positive integration then yields

\[
0<S/X\le E^2/(3h),\quad |S_\eta|/X\le2E^2/(3h),
\quad |\Pi|\le E^2/2,\quad |\Pi_\eta|\le E^2.
\]

These bounds include all partial angular corrections, the waiting interval
and the infinite power tail. They do not require pressure cancellation
inside a patch. The factors 1/h are retained as finite quantities.

Formula (4.16) now reduces exactly to

\[
N=4h\eta S/X-dS_\eta/X+4A_p\eta\Pi-d\Pi_\eta,
\quad d=1-\eta^2,\quad A_p=1/2+h.
\]

Therefore

\[
|N|\le E^2\left[10/3+2h+2/(3h)\right]<2E^2/h,
\qquad |w|=|N|/(EQ)\le4000E/h^2.
\]

This estimate deliberately keeps finite h in the energy contribution. In
the pure power tail, `S/(XE^2)=1/(4h)`, so `4h*S/(XE^2)=1`, not zero.
Dropping that term because h is small would destroy an order-one cancellation.

## 4. Absorb the pulse-end decay with a finite logarithmic inequality

At pulse end,

\[
E_{\rm end}(\eta)\le
\exp(-119.5T+.3-13-6.5/\lambda)\le e^{-6.5/\lambda}.
\]

All later E decrease. For h=lambda^2 the inequality
`exp(-6.5/lambda)<=h^4` is equivalent to
`8lambda log(1/lambda)<=6.5`. The left side decreases for T>=64 and is
checked outward at lambda=exp(-256). No exponential of -1/lambda is
underflowed and then called zero. Consequently

\[
|w|\le4000h^2\le4000e^{-1024}
<7.665\cdot10^{-442}<10^{-440}.
\]

The tiny bound is a consequence of the very conservative fixed parameter
hierarchy, not a numerically observed small residual or a closeness measure
for a full blowup solution.

## 5. Check the actual finite-radius cone

Here bs=0, so `vs=a=2-2l`, `Pc=p=ps1` and `Jc=p*w`. The strict shear
excess satisfies `a-2>=3h/2>0`, while a<=4. Both ratio inequalities hold:
the first margin is a>2, and the second is `2-(a-2)w^2>1.99`.

Because X>=Xend and L<=1,

\[
p=\frac{XQ}L\ge\frac{X_{\rm end}h}{2000},\qquad
X_{\rm end}=X_Re^{241T+2+13/\lambda}.
\]

Even discarding the favorable `13/lambda` term, `XR>=100` gives
`log(p)>=log(100)+233T+2-log(2000)>log(100)`. Thus p>100>a, and

\[
\frac{2(P_c-a)^2-(a-2)J_c^2}{p^2}
\ge2(1-4/100)^2-2\cdot10^{-880}>1.84>0.
\]

Lemma 4.5 proves the strict admissible cone on the stated finite interval.
The audit stores `log(3h/2)` separately; subtracting rounded a and two
cannot establish the strict excess.

## Numerical controls and their limits

The diagnostic patch uses the actual preceding A.11 solver at `Md=4`,
lambda=1e-4, h=exp(-2T), Tf=1000, with 32 interpolation integration panels.
This is **not the chosen Md64 family** and remains a failed global candidate.
Both 32/48 quadrature orders use 280 digits. The coefficient eta derivative
comes from the exact differentiated linear-plus-quadratic moment system.

Partial backward pressure includes `2ci beta+ci^2 beta^2` on each remaining
piece of a bump. A separate DOP853 integration evolves its forward pressure
increment in coefficient-scaled units. Its complete increment is zero to
quadrature accuracy, while its partial increment is demonstrably nonzero.
This rejects the tempting but incorrect rule “pressure-preserving means
pressure unchanged inside the patch.” Quadratic terms appear explicitly in
both formulas. The recorded positive omission budget bounds their entire
pressure contribution, including additions that finite precision rounds
away. Even 280-digit arithmetic cannot resolve every endpoint quadratic
correction against its linear part. Exact closure follows from the inherited
analytic nonlinear moment argument, not a small numerical residual.

Exterior controls independently integrate both slope ramps from (4.9),
check the logarithmic waiting identity, and retain the exact attenuation
factors E->h^6 E and XE^2->h^8 XE^2 on the steep hold. Terminal controls use
the Md4 member of the **new** h=exp(-8T) family, at 280/320 digits and eight
integration panels. A backward
ODE for Q/h is compared with the positive integral (A.16). An expm1 formula
preserves the finite-h change in Qp/h that ordinary arithmetic loses.

At terminal y=3, Q=0 exactly. The pure power tail has N=0 by cancellation
of two order-one normalized terms, but Pc=0 still fails Pc>2. This endpoint
failure is a mandatory negative control, not a failure of the finite interval
claim. No ratio w=0/0 is evaluated there.

Quadrature refinements and ODE errors are diagnostics, not certified
quadrature error bars. The continuum argument rests on sections 1-5 and
the explicitly inherited exact moment and pulse arguments.

The initial 32/48-order comparison with 16 interpolation panels missed the
patch relative tolerance: the eta=1/2 coefficient changed by about
1.24e-19 against a 1e-20 limit. Four terminal panels also missed the 1e-30
limit, with a 1.12e-24 relative change in the stable finite-h difference.
The protocol therefore doubles those panel counts; neither tolerance is
relaxed. Arithmetic precision alone did not resolve those integration errors.

## Reproduce and code guide

```bash
python research/navier_stokes_cascade/results/post_pulse_stress_audit/test_audit.py -v
python research/navier_stokes_cascade/results/post_pulse_stress_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/post_pulse_stress_audit/evidence.json \
  --output post_pulse_stress_reproduced.json
```

`bounds.py` rejects unsupported endpoints and parameter ranges, constructs
outward universal constants at T=64, checks each bound above, and records
the actual parameter logarithms separately. Its returned positive error
quantities are never replaced by zero.

In `diagnostics.py`, `PatchSchedule` caches only the angle-independent input
of the inherited solver; its generic pressure interface deliberately refuses
use. `discrepancy_jet` computes both the angular moment discrepancy and its
eta derivative. `Patch.__init__` differentiates the quadratic correction
system. `remaining` integrates the exact remaining portion of each edit;
`density` normalizes before conversion to binary64; `independent` evolves
the separate forward equation. `record` keeps both paths and the positive
quadratic omission budget.

`terminal_q_over_h` evaluates the positive normalized terminal integral.
`terminal_control` compares direct subtraction with the stable finite-h
difference and checks the pure power cancellation. `independent_terminal`
integrates backwards from the exact zero endpoint. `exterior_control`
compares both unit ramps, the waiting identity and the two hold factors.

`audit.py` hashes the implementation, protocol, derivation and inherited
evidence; applies the fixed gates; and reproduces the recorded results.
`test_audit.py` checks individual identities, forbidden endpoint promotion,
tiny-term preservation and evidence tampering. None of these scripts changes
the historical evidence.

## Next construction gate

Build the heat exterior and its compensating moments, then verify the
endpoint collar and preservation of the already checked inner data and
cones. A scalar heat identity alone is insufficient: its moment changes
must be restored. Five-moment axis-annulus attachment and the complete PDE
residual/smooth-force construction remain separate open tasks.
