# Q111: closed all-order chart majorants for the current radial request

Parent: `e6f9a799e6e58b629734e189a774801829d1e202` (Q110).
On 2 November 2026 its research CI 41542142416 and field CI 41542142363
were still running; Q109 field CI 41533237315 had passed. The inspected
remote head was Q110. This checkpoint supplies explicit closed majorants
for ONE of Q110's missing inputs: the same-hull source and output chart
transfers. The old finite recurrence formulas are inherited from Q87 and
core_weighted_reverse_chart; their new all-order closed domination and
current-state binding are proved below. No actual high-jet source, raw
kernel profile or improved A0 tail is certified. Historical evidence,
root, fields, accepted cover and width C_setup^-20 are unchanged.

## 1. Keep the actual fixed-label power chart and positive hull

The original charts satisfy U=R^d, 1/3<d<5/2, with four other normalized
coordinates unchanged and the label exponent fixed under differentiation.
Let 0<a<=R<=b be the SAME current norm hull. Then the exact rational bounds

    u0=min(1,a^3) <= U <= max(1,b^3)                         (1)

hold: for a<1 and d<3, a^d>=a^3; for a>=1, U>=1; the upper bound is
analogous. They are coordinate enclosures, not new cutoff supports or a
replacement kernel gap. `ChartBridge.plan` retrieves a,b from Q110's
current physical profile and d from that graph's original compact chart,
checking the original range. It never accepts an unrelated convenient
hull or silently moves the source between state objects.

The source input from Q110 is a weighted factorial-through-R coefficient.
The original own-pole Eedge must be paid before the raw source estimate.
All pointwise raw-jet chart inequalities below may also retain the SAME
physical weight at corresponding points; no quotient/weight derivative
is being taken. Original root, pressure, compensator, local K0 and full
mean/covariance memory remain in Q110's source construction. A coordinate
bound does not cancel any of those fields or change their epsilon grade.

## 2. Forward chart: dominate Q87's finite recurrence at every order

For H=b and any valid same-hull lower coordinate floor u0>0, Q87's
sufficient R-through to U-through factor is

    F_N=243 N! sum_(j=0)^N vj,
    v0=1, vn=sum_(j=1)^n aj v_(n-j),
    aj=H binom(j+2,2)/u0^j+4*1_(j=1).                       (2)

The radial increments bound R(U)=U^(1/d) for 0<1/d<3; the four linear
increments are the unchanged coordinates. The physical chart and native
input hypotheses of Q87 remain. We prove a bound for the sufficient
factor (2), not an actual native field supremum.

Let A(t)=4t+H[(1-t/u0)^(-3)-1]. Positive coefficients and the recurrence
give sum vj t^j=1/(1-A(t)). Since binom(j+2,2)<=3^j for j>=1, choose

    tau=min(1,u0)/(16(1+H)),  KU=1/tau.

Then 4tau<=1/4, and

    H sum binom(j+2,2)(tau/u0)^j
       <=3H/[16(1+H)-3] <3/16.

Hence A(tau)<7/16<1/2. Because tau<=1 and every coefficient is positive,

    sum_(j<=N) vj <= tau^(-N) sum_(j>=0) vj tau^j
                    <2 tau^(-N),
    F_N <=486 N! KU^N,
    KU=16(1+H)/min(1,u0).                                  (3)

Convergence here concerns only the known coordinate majorants; no
analyticity of the flat source r0 is assumed.

This proof holds for every fixed N without computing a length-N recurrence
or materializing N!. `forward_recipe` stores the scalar, factorial argument,
geometric base and power only. It is deliberately conservative and does
not claim sharp factorial growth. At N=100005 the same proof applies;
the original actual source inequality through that order remains missing.

## 3. Reverse chart: dominate the existing exponential recurrence

The frozen raw-U to R chart uses B=max(1,b^2,1/a) and coordinate increment
bounds

    bj=4*1_(j=1)+3B/a^(j-1), j>=1,
    A_R(t)=4t+3B t/(1-t/a),
    sum ej t^j=exp(A_R(t)),
    G_N=N! sum_(j=0)^N ej.                                  (4)

The original proof bounds both normalized-through and raw-through R
output from the raw U input; its N! padding is retained. Take

    sigma=min(1,a)/(16(1+B)),  KR=1/sigma.

As above, 4sigma<=1/4 and
3B sigma/(1-sigma/a)<=3B/[16(1+B)-1]<3/16.
Thus A_R(sigma)<1/2 and exp(A_R(sigma))<2 (the exponential series is
strictly bounded by the geometric series at 1/2). Positivity yields

    G_N <=2 N! KR^N,
    KR=16(1+max(1,b^2,1/a))/min(1,a).                         (5)

This closes the OUTPUT chart factor through output n. It does not replace
the original raw e=0 kernel profile, source-pole unweighting, output edge
floor or physical derivative/scale transfers outside this fixed chart.

## 4. Insert both factors without losing a factorial

Q110's original A0 inequality uses raw-U source seminorm N_J(Qcirc r0).
Our forward chart result (3) has factorial-through-U OUTPUT. To bound the
raw source one must STILL use Q110's J! conversion. Consequently a safe
source chart/raw chain is

    N_J(raw U source) <= J! F_J * R_through_source_J
                 <=486 (J!)^2 KU^J * R_through_source_J,        (6)

after the original source unweighting at its exact common pole. Qcirc's
factor two stays in the raw tail recipe, not in a new informal cancellation.
The extra factorial is conservative; the existing J! in the inherited raw
recipe must not be confused with the J! already present in (3).

For the current request J=100005 the plan therefore records TWO factorial
instructions with that argument. The raw operator retains all original
powers, L, N_n(k0 chi_U), N_J(f0), temporal inverse floor, output edge-floor
inverse and original geometry. Attach output factor (5) at n after its
raw-U output estimate, keeping the inherited 6^n raw-to-through padding
as well. Any overlapping conservative factors are retained, not cancelled
without proof. No 100005!, its square, giant power or actual jet array is
numerically constructed here.

`ChartBridge.recipe` first calls Q110's guarded recipe: current source
archive, complete bundle lengths, high-jet reserves, full common pole and
frozen width must still match. It then attaches these chart factors from
the identical state's hull. The low-output archive still fails that gate.
`plan` can specify the chart factors without pretending the full high-jet
archive or raw profiles have been supplied. Neither method marks a radial
tail or complete residual as certified.

For fixed a,b these factors are independent of S/epsilon. A family of
physical scales has uniform chart constants only if its NORMALIZED hulls
have common positive lower and finite upper bounds. Without those, (3)
and (5) explicitly expose the losses as a approaches zero or b grows.
This proof is not a uniform actual-pulse theorem, does not normalize away
physical scaling losses, and does not alter the accepted band or width.

## 5. Verification and next obligation

Exact rational tests compare (3) against the original Q87 recurrence and
(5) against the original exponential recurrence through moderate orders
on several positive hulls. Independent polynomial differentiation checks
an actual power-chart composition. Exact radical inequalities check (1).
A typed current graph verifies state/hull binding, the deferred 100005
request and the two distinct factorials; a moderate full request retains
Q110's raw profile/width/pole/reserve gates. Materialization has an explicit
small diagnostic cutoff; the analytical estimates have no order cutoff.
Frozen replay checks unchanged Q110 and all source/evidence bindings.
Passing these finite checks is not a continuum proof of native inputs.

Next prove the ACTUAL current source-through-J bound, raw e=0 profile
seminorms and full common-pole/archive inputs identified by Q110, now with
the explicit chart factors (3)--(6), or find a sharper radial cancellation.
Only then can the radial component be joined to Q108/Q109 and the original
inner/outer cone estimates with their remaining transfers. Native
coefficient and reserved-mean support contracts, uniform actual-pulse
size/derivative/support on every physical scale, complete viscous/cutoff/
interaction residual improvement, infinite realization and blow-up remain
unproved.
