# Actual axial-taper cutoff crossings and a necessary shape budget

Parent: `8c34b682cd1d1ee593f190ae46efabe65a647f1f`.
Status: `actual-axial-taper-crossing-and-necessary-cutoff-pullback-budget`.

The previous checkpoint ruled out a generic cancellation argument with a
manufactured local family. Here we use the RETAINED pre-modulation axial
taper itself. It crosses the entire variance-cutoff window at nondegenerate
angular points. An exact pullback recovers the unchanged cutoff from the
actual variance target. Consequently the parent's condition (S) implies a
specific NECESSARY shape inequality, in addition to the earlier SUFFICIENT
one. An angularly flat crossing cannot bypass this gap on the actual taper.

These are analytical consequences of the retained native profile, exact
earlier annular match and parameter bounds, not an independent existence
review. The actual cutoff is still not quantitatively specified: neither
(S) nor its failure is established. No primitive H bound is disproved;
(S) is a sufficient input to that bound, not claimed necessary for it.
No field, cutoff, root, support, width C^-20, frequency or pressure datum is
changed. Full pressure/axis/Cp/axial memory remains intact. FINAL global
jets, radial/state acceptance, physical pulses, complete covariance/means,
full viscous/cutoff/interaction residual and blow-up remain unproved.

Source: retained [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed pp.129--130 (A.5)--(A.9), and p.160 (C.9). PDF SHA256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The final theorem is not a premise. The native matching and parameter
sources are `annular_exterior_foundation`, `actual_fifth_native`,
`higher_jet_budget` and `axis_core_attachment`.

## 1. A genuine native slice with nonzero angular crossings

Use the original local log-radius y on the axial reduction, after the first
unit ramp. The SAME profile is

`E=P1*f(eta)*exp(-y/2)`, `U=k(y)*eta`, `f=(1+eta^2)^-1`,
`k(y)=4*(1-sigma(log(1+y)/Md))`, `Md=64`, `Td=exp(64)+10`.

It precedes the first reserved compact-repair interval, so it is part of
the compact PREMOD input. The earlier annular edits have ended and match
the ideal fields exactly there; they do not alter this formula. The first
ramp gives `exp(-1/2)*P<=P1<=exp(1/10)*P<2P`. Here the retained choice has
`log P=exp(64)+11`. At

`y0=exp(32)-1`, `k'(y0)=-(1/2)*exp(-32)`.

This is strictly inside the axial taper, before its final eleven flat units.
Writing Eamp=P1*exp(-y0/2), the actual shear definitions give

`a=1-2E_y/E=2`, `bs=2U_y/E=2k'(y0)*eta*(1+eta^2)/Eamp`,
`vs=a+bs^2/a=2+beta*eta^2*(1+eta^2)^2`,
`beta=2*k'(y0)^2/Eamp^2=(1/2)*exp(y0-64)/P1^2>0`.

No pressure coordinate is set to zero or recomputed from truncated memory.
This identity involves native fields only and does NOT replace actual p2
by the zero-pressure diagnostic of the preceding checkpoint.

Let A=C^4096 and delta=deltaL=exp(-A^20), unchanged. The inherited
logarithmic parameter certificate gives P1<A. Therefore

`exp(-65)*A^-2 < beta < 1`.

The lower bound uses y0>0 and log(2)<1. The upper bound uses
P1>=exp(-1/2)P and the actual log P. Set r=sqrt(delta/beta). Then

`log r < 65/2+log A-A^20/2`, hence `r<1e-6`.

This last inequality holds for every A>=2^16: use log A<=A,
log(1e6)<14, and A^20/2>2A>A+47. Only logarithms and finite universal
comparisons are evaluated; no actual A, beta, delta, G, H, N or physical
radius is materialized. Retained-parameter comparisons are checked with
outward intervals at 80 and 110 digits.

## 2. Pull the ENTIRE unchanged cutoff back from the actual target

Let u(eta)=eta+eta^3. On [0,1] this has derivative >=1 and inverse f_u.
For t in [0,1], define

`q(t)=r*sqrt((1+t)/8)`, `eta(t)=f_u(q(t))`.

Then 0<eta(t)<1 and

`vs(eta(t))=2+delta*(1+t)/8`.

Thus this one actual radial slice visits every scalar point of the cutoff
transition; reflection gives the negative-eta branch too. Its angular
gradient is nonzero, with

`sqrt(beta*delta)/(2sqrt(2)) <= vs_eta <= 8sqrt(beta*delta)`

on the positive branch. Indeed eta lies between q/2 and q, and
`vs_eta=2beta*eta*(1+eta^2)*(1+3eta^2)`.

Define the SAME actual shape `K(t)=zetaL(2+delta/8+(delta/8)t)` and
the actual variance target `s=zetaL(vs)*sqrt((2+delta/2-vs)/a)`.
On this slice a=2, so the exact identity is

`s(y0,eta(t))=(sqrt(delta)/4)*K(t)*sqrt(3-t)`.             (1)

It is the selected cutoff, not a replacement by the radial step sigma.
No assumption about its derivative size or monotonicity is used.

## 3. Quantitative inverse-composition bound through order four

Use normalized C4 norms `sum_(j=0..4) sup |D^j f|/j!`. For the inverse
of eta+eta^3, raw derivatives at eta=f_u(q) are

`f_u'=1/(1+3eta^2)`,
`f_u''=-6eta/(1+3eta^2)^3`,
`f_u'''=(90eta^2-6)/(1+3eta^2)^5`,
`f_u''''=(360eta-2160eta^3)/(1+3eta^2)^7`.

Their normalized absolute bounds on 0<=eta<=1 are respectively
1,3,16,105. Also |q|C4<r, by the four binomial coefficients of sqrt(1+t).
Finite normalized composition gives

`|eta(t)|C4 <= r/2+r+3r^2+16r^3+105r^4 <=126r<1`.

For any actual s on the full angular interval, the Taylor jet of composition
therefore obeys `|s(y0,eta(t))|C4 <= |s(y0,.)|C4`. This uses the complete
inverse jet, not an evaluation at a frozen angular point. Since
`|(3-t)^(-1/2)|C4<2`, (1) and the factorial Banach-algebra product yield

`|K|C4 <=8*delta^-1/2*|s(y0,.)|C4`.

Let `S=sum_(j=0..4) sup_(y,eta) |partial_eta^j s|/j!` on the same
compact PREMOD rectangle, exactly the uniform norm used by the parent, and
`Z=max(1,max_(j=1..4) sup |K^(j)|)` as in the prior checkpoint. Since
K(0)=1 and j!<=24, the justified actual-profile lower estimate is

`S >= (sqrt(delta)/192)*Z`.                            (2)

This is a necessary-budget extraction, not an upper bound proving (S).
It makes no assertion that the source-dependent loop or exact repair fails.

## 4. The precise remaining inequality

For the unchanged G=exp(A^32), condition (S), namely S<=G, NECESSARILY gives

`log Z <= A^32+A^20/2+log(192)`.                       (N)

The earlier sufficient implication remains

`log Z+12log(10)+20673log A+8A^20 <= A^32`.             (Z)

Neither inequality has been evaluated for the actual retained cutoff. A
violation of (N) would refute the specific S input on this actual taper;
satisfying (N) would NOT prove S. Conversely (Z) would prove S with the
retained generic native estimates. The gap between necessity and sufficiency
is explicit. The prior local obstruction remains historical evidence;
today's result additionally excludes angular flatness of all actual crossings.

## 5. Checks and next obligation

`pullback.py` checks the exact shear identity, inverse derivatives, positive
crossing, normalized composition constants and the necessary/sufficient
direction. Manufactured floating diagnostics recover a cutoff and its
fourth derivative through two independently computed inverse maps. They are
not interval proofs of actual fields or evaluations of the enormous hierarchy.
All historical source bindings and a fixed manifest of frozen records are
verified without modifying them. Scope-promotion and necessary-as-sufficient
mutations are rejected.

```sh
python -m unittest -v test_audit.py
python audit.py --verify-record evidence.json --output /tmp/actual-cutoff-pullback-fresh.json
```

Use a fresh output path. Next: identify/enclose the unchanged actual Z or
otherwise bound its complete target composition. Then join radial/state
cone and FINAL global derivative estimates before unconditional physical
pulse, full covariance/mean and residual conclusions are possible.
