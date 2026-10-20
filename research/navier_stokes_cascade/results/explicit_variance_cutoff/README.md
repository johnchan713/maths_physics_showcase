# Explicit source-admissible variance cutoff and the Z4/Z5 gates

Parent: `e1fbc36677e27239c2da686d75080c3e1d74c068`.
Status: `explicit-source-admissible-variance-cutoff-and-Z4-Z5-gates`.

This checkpoint corrects an inherited premise before using it.  Printed
equation (C.9) does **not** identify an already selected quantitative cutoff:
it says to let `zetaL` be any smooth nonnegative function, at most one, equal
to one below `2+deltaL/8` and zero above `2+deltaL/4`.  Thus earlier wording
such as “the actual retained cutoff” named an admissible class, not a frozen
shape whose derivatives could be recovered.  The historical records are not
rewritten.  We now make one explicit admissible source choice and freeze it.

For this choice, exact analytical bounds give

`Z4 < 3,000,000`, `Z5 < 1,000,000,000`.

They discharge both the earlier fourth target gate and the fifth gate after
combining with the parent's `R5<A^32`.  Hence the conditional same-frequency
compact modulation and exact five-moment repair estimates no longer have a
cutoff-shape premise.  This does not prove the remaining FINAL joint bounds,
the inner/outer physical-scale joins, uniform pulses, complete covariance and
mean corrections, the full residual, or blow-up.

Source: the retained [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
printed p.160, (C.9).  Its SHA-256 is
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The final theorem is not used as a premise.

## 1. The source-completion choice

Put

`phi(t)=exp(-1/t-1/(1-t))` for `0<t<1`, and `phi=0` otherwise.

Let `I=integral_0^1 phi` and define

`K(t)=1-I^(-1) integral_0^t phi`,

extended by one for `t<=0` and zero for `t>=1`.  Finally set

`zetaL(v)=K(8(v-2)/deltaL-1)`.

The bump is positive in the open interval and flat to every order at both
ends.  Therefore K is C-infinity, nonincreasing, takes values in `[0,1]`, and
has exactly the plateaus required in (C.9).  The scale `deltaL`, the transition
window, all radial widths, supports, roots, pressure gauges, moments and the
frequency hierarchy are unchanged.  What changes is that the previously
unspecified member of the admissible class is now fixed.

This does not contradict the historical compressed-cutoff obstruction.  That
obstruction proves that the range and plateau axioms alone give no uniform
derivative bound across the entire class.  Here a particular member is chosen
and bounded.

## 2. Exact first-five derivative enclosure

For `u=1/t`, `v=1/(1-t)`, write

`phi^(n)(t)=P_n(u,v) exp(-u-v)`.

Differentiation uses `u'=-u^2`, `v'=v^2`, so

`P_(n+1)=(-u^2 partial_u+v^2 partial_v)P_n+(u^2-v^2)P_n`.

For every nonnegative integer m, `x^m exp(-x)<=m!`; this follows directly
from `exp(x)>=x^m/m!`.  Summing the absolute monomials of P_n gives the exact
integer majorants

| n | terms in P_n | bound for `sup |phi^(n)|` |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 2 | 4 |
| 2 | 5 | 80 |
| 3 | 10 | 3,600 |
| 4 | 18 | 292,320 |

On `[1/3,2/3]`, `t(1-t)>=2/9`, hence
`I>=exp(-9/2)/3`.  Since `e<3`, `I^(-1)<3^6=729`.  For `1<=j<=5`,
`K^(j)=-phi^(j-1)/I`, and therefore the raw derivative bounds are

`729, 2,916, 58,320, 2,624,400, 213,101,280`.

Consequently the earlier raw quantities satisfy

`Z4<2,624,400 <3e6`, `Z5<213,101,280 <1e9`.

No numerical sampling establishes these suprema; the polynomial recurrence
and positive integral lower bound are continuum estimates.  Finite tests only
check the exact recurrence and bookkeeping.

## 3. Fourth and fifth gates

Let `A=C^4096>=2^16`, `G=exp(A^32)`, with the original
`deltaL=exp(-A^20)`.  Since `log 10<3` and `log A<A`,

`log Z4+12log10+20673log A+8A^20`
` <63+20673A+8A^20 <9A^20 <A^32`.

The middle inequalities hold already at `A=2^16` and strengthen with A.
Thus the sufficient fourth gate (Z) holds, so the actual target selected here
has the required angular C4 norm `S<=G`.  Combining this with the previous
moving-loop and joint-fourth checkpoints discharges their cutoff premise.

For the fifth gate, the parent proved `R5<A^32`.  The remaining left side is

`log Z5+32log A+10A^20+1100000G`.

The non-G terms are below `11A^20`.  Divide by G.  Since
`G=exp(A^32)>A^32>11A^20`, their quotient is below one, while
`G>A^32>1,100,001`.  Hence the displayed left side is strictly below `G^2`.
The parent's B5 implication now supplies the fifth primitive/frequency bound
at the same N and the same exact five-row repair, under its retained
native/regional/root premises.

## 4. Scope and next obligation

This is a correction plus an analytical construction choice, not a numerical
diagnostic and not an independent verification of all inherited foundations.
It closes the explicit Z4/Z5 obligations for the newly frozen source-admissible
cutoff.  It does not turn PREMOD pure-angular and one-radial bounds into FINAL
joint bounds.

Next join the complete compact pressure/flux and repaired-state bounds to the
inner and outer estimates at every retained physical scale.  Then verify the
global pressure/moment identities, coefficient and reserved-mean supports,
uniform actual pulse estimates, full covariance/axial-pressure mean
corrections, and the complete viscous/cutoff/interaction residual.

Run:

```sh
python test_audit.py -v
python audit.py --verify-record evidence.json --output /tmp/explicit-cutoff-fresh.json
```

The output path must be fresh; frozen evidence is never overwritten.
