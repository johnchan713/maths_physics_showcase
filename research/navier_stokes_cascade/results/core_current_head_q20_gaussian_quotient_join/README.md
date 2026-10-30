# Q20: whole-path Gaussian quotient jet conversion

Parent: `3f38f4e5a48d3d5d5f06b8559f36fadcc09bd755` (Q19).
This supplies the explicit conversion between two distinct norms used by
the source inverse. It does not certify actual source seminorms. All
original labels, endpoints, pressure gauge, roots and width remain fixed.

## 1. The missing conversion is triangular, not division by a small minimum

`core_canonical_entrance_jet_closure`, section 4, correctly bounds
J_j(chi*t)/P rather than J_j(t/P). Q19 retains full lifted-path f/P jets
as an inverse requirement. These quantities are not equal. On the SAME
closed lifted domain suppose P>0 and, for each multiindex I of total order
j<=R, the raw bounds are

    |D^I f|/P <= F_j,    |D^I log P| <= ell_j (j>=1).       (1)

The F_j may contain the retained epsilon order, half-edge weight, finite
poles and scale polynomials. Define

\[
 H_0=1,\qquad H_j=\sum_{i=1}^j {j-1\choose i-1}
                    \ell_i H_{j-i}.                       \tag{2}
\]

Then for either sign

\[
 |D^I P^{\pm1}|/P^{\pm1}\le H_j.                          \tag{3}
\]

Proof: write P=exp(g), fix one derivative of I, and differentiate
D_k exp(±g)=±exp(±g)D_k g another j-1 times. The sum of mixed binomial
coefficients at fixed total split i-1 is binom(j-1,i-1). Induction and
the triangle inequality give (2). It applies to arbitrary mixed jets,
not just repeated derivatives in a single coordinate.

Leibniz therefore gives the explicit conversion

\[
 |D^I(f/P)|\le Q_j,
 \qquad Q_j=\sum_{i=0}^j{j\choose i}H_i F_{j-i}.           \tag{4}
\]

Conversely the same formula bounds |D^I(Pu)|/P from bounds on D^I u.
The positive Gaussian cancels algebraically in every term. Neither ell_0,
1/min P nor any exponential in the pulse length appears. The conversion
uses only source jets through j, not additional source derivatives. If
F_j have a common epsilon/edge factor after finite pole embedding, that
factor is retained; the H_j supply only their own stated coefficient costs.

## 2. Bounds for the retained fixed-label Gaussian

In `physical_pulse_foundation`, P'/P=p=lambda-dref and P(L/2)=1.
The pulse length L=2r0/c_i is fixed during coefficient differentiation
within each label, and L<=C S. Thus

    log P(v,x)=integral_(L/2)^v p(w,x) dw, 0<=v<=L.        (5)

If p_j bounds all order-j mixed jets of the reference rate on this
entire domain, then a parameter-only derivative of log P has bound
(L/2)p_j. A derivative containing a path derivative has bound p_(j-1),
since partial_v log P=p. A sufficient common bound is therefore

\[
 \ell_j=(L/2)p_j+p_{j-1},\qquad 1\le j\le R.             \tag{6}
\]

The retained smooth reference-rate/frame estimates give finite scale
polynomials p_j at each fixed order, with constants depending on the
original frozen parameters. Equations (2),(4),(6) remain finite scale
polynomials. No smallness of high derivatives or new first band is needed
for this conversion. L and the label must not be differentiated as though
they were continuous physical variables. Moving-endpoint or moving-label
coordinates require their endpoint/chain terms before applying (6).
The existing fixed-label chart and its derivative pullbacks are the scope
here; no new choice of reference representative is made.

In particular the tiny value of P at the pulse endpoints does not cause
an exponential reciprocal loss in (4). Smooth zero extension of f/P at
those endpoints still requires the retained source/cutoff compatibility;
P is positive there at each fixed label, but uniform quotient bounds alone
do not assert vanishing traces.

## 3. Conditional propagation through the supported source grammar

Take the all-jet W bounds of Q19 on the WHOLE lifted label domain, with
the same P. The linear regular operator contains coefficient derivatives
of waves and full pressure. Apply (4) at its requested derivative orders,
retaining every derivative of P, to bound its normalized coefficients.
Mean-wave products use a Gaussian-free M factor times a W factor.
For two same-label waves f=P u and g=P v, (fg)/P=P uv; (2) and Leibniz
bound this without dividing a product by a Gaussian minimum. Its radial
weight is the full zeta, as in Q19. Nonzero different-label products vanish
only under the original closed support separation on the extended common
domain. Without that hypothesis, a quotient involving two different
Gaussian envelopes would need its own estimate.

This proves a conditional grammar closure for normalized source jets from
the corresponding whole-domain W/M jet estimates. It is not a derivation
of those input estimates merely from support, point values, epsilon orders
or the finite expression graph. Additive Gaussian/base exterior residual
tails remain outside the forward inverse. Raw source derivative orders
needed by the differential operator must still be requested; (4) introduces
zero extra derivatives beyond that existing request.

For the inverse pressure keep exactly

    pi=i(n.Kbg.t-nprime.t+n.f)/(k m |n|^2).

After normalization by P this uses both t/P and f/P with their complete
jets, alongside the old inverse-normal bounds. A purely normal source can
have nonzero pressure even when the projected velocity response vanishes.
Nothing here deletes that term, the moving-normal term, mean A0/A1
remainders, moving c_rho Pscalar, or the base Pi1/Pi2 pressure memory.

## 4. Limitation: a physical trace is not a whole-path estimate

A bound for a physical evaluation alone does not supply (1) on its lifted
domain. For an elementary counterexample on [0,1]^2, let the observed trace
be y=0 and let f_n(x,y)=n y. Every tangential x jet on that trace vanishes,
while sup |f_n|=n and the transverse derivative is n. Multiplying by a
positive Gaussian does not change the normalized counterexample.
This is a domain-quantifier warning, not an assertion that these are the
actual sourced fields or the actual physical evaluation map.

The precise remaining task is to establish, for the actual recomputed
sources, the whole-domain W bounds at each corrected Q14 derivative and
pole request, with a common band-uniform scale polynomial. Once they are
available, (2)--(6) supply the f/P jet requirement without a separate
Gaussian-division hypothesis. Root/frame/propagator premises, component
specific epsilon orders, full endpoint extensions and the residual
interactions remain required. No actual complete representative bound,
endpoint closure, unconditional force or blow-up theorem is established.

## 5. Evidence

Exact exponential-polynomial fixtures independently check quotient jets,
nonzero Gaussian derivative terms, mixed jets, and the absence of a
zeroth log cost. The finite polynomial recurrence and source-order count
are arithmetic checks, not continuum verification of the chosen huge field.
The Q19 evidence and inherited hash graph are reproduced unchanged.

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q20_gaussian_quotient_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q20_gaussian_quotient_join/audit.py --verify-record research/navier_stokes_cascade/results/core_current_head_q20_gaussian_quotient_join/evidence.json --output FRESH.json
```
