# Q154: physical diagonal closed-time and punctured-endpoint join

Parent: Q153, `a22f25e5f14763dc24c695233144e9be46a359fa`.
The live branch, Q153 evidence and CI status, Q7's endpoint obstruction, the
certified complete-potential endpoint theorem, endpoint foundation, current
endpoint jets and physical-pulse foundation were read first. Historical
evidence, the finalized physical family, roots, pressure gauge, supports,
accepted bands, diagonal cutoffs and `C^-20` width are unchanged.

Q7 correctly isolated a missing uniform closed-time bound at its then-current
head. Q149--Q153 now supply the missing same-family inputs: every actually
used band lies on one derivative-independent admitted tail; the physical
primary has all fixed one-sided jets; generated sources, compact maps and
tails have finite all-fixed-order ledgers; and Q153 binds them to one locally
finite complete-potential diagonal. This checkpoint records that later
correction and proves the punctured endpoint implication for that same
physical field. The singular origin remains excluded.

## 1. Keep the same representatives and pressure gauge

Put

\[
 \tau=1-t,\quad A=\tfrac12+h,\quad D=\tfrac12-h,
 \quad q=\tau+z^2q^{2h},\quad
 X={r^2\over2q},\quad \eta={z\over q^D}.               \tag{1}
\]

At the selected `h=1/200`,

\[
 L=1-2h\eta^2\ge1-2h={99\over100}                     \tag{2}
\]

on the closed terminal chart `T+eta^2=1`. No division by `1-eta^2` is
needed at `T=0`.

Use Q153's exact fixed-gauge representatives

\[
 \mathcal A=A_{\rm seed}+\sum_{j\ge1}\chi(a_jq)A_j,
 \qquad
 \mathcal B=B_{\rm seed}+\sum_{j\ge1}\chi(a_jq)b_je_\theta,
\]
\[
 u=\operatorname{curl}\mathcal A+\mathcal B,qquad
 p=p^{[0]}+\sum_{j\ge1}\chi(a_jq)(p^{[j]}-p^{[j-1]}). \tag{3}
\]

`A_seed` contains the background, wave and mean potentials. `B_seed`
contains the divergence-safe axisymmetric direct swirls. All cutoff-curl
terms stay in `curl A`; no compact potential is invented for direct swirl.
The full seed pressure, both exceptional same-base tails, changing
covariance pressure, axial velocity and moving pressure moment remain in
one gauge. No endpoint inversion, gauge reset or new cutoff sequence is
introduced.

## 2. The Q7 blocking inequality is now supplied conditionally

The exact Q7 premise was, for every punctured terminal compact `E` and every
fixed `alpha,s`,

\[
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^s(\mathcal A,\mathcal B,p)|<\infty. \tag{4}
\]

Interior smoothness and a residual flat at the origin do not imply (4).
The later physical chain supplies it as follows.

First use the complete accepted band range

\[
                     {1\over4}\le {q\over Q}\le4.       \tag{5}
\]

On `c<=q<=qmax`, active dyadic bands form the finite closed window

\[
 \max\{\ell_0,\lceil\log_2(1/(4q_{\max}))\rceil\}
 \le\ell\le\lfloor\log_2(4/c)\rfloor.                 \tag{6}
\]

For `c=2^-24`, `qmax=2^-20`, this is `ell=18,...,26`; the narrower
`1/2<=q/Q<=2` window would omit both boundary bands.

Let

\[
 q_e=\min(q_\star/16,2^{-\ell_0-2}),\qquad
 \tau\le q_e/16,\qquad |z|\le(q_e/16)^D.               \tag{7}
\]

The implicit equation gives `q<=3q_e/8`. Since `Q<=4q`,

\[
 4Q\le16q\le6q_e<q_\star,qquad
 Q\le\tfrac32q_e<2^{-\ell_0}.                          \tag{8}
\]

Thus every enlarged chart in the observation cylinder lies inside the same
admitted physical domain and beyond the retained first-band threshold. This
uses Q149's actual later-tail construction rather than silently certifying
an unspecified historical first band.

For `q>=c`, each slow box has finitely many mesh centers, including the
already selected one-sided `T=0` indices `-1,0,1`. If diagonal stage `j` is
active, `a_jc<1`; the retained `a_j>2^j/qstar` gives

\[
                         j<\log_2(q_\star/c).           \tag{9}
\]

Only finitely many bands, boxes and correction stages survive on the
compact. Q149 supplies fixed physical potential, pressure and complete-curl
jets uniformly on every such actual band. Q153 supplies every finite
generated-source quotient, smooth support extension, compact request and
tail bound on that same family. Finite products, constrained ODE variation,
complete curls, pressure/moment maps and compact inverses preserve one-sided
fixed-order bounds. Arbitrarily growing constants in stage or derivative
order are harmless because both are fixed before the finite set is chosen.
This proves (4), relative to those named continuum contracts.

The pulse-path endpoint and `t=1` are distinct interfaces: Q153's smooth
pulse cutoff supplies source support, while (4) uses the physical one-sided
`T=0` estimates and finite closed atlas. Neither is substituted for the
other.

## 3. Cover every punctured terminal compact

Let `|x|>=d>0` on a terminal compact in (7), and choose

\[
 b=\min(q_e/4,d^2/(16X_{\rm ext}),1/4),\qquad
 \zeta=\min(d/4,(b/16)^D),\qquad \tau\le b/16.          \tag{10}
\]

On `|z|>=zeta/2`, equation (1) gives a positive lower bound for `q`, so
(4) follows from the finite atlas. On `|z|<zeta`, one has `r>d/2`,
`q<=3b/8`, and

\[
                X\ge {d^2\over3b}>{16\over3}X_{\rm ext}. \tag{11}
\]

That whole equatorial neighborhood is the unchanged heat exterior. The two
regions overlap and cover positive-radius equatorial points and axis points
with `z!=0`; they deliberately do not cover `x=0`.

## 4. Complete exterior heat and pressure jets

In the heat exterior, the fixed-gauge representatives are `A=0`, direct
swirl `K e_theta`, and

\[
 K=C r^{-1-2h}H(4\tau/r^2),\qquad
 H(Z)={1\over\Gamma(1+h)}\int_0^\infty
 e^{-v}v^h(1+Zv)^{-h}\,dv,
\]
\[
 p=-\int_r^\infty K(\rho,\tau)^2{d\rho\over\rho}.      \tag{12}
\]

Set

\[
 A_j=(h)_j(1+h)_j,qquad
 S_j=\sum_{k=0}^j{j\choose k}A_kA_{j-k}.               \tag{13}
\]

Dominated differentiation yields

\[
 \partial_t^jK(r,1^-)=C4^jA_jr^{-1-2h-2j},
\]
\[
 \partial_t^jp(r,1^-)=-{C^24^jS_j\over2+4h+2j}
 r^{-2-4h-2j}.                                         \tag{14}
\]

The full convolution `S_j` is essential. Already

\[
                         S_2-2A_2=2A_1^2>0.            \tag{15}
\]

Radial derivatives retain both changing homogeneity and derivatives of
`Z=4tau/r^2`; Cartesian frame costs are bounded on positive-radius
compacts. These estimates complete (4) on the second region without
discarding pressure cross terms.

There is a sharp limitation. For `C!=0`, the ratio of successive swirl
time-Taylor coefficients is

\[
 {4\over r^2}{(h+j)(1+h+j)\over j+1}>{4j\over r^2}.    \tag{16}
\]

The pressure coefficients have the same factorial obstruction through
`S_j`. Both terminal time Taylor series have zero radius, even though every
one-sided smooth jet exists. A future extension must use smooth jet gluing;
analytic continuation or zero extension of the nonzero pressure jets is
invalid.

## 5. Compatible limits and remaining boundary

Using (4) at one extra time derivative gives

\[
 \|D_x^\alpha\partial_t^sG(t_2)-
 D_x^\alpha\partial_t^sG(t_1)\|_E
 \le M_{\alpha,s+1}|t_2-t_1|,qquad
 G=(\mathcal A,\mathcal B,p).                           \tag{17}
\]

Every fixed jet therefore has a uniform one-sided limit. Spatial coordinate
integration identifies its spatial derivatives and time integration
identifies successive time jets. For every fixed finite `k`,

\[
 \left\|D_x^\alpha\left[G(1-\tau)-
 \sum_{j=0}^k{(-\tau)^jT_j^G\over j!}\right]\right\|_E
 \le {M_{\alpha,k+1}\tau^{k+1}\over(k+1)!}.            \tag{18}
\]

Curl commutes with these finite limits. For a later localization request
with `m=a+j`, the safe finite ledger remains

\[
 \mathcal A:m+4,\qquad \mathcal B:m+3,\qquad p:m+2,
 \qquad \mathcal R:m+1.                                \tag{19}
\]

This closes the punctured endpoint implication for the Q153 physical field,
conditionally on the retained continuum providers. It does not regularize
the singular origin, construct a future solution, or prove an admissible
global force.

The next obligation is to replay the exact axisymmetric space-time
localization and zero initialization on this same field. That replay must
retain unit viscosity, every potential-cutoff derivative, the full
fixed-gauge pressure cutoff, both nonlinear cross orientations and the
cutoff-curl square before any smooth force extension, energy, comparison or
core lower-bound claim.

Passing tests checks endpoint arithmetic, complete pressure convolution and
frozen provenance. It does not independently prove the continuum estimates.
No unconditional pulse family, compact force, global residual theorem or
blow-up is claimed.

Reproduce with fresh output paths:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q154_physical_endpoint_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q154_physical_endpoint_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q154_physical_endpoint_join/evidence.json \
  --output q154_replay.json
```
