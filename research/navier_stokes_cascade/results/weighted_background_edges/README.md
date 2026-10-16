# All-order edge weights and a smooth relative background stress

Status: `all-order-weighted-edges-and-relative-background-stress-audited`.

The next estimate is stronger than a generic weighted perturbation:

\[
 \boxed{D^I\left[\frac{\widehat T-T_0}{|T_0|}\right]
                  =O_I(q^{2h})}\tag{1}
\]

on the CLOSED active annulus, with the quotient extended at both flat edges.
It follows from the ACTUAL first-stress edge factorization and common
interior supports, with additional finite seminorms in the already free
Borel cutoff choice. It does NOT follow just from a bound involving
zeta delta^-N. We also join the leading all-order edge bounds, derive
the weighted stress comparison, and preserve the actual background's
growth/oriented cone on one common small-q domain.

This is an analytical consequence of the explicitly cited, replayed
continuum lemmas, not an independent reproof of every inherited lemma or
a proof-assistant certificate. The numerical controls are not the enormous
actual profile. Uniform physical PULSE estimates, complete covariance
means, the complete residual improvement and blow-up remain unproved.

Parent: `ede56e1149c4ea3728f1376039755a951f6fd71f`. Source defining equations:
(4.26)--(4.27), (5.13), (5.22), (5.42)--(5.46), (6.23) in the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
PDF pages 44 and 60 were visually checked for the exact inner weight and
derivative convention. Its final theorem is not a premise. Every historical
source/evidence file and GOAL.md remains unchanged.

## 1. Fixed parameters, weight and all-order leading regularity

Keep Md64, the same h, t1=C^-20 and finite N=1+floor(Hhier^32).
Write ya=log(X/Xa), yb=log(Xb/X), and use EXACTLY

\[
 \zeta=\exp[-t_1^2/y_a^2-4/y_b^2],\qquad
 \delta=\min(1,y_a,y_b),\qquad X_a<X<X_b.\tag{2}
\]

Extend zeta by zero. Neither replacing t1^2 by one nor treating different
exponential constants as comparable is valid. If ca is replaced by ca',
the inner ratio is exp[(ca'-ca)/ya^2] times a smooth nonzero factor.
If ca'>ca, no polynomial upper loss can bound that ratio; if ca'<ca,
the required positive lower stress/zeta bound fails.

The [core contraction](../core_contraction_audit/README.md) supplies common
analytic-eta regularity and every fixed radial derivative of the exact
core. The [activation quotient](../inner_collar_activation/README.md)
propagates them through the flat activation, including its endpoint.
Subsequent reference stages, finite-N loops and primitives are smooth.
This assertion does not promote their low-order quantitative envelope to
an all-order bound uniform in N. For the variance root the continued
F_mu is strictly positive, also at zero variance; the phase density is
bounded away from zero. The finite-dimensional repair root has uniformly
invertible Jacobian B+DQ(c): the accepted contraction bounds
||B^-1 DQ(c)||<1/2. Differentiate its pointwise equation in eta;
at EVERY fixed order m the only highest derivative is multiplied by that
same Jacobian. All remaining terms use lower root derivatives and smooth
input jets. Induction and the smooth implicit-function theorem give a
smooth root. The C1 small root and this root coincide by uniqueness.
No derivative-order-dependent reselection of N is needed; constants may
contain arbitrary finite powers of the already fixed N. The same argument
applies to the exact joins. The Gamma-integral heat formula has integrable
majorants for each fixed derivative, including eta=+/-1.

Thus finalized leading fields and moment primitives are C-infinity on the
needed closed rectangles, using these defining-map lemmas. The
[regional join](../regional_cone_join/README.md) and local
[growth audit](../core_growth_activation/README.md) give positive shear,
growth and oriented cone gaps, including the limiting stress directions.
No numerical resolution of the actual finite N is being asserted.

## 2. Join every fixed edge derivative, retaining both exponential constants

On the inner collar, s=ya/t1, the exact factorization is

\[
 T_0=e^{-1/s^2} A_a(s,\eta),\qquad
 A_a=(1-\kappa_0)g(s)B(s,\eta),\quad |A_a|\ge c_a>0.\tag{3}
\]

Here A_a and all fixed jets are smooth on a smaller closed collar.
The positive lower norm follows from the core-growth activation estimate,
not an approximate stress-zero tolerance. The outer factor is

\[
 T_0=e^{-4/y_b^2}y_b^{-3}(b_0,y_b^6b_z),\qquad
 T_1=e^{-4/y_b^2}y_b^{-3}(b_1,0),\quad b_0>0,\tag{4}
\]

with smooth b0,bz,b1 and b1>=0, from the
[direct heat-edge audit](../outer_edge_stability/README.md), now identified
with the actual forward stresses by exact moments and common supports.

For E_c(y)=exp(-c/y^2), write E_c^(j)=E_c P_j(1/y). Then

\[
 P_0=1,\qquad P_{j+1}(z)=2cz^3P_j(z)-z^2P'_j(z),
                    \quad\deg P_j\le3j.\tag{5}
\]

For example P2=4c^2 z^6-6c z^4; repeatedly multiplying the logarithmic
derivative alone would miss the second term. Equations (3)--(5) bound
inner derivatives by exp(-t1^2/ya^2) ya^-3j and outer derivatives by
exp(-4/yb^2) yb^(-3-3j). Physical X derivatives replace logarithmic ones
by X^-j times a finite signed-Stirling combination. Since Xa>0 is fixed,
these costs are finite, not uniform as Xa or t1 vary.

The omitted exponential factor in (2) is smooth, positive and has smooth
reciprocal on the opposite edge collar. On the compact middle both
zeta and |T0| are positive; all fixed jets are bounded. Use overlapping
collars ya<=tc/16, yb<=1/2 and middle tc/32<=ya,
yb>=1/4 (shrinking collars if necessary). They cover the annulus.
For every total profile derivative order m, this proves

\[
 |T_0|\ge c\zeta,\qquad
 |\partial^I T_0|+|\partial^I T_1|
       \le C_I\zeta\delta^{-3m-3},\quad |I|\le m.\tag{6}
\]

T1 is identically zero on a fixed inner collar: the coefficient equations
are solved exactly there, before their cutoff. Every n>=2 stress has
common compact support K=[Xinner,Xv] inside (Xa,Xb), by the
[coefficient support audit](../coefficient_support_audit/README.md).
Thus Tn/zeta, extended by zero at the edges, has bounded jets of every
fixed order; these constants depend on n. Exact source identities and
the zero total moments are essential to these support assertions.

## 3. A stronger quotient, specific to this construction

Set W_n=T_n/|T_0| in the active interior. At the inner edge W1=0 on an
open collar. At the outer edge (4) gives the cancellation-safe formula

\[
 W_1=\left(\frac{b_1}{\sqrt{b_0^2+y_b^{12}b_z^2}},0\right),
 \quad W_1(X_b,\eta)=\left(\frac{2\eta^2X_b}{L^2},0\right).\tag{7}
\]

It is smooth jointly in X,eta, including eta=0. In the compact middle
|T0| never vanishes. For n>=2, Wn is smooth and supported in the SAME K.
Therefore EVERY Wn extends smoothly over the closed annulus. The large
factor Xb in (7) is retained; it is not a numerically small constant.

A generic weighted bound cannot substitute for (7): with leading angular
stress exp(-4/delta^2)delta^-3, the perturbation
-q^(2h)exp(-4/delta^2)delta^-4 reverses its sign when delta<q^(2h).
No fixed positive q prevents this. Exact supports and the matching
flat-factor powers, not small absolute errors, avoid that obstruction here.

## 4. One diagonal schedule for all finite jets and all smaller q

Use the same Stokes potential gauge and 0<=chi<=1, with chi=1 through
1/2 and zero from 1 onward. At step n>=2 add the finite seminorms

    all profile jets through order n of Wn and Tn/zeta

to the velocity/potential and physical requirements already in
[background jet foundation](../background_jet_foundation/README.md).
The stress itself is summed directly, without a potential curl:

\[
 \widehat T=T_0+\sum_{n\ge1}\chi(c_nq)q^{2nh}T_n.
\]

They are finite by sections 2--3. For order m<=n let M_n,m bound both
quotients and their profile jets, and let K_m bound the fixed Euler cutoff
jets through m. A sufficient constant is

\[
 C_{n,m}=\max\{1,2K_m(1+2nh)^m M_{n,m}\}.
\]

With sigma=c_n q, Euler jets of chi are fixed sigma jets, independent of
c_n. Choose, recursively, one finite log(c_n) at least

\[
 1+\max\left\{\log c_{n-1}+\log2,\log q_0^{-1}+\log2,
       \max_{m\le n}\frac{\log C_{n,m}+n\log2}{nh},
       \text{the existing physical/velocity requirements}\right\}.\tag{8}
\]

Here 0<q0<=1 is fixed. This is an admissible specification of free Borel
scales, NOT a numerical evaluation of the actual coefficient constants.
Euler/frozen conversions retain their Bell/factorial costs; the radial
velocity retains the extra potential cutoff derivative and M_n term.
Neither h,N nor any field support is changed. On the ENTIRE support
q<=c_n^-1, (8) implies C_n,m q^(nh)<=2^-n, hence

\[
 |D^I[\chi(c_nq)q^{2nh} W_n]|,
 |D^I[\chi(c_nq)q^{2nh} T_n/\zeta]|
           \le2^{-n}q^{nh},\quad |I|\le n.\tag{9}
\]

D consists of q partial_q at fixed X,eta, partial_X and partial_eta.
For fixed m split at J=max(m,2). The finite block n<=J is O_m(q^(2h));
the tail is bounded by 2^-J q^(h(J+1))<=2^-J q^(2h) for q<=1.
The sums are locally finite for q>0 because c_n grows at least geometrically.
Thus (1) holds on ONE domain 0<q<=q0 for every fixed derivative order,
with order-dependent constants. From (6), (9) and the product rule for
zeta, we also obtain the actual weighted comparison

\[
 |D^I(\widehat T-T_0)|\le C_I q^{2h}\zeta\delta^{-3m-3},
 \qquad |I|\le m.\tag{10}
\]

For the tail, derivatives of zeta cost at most delta^-3m; T1 determines
the harmless extra three. No uniform-in-n coefficient bound or convergent
uncut formal series is required.

## 5. A common background cone domain, but not a pulse/residual theorem

Let n0=T0/|T0| with its smooth edge extension, and R=(That-T0)/|T0|.
By (1), choose a SINGLE smaller qstar>0 with sup|R|<=1/2. Then

\[
 \widehat T=|T_0|(n_0+R),\quad
 \tfrac12|T_0|\le|\widehat T|\le\tfrac32|T_0|,\quad
 n_B=\frac{n_0+R}{|n_0+R|},\quad
 D^I(n_B-n_0)=O_I(q^{2h}).\tag{11}
\]

The normalized tangential velocity and first profile derivatives differ
from leading by O(q^(2h)), from the now supplied regularity/support and
cutoff premises of the background foundation. Radial velocity remains
O(q^h) in q^A normalization; it is not discarded.
On the leading closed annulus the smooth gap functions
F, |g|^2, -Ntheta, lambda^2, -n0 dot N and
1-c0^2(n0 dot K)^2/(n0 dot N)^2 have strictly positive minima.
Here N=g/|g|, K=(-Nz,Ntheta), c0=lambda/(2F Ntheta).
On a fixed compact neighborhood avoiding their denominator zeros they
have finite Lipschitz constants. Write rho for its radius, gamma for
the minimum gap and Lgap for a padded Lipschitz constant. If Cbounds
bounds the combined velocity/shear/direction perturbation by Cbounds*q^(2h),
it suffices to shrink qstar so

    Cbounds*qstar^(2h) < min(rho/2, gamma/(2 Lgap), 1/2).

Such a positive qstar exists at the SAME h,N. This is a value/first-jet
selection, not an order-by-order shrinking domain. Higher derivatives
change constants only. The summed background retains strict growth,
positive orientation and a strict cone margin, with |That|>=c*zeta/2.
No pressure equality on I4 is inferred: Pi1=K1-mu^2/(4X) and Pi2=K2
remain in all subsequent means, along with full axial velocity terms.

The next obligation is to apply this actual background certificate to the
physical pulse, carrier rounding, all-scale coordinate charts and source
supports, checking all derivative conventions. Complete-wave covariance
means and the full residual including viscosity/cutoffs/interactions are
still separate. The augmented background residual is not a proof of the
complete Navier--Stokes residual improvement or blow-up.

## 6. Reproduction and scope

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output NEW.json` here.
The parent support record is regenerated and compared, with 564 inherited
hashes checked. New controls cover exact flat-derivative polynomials,
both weight constants, all-smaller-q diagonal inequalities, finite-block
and tail exponents, factored heat-edge quotients including eta=0,+/-1,
and the generic weighted-perturbation counterexample. Manufactured
constants and heat-tail parameters are labelled; no actual giant profile
or actual Borel scales are numerically evaluated. The continuum argument
is sections 1--5, not the finite tests.
