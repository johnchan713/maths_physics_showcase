# Q114: current reduced-source polynomial degree and fixed log-height

Parent: `cf00be3b42b791c0d3dc4df4bfafcfcc16fc16b1` (Q113).
On 3 November 2026 the inspected remote was Q113; research CI 41563298464
and field CI 41563298469 were running. The frozen core-growth checkpoint
still concerns only the leading core and first collar. Q113 has settled
qualitative finiteness of the original fixed compact-shape factors. This
checkpoint isolates the remaining source's coefficient algebra: it gives
an explicit polynomial DEGREE bound and a fixed log-height bound without
expanding products of exponential coefficient sums. It does not verify
native field inequalities or evaluate source order 100005.

## 1. What the coefficient ring actually means

The existing core_symbolic_fresh_mean_join EP class represents

    P(S)=sum_k sum_e c_(k,e) exp(e) S^k, c_(k,e),e>=0,

with FINITE support. Its exponential exponents are fixed real constants,
not multiples of S. They can depend on the fixed derivative order, common
pole, original profiles and frozen width. They cannot depend on the running
band if P is to be a coefficient polynomial. Here S=ell^2 and
 epsilon=2^(-h ell); exp(S) cannot be absorbed into a polynomial in S.
Calling an arbitrary epsilon-dependent native norm a constant would not
supply an EP input. No such relabeling is allowed below.

Represent an upper bound by (D,H), meaning

    0<=P(S)<=exp(H)(1+S)^D, S>=0.                         (1)

Use a distinguished exact-zero symbol. A nonzero rational c is represented
by (0,c), since c<=exp(c); a term c exp(e) S^k by (k,c+e).
For two nonzero bounds the rules are

    addition:       (max(D1,D2), 1+max(H1,H2)),
    multiplication: (D1+D2, H1+H2).                       (2)

The sum rule follows from 2<=exp(1); the product rule is immediate.
Induction proves (1) for every FINITE positive expression. There is no
rounding or exponential evaluation in this certificate, and zero is
handled exactly. It is deliberately conservative, not an asymptotically
sharp log-height. EP's existing fixed exponential constants are not new
exponential growth in S. This is a sufficient coefficient compression,
not a proof that an arbitrary smooth field has polynomial norm growth.

## 2. Apply the algebra to the current Q104/Q110 reduced source

Let X,U,Z be the degrees of the three COMPLETE six-stage mean norm leaves,
A the degree of their common edge-absorption factor, W the degree of the
TOTAL old-plus-sourced covariance leaf, and b,V,G,c the degrees of the
four fixed-background leaves (degree labels, not their field values).
Put x=X+A, u=U+A, z=Z+A. Q104's exact positive reduced-source polynomial is

    Srr=2 b0 x0+x0^2+W0,
    Srz=b0 z0+G0 x0+x0 z0+W0,
    Stt=2 V0 u0+u0^2+W0,
    P_J=(J+1)(1+4 c0)x0+r_J Srr+(J+1)Srz+r_J Stt
                                      +(J+1)(J+2)x0,
    r_J=sum_(j=0)^J a^(-j-1), a=Rmin>0.                  (3)

Here suffixed symbols are the actual positive norm coefficients; lowercase
unsuffixed symbols above (3) are their degree labels. Under the original
fixed positive normalized hull, r_J is a fixed finite constant. For
nonzero leaves the polynomial degree is bounded by

    D_J=max(c+x, b+x, 2x, W, b+z, G+x, x+z, V+u, 2u).    (4)

The omitted plain x term is dominated by 2x, and the constant part of
1+4c0 is covered since c>=0. For vanishing leaves use (2), which may give
a smaller degree. Apply (2) to (3) to obtain a finite H_J as well. Thus

    P_J(S)<=exp(H_J)(1+S)^D_J.                           (5)

The standalone radial graph multiplier does not occur in (3) or (4).
That is an exact feature of the REDUCED source after Dr K0 is removed
under A0. Graph dependence inside the mean/covariance leaves and the
original chart transfers is retained; it is not set to zero. The full
pressure still contains K0=-Srr+epsilon D1 beta. All six means, both slow
histories, total wave covariance, local viscosity and moving compensator
remain in their original stages. No full-pressure order is promoted.

`reduced` implements (3) in the bound algebra. `SourceGrowth.source` also
compresses an existing Q110.source result, rechecks the identical Q113
state and delegates its original archive/pole/bundle checks. This route
is for supplied finite orders, not an attempt to materialize order 100005.
It never substitutes an abstract leaf coefficient for an actual norm.

## 3. What this discharges, and what is still missing

Conditional on the actual native leaf inequalities that Q104/Q110 use,
(5) gives a single explicit polynomial-growth envelope for the reduced
source's coefficient at its inherited grade B=1-kappa. It separates two
questions: positive polynomial algebra now has a finite structural
certificate; whether the actual current leaves obey those native uniform
inequalities on the whole requested proof domain remains unverified.

For the current output n=0 and kappa=1/100000, Q110 still requires
J=100005, fields through 100007, total covariance through 100016 and
profile jets through 100046, together with a common high pole covering
all nested inverses. To use (5) there one must supply valid band-independent
(D,H) leaves or legitimate EP leaves at precisely those orders. Smoothness
at each epsilon or increasing an archive maximum does not supply them.
The source's own-pole unweighting, whole-line/mixed traces, original raw
R/U chart conversions, fixed kernel profiles and temporal inverse-floor
polynomial must also retain their native proofs. If a hull or pole varies
with the band, its dependence must be bounded as a leaf; it cannot be
hidden in a supposedly fixed H_J.

Consequently this is a conditional source coefficient estimate, not an
actual A0 gain or radial cone closure. Numerical values of enormous fixed
coefficients are unnecessary for this analytical algebra lemma. The
remaining obstruction is the uniform native input inequality on the
consistent current high-jet graph, rather than EP's fixed exp(e) syntax.

## 4. Checks and next obligation

Exact tests compare (4) with the degree of Q104's independent EP expansion
on positive polynomial leaves, include nonconstant absorption, retain all
three mean components and covariance, and verify invariance under changes
to the standalone graph leaf. Moderate log-sum evaluations check (1) on
these expanded polynomials without overflowing exponentials. Typed current
source tests use manufactured leaves and label them as controls; they
verify state/archive routing, not canonical norm enclosures. The analytical
proof is the finite induction (2) and the expansion (3), not these tests.
Frozen replay preserves Q113 and all historical source bindings.

Next prove the actual high-order native leaf bounds and common-pole/trace
contracts, or find a sharper actual A0 cancellation, then join the radial,
angular/axial and inner/outer cones. Native coefficient/reserved-mean
supports, all-scale actual-pulse bounds, complete viscous/cutoff/interaction
residual improvement, infinite realization and blow-up remain unproved.
