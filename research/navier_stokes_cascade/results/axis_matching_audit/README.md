# A quantitative five-moment annulus, conditional on axis entry

This checkpoint constructs the five corrections needed to join an inner
profile to the completed reference exterior. It gives an explicit acceptance
ball for incoming moment functions and axial velocity, proves a small exact
root exists, and bounds the **strict relaxed** stress inequalities throughout
the joining region. The five directions remain uniformly independent after
normalization; a large matching radius does not damage the inverse bound.

**An actual regular axis reaching this acceptance ball is still unverified.**
The numerical incoming data below are manufactured fixtures, not an axis
solution. The activation collar, the conversion from relaxed to fully
admissible stress in Appendix C, the full PDE corrections and admissible
smooth forcing are also unverified. No blow-up result is claimed.

The source's order of choices matters here: the matching tolerance comes
before j0, sigma and Lambda in (B.40). The earlier axis pilot's parameters
cannot simply be carried forward without checking that tolerance.

Source: the supplied [manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Appendix B.1, B.7-B.9, especially (B.1)-(B.3) and (B.35)-(B.40), printed
pages 144-157. The moment matching page 156 was checked visually.
Frozen PDF SHA-256:
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
Parent: `6b995ec8160bb363588d371da41a579f0344cce6`.
The explicit constants, relative swirl bumps and acceptance ball below are
our specialization of that construction, not a verification of the full paper.

## 1. State the incoming conditions before solving the equations

Write P=Pstar, x=X/XR, y=log x, xc=exp(-6), f=(1+eta^2)^(-1).
The argument works for P>=16 and 0<h<=.01, hence includes the selected
Md64 exterior family. On the entry side y=-8, require the fields already
to have the form E=P f x^(1/10), U=4 eta+g(eta), with that form on a
neighborhood to its left. This is the endpoint form of (B.34).

Use the common pressure datum Pi0 fixed by the exterior. Let hats remove
the physical radius scales: M=XR Mhat, I=XR^(3/2) Ihat,
J=XR^(3/2) Jhat, S=XR Shat; Cp has no XR factor. The ideal primitives are

\[
\widehat M_0=4\eta x,\quad
\widehat I_0={5\sqrt2\over8}Pf x^{8/5},\quad
\widehat J_0=4\eta\widehat I_0,
\]
\[
\widehat S_0=16\eta^2x-{5\over12}P^2f^2x^{6/5},\qquad
C_{p,0}={5\over2}P^2f^2x^{1/5}.
\]

For actual-minus-ideal incoming moment differences define, in this order,

\[
r=\left(
{\Delta\widehat M\over x_c},
{\Delta\widehat J-4\eta\Delta\widehat I\over\sqrt2 Pf x_c^{8/5}},
{\Delta\widehat I\over\sqrt2 Pf x_c^{8/5}},
-{\Delta\widehat S-8\eta\Delta\widehat M\over P^2f^2x_c^{6/5}},
{\Delta C_p\over P^2f^2x_c^{1/5}}\right).
\]

Assume these are smooth functions on the closed eta interval and impose

\[
\max_i\|r_i\|_{C^1}\le\varepsilon_0=10^{-16},\qquad
\|g\|_{C^1}\le\varepsilon_0,\qquad
\|v\|_{C^1}=\sup|v|+\sup|v_\eta|.
\]

These are uniform hypotheses. Passing at several sampled angles does not
establish them. For a global axis attachment, the ledger must actually be
realized by a common smooth inner profile with the prescribed Pi0. This
checkpoint does not establish that realization.

`physical_to_rows` accepts all five moment differences with their angular
derivatives; it does not silently assign zero derivatives to bare values.
The derivatives of the combinations 4 eta I and 8 eta M are included.

## 2. Restore the axial field and solve five exact equations

On -8<y<-7 set U=4 eta+g[1-sigma(y+8)], leaving E ideal. All derivatives
match at the ends because sigma is flat. Put t=y+6 and let beta_c(t) be
`sigma'((t-c)/.08+.5)`. On 0<t<1 use

\[
U=4\eta+u_0\beta_{.3}+u_1\beta_{.7},\qquad
E=E_0(1+e_0\beta_{.15}+e_1\beta_{.5}+e_2\beta_{.85}).
\]

The five supports are mutually disjoint and lie strictly inside their patch.
The two U corrections are additive; the three E corrections are relative.
This is why their weights differ from the additive-E weights printed in B.8.

With columns c=(u0,u1,e0,e1,e2), the normalized exact map is

\[
Bc+q_\eta(c,c)=z.
\]

Here z is the negative incoming ledger after axial restoration. The two
diagonal linear blocks integrate beta against these exponential rates:

| Block | Row rates in t | Centers |
|---|---|---|
| U: M, J-4 eta I | 1, 8/5 | .3, .7 |
| E: I, -(S-8 eta M), Cp | 8/5, 6/5, 1/5 | .15, .5, .85 |

There are only two nonzero quadratic rows. Set
`zeta=exp(6/5)/(P^2 f^2)`. They are

\[
q_3=\tfrac12\int e^{6t/5}(\sum e_j\beta_j)^2dt
 -\zeta\int e^t(\sum u_j\beta_j)^2dt,\qquad
q_4=\tfrac12\int e^{t/5}(\sum e_j\beta_j)^2dt.
\]

The mixed U/E term in J is identically zero because its supports are disjoint.
The small U-square term in the energy row remains present, even at large P.

The inverse is controlled without evaluating a bump integral numerically.
For each block, B[k,j]=mass[k] node[k]^j. The nodes are exp(spacing*rate),
and positivity with exact bump mass .08 gives
`mass[k]>=.08 exp((first_center-.04)*rate[k])`. The coefficients of the
Lagrange polynomials for these distinct nodes bound the inverse rows.
Outward bounds are below 76 for the U block and 608 for the E block.
Use the common conservative bound **1000** on the maximum-row C1 norm.
No P, XR, eta or h appears in this linear inverse.

Since sigma'^2<=64 sigma', and
`||zeta||C1<=12 exp(6/5)/P^2<1`, the quadratic norm is below **100**.
The axial restoration changes its two linear rows by at most ||g||C1,
and its energy row by at most ||g||C1^2. Therefore

\[
\|z\|_{C^1}\le3\varepsilon_0,\qquad
8(1000)^2(100)(3\varepsilon_0)=2.4\,10^{-7}<1.
\]

The quadratic contraction gives a unique small exact root with

\[
\boxed{\max_j\|c_j\|_{C^1}\le C_b:=6000\varepsilon_0=6\,10^{-13}.}
\]

Uniform invertibility of the derivative of the map and implicit
differentiation give smooth coefficients; every higher fixed derivative
has its own finite bound. No division by 1-eta^2 occurs at the poles.
At the right endpoint the fields are ideal and all five moment functions
are exactly restored. The shared Pi0 then gives the same pressure, Qs and
Ns as the reference exterior; equality propagates along that exterior.

## 3. Bound the whole joining region, including bump interiors

The original pressure has the representation
`-Pi0/P^2=integral f^(2 theta) dmu(theta)` with positive measure on
0<=theta<=1. Its ideal part is (5/2) f^2. During the first unit transition,
E/P<=exp(.1). Afterward the stored scalar amplitude is at most
`2 exp(-.2) exp(-(y-1)/2)`: the interpolation factor is at most one,
and the terminal cutoff ratio is below two. Consequently

\[
\mu([0,1])\le\tfrac12(5+e^{.2}+4e^{-.4})<5,
\qquad |\Pi_0|/P^2<5,\quad |\Pi_{0\eta}|/P^2<10.
\]

The last derivative bound follows from |4 eta/(1+eta^2)|<=2.
All later pressure-preserving edits retain these bounds on the same datum.

For the ideal profile, direct substitution in (B.35) gives

\[
Q_0={9\over8}-{5h\over8}+2h\eta^2
 +{5\over4}{\eta^2(D+4d)\over1+\eta^2}>1.1,
\quad d=1-\eta^2,\ D=\tfrac12-h.
\]

Also E0/P>=exp(-.8)/2>.22 and |N0|/P^2<30. For the last bound, the
velocity terms contribute at most `(20+32(1+h))/P^2`, the datum terms
at most `4(.5+h)*5+10`, and the remaining radial terms at most
`[5+(25/3)h+25/3] exp(-1)`.

At every point of both patches, the cumulative normalized ledger satisfies
`||r||C1<=3 epsilon0+2 Cb+100 Cb^2<4 Cb`. Disjoint supports give
`||U-4 eta||C1, ||E/E0-1||C1<=64 Cb` and radial derivative bounds
`|U_t|, |(E/E0)_t|<=250000 Cb`. A direct endpoint/central-interval estimate
gives |sigma''|<20000, so these derivative bounds do not use a grid maximum.

The original moment formula now bounds the stress changes. For example,
`|Delta W|<=2 exp(2) ||r||C1`. With V the numerator of Q+W, the transformed
rows give `|Delta V|<=39 sqrt(2) P xc^(8/5) ||r||C1`; the ideal
`|V/(x Hhat)|` is below ten. The reciprocal denominator exceeds .99.
Substitution yields the following deliberately loose bounds; their
arithmetic is checked with outward intervals in `bounds.py`.

| Quantity | Uniform bound on -8<=y<=-5 |
|---|---|
| abs(Q-Q0) | 10^5 Cb |
| abs(N-N0)/P^2 | 10^4 Cb |
| abs(a-.8) | 10^6 Cb |
| E/P, Q | E/P>.2, Q>1 |
| abs(N)/P^2 | <32 |
| abs(P bs) | <3 million Cb |
| abs(w/P), w=N/(EQ) | <160 |

For the N estimate one uses
`||Delta S||C1<=16 xc ||r||C1+4 P^2 xc^(6/5) ||r||C1` and
`||Delta Cp||C1<=4 P^2 xc^(1/5) ||r||C1`, in addition to Delta W and
the axial field bound. The large pressure scale is removed explicitly:
the product bs*w equals (P bs)*(w/P).

Thus with `G=Q-bs N/(a E)` and `vs=a+bs^2/a`,

\[
G\ge1-{(3\,10^6 C_b)32\over .7\cdot.2}>.99,
\qquad
v_s\le.9+{(3\,10^6C_b)^2\over16^2\cdot.7}<1.
\]

Finally Pc=XR x G/L, L<=1, x>=exp(-8). The explicit choice

\[
\boxed{X_R\ge10000}\qquad\Longrightarrow\qquad
P_c\ge.99\,10000e^{-8}>3.32>2
\]

proves the strict relaxed condition throughout this annulus. The previously
sufficient exterior floor XR>=100 remains compatible, but is insufficient
for this join: even the ideal profile at eta=0, y=-8 has Pc<.038 when XR=100.
This rejects that radius for the stated joining criterion; it does not
contradict the earlier exterior bounds, which concerned other intervals.

## 4. Compatible axis data and the still-open entry obligation

This tolerance now supplies concrete axis targets. It is enough to require
the core normalized ledger at Xsep to have C1 norm <=epsilon0/4, and

\[
j_0\le\varepsilon_0/12,\quad
\|U_{\rm nat}(4/\Lambda)-U_*\|_{C^1}\le\varepsilon_0/12,\quad
\|G_i-U_{\rm nat}(4/\Lambda)\|_{C^1}\le\varepsilon_0/12.
\]

Then ||Gi-4 eta||C1<=epsilon0/4. The intervening constant axial offset
adds at most that amount and its square to the normalized core ledger
before y=-8, because the relevant positive powers are integrable at zero.
These conditions imply the entry hypotheses of section 1.

Choose **j0=1e-18**, **sigma=j0/100=1e-20**. These are axis data, not
a complete parameter choice for a solution. Their separation condition in
(B.2) can be made quantitative using the same pressure representation.
Write

\[
Z_*=-\eta P^2 K(\eta)+R(\eta),\qquad
K\ge1.25,\quad |R|\le40|\eta|+5j_0.
\]

The lower bound follows from the ideal part of the positive pressure measure;
the remainder bound follows directly from U*=4 eta+j0. Set
`delta_star=j0 P^2/100`. If |Zstar|<=delta_star and P>=16, then

\[
{|\eta|\over j_0}\le
{1/100+5/256\over1.25-40/256}<{1\over30}.
\]

It follows that Hstar>=.8 j0 on that region and
`chi=Hstar^2/(Hstar^2+sigma^2)>.9998>.99`. The unique zero of Hstar lies
between -j0/4 and -j0/5. There its pressure contribution dominates the
velocity remainder and Zstar>.2 j0 P^2>0. Also -Wstar>2.8 on the whole
parameter interval.

A common complex rectangle is available with
`rho=sigma/[4(.5+4+48+4*.05)]`. On it |L|>.9 and
`|Hstar +/- i sigma|>=3 sigma/4`. The same estimate used by the earlier
axis pilot gives `|zeta_star|<=M=4(1+8h)[1+5(8+.05)]/sigma^2`.
Once Lambda is fixed, `log C>=Lambda(2M+1)` ensures |g|<=exp(-Lambda)
on that rectangle. This is an amplitude bound, not a verification of
the nonlinear analytic contraction threshold for Lambda.

Increasing C later sets `XR=110(CP)^10` and
`log xsep=Tsh-10(log C+log P)`. Tsh must first be bounded independently
of C; xsep<exp(-8) and small core moments must then both be established.
A large radius alone neither bounds the core ledger nor removes Gi-4 eta.

**The next unresolved step is quantitative construction of the nonlinear
axis and its activation continuation at these new data.** No Lambda is
certified here. The native axial error, activation error and core ledger
have not been shown to satisfy the displayed budgets. The old j0=.05 lies
outside this sufficient entry budget; that does not prove every possible
join at the old parameters is impossible.

## 5. Numerical checks, retained failures and reproduction

`moments.py` evaluates the exact polynomial moment map with composite
positive Gauss-Legendre quadrature. The audit uses orders 32/48, eight
panels and 280/320 digits. It tests two smooth manufactured entry functions
at five angles and at two pressure scales. The large scale uses the actual
Md4 values P=exp(exp(4)+11), h=exp[-8(exp(4)+10)]; the historical Md4 outer
failure is retained. Neither fixture supplies an axis. The numerical datum
`-3.31462273001433255 P^2 f^2` lies inside the proved pressure envelope;
it is not represented as a recomputation of the complete outer pressure.

The exact existence statement uses the integral bounds and contraction above.
The numerical Newton root is a root of quadrature approximations. A separate
homogeneous expansion retains its quadratic grade and a positive truncation
budget: if a=1000 max(abs(z)) and r=4*1000*100*a<1/2, the remainder after
N grades is at most `a*r^N/(1-r)`. This budget does not certify quadrature
or floating-point rounding error. Refinement checks compare the changes
being made, including the stress perturbations, rather than allowing an
order-one background to hide errors in a tiny edit.

DOP853 independently integrates the original five moment densities with
linear and quadratic pieces in separate units. Tanh-sinh integration checks
the U-square and E-square pieces separately. The original radial source
equations are checked against derivatives of the moment primitives. The
implicit eta derivatives of the root are compared to independent re-solves.
The independent moderate-scale integrations use 100 digits, as in the
construction tests; their accuracy thresholds are unchanged. The 280/320-digit
refinements separately retain the much smaller finite-h terms at the Md4 scale.
Recorded stress extrema concern the declared 440 samples; the whole-interval
claim comes from section 3, conditional on its incoming hypotheses.

The audit retains these failures:

- A linear correction leaves a nonzero quadratic pressure discrepancy.
- Solving only four moment equations leaves the fifth pressure row unmatched.
- At eta=0, an odd fixture has zero correction values and nonzero correction
  derivatives. Freezing those derivatives corrupts Q after the apparent join.
- XR=100 fails the stated joining pressure criterion even for the ideal field.
- Shrinking j0 while leaving sigma=.002 can violate the axis alternative;
  sigma must be chosen consistently with j0.
- The first hand reduction of Q0 incorrectly used `-h eta^2/2` in place of
  `+2h eta^2`. Direct evaluation of the five primitives detected the gap
  `(5/2)h eta^2`. The corrected formula and a normalized regression retain
  this finite-h error; no acceptance threshold was relaxed.
- A geometrically large radius does not certify incoming core moments or
  repair an axial offset outside the chosen acceptance ball.

Run from the repository root:

```sh
python research/navier_stokes_cascade/results/axis_matching_audit/test_audit.py -v
python research/navier_stokes_cascade/results/axis_matching_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/axis_matching_audit/evidence.json \
  --output axis_matching_reproduced.json
```

The protocol and source hashes bind the evidence. All earlier audit sources
and their evidence remain unchanged. Passing this audit establishes the
declared conditional annulus checks, not actual axis entry, full stress
realization, a solution of the full PDE, or finite-time blow-up.
