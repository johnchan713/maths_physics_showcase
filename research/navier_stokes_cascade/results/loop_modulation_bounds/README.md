# Fixed-phase loop and modulation error bounds

Status: `fixed-phase-loop-and-modulation-errors-bounded`.

This checkpoint supplies the third estimate requested by the
[joined-stress review](../joined_stress_construction/README.md): first slow
derivatives of the mean-preserving loop, including its inverse phase map,
and constants for the resulting state and five-moment errors. It uses the
actual-profile input bounds of the
[compact envelope](../compact_jet_envelope/README.md). Those earlier analytic
arguments remain dependencies; reproducing their arithmetic does not independently
prove them. The new argument below is a mathematical estimate accompanied by
equation checks, not a formal certificate or independent peer review.

The proposed bounds **Cstate = D = H^16 are sufficient** for the modulation
alone, with the definitions below. The repair-to-state constant Ccorr and its
positive-field neighborhood remain unverified. Consequently the proposed
frequency is still not accepted, the full stress construction is unfinished,
and no Navier-Stokes solution, smooth forcing, or blow-up is claimed.

Parent: `bb72b484feb0437800ea0a8cf9a2c7e1aefab527`.
Source formulas: Appendix C, particularly (C.9)-(C.16), of the supplied
[manuscript](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
Its SHA-256 remains
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`;
printed page 160 was visually inspected. The quantitative derivative and
error budgets here extend the preceding repository arguments.

## 1. Fixed inputs and derivative order

Write y=log X. Keep every parameter, support endpoint, and axis pressure datum
from the preceding checkpoint. Put

\[
A=C^{4096}\ge 2^{16},\quad G=\exp(A^{32}),\quad H=\exp(G^4),
\quad d_0=(16A)^{-1},\quad \delta=\exp(-A^{20}).
\]

The letter H in this hierarchy is unrelated to the physical quantity
\(H_x=\sqrt{2X}E\). In this note a slow derivative means either partial_y or
partial_eta **before substituting the oscillating phase**.

The input envelope bounds E,U, all five moments and Pi through mixed order
two; the shear and integrated stress coordinates through order one; and the
listed positive reciprocals and repair normalization multipliers. Each such
norm is at most A. In particular a>=1/A, |ts|,|p1|,|p2|,vs<=A, where
ts=-bs/a and vs=a(1+ts^2). The radial region has A^-1<X<A and length less
than A. Also lambda^-1<A. None of these is inferred from the test fixtures.

We prove only **C0 errors of the four state coordinates** (a,b,p1,p2), and
**C1_eta errors of the moments**. The latter use first slow derivatives of
the loop, which use first derivatives of the input coordinates and hence the
already established second derivatives of the original moments. No extra
derivative has been silently assumed. A C1 error for p1,p2 is not asserted.

## 2. A nonsingular variance-root derivative

Use the exact variance root and smooth cutoff from the preceding stress note.
Let q(z,theta)=exp(z sin(theta))/I0(z), z=mu p, p=p2, and F(mu,p)=sqrt(V(mu,p))
for mu>=0, continued with its sign across zero. The positive variance series
establishes this smooth continuation at mu=0. The identities

\[
q_z=q(\sin\theta-g'),\quad
q_{zz}=q[(\sin\theta-g')^2-g''],\quad g=\log I_0
\]

give |q_z|<=2 exp(2|z|), |q_zz|<=5 exp(2|z|), since |g'|<=1 and
0<g''<=1. The integral identity (q(z)-1)/z=int_0^1 q_z(rz)dr avoids division
by a small p. It gives

\[
|t_\mu|\le2d_0e^{2|z|},\qquad
|t_p|\le3d_0\mu^2e^{2|z|},\qquad
|t_\theta|\le d_0\mu e^{2|z|}.
\]

The tilted density is at least exp(-2|z|). Minimizing its variance over the
subtracted constant therefore gives g''(z)>=exp(-2|z|)/2. Applying this on
the interval between z and 2z, and using I0(2z)/I0(z)^2>=1, yields
V_mu>=d0^2 mu exp(-4|z|). The inherited excursion bound sqrt(V)<=16 d0 mu
then gives

\[
F_\mu\ge(d_0/32)e^{-4|z|},\qquad
|F_p|\le3d_0\mu^2e^{2|z|}.
\]

The second estimate follows by differentiating the L2 norm of t-ts and
using the preceding bound on t_p. At mu=0 the exact values are
F_mu=d0/sqrt(2), F_p=0; division by V_mu there would be invalid.

Let zeta(vs) be the flat step which changes from one to zero between
2+delta/8 and 2+delta/4. Its derivative is at most 512/delta. On its support
set Delta=2+delta/2-vs>=delta/4 and

\[
s=\zeta(v_s)\sqrt{\Delta/a},\qquad F(\mu,p)=s.
\]

Extend s by zero off that support. Flatness makes the extension smooth.
Since Delta<3, for either slow coordinate x,

\[
|s_x|\le1024A^3/\delta,\qquad
\boxed{\mu_x=(s_x-F_p p_x)/F_\mu.}
\]

This covers both zero pressure and the zero-variance boundary. Merely
asserting smoothness of the root would not give the bound needed here.

## 3. The inverse phase map cannot be frozen

The previous uniform cap is

\[
Z=16(1+2048A^5)^2,\quad \mu_{max}=64A^{3/2}e^Z.
\]

For A>=2^16, log(64)+(3/2)log(A)<A; the difference increases with A.
Also Z<=16(2049)^2 A^10. Thus
log(mu_max)<[16(2049)^2+1]A^10<A^16, since
16(2049)^2+1<(2^16)^6. Adding log(A) for |mu p|, or using
|t|<=A+mu_max/A, leaves a logarithm below 2A^16<A^32.
Hence mu, |mu p|, |t|,
A, 16A, delta^-1 and all fixed constants used below are smaller than G.
The excursion estimate supplies |t|<G without exponentiating a crude bound
on q. These comparisons are symbolic; the huge quantities are not rounded
to infinity or materialized.

For G>=2^16, 10^6 G^10<exp(G). This one inequality absorbs the explicit
polynomial coefficients in the following estimates. Its logarithmic
difference is increasing on that interval. Also |v_x|<=4096 A/delta<G^3,
2<v<=A. The preceding root estimates give
|mu_x|<exp(7G), |t_x| at fixed theta <exp(10G), and
|t_theta|<exp(3G).

Define the phase lift with a fixed origin,

\[
w=\frac{a(1+t^2)}{2\pi v},\quad
\phi(x,\theta)=\int_0^\theta w(x,r)dr,\quad
\theta=\Theta(x,\phi).
\]

The exact variance equation gives int_0^(2pi) w=1. Thus Theta is a smooth
inverse circle map. Quantitatively G^-3<w<G^3 and

\[
w_x=w\left(\frac{a_x}{a}+\frac{2tt_x}{1+t^2}-\frac{v_x}{v}\right),
\quad \Theta_\phi=1/w,\quad
\boxed{\Theta_x=-\phi_x/w.}
\]

Here phi_x is taken at fixed theta; Theta_x is taken at fixed phi. On a
fundamental period, |phi_x|<=2pi sup|w_x|. Periodicity extends the same bounds
to every phase. It follows that |w_x|<exp(11G), |phi_x|<exp(12G), and
|Theta_x|<exp(13G). In particular,

\[
\boxed{\partial_x[t(x,\Theta(x,\phi))]
=t_x+t_\theta\Theta_x,\qquad |\partial_x t|_\phi<e^{17G}.}
\]

Differentiating v/(1+t^2) and -vt/(1+t^2) now gives first slow derivatives
of aL,bL below exp(18G): both scalar rational factors and their first
derivatives have absolute value at most one. Their slow J1 norms are below
exp(19G). Multiplication by E in the second periodic source costs at most A.
The zero-mean periodic primitive operator has sup norm at most two and
commutes with slow differentiation at fixed phase. Consequently the two
primitives script-A, script-B satisfy

\[
\|\mathcal A\|_{J_1},\ \|\mathcal B\|_{J_1},
\ |\mathcal A_\phi|,\ |\mathcal B_\phi|
<e^{22G}<e^{32G}<H.
\]

All bounds are uniform in phase and on the selected compact rectangle.
Both primitives vanish in the protected collars. Numerical checks explicitly
retain a failure when the inverse-phase derivative is omitted.

## 4. Transfer to fields, moments, and pressure coordinates

Use the inherited modulation E_N=E exp(script-A/N), U_N=U+script-B/N,
evaluated at phi=Ny. Assume only N>=4AH in the estimates of this section.
The angular C1 norm is |f|_sup+|f_eta|_sup; it is a product algebra.
The phase is independent of eta, so

\[
\|E_N-E\|_{C^1_\eta},\|U_N-U\|_{C^1_\eta}
\le d:=4AH/N\le1.
\]

Indeed ||exp(script-A/N)-1||C1<=2H/N. The modified fields have C1 norm
at most 2A, and E_N/E lies between 1/2 and 2. In particular the new inverse
E and inverse H_x have value bounds 2A. No comparison of radial derivatives
of the two fields is used.

Integrate the five ORIGINAL moment densities on the modulation support.
Product differences satisfy ||delta(UE)||C1<=3Ad and
||delta(U^2-E^2/2)||C1<=5Ad. Since the support length is below A,
sqrt(2X)<=2A, and 1/X<=A on that support, we obtain:

| Moment | C1_eta error at any radius through the repair patch |
|---|---:|
| M | A d |
| I | 2 A^2 d |
| J | 6 A^3 d |
| S | 5 A^2 d |
| Cp | (3/2) A^3 d |

Thus use **m=8 A^3 d=32 A^4 H/N** for each moment and Pi. The axis pressure
datum is unchanged. The Cp integration starts at the positive left support
endpoint, not at zero with a fictitious uniform lower bound on X.

For completeness, propagate these errors through the original (4.16).
Use H_x=sqrt(2X)E and let B be the numerator of the term divided by X H_x
in Qs. All four coefficients in B have magnitude at most one; hence
|B|<=4A and |delta B|<=4m. The input bounds include |W|,|Ns|<=A. Direct
subtraction, retaining the new denominators, gives

\[
|\Delta W|\le2Am,\quad |\Delta(H_x^{-1})|\le4A^3d,
\quad |\Delta(E^{-1})|\le2A^2d,
\]
\[
|\Delta Q_s|\le10A^2m+16A^5d,\quad
|\Delta N_s|\le12A^2m+Ad,
\]
\[
|\Delta p_1|\le10A^4m+16A^7d\le96A^7d,
\quad |\Delta p_2|\le24A^5m+4A^5d\le196A^8d.
\]

These are value estimates, not estimates on angular derivatives of p.
In Ns the pressure terms contribute at most 4m; dropping them is not allowed.

The exact shear identities from (C.13) give
|a_N-aL|<=2H/N and |b_N-bL|<=8AH/N. The exponential factor in b_N is
retained. Therefore a common four-coordinate state error constant is

\[
\boxed{C_{state,raw}=1024A^9H.}
\]

## 5. The transformed repair target keeps its lambda loss

Each of the five individual normalization multipliers has angular C1 norm
at most A, including the derivatives of e(eta)=E(rhop,eta). Only AFTER
normalization do we combine the M and J rows and divide by lambda. Thus the
transformed discrepancy has C1 norm at most

\[
\frac{2A}{\lambda}m
=\frac{64A^5H}{\lambda N}
\le\frac{64A^6H}{N}.
\]

The signs of the repair target do not alter its norm. We bound the difference
before subtraction: no nearly equal rounded M and J values are used.
A well-conditioned transformed matrix does not remove this input loss.

The inequalities 1024 A^9 H<H^16 and 64 A^6 H<H^16 prove that the inherited
proposal **Cstate=D=H^16** is sufficient for these two roles. They also hold
at the proposed N=1+floor(H^32), which exceeds the auxiliary condition
N>=4AH. This check does not certify that frequency: the complete acceptance
condition still needs the actual repair-to-state estimate and its neighborhood.

## 6. Review scope and next obligation

The exact chain-rule calculation and estimates close the fixed-phase loop
and modulation-error part of the current research plan, conditional on the
preceding actual-profile construction and envelope. Manufactured checks
exercise zero and nonzero pressure, root differentiation, an inverse phase
map, and differentiation of modulated fields. They are not evaluations of
the actual joined profile or numerical evidence of a singularity.

The next named obligation is to prove Ccorr and r for the actual first
reserved correction patch. This must propagate the five correcting bump
coefficients through field positivity, shears, cumulative moments and (4.16).
Only after that may the existing finite-frequency criterion be applied.
Even a completed leading stress construction would leave the full PDE
corrections, their convergence, admissible smooth forcing, and actual
finite-time divergence to be proved.

## Reproduce

With mpmath 1.3.0, from the repository root:

```sh
python research/navier_stokes_cascade/results/loop_modulation_bounds/test_audit.py -v
python research/navier_stokes_cascade/results/loop_modulation_bounds/audit.py \
  --verify-record research/navier_stokes_cascade/results/loop_modulation_bounds/evidence.json \
  --output loop_modulation_reproduced.json
```

`bounds.py` checks the scalar domination at two outward precisions and records
the derived error ledger. `review.py` independently differentiates the
manufactured implicit constructions, with refinement and deliberate omissions.
The evidence hashes the argument, code, and inherited inputs. Passing these
checks is narrower than proving every assertion in the manuscript.
