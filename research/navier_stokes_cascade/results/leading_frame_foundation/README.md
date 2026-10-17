# Current leading profile: uniform spectral floors and the actual moving frame

Parent: `889f2e5ba0837a3ea90abdacd7ae37d92f366b3d`.
Status: `current-leading-frame-floors-from-reviewed-core-and-same-root-regional-profile`.

This checkpoint removes the separate positive leading/frame FLOOR inputs
from the current primary-root and sourced-inverse estimates. It derives
them from the SAME reviewed exact-core, modulation, repair and compensated
exterior maps, not from the final manuscript theorem or a sampled huge
profile. It supplies explicit finite constants for all bands and labels.
It does not independently formalize the infinite-core contraction, close
all generated nonlinear covariance/mean inputs or prove a full residual.
No root, carrier rule, frequency, width `C^-20`, support or pressure gauge
is changed. Historical evidence is preserved.

The manuscript defining equations are (4.11), (4.20), (7.1)--(7.11), SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The live input is the [primary-covariance parent](../primary_covariance_foundation/README.md),
the [same-root pressure join](../global_pressure_moment_foundation/README.md)
and their reviewed regional/core estimates, rather than an independent
assumption of Theorem 4.6.

## 1. One closed-annulus certificate from the current defining maps

Let `Aenv=C^4096`, `G=exp(Aenv^32)` be the ALREADY selected hierarchy.
Use new names for the following bounds, not for new construction parameters:

\[
 \alpha={1\over2G^3},\qquad \nu={1\over2G^4},\qquad
 A_b=B_b=4A_{env},\qquad V_b=A_b+B_b^2/\alpha.       \tag{1}
\]

Here `a=-2 partial_y log F`, `b=bs`, `y=log X`, and `v=a+b^2/a`.
On the SAME closed annulus `[Xa,Xb] x [-1,1]`,

\[
 \alpha\le a\le A_b,\qquad |b|\le B_b,\qquad
 2+\nu\le v\le V_b.                               \tag{2}
\]

The regional source of EACH claim matters:

* Protected inner part: the actual activation has `a=k ar>C^-22`,
  `-b=k qr Er/E`, `0<k<=1`, and `Er/E=exp(-J)<=1` since `ar>0`.
  The core envelope gives `a<=10`, `|qr|<=CR` with `R=exp(B/100)<C`,
  and the reviewed protected growth is `v>2.29`.
* Compact loop plus exact repair: the four loop coordinates have magnitude
  at most `2Aenv`; the COMPLETE four-coordinate error is below
  `epsilon/2<1`. The current bounds give `a>=G^-3/2` and `v-2>=G^-4/2`,
  where each displayed `/2` means division by two, NOT a fractional exponent.
  Thus (1) is an unambiguous rewriting of the retained bounds.
* The protected right power has `a-2=2lambda>C^-1`. The SAME compensated
  exterior has `U=0`, hence `b=0`, `a-2>=h`, `v<120` through its heat edit.
  The terminal heat interval has `2+h<a<=2+2h`. Its endpoint is included
  by the smooth formula. Retained `lambda>C^-1`, `h=lambda^2` and `G>C`
  imply `nu<C^-2<h`, so no terminal growth floor is lost.

For C>=2 the hierarchy gives `log G=C^131072>32 log C`, hence
`G>C^32`. This also makes `alpha<C^-22`, as needed on the inner part.

These intervals are exactly those joined in the current global pressure
checkpoint. The modulation and both repairs remain the same roots; fields
alone would not identify their stresses or gauge.

Positivity of F does not come from dividing by a sampled velocity. The core
has `F=g Phi>0`; activation and modulation use positive exponentials.
The annular/compact additive repairs stay in their proved half-relative
positive neighborhoods; the heat multiplier and compensation are positive.
All fields join smoothly. The core anchor at Xa satisfies

\[
 f_a={e^{-2\Lambda M}\over4C}\le F(X_a,\eta)
 \le f^a={2e^{2\Lambda M}\over C}.                  \tag{3}
\]

The exact logarithmic ODE and (2), with `Yb=log(Xb/Xa)>0`, now give

\[
 F(X,\eta)=F(X_a,\eta)e^{-{1\over2}\int_{\log X_a}^{\log X}a\,dy},
 \qquad f_*=f_a e^{-A_bY_b/2}\le F\le f^a.          \tag{4}
\]

All constants are finite for the frozen profile. The floor can be extremely
small: (4) is a logarithmic analytical bound, not a number asserted to be
representable in binary64. No uniformity as C,h or the annulus change is
claimed. A mere positive function on an OPEN annulus would not imply (4).

## 2. Exact cancellation and uniform leading frame denominators

For `g=F(-a,b)`, `N=g/|g|`, `K=(-Nz,Ntheta)`, the identities are

\[
 |g|^2=F^2av,\quad -N_\theta=\sqrt{a/v},\quad
 \lambda_0^2=2aF^2(1-2/v),\quad c_0^2=(v-2)/2.       \tag{5}
\]

The last formula is obtained by dividing the defining growth equation by
`4F^2 Ntheta^2`; neither a pressure omission nor an approximation is used.
Consequently (2)--(4) give

\[
 g_*=f_*\sqrt{\alpha(2+\nu)}\le|g|\le f^a\sqrt{A_bV_b},
 \quad n_*=\sqrt{\alpha/V_b}\le-N_\theta,             \tag{6}
\]
\[
 \lambda_*^2={2\alpha f_*^2\nu\over V_b}
 \le\lambda_0^2\le2A_b(f^a)^2=(\lambda^*)^2,
 \quad {1\over2G^2}=\sqrt{\nu/2}\le|c_0|
 \le\sqrt{(V_b-2)/2}.                               \tag{7}
\]

These are the missing STRICT leading denominator inputs, supplied from
the reviewed same-profile value estimates. The current joined unit stress
direction retains its strict orientation/cone on the closed annulus even
where stress itself is zero. Its actual fixed tilt u and cone margin are
unchanged; no new tilt is selected from a possibly weaker coarse bound.

In a physical band `s=q/Q in [1/4,4]`, the leading F,g,lambda acquire the
SAME positive factor `s^(-1-h)`; N,K,c0 do not. Thus multiply their lower
bounds by `4^(-1-h)` and upper bounds by `4^(1+h)` before applying the
phase formulas. No q or derivative-order loss occurs in this conversion.

## 3. Carrier rounding and an explicit frame inverse

In this section use the chart bounds just described. Keep the EXISTING
tilt u, `1<=epsilon*k^2<=4`, and the actual definition of Bs. Then

\[
 B_*^2={\lambda_*\over4(1+u^2)^{3/2}}
 \le B_s^2\le{\lambda^*\over(1+u^2)^{3/2}}=(B^*)^2. \tag{8}
\]

The unchanged axial wave number satisfies
`pz=Bs(Kz-sigma*u*gz/(Ls*|g|^2))`. Since `Kz=Ntheta<=-nstar`,

\[
 L_s\ge {2u\over g_*n_*}\quad\Longrightarrow\quad
 p_z\le-B_s n_*/2<0.                                \tag{9}
\]

Only the angular carrier is rounded. Nearest NONZERO rounding costs at
most `1/k`, not `1/(2k)`. In the actual terminal shear `N=(-1,0)`,
`K=(0,-1)`, the unrounded `p/R0=sigma*Bs*u/(Ls*|g|)`.
The profile constants remain fixed while Ls grows, so `|p/R0|` tends to
zero including its `1/(k R0)` rounding error, while `pz=-Bs` stays nonzero.
The independent algebraic carrier control sets `R0=Bs=u=|g|=1`
(not a selected physical pulse). For `k>=Ls^2`, it has `|p|<=2/Ls`.
Thus angular periodicity does NOT provide a band-uniform lower bound for
`|p/R|`. It is the full tangential normal that must be bounded below.

Use the reviewed comparison `|n-Bs(sigma*(u/2+uv/Ls),K)|<=D/Ls`.
One finite value condition `Ls>=2D/Bstar` gives `|n_tan|>=Bstar/2` and
`|n|>=Bstar/2`; it also bounds `|sa|` by a fixed constant. Equation (9)
is a useful additional explicit axial-carrier floor, not a claim that the
moving axial normal cannot be affected by its retained epsilon*v term.

For the ACTUAL unit vectors Ka,Na and real sa, let
`U=[er-sa Ka,Na]`, `M=[[1,1],[c0 sqrt(1+sa^2),-c0 sqrt(1+sa^2)]]`.
The original left inverse is `Bleft=M^-1[er^T;Na^T]`. Exact algebra gives

\[
 B^\ell B=I_2,\quad n^TB=0,\quad
 \|B^\ell\|_2\le {1\over\sqrt2\min(1,|c_0|)}
 \le\nu^{-1/2}=\sqrt2G^2,                           \tag{10}
\]
\[
 \|B\|_2=\sqrt{2(1+s_a^2)}\max(1,|c_0|).             \tag{11}
\]

Indeed the two rows of M are orthogonal, with singular values sqrt(2)
and sqrt(2)*|c0|*sqrt(1+sa^2). The two rows er^T,Na^T are orthonormal.
Equation (10) concerns the SPECIFIED ambient left inverse, not an unnoticed
replacement by the Moore--Penrose inverse. The latter can be smaller.

## 4. One value domain, every fixed jet, and the existing inverse/root join

The finalized leading profiles and unit stress direction have all fixed
jets on the closed compact annulus by the reviewed core, activation,
finite-N modulation, same-root repair and outer Gamma constructions.
Only the denominator/value floors (6)--(8) are needed to choose a fixed
small enlargement; continuity preserves half their positive margins there.
For the cone, use the same smooth DIRECTION, not T0/zeta, which can be
unbounded at the outer edge. At eta=+/-1 and T=0 the relevant jets are
one-sided; no negative-time heat extension is assumed.

On that one enlargement, quotient differentiation gives every fixed jet
of N,K,c0 and the leading frame. Background differences O(epsilon^2),
full-interval rounding error `Ls/k`, and the phase comparison then give
the reviewed `E,d-dref=O(1/Ls)` with finite profile-dependent constants.
Once `|n_tan|>=Bstar/2` holds, derivatives of the ACTUAL B,Bleft have
polynomial S costs without a new smallness condition for each order.
No derivative of the frozen representative, integer rounding or sign is
taken. The fixed extension and one finite value-domain selection are
different from a derivative-dependent starting band.

The homogeneous Riccati inverse now uses the explicit
`mu=lambda_star/sqrt(1+(3u/2)^2)>0`. The old minimum length is at least
`25`, `8CE/mu`, and the retained covariance gap budget. Add the displayed
finite normal/carrier comparisons to its preliminary value acceptance,
if not already implied by its constants; do not claim a numerically
evaluated first giant integer or change the prescribed frequency/width.
Since `Ls>=2r0 ell^2`, all these conditions hold on one later half-line.
The same exponential majorant controls `S^2(epsilon+epsilon^2+1/k)`.

The sourced inverse and original primary roots therefore have their
leading/frame-floor inputs supplied by this current profile certificate.
They retain physical Haar costs, original positive branches, entrance
compatibility and full normal-source pressure. In particular the term
`(n.Kbg*t-n'.t+n.f)/(k*m*|n|^2)` is not erased by transversality.
This conclusion does not supply every generated source/complete covariance
target, fresh axial/radial mean or Pi1/Pi2 pressure update in the nonlinear
cycle. Those and every viscous/cutoff/interaction residual still need their
own estimates. The exact-core analytical dependencies have not been turned
into a proof-assistant certificate or an unconditional blow-up result.

## 5. Reproduction and controls

The checks use independent rational shear identities, actual frame/left
inverse products, both carrier signs, nonzero rounding and a disappearing
angular-component counterexample. A logarithmic manufactured profile shows
why underflow is not a zero floor; open-domain/zero-growth controls show
why mere pointwise positivity is insufficient. None evaluates the enormous
selected profile or proves a continuum quantifier by sampling.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
The output must be fresh. Analytical estimates are (1)--(11), conditional
only as stated on the reviewed exact-core/current global defining maps.
