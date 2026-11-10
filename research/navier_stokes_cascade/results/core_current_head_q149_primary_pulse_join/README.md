# Q149: rule-defined schedule correction and uniform primary pulses

Parent: Q148, `a69f150b0c1e87f2080538293be971434731c3be`.
The live branch, Q148 evidence and CI, the rule-defined background Borel
join, same-band gap join, canonical full-path ODE, fixed projector family,
compact/slow inverse domain and physical pulse transfer were read first.
Historical evidence, roots, profiles, supports, pressure datum, tilt,
frequency and the `C^-20` width remain unchanged.

Q148 correctly bound the finite positive-order coefficient/support theorem
and sharp reserved-pressure floor to the finalized Q147 profile. Its stated
next obstruction was too strong: the live repository does not use a
separately frozen numerical Borel scale list. The exact-real prescription
chooses each scale only after the corresponding finite continuum seminorms
exist. This checkpoint records that correction and joins the already audited
physical primary-pulse estimates to the same finalized profile.

## 1. Correction: the diagonal row follows from the retained rule

For every fixed coefficient order `n` and derivative order `m<=n`, Q148
provides a finite complete coefficient seminorm. The Borel provider includes
the implemented cutoff derivatives, physical rows, radial potential-curl
commutator and full pressure memory in a finite constant `C_(n,m)`.

Write `g_n=nh`, `b_n=log(a_n)`, `v_(n,m)=log C_(n,m)` and `ell=log 2`.
The retained exact-real rule chooses `b_n` strictly above

\[
 1+\max\left\{b_{n-1}+\ell,\ B+\ell,
 {v_{n,m}+n\ell\over g_n}:0\le m\le n,
 \hbox{all finite physical rows}\right\}.              \tag{1}
\]

Every entry is finite before this choice. Hence the whole sequence exists,
`a_n>2a_(n-1)`, and on `q<=a_n^-1`,

\[
 C_{n,m}q^{nh}\le2^{-n}e^{-nh}<2^{-n}.                 \tag{2}
\]

Equivalently, in the binary form used by Q120, a certified bound
`log2 C_(n,m)<=c_(n,m)` and a chosen `ell_n<=log2 a_n` satisfy

\[
 nh\,\ell_n-c_{n,m}-n\ge nh>0.                         \tag{3}
\]

Thus Q148's displayed diagonal inequality is a consequence of the actual
rule-defined construction, not an external unresolved numerical-scale
premise. Q120's warning remains valid for an unrelated pre-existing numeric
scale list; no such list is promoted here. This correction neither evaluates
the enormous scales nor changes them.

The locally finite Borel sum is therefore the same completed background:
it keeps every `Pi_n`, including the two exceptional mean-patch pressures,
and includes every cutoff-curl term. For any fixed jet its tail is bounded by

\[
 \sum_{n\ge N}2^{-n}(q^h)^n
 ={(q^h/2)^N\over1-q^h/2}\le2^{1-N}.                  \tag{4}
\]

No coefficient bound uniform in `n` and no convergence of the uncut formal
series is asserted or needed.

## 2. Same-profile primary-pulse join

The rule-defined background and same-band provider supply, on one
derivative-independent later tail:

1. positive background value, spectral, same-tilt cone and phase-normal
   gaps with the full `F,G` transport; and
2. the exact actual background, including pressure and potential-cutoff curl,
   required by the independently audited physical primary-pulse theorem.

That physical theorem supplies its own whole-path moving-frame/Riccati
existence, moving-normal and source-normal pressure, fixed potential/pressure
cutoffs, derivative transfer and smooth support extension. The separate
canonical full-path checkpoint explicitly records that its canonical family
has **not** been identified with this inherited physical primary family.
This checkpoint preserves that negative flag and does not merge the two.
The fixed-projector and compact/slow-inverse contracts are downstream
infrastructure for covariance and generated corrections, not premises used
to prove the primary estimates below.

Q147 supplies the accepted leading repair, global moments and unchanged
pressure gauge; Q148 supplies the positive-order support ledger on that same
profile. Consequently the independently audited physical primary-pulse
transfer applies without
selecting a new profile, root, band, cutoff or pressure representative.

Let `Q=2^-ell`, `S=ell^2`,

\[
 s_x={1\over2}+{h\over2},\qquad
 s_t=1+{3h\over2}.                                     \tag{5}
\]

For every actual used band `ell>=ell0` and every fixed spatial/time orders
`a,b`, the complete primary fields obey

\[
 |\nabla^a\partial_t^b u_{\rm primary}|
 \le C_{a,b}Q^{-(a+1)s_x-bs_t}S^B\sqrt\zeta\,\delta^{-B}P, \tag{6}
\]
\[
 |\nabla^a\partial_t^b p_{\rm primary}|
 \le C_{a,b}Q^{-(a+2)s_x-bs_t}S^B\sqrt\zeta\,\delta^{-B}P, \tag{7}
\]
\[
 |\nabla^a\partial_t^b {\cal A}_{\rm primary}|
 \le C_{a,b}Q^{-as_x-bs_t}S^B\sqrt\zeta\,\delta^{-B}P.    \tag{8}
\]

The constants and polynomial exponent may depend on the fixed derivative,
profiles, tilt and width, but not on the running band or pulse label. These
are uniform bounds, not absolute smallness as `Q` tends to zero.

The vector potential and pressure have flat cutoffs with closed supports

\[
 |\xi|\le r_0/2,\qquad |v/L_s-1/2|\le2/5,              \tag{9}
\]

inside the fixed radial shell and slow support. Curl does not enlarge that
support. Distinct enlarged labels have zero products, at most one actual
wave label is nonzero at a physical point, and the locally finite family has
smooth zero extensions. The one-sided `tau=0` background jets are retained;
there is no extension through `q=0`. Early bands are deliberately absent
from the working domain, so (6)--(9) cover every scale actually used by the
construction.

This closes the uniform **primary** pulse size, fixed physical derivative
and support obligation. A sourced inverse still requires its stated input
class and smooth endpoint support. No complete nonlinear generated-pulse
hierarchy is claimed.

## 3. First obstruction in the covariance step

The reserved background pressure remains

\[
 \Pi_1=K_1(\eta)-{\mu^2\over4X},\qquad \Pi_2=K_2(\eta), \tag{10}
\]

with all physical derivatives and cutoff terms. The complete wave also has
axial velocity and curl remainders. The known primary physical covariance
remainder carries the exact power

\[
                  Q^{-1-2h\kappa},\qquad \kappa=10^{-5}. \tag{11}
\]

Because `-1-2h*kappa<0`, (11) grows as `Q` decreases. It cannot be promoted
to a decaying error by enlarging a finite constant or by the pulse bounds
(6)--(8). The next checkpoint must derive the complete tensor covariance
and mean cancellation, including axial velocity, radial velocity, all curl
remainders, (10), moving pressure and every cutoff derivative. This is an
analytical obstruction, not a failed numerical tolerance.

Uniform primary pulses are now joined on the finalized profile. Complete
covariance/mean corrections, estimates for every nonlinear generated source,
the full viscosity/cutoff/all-interaction residual and blow-up remain
unproved.

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output fresh.json` here.
The record binds the Q148 evidence and each primitive provider without
editing any historical record.
