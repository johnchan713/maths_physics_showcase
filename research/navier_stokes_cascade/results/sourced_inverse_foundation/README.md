# Sourced pulse inverse: explicit existence, harmonic damping and jet bounds

Parent: `93103d22e03b2996ad8cf18ef31f77d090c22896`.
Status: `explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs`.

This checkpoint independently constructs the linear sourced pulse inverse
on the entire already accepted interval and supplies explicit harmonic
dependence, rather than retaining an unspecified finite sourced-inverse
existence lemma. It strengthens the old polynomial-length bound for high
harmonics. It does NOT repeat the completed primary-pulse checkpoint or
claim that all nonlinear targets have now been checked.

The current [coefficient/support foundation](../coefficient_support_foundation/README.md)
and the [physical-pulse foundation](../physical_pulse_foundation/README.md)
are reproduced. Their reviewed leading/frame inputs, exact-core and
background hypotheses remain explicit dependencies. The proof below is
an implication from those inputs, not unconditional construction of the
actual giant profile. All historical evidence is unchanged; parameters,
frequency, original roots, radius, width `C^-20` and supports stay fixed.

The original equations are (7.5), (7.10), (7.12)--(7.20) of the supplied
manuscript, SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The local retained text and PDF were consulted; its final theorem is not
an input. In particular, its finite-harmonic derivative statement is not
silently promoted to an infinite nonlinear construction.

## 1. The precise inputs and the full constrained equation

Write `L=Ls`, `epsilon=Q^h`, `S=ell^2`, `k=ceil(epsilon^-1/2)`.
The fixed accepted domain has `4Q<qstar`, `1/4<=q/Q<=4`, physical `T>=0`,
and `L` comparable to the fixed-radius multiple of S. Discrete data are
frozen before differentiation. The existing moving frame B has a left
inverse Bleft, spans the actual `n`-perpendicular plane, and obeys

\[
 B^{left}(A_nB-B')=\operatorname{diag}(\lambda,-\lambda)+E,
 \quad \lambda>0,\quad d=\epsilon k^2|n|^2\ge d_*>0,
 \quad \|E\|+|d-d_{ref}|\le H/L.                    \tag{1}
\]

Norms here are Euclidean operator norms; a change from another finite
matrix norm must be absorbed in H. B, Bleft, n, lambda, d and their fixed
coefficient jets have polynomial S bounds, with the same nonzero floors.
The positive reference envelope satisfies `P'/P=lambda-dref`, `P<=1`
and the previously reviewed Gaussian bounds. Its fixed normalized
parameter jets also have polynomial S bounds. L is fixed for a frozen
label. The proof only differentiates continuous parameters, not carriers,
floors, representatives or band labels.

For ANY smooth prescribed source f on its full lifted path, not
necessarily transverse and not necessarily small, the equations are

\[
 t'+Kt+m^2dt+ikm n\pi=-f,\quad n\cdot t=0,
 \quad t(0)=0,\quad m\in\mathbb Z\setminus\{0\},
\]
\[
 A_n=-K+\frac{n(n^TK-(n')^T)}{|n|^2},\quad
 t'=A_nt-m^2dt-\operatorname{Proj}_{n^\perp}f,
 \quad\pi=\frac{i(n\cdot Kt-n'\cdot t+n\cdot f)}{km|n|^2}.
                                                               \tag{2}
\]

For each finite m and compact parameter set the matrix and source in (2)
are continuous and bounded. Picard successive approximations converge
uniformly, since the jth difference is bounded by a constant times
`(M_m L)^j/j!`; no contraction smallness or source-size condition is
needed. The same estimate for differences proves uniqueness. Smooth
parameter dependence follows by the finite differentiated equations.
The constraint obeys `(n.t)'=-m^2d(n.t)`, so its zero initial datum is
preserved. Substitution of the displayed pressure cancels the entire
normal component. This proves finite sourced-inverse existence and
uniqueness FROM (1), rather than assuming a separate inverse lemma.
It does not prove the primary covariance roots or the core inputs of (1).

## 2. An explicit high-harmonic residence bound

In coordinates `t=Bz`, put `g=-Bleft Proj(f)` and `y=z/P`. With
`delta=m^2-1>=0`,

\[
 y'=[R-\delta d I]y+G,\quad
 R=\operatorname{diag}(0,-2\lambda)+E-(d-d_{ref})I,
 \quad G=g/P.                                               \tag{3}
\]

The logarithmic norm of R is at most H/L, since its displayed diagonal
part is negative semidefinite, also for complex amplitudes. Gronwall
therefore gives, for `0<=w<=v<=L`,

\[
 \|V_m(v,w)\|\le e^H\frac{P(v)}{P(w)}
      e^{-\delta d_*(v-w)},\quad
 V_m=e^{-\delta\int_w^v d}\,V_1.                            \tag{4}
\]

The scalar damping commutes with the entire noncommuting matrix
evolution. It is not necessary for E to commute with the growing frame.
If `|G|<=F0`, Duhamel yields

\[
 |y(v)|\le C_0F_0\ell_m,\quad C_0=e^H,\qquad
 \ell_m=\begin{cases}L,&|m|=1,\\
 \min\{L,1/[(m^2-1)d_*]\},&|m|\ge2.
 \end{cases}                                               \tag{5}
\]

The frame adds only fixed constants to the bound for t/P. This improves
the old `O(L)` source loss to `O((m^2-1)^-1)` in high modes. The strict
damping floor follows from the reviewed normal floor and
`1<=epsilon k^2<=4`; no new threshold depending on m is required.

## 3. All fixed parameter jets: damping derivatives must be charged

Let I be a continuous slow/transverse multi-index, with v held fixed.
Use ordinary derivatives and the usual multi-index binomial factors.
Let `F_I`, `R_J`, `d_J` bound the derivatives of G, R and d. Differentiating
(3), INCLUDING the `-delta D^J d` term, gives a sufficient recursion

\[
 Y_I\le C_0\ell_m\left[F_I+
 \sum_{0<J\le I}{I\choose J}(R_J+\delta d_J)Y_{I-J}\right]. \tag{6}
\]

All initial differentiated data vanish. For `|m|>=2`, `delta>=3`.
An explicit harmonic-independent bound for `delta Y_I` is

\[
 Z_0=C_0 F_0/d_*,\quad
 Z_I=C_0\left[\frac{F_I}{d_*}+
 \sum_{0<J\le I}{I\choose J}
 \left(\frac{R_J}{3d_*}+\frac{d_J}{d_*}\right)Z_{I-J}\right]. \tag{7}
\]

Indeed `delta ell_m<=1/d_*` and `ell_m<=1/(3d_*)`; induction in (6)
proves `Y_I<=Z_I/delta`. Thus ALL fixed slow/transverse parameter jets
retain the inverse-square high-harmonic gain. Omitting differentiated
damping is not justified, even though it has a favorable sign in values.
Uniformity in m requires the supplied normalized source jet bounds F_I
to be uniform in m; otherwise (6)--(7) retain their actual mode-dependent
input sizes. No decay of arbitrary nonlinear Fourier data is assumed.
One may enlarge polynomial S and edge-distance degrees in this finite
recursion, without changing epsilon order or the accepted domain.
For `|m|=1`, (6) supplies the ordinary finite polynomial-length bounds.

For a generic source this does NOT automatically give inverse-square
pulse-coordinate derivatives: constant normalized forcing gives
`y(0)=0`, `y'(0)=1` for every m. The missing condition is entrance-jet
compatibility, not a failure of damping. If G and all its pulse derivatives
vanish at entrance, the ODE gives `partial_v^r y(0)=0` for every r.
Differentiating (3) with respect to v now obeys EXACTLY (6), with no initial
boundary term. The same proof of (7) therefore covers every fixed mixed
parameter/pulse jet. The ACTUAL prescribed regular source's full smooth
temporal zero extension supplies precisely this stronger condition.
Only finitely many vanishing entrance jets are needed for any fixed order.
The constant-entrance control is not an admissible smoothly endpoint-
supported source and must not be used to deny this stronger conclusion.

Returning from `(xi,v)` to the common torus differentiates v with `O(S)`
factors. For entrance-compatible sources these retain the inverse-square
coefficient gain; for a general source pulse jets instead have sufficient
polynomial `1+m^2` costs. Physical differentiation of the oscillating
carrier still introduces powers of |m|. No infinite-harmonic physical
jet estimate or nonlinear summability is inferred from (7).

## 4. Pressure, support and the physical transfer

For high harmonics the fixed-parameter bounds implied by (2),(7) are

\[
 |D^I(t/P)|\le C_I/(m^2-1),\quad
 |D^I(\mathcal C/P)|\le C_I/[k|m|(m^2-1)],\quad
 |D^I(\pi/P)|\le C_I/(k|m|),\quad
 \mathcal C=\frac{i(n\mathbin\times t)}{km|n|^2}.           \tag{8}
\]

Constants include source bounds and fixed frame/edge costs. The pressure
estimate cannot be replaced by the potential estimate: take constant
`n=e3`, `f=e3 F(v)`, `t=0`. Then `pi=iF(v)/(km)`. F can be any nonzero
smooth compactly supported pulse. This exact admissible control retains
source-normal pressure and blocks an extra inverse-square pressure gain.
The moving-normal term in (2) is independently nonzero in the rational
fixture. The background radial flow and full Pi1/Pi2 are NOT set to zero
in supplying K and the same background balance.

If f vanishes on a full slow/transverse/shell parameter path, uniqueness
gives zero response there; differentiated zero data give smooth support
extension across those boundaries. Deck translations commute with the
path, so uniqueness gives common-torus descent without averaging distinct
preimages. Even a smooth source zero near both temporal endpoints can
leave a nonzero exit: the scalar response to a positive interior source
is a strictly positive damped integral. Gaussian smallness is not zero.
Apply the ORIGINAL cutoff to potential AND pressure before the complete
curl. Curls, graph derivatives and Cartesian frame derivatives do not
enlarge those closed supports. The uncut and cut support claims differ.

If the stipulated full-path source belongs to W_alpha, its inverse t is
W_alpha, pressure/potential W_(alpha+1/2), and complete curl remainder
W_(alpha+1/2-kappa). The actual complete cutoff curl is W_alpha.
For each FIXED finite harmonic set and physical derivative counts a,b,
the existing exact graph/frame transfer gives the respective Q powers

\[
 \delta u:\ -A+h\alpha-a s_x-b s_t,\quad
 \delta\pi:\ -2A+h(\alpha+1/2)-a s_x-b s_t,\quad
 \delta\mathcal A:\ 1/2-A+h(\alpha+1/2)-a s_x-b s_t,
                                                               \tag{9}
\]

where `A=1/2+h`, `sx=1/2+h/2`, `st=1+3h/2`. Edge weights remain
sqrt(zeta) with fixed inverse-distance losses; all fixed Gaussian tails
remain flat on the SAME domain by square completion. Sources designated
as additive flat tails are still not forwarded: a flat early input can
have a nonflat central response. No parameter, width or first band is
retuned. Constants are not uniform in derivative or nonlinear stage.

## 5. Checks, correction boundaries and next obligation

The new Fraction-only fixture checks 144 real principal-ODE components
for positive/negative modes, including moving normal, source normal and
complex reality. A separate noncommuting 2x2 numerical model checks (4)
and forced residence estimates; scalar integral diagnostics check fixed
parameter jets. These are manufactured controls, NOT actual-pulse
numerical evaluation or continuum supremum certificates. The analytical
proof is (1)--(9). The pulse-jet distinction is important: generic sources
need entrance compatibility; actual full smooth endpoint support supplies
it. Exact initial-jet recurrences test both cases. Pressure still has the
independent sharp obstruction described above.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
Output files must be fresh. Frozen source bindings and scope guards are
checked in addition to the exact and numerical controls.

Next audit the actual generated-source and complete-covariance mean
inputs against the current root/pressure/support foundations, rather than
reimporting a blanket collection of finite analytic lemmas. The original
positive primary roots, actual frame/background and retained exact core
remain dependencies; the new sourced inverse itself is now explicit from
them. Then verify all viscous/cutoff/interaction residual groups and the
single-field realization. No unconditional pulse-family theorem, complete
residual improvement or blow-up is established here.
