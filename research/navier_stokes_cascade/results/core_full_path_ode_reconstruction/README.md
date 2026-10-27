# Canonical full-path ODE on the reconstructed same-band tail

Parent: `e016f4a705ef8a04c89478b8707cd8a07273a112`.
Status: `canonical-full-path-ode-reconstructed-on-current-same-band-tail`.

The parent reconstructed the value, spectral, factored-cone and phase-normal
gaps on one accepted tail of the finalized leading profile and real Borel
background. This checkpoint reconstructs the third primitive contract on that
same tail: finite-interval existence, the Gaussian-adapted propagator, every
fixed mixed jet, the canonical entrance, complete pressure/potential formulas,
and cutoff/support placement.

It retains an essential distinction. The canonical primary is homogeneous,
so its source is exactly `f=0`. The adjacent zero-entry sourced inverse uses
the same operator and must retain the unprojected `n dot f` pressure source.
Inserting a source into the primary or projecting the source before pressure
reconstruction are both wrong. Historical evidence remains frozen. Exact
controls check algebra and provenance; they do not sample the selected field.

## 1. Constrained equation and complete pressure

Freeze a label and its full lifted interval `0<=v<=L`. Discrete carriers,
representatives, signs, length and entrance are fixed before differentiating
continuous parameters. For the full background phase put

\[
 n=\nabla_*\Phi,\qquad d=\epsilon k^2|n|^2.
\]

For every nonzero harmonic `m`, the exact constrained system is

\[
 t'+Kt+m^2dt+ikm n\pi=-f,\qquad n\cdot t=0.             \tag{1}
\]

Differentiating the constraint and taking the scalar product of (1) with
`n` gives

\[
 \pi={i\{n\cdot Kt-n'\cdot t+n\cdot f\}\over km|n|^2}. \tag{2}
\]

Substitution yields

\[
 t'=A_nt-m^2dt-\Pi_{n^\perp}f,
 \quad A_n=-K+{n(n^TK-(n')^T)\over|n|^2},              \tag{3}
\]

and `(n dot t)'=-m^2d(n dot t)`. Thus tangent entrance data remain tangent,
and (2)-(3) are exactly equivalent to (1). A purely normal source has zero
projected velocity forcing but nonzero pressure. This is the explicit
obstruction to deleting `n dot f`; the moving-normal term is equally required.

The homogeneous canonical primary takes `f=0`, `|m|=1`, and its original
reference-mixing entrance

\[
 t(0)=P(0)B_r(0)e_+.                                   \tag{4}
\]

It is not identified with the distinct inherited actual-mixing primary. Its
principal covariance roots are rebuilt from this canonical pulse.

## 2. Same-tail denominator and propagator estimates

The parent gives on every complete path

\[
 |n|^2\ge(31/32)^2\beta_{phase}>0.                     \tag{5}
\]

The first contract supplies every fixed background jet, so `n,n',K,d` and
their fixed derivatives are finite on the closed interval. The original
moving frame satisfies

\[
 B^\ell(A_nB-B')=\operatorname{diag}(\lambda,-\lambda)+E,
 \quad\|E\|+|d-d_{ref}|\le H/L.                        \tag{6}
\]

With `P'/P=lambda-dref`, the Hermitian energy estimate gives on every
subinterval

\[
 \|V_m(v,w)\|\le T_0{P(v)\over P(w)}
 e^{-(m^2-1)d_*(v-w)},\quad0\le w\le v\le L,           \tag{7}
\]

where `d_*>0`. The complete moving frame, both endpoint conversions and
the actual variable rate are retained. No numerical lower bound on `P` is
needed: its factors cancel between entrance/source terms and (7).

For the Gaussian-adapted ambient lift set `h=n/|n|^2`, `H_n=hn^T`,
`A=A_n-dI`, and

\[
 \widehat A=A+(d+P'/P)H_n.                              \tag{8}
\]

It agrees with `A` on tangent vectors and obeys

\[
 (n')^T+n^T\widehat A=(P'/P)n^T.                       \tag{9}
\]

This prevents a reciprocal Gaussian endpoint loss in normal extensions.

## 3. Existence and every fixed jet

Continuous finite-dimensional coefficients on `[0,L]` give Picard
differences bounded by `C(ML)^j/j!`. Their successive ratio is
`ML/(j+1)`, eventually at most `1/2` for every finite `M,L`; no contraction
smallness or source-size assumption is needed. The same estimate proves
uniqueness.

Differentiate (3) at a fixed mixed order. Each derivative satisfies a
linear equation whose forcing contains only finite coefficient jets and
lower solution jets. Duhamel with (7) therefore gives, inductively,

\[
 |D^It(v,x)|\le P(v,x)U_{|I|}(S),                       \tag{10}
\]

with a finite constant for every fixed order. The constants may be enormous
and are not uniform in the derivative order. They require no smaller
value-domain tail because existence and finite-jet induction use boundedness,
not a new perturbative smallness condition.

For sourced modes the residence factor is

\[
 \ell_m=\begin{cases}L,&|m|=1,\\
 \min\{L,[(m^2-1)d_*]^{-1}\},&|m|\ge2.
 \end{cases}                                            \tag{11}
\]

Differentiated damping terms `(m^2-1)D^Jd` remain in the jet recurrence.
Entrance-compatible sources extend (11) to fixed pulse-coordinate jets.
No infinite-harmonic summability or uniform nonlinear source bound is
asserted here.

Pressure (2) and potential

\[
 {\cal A}={i(n\times t)\over km|n|^2}                  \tag{12}
\]

inherit the same fixed-jet bounds. The local wave pressure is retained in
full. The distinct `Pi1/Pi2`
background-pressure memory and moving pressure moments are not inputs to this
isolated solve: they remain explicit later mean-ledger obligations and are not
reset or declared negligible here.

## 4. Cutoffs and support

The raw homogeneous solution is not compactly supported in `v`, and an
interior sourced response may remain nonzero at the outgoing endpoint. The
construction estimates the raw solution on the entire interval, then applies
the original temporal and transverse cutoffs to potential and full pressure:

\[
 \operatorname{curl}(\psi{\cal A})
 =\psi\operatorname{curl}{\cal A}+\nabla\psi\times{\cal A},
\]

\[
 \nabla(\psi\pi)=\psi\nabla\pi+\pi\nabla\psi.           \tag{13}
\]

Every cutoff derivative in (13) is retained. If a source vanishes on an
entire slow/transverse path, uniqueness gives zero response there and smooth
support extension across that boundary. Vanishing only at the temporal
entrance does not imply zero transverse support or zero raw exit.

## 5. Result and next obligation

Equations (1)-(13) reconstruct the third primitive contract on the current
same-band tail, relative to the frozen background/frame continuum providers.
They establish definitions and fixed-jet/support bounds for the primary and
adjacent sourced inverse, not complete residual smallness. No unconditional
pulse family, admissible force, exclusion or blow-up theorem follows.

The next primitive contract is the fixed compensator/projector family: prove
its actual supports and every fixed jet on this same domain. Then join the
compact and slow inverse weighted-domain estimates. Complete covariance,
fresh axial/full-pressure means, viscosity, cutoffs and all interactions
remain later obligations.
