# Leading core pressure: identification, endpoint recovery and retained memory

Parent: `0b16f1b18aa6277dc507e39fd7cecbd9a6e4016a`.
The parent CI was still running on 18 October 2026, with no failed step
reported at inspection. Historical fields, sources and frozen evidence stay
unchanged. No unconditional construction or blow-up result is claimed.

The requested dependency trace resolves a scope issue: the SAME-ROOT
leading pressure identity was already supplied by
`global_pressure_moment_foundation`, following `shape_join_pressure` and
`global_moment_stress`. Its retained core input is not an assumption that
the full corrected Navier-Stokes pressure equals centrifugal pressure.
This checkpoint explicitly joins that leading datum to the analytic core,
derives a quantitative angular inverse for a pressure-normalization defect,
and retains the nonzero higher-order background pressure memory.

The new analytical result is the inverse and its endpoint/range condition
in sections 3--4. A new exact control detects the unbounded homogeneous
mode that would invalidate the inference without endpoint regularity.
Passing its finite controls is not a continuum proof.

## 1. Define the reference pressure before selecting the radial scale

For the UNEDITED reference, use y=log(X/XR) and
E_id,sched(y,eta)=c(y) f(eta)^theta(y), f=(1+eta²)^-1, 0<=theta<=1.
The coefficients and logarithmic segment lengths are independent of eta
and XR. In particular the exact normalization is

\[
\Pi_{axis}(\eta)=-\frac12\int_{-\infty}^{\infty}E_{id,sched}(y,\eta)^2dy
             =-\int_0^\infty\frac{E_{id}(X,\eta)^2}{2X}dX.       \tag{1}
\]

The change of variable dX/X=dy cancels XR exactly, even when XR is later
chosen from the core attachment parameters. There is no pressure feedback
through that choice. Dropping 1/X in (1) would destroy this invariance.
P and the other fixed schedule parameters are not varied in this argument.

The pushforward of c(y)²dy/(2P²) under y->theta(y) is a positive measure
mu on [0,1]. The exact inner power contributes (5/2)delta_1. The reference
mass estimate in `core_growth_activation`, derived before any quadrature,
gives mu([0,1])<5. Therefore the analytic core receives precisely

\[
\Pi_{axis}=-P^2\int_0^1(1+\eta^2)^{-2\theta}d\mu(\theta),
\qquad \mu\ge(5/2)\delta_1,\quad\mu([0,1])<5.                 \tag{2}
\]

This is the defining measure, not the discrete measure used by diagnostic
quadrature. The common logarithm and holomorphic bounds used by
`core_contraction_audit` thus concern this very datum. The infinite-space
core construction and its original-equation recovery remain the retained
local theorem; finite Taylor coefficients are not substituted for its limit.

Source: retained manuscript Lemma A.5, (A.21)--(A.23), Appendix B and
(4.25), SHA256
`0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f`.
The supplied manuscript's final theorem is not a premise.

## 2. The completed leading profile has the SAME axis normalization

Let E0 denote the leading swirl after the accepted angular, axial, heat,
annular and compact repairs. These are the accepted leading roots, not the
different covariance roots of the subsequently rebuilt physical primary.
The same-root transformed five-row systems have the original zero set
(M,R,J,S,Cp), with the retained lambda>0. The complete edit ledger gives

\[
\Pi_{axis}+\int_0^\infty E_0^2/(2X)dX=0,
\qquad
\Pi_0(X,\eta)=\Pi_{axis}+\int_0^X E_0^2/(2x)dx
            =-\int_X^\infty E_0^2/(2x)dx.                     \tag{3}
\]

The heat edit and its compensation can have nonzero memory between their
supports. Only their TOTAL increments cancel. Neither the core gauge nor
incoming pressure is reset. Later physical-primary/covariance rebuilding
does not edit E0 or these leading roots. Thus (1)--(3) identify the core's
datum with the completed LEADING backward pressure in the same gauge.
The earlier local scope warning is not an additional unsolved leading
pressure equation after this trace. It also does not prove every hypothesis
of the accepted leading maps or the later corrected construction.

All fixed eta derivatives in (3) are meaningful, including one-sided
endpoint derivatives. Near the axis E0²/(2X) is bounded and smooth. On
finite radial intervals the accepted maps supply smooth parameter jets.
In the unchanged exact heat tail, the Gamma estimates of
`annular_exterior_foundation` give, for each fixed k, an integrable bound
|partial_eta^k(E0²/(2X))|<=C_k X^(-2-2h), h>0 fixed. Explicitly,

\[
\int_B^\infty X^{-2-2h}dX=\frac{B^{-1-2h}}{1+2h},\quad B>0.   \tag{4}
\]

Dominated differentiation proves each fixed jet of (3). This is not a
claim of complex analyticity of the heat endpoint; those heat Taylor
series can have zero radius. There is no uniform-in-k constant or exchange
of h->0 with the other renormalized moment integrals.

## 3. New quantitative recovery from the axial pressure-stress defect

For ANY candidate leading swirl with the same forward datum, write its
functional normalization defect as

\[
\delta(\eta)=\Pi_{axis}+\int_0^\infty E_0^2/(2X)dX
            =\Pi_{forward}(X,\eta)-\Pi_{backward}(X,\eta).      \tag{5}
\]

The last difference is independent of X. In the preserved leading stress
identity, its exact contribution is
Delta(rT_z)=X[4A eta delta-d delta']/L, where the other moment contributions
are held fixed or subtracted, A=1/2+h,
d=1-eta² and L=1-2h eta². Keep the derivative and the factor X. Define
g=K_A delta=4A eta delta-d delta' and Gamma=g/L. The following lemma
is independent of any fluid construction.

For A>=1/2 and a bounded smooth delta on [-1,1], put p=2A. Then

\[
(d^p\delta)'=-d^{p-1}g,\qquad
\int_{-1}^1d^{p-1}g\,d\eta=0,                                \tag{6}
\]
\[
\delta(\eta)=d(\eta)^{-p}
 \begin{cases}\displaystyle\int_\eta^1(1-s^2)^{p-1}g(s)ds,&0\le\eta<1,\\
 \displaystyle-\int_{-1}^\eta(1-s^2)^{p-1}g(s)ds,&-1<\eta\le0.
 \end{cases}                                                 \tag{7}
\]

The endpoint products d^p delta vanish because delta is bounded and p>0.
This proves (6)--(7); the second relation in (6) is the compatibility
condition for the two endpoint inversions. Conversely, any smooth g with
that weighted integral zero defines a unique smooth delta by (7).
For the right endpoint set s=eta+(1-eta)t. Formula (7) becomes

\[
\delta(\eta)=(1+\eta)^{-p}\int_0^1(1-t)^{p-1}
 [1+\eta+(1-\eta)t]^{p-1}g(\eta+(1-\eta)t)dt.                  \tag{8}
\]

For 0<=eta<=1 every fixed differentiated integrand is dominated by a
constant times (1-t)^(p-1), which is integrable. Reflection handles the
left endpoint. Compatibility matches the values at zero; the ordinary
first-order equation there matches every higher derivative. In particular
delta(1)=g(1)/(4A), delta(-1)=-g(-1)/(4A).

For eta>=0, use 1+s<=2 in (7). Reflection gives the other half, so
||delta||infinity<=2^(p-1)/p * ||g||infinity. In the actual core range
1/2<=A<=.51, convexity gives 2^(p-1)<=p (the chord on 0<=p-1<=1).
Consequently

\[
\boxed{\|\delta\|_\infty\le\|K_A\delta\|_\infty
                        \le\|\Gamma\|_\infty.}                \tag{9}
\]

The last inequality uses 0<L<=1. In particular zero axial pressure-stress
defect on the WHOLE angular interval implies exact pressure normalization.
A single zero value, a small norm, or a finite sample does not imply it.

## 4. All fixed angular jets and the indispensable endpoint condition

For k>=1 exact differentiation, including both derivatives of d, gives

\[
K_{A+k/2}(\delta^{(k)})
 =g^{(k)}-k(4A+k-1)\delta^{(k-1)}.                            \tag{10}
\]

The same inverse proof has constant 2^(2A+k-1)/(2A+k)<=2^k on the original
A range. Thus with B0=G0 and

\[
G_k=\|\Gamma^{(k)}\|_\infty+4hk\|\Gamma^{(k-1)}\|_\infty
                    +2hk(k-1)\|\Gamma^{(k-2)}\|_\infty,
\quad B_k=2^k[G_k+k(4A+k-1)B_{k-1}],                         \tag{11}
\]

negative-index norms zero, induction gives ||delta^(k)||infinity<=Bk.
Every fixed bound is finite. These are not uniform-in-order or physical
scale bounds and do not make a nonzero defect exactly zero.

Endpoint regularity is ESSENTIAL. If it is removed, delta=d^(-2A) is a
nonzero homogeneous solution K_A delta=0 on (-1,1). Its weighted endpoint
trace is one, so (6)--(7) cannot be invoked. For A=1/2 this is the rational
fixture 1/(1-eta²), unbounded as eta tends to either endpoint. Similarly
g=1 fails the compatibility integral in (6), so there is no bounded
solution matching both endpoints. These are exact ODE controls, not fluid
counterexamples. The leading construction's smooth closed angular traces
supply the needed endpoint condition relative to its accepted map inputs.

## 5. Leading identification must NOT erase corrected pressure memory

On the actual reserved I4 patch the coefficient foundation gives, exactly,

\[
q^{2A}p_B=\Pi_0+\chi(c_1q)q^{2h}[K_1(\eta)-\mu^2/(4X)]
                         +\chi(c_2q)q^{4h}K_2(\eta),\quad\mu>0. \tag{12}
\]

Every positive-order background velocity vanishes there, but the first two
pressure orders remain. At fixed q,eta the radial derivative of the
background pressure difference from Pi0 is
chi(c1 q)q^(2h) mu²/(4X²), strictly positive when its cutoff is one.
Thus even the pre-wave background pressure is not generally the leading
centrifugal pressure. K2 also retains later-source memory despite its zero
local radial derivative. This is a statement about the same background
construction, not a claim that a manufactured K2 value is its actual value
or a substitution for the complete canonical pressure after waves/means.
All those pressure terms remain in the canonical force and comparison.

## 6. Result, limits and the next concrete acceptance gate

The leading pressure-family/global-backward-pressure compatibility follows
from the actual definition (1), the same-root ledger (3), fixed-order tail
control and the regular angular gauge. It was already implicit in the
current leading checkpoint, not a new full-fluid existence theorem. The
new inverse (7)--(11) makes its defect criterion quantitative and detects
the exact endpoint obstruction. Historical scope statements are preserved
and interpreted at their original checkpoint.

The next audit should address actual retained value/defining-map premises,
starting with the canonical primary's common root-value acceptance gate:
theta=12 C_H/[j_a^- min(Amin,u) L] <=1/(2kappa+1),
kappa=24 j_a^+/(gamma j_a^-), together with its common normal/Riccati/cone
thresholds on the original accepted domain. The diagnostic budgets use
manufactured inputs; they do not by themselves verify those constants at
the actual first accepted band. Establish the inequality and the permitted
domain choices analytically, or record the exact obstruction, before an
unconditional pulse or complete residual claim. No acceptance failure is
asserted here. Further repetition of completed comparison implications is
not a substitute for this check.

The Fraction-only controls check the differentiated operator, weighted
compatibility, two endpoint inversions, singular mode, radial-scale
Jacobian, stress-defect scaling and higher pressure memory. They are finite
algebraic fixtures, not the actual infinite core or global fluid. The
analytical argument above carries the continuum quantifiers.

```sh
python3 research/navier_stokes_cascade/results/core_pressure_identification/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_pressure_identification/audit.py \
 --verify-record research/navier_stokes_cascade/results/core_pressure_identification/evidence.json \
 --output /tmp/core_pressure_identification_fresh.json
```
