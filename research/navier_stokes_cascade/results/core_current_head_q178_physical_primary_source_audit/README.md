# Q178: physical primary-pulse source audit

Parent: Q177, `9cd8abaa2c69e8b6108dcafdcf22f9e80b392563`.
The live branch, Q177, Q148/Q149, `physical_pulse_transfer`,
`physical_pulse_foundation`, frozen evidence and current CI were read first.
Historical evidence, profiles, roots, pressure gauge, frequency, supports and
the `C^-20` width are unchanged.

Q177's next obligation was to audit the physical-pulse theorem and its Q149
same-profile join before moving into covariance means. This checkpoint closes
that source audit. It does not reprove the retained exact-core, existence,
inverse or gap lemmas on which the physical theorem explicitly depends.

## 1. What the sources establish

Q149 correctly replaces Q148's proposed external numeric diagonal gate by
the rule-defined Borel schedule. For each fixed coefficient order `n`, all
finitely many retained seminorms `C_(n,m)` exist before `a_n` is chosen. The
choice therefore enforces

\[
 C_{n,m}q^{nh}\le 2^{-n}e^{-nh}<2^{-n}
 \qquad(q\le a_n^{-1},\ 0\le m\le n).                 \tag{1}
\]

There is no separately frozen numerical scale list to certify. For each
fixed jet, the same completed background has the exact tail

\[
 \sum_{n\ge N}2^{-n}(q^h)^n
 ={(q^h/2)^N\over1-q^h/2}\le2^{1-N}.                 \tag{2}
\]

On this same finalized profile, let `Q=2^-ell`, `S=ell^2`,

\[
 s_x={1\over2}+{h\over2},\qquad s_t=1+{3h\over2}.
\]

The audited physical theorem gives one later band `ell0` which works for
every later used scale and every fixed physical derivative order. For fixed
`a,b`, the complete primary fields satisfy

\[
 |\nabla^a\partial_t^b u_p|
 \le C_{a,b}Q^{-((a+1)s_x+bs_t)}S^B\sqrt\zeta\,\delta^{-B}P, \tag{3}
\]
\[
 |\nabla^a\partial_t^b p_p|
 \le C_{a,b}Q^{-((a+2)s_x+bs_t)}S^B\sqrt\zeta\,\delta^{-B}P, \tag{4}
\]
\[
 |\nabla^a\partial_t^b {\cal A}_p|
 \le C_{a,b}Q^{-(as_x+bs_t)}S^B\sqrt\zeta\,\delta^{-B}P.    \tag{5}
\]

These are uniform scale-dependent bounds, not absolute smallness. They keep
the complete curl, axial and radial velocity, graph coefficients, moving
frame, cylindrical terms, pressure cutoffs and the pressure formula

\[
 p_m={i\over km|n|^2}
 \bigl(n\mathbin{\cdot}K_{bg}t_m-n'\mathbin{\cdot}t_m
       +n\mathbin{\cdot}f_m\bigr).                    \tag{6}
\]

Thus the moving-normal and source-normal pressure terms are present. The
closed supports `|xi|<=r0/2` and `|v/L_s-1/2|<=2/5` are preserved by curl;
the zero extension is smooth on the accepted one-sided physical domain.

The lower bound `cP<=x<=CP` applies to the uncut radial ODE amplitude. It is
not a positive pointwise lower bound for the oscillating, cutoff velocity.
This distinction is explicit in the foundation evidence and is preserved
here.

## 2. Exact input boundary

The sourced inverse is a mapping theorem for a stipulated input class with
explicit bounds, harmonic support and smooth endpoint extension. The primary
bounds (3)--(5) do not by themselves prove that every generated nonlinear
source belongs to that class. In particular, a positive Gaussian envelope at
an endpoint is not a zero-extension theorem. Bounds and full endpoint
supports for all generated sources remain an analytical obligation.

This is the precise conditional result: uniform values, every fixed physical
derivative and closed supports are established for the actual **primary**
pulse family on all scales used by the construction, relative to the retained
exact-core and finite lemmas. No complete generated-pulse hierarchy is
established.

## 3. First covariance obstruction

The reserved mean pressure retains

\[
 \Pi_1=K_1(\eta)-{\mu^2\over4X},\qquad \Pi_2=K_2(\eta), \tag{7}
\]

and the complete wave retains axial velocity, radial velocity, curl
remainders, moving pressure and cutoff derivatives. Before their exact mean
cancellation, the known primary covariance remainder has the power

\[
                 Q^{-1-2h\kappa},\qquad \kappa=10^{-5}. \tag{8}
\]

For `0<Q<1`, the exponent `-1-2h*kappa` is negative. Equivalently, at
`Q=2^-ell`, (8) equals `2^((1+2h*kappa)ell)` and grows strictly with `ell`.
No finite constant in (3)--(5) turns it into absolute decay. The next step
must audit the complete covariance and mean identities, beginning with Q150,
and retain every tensor entry, both velocity components, curl remainders,
(7), moving pressure and cutoff terms.

The exact arithmetic tests check (1)--(5), (8), provider flags and frozen
provenance. The inherited continuum theorems are analytical premises; test
success is not an independent continuum proof. Complete covariance means,
all generated source inputs, the viscosity/cutoff/all-interaction residual,
force exclusion and blow-up remain unproved.

Run `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output fresh.json` here.
