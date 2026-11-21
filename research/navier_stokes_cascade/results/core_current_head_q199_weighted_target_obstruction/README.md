# Q199: exact-moment weighted current-target obstruction

Parent: Q198, `10f9e41ccbcca8538dab3eb00674321085b50718`.
The live branch, Q198, both active CI runs, the current two-edge target
budget, generic signed covariance provider, finite support closure, original
root budget and Q184 logarithmic correction were read first. Historical
evidence, fields, roots, pressure, supports, accepted scales and the
unchanged `C^-20` activation width are preserved. Q198 CI had completed 152
steps without failure and was running step 153 when this checkpoint began.

Q198 correctly identified the missing weighted current-target inequality,
but its simple constant example did not respect two side conditions already
available in the construction: exact projected moment cancellation and
smooth compact interior support. This checkpoint corrects that weakness. It
shows that even those side conditions, together with small unweighted jets
of every fixed order, do not imply the required flat-edge estimate. This is
an obstruction to an inference from the present evidence, not a
counterexample to the intended canonical field.

## 1. The exact missing inequality

At the Q198 rates

\[
 B=1-\kappa,\qquad C=\tfrac32-\kappa,
 \qquad D=\tfrac32-2\kappa,\qquad C_0=2-2\kappa,     \tag{1}
\]

the generic covariance and projected-target algebra has conditional gain
`C0-C=1/2-kappa`. To feed that target through the original reciprocal-root
map, however, each angular/axial current row must satisfy, for every fixed
mixed jet `n`,

\[
 \sup_{\rm chart}\epsilon^{-C}\zeta^{-1}\delta^{p_n}
 J_n\!\left(R^{(e+1)/2}E_e\right)
 \le \widehat E_{e,n}(S),\qquad e\in\{1,2\},        \tag{2}
\]

on every accepted band. The exact moment projection, whole moving pressure
product, reciprocal-root derivatives, nonlinear chart conversion and the
expanded pressure/source pole remain part of (2).

The frozen dependency graph explicitly declines to certify (2): the target
budget marks actual weighted current rows false; the generic signed provider
marks its actual weighted root/chart/support inputs false; the support
closure marks the original weighted input contract false; and the root
budget marks the actual current target false. Their algebraic gates pass,
but no evidence path turns those conditional recipes into an actual-field
certificate.

## 2. A two-bump source with exact moment cancellation

Let `phi` be a nonnegative smooth bump supported in `(0,1/4)` with mass
`m>0`. For small `eta>0`, define

\[
 g_\eta(x)=e^{-a/(4\eta^2)}
 \left[\phi\!\left({x-\eta\over\eta}\right)
 -\phi\!\left({x-7\eta/4\over\eta}\right)\right].  \tag{3}
\]

The first bump is supported in `(eta,5eta/4)` and the second in
`(7eta/4,2eta)`. They are disjoint, smoothly compactly supported in the
interior, and have exactly opposite masses:

\[
 \int g_\eta=\eta m e^{-a/(4\eta^2)}
              -\eta m e^{-a/(4\eta^2)}=0.            \tag{4}
\]

At `x=3eta/2`, between the two supports, the oriented primitive contains
the whole first mass,

\[
 F_\eta(3\eta/2)=\eta m e^{-a/(4\eta^2)}.            \tag{5}
\]

The smooth bounded radial factor `R^{-e}` does not alter the edge order, and
the sources can be converted to the required radial moment by multiplying
by `2R^{-e-1}`.

## 3. Every fixed unweighted jet is small, but every fixed edge pole fails

Every fixed derivative order `j` has an unweighted bound

\[
 \|D^jg_\eta\|_\infty
 \le C_j e^{-a/(4\eta^2)}\eta^{-j}\longrightarrow0. \tag{6}
\]

Thus this family defeats neither smoothness nor any prescribed finite list
of unweighted seminorms. At the observation point, the retained flat weight
gives the exponential gap

\[
 {F_\eta\over\zeta(3\eta/2)}
 \ge \eta m\exp\!\left({7a\over36\eta^2}\right).    \tag{7}
\]

For any fixed inverse-distance pole `d`, choose an integer
`k>(d+1)/2`. The elementary bound `exp(y)>=y^k/k!` yields

\[
 \eta^d{F_\eta\over\zeta(3\eta/2)}
 \ge c_{a,m,k}\eta^{d+1-2k}\longrightarrow\infty.  \tag{8}
\]

Therefore no fixed polynomial edge pole and no constant independent of
`eta` follows from smooth compact support, (4), and all fixed unweighted
jet bounds. A finite-regularity polynomial bump is used only to verify the
positive mass and rational geometry in the executable diagnostics; the
analytical obstruction uses a smooth bump.

## 4. What this corrects and what it does not prove

Q198's constant function already showed that an unweighted bound alone does
not imply division by `zeta`, but it lacked exact zero moment and compact
interior support. Equations (3)--(8) close that logical gap. They do not show
that the actual canonical current rows fail (2). Those rows have additional
structure from complete pressure, fast/slow maps, wave factors and support
identities; that structure must be used directly to prove (2).

The conditional Q198 orders remain unchanged. In particular, the projected
target is still conditionally `C0`, while Haar pressure/scalar changes remain
`C`, pointwise pressure remains `D`, and the physical zero-jet regular
residual remains

\[
 -2A-\tfrac12+hD=-{15025001\over10000000}<0.        \tag{9}
\]

Passing the exact tests does not supply a continuum enclosure or physical
decay.

## 5. Remaining boundary

The current evidence cannot justify constructing the actual next signed
pulse, so advancing to post-signed fast/slow means as though (2) were proved
would be circular. The next obligation is to derive (2) from the complete
current pressure/fast/slow equations and original weighted wave inputs on
every accepted band, retaining the expanded pressure pole and whole-path
support. Only then may the original-root signed correction be declared
actual and the fast/slow mean relaxation recomputed.

No unconditional pulse bound, complete residual improvement, admissible
force, comparison exclusion, infinite realization or blow-up is claimed.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q199_weighted_target_obstruction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q199_weighted_target_obstruction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q199_weighted_target_obstruction/evidence.json \
  --output q199-fresh.json
```
