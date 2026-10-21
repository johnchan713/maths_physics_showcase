# One accepted primary tail and quantitative all-band cutoff flatness

Parent: `9f62f1c849779bf3516253aa42b4d8e43bd571cd`.
Status: `same-profile-primary-tail-value-gates-and-dyadic-edge-flatness-rechecked`.

The latest coefficient checkpoint preserves the same repaired leading field,
finite coefficient supports, nonzero leading radial flow, and both exceptional
reserved pressure orders. Existing work already supplies primary pulses on a
separately certified tail, relative to the named analytical foundations. This
checkpoint assembles its finite acceptance conditions into one integer selector
and makes the cutoff-error flatness constant explicit after edge absorption
and summation over every active band. It does not repeat the primary ODE or
claim verification of the complete nonlinear pulse family.

The selected coefficient cutoff, tilt `u^2=2^19 G^22`, raw deficit
`d0=1/(2^19 G^22)`, pulse radius `r0=1/6057072`, width `C^-20`, roots,
carriers, gauges, and historical evidence remain as recorded. The new tail is
a certified observation/definition domain; it is not identification with an
unknown historical first band or an earlier localized field and force.

## 1. Actual finite budgets, with one threshold before differentiation

Write `Q=2^-ell`, `epsilon=Q^h`, `S=ell^2`, and
`k=ceil(epsilon^-1/2)`, where `0<h<1/100`. Use the positive value budgets
from the retained leading/frame/cutoff lemmas: `Ccoeff, rb, omega, Lpsi,
Cmesh, Cpert, Dbox, dchart`. Set

\[
 \gamma=\omega^2d_0/(9\kappa_u),\qquad
 B_0=\min\{r_b/8,\omega^2d_0/(16L_\Psi),\gamma/(128C_{pert})\},
 \quad m=\lceil\log_2(C_{coeff}/B_0)\rceil_+.             \tag{1}
\]

Here `Ccoeff=Cstar(729,2916)` retains the actual selected cutoff and the
finalized leading/first/second coefficient inputs. Let `Lreq` be the finite
maximum of the normal, moving-frame, Riccati, covariance-cone and any requested
root-comparison lengths in the pulse/value foundations. They include
`25`, `8CE/mu`, and the positive Haar-column/cone thresholds. A comparison
length may be retained although existence of the canonical roots does not
require comparison with historical roots.

Let `qstar<=1` be a valid background domain endpoint. Pad the slow-coordinate
error by a finite value constant `Cslow>=1`, with required tolerance
`0<tau<=1`. No manufactured fixture supplies these actual constants. For the
fixed parameters the inherited lemmas give finite constants; the formulas
below are real-variable estimates, with exact rational outward bounds usable
in the implementation. Define

\[
 \ell_* = \max\left\{
 16,\ 2+\left\lceil\frac m{2h}\right\rceil,
 \left\lceil\max\{1,16C_{mesh}/\gamma,D_{box}/d_{chart}\}^{1/6}\right\rceil,
 \left\lceil\sqrt{L_{req}/(2r_0)}\right\rceil,
 3+\lceil\log_2(1/q_*)\rceil_+,
 \left\lceil\left(\frac{3\,8!\,4^8C_{slow}}{\tau h^8}\right)^{1/4}\right\rceil
 \right\},\qquad q_*=q_{star}.                         \tag{2}
\]

All gates hold for EVERY integer `ell>=ell_*`. In particular
`floor(2h(ell-2))>=m` proves the mixed value-error gates from (1), while
`Cmesh ell^-6<=gamma/16` supplies the positive current cone budget.
The length bound keeps the true radius factor `Ls>=2r0 ell^2`.
The background condition gives `4Q<=qstar/2<qstar` and the mesh/chart row
keeps the independent enlargement constraint.

The slow gate is an all-half-line estimate, rather than a sampled endpoint:
`k^-1<=sqrt(epsilon)`, `log(2)>=1/2`, and the eighth exponential-series
term give

\[
 C_{slow}\ell^4(\epsilon+\epsilon^2+k^{-1})
 \le3C_{slow}\ell^4e^{-h\ell/4}
 \le\frac{3\,8!\,4^8C_{slow}}{h^8\ell^4}\le\tau.       \tag{3}
\]

The selector uses integer ceilings, roots and the required binary exponent;
it never constructs `2^ell`, `2^floor(2h(ell-2))` or tiny `q`. This matters
even in the rational control `h=2^-200`, where the scale exponent itself is
enormous. No actual giant field or norm enclosure is numerically evaluated.
There is no derivative order in (1)--(3). Fixed higher jets change constants,
not this common accepted tail.

## 2. Same profile, common physical domain and primary estimates

Define

\[
 q_v=\min\{q_{star}/64,\ 2^{-\ell_*-2}\}.                \tag{4}
\]

At `0<q<=qv`, an active enlarged chart has `q/Q` in `[1/4,4]`. Hence
`Q<=4q<=2^-ell_*` and `ell>=ell_*`. Its whole enlargement satisfies
`q'<=4Q<=16q<=qstar/4`. The background, source memory and value gates apply
throughout that chart. This supplies the required common domain without
assuming acceptance of earlier bands. Discrete carriers, representatives,
signs, covers and copies are frozen before differentiating.

On this domain the retained frame and Riccati argument gives the positive
uncut radial amplitude `cP<=x<=CP`; the positive covariance roots retain their
actual transverse Haar area. With the inherited complete curl, graph and
Cartesian-frame terms, the primary physical losses are

\[
 s_x=1/2+h/2,\qquad s_t=1+3h/2.                          \tag{5}
\]

For each fixed spatial/time derivative pair `(a,b)`, uniformly over accepted
bands, labels and copies, the velocity, wave pressure, and fixed potential
have respectively the powers

\[
 Q^{-(a+1)s_x-bs_t},\qquad Q^{-(a+2)s_x-bs_t},
 \qquad Q^{-as_x-bs_t},                                \tag{6}
\]

times the retained polynomial `S` and `sqrt(zeta) delta^-J P` factors.
These are conditional transfers of the reviewed primary theorem, not a new
construction of all sourced corrections. No positive pointwise lower bound
for an oscillating/cutoff velocity follows from `x>0`.

The potential and pressure cutoffs precede the complete curl. Derivatives do
not enlarge their closed supports: `|xi|<=r0/2`,
`|v/Ls-1/2|<=2/5`, and the same radial shell. The source-normal and
moving-normal pressure terms, leading radial velocity and reserved Pi1/Pi2
remain. The latter keep the derivative and pressure-error formulas of the
parent checkpoint; leading pressure alone does not replace them.

## 3. An explicit retained edge factor

Put `ca=t1^2/2>0`, `cb=2`, so
`sqrt(zeta)=exp(-ca/ya^2-cb/yb^2)` with `t1=C^-20`. For a fixed integer
edge power `J>=0`, let `j=ceil(J/2)` and

\[
 E_J=3+j!\{(2/c_a)^j+(2/c_b)^j\}.                       \tag{7}
\]

Since `delta^-J<=1+ya^-J+yb^-J`, `y^-J<=1+y^-2j`, and
`x^j exp(-x)<=j!`, splitting each edge exponential in half proves

\[
 \delta^{-J}\sqrt\zeta
 \le E_J\exp[-c_a/(2y_a^2)-c_b/(2y_b^2)].               \tag{8}
\]

Thus one half-edge factor absorbs every fixed finite polynomial edge loss
and leaves a flat factor for zero extension at either shell boundary.
For derivatives apply (8) at their actual finite enlarged edge exponent.
The narrow inner factor is explicit: `(2/ca)^j=(4/t1^2)^j`; it cannot be
silently replaced by an order-one constant uniform in the setup width.

## 4. Cutoff errors summed over all active physical bands

For an error ALREADY proved by the pulse/source lemmas to obey the full-path
envelope AND to retain the closed enlarged label supports, the cutoff region
gives, after any fixed physical derivative,

\[
 |e_\ell|\le C\ell^B Q^{-\Lambda}
               \delta^{-J}\sqrt\zeta\,e^{-a_0\ell^2},
 \qquad a_0=c_G r_0/25>0.                              \tag{9}
\]

Polynomial `S` losses have been included in the integer `B`. For each desired
integer flatness order `K>=0`, set `D=B+max(0,Lambda+K)`. On the SAME whole
domain `ell>=ell_*>=1`, `log ell<=ell` and `log2<1` give

\[
 \ell^B Q^{-(\Lambda+K)}e^{-a_0\ell^2}
 \le e^{D\ell-a_0\ell^2}\le e^{D^2/(4a_0)}.             \tag{10}
\]

The exact identity behind (10) is a nonnegative completed square; it imposes
no new derivative-dependent starting band. Combining (8)--(10) bounds each
error by `C E_J exp(D^2/(4a0)) Q^K`, with half-edge flatness still present.

At fixed `q>0`, the enlarged membership condition implies
`q/4<=Q<=4q`. Such a dyadic interval contains at most FIVE integer bands,
including its closed endpoints. The geometry's separated auxiliary supports
give at most one actual label at a point; even allowing one contribution
per enlarged band is a valid conservative bound. Consequently

\[
 \sum_{\ell\ {\rm active}}|e_\ell|
 \le5\,4^K C E_J e^{D^2/(4a_0)}q^K.                    \tag{11}
\]

There is no number-of-bands growth near `q=0`. Derivatives retain the same
closed supports and satisfy (11) with their own finite `B,J,Lambda,C`.
This proves all-fixed-jet flatness of this additive cutoff-error ledger
on a common tail, conditional on (9) and the stated support hypothesis. It
does not turn an arbitrary source
into an envelope-controlled source. Time-boundary jets remain one-sided
at `T=0`; no global temporal zero extension is inferred.

## 5. Exact frontier and reproduction

The unclosed input is now specific: every generated full harmonic source
must satisfy its full lifted-path envelope, smooth endpoint extension,
edge/jet bounds and required mode moments. A finite sourced inverse mapping
does not prove that membership. Complete-wave covariance means must also
retain axial velocity, curl remainders, signed squares and both exceptional
pressure orders. Next independently join that complete covariance/mean map
on this common tail, then audit generated sources and the complete
viscosity/cutoff/interaction residual. Old force identification and blow-up
remain unproved.

Reproduce here with `python3 test_audit.py -v` and
`python3 audit.py --verify-record evidence.json --output fresh.json`.
Exact rational controls exercise coupled gate failures, huge-exponent-safe
selection, fixed dyadic overlap, square completion and retained width costs.
They are not continuum norm computations. Historical records are preserved.

Local validation on 2026-10-21: 209 tests passed across this checkpoint (24),
physical-pulse foundation (37), sourced-inverse foundation (34), explicit
tilt/tail (12), auxiliary geometry (33), coefficient-support recheck (24),
and canonical value-domain selection (45). Their frozen evidence records
were verified by those suites. The numerical dependency suites used isolated
mpmath 1.3.0 dependencies; test counts do not certify analytical premises.
