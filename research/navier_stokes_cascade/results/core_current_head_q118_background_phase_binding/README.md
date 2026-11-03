# Q118: current full-background phase, rate and linear-coefficient binding

Parent: `21d2c754f19eeab8888209321d7423aeeedffeab` (Q117).
On 3 November 2026 the inspected remote was Q117; research CI 41587075436
and field CI 41587075430 were running. Q116's two runs were also running.
The latest notes, frozen evidence and the local leading-core/first-collar
scope at width C_setup^-20 were read. This checkpoint connects EXISTING
full-background phase formulas to Q93--Q117's current source graph. It
retains pressure memory and pays the coordinate conversions explicitly.
Historical frozen evidence, roots, fields, profiles and selectors are unchanged.

## 1. Reuse the actual full-background coefficient recipe

core_phase_coefficient_budget already derives the frozen phase and rates
from FULL F,G,btilde, the original path/carrier bounds and leading value
gaps. core_background_finite_head_budget supplies its full slow-jet input
recipe, including the finite finalized coefficient head, geometric tail,
and Borel Euler-cutoff derivatives. Its tuple retains both Pi1/Pi2 pressure
histories, moment repairs and the radial cutoff-curl commutator. A leading
background norm alone cannot replace this input.

For reference, at frozen carriers p,pz,x0 and path v,

    H=pF+pzG,
    n=(x0-v H_R,p/R,pz-epsilon v H_Z),
    n'=(-H_R,0,-epsilon H_Z),
    K=[[0,-2F,0],[2F+R F_R,0,0],[G_R,0,0]],
    Eik/epsilon=v(H_T-G H_Z)+btilde*n_r.                (1)

The old positive recipes include every displayed term. Their Eik, bR and
Bz bounds retain the original epsilon factor before any coefficient
estimate; the adapter changes no epsilon grade. The reference rate is

    p_ref=P'/P=lambda0/sqrt(1+s_sigma(v)^2)
              -lambda0(1+s_sigma(v)^2)/(1+u^2)^(3/2).  (2)

The old finite Taylor-composition formula for (2), not a constant lambda0,
supplies EVERY requested rate jet. In particular the third growth jet is
not set to zero. No new phase, carrier rounding or path is selected.

`BackgroundPhase` accepts the old typed BackgroundBudget and the current
Postslow graph. It requires the exact original closed R hull and uses
that budget's `phase_inputs`. Its lower contract still requires that the
budget actually bound the SAME physical full background, with the original
h, coefficient profiles, Borel schedule, common-lift path and value bounds.
Object identity and matching hulls do not prove any of those inequalities.
Leading lambda/shear/tilt bounds and the normal VALUE floor remain actual
native premises; upper jets cannot establish a normal lower bound.

## 2. Coordinate audit: raw R bounds are not raw common-U bounds

The old phase formulas use ordinary raw R-coordinate jets: R_1=1,
R_j=0 for j>=2 and (1/R)_j=j! Rmin^(-j-1). Q93's native phase slot,
however, is in the common U chart before its return-to-R conversion.
Simply passing those R bounds as U bounds would leave a missing transfer.
This checkpoint supplies the transfer; it does not declare that the old
standalone phase-budget formulas were wrong.

Fix a finite common ceiling N BEFORE any norm evaluation. On the same
positive hull a<=R<=b, with the original 1/3<d_radial<5/2, put

    umin=min(1,a^3), B=16(1+b)/min(1,umin),
    A_kind,N=sum_(j=0)^N (5^j/j!) raw_R_kind,j.

Q111 proves a factorial-through-R to factorial-through-U bound
486 j! B^j. A raw common-U jet of exact order j costs another j!.
For j<=N, positivity and elementary inequalities give

    j! [486 j! B^j] A_kind,j
       <= exp(486+2N^2+NB) A_kind,N.                   (3)

Indeed j!<=exp(j^2), B^j<=exp(jB) for B>0, and 486<=exp(486).
The single fixed exponent in (3) is independent of the running S/epsilon.
It is stored symbolically as an EP coefficient; no exponential or giant
factorial is evaluated to implement this forward factor. All raw jets
through N use this same conservative bound. It is not an order-uniform
estimate, an analytic-shape assumption or a smaller physical scale.

The adapter supplies (3) for normal, matrix, moving normal, reference rate,
and the linear phase coefficients. Q93/Q96 then perform their EXISTING
common-U to factorial-R and raw-path conversions at the requested order.
Thus the path is R -> U -> R -> lifted path. It is conservative, but each
leg has its own proof. No U-to-R factor is mistakenly used as its inverse.
The native five-coordinate sum is retained; Q95's six-variable through
norm, if requested later, still has its own six-variable count.

## 3. Derive the three fixed linear slots at their actual ceiling

Q94 source order j calls Q92's full linear envelope with coefficient
reserve N_j=j+2. For kind in Eik,bR,Bz, supply

    C_kind,j = 2 Tchart_(N_j) 16*188^(N_j) (N_j)!
          * [sum_(k=0)^(N_j) (5^k/k!) common_raw_kind,k]
                                                 evaluated at (25/16)S. (4)

This is exactly Q93's retained same-point return conversion. The forward
factor (3) has already produced the common raw bounds; (4) is applied ONCE
on the return. Background value and pressure derivatives are never
replaced by leading values. The remaining psi/cutoff coefficient is still
its original native R through-(j+2) bound, supplied separately. It is not
derived from phase smoothness.

For Q117 complete output n, the response ceiling is R=n+1. Its full source
table runs through R, so the largest phase reserve in (4) is n+3. The old
phase formulas need full background through one extra order, and the
background recipe needs one more Euler Borel-cutoff jet:

| This adapter's request | Order |
|---|---:|
| Current Q95 source ceiling | n+1 |
| Fixed common phase ceiling | n+3 |
| Full background slow jets | n+4 |
| Borel Euler-cutoff jets | n+5 |
| Existing source mean-factory reserve | n+3 |

At background order m its finite profile head is 1<=ell<max(2,m), with
profile jets through m. These requests are not the total primitive DAG:
root, target, regular/cut pulse, normal floor, frame and support providers
retain their separate nested obligations. No giant requested table is
materialized and no low-order archive is relabelled as high-order.

## 4. Bind the operator rate to the identical current response

`BackgroundResponse` subclasses Q117.CovariantResponse on that SAME module
family. It requires the original BackgroundPhase object on the current
wave factory. Caller-supplied phase-linear and rate tables must be None;
they are derived internally. For raw path rate order j the inherited
Q96.phase method is applied to `profile_rate`, just as for normal/matrix
jets, including the native factorial sum, state chart, S rescaling and
raw path factors. No second path conversion is applied to its result.
Insufficient ceilings, changed states/hulls and missing entrance contracts
are rejected before the large background/source expansion.

Q117 then retains its derived damping floor, covariant solve, both Gaussian
Bell-pole additions, packet evaluation, full pressure, complete cut curl
and original cutoff residual obligations. Q97 still requests both m=1 and
m=2 pairs, both cross orientations and the new square. Its covariance pole
and moment-reserve guards are unchanged. Neither local wave pressure nor
full-background norms reset current compact pressure, axial A1 velocity,
moving scalar compensation or K0 in the radial residual.

Log-Gaussian quotient/restoration jets, psi, physical cutoff jets, normal
value floor, tangent/length and regular/cut pulse inputs remain their
original providers. The rate derived here is UNWEIGHTED under the old
background contract; a poled coefficient cannot silently be substituted
for it. No actual whole-path source inequality follows merely from the
new typed binding or its metadata checks.

## 5. Verification and frontier

Exact checks compare derived coefficient tables with the retained phase
formulas and every coordinate factor. Trace controls reach the actual
background order 4 and Borel-cutoff order 5 at output zero, retaining
coefficient-one and the rest of the finite head. The varying-growth recipe
is exercised through third order. A manufactured current graph verifies
complete response, derived rate jets, both covariance carriers and the
moment guard. These are routing/algebra controls with manufactured lower
leaves, not continuum field enclosures. The analytical contributions are
the basis correction (3), finite request composition and application of
the old full-background identities (1)--(2) to this current state.
Frozen replay preserves Q117 and all historical bindings.

Next prove the actual full-background/profile/path and remaining pulse/
cutoff/source inequalities on their common whole proof domain, then rebuild
the high-pole covariance/mean/radial graph. Source order 100005, covariance
order 100016 and profile order 100046 remain unverified. No radial A0 gain,
full cone closure, unconditional all-scale pulse estimate, complete
viscous/cutoff/interaction residual improvement, infinite realization or
blow-up is established.
