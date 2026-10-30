# Q26: the current later-tail source domain is closed under the retained maps

Parent: `8f5b46e6d0f2df1eb27ac2ca84c1cec8400b2448` (Q25).
This checkpoint identifies the domain of the LOCAL source chain. Q2 cites
`core_primary_tail_recheck`'s certified observation/definition tail; Q6 states
its source/residual claim on the accepted later-tail chart. Neither statement
identifies that domain with an arbitrary historical first band or Q8's
localized field. Q11's root-parameter correction also remains in force.

The result here is domain closure, not a source-norm theorem: on the existing
certified tail, each label entering the local source uses Q25's direct-cone
bounds, and the actual nonlocal maps do not import earlier bands. This is an
analytical consequence of the fixed-q geometry and the exact map formulas.

## 1. Bind the exact domain, labels and whole background chart

Retain the SAME exact root h, field, profiles, selector constants and integer
ell_* from `core_primary_tail_recheck` (2). Set, as in its (4),

    qv=min(qstar/64, 2^(-ell_*-2)), Dv={0<q(z,t)<=qv}.   (1)

q is independent of the physical radius and of auxiliary torus coordinates.
For any label whose enlarged dyadic support meets a point in Dv,
q/Q in [1/4,4], Q=2^-ell, hence

    Q<=4q<=2^-ell_*, ell>=ell_*.
    q'<=4Q<=16q<=qstar/4<qstar                         (2)

throughout that label's enlarged slow chart. Thus its background estimates
have their whole-chart domain, and its original canonical length gate
applies. Q25 then supplies its column/root lower constants on that chart,
relative to the retained actual provider hypotheses. This uses the
rule-defined certified selector, not an unevaluated historical numerical
floor. No new selector or small-q threshold is chosen.

The factor 16 in (2) is a bound for checking ONE background chart. It is
not a rule that moves source evaluation from q to 16q at every recursion.
An enlarged slow chart may extend outside Dv; (2) certifies its background
inputs, not the eligibility of every other label that might be active at
an arbitrary exterior slow point. Local source estimates use intersections
with Dv, or slow neighborhoods contained in its interior, with closed
boundary jets obtained by continuity. They do not assert recursive source
bounds on the entire unrestricted slow box outside (1).

## 2. Every nonlocal source operation stays on the same q fiber

The following are exact geometric facts, not support-separation shortcuts.

| Operation | Variables it integrates/translates | Variables fixed |
|---|---|---|
| Angular/auxiliary Haar projection | angle, torus | physical z,t, hence q |
| N_t inverse on zero-Haar data | torus Fourier modes | z,t,q |
| Shifted T_e and A_e | radial characteristic U and torus shift M a v_r | z,t,q |
| Radial moments, pressure primitive and compensators | radius | z,t,q |
| Fixed five-profile slow response | radial profiles and radial moments | z,t,q |
| Lifted homogeneous/forced pulse ODE | auxiliary pulse time v | slow coordinates, z,t,q |

For example the actual compact integrand is F(U+a,H+Ma v_r), not
F(U+a,H). The shift stays in the q fiber in either case, but only the
shifted formula is the construction. The full radial hull and all required
background regions must remain available. Restricting the radial interval
to a single mesh box would invalidate this argument and the compact maps.
The slow potential can fill gaps between bumps; it remains at the same z,t.

The common atlas freezes the band set on a slow neighborhood independently
of radial position. Radial integration may gather O(S^3) labels, as in Q17,
but every gathered label at that fixed z,t still satisfies (2). Deck lifts
are evaluated separately. No average of preimages or discarded pressure
memory is used to obtain this domain statement.

## 3. Closure through arbitrary finite stages

Let R_D denote restriction to complete q fibers in the INTERIOR of Dv.
For each nonlocal operator K in the table, its formula gives

    R_D K f = K_D R_D f,                               (3)

where K_D uses the same full fiber and original coefficients. Equation (3)
means two globally defined inputs that agree on Dv give identical outputs
there. It also defines the restricted operator without knowing input data
at q>qv. Products, sums and smooth coefficient multiplication are local.
Continuous derivatives commute with restriction to an open domain, since
they use its germs. The whole moving pressure product derivative must still
be taken; restriction is not permission to freeze coefficients.

Induction over any finite source-expression tree now proves domain closure.
No factor 16 is compounded with stage, derivative order or number of compact
integrations. Higher derivative requests change seminorm requirements, not
the q fiber or accepted-label floor. Smooth closed-boundary jets at q=qv
are limits from inside the same charts; they are not distributional
zero-extension derivatives.

This discharges Q25's domain-membership question for the CURRENT LOCAL
later-tail chain, relative to its explicit rule-defined domain (1). It does
not discharge Q24 margins for a different earlier domain. Q23's finite-head
upper-propagator theorem remains valid but is not needed to justify labels
inside (1), which already satisfy the direct-cone selector.

## 4. Restriction is not multiplication by a physical cutoff

No velocity, pressure or force is redefined here. If one instead multiplies
by a cutoff chi(q), derivatives produce commutators:

    partial_t(chi(q)f)=chi partial_t f+chi'(q)q_t f.      (4)

Moreover the divergence-free localized velocity is curl(chi A)+chi B,
not merely chi u. The full pressure cutoff, viscosity terms and both
nonlinear cross orientations remain. Thus (3) does not extend the local
proof to Q8's historical localization or its endpoint cylinder. Q11 already
blocks identifying the h=1/200 records with the exact-root field.

For a separately specified historical q-domain, the sufficient membership
condition is q_old<=qv, or equivalently the appropriate minimum-label and
background-chart inequalities. No such comparison is established here.
The distinct question of how a same-parameter localization is chosen remains
an obligation before any global force or endpoint conclusion.

## 5. Result and next obligation

Established: accepted-label coverage, full background-chart inclusion and
fixed-q closure of the local source grammar through every finite stage.
Inherited premises remain: same-profile background/selector validity,
full radial supports, and all analytical norm estimates for the operators.
Domain closure alone does not prove the weighted source inequality.

Next compose the first signed-prefix source recipe on Dv with Q19's typed
weights, Q20--Q21's whole-path conversion, Q23's upper bounds and Q25's
canonical root constants, using the actual complete coefficient inputs.
Retain all six covariance entries, full Pi1/Pi2 memory, axial A1 velocity,
A0 pressure errors and the differentiated c_rho P moment. No unconditional
pulse bound, complete residual improvement, force or blow-up is claimed.

The exact tests verify dyadic boundary cases, fixed-q nonlocal composition,
radial support filling and nonzero cutoff commutators. They are manufactured
controls, not evaluations of actual selected fields. Frozen evidence and
source bindings are replayed unchanged. Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.
