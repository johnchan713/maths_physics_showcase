# First-real physical residual: the compact/regular stage interface

Parent `82e244169720ae94633358a66afbadeea2c36fa0` conditionally closes the
radial and axial compact tails of the **post-signed, pre-sourced** first-real
state at every fixed physical jet and requested power.  The live branch also
contains a sharper first-sourced-wave foundation, where the regular nonzero
residual improves from

    B=1-kappa  to  D=3/2-2kappa.

Those two statements cannot be spliced without changing the state.  The
first sourced wave changes the complete six-entry covariance, local mean
pressure, axial pressure moment and subsequent fast/slow means.  Therefore it
also changes both compact inputs.  This checkpoint records the exact stage
ledger and the resulting obstruction.  It preserves every historical record;
no field, cutoff, root, width, pressure gauge or accepted band is changed.

The conclusion is conditional on the named primitive accepted-tail and pulse
contracts.  Exact controls below distinguish algebraic identities from those
continuum premises.  They do not constitute an actual giant-jet enclosure.

## 1. One physical residual must use one state

Let `R` denote the first-real state closed by the parent.  With
`Qcirc=I-Haar`, its exact residual allocation is

    Rphysical_R = Rregular_R + Tcompact_R
                    + TGaussian_R + Tbase_exterior,
    Tcompact_R = (-A0 Qcirc gr_R, 0, +A1 Dz Qcirc Phi_R).       (1)

The parent proves, from the same-prefix all-fixed-order dependency graph,

    |Dphysical^m Tcompact_R| <= C_(m,N) q^N                  (2)

for every fixed `m,N`.  The Gaussian and assigned base/exterior tails retain
their separate conditional flatness contracts.  The regular nonzero source
at this state is only `W_B`: the complete original primary linear,
viscous, pressure, wave-wave and actual mean-wave terms are retained, and the
original primary contribution is limiting.

Now apply the first sourced correction `v`, with its full forced pressure,
cut potential and complete curl.  Call the changed state `S`.  Before any
later relabelling,

    W_S = W_R + Delta W,
    Delta W = Bil(w_R,v) + Cov(v),                            (3)

where all six cylindrical covariance entries, both cross orientations and
every matching harmonic are kept.  The radial pressure source changes by

    Delta g_p = -D1 Delta W_rr - Dz Delta W_rz
                    + Delta W_thetatheta/R,
    Delta P = int Haar(Delta g_p),
    Delta p_m = T0(Delta g_p-rho Delta P).                   (4)

After including the actual induced mean changes, write their full radial and
axial compact inputs as `Delta gr=gr_S-gr_R` and
`Delta Phi=Phi_S-Phi_R`.  Exact linearity gives

    Tcompact_S-Tcompact_R
       = (-A0 Qcirc Delta gr, 0, +A1 Dz Qcirc Delta Phi).     (5)

Neither term in (5) vanishes structurally.  In particular
`Delta Phi` retains `Delta W_zz`, `Delta p_m`, the axial-velocity changes and
their squares/crosses.  Equation (4) retains the moving `rho Delta P`
compensator and the same old background pressure memory.  Angular averaging
does not remove these mean quantities.  A manufactured compact covariance
fixture in the controls makes both `Qcirc Delta gr` and
`Qcirc Delta Phi` nonzero.

Consequently (2) is a bound for `Tcompact_R`, not `Tcompact_S`.  Reusing it
after (3) would omit changed covariance and pressure terms.  Conversely, the
improved order `D` is obtained only after the sourced cancellation and belongs
to `Rregular_S`, not `Rregular_R`.  A complete physical estimate may use the
pair `(Rregular_R,Tcompact_R)` or the pair
`(Rregular_S,Tcompact_S)`; it may not use the regular term from one row and
the compact term from the other.

The earlier `core_physical_tail_join` gives a conditional post-sourced tail
interface using explicit high-jet recipes.  The new initial acyclic closure
does not independently discharge that changed-state recipe: its dependency
graph ends before `v`, (3), (4), and the postwave mean refresh.  The next
proof obligation is therefore to carry the changed-state graph through those
operations and certify its primitive contracts on the same domain.

## 2. Exact physical inequalities on the two synchronized rows

Use the retained physical normalization

    Lambda=3/2+2h,    epsilon=Q^h,    0<h<1/100.            (6)

For normalized regular order `a`, the established fixed-jet dictionary gives

    |Dphysical^m Rregular| <= P_m(S+K) q^(h a-Lambda-2m).   (7)

At state `R`, `a=B=1-kappa`, so

    hB-Lambda-2m = -3/2-h(1+kappa)-2m < 0.                 (8)

At state `S`, the complete sourced linear remainder, both old/new wave
crosses and actual mean/new-wave interactions have order
`D=3/2-2kappa`; the square is stronger.  Thus

    hD-Lambda-2m = -3/2-h(1/2+2kappa)-2m < 0.              (9)

Both equalities retain viscosity, cutoffs, forced pressure and all
interactions in the regular operator.  They show that neither synchronized
finite row yields absolute physical decay, even at `m=N=0`.  This is a
failure of the available upper bounds, not a lower bound for the residual.

After reserving `eta>0` to absorb the polynomial in (7), a later complete
regular order `a` can imply a requested `q^N` bound only if

    h a > N+Lambda+2m+eta.                                 (10)

For the illustrative retained values `h=1/200`, `kappa=1/100000`, the
zero-jet exponents are exactly

    R state: -30100001/20000000,
    S state: -15025001/10000000.

The sourced solve improves the exponent by
`h(D-B)=h(1/2-kappa)>0`, but both remain negative.  Flat compact tails and
the sourced regular gain therefore solve different parts of the ledger and
do not establish a complete residual improvement.

## 3. What is established and what remains

This checkpoint establishes the exact state-indexed allocation (1)--(5),
the two physical exponents (8)--(9), and the necessary sufficient condition
(10).  It corrects any inference that the latest initial compact-tail DAG
automatically supplies the post-sourced compact inputs.  It does not revoke
the older conditional post-sourced recipes.

Next construct the post-sourced acyclic high-jet graph with the full complex
covariance (3), pressure/moments (4), actual axial velocity and refreshed
fast/slow means, and close its `A0/A1` requests from the primitive same-domain
contracts.  Then the later sourced families require the same synchronized
join, together with uniform physical-scale support/derivative bounds,
localization cutoffs, a controlled diagonal and forcing estimates.  No
unconditional pulse family, complete residual decay, infinite realization or
blow-up is asserted.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
