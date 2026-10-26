# Sharper initial target rate restores a strict first-real margin

Parent `5e6063a30f0dba6fe98f8c12f2f388fa8f0ec93f` excluded an internal
graph-frequency or radial-viscosity cancellation for a nonzero actual
selected covariance profile.  It retained the historical target order

    alpha=2-3kappa

and therefore left the first-real final row at the equality order
`nu=3-6kappa`.  Rechecking the frozen input chain reveals a rate downgrade:
`core_initial_haar_row_budget` already proves both actual pre-signed Haar
rows at

    beta=2-2kappa=alpha+kappa.                              (1)

`core_symbolic_initial_haar_join` carries that stronger rate into the
projected source, but then deliberately lowers it to `alpha` before the
two-edge target and original signed root.  Those subsequent maps are linear
and preserve epsilon order.  Retaining (1) supplies the missing strict
`kappa` margin.

This checkpoint corrects the conditional rate propagation.  It does not
supply the still-missing same-state high-order profile inequalities or prove
the complete residual.

## 1. The actual projected target has order beta

At the same postslow, pre-signed state, the frozen weighted Haar estimate is

    ||bar E_theta||_n+||bar E_z||_n
      <=epsilon^beta zeta delta^(-2p) P_n(S).               (2)

It retains the complete primary covariance, zero-Haar fast correlations,
sharp slow fields, actual `gd-A1 gd`, full compact pressure and moving
compensator.  Its limiting terms are the selected real-primary curl-square
errors.  The full auxiliary-dependent rows remain only order `alpha`; (2)
is specific to Haar averaging and the initial real prefix.

The two projector moment bumps have order three.  Thus for

    h_2=bar E_theta-sigma_2 epsilon partial_Z J_theta,
    h_1=bar E_z-sigma_1 epsilon partial_Z(J_z+c_rho P),

both `h_e` have order `min(beta,3)=beta`.  The whole moving product
`c_rho P` remains.  The exact zero moments, logarithmic pullback, chosen
rho-squared projector, two-edge primitive, reciprocal edge weight and
native-chart return multiply by fixed coefficient families with polynomial
`S` losses but no epsilon loss.  Therefore

    Xi_e=-T_e h_e in M_beta.                               (3)

The original positive-root signed inverse is linear in `Xi`.  No updated
positive root or amplitude sign condition is introduced.  The finite-loss
profile audit of the parent also shows that the auxiliary cutoff/operator
cannot independently change this epsilon order; seven extra auxiliary jets
suffice to transfer a hypothetical profile gain back to the target.

Equation (3) is conditional on the same weighted fixed-map, root, support
and profile premises as the frozen chain.  It corrects only the gratuitous
`beta -> alpha` downgrade.

## 2. Recompute the complete first-real covariance rates

Apply the general phase-locked real formulas with input `beta`, rather than
the old `alpha`.  The signed principal velocity and complete-curl remainder
have orders

    beta-1/2=3/2-2kappa,
    beta-kappa=2-3kappa.                                  (4)

All signed squares and curl terms remain.  The three complete covariance
error orders become

    beta+1-2kappa=3-4kappa,
    2beta-1       =3-4kappa,
    2beta-2kappa  =4-6kappa.                              (5)

Hence the real pair remainder is `3-4kappa`; the full covariance change is
only order `beta`.  Before Haar averaging, its selected tangential response
pays one radial graph loss and has order

    beta-kappa=2-3kappa=alpha.                             (6)

This is a gain of `kappa` over the previously used
`D=alpha-kappa=2-4kappa` response.  Other covariance entries and the whole
auxiliary-dependent pressure remain in the complete response; (6) does not
drop pressure memory or replace the full tensor by its selected entries.

## 3. Rejoin the fresh fast and final mean rows

With (6), the four current fast-source groups have orders

| Source group | Corrected order |
|---|---:|
| original fast difference | `alpha=2-3kappa` |
| original slow difference | `beta=2-2kappa` |
| complete signed covariance response | `alpha=2-3kappa` |
| corrected original A1 complement | `2` |

Thus the complete source is order `alpha`, a strict `kappa` improvement over
the old `D`.  The fresh-fast difference groups are then

| Group | Corrected order |
|---|---:|
| slow time | `3-3kappa` |
| base radial flux | `3-4kappa` |
| radial viscosity | `3-5kappa` |
| old-mean transport | `4-5kappa` |
| changed axial pressure | `3-3kappa` |
| new-mean square | `5-7kappa` |

Radial viscosity remains the limiting fast group, but it is now strictly
above the inherited required order `nu=3-6kappa`.  The post-signed Haar row
from (5) is `3-4kappa`.  The sharp slow difference is at least
`3-3kappa`.  The retained raw `A1` estimate was deliberately integrated to
order `3-5kappa`; it ties the improved viscosity without dropping below it.
Consequently the complete conditional first-real final row satisfies

    E_final in M_(3-5kappa),
    (3-5kappa)-nu=kappa>0.                                (7)

Equation (7) resolves the previous **rate** obstruction and supplies strict
absorption relative to `nu`.  It does not materialize the symbolic raw A1
or same-state `Phi` inputs: their roughly two-hundred-thousand-order field,
profile and support inequalities remain lower obligations.  Passing the
rate ledger cannot replace them.

## 4. Scope and next obligation

The correction is confined to the first real-primary signed prefix.  Later
complex waves do not have the same sharp Haar quadrature, and their generic
target/covariance rates are unchanged.  No full pointwise row is promoted to
the Haar rate.  Complete pressure reconstruction, signed squares, viscosity,
cutoffs, old/new interactions and the additive flat ledger remain in their
existing conditional formulas.

Next supply the actual same-state high-jet and whole-support bounds required
by the `A1` recipe and the corrected beta-order signed target, then re-evaluate
the radial `A0` and remaining complete residual interfaces.  Global moments,
uniform generated-source supports, all later interactions, infinite
realization and blow-up remain unproved.

Exact controls evaluate the corrected rational ledger, including every
source and final-row group.  Seven frozen records and their source digests
are replayed.  The controls establish algebraic consistency and provenance;
they are not continuum suprema or a proof of the missing high jets.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
