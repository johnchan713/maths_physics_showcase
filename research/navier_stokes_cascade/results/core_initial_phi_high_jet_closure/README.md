# Conditional same-state Phi high-jet closure from the all-fixed-order DAG

Parent `be781606f848d8d7ce3389f54288a25c6156a756` computes the complete
finite request graph for the corrected first-real `A1` tail.  At output
order zero it reaches the safe canonical covariance index 200018 and leaves
the quantitative same-state inequality

    N_199994(Qcirc Phi_postsigned) <= epsilon P_199994(S)       (1)

as the immediate blocker.  This checkpoint audits whether (1) is a genuinely
new independent hypothesis.  It is not: the retained first-real construction
is a finite acyclic composition of all-fixed-order positive norm maps.  For
every fixed requested order, those maps imply (1) from their named primitive
accepted-tail contracts.  No numerical list of 200,018 derivatives is
required.

This is a conditional implication.  It does not certify the primitive
background, phase, normal, value/cone and canonical ODE contracts themselves.
No field, root, profile, support, pressure gauge, width or accepted band is
changed, and historical frozen evidence is preserved.

## 1. The post-signed state is an acyclic finite prefix

The first real prefix has the following dependency order:

1. fixed analytic background/Borel coefficients, same-band phase and normal
   data, canonical homogeneous ODE/path cutoffs, fixed `rho`, and the original
   compact/slow inverse profiles;
2. the original positive primary root and complete curl;
3. the complete primary covariance;
4. initial fast means, then initial slow means and complete compensated
   pressure/scalars;
5. the complete initial Haar rows;
6. the rho-squared projected target, original-root signed wave and its full
   six-entry covariance;
7. the same post-signed pressure/rows, `Qcirc Psi`, and finally

       Qcirc Phi = -2 G_R Qcirc Psi
                   +Qcirc(gamma^2+Wzz)+Qcirc pm.             (2)

The desired new `A1` output does not occur in this graph.  Equation (2) uses
the post-signed, pre-new-fast pressure and retains local pressure, the moving
compensator, Pi1/Pi2 memory and every covariance component.  Consequently
there is no circular estimate in deriving its input.

## 2. Why “all fixed orders” covers 200018

Each primitive contract quantifies over every finite derivative order.  Its
constants may depend on that order and on the one fixed smooth profile.  Such
constants are degree-zero positive polynomials in the band variable `S`.
Finite Leibniz sums, chart-composition recurrences, inverse/root recurrences,
ODE variation, Haar integration, compact primitives and pressure/moment maps
preserve finiteness and the positive-polynomial class.  Denominators use the
same retained positive normal/root/value gaps.

The signed product keeps half of the flat edge weight; the complete cut curl,
covariance, compact means and pressure have smooth whole-hull extension by the
existing finite-support theorem.  Hence the induction applies at any one
finite request, including the parent's safe ceiling 200018.  Giant
factorials and profile seminorms can remain symbolic.  Their size does not
change an epsilon exponent because the derivative order is fixed before the
band tends to infinity.

This argument would fail for an order growing with the band, an infinite
stage with no controlled diagonal, or an unproved denominator floor.  None of
those is asserted here.

## 3. Propagate the rate without dropping pressure memory

The already derived same-state groups in (2) have orders

| Group | Order |
|---|---:|
| `Qcirc Psi` | 1 |
| `gamma^2` | `2-2 kappa` |
| complete `Wzz` | 1 |
| `Qcirc pm` | 1 |

Therefore `Qcirc Phi` has order one at every fixed requested jet.  The
explicit `Dz=epsilon partial_Z` gives

    N_199993(Dz Qcirc Phi_postsigned)
       <= epsilon^2 Q_199993(S).                           (3)

The compact error estimate with the least strict count `p=99995` adds
`p kappa`, so the retained final-row tail has order

    2+p kappa = 3-5 kappa,

with strict `kappa` margin over `nu=3-6 kappa`.  No new target rate or
pressure truncation is used.

## 4. What is closed and what remains

Relative to the five primitive accepted-tail contracts, (1)-(3) are now a
proved finite-order consequence rather than an additional “actual giant
jet list” premise.  Qualitative profile/support closure and the exact
order-200018 padding from the parent are retained.

The primitive contracts are still not independently certified here.  They
include the actual analytic background/Borel identification, same-band
phase/normal and value/cone gaps, canonical homogeneous ODE with full-path
cutoffs, the fixed rho/projector data, and the exact compact/slow inverse
maps on their common domains.  A failure of any one of those premises would
invalidate the induction.  Passing tests does not establish their continuum
suprema.

The next residual obligation is the radial `A0 Qcirc gr_final` tail with its
larger graph, followed by global moment/support interfaces and the complete
viscosity, cutoff and interaction ledger.  Infinite realization and blow-up
remain unproved.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
