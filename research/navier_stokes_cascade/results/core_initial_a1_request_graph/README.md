# Corrected first-real A1 finite request graph and exact remaining inequalities

Parent `50cd291d22ab147bc3a89186b142496fbda5e45f` restores the sharper
initial target order `beta=2-2 kappa` and consequently gives the conditional
first-real final row order `3-5 kappa`, with strict `kappa` margin over
`nu=3-6 kappa`.  The remaining named input was described only as a roughly
two-hundred-thousand-order same-state `Phi`/profile requirement.  This
checkpoint computes its complete finite padding, corrects its relationship
to the older `400009` request, and separates support/profile regularity from
the quantitative field inequality that is still missing.

No field, root, projector, support, width, gauge, accepted band or frozen
record is changed.  The result is a dependency theorem and a precise
obstruction, not a numerical evaluation of the giant constants.

## 1. The immediate A1 count and raw source orders

At final-row output order `n`, retain the exact post-signed identity

    A1 Qcirc Ez = A1 Dz Qcirc Phi,
    Phi = -2 G_R Psi + gamma^2 + Wzz + pm.                 (1)

The same-state pressure `pm` in (1) is the complete post-signed,
pre-new-fast pressure.  Postfast or postslow pressure cannot be substituted.
The source has order two.  With `kappa=10^-5`, the least integer satisfying

    2 + p kappa > nu = 3-6 kappa                           (2)

is `p=99995`; `p=99994` gives equality.  The retained compact-kernel rule
therefore requests

    J_DzPhi = n+2p+3 = n+199993,
    J_Phi   = J_DzPhi+1 = n+199994.                       (3)

It produces order `2+p kappa=3-5 kappa`, preserving the strict `kappa`
margin.  This is the corrected A1 tail used in the conditional final row.

## 2. Pad the whole same-choice graph once

The axial-potential normal form reads current fields through `J_Phi+6` and
its conservative pressure-field recipe through `J_Phi+8`.  Put

    Nstate = J_Phi+8 = n+200002.                          (4)

The existing same-choice first-real covariance interface needs its signed
target through `Nstate+3` and the fixed moving projector/profile interface
through `Nstate+4`.  The target source is the complete initial Haar state.
Its projected moment reads `SharpInitial.output(j+1)`, so a Haar ledger
through `j` requires one additional lower scalar slot.  This off-by-one
dependency must be padded even though the historical constructor guard only
compared equal declared maxima.

For `n=0`, a sufficient graph is therefore

| Quantity | Required order |
|---|---:|
| `Dz Qcirc Phi_postsigned` | 199993 |
| `Qcirc Phi_postsigned` | 199994 |
| conservative same-state fields/pressure | 200002 |
| signed target and initial Haar source | 200005 |
| fixed rho-squared projector/profile | 200006 |
| lower sharp-initial scalar ledger | 200006 |
| initial field slots used by the Haar row | 200008 |
| safe canonical covariance ledger | 200018 |

These indices are finite and shift by exactly `n` for output order `n`.
They need not be materialized as arrays, factorials or sampled derivatives;
the positive recurrences and all-fixed-order estimates are symbolic.

The earlier `400009` signed-target index came from a graph that requested
both the radial `A0` physical residual tail and the axial `A1` mean-row tail.
The immediate corrected first-real final mean row contains the A1 tail but
not that radial A0 residual term.  Importing the A0 graph here doubles the
current request without proving anything about (1).  The radial A0 graph is
still required when the complete physical residual is assembled; it has not
been removed or improved.

For the illustrative inherited values `old_pole=1`, target allowance 12 and
retained wave pole 30, (4) gives the sufficient common pole

    q=max(30,6,6 Nstate+9,6 Nstate+18)=1200030,           (5)

with `KT1,KA1` on `2q` and `Kp,Kmom` on both `2q` and `4q`.  Equation (5)
is a symbolic domain requirement, not a claim that the lower input bounds
hold there.

## 3. What fixed smooth profiles and support already give

The explicit rho-squared projector has a band-uniform finite-jet recipe at
every fixed order from the one fixed smooth compact `rhohat`.  The retained
finite-support theorem propagates weighted flatness through the signed
product, complete cut curl/pressure, covariance, compact means, moving
pressure reconstruction and global chart representatives.  Consequently,
conditional on the weighted input contract, every finite index in the table
has smooth whole-hull zero extension and matching transported seam traces.
No list of 200,018 separate trace tests is needed.

This qualitative conclusion does not supply the epsilon-weighted field
sizes.  Smoothness of a fixed profile makes its requested seminorm finite;
it does not prove that the recomputed current fields obey a band-uniform
power of epsilon.  Likewise, passing finite controls does not evaluate a
continuum supremum.

## 4. The precise unresolved inequalities

A direct sufficient input for the A1 tail is

    N_199994(Qcirc Phi_postsigned)
        <= epsilon P_199994(S),                           (6)

equivalently, after retaining the explicit `Dz=epsilon partial_Z`,

    N_199993(Dz Qcirc Phi_postsigned)
        <= epsilon^2 Q_199993(S).                         (7)

Here `P_j,Q_j` must be finite positive polynomials on the same accepted
bands and the same weighted-to-raw support conversion.  Equations (6)-(7)
must come from the complete current fields and pressure in (1); they cannot
be supplied by the desired A1 output itself.

A sufficient lower-level formulation is that every actual initial
root/covariance, fast/slow mean, pressure and signed-target input through the
safe order 200018 obey its recorded epsilon, `zeta` and edge-pole inequality,
with the identical positive root, rho-squared projector, whole compact hull,
endpoint traces and pressure memory.  Existing symbolic norm maps show that
these inequalities imply (6).  They do not independently prove the earliest
canonical/background/weighted input bounds at those orders.

Thus the support/profile part of the immediate obligation is closed
conditionally for every required finite order.  The band-uniform weighted
field inequality (6), or the displayed sufficient lower-level family, is the
remaining blocker.  The radial A0 tail, global moment/support interfaces,
uniform physical pulse estimates, all remaining cutoffs and interactions,
the complete residual, infinite realization and blow-up remain open.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
