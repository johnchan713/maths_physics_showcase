# Conditional initial final row and the sharp viscosity obstruction

Parent `3a966043f2b68ad1ab8f6b190c409df5734fc9da` recomputed the complete
postfast pressure and scalar defects, applied the current sharp slow map,
and derived the least strict request for the new `A1` remainder. This
checkpoint turns that request into the original raw compact coefficient
recipe, assembles the complete conditional final mean rows and pressure,
and identifies the single zero-margin group in the fast difference.

The result is still conditional on the inherited weighted fields, profile
seminorms, whole-support trace identities, and pressure/operator bounds.
The exact and manufactured controls do not evaluate the required
two-hundred-thousand-order continuum jets.

## 1. Assemble the complete final rows and pressure

Let `H=(Htheta,Hz)` be the complete post-signed Haar rows, `F` the new
actual fast increment, and `S` the refreshed sharp slow increment. With
`Uold` the complete original fast-plus-slow mean, define the existing full
row-difference functional at the two successive states by

    Rf = Difference(Uold,F),
    Rs = Difference(Uold+F,S).

Both differences retain every old/new cross and square, cylindrical and
graph derivative, vector viscosity term, and the exact reconstructed
pressure change `T0(dgr-rho dP)`. Its positive norm envelope is
`Kp(Gr+rho Pd)`. The exact temporal cancellation gives

    Efinal = H + (0,A1 Qcirc Ez_postsigned) + Rf + Rs.     (1)

No extra `Qcirc` factor belongs on the full `Rf` or `Rs`. The second
difference must use `Uold+F`, including the actual axial velocity
`gd-A1 gd`; using `Uold` again drops nonzero interactions.

The pressure after the slow step is likewise the sequential sum

    pm_final = pm_prefast + dpm_initial-fast + dpm_initial-slow
             + dpm_signed-wave + dpm_new-fast + dpm_sharp-slow. (2)

Equation (2) preserves all earlier local, nonlocal, moving-compensator and
axial pressure memory. The additive gauge remains a separate exact ledger.

## 2. Materialize the symbolic raw A1 recipe

The `A1` term in (1) uses the post-signed input state that generated `F`:

    A1 Qcirc Ez_postsigned = A1 Dz Qcirc Phi_postsigned,
    Phi_postsigned = -2 G_R Psi + gamma^2 + Wzz + pm.      (3)

Postfast or postslow pressure cannot be substituted in (3). For
`kappa=10^-5`, source order `2`, and target `nu=3-6kappa`, the least strict
integration count is `p=99995`. The original raw kernel estimate therefore
uses

    J_DzPhi = n+2p+3 = n+199993,
    J_Phi   = J_DzPhi+1 = n+199994,                        (4)

and reaches `2+p kappa=3-5kappa`, with margin `kappa` over `nu`.
The executable recipe records the exact factors

    2*6^n*8,  2^(2n+J+p+3), 3^p, J!, L,
    N_n(k_1 chi_U), N_J(f_1), Eedge(source pole,width),
    the output edge floor, and the frozen temporal floor.

This materializes a finite symbolic instruction, not its actual input
inequality. The same-state potential dependency still reaches current
fields through `J_Phi+6` and the conservative pressure field recipe through
`J_Phi+8`. Those field/profile/support bounds remain absent.

## 3. Locate the exact equality in the fast difference

Write

    B=1-kappa, alpha=2-3kappa,
    D=alpha-kappa=2-4kappa, nu=3-6kappa.

The complete source groups for the new fast field have orders

| Source group | Order | Margin over `D` |
|---|---:|---:|
| original fast difference | `alpha` | `kappa` |
| original slow difference | `2-2kappa` | `2kappa` |
| complete signed-covariance response | `D` | `0` |
| corrected original `A1` complement | `2` | `4kappa` |

Thus the signed-covariance response, whose pointwise row pays the retained
graph loss, is the unique source-rate bottleneck in this positive-majorant
route. The corrected original `A1` term is no longer the bottleneck.

For the full new-fast row difference the retained groups are

| Fast-difference group | Order | Margin over `nu` |
|---|---:|---:|
| slow time | `D+1` | `2kappa` |
| base radial flux | `D+1-kappa` | `kappa` |
| radial viscosity | `D+1-2kappa` | `0` |
| old-mean transport | `B+D+1-kappa` | `1` |
| changed axial pressure | `D+1` | `2kappa` |
| new-mean square | `2D+1-kappa` | `2-3kappa` |

The radial viscosity term is the unique equality. Every sharp-slow
difference group has margin at least `3kappa`.

Within these inherited estimates, let `delta_source` be a proved gain over
`D` for the complete fast source and `delta_viscosity` a proved additional
gain or exact cancellation in the radial viscosity contribution. The
limiting order becomes

    nu + delta_source + delta_viscosity.

Strict absorption therefore requires the precise inequality

    delta_source + delta_viscosity > 0.                    (5)

Neither gain is present in the current evidence. Generic positive norm
bounds at equality cannot supply (5). This does not prove that another
argument is impossible: a sharper signed-covariance estimate or an exact
structural cancellation could resolve it.

## 4. Scope and next obligation

Equations (1)--(5), the six-stage pressure sum, and the raw `A1` coefficient
recipe are established conditionally. The exact controls also reproduce the
independent frozen first-real rate table and confirm that all pressure and
viscosity components remain nonzero in manufactured data.

The next obligation is now twofold: validate the actual same-state `Phi`
high jets and whole-support profile bounds in (4), and prove one positive
gain in (5), or an exact cancellation replacing that estimate. Only then
can the conditional `nu` final row be promoted to a strictly absorbable
complete residual contribution. The radial `A0` term, global moment/support
interfaces, uniform physical-scale estimates, all remaining cutoffs and
interactions, infinite realization, and blow-up remain open.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
