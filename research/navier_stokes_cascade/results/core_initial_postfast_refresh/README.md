# Complete initial postfast pressure/scalar refresh and sharp slow input

Parent `87fdc99cfa64dcd36697ff73d35856607ae39efa` replaced the unnecessary
normalized `200001`-jet estimate for the original prefast `A1` remainder by
the exact compact quotient and fed the resulting complete order-`D` source
to the original fast map. This checkpoint performs the next required
bookkeeping: it recomputes the pointwise pressure and scalar defects after
that new fast field, applies the original sharp slow five-map, and derives
the precise remaining request for the new `A1` time error.

All statements are conditional coefficient estimates. The inherited actual
weighted fields, fixed profiles, whole-support identities, and operator
bounds remain premises. No finite symbolic test below is a continuum norm
certificate.

## 1. Complete postfast pressure memory

Let `pm0` be the full original prefast pressure, including the compensated
radial moment, and let `dpm_f0`, `dpm_s0`, `dpm_w`, and `dpm_f1` denote,
respectively, the original fast difference, original sharp slow difference,
complete signed-wave covariance response, and new corrected-fast difference.
The current pointwise pressure envelope is the sequential sum

    pm_postfast = pm0+dpm_f0+dpm_s0+dpm_w+dpm_f1.          (1)

For the last term, the old field is the complete original fast-plus-slow
mean and the increment is the actual new `(x,u,z)` fast field, whose axial
entry is `gd-A1 gd`. Each difference uses the full `T0(gr+rho P)` pressure
program. Thus (1) keeps local, nonlocal, axial, moving-compensator, and
signed-covariance pressure memory. The additive pressure gauge remains a
separate exact ledger and is not estimated by (1).

## 2. Refresh the scalars once

The same-choice pressure recomputation already returns

    updated = original postslow scalars + complete signed-wave change.

Consequently the postfast defects are exactly

    current = updated + fast_change(old_full,new_fast).    (2)

The full-old/new-fast quadratic and cross terms in `fast_change` are kept,
including the distinct axial velocity. Adding `wave_change` again would
double count every signed-wave scalar contribution; the controls compare
that incorrect expression directly with (2).

The original moving-hull weighted five-map now gives

    u_s = floor H (Kp P + Kz Jz),
    z_s = floor H Kg Jtheta,
    x_s = floor (n+1) H_(n+1) Kpsi_(n+1) Jtheta_(n+1).    (3)

Hence, with `kappa=10^-5`,

    alpha=2-3 kappa,  D=2-4 kappa,  nu=3-6 kappa,

the refreshed scalar and tangential/axial slow orders are `alpha`, while
the radial slow field has order `alpha+1`. The resulting scalar reset has
order `alpha+1`, conditionally on the existing exact reserved-hull solve.

## 3. The new exact A1 error is still a high-jet obstruction

The exact `A1` remainder left by the new fast cancellation must be evaluated
on the post-signed input state that generated that fast field. Its complete
pressure is retained; neither the refreshed postfast pressure (1) nor a
later postslow pressure may be substituted into this source identity:

    A1 Qcirc Ez_postsigned = A1 Dz Qcirc Phi_postsigned,
    Phi_postsigned = -2 G_R Psi + gamma^2 + Wzz + pm.      (4)

The source in (4) has normalized order `2`. Direct order two is below
`nu`. Every allowed compact integration gains `kappa`, and the inequality

    2+p kappa > 3-6 kappa                                (5)

has least integer solution `p=99995`; `p=99994` gives equality and is not
a strict bound. The potential route therefore requests

    J = n+2p+4 = n+199994,

with margin exactly `kappa`. This is seven raw orders smaller than the old
generic full-row request `n+200001`, but it is still unmaterialized. At this
`J`, the complete dependency ledger reaches current fields `J+6`, current
`Ez` and `Etheta` `J+4`, coarse pressure `J+5`, and its conservative field
recipe `J+8`. The improvement in the raw potential index must not be
misreported as an actual all-scale field enclosure.

There is a second sharp obstruction. The inherited final fast-difference
order satisfies

    D+1-2 kappa = nu,

so it has zero strict margin. Thus (1)--(5) do not establish a complete
residual improvement even if the named high jets are later supplied.

## 4. Scope and next obligation

Established analytically here are the five-term pressure decomposition,
the no-double-count scalar refresh, the sharp slow coefficient formulas,
and the least strict `A1` integration count. The executable checks use exact
positive polynomials and manufactured finite data; they verify identities,
dependency use, and guards only.

The next obligation is to materialize or sharpen the complete same-state
post-signed `Phi` high-jet/support estimates in (4), retain the exact `A1`
term, and assemble the final row including the sharp slow difference. The
zero-margin fast difference must then be resolved before any strict complete
residual improvement. Later `A0`, global pressure/moment support, uniform
physical pulse estimates, viscosity/cutoff/interactions, infinite
realization, and blow-up remain open.

Reproduce the checkpoint with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
