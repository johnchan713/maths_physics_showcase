# Q108: full current angular row from the incoming cancellation and sourced Haar

Parent: `4bd4b8fc98fdd76ec31eecccd639f6c44381d412` (Q107).
On 2 November 2026 its research CI 41514853426 and field CI 41514853252
were still running; Q106 field CI 41504983092 had passed. The inspected
remote head was Q107. This checkpoint composes the ALREADY established
Q87 angular cancellation with Q107 and Q105 on one identical current
state. It supplies a conditional coefficient for the FULL angular row,
including its Haar part. It neither supplies the old axial Haar row nor
closes the stronger mean-row target. Frozen evidence, fields, root,
accepted cover and width C_setup^-20 are unchanged.

## 1. Recover the actual incoming angular coefficient without a new callback

Q107 writes Hold for the completed first-real incoming row before the
sourced wave. For its angular component Q87 already gives a coefficient:
there is no need to replace Hold by an absolute velocity norm or to
postulate a new fine-Haar norm. Use exactly the Q86.FinalMeans instance
held by Q105's current fresh-means construction. Recheck its corrected
Q47 postslow type, same pressure state and original root identity through
Q86.FinalMeans. All four old mean stages and their full correlations are
thereby retained.

Q87's exact identity is

    Eold_theta = Haar(E_signed_theta)
                  +Rold_fast_theta+Rold_slow_theta.              (1)

The first term contains the actual initial Haar row from
state.target.source.haar.rows(n), including its assigned-stress
cancellation, PLUS the complete signed covariance change from Q85 on
that same state. The latter is explicitly requested with target_rate=2;
this remains Q51's original physical target premise, NOT a rate upgrade
inferred from a dictionary or a passing test. There is no angular mean
pressure gradient. Full pressure and the axial gd-A1 gd velocity remain
in the fields used to evaluate both differences in (1).

At N=n+3, retrieve old_full(N), fast_unweighted(N), slow(N,False), and
fixed multipliers from that same mean instance. Call the frozen Q87
angular_difference polynomial for (O,F) and (O+F,S). Retrieve the initial
Haar coefficient at pole 2h.pole and the signed covariance at its OWN
returned pole. Unweight EACH using its original edge factor before adding:

    Hi = Eedge(2h.pole) InitialHaar_n,
    Hs = [(n+1)(graph+1)+2r] Eedge(p_signed) signed_(n+1),
    Aold_n = 243 [Hi+Hs+Rold_fast_n+Rold_slow_n].                 (2)

This is Q87's existing positive coefficient program, evaluated on the
current Q105-held instance. The factor 243 is the original safe raw-max/
factorial-through padding, not a newly claimed cancellation. `old` uses
only that program and original lower constructors. Q87's separate loaded
module is used only for its pure angular_difference polynomial; class
and state checks use the identical native Q86 module inside Q105's own
factory graph. No copied/recreated physical state is substituted.

On the inherited Q51/Q85/Q86/Q87 all-fixed-order input and support
contracts, Aold has order C0=2-2kappa. This is a conditional normalized
coefficient claim; actual continuum provider inequalities and physical
uniformity have not been evaluated here. Haar contraction gives
||Hold_theta||<=Aold, with no extra Qcirc factor. Bounding the full old
angular row is sufficient and conservatively retains any old nonzero modes.

## 2. Add the current sourced Haar increment exactly once

At Q97's post_first_sourced_wave_fixed_means stage the old means, stress,
base and moving compensator are fixed. Q107 derives

    delta Htheta=(partial_R+2/R) Haar delta Wrt
                              +epsilon partial_Z Haar delta Wzt. (3)

Both cross orientations and the sourced square remain in bundle.change.
There is no graph epsilon loss after Haar averaging. Invoke the SAME
current Q98 pressure bound through n+1, retaining its state/stage/pole and
covariance n+2 guards, then Q107's pure envelope. The pressure part of
that call is not needed in (3); it preserves the original complete
changed-pressure contract for this state. No new pressure or source
callback is supplied. Old pressure and stresses remain inside (1), and
no wave change is counted in (1) a second time.

The incoming pre-new-fast angular Haar coefficient is therefore

    Apre_n = Aold_n + Eedge(p_cov) delta_Htheta_n.                (4)

Every weight is removed at its OWN exact pole on the original whole hull.
All returned coefficients here are UNWEIGHTED factorial-through-R norms.
This avoids inferring monotonicity of operator norms or promoting an
unweighted result to an unrelated weighted chart. Closed physical cover,
whole potential hull, seam/mixed-trace identities and native transfers
remain their original conditional obligations.

The incoming Haar order from (4) is C=3/2-kappa, since Q107's sourced
increment is at C and the old contribution is at C0>C. This is the current
pre-fast stage, not the original primary prefast row or final postslow
scalar. The Q103 scalar reset cannot replace (4).

## 3. Full final angular row, including the Haar part

Q105's exact cancellation with the NEW fast and slow fields gives

    Efinal_theta = Hold_theta+delta Htheta+RFtheta+RStheta.       (5)

RF and RS are the FULL differences, not their Qcirc projections; their
Haar parts must also remain. They contain all old/new mean products,
axial velocities, time, viscosity, cutoffs and the full pressure
construction upstream. `final` calls the original Q105 bound and takes
its two actual angular difference coefficients, unweighted at their
returned common pole p_final. The result is

    Afinal_n = Apre_n + Eedge(p_final)(RFtheta_n+RStheta_n).       (6)

There is no factor two on the full differences in (5)--(6). Using Q105's
Qtheta coefficient alone would omit (4); adding Qtheta to (6) would
unnecessarily double-count the differences. Metadata guards require the
same source state and sufficient returned output order. The original
projected construction is not mutated by this adapter.

For 0<kappa<1/4 the four available grades are

    old: 2-2kappa; sourced Haar: 3/2-kappa;
    RFtheta: 5/2-4kappa; RStheta: 5/2-3kappa.

Their minimum is C=3/2-kappa throughout that range. The positive margins
over C for the other three terms are 1/2-kappa, 1-3kappa and 1-2kappa.
Thus this is a FULL angular-row coefficient at C, not just a projected
angular estimate. It still falls short of the old stronger conditional
mean-row target C0 by 1/2-kappa. This is an upper-bound gap, not a lower
bound proving the canonical field cannot meet C0. Extra cancellation in
(3), or a stronger sourced covariance/response estimate, is still needed.
No claim is made that the full residual has improved to C.

## 4. Requests, verification and remaining scope

The old program requests mean fields through n+3 and the old signed
covariance through n+1, with Q85's original target/profile request graph.
The sourced program requests Q98 pressure through n+1 and covariance
through n+2. The final join retains Q105's larger row/covariance reserves
n+10,n+12, its pressure output n+1 and all earlier high-jet/operator/shape
requirements. Calling reset.check(n+1,bundle) conservatively enforces that
whole final graph even when only `bound` is requested. No giant actual
jet is materialized or declared certified.

Tests independently check the full angular cancellation, absence of an
extra Q factor, exact own-pole unweighting arithmetic and the grade
minimum across the retained kappa range. A typed current factory is built
through Q97--Q105. Manufactured lower coefficient leaves exercise the
old same-state requests and the REAL Q98 changed-pressure guards; they
are not actual canonical norm enclosures. Frozen replay recursively
checks Q107 and all unchanged source/evidence bindings. Passing these
finite checks is not a continuum proof.

Next derive and bind Hold_z with its complete averaged shifted pressure,
all assigned-stress and mean correlations, and its native coefficient
contracts; Q87's pressure-free angular argument must NOT be copied to
that slot. Then combine it with Q107's axial increment and Q105's current
A1 tail/full differences. Separately close the sourced target gap in (3)
and the radial A0/scalar gates before joining all cone components to the
inner/outer estimates. Actual pulse size/derivative/support estimates on
every physical scale, reserved-mean support contracts, complete viscous/
cutoff/interaction residual improvement, infinite realization and blow-up
remain unproved.
