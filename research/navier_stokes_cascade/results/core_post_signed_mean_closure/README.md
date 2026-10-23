# Fresh means after the generic signed wave: full regular mean-row closure

Parent: `813240e752d550d4ed0cbcc81c62f34e375d3836`.
The parent joins the generic signed covariance, full pressure changes
and current projected targets. At that stage Haar tangential rows have
order C0=2-2kappa, but the complementary rows can have order D=3/2-2kappa
and scalar/radial Haar defects can have order C=3/2-kappa. This checkpoint
derives the prescribed NEXT fast/slow corrections and joins all three
full REGULAR mean rows at C0, conditional on the original input contracts.

The broad conditional finite-cycle route was already recorded in
`current_target_cycle_bridge` and `finite_induction_foundation`. The new
result is the quantitative post-signed mean block: exact current row
differences, all three fresh scalar recipes, actual A1 velocities, A0/A1
raw source requests, and one common operator domain covering the ENTIRE
new finite request graph. The low-output signed pole alone is insufficient.

Historical evidence, actual fields, original positive roots, width
C_setup^-20, cutoff profiles, carriers and accepted band are unchanged.
These are conditional analytical estimates with finite exact diagnostics.
Original weighted inputs, enormous jets, supports and closed-hull
background contracts are not independently validated. No complete
nonzero-angular/physical residual, absolute decay, infinite field or
blow-up follows from this mean block.

## 1. Derive the post-signed fast source from the SAME current state

Put B=1-kappa, C=B+1/2, D=C-kappa and C0=2B, with 0<kappa<1/4.
Use Qcirc=I-Haar for the auxiliary torus complement; this is distinct
from the radial moment projectors Qe in the parent's target construction.
Means are still the previous current means, at orders B (tangential)
and B+1 (radial). The complete wave now includes the new generic signed
velocity, both quadratures, all old harmonics and the complete cut curl.

The current tangential rows are the old FINAL rows plus the row response
to the parent's complete covariance change. The former have order C.
The latter include all six entries and the SAME reconstructed pressure;
their graph derivative loses kappa, giving D. Thus

    F=(Qcirc Etheta_postSigned,Qcirc Ez_postSigned) in M_D. (1)

This is a consequence of those row differences, not a new fine-row input
callback. `post_signed_rows` records the old `FinalMeanRowClosure` DAG
and the generic `covariance_budget -> covariance_rows` program. It does
not replace an old raw B norm by D, transplant another mean stage, or
forward a compact/Gaussian additive residual tail through an inverse.

Apply the SAME frozen temporal inverse and compact T1:

    du=-ci^-1 Nt^-1 Ftheta, gd=-ci^-1 Nt^-1 Fz,
    psi=T1 gd, dbeta=-epsilon partial_Z psi,
    dgamma=D1 psi=gd-A1 gd.                            (2)

Actual A1 velocity remains in every later product. The original exact
time identity is

    ci Nt(du,dgamma)=-(Ftheta,Fz)+(0,A1 Fz).           (3)

It has a positive A1 residual sign and no additional temporal inverse
cost on A1 Fz. All three increments have zero auxiliary Haar, divergence
is exact, and both weighted Haar velocity moments are preserved.

For a field output n the actual fast norm functional is

    u_n <= J_n Gamma_n I_n (n+4)! Ftheta_(n+4),
    gamma_n <= J_n Gamma_n I_n (n+4)! (1+KA1_n)Fz_(n+4),
    beta_n <= J_n Gamma_(n+1) I_n (n+5)! KT1_(n+1)Fz_(n+5),
    J_n=binom(n+5,5), Gamma_n=24*2^(n+4).             (4)

Here I_n bounds the retained temporal reciprocal, and KT1/KA1 are the
ORIGINAL weighted raw-profile recipes at the new source pole. The
factorials and five-variable raw/normalized conversion are explicit.
The field orders in (4) are D,D,D+1. No extra graph loss is charged to
dgamma, because (2) uses the exact gd-A1 gd identity.

## 2. Recompute all three scalars AFTER this actual fast increment

The parent showed that a generic signed update can change P,Jtheta,Jz
at order C, even though its projected targets are C0. Keep those changes.
The complete postfast scalar polynomial is

    S_postfast = S_previous_postslow
                  + wave_change(complete signed covariance)
                  + fast_change(FULL old means, actual new fast). (5)

`postfast_scalar_envelope` implements this sum in all three channels.
Old/fast correlations generally survive Haar averaging; the old means
cannot be replaced by their Haar factors in the fast term. The orders
of the three groups are C+1, C, and min(B+D,2D), respectively, so (5)
has order C. Pressure P and axial velocity covariance in Jz are not reset.

Feed these FRESH scalars at n and n+1 into the SAME fixed five-map
response `response`, giving slow tangential fields at C and radial field
at C+1. The response uses current P,Jz for angular velocity and current
Jtheta for axial velocity/potential. It is supported on the whole original
reserved POTENTIAL HULL, including gaps between the individual bumps.
The retained base contract G=0 is required on that WHOLE hull. Nonconstant
G there can leave order-C linear scalar leaks despite the two velocity
moments; the exact negative control detects them. No actual G/support
enclosure is supplied by this checkpoint.

The old Haar field before this slow step is exactly

    previous old_Haar + previous latest slow increment.   (6)

The signed wave changes no means, and the new complete fast increment
has zero Haar. The parent already derived both terms of (6) from the
initial seed and first signed/fresh scalar chain. Both retained slow
orders exceed one. Thus (6) still has tangential Haar order 1 and radial
Haar order 2. For the auxiliary-independent NEW slow increment the exact
`haar_slow_reset` functional gives

    (Pnew,Jtheta_new,Jz_new) in S_(C+1),
    epsilon partial_Z(Jz_new+c_rho Pnew) in S_(C+2).    (7)

The derivative of c_rho is retained. Using stale presigned/prefast
scalars in the five-map would not give this cancellation.

## 3. Full pointwise tangential differences, then exact fast cancellation

Let Rf be the complete change in the tangential rows EXCEPT the principal
ci Nt of the fast tangential increment. Let Rs be the change from the
new slow increment, whose ci Nt is exactly zero. Both are the existing
`fast.difference` functional applied to the actual three-component fields,
at through-(n+3) input order. It contains the full cylindrical quadratic
flux differences, slow time, vector viscosity, and

    dgr = -t* dbeta-D1(2b dbeta+2beta dbeta+dbeta^2)
            -Dz(b dgamma+G dbeta+beta dgamma+dbeta gamma
                          +dbeta dgamma)
            +(2V du+2v du+du^2)/R
            +epsilon(Delta0-R^-2)dbeta,
    dP=int Haar(dgr), dp=T0(dgr-rho dP).              (8)

Every later derivative acts on the whole dp, including the moving rho
and current dP. Actual A1 velocities stay inside (8) for BOTH updates;
the second difference uses the old means PLUS the new fast increment.
No independent fine-pressure or fine-row callback is assumed.

For old full means B/B+1 and fast D/D+1, the sufficient remainder order
is D+1-2kappa; for slow C/C+1 it is C+1-2kappa. These conservative bounds
retain the worst full graph-viscous loss. Both exceed C0, with margins
1/2-2kappa and 1/2-kappa. The pressure fields themselves are not assigned
these tangential remainder orders.

Combining (3) with the exact row differences gives the FULL identity

    Efinal = Haar(EpostSigned)+(0,A1 Fz)+Rf+Rs.        (9)

The full Rf and Rs appear, not merely their Qcirc projections. Their Haar
parts are needed as well. `final_row_envelope` implements (9) without
adding an unnecessary Qcirc factor to these full differences. The parent
supplies Haar(EpostSigned) at C0. Thus a strict C0 bound for the exact A1
remainder completes both full regular tangential rows.

## 4. Original raw compact estimates, with strict epsilon margins

The original nonzero-mode integration-by-parts recipe gives A_e applied
to a source at order a a bound at a+p*kappa, using source derivatives
through J=n+2p+3. The raw N_n(k_e chi_U), N_J(f_e), J!, original width,
whole-line smooth extension, output weight floor, and frozen graph
M^-p temporal floor are all retained. An arbitrary unweighted inverse
norm cannot replace those factors.

For A1 Fz choose the least positive integer p1 with D+p1*kappa>C0.
For the radial source choose the least positive p0 with B+p0*kappa>C0.
At kappa=1/100000 these are

    p1=50001, A1 raw source J1=n+100005,
    p0=100000, A0 raw source J0=n+200003.               (10)

Both final powers are 2-kappa, strictly above C0. These counts are
sufficient for this specific raw-kernel rule, not universal lower bounds
on every possible proof. They are NOT evaluations or certifications of
those actual high derivatives.

Rebuild the final mean pressure from the FULL current radial source gr,
including all updated mean fields and complete unchanged post-signed
wave covariance. Its coarse pointwise source order remains B. The exact
unit-mass compensation and pressure identity give

    Er_final=-rho Pnew-A0(gr-rho Pnew),
    int Haar(gr-rho Pnew)=0, A0 source=A0 Qcirc gr.     (11)

The compensator has order C+1 by (7), and the nonzero-mode A0 estimate
in (10) is stronger than C0. Hence the full regular radial row is also
C0. This does not promote the full pressure field or discard its old
wave/base/Pi1/Pi2 memory. It controls a regular mean row, not an entire
physical residual.

## 5. One finite pole for the ENTIRE new request graph

The parent's signed pressure call needed source domain 2q. That is
insufficient for the new postfast quadratic sources. With q also at
least 6p to cover the OLD full row domain 12p, a valid new ledger is

| Quantity | Pole |
|---|---:|
| Complete signed covariance / prefast rows | 2q |
| Actual new fast fields | 2q |
| Postfast full sources and final gr | 4q |
| New slow fields on the reserved closed hull | 0 |
| KT1, KA1 proof domain | 2q |
| Kp, Kmom proof domains | 2q and 4q |

Use the original `shift_recipe(q/3,L,width,gap)` for these domains.
Its constants give a common positive-polynomial argument S+K_mean on
the SAME S>=1 domain, without a new band gate or changed physical field.
Embed lower-pole inputs before products, for example by the sufficient
max(1,L)^ceil(4q) factor. Original raw profile and cutoff geometry remain.
The ledger uses the retained unweighted fixed base/graph/chart coefficient
majorants. If another fixed multiplier has an edge pole, add that pole
at its actual products before choosing the common domains.

Crucially q must cover ALL requested field jets, not just final output n.
`dependency_trace` follows the new field/scalar/fast calls symbolically.
The slow radial field reads postfast scalars one order higher; their fast
change then reads actual fast fields and the prefast row up to seven
further orders. At the A0 source J0 this yields

    highest new prefast row: J0+10,
    highest new covariance coefficient: J0+13,
    highest NEXT signed current-target input: J0+14.   (12)

At n=0 and kappa=1/100000 the last index is 200017. Older frozen DAGs
are explicit references, NOT expanded or truncated; their nested original
A1/initial-root requests can be still larger. Thus 200017 is a NEW-LAYER
request count, not a total original-profile derivative bound.

Choose the original finite p for every nested old input before using the
old row pole 12p. With that fixed common pole, the original TargetBudget
quotient allowance can be padded to d=12p through all these finite jets.
For N equal to the largest new covariance input, a sufficient common pole is

    q=max(q_previous_signed,6p,6N+9,6N+12p+6).          (13)

For illustration ONLY, p=1 gives N=200016 and q=1200114 at n=0.
The old low-output q=30 would not cover this new request graph.
`common_pole` computes the integer indices in (12)--(13);
the giant kernel/edge constants are retained as DEFERRED arithmetic
instructions. No million-order factorial or selected field is allocated.
The actual validity of the old input bounds at the chosen p remains an
explicit missing inequality. Changing an estimate wrapper to that finite
p changes constants, not roots, pulses or the accepted domain.

## 6. Scope and next obligation

The positive functionals, exact full row differences, strict rate margins,
raw compact requests, support/moment preservation and conditional operator
domain expansions join the mean block at C0. Exact diagnostics use finite
Laurent/compact polynomial fields. The pointwise difference test allows
arbitrary COMPLETE pressure outputs, so it checks an algebraic identity;
separate original compact-operator controls test the inverse identity.
These manufactured tests are not a canonical continuum solve.

No actual high-jet/support/background enclosure or complete next-cycle
nonzero residual is asserted. The new mean feedback alone has orders D
and C when inserted against the complete wave. `mean_feedback` retains
the existing harmonic-aware coefficient 2[d*graph+d+R^-1]+2H*normal,
with H=2 for the unchanged input harmonics. Phase advection costs a true
epsilon^-1/2; the complete wave contributes epsilon^1/2, leaving exactly
the increment order. This bound uses the derived new fast/slow fields,
not candidate improved orders assigned to coarse increments. The full new signed
linear/pressure/multiharmonic residual and all additive physical tails
must still be joined on these same expanded domains. Next finish that
join, validate the actual lower inputs, and only then advance the full
finite canonical seed/induction. Sufficient physical regular orders,
controlled infinite diagonal and forcing regularity remain open before
any unconditional pulse or blow-up claim.
