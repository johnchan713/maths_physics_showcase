# Recover the calibrated normal floor from the original length proof

Parent: `0ddb0ef928e3e51b76173f847a18b150e3687cba`.
Status: `reference-scaled-normal-floor-on-original-frame-tail`.

This checkpoint corrects the recent progress frontier. The native first-jet
and directional estimates remain valid, but their separate retained-power
normal-margin inequality is NOT a prerequisite for the original full-frame
accepted-tail route. The earlier core_normal_length_closure, sourced-normal
transfer and tangent-propagator closure already supply that route. Reusing
the FULL original error estimate, rather than only its coarse conclusion,
also proves the requested calibrated phase-beta floor. We derive a sharper
constant below and pass it through the unchanged sourced-inverse recipes.

The result is conditional on valid ORIGINAL background/frame primitive
bounds, physical domains and reference calibration. It does not verify
their actual numerical enclosures, an unrelated historical band, generated
sources or the full residual. No new band, norm bound, length, cutoff,
carrier, tilt, root, support, pressure datum or width C^-20 is selected.
Historical records, including their narrower conclusions, remain frozen.

## 1. The retained proof contains more than its coarse uniform conclusion

Keep the original reference magnitude bmin<=Bs<=bmax, unit vector K,
actual pulse length L, and positive constant Dt from
[core_normal_length_closure](../core_normal_length_closure/README.md),
(1)--(3). Its SAME complete tangential-normal estimate is

    |n_tan-Bs K| <= Dt/L.                                   (1)

It includes both tilt terms, nearest-NONZERO rounding, box displacement
and epsilon*L*(p F_Z+pz G_Z). The original full-background Mz controls
both axial derivatives; G_Z is not dropped. The full radial cutoff curl
and pressure/moment data remain governed by their existing definitions.
The original selector already supplies the whole-path slow gate used to
prove (1). No leading-only replacement of H_Z is made.

If L>=Da=2Dt/bmin, (1) gives

    |n_tan| >= Bs-bmin/2 >= Bs/2.                          (2)

The last step uses bmin<=Bs. Keeping only the weaker conclusion
|n|>=bmin/2 would discard this reference-relative information.
Two unrelated statements |n|>=bmin/2 and Bs^2>=beta do NOT by themselves
imply |n|>=sqrt(beta)/2. Equations (1)--(2) DO imply it. This is why a
numerically coarse bmin^2<beta is not an additional obstruction here.

The warning in the older normal-length checkpoint about inferring a
stronger inverse constant from its coarse conclusion alone remains valid.
The new result instead reuses its stronger full error estimate and the
already retained carrier calibration. It changes neither that proof nor
its frozen controls.

## 2. Use the slack already present in the original Riccati length

Write sr=3u/2, mu=lambda_min/(1+sr), and LR=8CE/mu, with the ORIGINAL
CE and primitive bounds from frame_column_foundation. The earlier
normal-length proof establishes, before any ODE denominator is used,

    CE>=2lambda_min Da,
    LR>=16(1+sr)Da.                                       (3)

This uses the physical reference consistency
lambda_min^2<=2 Fmax gmax, itself a consequence of the growth formula.
It is not inferred from arbitrary positive diagnostic inputs.

On the rule-defined accepted tail the SAME selector guarantees

    L>=2r0 ell^2>=Lreq>=max(25,LR).

Define the uniform, scale-independent number

    theta=Dt/(Lreq*bmin)
       <=1/[32(1+sr)]<=1/32.                              (4)

For every label and every point of its whole accepted path,

    |n_tan-Bs K| <=theta*bmin<=theta*Bs,
    |n|>=|n_tan|>=(1-theta)Bs>=31Bs/32.                   (5)

This uses the existing length's slack. There is no larger Lreq or later
tail selection and no substitution of new conservative Mz or Db into an
old length. A supplied length with unknown provenance still has to pass
the ORIGINAL semantic check Lreq>=max(25,8CE/mu), with those same inputs.

## 3. The original reference calibration supplies phase beta

For the actual fixed tilt u_actual, frozen label and ORIGINAL rounded
k=ceil(epsilon^-1/2), 0<epsilon<=1, the defining formula is

    Bs^2=lambda0/[epsilon*k^2*(1+u_actual^2)^(3/2)],
    1<=epsilon*k^2<=4.                                   (6)

This is the reference magnitude from actual_band_pulse_bounds §3 and
frame_column_foundation §1, not an extra fitted amplitude. Let the
original valid primitive values satisfy lambda0>=lambda_min>0 and
0<=u_actual<=u. Then

    Bs^2>=beta_phase:=lambda_min/[4(1+u^2)^2]>0,            (7)

since (1+u^2)^2>=(1+u_actual^2)^(3/2). In particular the exponent 2
is a conservative bound, not a change of the original 3/2 calibration.
The factor 4 in (6)--(7) and the complete integer-rounding rule remain.

Combining (5) and (7) gives the uniform bound

    |n|^-2 <= qref=1/[(1-theta)^2 beta_phase]
            <=1024/(961 beta_phase)<4/beta_phase.          (8)

Even the weaker normal length Da, without the extra Riccati slack,
recovers 4/beta_phase from (2). Thus the recent native sufficient margin
4[A/S+B/S^2+C/S^3]^2<=beta_phase is optional on THIS original full-frame
route. Its validity and usefulness as a separate LOW-ORDER test remain;
that route need not assume the full original Db/frame estimates used here.
Neither route certifies an arbitrary separately fixed old numerical band.

An exact manufactured-budget control makes this distinction explicit: a
later deliberately loose native envelope fails its optional margin while
the ORIGINAL frame certificate passes at the SAME length and band. That
envelope is not substituted for the original valid Mz or Db. This is a
logical comparison of sufficient estimates, not a Navier--Stokes field
counterexample or an enclosure of the actual enormous selected constants.

## 4. Keep the original sourced inverse and its full pressure

The existing original normal bound qold=4/bmin^2 remains valid. Use

    quse=min(qold,qref).                                  (9)

The original phase coefficients, CE, Cd, endpoint frame factors, tangent
propagator T(S), Gaussian ratio P(v)/P(w), source and length polynomial
are held FIXED. Only the scalar reciprocal upper bound in ForcedBudget
is replaced by (9). The tangent theorem remains the existing
[core_tangent_propagator_closure](../core_tangent_propagator_closure/README.md);
we do not rederive its constants from the smaller quse.

Every fixed-order reciprocal, lift, adapted generator, projected-source,
ambient, velocity, pressure and potential majorant is a polynomial with
nonnegative coefficients in q and S at those fixed other inputs. This
was proved in core_sourced_normal_transfer. Induction in the unchanged
positive recurrences therefore gives COEFFICIENTWISE

    B_kind,j(S;quse)<=B_kind,j(S;qold), every S>=0.           (10)

The same convolution with the original cutoff jets preserves (10).
There are no added powers of S, inverse Gaussian costs, support losses
or physical q gains. This is a monotonicity statement about upper-bound
recipes, not a change of the actual solution or an improved residual order.
It does not compare feedback sources whose own budgets change with q.

The normalized ODE pressure still contains

    i/(km) * [ (n^T Kmat-nprime^T)t+n^T f ]/|n|^2.         (11)

The unprojected n^T f and moving-normal term remain. For purely normal
forcing Pi f=0, zero-entry velocity can vanish while (11) is nonzero.
All original potential/pressure cutoff derivatives, entrance compatibility,
whole-path support requirements, global Pi1/Pi2 memory and fresh axial
and pressure means remain separate from this denominator estimate.

## 5. Corrected frontier and reproduction

Removed extra prerequisite for this original full-frame accepted-tail
route: a separate native/Ctan/Cphi comparison against the same Lreq,
including a comparison bmin^2>=beta_phase. Added analytical assumptions:
NONE; (1), (3), (6) are retained original estimates/definitions. Actual
validity of the original primitive bounds, complete path jets and physical
domain is retained, not proved by rational controls. In particular the new
leading-first-jet estimates alone do not supply Db or every phase jet.

Next: establish the current generated full-path source and weighted
target/freezing/root/cone inequalities, then close complete covariance and
fresh axial/full-pressure means and every viscosity/cutoff/interaction
residual term. Do not reopen the already derived original normal and
tangent-propagator implications as independent inputs. Unconditional
pulse bounds, residual decay, an infinite field and blow-up remain unproved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
Use a fresh output. Exact algebra controls distinguish the full-error
argument from the invalid coarse-floor shortcut, check original rounding,
retained length slack and coefficientwise inverse/cutoff monotonicity.
Manufactured data are not samples or enclosures of the selected field.
