# Original prefast quotient correction and conditional whole-source fast inverse

Parent `8431337a981afb25fd228070e85c86a057519dd7` was verified remotely on
25 October 2026. Its GitHub Actions run 40115447738 completed successfully,
including the initial-cancellation checkpoint. Historical frozen evidence,
the original fields, roots, carriers, supports, width and pressure gauges are
unchanged.

The parent safely retained the original axial `A1` fast-time remainder and
bounded it by the generic oscillatory-kernel route. For the normalized first
source this produced a conditional request for raw source order `200001`.
That request is sufficient, but it is not necessary for this particular
prefast row. Applying the already proved compact-source quotient identity
before estimating gives a direct order-two bound from low finite jets. This
checkpoint records that correction, joins it to the parent's complete regular
source, and substitutes the resulting whole source into the original
temporal/compact fast-map coefficient program.

The result remains conditional on the original weighted covariance, fixed
profile, whole-support, pressure and moment inequalities. The exact identities
and positive coefficient programs do not validate those continuum inputs.

## 1. Keep the local pressure complement before applying A1

At the ORIGINAL prefast state the perturbation means vanish. Write the six
primary covariance entries as `W`, and set

    S=Wrr,  H=Wrz,  C=Wtheta,theta.

The complete radial mean source and pressure are

    gr=-D_r S-S/R-Dz H+C/R,
    P=int Haar(gr),
    pm=T0(gr-rho P).

Define

    r0=-S/R-Dz H+C/R,       K0=-S.

Then, exactly,

    gr=r0+D_r K0.

On the SAME whole compact potential hull, including the zero extensions and
endpoint traces, the compact identities give

    A0 D_r K0=0,   T0 D_r K0=K0,
    P=int Haar(r0),
    pm=K0+T0(r0-rho P).                                    (1)

This is a pressure decomposition, not deletion of pressure. In particular
the local term `-Wrr` must remain; omitting it changes the axial force. The
full pressure still contains the moving compensator `rho P`. Only after the
auxiliary projection is taken does `Qcirc(rho P)=0`, because `rho P` is
auxiliary independent and `T0` preserves auxiliary Fourier modes.

The original axial row is

    Ez=D1(Wrz-Sigma_z)+Dz(Wzz+pm).

`Sigma_z` is auxiliary independent. Under the retained whole-support
contract, `A1 D1=0` on `Qcirc Wrz`. Therefore (1) yields

    A1 Qcirc Ez
      =A1 Dz Qcirc Phi0,
    Phi0=Wzz-Wrr+T0(r0-rho P).                            (2)

Every viscosity, cutoff and pressure term represented in the original
prefast covariance/pressure construction remains accounted for. Equation
(2) concerns only the exact `A1` remainder created by the fast-time
cancellation. It does not replace the pointwise row or the pressure field.

## 2. Direct weighted bound and corrected derivative ledger

Let `W_N` bound all six original covariance entries through normalized order
`N=n+3`, let `r_N` bound `1/R`, and put `d=n+4`. A sufficient positive
program is

    R0=(d+2 r_N) W_N,
    P0=Kmom R0,
    Pfull=KT0(R0+rho P0),
    Qpressure=2 KT0 R0,
    QPhi0=4 W_N+Qpressure,
    Tail_n=d KA1_n QPhi0.                                 (3)

`Pfull` is recorded to show that the complete compensated pressure remains;
it is not substituted into the sharper projected channel. The factors two
in (3) are the norm of `Qcirc`. `KA1_n` is read from the ORIGINAL compact
profile recipe at the same supplied pole and axial support, rather than from
an arbitrary operator norm.

The covariance has order one. Formula (2) retains the explicit `epsilon` in
`Dz`, so the projected tail has order two:

    ||A1 Qcirc Ez_prefast||_n
       <=epsilon^2 Tail_n(S).                             (4)

The parent's complete regular current source has order

    D=2-4 kappa.

Thus (4) has strict margin `4 kappa`, and the whole current source has order
`D`. No nonstationary integration is needed for this normalized source. At
fast output `m`, the direct tail reads primary covariance only through
`m+7` and `m+8` for the two source calls. The original fast/slow difference
program already reaches old primary order `m+18`, so that remains the largest
old-primary request. The previous conditional raw order `m+200009` is not a
requirement of this normalized route.

This corrects the parent's derivative ledger, not its exact cancellation or
its generic high-jet estimate. The archived `200001` route remains a valid
sufficient estimate for an arbitrary zero-Haar row lacking (2), and it may
still matter in distinct physical-transfer arguments. It cannot be presented
as the necessary next input for the source treated here.

## 3. Feed the complete source into the original fast map

Let `(Ftheta_n,Fz_n)` be the parent's regular source, containing the original
fast difference, sharp slow difference and complete same-choice covariance
response, all with full pressure and viscosity. Define

    Htheta_n=Ftheta_n,
    Hz_n=Fz_n+Tail_n.                                     (5)

The projection norm is already included in both terms of (5). It must not be
charged again. For the original fixed-label temporal and compact maps, with
`I_m>=c_i^-1/S`, the conditional coefficient program at output `m` is

    Gm=binom(m+5,5) 24 2^(m+4) I_m S,
    Um=(m+4)! Gm [Htheta_(m+4)+(1+KA1_m)Hz_(m+4)],
    Xm=(m+5)! binom(m+5,5) 24 2^(m+5)
                      I_m S KT1_(m+1) Hz_(m+5).            (6)

The new fast velocity bounds are `(x,u,z)=(Xm,Um,Um)`. The axial velocity is
the actual `gd-A1 gd`; (6) never replaces it by `gd`. Conditional orders are
`D+1` for the radial field and `D` for the angular and axial fields.

The word “conditional” is essential. Equations (3), (5) and (6) materialize
finite positive-polynomial norm recipes, subject to the original actual
profile, support and field enclosures. They do not sample or prove the
selected enormous canonical fields. The short manufactured fixture used by
the tests deliberately cannot evaluate a full new fast output; it checks the
formula independently and verifies that an uncovered ledger is rejected.

## 4. What is established and what remains

Analytically established here:

* the exact original prefast quotient identity (2), including the local
  `-Wrr` pressure memory;
* the conditional direct order-two `A1` bound (3)--(4);
* the correction of the normalized `200001`-jet requirement;
* the complete conditional source (5) at order `D`; and
* the original temporal/compact fast-map coefficients (6), including the
  actual axial correction and all source channels.

Still required are actual continuum enclosures of every named weighted
primary/profile/support input and the common whole-hull trace identities.
The immediate next calculation is to recompute the complete pressure and
scalar defects after adding (6), retain its new exact `A1` time error, and
only then apply the sharp slow response. The first-real zero-margin issue,
later `A0` high-jet inputs, global pressure/moment support, uniform physical
pulse estimates, complete viscosity/cutoff/interaction residual, infinite
realization and blow-up all remain open.

The exact controls check the radial quotient, local pressure sign,
auxiliary projection, positive coefficient formulas, corrected request
orders and the fast-map substitution. Their symbolic and manufactured data
are diagnostics, not continuum suprema. Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
