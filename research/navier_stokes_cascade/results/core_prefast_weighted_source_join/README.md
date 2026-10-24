# Complete prefast weighted source from the charted primary

Parent: `0632e406f3244f30273bc04519ccd0b143507bd3`.
On 24 October 2026 its symbolic compact-map CI step passed; the full
run was still in progress without failed steps. The moving reserved-hull
weight comparison is ALREADY proved in `core_reserved_geometry_recheck`;
it is not a new result here. This checkpoint joins the explicit primary
covariance to the full PREFAST mean rows, supplying the first weighted
zero-Haar fast-source inputs relative to the remaining full background
stress and fixed pressure/moment map bounds. All historical evidence and
original bands, widths, supports, carriers and pressure gauges remain.

## 1. Use the full auxiliary-dependent covariance

Before the first fast correction the perturbation mean velocities are
zero. The base (b,V,G) is the SAME complete base balance, including its
pressure and additive flat residual split. Angular averaging kills every
linear primary-wave term, including time, cutoff and viscosity terms,
because the original wave has only nonzero angular modes and the base
operators are angular independent. This does not kill quadratic products
or auxiliary derivatives. The complete angular covariance W remains a
function of the auxiliary torus; it is NOT its Haar average.

For output n put N=n+3. The derived charted-primary recipe gives a
through-order norm W_N in epsilon*zeta*delta^-p, p>=12N+18, for ALL
six entries. Pick a common p also covering the full background stress
Sigma and weighted compensator rho. The fixed pressure operator Kp must
act at that p on the SAME whole integration hull. These lower inputs
remain assumptions of this join; they are not current-row bounds.

Let l=Rmin, B=max(1,Rmax^2,l^-1). The normalized through-order coefficient
bounds for the graph direction and reciprocal radius can be chosen as

    g_N=2+6B sum_(j=0)^N (j+1) l^-j,
    r_N=sum_(j=0)^N l^(-j-1).                          (1)

The complete-curl checkpoint gives normalized radial-power jets
3B (j+1) l^-j for d_r R^(d_r-1). The two components of vr have l1 norm<2;
the extra constant2 covers the ordinary radial derivative. The retained
M<=epsilon^-kappa supplies the graph cost, so Dr through output n is
bounded by epsilon^-kappa*(n+1)g_N times the through-(n+1) input.
Axial derivatives cost epsilon*(n+1). No physical scale is reselected.

## 2. Whole radial pressure before zero-Haar projection

The exact regular rows at this prefix are

    gr=-D1 Wrr-Dz Wrz+Wtheta,theta/R,
    P=int bar(gr), pm=T0(gr-rho P),
    Etheta=D2(Wrtheta-Sigma_theta)+Dz Wtheta,z,
    Ez=D1(Wrz-Sigma_z)+Dz(Wzz+pm).                    (2)

These are the historical full conservative mean equations with x=u=z=0,
not new equations. Kp includes the original shifted primitive, cutoff,
multiplier and chart constants; Kmom bounds the WHOLE radial/Haar moment.
Neither P nor rho P is discarded when taking a zero-Haar source. Their
moving slow derivatives remain in pm. Pi1/Pi2 cancel only against the
same complete base subtraction. Additive flat residuals stay in their
separate exact ledger and are not assigned a shell weight by assertion.

Set d=n+4, dr=d g_N, w=W_N and sigma=Sigma_N. Sufficient polynomials are

    Gr=(dr+r_N)w+d w+r_N w,
    PiP=Kmom Gr, Pm=Kp(Gr+Rho PiP),
    Et=(dr+2r_N)(w+sigma)+d w,
    Ez=(dr+r_N)(w+sigma)+d(w+Pm),
    Jtheta=Kmom w, Jz=Kmom(w+Gr).                     (3)

They reproduce EXACTLY the existing `full_mean` polynomials at zero
perturbation means. With mu=1-kappa, (2) gives

    |gr|_(n+1), |pm|_(n+1), |Etheta|_n, |Ez|_n
        <= epsilon^mu zeta delta^-p times (3),         (4)

using the appropriate polynomial for each field. N=n+3 safely covers
the pressure source derivative and all coefficient products. P itself
has the unweighted moment bound epsilon^mu PiP. There is no claimed
shell weight for a scalar radial integral. Rho has its own weighted
norm, so rho P has pole p after padding. Every weighted operation in
this prefast calculation is LINEAR in W or Sigma; no product-pole 2p
is needed here. Postfast quadratic mean products DO require the older
enlarged pole ledger. Formula (4) is a conservative size estimate.

The selected principal fluxes match only AFTER auxiliary Haar averaging.
Thus one cannot replace w+sigma in (3) by the fine Haar curl-square
error before the fast projection. The auxiliary fluctuations still have
order epsilon; differentiating them in Dr costs epsilon^-kappa. The
historical sharper initial Haar-row order 2-2kappa does not apply here.

## 3. Supply both raw fast-source orders without stale reuse

Since zeta is independent of the auxiliary torus, averaging is a
contraction in its weighted coefficient norm, and f^circ=f-bar(f)
has norm at most 2|f|. For fast output m the raw input coefficients are

    Et_raw=2(m+4)! Et_(m+4),
    Ez_raw=2(m+4)! Ez_(m+4),
    Ez_radial_raw=2(m+5)! Ez_(m+5).                    (5)

The factorial converts the normalized coefficient norm to the raw maximum;
it is not the five-variable output conversion J_m. These three inputs
now come from the charted primary plus (3), instead of independent current
row callbacks. Use ONE common pole p>=12(m+8)+18 for both calls, covering
all remaining lower inputs as well. The covariance order needed is m+8;
its internal original-root/background derivative requirements remain.

Substitute (5) into the parent's symbolic weighted fast-map formula with
KT at m+1 and KA at m. It gives tangential order mu, radial order mu+1,
the same pole p, one explicit S factor and the actual axial gd-A1 gd.
This discharges that formula's source slot for this FIRST prefast prefix,
conditional on the named lower hypotheses. The current implementation
returns positive scale polynomials; giant exponential operator constants
remain symbolic as in the parent, not silently cast to rationals.

## 4. Evidence and next obligation

Exact controls compare every polynomial in (3) with the historical full
mean routine, check the graph-power coefficients, full pressure memory,
factorials and distinct input orders. An explicit cosine covariance has
nonzero zero-Haar fluctuation despite zero Haar average, detecting the
invalid promotion of the fine Haar estimate. These are manufactured
finite diagnostics, not continuum suprema or actual giant field sampling.

Next recompute the COMPLETE pressure and scalar defects AFTER adding the
actual first fast velocities, retaining A1, all mean squares and all six
covariances. Feed those fresh scalars into the existing slow response on
its entire moving potential hull using the already proved weight floor.
Do not use the prefast scalars (3) for that later step. The full stress,
fixed pressure/moment operators, common-hull and phase/value validity
remain lower premises. Global moments, later sourced waves, all physical
scale jets and the complete viscosity/cutoff/interaction residual remain
open. No unconditional pulse, full residual improvement or blow-up is
proved.
