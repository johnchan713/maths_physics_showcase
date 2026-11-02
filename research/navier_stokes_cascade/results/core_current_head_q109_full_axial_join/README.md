# Q109: full current axial mean row with complete averaged pressure

Parent: `7a19b44abd921d1d43f1133ce42afd2864440b97` (Q108).
On 2 November 2026 its research CI 41523775785 and field CI 41523775888
were still running; Q107 field CI 41514853252 had passed. The inspected
remote head was Q108. This checkpoint supplies the incoming axial Haar
coefficient and composes the full current axial mean row on Q108's exact
state. These are conditional normalized estimates, not actual continuum
suprema or complete Navier--Stokes residual bounds. Historical evidence,
fields, root, accepted cover and width C_setup^-20 remain unchanged.

## 1. Take Haar of the OLD complete first-real identity

Write Einit for the initial completed mean row, delta Ws for the first
REAL signed covariance change, and R1,R2 for the subsequent first-real
fast and slow FULL mean differences. Q105/Q108's old finite-cycle
identity gives

    Hold_z = Haar(Einit_z + delta Es_z) + Haar(R1z+R2z),
    delta Es_z = D1 delta Ws_rz + Dz(delta Ws_zz+delta pm_s),      (1)

where D1=Dr+1/R and Dz=epsilon partial_Z. The old compact remainder has
zero Haar. Indeed the SAME original A1 is

    A1 f = R^-1 chi_R int_a^b r f(r,y+shift(r,R)) dr,
    Haar(A1 Qcirc f)=R^-1 chi_R int_a^b r Haar(Qcirc f) dr=0.      (2)

Translation invariance of the entire auxiliary torus, chi independent
of auxiliary variables, and the original smooth whole-hull integration
contract justify (2). This does not commute Dz through a moving A1,
discard its nonzero modes, or assert that A1 Qcirc f itself is zero.
Most importantly, the actual old axial fast VELOCITY still includes
its gd-A1 gd correction in every product used for R1,R2. Only the Haar
of the separate old residual remainder in the full row identity vanishes.
The finite Fourier control explicitly detects the change in averaged
velocity squares if that correction is deleted.

Q106's full Haar equation retains averaged full-mean products, axial
time and both viscous terms, assigned stress and the complete averaged
pressure. Equation (1) is an alternative coefficient route through its
exact cancellation, not permission to replace those products by products
of Haar means or to copy Q87's pressure-free angular polynomial to z.

## 2. Keep the complete OLD signed pressure change

At the fixed-old-mean signed-wave stage let Dij=Haar(delta Ws_ij). Then

    Hs=(Dtt-Drr)/R-Dz Drz,
    delta Ps=int_a^b Hs dr,
    delta pbar_s=-Drr+L0(Hs-rho delta Ps).                       (3)

L0 is the SAME characteristic T0 after averaging. Unit mass of the SAME
moving rho and the original covariance endpoint/mixed-trace identities
allow L0 to reduce to the left radial integral, as in Q107. They do not
remove -Drr, rho_Z delta Ps or rho partial_Z delta Ps. The initial full
pressure, including earlier base/Pi1/Pi2 subtraction and mean correlations,
is already in InitialHaar.rows(n). It is neither set to zero nor added a
second time in (3). Unmatched seams require Q98's defect terms.

To bound (3), `old` uses Q85's complete SIGNED change, not its principal
part and not total wave covariance or signed-plus-oldcurl. The old curl
already occurs in the initial Haar row. The request target_rate=2 is
Q51's original physical premise on the identical Q37 root/state; it is
not certified by requesting that number. For each j<=n+2 retrieve the
original normalized signed coefficient Sj at its returned pole pj. Embed
coefficients upward into the current Q101 pressure archive's common pole
p, via B(pj,p)=max(1,Lhi)^(p-pj), rejecting pj>p. Form raw sufficient jets

    Dj = j! B(pj,p) Sj.                                        (4)

The safe factorial padding dominates each raw derivative because the
through-j norm contains 5^j/j! times it. Actual chart/root/curl/native
contracts are still required. No operator norm monotonicity is inferred
from this coefficient embedding.

Apply the frozen Q98.envelope at output n+1 to these jets through n+2,
using the EXACT Q101 pressure(n+1,p) archive, its original KT/KA, the
whole positive hull, edge(p) and the original weighted moving rho through
n+1 (embedded from pole zero with safe 243 padding). This reconstructs
local plus primitive pressure. Then apply Q107's pure Haar row envelope
to obtain a sufficient signed Haar axial coefficient Zs_n at p. The
full pressure norm is a safe bound for its Haar; no new fine pressure
callback or unjustified weighted ordinary primitive is introduced.

## 3. Bind the complete incoming axial Haar coefficient

Use exactly the Q108.FullAngular-held Q105/Q86 state. The initial input
is the original state.target.source.haar.rows(n)['Ez'] at pole 2h.pole.
It is the complete initial Haar row program, including pressure; using
its Etheta slot here would be wrong. Take R1z and R2z from Q101.bound's
fast_difference and slow_difference on that same state and common pole p.
These are OLD first-real differences, not the later sourced mean updates.

Q101.bound also computes current sourced rows and an old compact tail
for its other outputs. Those values are not added here. In particular
its E_z total and post_signed_A1 are NOT inputs to Hold_z. This keeps the
signed wave change and the later sourced wave change separate and avoids
double-counting the latter. The conservative upstream calls/reserves and
the actual axial velocity correction remain intact.

A sufficient unweighted factorial-through-R coefficient is

    Aold_z = 243 [ Eedge(2h.pole) InitialHaar_z
                            +Eedge(p)(Zs_n+R1z_n+R2z_n) ].      (5)

Every weight is removed at its own exact pole. The extra positive 243
padding is conservative, not a physical cancellation. Under the original
native/trace/support contracts the initial and signed Haar terms have
sufficient order C0=2-2kappa; Q85's stronger target premise is retained.
The old first-real fast/slow differences have orders 3-6kappa and
3-5kappa, both greater than C0 for 0<kappa<1/4. Thus (5) gives Hold_z at
C0 without claiming a C0 bound for the entire OLD nonzero-Haar axial row.

## 4. Join the current sourced increment, tail and full fresh differences

Q107 supplies delta Hz at C=3/2-kappa from bundle.change and the SAME
current Q98 pressure bound. Add it once, unweighted at p_cov, to (5).
Q105 then gives the complete current axial mean row

    Efinal_z = Hold_z+delta Hz+A1 Qcirc Epre_z+RFz+RSz.           (6)

The CURRENT A1 tail in (6) is kept: it acts on the current total-wave
pre-new-fast row, not on the earlier first-real row. Its Haar is zero,
but (6) is now the FULL row, including nonzero auxiliary modes. Discarding
it on the basis of (2) would be an error. `final` invokes the original
Q105 bound and uses its pre_fast_A1, FULL fast_difference and FULL
slow_difference, all unweighted at their actual common returned pole.
No extra Qcirc factor is attached to those full differences.

The available grades are respectively C0, C, 2-kappa, 5/2-4kappa and
5/2-3kappa. Their minimum is C throughout 0<kappa<1/4. Along with Q108
this supplies conditional full angular and axial MEAN-row coefficients
at C, including their Haar parts. Nonzero physical angular wave residuals
are separate. This does not close the stronger target C0: the available
upper-bound gap remains 1/2-kappa. Nor does it close the radial A0/scalar
or compact modulation/repair cone gates. No canonical lower bound or
unconditional pulse/residual improvement follows.

## 5. Requests and checks

The additional old signed jet requests run through n+2 (Q85 target through
n+5 and its original deeper profile graph), with the old pressure archive
through n+1 and original initial Haar row n. Q101's old differences need
fields through n+3. The overall Q105 row/covariance reserves n+10,n+12
and all previous mixed cutoff/operator/high-jet and physical-cover
contracts remain. Native providers and large actual jets are not certified.

Exact controls prove the Haar identity for a moving characteristic
operator while retaining a nonzero A1 velocity/correlation contribution,
check the surviving local axial pressure force, and independently check
the factorial padding and grade minimum. A typed current factory with
manufactured lower leaves checks initial axial versus angular slots,
all signed jet calls, exact pressure pole, current Q98 stage/order guards,
and the distinction between the old zero-Haar remainder and the retained
current full tail. Frozen replay binds all sources to unchanged Q108
history. Tests are finite diagnostics, not a continuum proof.

Next audit the radial A0/scalar component on this same completed state,
and the special sourced cancellations needed to improve C to the required
cone/mean target, before joining with the inner/outer estimates. Native
coefficient and reserved-mean support contracts, uniform actual-pulse
size/derivative/support estimates across all physical scales, complete
viscous/cutoff/interaction residual improvement, infinite realization and
blow-up remain unproved.
