# Actual-profile gate for the initial source-to-viscosity equality

Parent `c7519993b5e45ec4e6f3fe08ebfdc32280d91b9b` proved that the inherited
size, support and zero-Haar premises permit simultaneous equality in the
first-real covariance source and radial-viscosity estimates.  It correctly
left open whether the specific constructed pulse has an exact cancellation.
This checkpoint audits that possibility against the actual auxiliary
profile and irrational graph directions.

The result is conditional but sharper.  Whenever one selected signed target
component is nonzero, its actual compactly supported covariance profile has
no internal graph or radial-viscosity cancellation.  Auxiliary independence
of the target and signed **multiplier** does not make the physical pulse
profile auxiliary-independent.  This excludes two proposed exact escape
routes.  It does not give an epsilon-uniform lower bound, so the strict
residual gain remains unproved.

## 1. Keep the physical profile after applying the signed inverse

Use the same fixed positive roots, cutoffs, slow partition and real
phase-locked columns as `real_signed_covariance_update`.  For either selected
entry `e` in `{r theta,r z}`, write the principal cross covariance on one
common auxiliary torus as

    Lambda_e(y)=sum_beta eta_beta^2 Lambda_beta,e(y).       (1)

Its Haar average is the auxiliary-independent target `Xi_e`.  Each physical
`Lambda_beta,e` still contains the squared transverse and temporal pulse
cutoffs, schematically

    [chi_g(xi) psi(v/L)]^2

together with the ODE columns and the fixed signed multiplier.  The latter is
auxiliary-independent; the displayed profile is not.  The inherited
construction puts `chi_g` strictly inside the transverse rectangle, and the
slow partition does not change auxiliary support.  Therefore, at a fixed
slow point where `Xi_e != 0`, (1) is nonzero and vanishes on a nonempty
transverse open set.  It is nonconstant.

This statement uses the exact Haar identity and the actual support contract.
It does not assert `Xi_e != 0` everywhere.  If both selected targets vanish,
or have a better epsilon order, that is a separate branch which may remove
or improve the limiting source.

## 2. The irrational Fourier multiplier is injective

On the original common torus let

    L=v_r dot partial_y,     N_t=v_t dot partial_y,
    v_r=(1,1-sqrt(2)),       v_t=(sqrt(2)-1,1).             (2)

For a nonzero integer Fourier mode `k=(k1,k2)`, set `a=k1+k2`.  Its two
divisors satisfy the exact identity

    (v_r dot k)(v_t dot k)
      =(sqrt(2)-1)[a^2-2 k2^2].                            (3)

The bracket is a nonzero integer: equality would make `sqrt(2)=a/k2`,
apart from the already excluded zero mode.  Hence neither divisor vanishes.
The zero-Haar temporal inverse is defined mode by mode, and the leading
source-to-viscosity map has multiplier

    L^2 N_t^{-1}L: q_k -> (2 pi i v_r dot k)^3
                             /(2 pi i v_t dot k) q_k.      (4)

It is nonzero for every `k != 0`.  Its kernel on smooth periodic functions
therefore consists exactly of constants.  By Section 1, a nonzero selected
actual profile is not constant.  Consequently

    L Lambda_e != 0,
    L^2 N_t^{-1}L Lambda_e != 0.                           (5)

Thus neither radial-frequency orthogonality nor vanishing of the active
principal-cross graph coefficient can cancel the actual leading term when
`Xi_e != 0`.  Compact support, rather than a numerical sample of the pulse,
is what rules out a constant profile.

## 3. Exact rate ledger and the remaining inequality

Put

    alpha=2-3 kappa,  D=alpha-kappa,  nu=D+1-2 kappa=3-6 kappa.

The selected tangential row contains the leading graph source

    M L Lambda_e,                                           (6)

the fresh fast inverse contains

    -c_i^{-1} M N_t^{-1}L Lambda_e,                         (7)

and its radial vector-viscosity difference contains

    epsilon c_i^{-1} M^3 L^2N_t^{-1}L Lambda_e.              (8)

Here `N_t^{-1}` is applied to the zero-Haar function `L Lambda_e`.
Equivalently one may write `L^3N_t^{-1}Qcirc Lambda_e`; it is not an
inverse on the generally nonzero Haar mode of `Lambda_e`.  Fixed nonzero
factors from `d_r R^(d_r-1)` multiply (6)--(8), while differentiating those
factors produces only the already higher, lower-graph-power terms.

Lower graph-power connection terms retain the positive margins already
listed by the parent.  Equations (5)--(8) exclude exact annihilation inside
the principal selected profile.  They do **not** show that its normalized
coefficient stays bounded below as the physical band tends to infinity.
Nonzero functions can still be smaller than every presently recorded upper
envelope.

A sufficient quantitative route to strict closure is now explicit: for
some `delta>0` and every required fixed jet order `n`, prove on every
accepted band

    N_n(L^2 N_t^{-1}L Lambda_e)
       <= epsilon^(alpha+delta) P_n(S)                      (9)

for both selected components, with finite positive polynomials `P_n` and
the same-state roots, cutoffs and supports.  Equivalently, a corresponding
positive-order gain for every active `Xi_e` may be propagated through the
fixed signed inverse.  Current evidence supplies only order `alpha`, so (9)
is not established.  A new complete same-order row identity could also
replace (9), but the inherited positive rate ledger has no other group at
`nu`.

The actual-profile audit therefore corrects the interpretation of
"auxiliary-independent" and narrows the obstruction: the two simple exact
cancellations are unavailable on a nonzero target branch; only quantitative
extra smallness (or a new complete-row identity) can give the strict margin.
The same-state two-hundred-thousand-order `Phi`/profile bounds, radial `A0`
route, global moments and supports, remaining cutoff/interactions, complete
residual, infinite realization and blow-up remain open.

## 4. Evidence

Exact controls verify (3) on representative near-Pell modes and check
injectivity on a real finite zero-Haar Fourier fixture.  Logical controls
separate the nonzero-target branch, strict transverse support and the absent
quantitative epsilon gain.  Five frozen dependency records and their source
digests are replayed.  These checks support the algebra and provenance; they
are not continuum profile bounds or numerical proof of (9).

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
