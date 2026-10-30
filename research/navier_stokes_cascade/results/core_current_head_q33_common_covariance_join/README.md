# Q33: common-chart covariance and initial selected-flux defect

Parent: `9127cae3f55522b11339923754ec00b62e1b47b1` (Q32).
This composes the canonical column/wave norm recipe with Q32's original
leading target and the weighted common-atlas covariance estimate. It
supplies the W and D inputs of the INITIAL Haar-row recipe from those
lower inputs. It does not claim that all canonical ODE/value or mean/source
premises have been independently proved, or that a later complex sourced
wave has the same covariance cancellation.

## 1. The edge weight is exactly invariant under native-band transfer

On the existing common neighborhood choose its coarsest band Q*, as in
`core_primary_aggregate_atlas`. For a native band Qi, put r=Q*/Qi.
The frozen coordinate changes are

    Ri=sqrt(r) R*, Zi=r^D Z*, Ti=r T*, si=r s*, D=1/2-h.

Consequently Xi=Ri^2/(2si)=X* and eta_i=Zi si^-D=eta_*.
The SAME ya,yb,delta and zeta therefore occur in all active native charts.
There is no comparison between different flat exponentials and no inverse
weight floor in the transfer. For U=R^d the map is Ui=r^(d/2)U*.
The already certified derivative factor 188^n covers this fixed linear
map and the common-torus map; its higher derivatives vanish. Discrete
labels remain frozen on the full neighborhood and radial integration hull.

The common scale has epsilon_i<=epsilon*, ki^-2<=epsilon*, and
Si<=(25/16)S*. These are the matched coarsest-chart bounds, not unrestricted
claims for a finer reference. Full coefficient cutoffs and their derivatives
remain in the native wave/curl estimates.

## 2. Bind the canonical target before aggregating waves

`canonical_root` first uses the original `root_join` factory to bind H and
W to the SAME canonical entrance/column recipe. It then supplies Q32's
factorial target bound for s^(-1-h)T0/zeta, converted into that factory's
positive-polynomial representation. K0,beta,R,mass,width, loghalf and all
complete-curl operator bounds are retained. This replaces norm callbacks,
not the physical root equation, entrance or pressure gauge.

Q32's geometry data must enclose the actual enlarged charts used here.
In particular its optional U chart must be the SAME fixed d and positive
Umin as the native root/compact-map chart. Finite manufactured inputs do
not establish those enclosures or canonical value floors.

For normalized root amplitude and complete-curl polynomials A_j,C_j, put

    a_n(S)=16*188^n*n! sum_(j<=n) A_j(25*S/16),
    b_n(S)=16*188^n*n! sum_(j<=n) C_j(25*S/16).          (1)

The factor 16 includes the existing physical normalization bound 8 and
vector norm padding 2. These are exactly the weighted coefficients in
`core_initial_haar_row_budget.weighted_covariance`; NO Eedge factor is
applied. Since delta<=1, their common half-weight pole can be padded to
6n+9. Closed enlarged label separation bounds the pointwise label norm;
normalized Haar measure introduces no deck-copy count. These are existing
support hypotheses, not consequences of frequency separation alone.

## 3. Full covariance versus the two selected entries

For the original REAL primary the exact angular covariance is

    W_ab=(epsilon/2) sum_labels(A_a A_b+k^-2 B_a B_b).

The sine/cosine cross vanishes; the complete curl square does not. Leibniz
and (1), with the weight still attached, give the raw through-order-n bounds

    W_n=2^(n-1)(a_n^2+b_n^2),
    Dcurl_n=2^(n-1)b_n^2, p_n=12n+18.                 (2)

Every one of the six symmetric covariance entries is bounded by
epsilon zeta delta^-p_n W_n. Only the two selected principal Haar entries
match the leading stress. The existing squared-partition assembly proves
this equality before any derivatives are bounded; neither signs nor
partitions are counted twice.

The complete selected-flux defect has two parts:

    bar W_re-Sigma_norm,e
      = bar Wcurl,re-(Sigma_norm-Sigma_norm,0)_e.       (3)

In the SAME common chart Q32 bounds the second part by
epsilon^3 zeta delta^(-3n-3) H_n, where

    H_n=4*Cnative_n*S_n.                               (4)

Here S_n is Q31's higher-stress constant, not the band variable S. This
axisymmetric background is evaluated directly with Q* in its enlarged
chart. It is ONE field, not a sum of localized wave copies. Therefore (4)
does not acquire the wave's 188^n, 25/16 rescaling or label-count factors.
Those factors remain in (1); omitting them from waves would be invalid.

For 0<kappa<1/4 and epsilon<=1, epsilon^3<=epsilon^(2-2kappa).
Also 3n+3<=p_n. Thus a sufficient INITIAL defect input is

    D_n=Dcurl_n+H_n,
    |D^I(bar W_re-Sigma_norm,e)|
      <=epsilon^(2-2kappa) zeta delta^-p_n D_n, |I|<=n. (5)

The sign in (3) matters as an identity; the norm bound (5) sums magnitudes.
There is no lower bound on the defect and no inference that stress or
pressure vanishes on an auxiliary support gap. A purely principal estimate
with D=0 would lose both terms in (3).

## 4. Exact interface and remaining work

`common_inputs` returns (2), (4), (5) in the existing Haar-row polynomial
type, along with their orders and common pole. `PrimaryProxy` copies exact
polynomial coefficients between historical helper module instances; it
does not evaluate polynomials or fit a numerical bound. Independent W,D
callbacks are removed relative to the listed lower contracts.

The existing initial row of output order j requests fields through j+3.
Use this adapter at that requested order, with its corresponding common
pole, before adding fresh means or signed targets. Complete curls request
one higher amplitude/root jet. This is not a claim to replace Q14's fully
expanded primitive ceiling by j+4: the mean/target/pressure routes request
additional derivatives. All distinct compact-map pole constants remain.

Next compose these initial W,D inputs with actual fast/slow mean and
pressure bounds, then the signed target and first complete source. The
canonical full-path ODE/value gates, remaining background/graph/curl input
identifications and uniform physical support estimates must still be
checked on Dv. Pi1/Pi2, moving pressure moments, full axial velocity and
the other four covariance entries remain in their original ledgers.
No unconditional full source, residual improvement or blow-up is proved.

Tests check invariant profile coordinates, exact real covariance including
cross cancellation and curl square, the coarsest-order inequalities,
polynomial transfer into the unchanged row recipe and request counts.
They are manufactured controls, not continuum supremum certificates.
Q32 and its transitive frozen sources are replayed unchanged. On 30 October
2026 Q32 focused CI 41093321841 had passed; broad CI 41093322021 was running.
Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` here.
