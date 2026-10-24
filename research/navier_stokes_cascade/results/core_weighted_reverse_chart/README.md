# Same-weight reverse radial chart and complete initial row interface

Parent: `e30f082183bd3e2ee04c60b6f67189e524616e43`.
On 24 October 2026 its CI was in progress without failed steps; the
weighted covariance/row-gate step had passed. The parent isolated the
missing U-to-R tensor-jet inequality. Here it is derived on the SAME
positive radial hull and inserted into the existing complete initial
mean-row recipe. Historical evidence and all original roots, cutoffs,
supports, physical bands, widths and pressure gauges are preserved.

## 1. Coordinate jets on the original hull

Let 0<Rmin<=R<=Rmax be the original common-band hull and U=R^d,
1/3<d<5/2. The other four coefficient coordinates are unchanged.
The same physical cylindrical components are re-expressed as functions;
this is not a change to coordinate-basis vector or tensor components.
Put B=max(1,Rmax^2,Rmin^-1). Then

    |d R^(d-1)| <= 3B,
    |partial_R^j U|/j! <= u_j=3B Rmin^(-(j-1)), j>=1. (1)

Indeed R^(d-1)<=B separately on R<=1 and R>=1, while each additional
factor |d-1-m|<=m+2, m>=0. The product of these j-1 factors is <=j!.
Thus (1) follows from the SAME radial-power derivative bound used in
the complete-curl checkpoint, now for the reverse coordinate increment.
It is not the earlier R(U) estimate used backwards. The endpoints of
the displayed d interval also obey these safe bounds.

Define positive formal series (only the finite jet through N is used)

    a(t)=4t+sum_(j>=1) u_j t^j,
    e_0=1, e_n=(1/n)sum_(j=1)^n j a_j e_(n-j),
    Tr_N=N! sum_(n=0)^N e_n.                          (2)

The recurrence is exactly [t^n]exp(a(t)), not an analytic convergence
assumption on the physical field. These are finite positive rational
constants for rational radius enclosures; for arbitrary fixed positive
enclosures the same proof gives finite real constants. Rational outer
enclosures can be chosen once without changing a physical support.
They are degree-zero positive polynomials in S, independent of band
number, epsilon and the edge position. They are not uniform in N or
as Rmin tends to zero. The four identity coordinates account for 4t.

## 2. Pointwise raw and factorial norms, with the weight retained

For any scalar component F in the five U coordinates at the corresponding
point, suppose the RAW maximum through N satisfies

    max_(|alpha|<=N)|partial_U^alpha F|
        <= epsilon^b zeta delta^-p K_N(S).             (3)

Let M be the right side. The finite multivariate Taylor formula for
composition has, at outer order k, the majorant

    M sum_(|alpha|=k) (Delta Phi)^alpha/alpha!
        = M (sum_i Delta Phi_i)^k/k!.

Replacing the positive-degree coordinate jets by a(t) bounds the SUM
of absolute normalized mixed R jets of exact degree n by M e_n.
Consequently both the through-order factorial norm and the SUM of raw
mixed derivatives satisfy

    sum_(|alpha|<=N) |partial_R^alpha F|/alpha!
        <= M sum_(n<=N)e_n <= M Tr_N,
    sum_(|alpha|<=N) |partial_R^alpha F| <= M Tr_N.     (4)

The last step uses alpha!<=N!. This intentionally dominates the
parent's raw maximum and the historical row's submultiplicative
factorial norm, avoiding an identification of those conventions.
No additional derivative of F beyond N is requested. Formula (4)
applies componentwise and to component l1 sums by linearity.

The weights in (3) are the ACTUAL weights at the same physical point.
The parent's exact X=r^2/(2q), eta=z/q^D equalities identify them in
every band and after this coordinate composition. The differentiated
F already contains all derivatives of its flat factors. The Taylor
argument never differentiates F/zeta, bounds zeta^-1, or compares flat
exponentials at different points. Hence epsilon^b, zeta and p are
unchanged. On the flat boundary, smooth zero extension and the positive
radial diffeomorphism preserve derivative supports. All earlier enlarged
chart, common-index and closed-support separation premises remain.

## 3. Close the tensor slot, preserve the other row contracts

For the parent's U-chart covariance polynomials W_N,Dcurl_N and pole
p_N=12N+18, (4) supplies the previously REQUIRED transfer input:

    W_R,N=Tr_N W_N, D_R,N=Tr_N Dcurl_N+S_R,N.          (5)

S_R,N is the existing weighted higher-stress input in the ROW R chart
at epsilon order3 and pole at most p_N. It is added AFTER the tensor
transfer, not multiplied again. All six covariance entries, the original
real angular quadratures and selected principal-flux cancellation remain
exactly as in the parent. The tensor error order stays 2-2kappa.

`from_original` obtains Rmin,Rmax from the same returned curl operator;
it constructs Tr rather than accepting a new independent norm callback.
`ChartedPrimary.initial_rows(n,data_R,higher_stress_R)` inserts (5) at
N=n+3 into the historical row routine. Other fields and operators are
STILL inputs satisfying that routine's R-chart contracts. A caller with
only U-chart raw scalar-field bounds may apply (4) to them too, but an
operator norm must be derived for its input/output norms, and a compact
response must retain its actual potential hull. A weight identity is
not by itself a proof of either assertion.

In particular b,V,G,Ft,Fr,U,Zg,X,PiP,Rho,Kp,r retain their historical
orders and poles. Choose a common finite p covering ALL those fields,
not just the covariance pole p_N. The weighted fast map needs its higher
source jets and actual plateau gap; the slow response needs the whole
reserved hull weight floor. The pressure operator Kp is evaluated at
the doubled product pole 2p. The resulting existing row estimate is
epsilon^(2-2kappa) zeta delta^(-2p), conditional on these inputs.
The axial velocity gd-A1 gd, all covariance correlations, Pi1/Pi2 and
the WHOLE moving pressure moment PiP remain in that contract. No new
mean or residual rate improvement is claimed by coordinate conversion.

## 4. Evidence and next obligation

Independent finite Taylor expansions check (2). Exact derivatives of
fractional powers check (1); mixed five-variable derivatives including
a flat exponential check both norms in (4). Failure controls detect a
unit transfer, missing nonlinear chain terms, omitted identity variables
and missing raw/factorial conversion. Exact substitutions reproduce the
full historical rows and detect deletion of the pressure moment. Five
historical evidence records and all declared sources are rebound.
These manufactured controls support the derivation, not continuum
suprema or evaluation of the selected gigantic field.

Next supply the linked weighted fast/slow and full-pressure/profile
operator inputs on this SAME R hull, with the n+4/n+5 source jets,
actual cutoff plateaus, reserved-mean potential supports and common
pole ledger. Then compose the existing signed-source and subsequent
complete-covariance mean updates. Full background and phase validity,
value/cone acceptance, global moment identities, all physical-scale
jets, viscosity, cutoff tails and nonlinear interactions in the COMPLETE
residual remain their original obligations. Unconditional pulse bounds,
complete residual improvement, infinite realization and blow-up remain
unproved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
