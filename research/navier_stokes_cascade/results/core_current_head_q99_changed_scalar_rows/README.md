# Q99: current postslow scalars and integrated sourced-wave rows

Parent: `e8e3bc333904269dd573b58f279e13569f1d38cc` (Q98).
On 2 November 2026 its CI runs 41426008465 and 41426008526 were still
running. The earlier Q97 field-comparison run 41418521802 had passed.
This checkpoint composes the frozen scalar calculus with Q97/Q98's
CURRENT changed covariance and shifted pressure. The incoming scalar
state is the COMPLETED first-real postslow state, not its earlier signed
covariance state. It derives conditional scalar and integrated-row bounds;
it does not yet construct fresh fast or slow means.

## 1. Exact state addition and complete axial pressure product

Hold every incoming final mean velocity fixed. Let Delta W be Q97's full
complex covariance change and delta pm be Q98's full shifted compact
pressure increment. For the original unit-mass moving rho and
c_rho=(int R^2 rho dR)/2, the frozen global identities give

    Theta_changed=Theta_POSTSLOW+delta Theta_wave,
    delta Theta_wave=(delta P,delta Jtheta,delta Jz),
    delta Jtheta=int R^2 Haar(Delta Wtheta,z) dR,
    delta Jz=int R Haar(Delta Wzz) dR
                               -int R^2 Haar(delta gr) dR/2,
    int R Haar(delta pm) dR=-int R^2 Haar(delta gr) dR/2
                                               +c_rho delta P. (1)

These are signed identities. The implementation adds their positive norm
majorants once; it does not claim equality of actual fields and majorants.
The old scalars are obtained from Q86 FinalMeans.postslow.scalars, whose
Q41 reset retains previous Haar means and the actual slow response.
Reusing SignedCovariance.pressure.updated_scalars here would return the
earlier pre-fast/pre-slow scalar stage. No old scalar is reset to zero,
and Q97's wave change must not be counted a second time in later updates.

All supports refer to the full fixed positive hull and every transported
fiber. The whole compact mixed-trace and pressure contracts from Q97/Q98
remain necessary; Q98's explicit seam corrections apply if they fail.
No inference is made from samples or local reserved bumps.

## 2. Integrated rows are stronger than the full pointwise row estimate

Write D_e=Dr+e/R and Dz=epsilon partial_Z. With unchanged mean velocities,
the COMPLETE zero-angular momentum row changes are

    delta Etheta=D2 Delta Wr,theta+Dz Delta Wz,theta,
    delta Ez=D1 Delta Wrz+Dz(Delta Wzz+delta pm).          (2)

Every time, viscous and mean-transport term has zero CHANGE because its
mean velocity is fixed, not because those terms vanished in the old row.
The sourced harmonic pressure is nonzero angular frequency and is not
a replacement for delta pm. Haar removes the graph derivative, then
integration by parts with compact boundary traces gives

    int R^2 Haar(delta Etheta) dR=epsilon partial_Z delta Jtheta,
    int R Haar(delta Ez) dR
                  =epsilon partial_Z(delta Jz+c_rho delta P). (3)

The entire c_rho P product is differentiated. In particular its c_rho,Z P
term survives. No slow derivative is commuted through the moving T0.
Old pressure, Pi1/Pi2 and all previous mean/covariance memory remain
additive in (1)--(3). These are integrated rows, not pointwise Haar rows.

For 0<kappa<1/4, put C=3/2-kappa and alpha=2-3kappa. Incoming Q41
postslow scalars have order alpha+1 and their differentiated integrated
rows order alpha+2. Q97 supplies scalar changes of order C and the whole
axial differentiated product of order C+1. Thus sufficient new orders are

    scalars: min(alpha+1,C)=C,
    integrated rows: min(alpha+2,C+1)=C+1.                (4)

These are conditional normalized orders before physical transfer; all
native/support and exact reset premises remain. A scalar bound of order
C does not certify a fresh mean update or an improved full residual.

## 3. Coefficients, requests and the precise next obstruction

For raw Q97 moment coefficients m_j define the conservative normalized
through-n conversion N_n(m)=sum_(j<=n) 5^j m_j/j!. The five-variable count
also bounds these two-slow-variable functions. At integrated output n,
use old postslow scalar norms O_k through n+1. Sufficient coefficients are

    Theta_k through n+1 <= O_k+N_(n+1)(delta Theta_k),
    integrated theta through n / epsilon
       <=(n+1)O_Jtheta+sum_(j<=n) 5^j delta Jtheta_(j+1)/j!,
    integrated axial through n / epsilon
       <=old_postslow.axial(n)+N_n(delta whole_axial).    (5)

The last Q97 array already differentiates the WHOLE pressure product.
Q41 axial(n) similarly bounds the whole incoming product. Equations (5)
drop only beneficial differences of epsilon grades when adding old and
new bounds. They request covariance through n+2, Q96 source through n+3,
Q94 mean factory through n+5, conservative target/profile n+19,n+32,
plus the inherited postslow/compact request graph. Q98 is also checked
at output n and the exact current covariance pole; (3)'s moment route
does not assert a new pointwise pressure derivative estimate at n+1.
All constants are symbolic positive exponential-polynomials, not enclosed
numerical continuum seminorms.

The next fast inverse requires the FULL auxiliary nonzero-Haar rows,
with whole shifted pressure, at its requested higher jets. In particular
one must derive same-domain coefficients for

    ||Qcirc Etheta_changed||_j, ||Qcirc Ez_changed||_j,
    Qcirc=I-Haar_aux,                                   (6)

and the exact-pole KT1/KA1/time inverse and compact-tail contracts. Equations
(3)--(5) do not bound (6). For example a compact covariance comparison
Delta Wr,theta=h(R) cos(y), with all other entries zero, has zero Haar
moments and zero delta pm but nonzero D2 Delta Wr,theta. This is a linear
tensor comparison, not an assertion of canonical covariance realizability
or a counterexample to the actual iteration. It demonstrates the missing
information in a scalar-only inference. The direct graph estimate for
changed pointwise rows has only order C-kappa before further structure;
no fine fast-field rate is assigned from scalar cancellation alone.

Next derive (6) by adding the sourced covariance/pressure row changes to
the actual incoming final rows, then apply the original fast inverse,
recompute its FULL correlations and fresh slow response. Q39's first-real
constructor cannot simply be reused with stale old means or scalars.
Actual trace/native contracts, all-scale pulse estimates, complete residual
improvement, infinite realization and blow-up remain unproved. No original
field, root, width C_setup^-20, accepted band or historical evidence changes.

Tests verify the two integrated identities including axial pressure,
detect loss of the moving product, retain the old scalar stage and expose
a zero-Haar comparison invisible to the moments. Manufactured norm leaves
test the state interface; none certify actual continuum data. Reproduce
with `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`.
