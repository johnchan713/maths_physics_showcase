# Same-choice recomputation after the rho-squared projector realization

Parent `8076cd797ce1c04b367b1f735764bb728b74d208`, verified on the remote
branch on 24 October 2026. Its new projector-realization CI step passed;
run 40100438812 remained running without reported failures. Historical
frozen records, old signed-target objects and unrelated work are preserved.

The parent explicitly CHOSE new projectors

    sigma_e=R^-e sqrt(s) rho^2/A2, A2=int rhohat^2, e=1,2,

and derived their uniform finite-jet bounds from the retained fixed rho
profile. It did not install that choice in a current signed target.
This checkpoint recomputes the CONDITIONAL norm chain on that choice,
rather than copying coefficients belonging to unspecified old sigma_e.
It evaluates the original signed-amplitude/root/curl recipes in the
positive symbolic coefficient ring. It does not evaluate continuum
fields or certify the original lower profile/initial-field enclosures.

## 1. Bind the new projected sources to the same pre-signed state

The complete initial Haar state is independent of this FIRST signed
projector choice. Require precisely the same primary chart and positive
primary root object as the pressure-profile geometry container. The
original pressure profile, exact mass/support and initial field/row
identifications remain mathematical premises, not consequences of type
or identity checks. In particular, initial Haar pressure and its whole
c_rho norm must concern the same retained rho profile.

Create a NEW ChosenSource from that initial Haar state and the explicit
rho-squared family. The actual source definition is now fixed as

    h_e=E_e-sigma_e M_e,
    M_2=epsilon partial_Z Jtheta,
    M_1=epsilon partial_Z(Jz+c_rho P),
    Xi_e=-T_e h_e.

These exact moments give int R^e h_e=0 by the parent's normalization.
No scalar moment is reset to zero. All old mean velocities, complete
covariance and Pi1/Pi2 pressure memory remain inside the incoming E_e
and scalar defects. The full moving product is differentiated.

For each n the recomputed logarithmic source coefficient is

    H_e,n=PB_e(n) Row_e,n + B_log,e(n) Moment_e,n,
    B_log,e(n)=exp(E) PB_e(n) NativeSigma_e(n).

PB is the previously derived same-weight Haar pullback. NativeSigma is
the parent's rho-squared derivative recipe. The weight floor and pullback
are inserted ONCE. This removes an independent arbitrary sigma-norm
callback for this explicit choice; fixed rho-profile and complete initial
Haar inputs still require their stated uniform enclosures.

## 2. Rebuild the signed root and full covariance, rather than relabeling

From the new source, construct a NEW original SignedTarget, then a NEW
RealCovariance. The primary positive root is unchanged; the target-dependent
signed amplitude is freshly computed from the new H_e,n through the
two-edge primitive, reciprocal weight, native-chart composition and
original signed-amplitude formula. The complete curl is freshly computed
from that new amplitude, including its next derivative. No old signed
amplitude/curl norm is substituted, and no fine target callback is used.

All six covariance entries retain the same exact first-real tensor
identity. Its upper bound is rebuilt using the new principal and full
curl coefficients, keeping the principal/signed cross, both curl crosses,
signed square and curl square. The same-weight reverse chart transfer
is applied once. This remains the FIRST real phase-locked update; it
does not extend the real-quadrature cancellation to later complex waves.

The algebraic rates remain alpha=2-3kappa for the complete change and
nu=3-6kappa for its real error, conditional on all original value/support
and input-class premises. New coefficients do not independently verify
the required signed-root value/cone bound or actual coefficient supports.

## 3. Native pressure/projector coefficients on the new domain

For pressure output n, the new covariance requests order n+2, and its
curl requests target/source jets through n+3. Initial projected moments
at that last order also request whole c_rho jets through n+4. The
constructor checks this finite interface before evaluation. The finite
MovingCompensator interface conservatively receives fixed-profile data
through n+4; analytically this last c_rho derivative uses the fixed H_0
moment and coordinate jets, not a new fourth derivative of rhohat.

Rebuild the moving pressure-compensator container around the NEW
covariance state with the SAME profile bounds, hull and width. At the
new covariance pole, use the previously proved native Haar pressure
constant and both T_e multiplier constants, plus the derived Kmom,
rho and c_rho bounds. The pressure callback's common projector bound is

    B_native(n)=WeightedSigma_1(n)+WeightedSigma_2(n).

It is NOT B_log: this functional acts on native R-chart projected rows.
Using the logarithmic bump there would mix coefficient charts. Each
channel separately is dominated by this common positive sum.

The full pressure_target_budget is then reevaluated on all NEW complete
covariance coefficients. Add changes to the OLD current scalar defects;
do not replace them by wave changes. Keep the complete next-order
epsilon partial_Z(Jz+c_rho P), including c_rho,Z P. The historical
source, root and pressure functions are unchanged and their evidence
is replayed, so earlier diagnostics remain reproducible.

## 4. Identity guards and scope

The constructor records the chosen profile-bound data, geometry and
original root. Public covariance/pressure operations reject changed
profile data, changed pressure-container data, swapped source/family/
target identities or substituted primary roots. These checks prevent
reuse of cached coefficients after routine object substitution. They
do not prove that supplied norm bounds enclose a physical field.

This closes the missing CONDITIONAL same-choice norm recomputation,
not the original continuum input inequalities. In particular the
initial Haar bounds, actual fixed rhohat jets, root/value/cone and
support identifications, whole reserved potential hull and all nested
high-order source data remain premises. Marked manufactured LOWER Haar
norms are used only in diagnostics; the signed amplitude/root/curl and
complete pressure recipes themselves are evaluated without mocked
signed coefficients. Historical targets are not asserted equal to these
new targets, and no actual continuum solve is claimed.

Next apply the fresh fast/slow mean construction to this SAME recomputed
state and its CURRENT defects, retaining actual A1 velocities and
the full required domain. The first-real request ledger still reaches
target order 400009 at output zero, before nested old dependencies.
This finite low-order reproduction does not supply those high-jet
enclosures or a strict gain beyond the fast-difference boundary.
Uniform actual pulses, complete viscous/cutoff/interacting residual
closure, infinite realization and blow-up remain unproved.

## Evidence

Exact controls compare the recomputed source, original SignedTarget
amplitude, all covariance summands and full pressure/scalar functionals.
They distinguish native and logarithmic projector bounds and retain
old scalar additions plus whole next-order pressure memory. Negative
controls reject insufficient ledgers and changed family/root/pressure
identities. Frozen records are hashed and replayed unchanged. These
finite algebraic checks are not a continuum proof.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
