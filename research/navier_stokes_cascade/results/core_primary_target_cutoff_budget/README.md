# Fixed-profile leading target, flat logarithm and complete cutoff budgets

Parent: `dd73978319d6f96c410964552b1fdddfa8e1437a`.
Status: `conditional-fixed-profile-primary-target-and-cutoff-budget`.
On 24 October 2026 the parent's CI run 40036876326 was in progress without
failed steps; its complete curl-operator step had passed. Historical frozen
records, original target/roots, profiles, width, geometry and accepted bands
are unchanged. This supplies the remaining target/loghalf/cutoff norm recipes
of the quantitative primary chain from fixed profile seminorms and chart
identities. No selected giant profile seminorm is numerically evaluated.

## 1. The actual target has a smooth cubic-pole numerator

Distinguish shell length Le=log(Xb/Xa) from pulse length Ls. Set
x=log(X/Xa), yb=Le-x, a=t1^2, b=4 and retain EXACTLY

    zeta=exp(-a/x^2-b/yb^2), delta=min(1,x,yb).

The retained all-order edge foundations give the ORIGINAL finalized leading
stress factorizations

    inner: T0=exp(-a/x^2) Aa(x/t1,eta),
    outer: T0=exp(-b/yb^2) yb^-3 (b0,yb^6 bz), b0>0.    (1)

The stated factors have smooth fixed jets on the closed respective collars;
the actual compact modulation/repair and full incoming pressure memory are
part of those defining-map premises. On the middle T0 and zeta are smooth
with zeta>0. Consequently

    Bprof=yb^3 T0/zeta                                 (2)

is smooth on the WHOLE closed native shell, including both endpoints and
eta in [-1,1]. At the inner edge it equals
yb^3 exp(b/yb^2) Aa; at the outer edge it equals
exp(a/x^2)(b0,yb^6 bz). Those expressions agree with (2) on their collars.
This does NOT say T0/zeta is bounded: its angular component retains the
cubic pole. It neither changes the target to the summed stress nor discards
the exceptional pressure terms used to construct the finalized T0.

Let M_j bound the factorial exact-degree two-variable norm J_j(Bprof) on
that compact shell. Each M_j is a finite FIXED profile seminorm, by (1)--(2)
and the inherited all-order smooth defining maps. No uniformity in frequency,
width tending to zero or derivative order is claimed. Just as in
`core_fixed_profile_tame`, numerical values are needed for a numerical norm
certificate, not for a conditional fixed-profile asymptotic recipe.

Since (partial_x^k/k!)yb^-3=binom(k+2,2)yb^(-k-3), define

    F_j=sum_(i=0)^j M_i binom(j-i+2,2).
    J_j(T0/zeta)<=F_j delta^(-j-3).                    (3)

This is sharper than obtaining the same quotient by separate weighted
derivatives of T0 and zeta^-1. It follows from the actual factorizations;
generic smooth zero extension without them would not prove (2) or (3).

## 2. Nonlinear native-to-band chart and the original scale factor

Use the SAME accepted physical branch
s-Z^2 s^(2h)=T, s=q/Q in [1/4,4], eta=Z s^(-D), D=1/2-h,
R>=Rmin>0. The fixed-profile chart foundation proves, for i+j>=1,

    |partial_Z^i partial_T^j log s|/(i!j!)
         <=(1/3)32^i 64^j.

Write U_j=(1/3)sum_(i=0)^j32^i64^(j-i) for j>=1 and U_0=0.
These bounds extend to one-sided physical time endpoints; only the coordinate
is analytically estimated, not the smooth profile at complex arguments.
The profile map is x=2log R-log(2Xa)-log s, eta=Z s^-D.
Positive normalized exponential recurrences give E_0=2, B_0=16 and

    E_j=(1/j)sum_(i=1)^j i(U_i/2) E_(j-i),
    B_j=(1/j)sum_(i=1)^j i(2U_i) B_(j-i),
    C_j=2/(j Rmin^j)+U_j+2E_j+E_(j-1), j>=1.           (4)

Here E bounds s^-D, B bounds s^-(1+h), since |Z|<=2, s^-D<=2,
s^-(1+h)<=16 and 1+h<2. C bounds the SUM of normalized increments of
both profile coordinates. All higher chart derivatives are retained.
For C(z)=sum_(j>=1) C_j z^j set

    Q_n=sum_(k=0)^n F_k [z^n]C(z)^k,
    T_n=sum_(i=0)^n B_i Q_(n-i).                        (5)

Finite multivariate Taylor composition and (3) prove

    J_n(s^-(1+h) T0/zeta)<=T_n delta^(-n-3).            (6)

Indeed each outer order k<=n can be padded to pole n+3. The target here
is exactly Tbeta=(Q/q)^(A+1/2)T0=s^-(1+h)T0. Its physical covariance
normalization, squared partitions and both sign columns remain those of
`core_localized_stress_assembly`. The root budget's existing target allowance
delta^(-3n-3) therefore accepts T_n without changing a root or a value gap.
These are constants in S, permitted as degree-zero positive scale polynomials.
The existing full R-to-U operator/chart transport is applied AFTER this
R-chart recipe; no differently normalized jets are silently reused.

## 3. Flat logarithm with both original constants

Let ell=log sqrt(zeta)=-a/(2x^2)-b/(2yb^2). For k>=1 its pure normalized
profile derivative is bounded by
(a+b)(k+1)delta^(-k-2)/2. Thus for n>=1,

    ell_n=sum_(k=1)^n [(a+b)(k+1)/2] [z^n]C(z)^k,
    J_n(ell)<=ell_n delta^(-3n).                        (7)

The parent exponential recurrence for sqrt(zeta) can use (7) directly.
`loghalf(0)=0` is an UNUSED norm slot in that recurrence, never a bound for
the unbounded logarithm itself. The width-squared a=t1^2 and outer b=4
are retained; neither exponential is replaced by a comparable weight.

## 4. Complete fixed cutoffs, with the slow partition outside H

The original normalized bump translates are smooth: their periodic squared
denominator has a strictly positive floor. For each fixed derivative order,
their normalized one-dimensional profile seminorms Kd_j (dyadic) and Km_j
(mesh) are finite uniformly in translates. Let Ktr_j, Ktm_j be the same
seminorms for the fixed transverse and relative temporal cutoffs. These
are properties of the already fixed profiles, not free new field estimates.

The dyadic argument is log2(q)+ell=log2(s), with integer ell frozen and
1/log 2<2. Replace C in (5) by the increment series 2sum U_j z^j to
compose Kd. The original mesh width is S^-3, so its three-fold tensor
product has normalized degree-n bound

    Mesh_n=S^(3n) sum_(i+j+k=n) Km_i Km_j Km_k.          (8)

No derivative of an integer mesh label is taken. In the native auxiliary
rectangle, the dual covectors lambda_r,lambda_t of (vr,vt) both have l1
norm below 2. The transverse argument is lambda_r(Y-center)/r0, while

    v/Ls=1/2+lambda_t(Y-center)/(ci Ls), ci Ls=2r0.

Therefore sufficient normalized cutoff factors are
Ktr_j(2/r0)^j and Ktm_j(1/r0)^j. The relative temporal cutoff has NO
extra S loss; its exact ci Ls cancellation must precede estimating the
separate factors. All fixed radius costs remain. Convolve these two factors
with the dyadic and mesh factors for Pi_n(S), and use n!Pi_n as a raw
maximum bound for the parent's homogeneous-wave product rule. Each Pi_n
has polynomial degree at most 3n. Physical graph curls and common covering
factors are still the existing operators/atlas conversions, not extra
ordinary cutoff derivatives in this native chart.

The complete localized coefficient has the slow partition eta_beta ONCE.
The factory includes its jet bound in the wave norm slot, while H remains
the ORIGINAL unrooted chi*psi Haar column with eta_beta outside H. No slow
partition is divided out, absorbed into the root equation, or multiplied
twice. The unlocalized chi*psi bound is also available explicitly. The exact
principal assembly still uses sum eta_beta^2=1 with both signs already in H.
All slow/transverse cutoff commutators remain in the regular complete curl
and residual; they are not declared Gaussian-flat.

## 5. Joint recipe, checks and scope

`from_original` derives T_n, ell_n and the full fixed-cutoff jets, then joins
the same canonical H/W and full D recipes on the ORIGINAL accepted tail.
It optionally transports every root input to U consistently. K0, beta,
Rroot, mass, width, physical roots and value/cone conditions stay fixed.
The complete-wave covariance/pressure maps are not replaced by the two
principal columns. The whole axial-moment input ledger retains background
through n+13 and fixed profile/cutoff seminorms through n+12 by this route.

Exact manufactured edge jets and fourth-order implicit chart expansions
test (3)--(7), including eta endpoints, nonlinear chain terms and the
original scale prefactor. Independent polynomial cutoff expansions test
(8), the dual-covector factors and ci Ls cancellation. Original-gate joins
check all five root inputs and complete-curl/signed norm recipes. These
finite controls do not evaluate the actual giant profile or its seminorms.
Eight historical frozen records and their declared source hashes are bound.

Next substitute these explicit primary inputs into the existing complete
covariance/mean and scalar/source joins on the SAME weighted physical hull;
audit which linked input, value-domain or full residual inequality still
remains. Initial weighted Haar-row and later conditional mean closures
already proved in their checkpoints are not reopened here. Validity of the
original full background/phase jet polynomials, defining-map/profile lemmas
and original physical value/cone neighborhood remains conditional. All six
covariance entries, fresh axial velocity, moving pressure moments, Pi1/Pi2,
reserved supports, viscosity, cutoffs and interactions remain in the full
ledger. There is no unconditional pulse, complete physical residual-decay,
infinite realization or blow-up theorem from these norm substitutions.

Reproduce with `python test_audit.py -v`, then
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
Install the declared SymPy requirement if needed. The analytical fixed-order
continuum implications are (1)--(8), not the passing finite controls.
