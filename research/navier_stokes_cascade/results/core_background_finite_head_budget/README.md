# Explicit finite-head background norms and the unchanged normal-floor join

Parent: `84a99f5292c1f63037275db47cad301e9717cb8c`.
The preceding checkpoint reduces frozen-phase coefficients and the inverse's
normal floor to full-background jets and a comparison constant Cphi.
The existing background jet foundation proves these jets finite but permits
an unspecified enlargement on a remaining compact q interval. This checkpoint
replaces that enlargement by an explicit finite coefficient head plus a
geometric tail, composes it with the actual annular chart, and supplies Cphi.
It reduces the remaining selector join to inequality (12) below.

All formulas are relative to the actual leading-profile regularity/gaps,
finalized coefficient seminorms and the retained normalized Borel diagonal.
The code evaluates rational MANUFACTURED inputs only; actual seminorms and
their selector substitution have not been evaluated. No profile, carrier,
cutoff schedule, reserved support, accepted band or width C^-20 is changed.
Historical evidence remains frozen. No full residual or blow-up claim follows.

## 1. Retain the complete coefficient problem and radial curl term

Use E=q partial_q at fixed X,eta, A=1/2+h, D=1/2-h and 0<h<1/100.
On the existing enlarged positive annulus, write the normalized velocity as

    Ecal=q^A utheta = E0+sum_n q^(2nh) chi(a_n q) E_n,
    Ucal=q^A uz     = U0+sum_n q^(2nh) chi(a_n q) U_n,
    Wcal=q^(1/2) ur = W0+sum_n q^(2nh)
                         [chi(a_n q) W_n^r-(E chi)(a_n q) J_n^r], (1)
    W_n^r=V_n/Rprof, J_n^r=2 eta M_n/(Lprof Rprof).

Here Lprof is the fixed profile-coordinate denominator; it is distinct from
the pulse path length L below. These are the particular Stokes potentials
and the full globally finalized coefficients of the inherited construction.
The coefficients include forward radial/axial reconstruction, viscosity,
all lower products, full Pi1/Pi2 pressure memory and the five moment repairs.
They are not compact pressure cutoffs or local-core coefficients. In
particular M_n and V_n may retain memory on the reserved annulus even where
U_n=0. No new pressure or moment identity is inferred from a velocity norm.

Let M_(n,m) bound every X,eta derivative through total order m of the
COMPLETE tuple (E_n,U_n,W_n^r,J_n^r), and let K_j>=1 bound all Euler jets
through j of the actual Borel chi. These are distinct from the radial
coefficient cutoff kappa(X). The background foundation's exact product rule
and its radial commutator give

    C_(n,m)=max(1,2 K_(m+1)(1+2nh)^m M_(n,m)),
    |E^j partial_X^a partial_eta^b increment_n|
       <= C_(n,m) q^(2nh), j+a+b<=m.                         (2)

The extra cutoff derivative m+1 is essential. The input convention is Euler;
q^j partial_q^j has a Stirling/Bell conversion cost and cannot be substituted
with the same constants. The API rejects that convention.

Assume the SAME retained schedule satisfies

    C_(n,m) a_n^(-nh)<=2^-n, m<=n,                          (3)

alongside every original physical/weighted-stress schedule requirement.
The foundation gives an admissible symbolic construction, not evaluated
numerical a_n. This checkpoint neither chooses another schedule nor infers
(3) from an unconverted frozen-derivative diagonal. Wherever a summand or
its cutoff derivatives are nonzero, q<=a_n^-1, so (2)--(3) imply
2^-n q^(nh) there; all cutoff jets vanish outside. Thus the estimate covers
the entire interval, including transition collars.

## 2. A finite explicit constant on every 0<q<=qstar<=1

For each fixed m take N=max(2,m), and define

    Bcorr_m=sum_(1<=n<N) C_(n,m)+2^(1-N).                    (4)

For the finite head, q^(2nh)<=q^(2h). For n>=N, use (3) and

    sum_(n>=N) 2^-n q^((n-2)h)
      =2^-N (q^h)^(N-2)/(1-q^h/2)<=2^(1-N).                 (5)

Hence the normalized correction tuple has every Euler/X/eta jet through m
bounded by Bcorr_m q^(2h) on the WHOLE existing q interval. The same bounds
justify termwise differentiation at each fixed order by uniform summation.
Coefficient ONE must remain in the finite head even for m=0,1: starting
the diagonal tail at n=1 gives q^-h after division by q^(2h), which is
unbounded. This is why N=max(2,m), not N=m or max(1,m).

Let Lead_m bound the leading tuple jets through m in the same norm.
Then every normalized full-background jet through m is bounded explicitly
by Lead_m+Bcorr_m. There is no unspecified compact-interval padding or
order-dependent shrinkage of qstar. Inputs for higher fixed orders are
finite, though their values may be very large.

## 3. Finite real jet composition into the retained slow chart

Use the original chart s=q/Q in [1/4,4],

    s-Z^2 s^(2h)=T, X=R^2/(2s), eta=Z s^-D,
    Rmin<=R<=Rmax, Rmin>0, |eta|<=1.                        (6)

All bounds refer to the same enlarged domain required by the parent; the
line segments used below must lie there. The coordinate majorant in
`core_fixed_profile_tame` gives the zero-constant log jet majorant

    Y(z,t)=(1/3)[(1-32z)^-1(1-64t)^-1-1].

Let Y_m be the sum of its coefficients through degree m, and let
Exp_m(a)=sum_(k=0)^m a^k/k!. Set

    Rpoly=Rmax^2+1_(m>=1)2Rmax+1_(m>=2),
    Xdelta=2[Rpoly Exp_m(Y_m)-Rmax^2],
    etadelta=(1+4*1_(m>=1))Exp_m(Y_m/2)-1,
    Comp_m=sum_(i+a+b<=m) Y_m^i Xdelta^a etadelta^b/(i!a!b!),
    Recip_m=sum_(j=0)^m Rmin^(-j-1),
    CF_m=m! 4 Exp_m(Y_m) Recip_m Comp_m,
    CG_m=m! 4 Exp_m(Y_m) Comp_m,
    CB_m=m! 2 Exp_m(Y_m/2) Comp_m.                           (7)

To see (7), write ell=log(s/s0). Then X=X0 exp(-ell) with
R replaced by R0+delta_R, and
eta=eta0 exp(-D ell)+s0^-D delta_Z exp(-D ell).
Use 1/(2s0)<=2, |eta0|<=1, s0^-D<=4, D<=1/2. The displayed Xdelta
and etadelta bound the zero-constant coefficient sums. Compose the
profile's finite Taylor jet in (ell,X-X0,eta-eta0), then multiply by
s^-A/R, s^-A or s^-1/2 respectively. Their base sizes are at most
4/Rmin,4,2 and their log exponents have absolute values <=1,1,1/2.
Multiplication and the finite Taylor sum give normalized coefficients;
m! bounds conversion of ANY mixed coefficient of degree <=m to its raw
jet. Only the coordinates use analyticity. The actual profiles need only
C^m real regularity, with boundary jets interpreted by continuity.

For F=Q^A utheta/R, G=Q^A uz and btilde=Q^(1/2)ur, sufficient RAW slow
jet bounds through m are therefore

    F_m=CF_m(Lead_m+Bcorr_m),
    G_m=CG_m(Lead_m+Bcorr_m),
    btilde_m=CB_m(Lead_m+Bcorr_m).                           (8)

A shared leading F/G bound is Jlead_m=(CF_m+CG_m)Lead_m;
a shared full bound is Jfull_m=(CF_m+CG_m)(Lead_m+Bcorr_m).
Since q^(2h)=epsilon^2 s^(2h)<=4epsilon^2, the difference from the
leading F/G is bounded by epsilon^2 Jerr_m, with
Jerr_m=4(CF_m+CG_m)Bcorr_m. These estimates retain the radial contribution
in (1) for btilde and hence for the full parent's eikonal recipe.

## 4. Explicit Cphi from the actual frozen reference

Retain U>=u, lambda0<=lambda_max, |g0|>=gmin>0 and
Llow S<=L<=Lhigh S. The slow box has max-norm diameter Dbox S^-3 in
(R,Z,T), and the enlarged chart contains its reference-to-point segments.
Write B=max(1,lambda_max), qz=B[1+U/(Llow gmin)], p0=Rmax qz,
pplus=p0+1. These are outward bounds, not replacement carriers.

At the frozen representative R0 the ACTUAL unrounded carrier obeys

    (ptilde/R0,pz)=Bs[Kunit-sigma*u*g0/(L|g0|^2)],
    ptilde F0_R(rep)+pz G0_R(rep)=-sigma Bs u/L,
    x0=sigma Bs u/2, |p-ptilde|<=1/k.                       (9)

The nearest NONZERO rounding requires a unit error, not half a unit.
Use (8), the leading second-derivative bound on each box segment, and the
parent's exact full normal. Define seven nonnegative coefficients

    ctilt=2 B U/(Llow gmin),
    crbox=Lhigh(pplus+qz)3 Dbox Jlead_2,
    crbg=Lhigh(pplus+qz)Jerr_1,
    crround=Lhigh Jlead_1,
    ctheta_round=1/Rmin,
    ctheta_box=p0 Dbox/Rmin^2,
    cztransport=Lhigh(pplus+qz)Jfull_1,
    Cphi=sum of these seven coefficients.                  (10)

The first term bounds BOTH tangential tilt shifts by ctilt/S. The radial
normal errors are at most crbox/S^2+crbg S epsilon^2+crround S/k.
The angular errors are at most ctheta_round/k+ctheta_box/S^3.
The full axial transport is at most cztransport S epsilon. Thus, for S>=1,

    |n-Bs(s_sigma,Kunit)|
      <= Cphi[S^-1+S^-2+S(epsilon+epsilon^2+k^-1)].          (11)

This is a component-sum bound, hence also a Euclidean bound. No angular
component lower bound is inferred. Cphi uses only leading tuple jets
through TWO, finalized coefficient ONE tuple jets through ONE, actual
Borel cutoff Euler jets through TWO, and the normalized schedule (3),
besides the fixed shell, box, growth and shear data. The small coefficient
head does not mean the other coefficients or their pressure disappear:
they are controlled by the whole tail (5).

## 5. Join to the existing selector: a remaining explicit inequality

The parent requires S^2(epsilon+epsilon^2+k^-1)<=1 and
36 Cphi^2<=beta S^2, beta=lambda_min/[4(1+U^2)^2].
The SAME `core_primary_tail_recheck` selector already imposes

    Cslow ell^4(epsilon+epsilon^2+k^-1)<=tau,
    ell^2>=Lreq/(2 r0), S=ell^2, Cslow>=1, tau<=1.

Its first gate implies the first parent gate on the entire accepted half
line. Its unchanged length gate implies the second if

    beta Lreq^2 >= 144 r0^2 Cphi^2.                         (12)

This is a compatibility test on the EXISTING Lreq, not authorization to
increase Lreq or select another band. Failure is reported explicitly.
All original selector gates and all earlier domain premises must still
hold. Once actual inputs establish (12), the parent's triangle inequality
gives |n_tan|^2>=beta/4 and |n|^-2<=4/beta on the same accepted domain.
Until then the floor is conditional; a passing manufactured selector
fixture does not discharge the actual substitution.

`phase_inputs` transfers (8) as constant positive polynomials into the
parent's raw coefficient recipes. The latter still retain varying growth,
full axial transport, radial flow and unprojected source/pressure terms.
The sourced inverse additionally requires its original tangent propagator,
source and path-coordinate derivative polynomials; these are not derived
by (12). The complete signed-prefix covariance, fresh pressure/moment
means and complete physical residual remain their separate obligations.

## 6. Evidence and next obligation

Reproduce using `python test_audit.py -v` and `python audit.py
--verify-record evidence.json --output /path/to/fresh.json`.
The new controls use exact rational finite jets, exact geometric sums,
manufactured normal errors and the original selector fixture. The audit
binds six historical evidence records and all their declared sources.
Inherited background numerical diagnostics remain numerical diagnostics,
not continuum proofs of (1)--(12) or leading-profile existence.

Next: supply certified leading through-two jets, finalized coefficient-one
through-one seminorms, cutoff Euler norms and the existing normalized
schedule/gap/geometry data; evaluate (10)--(12) against the unchanged
accepted selector. Then discharge the homogeneous tangent propagator and
remaining actual-root/support/fixed-map premises before promoting the
first sourced complete-wave/mean estimate. No unconditional actual pulse
bounds, complete residual improvement or blow-up have been proved.
