# Q79: explicit full-background and ambient value bridge

Parent: `8fc7b66adfca0b6a269e21e648bfa90a3c85bac5` (Q78).
The leading-field derivative bounds are now available through the terminal
stress edge. This checkpoint joins their low orders to the full-background
comparison and makes Q23's finite-head **ambient** bound explicit, without
an unevaluated maximum over its labels. It also reduces the ambient VALUE
extension's full-background request from the generic second-jet recipe to
first jets by using the actual affine normal. It does not reprove Q23's
existence conclusion or infer column positivity from an upper ODE bound.

The statements below are conditional on the same finalized coefficient
seminorms, normalized Borel schedule, leading spectral gaps and chart
geometry as the cited providers. No physical band, length, root, pressure
gauge, support or C^-20 switch width is changed. The new threshold is used
only to partition an estimate of the existing family. Historical evidence
is unchanged; executable rational examples are manufactured controls.

## 1. Low-order data on one fixed positive physical hull

Keep R in [r,H], 0<r<=H, s in [1/4,4], |eta|<=1 and q<=qstar<=1.
The hull must include the full radial fibers used by Q26, not just one
mesh box. Its leading X range lies in [r^2/8,2H^2]. Use Q63--Q78's
modulated, exactly repaired and exterior bounds on the stress shell;
where this larger hull extends beyond it, retain the inherited analytic
core and real heat tail. The latter has one-sided eta jets at its poles,
not a claimed complex analytic continuation. These finite profile bounds
and the finite real chart composition of
`core_background_finite_head_budget` give constants J0,J1,J2 for the
leading F,G slow jets. No pressure trace is reset in obtaining them.

Let M_(n,m) be the COMPLETE coefficient tuple seminorm, K_(m+1) the
actual Borel cutoff Euler seminorm, and keep exactly its retained schedule

    C_(n,m)=max(1,2 K_(m+1)(1+2nh)^m M_(n,m)),
    C_(n,m) a_n^(-nh)<=2^-n, m<=n.                    (1)

These are assumptions about the globally finalized coefficients, including
viscosity, radial curl commutators and full pressure/moment memory. Bounds
for the leading E0,U0 alone do not prove (1). Put c1=C_(1,1),
B1=c1+1/2 and Jchart=CF_1+CG_1 in the inherited chart convention.
Its finite-head theorem gives, simultaneously on every original band,

    ||F-F0,G-G0||_(C^1) <= epsilon^2 E1,
    E1=4 Jchart B1,
    ||F,G||_(C^1) <= Jfull1=Jchart(Lead1+B1).          (2)

Here the norms are component-sum upper bounds for all raw slow derivatives
through the indicated order. Lead1 need only bound leading E0,U0 when
using the F,G entries of the recipe. No conclusion about its unused
btilde entry is licensed by pair-only leading bounds;
F0,G0 in (2) mean leading F,G at the SAME point. For leading values at a
frozen representative and max-coordinate box radius Dbox S^-3, Q22 gives

    Db*=3 Dbox J2+E1,       Mz*=Jfull1.               (3)

Thus (3) uses leading jets through two, finalized coefficient ONE through
one, cutoff Euler jets through TWO, and the whole schedule (1). For higher
fixed output n, Q61 still requests corrections through n+6 and the finite
head n_coefficient<max(2,n+6); this value-order bridge does not erase them.
It does not by itself bound the radial background needed by the eikonal.

## 2. Explicit tail floor and tangent estimate, without reselection

Insert (3) as ESTIMATE constants in `frame_column_foundation`, retaining
its fixed leading R,F,g,c,lambda,Bs floors and upper bounds. Denote the
result by Dr,Dt,Drprime,Dtprime,Ds,CE*,Cd*,S0. Let

    L*=max(1,2Dt/bmin,Ds),
    ell0=ceil(max(16,360(4/h)^5,L*/(2r0))).             (4)

Q23's exponential-versus-power proof gives the original slow gate and
L>=L* for every existing ell>=ell0, since S=ell^2 and
2r0 S<=L<=S. No Riccati gate is needed for this upper estimate.
The actual full and tangential normal lengths are at least nu_t=bmin/2.
With

    Nt=bmax(1+3U/2)+Dr+Dt,
    Mt=Drprime+Dtprime,
    Tt=Kframe exp(CE*+Cd*),
    Kframe=(1+S0)max(1,cmax)/min(1,cmin),              (5)

we have |n|<=Nt, |nprime|<=Mt/L<=Mt, and the physical tangent two-time
bound Tt P(v)/P(w). These are the exact inherited frame/energy conclusions
with new analytical majorants, not a claim that the original numerical
selector used these constants. All fixed spectral hypotheses remain.

## 3. The finite head has an explicit ambient energy bound

On ell<ell0, nearest NONZERO rounding gives |p_theta|>=1/k. Define

    kh=1+2^ceil(h ell0/2),  Lh=ell0^2,
    nu_h=1/(kh H).                                    (6)

Indeed k=ceil(2^(h ell/2))<=kh. Keep the giant power in (6) as a symbolic
finite constant; it need not be allocated. Let qtheta,qz,qrad be the
carrier upper bounds in `core_phase_coefficient_budget` (1), using a valid
common Lmin (one suffices when all retained L>=1). From (2), put

    M=2(qtheta+qz)Jfull1,
    K=(H+5)Jfull1,
    Nh=qrad+qtheta/r+qz+Lh M,
    Dref=lambda_max(1+9U^2/4).                         (7)

The component sums in (7) bound |nprime|, ||Kbg|| and |n| on every
head label: n=(x0-v Hphase_R,p_theta/R,pz-epsilon v Hphase_Z),
Hphase=p_theta F+pz G, epsilon<=1. This does not sample any labels.

For N=n.n, Pi=n n^T/N, m=nprime, d=epsilon k^2 N and p=P'/P,
the canonical ambient generator is exactly

    Ahat=-(I-Pi)Kbg-n m^T/N-d(I-Pi)+p Pi,
    Ahat-pI=-(I-Pi)Kbg-n m^T/N-(d+p)(I-Pi).           (8)

Its normal component satisfies (n.t)'=p(n.t). Since p=lambda-dref,
lambda>=0 and d>=0, the last term in the Hermitian energy estimate is
bounded ABOVE by dref|t|^2, not by (d+|p|)|t|^2. Orthogonal projection
has norm one. Therefore the logarithmic norm of Ahat-pI is at most

    K+M/nu_h+Dref,   Th=exp[Lh(K+M/nu_h+Dref)].        (9)

Gronwall gives ||Vhat(v,w)||<=Th P(v)/P(w) on EVERY head subinterval.
The same estimate applies to complex vectors by the real Hermitian
pairing. Actual damping is retained; its large positive size costs
nothing in (9). Dropping the moving-normal term would make (8) false.
No numerical minimum of P and no d0 or second background jet is required.

## 4. Exact affine lift on the tail

Write hvec=n/N. The actual path derivative is

    hvecprime=(I-2Pi)m/N,
    |hvecprime|=|m|/N,
    A0 hvec-hvecprime=-(I-Pi)(Kbg n+m)/N,             (10)

where A0=-(I-Pi)Kbg-n m^T/N. The first matrix in (10) is an orthogonal
reflection. Thus a generic total-jet bound for hvecprime is unnecessary:
the affine path only uses first background derivatives, even though slow
derivatives of the ambient generator still require higher background jets.

For an arbitrary ambient solution set a=n.t and z=t-hvec a. Then

    a(v)=P(v)a(w)/P(w),
    z'= (A0-dI)z+(A0 hvec-hvecprime)a,
    z(w)=(I-Pi(w))t(w).                               (11)

Both z and its forcing are tangent. Apply (5) to (11); (10) bounds the
forcing coefficient by Qt=K/nu_t+Mt/nu_t^2. The orthogonal initial
projection has norm one. It follows that, with S=ell^2,

    ||Vhat(v,w)|| <= At(S) P(v)/P(w),
    At(S)=Tt[1+S Qt Nt]+Nt/nu_t.                      (12)

This is a genuine TWO-time estimate for arbitrary ambient initial data,
not just the canonical growing vector. Combining (9),(12) gives on all
original accepted bands the explicit polynomial

    A(S)=Th+Tt+Nt/nu_t+S Tt Qt Nt.                    (13)

The uniform normal and tangential squared reciprocal is
max(nu_h^-2,nu_t^-2). This may be much worse than the original normal
constant: it does not certify that old constant or its selector inequality.
The actual damping floor is
min(nu_h^2,nu_t^2), since epsilon k^2>=1. Omit head constants when that
set is empty. Finite constants in (6),(9) may be enormous but are fixed,
independent of the running S and of the derivative order subsequently
requested. Exponentiating a polynomial in the running S on the infinite
tail would not be an acceptable substitute for (12).

## 5. What this supplies, and the precise remaining input

Equations (2),(6)--(13) supply Q61's low-order background difference,
normal/entrance value denominators and ambient value-polynomial recipes,
relative to the stated coefficient, schedule, geometry and spectral data.
They do not imply Q61's normalized column inverse K0, positive square
root, or any source/support claim. Q25's direct canonical cone proof and
Q26's fixed-q domain join remain the route on their original certified
domain; a finite-head ODE upper bound cannot extend a cone margin.

The remaining coefficient input can be stated without a hidden supremum:
on the full profile hull prove the stated M_(1,1) bounds for the globally
finalized tuple (E1,U1,W1^r,J1^r), and prove (1) for the SAME schedule
at each required order. Until that and the leading value gaps are supplied,
(2),(13) are conditional majorants, not an actual all-scale pulse
certificate. Derivative bounds for E0,U0 do not imply either a positive
spectral gap or these positive-order coefficient estimates. No replacement
pressure cutoff or local-core coefficient is allowed in this substitution.

Next quantify that finalized coefficient-one tuple and its global moment/
reserved-support identities on this hull, then join the original column
and root domain. Higher n+6 inputs, complete-wave means with axial velocity
and pressure, radial pressure tails, viscosity/cutoff/interactions and the
complete residual remain open. No unconditional pulse or blow-up claim.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. Exact controls check the moving-normal
lift independently, the damping-aware energy bound, finite-head rounding
and the inherited low-order request set. They are not a continuum proof.
