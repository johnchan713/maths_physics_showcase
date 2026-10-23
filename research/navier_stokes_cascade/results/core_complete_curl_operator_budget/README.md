# Complete curl coefficient budget, including radial-chart conjugation

Parent: `95da5280b58140f76fd390b07e27940d6a455d5a`.
Status: `conditional-complete-curl-operator-and-radial-chart-budget`.
The parent CI run 40030983135 was in progress without failed steps on
23 October 2026; its canonical Haar-column step had passed. Historical
evidence, physical fields, roots, radius, widths and accepted bands remain
unchanged. This derives the independent D_j input still retained in the
root/complete-curl budget. It does not prove the full physical residual.

## 1. Exact cylindrical curl before estimating coefficients

In the same frozen native band use coordinates (R,Z,T,Y1,Y2) and the
complete normal n, N=|n|^2, h=n/N. The primary coefficient A includes its
original roots and every slow/transverse/temporal cutoff, and n.A=0.
The real potential is -C sin(k Phi)/k with C=h cross A. Coefficients
are theta-independent, but the cylindrical basis connection remains:

    B=curl_* C
     =er cross Dr C + ez cross epsilon partial_Z C + ez C_theta/R,
    Dr=partial_R+M g(R)(vr1 partial_Y1+vr2 partial_Y2),
    g(R)=d_r R^(d_r-1), vr=(1,-(sqrt(2)-1)).             (1)

In components this is (-epsilon C_theta,Z,
epsilon C_r,Z-Dr C_z, Dr C_theta+C_theta/R). A component inverse n_theta^-1
is never used. There is no T derivative in a spatial curl, although all
mixed T derivatives of its coefficients remain in the jet bounds.

Set eta=epsilon^kappa. These frozen band constants have zero continuous
derivatives. The retained scale identity M<=epsilon^-kappa S^-rho gives
eta M<=1; also eta<=1 and eta epsilon<=1. Thus no negative epsilon power
remains in the normalized operator eta B. Define vectors

    a_R=eta er, a_Z=eta epsilon ez, a_T=0,
    a_Ymu=eta M g vrmu er.

Let cross(q) denote the matrix q cross. Expanding (1) EXACTLY gives

    eta B=D0 A+sum_mu Dmu partial_mu A,
    Dmu=cross(a_mu)cross(h),
    D0=sum_mu cross(a_mu)cross(partial_mu h)
            +(eta/R)ez etheta^T cross(h).              (2)

The original moving normal and full background are in h and all its
derivatives. The graph acts on the ENTIRE h cross A, not just on A.
The original graph operator is a coefficient times a derivative, so it
does not generate a spurious derivative of g at order zero. Derivatives
of Dmu and D0 DO differentiate g; those terms are retained below.

## 2. Explicit positive raw mixed-jet coefficients

Use the parent's full-normal lift recurrence h_j, bounding the Euclidean
norm of every raw mixed derivative of total order j. It is derived from
the complete PhaseBudget normal and the reference-scaled normal reciprocal,
not an independent potential input. On the original closed positive hull
Rmin<=R<=Rmax, set

    Bpow=max(1,Rmax^2,Rmin^-1),
    g_j=3 Bpow (j+1)! Rmin^-j,
    r_j=j! Rmin^(-j-1),
    V_R,j=2*1_(j=0)+2 g_j.                             (3)

The original auxiliary geometry gives 1/3<d_r<5/2. Indeed rho>1/4,
rho<1, 0<hconstruction<1/100 and kappa=1/100000 imply
1/2-2/(100*100000)>1/3 and d_r<2(1+1/100)<5/2.
For alpha=d_r-1 in (-1,3/2), R^alpha<=Bpow and
product_(k=0)^(j-1)|alpha-k|<=(j+1)!. Thus (3) bounds every derivative
of g. Only pure R derivatives of g and 1/R are nonzero. Since
|vr1|+|vr2|<2, V_R,j bounds the SUM of all derivative-direction scalar
coefficient jets, including radial, axial and both graph directions.

For induced matrix infinity norm, ||cross(h)||<=2|h|_2, whereas cross(er)
and cross(ez) have norm one. Raw convolution is
(a*b)_j=sum_i binom(j,i)a_i b_(j-i). By differentiating (2),

    Draw_R,j=2[(V_R*h)_j+(V_R*h_shift)_j+(r*h)_j],
    h_shift,j=h_(j+1),
    D_R,j=(5^j/j!) Draw_R,j.                           (4)

D_R,j bounds the SUM of the factorial exact-degree J_j norms of D0 and
ALL first-order matrices. The factor 5^j/j! converts raw mixed maxima
in the same five coefficient variables. The shifted term is necessary:
D0 contains first derivatives of h. Every graph, cylindrical connection,
normal reciprocal and coefficient-derivative term is present in (4).
All quantities are finite positive scale polynomials at each fixed order.
Their constants may depend on the fixed radial hull and original profiles;
there is no assertion of uniformity at the axis Rmin=0.

## 3. Re-express in U=R^d_r consistently

The atlas also uses U=R^d_r. This is a coordinate re-expression of the SAME
field on a positive radial hull, not a different scale or pulse. One must
first use partial_R=g partial_U. Thus the first-order matrix D_R in (2)
becomes g D_R; the zeroth-order D0 does not change. BEFORE composition
with R(U), replace V_R in the FIRST term of (4) by

    V_U,j=1_(j=0)+3 g_j.                               (5)

This bounds g a_R, a_Z and both a_Y directions. The D0 term continues
to use V_R and h_shift. Reusing D_R unchanged in the U chart is wrong,
even if the coefficient functions have already been composed with R(U).

For a quantitative finite-jet composition recipe, put

    Ufloor=min(1,Rmin^3),
    L_j=Rmax (j+2)!/(2 Ufloor^j), j>=1,
    a_j=1_(j=1)+L_j/j!, a_0=0,
    Comp_n(F)=sum_(k=0)^n F_k [z^n](sum_(j=1)^n a_j z^j)^k. (6)

Because 1/d_r in (2/5,3), U>=Ufloor and R<=Rmax, the falling-factorial
formula bounds |partial_U^j R(U)| by L_j. The four identity coordinates
have only their first increment. Every coordinate increment has normalized
degree-j sum bounded by a_j; finite multivariate Taylor/Leibniz therefore
gives J_n(F composed with R(U))<=Comp_n(F). This uses smooth finite jets,
not analyticity of the background or convergence of an infinite series.

Apply (6) to the normalized precomposition version of (4) with (5) to get
D_U,n. Also compose EVERY H, target, wave and loghalf input jet of the
root budget. The value gaps K0, beta and Rroot are unchanged by this
positive coordinate map, as are mass, width and actual roots. Existing
edge losses are nondecreasing in outer derivative order k<=n; padding
them to order n preserves their original form and the same pulled-back
zeta, delta and Gaussian weights. No second edge unweighting is applied.
The code supports explicit R and U norm recipes and rejects mixed chart
names. Auxiliary common-cover factors remain the existing atlas conversion;
they are not silently charged twice here.

The signed-target jet callback xi is supplied SEPARATELY to `RootBudget.signed`.
It must also be transported with `operator.transport(xi,'U')` before that
call. Recharting only the stored root callbacks does not rechart a later
external target automatically. The same rule applies to any later external
source/target input; its established edge powers can be padded to the
requested outer order as above. The API makes this extra obligation explicit.

## 4. Join and derivative ledger

`from_original` reuses the preceding canonical Haar-column/entrance factory,
its original normal/tangent/length and fixed-positive-tilt Gaussian gates.
It derives D_j from that SAME full normal lift and then supplies H, W and
D to the existing root/curl recipe. Only target/loghalf, fixed cutoffs and
the original value/cone hypotheses remain their retained lower inputs.
The independent graph/curl-operator callback is removed relative to those
premises. No physical root or primitive is rebuilt.

The parent's complete-curl polynomial consequently remains

    C_n=sum_(j=0)^n D_j[A_(n-j)+(n-j+1)A_(n-j+1)].       (7)

Its Gaussian and shell-edge losses, common-atlas normalization, scalar
moment and full pressure chain are unchanged. D_j needs normal jets through
j+1. The PhaseBudget normal at order q uses full background through q+1,
so D_j needs full background through j+2. The whole moving axial-pressure
moment of output order n uses D through n+11 and amplitude inputs H/W
through n+12, hence full background through n+13 by this route. Chart
composition uses only outer orders at most the requested order; it adds
no derivative-order loss. This is a sufficient ledger, not a sharpness claim.

## 5. Evidence and remaining obligations

Exact manufactured cylindrical/graph fields verify (1)--(2), the complete
divergence-of-curl identity and the nonzero defects from dropping graph,
axial, connection or potential-derivative terms. Independent rational
whole-box enclosures test the SUM of every matrix jet through order2.
An exact U=R^2 re-expression detects missing operator conjugation. Finite
signed Taylor controls test (6) through order4. Original-tail controls
join the same H/W/D recipes through the inverse, root, signed map and full
curl recurrence in both charts. These are finite manufactured controls,
not actual background suprema or numerical evaluation of the giant band.
Eight historical frozen records and their declared sources are preserved.

Next justify the remaining weighted target/loghalf and cutoff jet inputs
from the actual fixed stress/chart profiles, retaining the original value
and cone domain, or identify their precise blocking inequality. Full
background/phase primitive validity still remains a premise; deriving an
operator norm recipe is not evaluation of those suprema. Complete-wave
covariance means retain all six tensor entries, axial velocity, moving
pressure moments, Pi1/Pi2 memory and reserved supports. Forced pressure
still includes its unprojected n.f term; this curl calculation does not
change it. Every viscosity/cutoff/interaction residual and infinite
realization must still be joined before unconditional pulse bounds,
complete physical residual improvement or blow-up can be asserted.

Reproduce with `python test_audit.py -v`, then
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
Install the declared SymPy requirement if necessary. The continuum
conditional argument is (1)--(7), not the passing finite tests.
