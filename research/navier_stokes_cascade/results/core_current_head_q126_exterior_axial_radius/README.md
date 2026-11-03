# Q126: a fixed angular radius for the same exterior axial pulse

Parent: `86477ffed871c9487bbd3fec8bf4230d537a380b` (Q125).
On 3 November 2026 the freshly fetched branch agreed with this parent.
Q125 research CI 41649355593 and field CI 41649355674 were in progress;
Q124 field CI 41641232330 and Q123 field CI 41631558011 had passed.
The latest notes, frozen evidence, Q77/Q78, global_axial_angular_bound and
root_foundation were read first. Historical evidence, parameters, real roots,
frequency, supports, pressure datum and scales are unchanged. The original
core_growth_activation result remains local leading core/first collar with
width C_setup^-20.

Q77 supplied all-order recurrences for the exterior axial strength and the
PRE-HEAT future energy. Here we replace their high-angular-order inputs by
a fixed complex-eta disk and a closed geometric estimate for the SAME real
pulse. This is conditional on the retained exact stage formulas, value
bounds and matching identities; no sampled giant profile is substituted.
The construction concerns an exterior SLOW PROFILE pulse, not a canonical
physical high-frequency wave. The Gamma heat field is not extended to
negative or complex Z by this argument.

## 1. Parameter caps and the angular reference repair

Keep T>=64, lambda=exp(-4T), beta=1-lambda, Tu=30 log(1/lambda).
Monotonicity of T exp(-4T) and e>2 give the convenient universal caps

    lambda<=L:=2^-256, T lambda<=2^-250,
    30lambda log(1/lambda)=120Tlambda<1.             (1)

At every real eta0 in [-1,1], use complex disks |eta-eta0|<1/2. Then
Re(1+eta^2)>=3/4, |1+eta^2|<=13/4 and the original principal branch
of a=log(2/(1+eta^2)) is holomorphic with |a|<=4. The exact Q77 angular
target has the disk bound

    Z=5 exp(4) exp[-beta(Tu-3)]/beta
       <=10e^8 lambda^30 <2^15 lambda^30.           (2)

The last inequality follows from e<11/4 and 10(11/4)^8<2^15.
A rational proof of e<11/4 is the sum through order four, 65/24, plus
the tail bound (1/120)/(1-1/6)=1/100. The original interpolation step is
evaluated only at real radial arguments. No complex radial cutoff is used.

The angular equation has the SAME constant B,Q, inverse bound <=5 and
bilinear bound <=32. These estimates also hold for complex coefficients:
B is bounded by absolute row sums, and Q's bound is the absolute positive
bump-integral estimate in root_foundation/pulse_moment_audit. The retained
real root has value norm <=7000lambda^30. Put

    r=2^19 lambda^30 <=1/1024.

On the complex value ball ||u||<=r, the map -B^-1 z-B^-1 Q(u,u) has
image norm <=5*2^15lambda^30+160r^2<r and Lipschitz constant
320r<1/2. Its iterates converge uniformly on the disk, so the fixed point
is holomorphic there. On the real interval the accepted root lies inside
this ball and solves the same equation; uniqueness identifies it exactly.
No new root or high-order smallness condition is selected. Consequently

    |u_i^(j)(eta0)|/j! <= r 2^j.                    (3)

This replaces the angular reference-coefficient recurrence only in this
exterior stage. It says nothing about parameter analyticity of compact
variance cutoffs or the compact leading repair.

## 2. Complex bound for the complete PRE-HEAT future mass

Use Q77's exact positive measure at eta=0 AFTER removing the two angular
bumps. Its total mass is <=512T. The unedited normalized density at any
eta has relative factor (1+eta^2)^(2(1-theta)), 0<=theta<=1. On the same
disk its modulus is <=16, using the original logarithm branch. Restoring
the two relative bumps, each with value <=64, multiplies it by
(1+sum u_i beta_i)^2, of modulus at most (1+128r)^2<2. Therefore

    |R(eta)| <=2^14 T                                (4)

for the COMPLETE pre-heat future energy, including interpolation, angular
corrections, release, wait and infinite tail. Fixed eta-independent stage
endpoints and the finite dominating positive measure justify holomorphic
parameter integration. No infinite interval is bounded by length times
a nondecaying pointwise constant, and no bump square is omitted.

The earlier heat compensation and terminal heat edit have exact canceling
S changes as functions of real eta and are independent of trial strength.
Thus the PRE-HEAT scalar energy equation still selects the actual strength
after both edits. This identity does not make the actual partial pressure,
actual post-compensation S trace, or heat-modified future mass equal to its
pre-heat counterpart. It does not give a complex Gamma integral at Z=0.
The complex function R in (4) is the analytic extension of the pre-heat
reference used in the scalar root equation only.

## 3. Same positive axial strength, with all affine end energy retained

Write the exact Q77 equation as

    F(A,eta)=K_b A^2-C0+lambda s0(eta)
       -(lambda/2)exp(-26)R(eta)
       +lambda sum_i tau_i (v_i A+b_i p(eta))^2=0,
    p=eta(1+eta^2), C0=(1-exp(-26))/4.               (5)

Keep 0<K_b<.246, 0<tau_i<=20, |v_i|,|b_i|<=lambda^20, and the exact
pre-heat s0=s_*+s_1 p^2 with |s_*|<=500T, |s_1|<=1. These constants
and R are independent of trial A. The accepted real strength A0(eta0) is
in [.9,1.2] and has J0=F_A(A0,eta0)>.35, so J0^-1<3.
On the disk |p|<5. For |A|<=1.3 the end coefficients have modulus
<7lambda^20, and |s0|<=501T. The small terms g in (5), excluding
K_b A^2-C0, obey

    |g(A,eta)| <= (501+8192)Tlambda+1960lambda^41
                    <2^14 Tlambda.                 (6)

Here 1960lambda^40<T follows from (1) and T>=64. Since F(A0,eta0)=0,
|F(A0,eta)|<=2^15 Tlambda<1/1000.
The coefficient of A^2 is the constant
K=K_b+lambda sum tau_i v_i^2<1/4. Its linear coefficient is
B(eta)=2lambda sum tau_i v_i b_i p. The ODD cross term is retained;
|B(eta)-B(eta0)|<=560lambda^41<1024lambda^41<1/1000.

For A=A0+z solve the same equation by

    z=-[F(A0,eta)+(B(eta)-B(eta0))z+Kz^2]/J0.        (7)

On |z|<=1/10 the image has modulus at most
3[.001+.001*.1+.25*.01]=.0108<.1, and the Lipschitz constant is
at most 3[.001+2*.25*.1]=.153<1/2. Uniform complex contraction gives a
holomorphic branch on the same eta disk, with |A|<=1.3<3/2.
For real eta this branch is real and positive. The leading quadratic
coefficient K>0 and F(0,eta)<0, since (6) is less than C0>1/8.
Thus (5) has exactly one positive real root. The branch is therefore the
SAME accepted strength, regardless of whether the unknown accepted root
was initially known to lie in the local radius-1/10 ball. The local
extensions agree wherever their real intervals overlap. Cauchy now gives

    |A^(j)(eta0)|/j! <=(3/2)2^j.                    (8)

No unknown higher derivatives of A, R or the end corrections enter this
proof. Neither B nor the tiny quadratic end energy is set to zero.

## 4. Closed exterior U norm and the current first-discrepancy interface

On the same REAL radial pulse interval [Xp,Xp exp(13/lambda)],

    U=e_b f exp[-(.5+lambda)y]
              [A R0(lambda y)+sum_i (v_i A+b_i p) beta_i(y)].

The retained e_b<=1, |f|<=2 on the angular disk, real R0<=13 and
beta_i<=64 imply

    |U(eta)| <=39+1792lambda^20 <40.                 (9)

Hence every fixed factorial angular norm on this whole interval obeys

    sup_X |U0(X,.)|_m <=40 sum_(j=0)^m 2^j
                         =40(2^(m+1)-1).           (10)

The binary budget is at most m+7. This is a closed estimate for the same
EXTERIOR axial region, not a newly introduced arbitrary profile seminorm.
The sharper old third-order estimate can still be used at its original
order. The larger new all-order upper bound does not itself recertify any
old numerical selector or retained scale.

The old global radial cover is axis/core up to modulation entry, compact
modulation/repair, the reference taper with |U|_m<=8, zero axial power
stages, this pulse, and U=0 on all later post-pulse/heat stages. Thus, given
valid axis and compact-region bounds at order m, their maximum with 8 and
(10) supplies the SAME whole-support leading U norm. These two remaining
regional bounds and matching/support identities are NOT proved by this
maximum operation. global_U exposes that limitation explicitly.
current_first passes this whole-support bound at m+2 to Q125's first
incoming-discrepancy estimate, retaining its separate leading F0 order m
input and same first-coefficient Cauchy seed. No pressure or moment trace
is reset on the zero-U stages. At angular order 100005 the diagnostic
requests global U order 100007 using integer budgets only. Its supplied
axis/compact and first-source constants are manufactured, not evaluations
of the actual selected profile or a retained a_n.

## 5. Validation and remaining obligation

Exact rational checks certify each displayed sufficient value inequality
under (1); they do not certify the retained hypotheses themselves. A
manufactured scalar energy polynomial checks higher same-root jets and
the nonzero odd affine-energy term against Q77's positive recurrence.
A finite positive future-mass control checks angular jets with the same
normalizations. Tests check binary budgets, global cover assembly and the
Q125 join, with explicit identity/scope rejection. Historical root, axial,
and energy tests are replayed separately. Passing tests is not a continuum
proof. The new audit binds this derivation and replays Q125's frozen
source/evidence closure; independent output must match byte-for-byte.

Next bound the actual axis and compact leading angular norms at the needed
orders, and the required global F0 and later positive-coefficient inputs.
Then verify the existing scale, native pointwise high-pole source and
common-pole/trace contracts before the compact-modulation/exact-repair cone
join. Full-wave covariance means including axial velocity and pressure,
uniform actual physical-pulse size/derivative/support bounds and the full
viscous/cutoff/interaction residual remain open. No unconditional pulse
bound, complete residual improvement or blow-up follows.
