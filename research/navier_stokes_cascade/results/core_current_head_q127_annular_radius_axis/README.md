# Q127: same earlier annular root radius and the leading axis angular cover

Parent: `aa596da2976852c6d8a8783d3a97ff702d94c4b2` (Q126).
On 3 November 2026 the freshly fetched branch agreed with this parent.
Q126 research CI 41659390712 and field CI 41659390704 were in progress;
Q125 field CI 41649355674 and Q124 field CI 41641232330 had passed.
Current notes, frozen evidence, Q70/Q72/Q73, annular_exterior_foundation,
and the original whole-axis bound were inspected first. Historical records,
real roots, supports, frequency, pressure datum, physical scales and the
C_setup^-20 width are unchanged. The core_growth_activation conclusion
remains local leading core/first activation collar.

Q126 removed the exterior axial high-jet recurrence. Here we derive a
fixed angular radius for the EARLIER eta-dependent annular join, retain
its varying quadratic terms, and propagate the existing leading axis cap
to every fixed angular order. The annular result applies to the PREMOD
profile. The later compact modulation and repair still need their separate
real-jet estimates. No analytic radius is silently passed through a flat
variance cutoff depending on eta.

## 1. Inputs already supplied before the earlier join

At radial order ZERO take Q73's positive H=b/4 and cap Z for the actual
five-component target on D_H={dist(eta,[-1,1])<H}; H<=1. Q73 constructs
Z=500D exp(12) from the actual prefix and ideal moment comparison at
Xc=R exp(-6), BEFORE the unknown earlier correction. It includes all five
axis/core traces, the common pressure datum and full incoming memory.
There is no use of a later exact match to define that target. These
order-zero complex caps remain valid when the ANGULAR request is raised;
no radial derivative norm or complex radial cutoff is required below.
The selected parameters may make H extremely small and Z extremely large.

The earlier equation is

    B c + Q_eta(c,c)=z(eta).                         (1)

B is constant. The only varying coefficient is proportional to
(1+eta^2)^2, so Q is a polynomial of degree four. The inherited absolute
bump-integral estimate bounds every normalized coefficient derivative by
100 at every real center. This also bounds the complex multilinear maps
by the same absolute-coefficient argument. The accepted real root has

    ||c||_infinity <= d=6e-13,
    ||[B+2Q_eta(c,.)]^-1||_infinity <= beta=2000.     (2)

The latter follows from the original inverse bound 1000 and the same small
real root, not a new high-order contraction. P>=16 in the retained annulus
is part of that estimate. These are the eta-dependent EARLIER coefficients,
not the constant-Q coefficients of the later compact repair.

## 2. A complex radius retaining all coefficient variations

Fix a real eta0 and c0=c(eta0). For |w|<=s<=1/2 the polynomial expansion
gives

    ||Q_(eta0+w)|| <=100 sum_(i=0)^4 s^i <=200=:q,
    ||Q_(eta0+w)-Q_eta0|| <=200s.

Cauchy on the target disk gives, for s<=H/2,

    ||z(eta0+w)-z(eta0)|| <=2Zs/H.

Set

    a=1/(16 beta q), S=2Z/H+q d^2,
    sigma=min(H/2,1/2,a/(8 beta S)).                 (3)

S>0 even when Z=0 because d is the fixed positive upper bound. This is an
estimating radius in eta, not a physical pulse/band selector. For c=c0+u,
J0=B+2Q_eta0(c0,.) and Delta Q=Q_(eta0+w)-Q_eta0, the exact equation is

    u=J0^-1[Delta z-Delta Q(c0,c0)
                         -2Delta Q(c0,u)-Q_(eta0+w)(u,u)]. (4)

The Delta Q terms MUST remain. They include the differentiated quadratic
coefficients multiplying the undifferentiated root and its high jets.
They cannot be replaced by the constant-Q recurrence of another repair.

On ||u||<=a, |w|<=sigma, the constant part in (4) is <=beta S sigma<=a/8.
The linear perturbation has norm <=ell=2 beta q sigma d<=beta q d
=2.4e-7<1/8. The quadratic image is <=beta q a^2=a/16. Thus the image
is inside a/2 and the Lipschitz constant is

    ell+2 beta q a <1/4.

Uniform complex contraction proves a holomorphic branch on the disk.
For real eta in that disk, the ALREADY accepted root differs from c0 by
at most 2d<a. It solves the same equation in the same contraction ball,
so it is precisely this branch. Overlapping real intervals identify the
local extensions with the original root. There is no new root, no change
to N and no smallness requirement in every higher angular norm. Cauchy
therefore gives

    ||partial_eta^k c(eta0)||/k! <=a sigma^-k, k>=1, (5)

while its real zeroth-order bound remains d. The looser complex value
cap d+a must not replace the old real smallness bound in acceptance tests.
The old sharper C6 estimates also remain valid at their original orders.

## 3. PREMOD fields, inverse swirl and pressure scope

Use Q73's actual continuation cap ua at radial order zero. On [Xa,Xi]
U has modulus <=ua on the original larger strip. On the shape U=Gi;
on the fixed restoration U=4eta+g(1-sigma_radial), |g|<=ua+8, so
|U|<=ua+16. Here sigma_radial is evaluated only at its REAL radial
argument; it is distinct from the angular radius sigma in (3).
On the earlier join the same native field is

    U=4eta+sum_(two U bumps) c_l(eta) beta_l(y).

The actual width-.08 bumps have value <=64 and no unit-mass prefactor.
Thus on the complex disk |U|<=8+128(d+a)<16. The ideal exit after the
bumps has |U|<=8. Therefore on the WHOLE PREMOD interval [Xa,R],

    sup_X |U_PREMOD|_m <=(ua+16) sum_(j=0)^m sigma^-j. (6)

The preceding intervals have analytic eta caps on larger strips, so using
the smaller uniform sigma is legitimate. This proves a closed geometric
alternative to Q73's high-angular root recurrence, relative to its explicit
order-zero prefix constants. It is not a radial analyticity claim.

There is also a useful nonzero-swirl consequence on -6<=y=log(X/R)<=0.
There E=P f exp(y/10)(1+sum_(three E bumps)c_l beta_l), f=(1+eta^2)^-1.
On the disk |f|<=2 and |f^-1|<=13/4. The correction is <=192(d+a)<1/2;
exp(y/10)>=1-3/5=2/5. Consequently

    |E|<=4P, |E^-1|<=20/P                            (7)

on this earlier annular patch and ideal exit. This does not extrapolate
(7) back through the entire shape interval, where Q73 has its own caps.
Nor does it set any local pressure to zero. The target retains all five
prefix moments and the exact endpoint match remains a functional identity;
partial pressure inside the bumps requires its genuine integral and both
pressure terms in the forward state map, as in Q67/Q70.

## 4. Deferred budgets and the same leading axis cover

If H^-1<=2^b_H and Z<=2^b_Z with nonnegative integer budgets, then
S<=2^(2+b_Z+b_H) and 8 beta/a=102400000000<2^37. Hence

    sigma^-1 <=2^(39+b_H+b_Z).                       (8)

A binary bound for (6), given ua<=2^b_u, is

    1+max(b_u,4)+ceil(log2(m+1))+m(39+b_H+b_Z).       (9)

No factorial, high power, target polynomial or high-order root solve is
evaluated by this budget. Actual H,Z,ua must be bound by the same Q73
order-zero construction; the example inputs in tests are manufactured.

For the LEADING axis region up to the original flat modulation entry,
global_axial_angular_bound (4) already supplies its original pointwise
complex cap K=K(ell_axis), with the same delta=rho_axis/(800S) Cauchy
radius, including X=0. Its identical proof at every fixed order gives

    K_axis,m <=K sum_(j=0)^m delta^-j.               (10)

This is an extension of its displayed third-order angular estimate, not
an extension of the radial cutoff/collar or a change to delta. Binary
budgets cost b_K+ceil(log2(m+1))+m ceil(log2(1/delta)). The helper requires
this LEADING axis identity; it is not Q125's first-coefficient seed bound.

Combining (10) with Q126 now bounds the whole leading U by the maximum of
the axis norm, a separately proved FINAL compact modulation/repair norm,
8 for the reference taper, and 40(2^(m+1)-1) for the exterior axial pulse;
later axial fields vanish. The API explicitly rejects its PREMOD output
as that final compact input. The final loop contains a flat variance cutoff
composed with eta-dependent data; analytic PREMOD inputs alone do not
establish a complex extension after that operation. Its finite real
composition estimates Q65/Q66/Q64 and same-root repair remain necessary.
This does not assert failure or nonanalyticity of the actual final field.

## 5. Checks and next obligation

Exact rational checks verify every value-ball inequality, inverse-radius
budget and moderate geometric enclosure. Scalar varying-quadratic controls
check nonzero higher root jets and the missing-Q-derivative countercontrol,
including Q_eta derivatives acting on c0^2. Same-cover axis and global
assembly checks retain the stage and order conventions. These are algebra
and provenance checks, not numerical selected-profile evaluation or a
continuum proof by testing. The audit replays Q126 and its frozen source
closure; fresh output must match byte-for-byte.

Next propagate these actual PREMOD field/state caps through the real
variance-cutoff, signed-root, moving-phase and fixed-frequency modulation
maps, then the SAME later compact five-row repair with complete pressure
memory. This is the remaining final compact U input; global F0 and later
positive-coefficient norms and all original scale/source gates also remain.
The native order-100005 pointwise source/common-pole and full cone join,
complete-wave covariance means including axial velocity and pressure,
uniform actual physical pulses, and the full viscous/cutoff/interaction
residual remain unproved. No unconditional pulse bound, residual improvement
or blow-up follows.
