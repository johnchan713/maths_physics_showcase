# Regional weighted-value gate from finite first-coefficient seminorms

Parent: `0c7ff52961d301880d0078d462ea8fe402328855`.
Status: `finite-seminorm-weighted-input-and-regional-background-gap-join`.
On 18 October 2026 the parent's new background-gap CI step passed; the full
workflow remained running. Historical source and frozen evidence are unchanged.

The parent gives an explicit sufficient error radius r_* but retains a
lumped bound Cval*q^(2h). This checkpoint DERIVES that constant from the
actual shared Borel-budget specification and finite first-coefficient
seminorms. It also joins the already proved sign-stable outer interval,
so a large harmless outer stress correction is not mistaken for cone failure.
The new result is this quantitative input/regional JOIN, not a new derivation
of the heat-tail sign theorem or another canonical initial-seed table.

The proof is relative to the reviewed actual leading profile, common
coefficient supports, exact pressure/moment equations, smooth W1 edge
factorization and admissible common cutoff schedule. It does not independently
prove every primitive core/regional small-root estimate, evaluate huge
coefficient constants or certify the historical first band. No unconditional
pulse, complete residual improvement, single-field/force or blow-up follows.

## 1. Use ONE schedule for velocity AND the prescribed stress

Retain h>0, h<1/100, the frozen N, positive lambda, width C_setup^-20,
all roots/carriers/supports and pressure gauge. At fixed normalized q,X,eta,
the actual tangential increments and prescribed stress are

    e=sum_(n>=1) chi(c_n q) q^(2nh) E_n,
    u=sum_(n>=1) chi(c_n q) q^(2nh) U_n,
    That=T0+sum_(n>=1) chi(c_n q) q^(2nh) T_n.                  (1)

The common chi satisfies 0<=chi<=1, is one through 1/2 and zero from 1.
Use the SAME c_n in these fields, not independently selected stress and
velocity schedules. The stress is the prescribed Borel-summed stress in
the augmented-background equation, NOT an assertion of a full wave
covariance identity or cancellation of the bare momentum residual.

Write W_n=T_n/|T0| with the exact smooth factored extension at stress edges.
For every n>=2 the existing combined cutoff prescription includes a finite
C_n bounding BOTH (E_n,partial_X E_n,partial_X U_n) in scalar sup norm
and W_n in Euclidean sup norm. Its stronger all-jet requirements also stay.
If a recorded vector seminorm uses componentwise maxima, the existing
factor 2 in its C_n,0 covers the sqrt(2) norm conversion; it is not omitted.
The retained full-support condition is

    C_n q^(nh)<=2^-n whenever q<=c_n^-1.                      (2)

It follows from taking the maximum of the velocity and relative-stress
seminorm requirements at the same stage; see `background_jet_foundation`,
(10)--(11), and `weighted_background_edges`, (8)--(9). Each coefficient
seminorm is finite by the common-collar construction and exact support/edge
proofs. This uses the actual admissible definition, not numerical evaluation
of its c_n. For a dyadic subbudget, if p_n>=log2(max(1,Cvel_n,Cstress_n)),
then any integer b_n with

    b_n>=b_(n-1)+1,  b_n*n*h>=p_n+n                         (3)

gives (2) for c_n=2^b_n. Add strict padding and EVERY earlier physical,
potential and higher-jet requirement to (3). This is not a proposal to
replace the actual cutoff scales by the diagnostic dyadic examples.

If q>c_n^-1 the term is zero. On its WHOLE support (2) gives a bound
2^-n q^(nh). The cutoff is independent of X,eta, so the same statement
holds for the displayed radial profile jets. With z=q^h in (0,1],

    sum_(n>=2) 2^-n q^(nh)
       = z^2/[4(1-z/2)] <= z^2/2 = q^(2h)/2.                (4)

The actual sums are locally finite for q>0; the infinite positive majorant
is only a bound. Equation (4) is a whole-small-q estimate, not an endpoint
sample or a convergent uncut coefficient-series claim.

No first-cutoff-one assumption is needed: |chi(c1 q)|<=1 suffices.
This value/first-radial-jet bound therefore avoids a hidden finite-block
constant on the rest of the q interval. Higher jets retain the previous
finite-block/tail proof and order-dependent constants; they are not all
assigned the same first-coefficient bound.

## 2. The explicit finite-seminorm constant

On a closed profile annulus K=[Xlo,Xhi], Xlo>0, put

    Lr=sqrt(2Xlo), Hr=sqrt(2Xhi),
    A_E=sup_K |E1|, A_EX=sup_K |partial_X E1|,
    A_UX=sup_K |partial_X U1|, A_W=sup_K |W1|.

All suprema also include eta in [-1,1], using the actual smooth endpoint
extensions. These are finite functionals of the finalized FIRST coefficient
and stress, not diagnostic replacement values. Define

    B_E=A_E+1/2, B_EX=A_EX+1/2,
    B_UX=A_UX+1/2, B_W=A_W+1/2.                              (5)

Equations (1)--(4) give bounds B_E q^(2h), B_EX q^(2h),
B_UX q^(2h) for e,e_X,u_X, respectively, and

    |(That-T0)/|T0|| <= B_W q^(2h).                         (6)

The relative stress bound in (6) is INDEPENDENT input to the shared budget,
not a consequence of the local velocity seminorm. Full pressure sources
and moments determine T_n; in particular vanishing positive velocity on I4
does not erase Pi1=K1-mu^2/(4X) or Pi2=K2. Neither pressure memory is reset.
A zero-velocity/arbitrary-stress sequence is an abstract budget countercontrol,
not asserted to be the actual PDE coefficient family. It shows why (2)
cannot be justified from velocity-only cutoff requirements.

In the actual physical band let s=q/Q in [1/4,4], R=sqrt(s)*sqrt(2X)
and A=1/2+h. At fixed Z,T, s and eta are independent of R: q depends
only on physical z,t. Thus the EXACT pullbacks are

    Delta F=s^(-A)*e/R,
    Delta g_theta=R partial_R Delta F
                 =s^(-A)*[(R/s)e_X-e/R],
    Delta g_z=partial_R(s^(-A)*u)=s^(-A)*(R/s)u_X.            (7)

The -e/R term must remain. Equivalently their common prefactor in profile
radius sqrt(2X) is s^(-1-h). For h<1/100<1/8,
s^(-1-h)<=4^(1+h)<5: the last inequality follows from 5^8>4^9.
Consequently a fully explicit sufficient constant for the parent's input is

    C_K=max(1, 5 B_E/Lr,
                   5[Hr(B_EX+B_UX)+B_E/Lr], B_W),
    max(|Delta F|, |Delta g|, |(That-T0)/|T0||)
                         <= C_K q^(2h) on K.                 (8)

We use the component sum to bound the Euclidean norm of Delta g; no square
root is underestimated. The stress quotient loses no chart factor: the
same positive band normalization multiplies numerator and leading stress.

Equation (8) discharges Cval as a SEPARATE black-box bound given the
specified coefficient/schedule inputs. It does NOT evaluate A_E,A_EX,A_UX,
A_W for the huge actual profile. The nonzero radial base flow is not set
to zero: its Stokes-potential cutoff commutator, including the extra cutoff
derivative, remains in the full velocity specification and regular operator.
It is absent from (7) only because the parent's value-gap variables are
F and the tangential shear g, not the full momentum operator.

## 3. Use the compact gate and the already sign-stable outer region

Take Xc at the beginning of the reviewed terminal heat interval (y=1/2
in `outer_edge_stability`), with the actual common positive-order support
ending to its left. Cover the active annulus by

    Kc=[Xa,Xc],  Ko=[Xc,Xb].

Use (8) with K=Kc and the parent's explicit radius r_Kc, evaluated using
valid leading gap bounds on Kc (the old global bounds are also sufficient).
The compact value gate is now

    C_Kc q^(2h)<=r_Kc.                                      (9)

On Ko, exact common supports and forward/backward moment identification
give the already reviewed identities

    actual velocity = leading heat-tail velocity,
    That=T0+tau(q) T1, tau(q)=q^(2h)chi(c1q)>=0,
    T0theta>0, T1theta>=0, T1z=0, N=(-1,0), K=(0,-1).        (10)

The current shear growth is therefore unchanged. For the SAME fixed u,
adding tau T1 increases positive angular stress and leaves axial stress
unchanged. Its oriented unit angular component does not decrease and

    ((1+u^2)c^2 That_z^2)/(u^2 That_theta^2)
             <= ((1+u^2)c^2 T0z^2)/(u^2 T0theta^2).           (11)

Thus the already accepted leading same-tilt cone is preserved on Ko
WITHOUT requiring tau |W1| to be small. This does not rederive the heat
sign theorem; it uses its actual sign/support identities in the new
compact-plus-outer quantitative cover. A signed cutoff would invalidate
the argument. Common closed-edge directions and fixed jets remain those
of the original factorizations.

The distinction can matter greatly. The EXACT outer edge value is

    W1(Xb,+/-1)=(2Xb/(1-2h)^2,0).                           (12)

Hence any global A_W is at least 2Xb/(1-2h)^2. When the first cutoff is
one, a global small-relative-error certificate would require
q^(2h)*2Xb/(1-2h)^2<=r_* just at this edge. That is necessary for THAT
small-error test, NOT for actual cone preservation: (10)--(11) can hold
when the test fails. The parent's global test was already stated only
as sufficient. Here we avoid imposing its potentially severe harmless
outer-edge cost. Compact seminorms can STILL depend on Xb through actual
incoming pressure/moment memory; no Xb-uniform bound is claimed.

This refinement concerns the growth/orientation/CONE interface only. If
another lemma needs the global magnitude upper bound |That|<=3|T0|/2,
it must still retain its own global small-error gate. Formula (8) also
applies on the whole annulus and supplies C_full for such a requirement;
the compact gate alone does NOT prove C_full*q^(2h)<=1/2. No enlarged
domain for every pulse, root, derivative or residual map is certified here.

## 4. Domain, historical-instance and remaining analytical boundary

Given the reviewed inputs, (9)--(11) preserve background growth and the
same-tilt cone on the full closed active shell. A sufficient full-chart
band requirement is

    ell >= 2+log2(C_Kc/r_Kc)/(2h), Q=2^-ell, q<=4Q,           (13)

rounded upward and combined with EVERY previous geometry, normal, Riccati,
root, cone and corrected-freezing gate. The leading quantities, actual
profile, h,N, width, roots, pressure gauge, Borel specification and frozen
pulse frames are unchanged. This is a sufficient later-tail selection;
it does not certify a historically fixed first band or transplant the
constant-reference shear frame into a moving physical pulse.

The numerical first-coefficient seminorms and primitive regional accepted-root
inequalities have not been evaluated or independently closed here. Passing
the tests does not supply those continuum inputs. A later distinct realized
field/force still needs its own controlled schedules and complete residual
join; historical fields and forces are not silently recertified.
Pi1/Pi2, actual A0/A1 velocities, viscosity, axial pressure/velocity and
both interaction orientations remain in the full residual. Spatial/diagonal
cutoff derivatives are not reclassified as temporal Gaussian tails.

Next obtain justified quantitative first-coefficient and compact leading
gap bounds from the actual retained core/regional defining maps, then check
their combined later-tail gate. Do not repeat the heat-sign or canonical
seed proofs, or claim that this regional value estimate is an unconditional
pulse, complete residual, single-field/force or blow-up theorem.

## 5. Reproduction and interpretation

Exact rational controls verify the full geometric tail, combined dyadic
subbudgets, nonunit chart factors and the inverse-radius derivative in
(7). Outer algebra controls keep the same u, including very large positive
corrections, and detect sign failure for a negative cutoff. They are
manufactured scalar/profile/sequence controls, not actual coefficient
norm evaluations, heat simulations or continuum supremum certificates.
An additional composed control uses compatible local radial jets, the new
C_K and the parent's explicit radius, then independently checks its spectral
and same-tilt conclusions with exact squared norms. It is a local fixture,
not evidence that the giant profile meets those bounds on a whole annulus.
Historical source bindings are checked recursively; fresh relevant test
and evidence replays are separate commands with original comparators.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json`.
Use new output paths; preserve all historical evidence.
