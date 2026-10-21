# Fixed-profile tame bounds and correction of the compensator data gate

Parent: `77a8b462019590178972b029ece632eac80e2313`.

The previous checkpoint correctly proves that mass, support and second
moment do not give profile-independent pressure-shape bounds. Its conclusion
that a numerical profile or numerical seminorms must be supplied before
continuing an asymptotic argument was too strong. The construction fixes one
smooth compact physical profile, and its estimates may depend on that profile.
Every fixed-order seminorm is then finite. An explicit profile is required
to evaluate numerical constants or certify an already specified finite band;
it is not required for existence of profile-dependent tame constants.

This checkpoint supplies those constants at every fixed radial and slow
derivative order and joins the entire moving-density pressure commutator.
It does not choose a new profile or change the accepted tail. Historical
notes and frozen evidence remain unchanged; this note corrects the
interpretation of the preceding gate.

## 1. Quantifiers and the actual remaining obligation

For the SAME fixed smooth compact rhohat and reference rhohat_*, set

    hhat=rhohat_*-rhohat,
    Dhat(x)=integral_left^x hhat(y)dy, E=x partial_x.

Both densities have unit mass. Thus Dhat is smoothly compactly supported
in their convex hull [l,u], 0<l<u, inside the reserved interior. In particular

    B0=||Dhat||_infinity <= (1+||rhohat||_L1)/2 < infinity,
    ||E^k Dhat||_infinity < infinity for every fixed k.       (1)

The right-hand side depends on the selected profile; signed unit mass is
not an L1 bound. The previous counterfamilies remain valid, but vary the
profile. They do not contradict

    for each fixed profile, for each finite order,
    there is a finite constant uniform in the physical band. (2)

No uniformity over arbitrary profiles, derivative orders, tilt tending to
zero, or cycle length is asserted. The genuine unresolved input is still
the ACTUAL linked target/field bounds at the necessary higher orders and
their complete residual cancellations. Merely declaring target smoothness
does not provide uniform band or stage estimates.

## 2. A finite-jet composition majorant without analytic profiles

Use the accepted chart 1/4<=s<=4 and the inherited coordinate contraction.
At any real center s0 write s=s0(1+t). On the normalized complex bidisk
of radius 1/16 the coordinate satisfies |t|<=1/4. Therefore

    u=log(s/s0)=log(1+t), |u|<=sum_(k>=1)(1/4)^k/k<=1/3.

Cauchy's formula and the same conversions to physical Z,T increments give

    |partial_Z^i partial_T^j u|/(i!j!)
       <=(1/3)32^i 64^j, i+j>=1; u(center)=0.             (3)

Only the coordinate is extended to complex arguments. For a smooth real
profile f, b=0 or 1, and nonnegative radial order r, consider

    V_b=s^(-b/2) f(R/sqrt(s)), c=b+r,
    M_(b,r,k)=||(E+c)^k partial_x^r f||_infinity.

With x0=R/sqrt(s0), radial differentiation and the real chain rule yield
the finite jet identity through any total slow order m:

    partial_R^r V_b
       = s0^(-c/2) sum_(k=0)^m (-u/2)^k/k!
                    (E+c)^k partial_x^r f(x0) + O(jet order m+1). (4)

This is a statement about derivatives at the center, not convergence of an
infinite series for f. Ordinary finite Taylor's theorem for
v -> exp(-cv/2) partial_x^r f(x0 exp(-v/2)) proves it for C^m profiles.
It also applies at physical endpoints using their smooth one-sided limits.

Introduce a formal series with nonnegative rational coefficients

    A(z,t)=(1/3)[1/((1-32z)(1-64t))-1],
    a_(m,k)=sum_(i+j<=m) [z^i t^j] A(z,t)^k,
    G_(b,r,m)=2^(b+r) sum_(k=0)^m M_(b,r,k) a_(m,k)/(2^k k!). (5)

Here A^0=1, a_(m,0)=1. Coefficientwise multiplication of the bounds (3),
followed by (4), proves for the normalized unweighted slow norm

    |partial_R^r V_b|_m
      :=sum_(i+j<=m) sup |partial_Z^i partial_T^j partial_R^r V_b|/(i!j!)
      <= G_(b,r,m).                                      (6)

Every quantity on the right is finite and independent of the band index.
Use b=0,f=Dhat for the mismatch primitive, and b=1,f=rhohat or rhohat_*
for either complete moving density. The generator E+b+r and the prefactor
2^(b+r) are essential for density and mixed radial derivatives.

All these fields and their derivatives have the actual moving hull
[sqrt(s)l,sqrt(s)u]. If that whole hull is strictly inside the same active
annulus, the inherited weight comparison multiplies (6) by exp(E_geometry).
The union interval [l/2,2u] is used for integration, not for an asserted
Cartesian weight floor. Supports reaching the flat annular edges need the
separate edge estimates already recorded in the repository.

## 3. Join to the mismatch and every moving-density commutator

Write F=KJ and d_m=G_(0,0,m)(Dhat), r_m=G_(1,0,m)(rhohat_*).
Let f_m=|F|_m and q_j(m)=(m+1)...(m+j), q_0=1. The exact parent source is

    g_v-rho P_v = D_rho,R [epsilon^2 F_TZ+epsilon^4 F_ZZZ]+C,
    C=epsilon^2[rho_*,Z F_T+rho_*,T F_Z+rho_*,TZ F]
      +epsilon^4[3rho_*,Z F_ZZ+3rho_*,ZZ F_Z+rho_*,ZZZ F]. (7)

The axial-viscous mismatch pressure contribution satisfies

    |epsilon^5 partial_Z[D_rho partial_Z^3(KJ)]|_m
      <=epsilon^5 q_4(m) d_(m+1) |K|_(m+4)|J|_(m+4).     (8)

The slow-time mismatch pressure gradient is also retained:

    |epsilon^3 partial_Z[D_rho partial_T partial_Z(KJ)]|_m
      <=epsilon^3 q_3(m) d_(m+1) |K|_(m+3)|J|_(m+3).     (9)

Leibniz, submultiplicativity and derivative shifts give ALL six terms in C
the finite bound

    |C|_m <= c_m :=
      epsilon^2[2(m+1)^2 r_(m+1) f_(m+1)+q_2(m) r_(m+2) f_m]
      +epsilon^4[3(m+1)q_2(m) r_(m+1) f_(m+2)
                 +3(m+1)q_2(m) r_(m+2) f_(m+1)
                 +q_3(m) r_(m+3) f_m].                  (10)

Each radial integral in C vanishes, including after slow differentiation,
because every displayed derivative of the unit-mass rho_* has zero total.
The ordinary primitive is therefore compact. On the common integration
hull of length L_hull=2u-l/2, the half-L1 estimate implies

    |epsilon partial_Z I0 C|_m
       <=epsilon(m+1)(L_hull/2)c_(m+1).                  (11)

Thus fixed-profile smoothness supplies uniform coefficients also for the
entire moving-density commutator. It does not remove that commutator or
the persistent radial principal coefficient -epsilon^4 K rho. That radial
coefficient has a finite upper norm by (6), and nonzero radial L1 size by
the parent identity. These are statements about the isolated repair terms.
The complete covariance, axial velocity, Pi1/Pi2 memory and all other
transport, viscosity, cutoff and nonlinear terms remain required.

## 4. Strict exponent margins versus numerical finite-band certification

If actual linked targets satisfy |J|_(m+4)<=A_m epsilon^alpha and
|K|_(m+4)<=K_m uniformly on the accepted bands, (8) is bounded by

    C_m epsilon^(alpha+5),
    C_m=q_4(m) d_(m+1) K_m A_m < infinity.                (12)

No numerical value for C_m is needed for this asymptotic order statement.
For a smaller requested gain delta<5 and tolerance tau>0, its sufficient
finite-band inequality is

    C_m epsilon^(5-delta)<=tau.                          (13)

For finitely many orders and conditions with positive exponent margins,
there exists one positive epsilon threshold satisfying them all. This is
an existence observation; no new tail is selected or certified here.
An already fixed first band requires checking its actual epsilon and
constants. Zero exponent margin instead requires C_m<=tau, and small
epsilon alone cannot supply it. There need not be one threshold for
infinitely many derivative orders or an infinite correction cycle.

The analogous statement for (9) has base gain 3, and (11) retains the
orders and input losses displayed in (10). If the actual high target jets
grow with frequency or band, their losses must first enter A_m and the
exponent; a nominal epsilon power does not absorb unspecified target jets.
The principal same-order derivative-loss obstruction remains intact.

Next establish those actual higher target and current-field bounds with
the linked cancellation mechanism, then join the full residual on the
existing common domain. No unconditional pulse bound, residual improvement,
infinite corrected field or blow-up theorem is claimed.

## 5. Evidence and reproduction

Rational formal bivariate powers independently check the log-scale jet and
the Euler-generator composition against direct power substitution in the
implicit chart. Mixed radial orders and density weights detect omission
of the generator shift. Exact polynomial product identities check every
term of (7), its zero radial mass and the compact commutator pressure.
Strict/zero-margin controls distinguish fixed constants from missing
band-dependent bounds. Historical records are hash-bound and unchanged.
Finite checks do not prove the continuum estimates; sections 1--4 do.

```sh
python3 research/navier_stokes_cascade/results/core_fixed_profile_tame/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_fixed_profile_tame/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_fixed_profile_tame/evidence.json \
  --output /tmp/core_fixed_profile_tame_fresh.json
```
