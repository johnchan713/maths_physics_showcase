# Current scalar target, projected target and the higher-jet pressure join

Parent: `a421ceab5aa765b6c08f8025bf3aa2ef331a3527`.
The preceding checkpoint supplies finite profile-dependent bounds at each
fixed derivative order. Here we join those bounds to the CURRENT angular
scalar defect Jtheta used by the axial slow inverse, distinguishing it from
the projected signed target Xi_theta used by the signed wave update.

The exact current identities and the conditional cycle already retain both
objects. This checkpoint does not replace that cycle or independently prove
its incoming field hypotheses. Its new estimate resolves which current data
control the higher scalar-target jets in the pressure bounds. Historical
evidence, profiles, fixed tilt, width and accepted domain are preserved.

## 1. The scalar channel survives projection

Use the current full equations, both zero velocity moments, and complete
Haar-averaged covariance. With D2=partial_R+2/R, define

    J=Jtheta=integral R^2 [G v0+V g0+bar(gamma v)+bar(Wtheta,z)] dR,
    M=integral R^2 bar(Etheta) dR=epsilon partial_Z J,
    Q2 H=H-sigma_theta integral R^2 H, integral R^2 sigma_theta=1,
    Xi=-I2 Q2 bar(Etheta), I2 H=R^-2 integral_left^R r^2 H dr.

The two exact channels reconstruct the current angular equation:

    bar(Etheta)=-D2 Xi+sigma_theta M.                       (1)

The scalar part has been removed by Q2. Therefore Xi alone cannot control
M: for any smooth scalar m, the compact source E=sigma_theta m has Xi=0
and weighted integral m. This is an operator countercontrol, not a claim
that these sources are actual canonical pulse states. An estimate for the
projected target cannot supply the discarded scalar channel without another
equation or induction hypothesis.

The retained current cycle explicitly includes both projected target bounds
and scalar defects in its input bundle. It is consistent with (1). In
particular, substituting its improved Xi order for the scalar J order would
be an extra unjustified inference. The pressure slow inverse uses J, not Xi.

## 2. Current-field control of the zero-Z scalar jets

Let the fixed positive integration interval I=[a,d] contain the supports of
ALL terms in J, including the whole current wave covariance. A wider wave
support must not be assigned to the reserved slow potential hull. Write
L=d-a, M2=L max(a^2,d^2). Let |.|_n be the normalized mixed Z,T norm of the
previous checkpoint and let

    <f>_n=sum_(j<=n) sup |partial_T^j f|/j!,
    A_n=<J>_n.

The supremum includes the full slow domain and the actual radial/auxiliary
variables when present. Differentiation commutes with integration of smooth
zero extensions. Product differentiation and Haar contraction give

    A_n <= F_n := M2 [<G>_n <v0>_n + <V>_n <g0>_n
                       +<gamma>_n <v>_n + <Wtheta,z>_n].   (2)

The base factors use unweighted norms. This is a sufficient raw bound; it
may lose cancellations. Every product is averaged after multiplication.
In particular zero individual Haar means do not remove bar(gamma v).
No radial or Z derivative is needed to estimate the pure-T channel (2).
If a stronger scalar order relies on cancellation between the displayed
terms, its actual current value or exact update difference must be bounded.

## 3. Recover all remaining scalar jets from the actual integrated residual

Epsilon is constant within the Q-chart. For i>=1,

    partial_Z^i partial_T^j J
       =epsilon^-1 partial_Z^(i-1) partial_T^j M.           (3)

Split the normalized mixed norm into i=0 and i>=1. Since 1/i<=1,

    |J|_n <= A_n+epsilon^-1 |M|_(n-1)
           <= F_n+epsilon^-1 |M|_(n-1), n>=1.             (4)

This is a local identity and norm estimate; it introduces no slow-domain
rectangle, reference trace, new boundary condition or first-band choice.
For n=0 use |J|_0=A_0<=F_0. A trace-only reconstruction would require an
additional integration-domain argument; it is unnecessary for (4).

The epsilon^-1 cost matters. It is not absorbed into a fixed profile
constant. For the usual consistent scalar input J in S_C, one has
A_n=O(epsilon^C) and M in S_(C+1) at every required fixed order. Then (4)
recovers J in S_C. A projected Xi bound, even Xi=0, supplies neither of
those two scalar inputs by itself.

## 4. Quantitative join to both pressure mismatch terms

Let d_m be the proved fixed-profile bound for D_rho=I0(rho_*-rho), and
k_n=|K|_n for the fixed slow inverse coefficient F=KJ. Put
q_j(m)=(m+1)...(m+j). The previous checkpoint's axial-viscous mismatch
pressure gradient now has the sufficient CURRENT-data bound

    |epsilon^5 partial_Z[D_rho partial_Z^3(KJ)]|_m
      <=q_4(m) d_(m+1) k_(m+4)
         [epsilon^5 F_(m+4)+epsilon^4 |M|_(m+3)].          (5)

The slow-time mismatch term is retained as well:

    |epsilon^3 partial_Z[D_rho partial_T partial_Z(KJ)]|_m
      <=q_3(m) d_(m+1) k_(m+3)
         [epsilon^3 F_(m+3)+epsilon^2 |M|_(m+2)].          (6)

For sharper estimates A can replace F. All six moving-density commutators
remain bounded by the preceding checkpoint's formula after substituting
f_n=|KJ|_n<=k_n(F_n+epsilon^-1|M|_(n-1)). No commutator, full pressure
memory or persistent radial coefficient is removed by this substitution.

Thus a concrete sufficient next inequality for proposed tolerances T_m is
that the right sides of (5),(6), together with the retained commutator and
other residual budgets, fit T_m on the same accepted bands. Equations
(2)--(6) supply the missing scalar-input routing. They do not certify that
the actual field norms or scalar cancellation bounds meet those tolerances.

## 5. Fresh fast data and order checking

For a zero-Haar fast increment (db,dv,dg), the base-linear change in J
vanishes, but exactly

    delta J=integral R^2 bar(gamma dv+dg v+dg dv) dR.        (7)

Old/new correlations survive. For old tangential order mu and fast order H,
the raw scalar change is only O(epsilon^min(mu+H,2H)). Its integrated residual
change epsilon partial_Z delta J has one extra epsilon under the retained
fixed-jet hypotheses. This distinction agrees with the earlier conditional
cycle; an extra epsilon must not be assigned to delta J itself.

A smooth compact control makes (7) nonzero even when both old and fast
Haar means vanish. Choose compact psi, f=(partial_R+1/R)psi and a constant
carrier c. With the commuting graph derivatives D_R=partial_R+c partial_y
and D_Z=epsilon partial_Z, use the complete curl of psi a A cos(y):

    beta=-epsilon psi partial_Z(a A)cos(y),
    gamma=a A[f cos(y)-c psi sin(y)], v=gamma/R.

Use the same construction with amplitude d B for the increment. Both
weighted HAAR velocity moments vanish, and divergence is exactly zero with
the graph derivatives. Holding the complete wave fixed, the scalar change is

    delta J=[a d A B+(d^2/2)B^2] integral R(f^2+c^2 psi^2)dR. (8)

The integral is positive for nonzero psi. Slow amplitudes can give nonzero
fourth Z jets as well. This is an admissible moment/divergence control,
not the prescribed fast inverse or an actual canonical cycle. It shows why
the slow solve must read the recomputed postfast scalar data.

Next prove the actual pure-T scalar bounds/cancellations and integrated
residual jets in (5),(6), alongside the projected target recurrence and the
full current field estimates. Join all generated-source supports, full
covariance, viscosity, cutoffs and interactions before claiming complete
residual improvement. No unconditional pulse or blow-up theorem is proved.

## 6. Validation

Exact rational current fields independently check the complete angular
equation, scalar moment, compact projection and reconstruction (1). Separate
finite-rank controls show its scalar blindness. Coefficient majorants check
(2),(4)--(6); direct derivatives verify (3) and the product jet routing.
Zero-Haar curl controls verify (7),(8), both moments, divergence and surviving
fourth jets. Fixtures have finite endpoint regularity and are not the actual
smooth global construction. Continuum conclusions follow from the displayed
identities and estimates, not from test counts. Historical evidence is intact.

```sh
python3 research/navier_stokes_cascade/results/core_scalar_target_pressure_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_scalar_target_pressure_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_scalar_target_pressure_join/evidence.json \
  --output /tmp/core_scalar_target_pressure_join_fresh.json
```
