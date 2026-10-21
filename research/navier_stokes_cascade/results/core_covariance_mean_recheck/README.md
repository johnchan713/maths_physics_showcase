# Complete covariance mean defects and quadratic conditioning

Parent: `c9e9fc845bf3b762a3f56766c3b8c885181d8c99`. This checkpoint continues
`core_primary_tail_recheck`; it does not replace any frozen evidence. The
initial mean identities already derived in `complete_covariance_initial_means`
are retained. The new analytical result is a pair of projected energy lower
bounds for the **recomputed nonlinear defects**, with an explicit axial
conditioning cost of order lambda^{-2} on a fixed reserved interval. This is
a necessary compatibility budget, not a failure theorem for the actual pulse.

## 1. State, supports and complete pressure

Fix a slow point (Z,T) in one Q-chart and a reserved hull I=[r_-,r_+],
0<r_-<r_+<infinity. All slow increments are smooth and compactly supported
inside this hull; all moment integrals below are over (0,infinity). When the
hull moves with Z,T, differentiate the compactly extended integrands, not a
truncated pressure integral. At the actual reserved patch

    G=0,  V(R)=a_eff R^(-1-2 lambda),  a_eff>0,  lambda>0,

where a_eff=a(eta)s^(lambda-h) is independent of R at the fixed slow point.
The radial background b=O(epsilon) is retained. Keep the complete real wave
covariance W fixed during this initial slow update. In particular, W contains
all six entries, axial velocity, principal/curl cross terms and curl squares.
The incoming defects P_in,Jtheta_in,Jz_in are computed from that full current
state after the initial fast correction. Its three Haar means must be exactly
zero. Only then do its products with the auxiliary-independent slow increments
(B,u,g) have zero Haar average. No pointwise cancellation is asserted.

The inherited five-map supplies

    D1 B + epsilon g_Z = 0,
    integral R^2 u = integral R g = 0,
    integral 2 V u/R = -P_in,
    integral R^2 V g = -Jtheta_in,
    integral -R V u = -Jz_in,               D1=partial_R+1/R.

These rows cancel linear changes only. The source pressure includes both
Pi1=K1(eta)-mu^2/(4X) and Pi2=K2(eta), with their distinct physical powers
and complete cutoff derivatives. Their actual axial gradients remain in the
incoming state and linearized balance. They are not deleted by the fact that
one plateau contribution has zero axial derivative. Pressure is recomputed
from the entire radial source, including its compact primitive and moment
compensation, after the update. At G=0 the complete defect is
Jz=integral R overline(gamma^2+W_zz) - (integral R^2 overline(g_r))/2.
The actual axial residual moment is epsilon partial_Z(Jz+c_rho P), where
c_rho=(integral R^2 rho)/2 for the unit-mass physical compensator rho.
The WHOLE moving product is differentiated, including (partial_Z c_rho)P.
Thus pressure is not an independent adjustable cancellation below; the
defect budgets also do not by themselves bound this differentiated residual.

## 2. Exact dynamic remainder

Write S=2bB+B^2, H=bg+GB+Bg (G=0 here), tau=epsilon partial_T
+epsilon^3 partial_ZZ. Subtract the incoming state and the linear V u term
from the Haar-averaged radial source. The exact remaining source is

    Rg = epsilon B_T - (partial_R+1/R)S - epsilon H_Z + u^2/R
         +epsilon (B_RR+B_R/R+epsilon^2 B_ZZ-B/R^2).

This includes radial vector viscosity, axial viscosity, radial transport and
axial transport. Define

    E_u=integral u^2/R, F_u=integral R u^2, E_g=integral R g^2,
    C_P=tau integral B - integral S/R - epsilon partial_Z integral H,
    C_Q=tau integral R^2 B + integral R S
        -epsilon partial_Z integral R^2 H.

Integration by parts gives integral Rg=E_u+C_P and integral R^2 Rg=F_u+C_Q.
Indeed integral (B_RR+B_R/R-B/R^2)=0, and its R^2-weighted integral is also
zero: the three interior contributions are 2 integral B, -integral B and
-integral B. Boundary terms vanish by compact support. Axial viscosity remains
in tau. Using all five linear rows, the complete new defects are exactly

    P_new      = E_u+C_P,
    Jtheta_new = integral R^2 g u,
    Jz_new     = E_g-(F_u+C_Q)/2.                         (1)

The new axial residual moment therefore contains the derivative of
E_g-(F_u+C_Q)/2+c_rho(E_u+C_P), with the full compensator product retained.

These identities hold under the initial-zero-Haar and frozen-W conditions in
section 1. They must not be applied unchanged at later stages with nonzero
old means or a changed wave: their additional cross terms must be recomputed.

## 3. Projected angular lower bound

In L^2(I,dR/R), the zero angular moment is <u,R^3>=0 and the pressure row is
<u,2V>=-P_in. Orthogonally project 2V off R^3 and apply Cauchy--Schwarz:

    D_P = 4 [integral_I V^2/R
             -(integral_I R^2 V)^2/(integral_I R^5)] > 0,
    E_u >= P_in^2/D_P,
    P_new-C_P >= P_in^2/D_P.                            (2)

Strict positivity follows because V=a_eff R^(-1-2 lambda) is not proportional
to R^3 on a nondegenerate interval for lambda>=0. This is the optimal bound in
the Hilbert relaxation with these two constraints. The actual five-map has
an additional angular constraint and prescribed smooth bumps; equality is not
claimed for it. In a stationary, Z-independent slow update the divergence row
forces B=c/R, and compact support forces c=0. Hence C_P=C_Q=0 and a nonzero
P_in leaves strictly positive P_new. For the evolving actual state, C_P can
have either sign; (2) does not rule out dynamic cancellation.

## 4. Axial flux costs quadratic inverse tilt

In L^2(I,R dR), <g,1>=0 and <g,RV>=-Jtheta_in. Projection now gives

    I0=integral_I R,
    D_z=integral_I R (RV)^2 -(integral_I R^2 V)^2/I0 > 0,
    E_g >= Jtheta_in^2/D_z,
    Jz_new+(F_u+C_Q)/2 >= Jtheta_in^2/D_z.              (3)

For V as above this is

    D_z/a_eff^2 = integral_I R^(1-4 lambda)
                 -(integral_I R^(1-2 lambda))^2/I0.

To quantify it without subtracting nearly equal floating-point moments, put
y=log R, ybar=(integral_I R log R)/I0,
L=integral_I R(log R-ybar)^2>0. For 0<lambda<=lambda_bar, define

    c=min(1,r_-^(-2 lambda_bar),r_+^(-2 lambda_bar)),
    C=max(1,r_-^(-2 lambda_bar),r_+^(-2 lambda_bar)).

The pairwise variance identity and the mean value theorem give

    D_z/a_eff^2 = (1/(2 I0)) integral_I integral_I
       R S [exp(-2 lambda log R)-exp(-2 lambda log S)]^2 dR dS,
    4 c^2 lambda^2 L <= D_z/a_eff^2 <= 4 C^2 lambda^2 L.       (4)

Here the symbol S in the double integral is a radial integration variable,
unrelated to the transport expression in section 2. The derivative of the
exponential in log R has magnitude between 2 lambda c and 2 lambda C.
The corresponding pairwise identity for log R proves (4). Dominated
convergence on this fixed finite interval also proves

    lim_{lambda->0} D_z/(a_eff^2 lambda^2)=4L.

Consequently the explicit necessary lower budget is

    Jz_new+(F_u+C_Q)/2
       >= Jtheta_in^2/(4 a_eff^2 C^2 lambda^2 L).          (5)

At lambda=0, RV is constant and the zero-flux row forces Jtheta_in=0.
For every fixed positive lambda the denominator is positive; no retuning of
the construction's tilt is proposed. The constants depend on the actual
reserved hull and a_eff. The limit is on a fixed hull; it does not establish
a uniform constant on moving hulls or across physical scales. This lifts the
inherited inverse's lambda^{-1} amplitude conditioning to a lambda^{-2}
nonlinear axial-energy budget. It does not assert that the actual target is
bounded away from zero or that (F_u+C_Q)/2 cannot cancel it.

There is also a joint budget without the unknown F_u. Compact support gives
F_u<=r_+^2 E_u=r_+^2(P_new-C_P). Therefore

    Jz_new + [r_+^2(P_new-C_P)+C_Q]/2 >= Jtheta_in^2/D_z.   (6)

In particular, exact simultaneous P_new=Jz_new=0 requires

    -C_P >= P_in^2/D_P,
    C_Q-r_+^2 C_P >= 2 Jtheta_in^2/D_z.                    (7)

This is a necessary two-defect cancellation gate, not sufficient closure of
all three defects. It allows dynamic cancellation explicitly. A stationary
Z-independent update has C_P=C_Q=0, so (7) excludes simultaneous exact
cancellation when either P_in or Jtheta_in is nonzero. It does not exclude
the actual evolving construction.

## 5. Checkpoint boundary and next obligation

The preceding common primary tail and support certificates remain unchanged.
Equations (1)--(7) are valid on each such chart if its inherited admissibility
conditions hold; they require no new independent tail choice. They provide
neither source admissibility for every subsequent harmonic nor a uniform
all-scale bound on these moment constants. Even if incoming targets have
size O(epsilon^alpha), their quadratic energy has a fixed-parameter cost
that can include lambda^{-2}; the linear five-map alone does not imply a
later improved signed-target exponent. The raw-state ceiling and cancellation
gate in `current_target_moment_audit` remain open obligations.

Next: bound the actual C_P,C_Q and g*u, plus all later old-mean and changing-W
cross terms, in the same weighted source norms on the common tail, including
both pressure tails. Join those estimates to the signed current-target cycle
before asserting any complete residual gain. Uniform generated-source
supports/jets, viscosity and cutoff interactions, and a full residual bound
are still unproved. No unconditional pulse bound or blow-up is claimed.

## 6. Reproducible evidence

`proof_checks.py` integrates a compact polynomial fixture directly through
the complete averaged radial operator and independently through C_P,C_Q. It
keeps nonzero time, axial transport, radial viscosity and axial viscosity;
omission controls detect their loss. It checks exact complete-covariance and
both-tail controls from the initial-mean source. An exact continuum moment
example at lambda=1/4 on [1,16] has D_z=99/85. Finite weighted-space controls
verify projection equality, the strict energy gap without projection, and
arbitrarily small nonzero axial denominators without roundoff. A joint
finite-space control admits simultaneous pressure/axial cancellation only
with the retained dynamic terms and verifies the eliminated budget (6).
These fixtures
are not the actual C-infinity profile, continuum suprema or a proof of
all-scale admissibility. Sections 2--4 supply the analytical arguments.

Run from the repository root (standard library only for this checkpoint):

```sh
python3 research/navier_stokes_cascade/results/core_covariance_mean_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_covariance_mean_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_covariance_mean_recheck/evidence.json \
  --output /tmp/core_covariance_mean_recheck_fresh.json
```

Six historical records and all of their declared source hashes are bound.
The CLI refuses to overwrite any output. CI repeats the tests and frozen-record
reproduction; passing them does not close the analytical obligations above.
