# Pointwise pressure-shape gate beyond the scalar moment cancellation

Parent: `bbac37527263fbd9bc91e277c1e56bcc0aacceef`. The previous checkpoint
identified an exact second-moment condition cancelling the principal fourth
target derivative in the INTEGRATED axial residual. Here we derive the
corresponding pointwise-in-radius coefficient in the Haar-averaged pressure
and axial equation. The scalar moment condition is insufficient to remove it.

The density rho_* below is a comparison object derived from the existing
axial potential, not a replacement for the actual fixed compensator rho.
No actual profile, pressure choice, tilt, width or common tail is changed.
All statements about high target jets concern the isolated slow-repair
linear response, with the same scope as the parent. Actual canonical targets
are linked to their incoming fields; no failure of that sequence is proved.

## 1. A fixed positive reference density from the prescribed potential

Use the parent's phi0=eta0, phi1=g^-1 eta0(x/g), g>1, with ordered disjoint
nonnegative smooth bumps. Define the TARGET-INDEPENDENT functions

    f=phi1-g phi0,
    A(x)=integral_left^x y f(y) dy,
    psi_f(x)=A(x)/x,
    D=g mu_1 log g,  mu_1=integral x eta0(x) dx,
    rhohat_*(x)=-psi_f(x)/D.                                (1)

The first bump makes A nonpositive, the gap keeps A=-g mu_1<0, and the
second bump raises it monotonically back to zero. Thus rhohat_*>=0,
smoothly compactly supported in the same potential hull. The parent's
logarithmic identity gives integral rhohat_*=1. In the gap it is exactly
1/(x log g), so it has a positive floor on every smaller closed gap.
This definition works also when the actual target vanishes: no target is
divided by. It is independent of lambda and a(eta), since those factors
cancel on normalization.

Transfer it to the same Q-chart:

    rho_*(R,Z,T)=s^-1/2 rhohat_*(R/sqrt(s)),
    F=K Jtheta,  K=s^h mu_1 log g/[a mu_(1-2lambda)(1-g^(-2lambda))].

The exact axial inverse then factorizes as

    Psi_R=-F rho_*,  B=-epsilon partial_Z Psi_R
                          =epsilon partial_Z(F rho_*).    (2)

These identities include the full moving shape and physical normalization.
The previous coordinate/support certificate supplies their fixed-jet chart
geometry costs on every accepted band; the higher target and a^-1 bounds
remain separate input obligations.

## 2. Exact pressure source and its whole moving-shape commutator

Isolate the slow-time plus axial-viscous radial source
g_v=epsilon B_T+epsilon^3 B_ZZ from the full linearized equations. Its total is

    P_v=epsilon^2 F_TZ+epsilon^4 F_ZZZ =: U(F).

For the ACTUAL rho, reconstruct its Haar pressure by the ordinary compact
primitive p_v=I0(g_v-rho P_v), I0 h=integral_left^R h. Both densities have
unit mass, with smooth compact extension. Product differentiation gives EXACTLY

    g_v-rho P_v = (rho_*-rho) U(F)
      +epsilon^2[rho_*,Z F_T+rho_*,T F_Z+rho_*,TZ F]
      +epsilon^4[3rho_*,Z F_ZZ+3rho_*,ZZ F_Z+rho_*,ZZZ F].    (3)

Every derivative acts on the complete moving density. Other radial vector
viscosity, time, transport, old means and squares stay in the full equation;
(3) is an identified component, not a declaration that they vanish. Complete
wave covariance is held fixed during this mean operation, so it cancels in
the difference while remaining in both incoming/outgoing equations. The
actual A0 pressure has zero Haar mean and is not removed from the full field.
Both Pi1/Pi2 background tails and their physical cutoff derivatives are retained.

Put

    D_rho=I0(rho_*-rho).

This is compact at both ends because the density difference has zero total.
The principal third target-Z-derivative coefficient in p_v is
epsilon^4 K D_rho, and the principal fourth coefficient in the pointwise
Haar-averaged axial equation (+epsilon partial_Z p_v) is

    epsilon^5 K D_rho(R,Z,T).                              (4)

All remaining terms in (3) contain at most two Z derivatives of F before
the axial gradient. The other linear mean-equation groups have lower
target-Z-derivative order as in the parent. Therefore this principal
coefficient vanishes for all R if and only if rho=rho_* as functions, since
K>0 and partial_R D_rho=rho_*-rho. Matching one scalar moment is weaker.

Indeed the exact integration-by-parts identity is

    integral R D_rho = c_rho-c_*,
    c_*=(1/2)integral R^2 rho_*=K2/(2K).                    (5)

Thus (4) integrates to precisely the parent's coefficient
epsilon^5(c_rho K-K2/2). Its scalar cancellation says only
integral R D_rho=0; it does not imply D_rho=0. The same distinction matters
for the compact signed axial target: after that scalar cancellation,
the fourth-jet pressure response is -epsilon^5 K I1(D_rho), with
I1 h=R^-1 integral_left^R r h(r) dr. A nonzero D_rho can survive this
zero-moment primitive too.

## 3. A positive smooth countercontrol with identical scalar moments

Choose nonzero C-infinity w compactly inside the positive gap of rhohat_*.
Let h=w''' and choose a sufficiently small nonzero constant delta. Then

    rhohat_alt=rhohat_*+delta w'''

is still nonnegative with unit mass. Integration by parts proves
integral w'''=integral x^2 w'''=0, so its second moment ALSO agrees exactly
with rhohat_*. Nevertheless at s=1,

    D_alt=-delta w'',  integral x D_alt=0,
    I1(D_alt)=-delta(x w'-w)/x.

Both functions are nonzero for nonzero compact w. Consequently the scalar
fourth-derivative coefficient cancels while the pointwise axial and signed
target coefficients (4) survive. The construction remains valid after the
actual dilation: D_alt(R,Z,T)=-delta w''(R/sqrt(s)). This countercontrol
demonstrates insufficiency of mass, second moment, positivity and common
support alone. It is not asserted to be the actual chosen compensator.

For an explicit positivity budget let U_gap be an upper bound for x on
the perturbation support. Use log g<=g-1 and
rhohat_*>=1/[U_gap(g-1)]. If M>=||w'''||_infinity, taking
0<|delta|<=1/[2U_gap(g-1)max(1,M)] keeps at least half the reference density
on the perturbation support. No small tilt or pressure-field approximation
is used to establish this admissibility.

## 4. A shape match moves the highest jet; it does not remove every residual

If the actual rho happened to equal rho_*, the mismatch term in (3) would
vanish, but all moving-density commutators there would remain. Moreover the
Haar radial residual identity is always partial_R p_v-g_v=-rho P_v. Its
principal third target-jet coefficient is

    -epsilon^4 K rho(R,Z,T).                               (6)

No unit-mass compensator can make this coefficient identically zero.
Its radial L1 norm is at least epsilon^4 K, because integral |rho|>=1.
If its support has length Lrho, its supremum is at least epsilon^4 K/Lrho.
These are bounds for the isolated linear-response COEFFICIENT, not numerical
lower bounds for the actual full residual. They prohibit claiming a full
mean residual gain from a compensator shape choice alone; actual target-jet
bounds or linked-cycle cancellations remain necessary.

## 5. Quantitative pressure-shape norm input

For any compact density difference h=rho_*-rho of zero total,

    ||I0 h||_infinity <= (1/2)||h||_L1.                     (7)

Its positive and negative masses are equal, and any cumulative integral
lies between their opposite extrema. The constant 1/2 is sharp for separated
positive/negative bumps. All slow derivatives of h also have zero total,
so the same bound applies to each differentiated primitive. Radial
derivatives instead satisfy partial_R^k D_rho=partial_R^(k-1)h, k>=1.

Let H_m=sum_(i+j<=m) ||partial_Z^i partial_T^j h||_L1/(i!j!) uniformly on
the slow domain. For the parent's unweighted normalized slow supremum norms,
|D_rho|_m<=H_m/2. In particular the mismatch part of the axial pressure
gradient obeys the sufficient bound

    |epsilon^5 partial_Z[D_rho partial_Z^3(KJ)]|_m
      <=epsilon^5 (m+1)(m+2)(m+3)(m+4)
          (H_(m+1)/2)|K|_(m+4)|J|_(m+4).                  (8)

This includes the derivative of D_rho; (3)'s other terms need their own
retained lower-jet bounds. Fixed-profile density jets and their uniform
coordinate pullbacks give band-uniform H_m if those profile bounds are
supplied. The actual D_rho or its norm must be controlled for pointwise
pressure, rather than inferred from (5).

Next identify/bound the actual compensator's full shape mismatch and prove
the actual higher target-jet budgets in (8) and the radial term (6), or
show their cancellation in the complete linked cycle. Generated-source
supports, all viscosity/cutoff interactions and full residual closure are
still unproved. No unconditional pulse or blow-up theorem is claimed.

## 6. Evidence

Exact radial polynomial/logarithmic integrals normalize the prescribed
reference density. An explicit compact perturbation preserves both moments,
has a rational positivity certificate, and gives nonzero compact pressure
and signed-target primitives. Independent plus/minus evaluations of the
full averaged radial operator isolate the high-jet source before actual
pressure compensation. Separate moving-density algebra reconstructs the
entire compensated pressure source and detects omitted commutators. The
finite fixtures have finite endpoint regularity; the smooth construction
and estimates above use C-infinity profiles. Diagnostics are not continuum
suprema or actual canonical fields. Historical evidence is unchanged.

```sh
python3 research/navier_stokes_cascade/results/core_pressure_shape_recheck/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_pressure_shape_recheck/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_pressure_shape_recheck/evidence.json \
  --output /tmp/core_pressure_shape_recheck_fresh.json
```
