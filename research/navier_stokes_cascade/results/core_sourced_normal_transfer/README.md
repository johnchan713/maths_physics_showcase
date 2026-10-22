# Original normal floor transferred through the sourced inverse

Parent: `772380fceca23a1a99d12749b238c494c4e88897`.
The original accepted-tail normal implication now feeds the existing first
sourced inverse with **q = 4/bmin^2**, without the stronger optional
phase-beta comparison gate. At every fixed derivative order this changes
finite constants only: it adds no powers of S, inverse Gaussian, or new
support loss. The source, coefficient-jet and tangent-propagator bounds
remain conditional inputs. This does not close the complete pulse or residual.

## 1. Join the existing estimates without changing a frozen field

`core_normal_length_closure` proves |n|^-2 <= q on the original rule-defined
accepted primary tail, conditional on the validity of its ORIGINAL primitive
bounds. Its original Riccati length already dominates the normal threshold.
Recompute that certificate on the supplied old band and retain its original
length semantics. Do not select a new band, enlarge its old length budget,
or replace its primitive inputs by larger later estimates.

The positive recipes in `core_phase_coefficient_budget` for full n, nprime,
K and P'/P use no inverse-normal floor. Consequently they remain valid
with the original denominator. Use the original [Rmin,Rmax], u, lambda_max,
L>=1, and the already retained carrier bounds

    |p| <= Pmax = Rmax*bmax*(1+u/gmin)+1,
    |pz| <= Zmax = bmax*(1+u/gmin),
    |x0| <= bmax*u/2.

These follow from the same frozen carriers and rounding rule; using L>=1
is conservative. There is no need to compare this loose bmax with a
separately loose sqrt(lambda_max) bound. Full F,G,btilde,path jet bounds
must still refer to the same physical neighborhood and complete background.
In particular G_Z and full H_Z=pF_Z+pzG_Z are retained. No pressure memory
endpoint is truncated to a chart prefix.

`forced_from_original` constructs this PhaseBudget and passes its callbacks
to the UNCHANGED `ForcedBudget`. It recomputes the old certificate, rather
than trusting a supplied success flag. Its rational q is exactly 4/bmin^2.
The older `normal_gate`/`PhaseBudget.forced` route and frozen evidence remain
historical sufficient estimates; their stronger 4/beta_phase is not inferred.

## 2. All-order dependence on the reciprocal constant

Hold all other input polynomials fixed: normal jets N_j, moving normal,
K, source f_j, profile rate, tangent T(S), length L(S), and cutoff jets.
Regard the scalar upper bound q>0 as a formal variable. All recurrences
are positive polynomials in q with positive-polynomial coefficients in S.
This is algebra on majorants, not a claim that the exact solution is a
polynomial in q. Write R_j for the reciprocal-normal-square jet bound.

    R_0=q,
    R_j=q sum_(i=1..j) binom(j,i) (N*N)_i R_(j-i).

Thus deg_q R_j <= j+1 by induction. Raw Leibniz convolution adds degrees;
it does not introduce new denominators. The remaining bounds are:

| Budget | Upper degree in q |
| --- | ---: |
| normal square, pressure row n^T K - nprime^T, damping | 0 |
| reciprocal R_j, lift h_j, undamped generator, normal projector, adapted generator at order j | j+1 |
| orthogonal source projection at order 0 | 0 |
| differentiated source projection at order j>=1 | j+1 |
| ambient propagator majorant | 2 |
| velocity jet V_j | 4j+2 |
| pressure and potential jets | 4j+3 |

For the ambient propagator, Q=A0*h-hprime has degree at most two and

    Kamb = T*(1+L*Q*N_0) + h_0*N_0

has degree at most two. This uses the Gaussian-adapted ambient extension
already derived in the sourced inverse; a naive normal extension cannot
be substituted. The required tangent estimate remains

    ||V_T(s,w)|| <= T(S) P(s)/P(w),   0<=w<=s<=L_actual,
    L_actual <= L(S), T(S)>=1.                              (1)

For velocity, V_0 has degree at most two. For j>=1 the parameter term is

    Kamb*L*(projected_f_j + sum_(i=1..j) binom(j,i) Ahat_i V_(j-i)).

Its largest possible degree, using the induction hypothesis, is

    2+(i+1)+[4(j-i)+2] = 4j+5-3i <= 4j+2,   i>=1.

The pulse derivative terms have degree at most 4j-1, and the direct
parameter source term at most j+3 <=4j+2. Pressure and potential have
degree at most 4j+3, by convolution with a reciprocal or lift jet.
The ambient factor is charged at EVERY differentiated inverse step;
omitting this repeated cost would undercount the q degree.

With D_velocity(j)=4j+2, D_pressure(j)=D_potential(j)=4j+3, let

    r=max(1,q_new/q_old).

Positivity gives the coefficientwise estimate, for every S>=0,

    B_kind,j(S;q_new) <= r^D_kind(j) B_kind,j(S;q_old).        (2)

This compares algebraic majorants at identical other inputs. If an old q
is unproved on the actual field, its old output is merely a comparison
polynomial; only the new q supplied by the original certificate is used
as a valid inverse bound. The q values must be uniform constants across
the family for the 'no new S powers' conclusion. Apply (2) termwise after
the SAME cutoff: lower derivative degrees are no larger than D_kind(j).
This does not bound a feedback source that itself changes with q; such
source dependence must first be included in its own budget. No complete
covariance/mean iteration is silently closed by (2).

## 3. Full pressure and support do not disappear

The normalized ODE pressure contains BOTH

    i/ m * [(n^T K - nprime^T)t/|n|^2 + n^T f/|n|^2].       (3)

An overall 1/k is tracked by the existing physical normalization. The
unprojected n.f term survives although the velocity equation uses Pi*f.
For n=e3, K=nprime=0 and f=g e3, zero-entry velocity is exactly zero but
pressure is i*g/m (before 1/k). Projecting the pressure source would erase
this term incorrectly. The new bridge uses the existing full pressure
formula, moving-normal term and SAME potential/pressure cutoff derivatives.
Pi1/Pi2 global pressure memory and axial mean terms remain separate complete
mean-map obligations. Orthogonal projection also does not prove support:
the source must vanish on the whole relevant parameter path to infer a
zero uncut response there. Endpoint source values alone are insufficient.

## 4. Progress ledger, validation and next obligation

Removed: the separate Cphi/phase-beta value gate from the first sourced
inverse when using the original accepted-tail normal certificate.
Added analytical assumptions: **none**. Retained: valid original primitive
inputs and accepted-tail membership; actual full-background/phase/source
jets on one domain; (1); original supports and fixed cutoffs. The bridge
does not identify an unknown historical numeric band with that tail.

Exact formal-q controls reuse every positive recurrence in the inherited
inverse, attaining the output degrees through order three, and
compare with rational substitution at two S values. Independent controls
check coefficientwise dilation, pure-normal-source pressure, moving normal,
and cutoff products. These are finite manufactured regression controls;
the all-order conclusion follows from the induction above. They are not
continuum enclosures. Three historical records and their declared sources
are bound; historical files are unchanged.

Next: establish (1) for the SAME canonical entrance and source domains,
then propagate the generated complete covariance through fresh axial and
full-memory pressure means. Full physical-scale derivative/support bounds
and viscosity/cutoff/all-interaction residual improvement remain open.

```sh
python3 research/navier_stokes_cascade/results/core_sourced_normal_transfer/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_sourced_normal_transfer/audit.py \
 --verify-record research/navier_stokes_cascade/results/core_sourced_normal_transfer/evidence.json \
 --output /tmp/core_sourced_normal_transfer_fresh.json
```
