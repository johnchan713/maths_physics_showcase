# Finite current scalar-cycle budget with scale losses retained

Parent: `b2d5274a30b1fcf8fef6982bf3e051d2805682eb`.
The last checkpoint routed higher angular scalar jets into pressure through
pure-time flux and integrated-residual bounds. The conditional cycle already
states the scalar order C+1 after a fresh slow solve. The new result here is
an explicit, noncircular norm recurrence for that scalar through BOTH wave
updates, fast means and fresh slow means, followed by its substitution into
the higher-jet pressure inequalities. Polynomial scale losses are retained.

This discharges the angular scalar routing relative to the specified current
map/seed bounds. It does not independently prove those operator inputs or
close the other scalar defects, full pointwise residual or infinite cycle.
Historical records, profiles, tilt, width and the common tail are unchanged.

## 1. Exact stage structure preserves the sharper old Haar rate

Use Q=2^-ell, epsilon=Q^h, S=ell^2>=1 and kappa=10^-5. Fix any finite
cycle index and any finite normalized Z,T derivative order n. Constants
below may depend on that index and order, fixed inverse tilt and profiles.
The norms are the unweighted mixed slow coefficient norms of the parent.
Finite S-power losses can be majorized by nonnegative integer powers since
S>=1. Passing from the retained edge-weighted field classes to these input
norms still uses their source/edge estimates; none is inferred from support
alone.

The canonical initial mean input has Haar tangential velocities of order
one and Haar radial velocity of order two. Every actual fast increment,
including its curl remainder, has zero Haar mean. Every wave update leaves
mean velocities fixed. Every later slow tangential increment has order
C>1, with radial increment of order C+1. The prefix slow orders
2-3kappa and 3/2-kappa and all later C>=C0=2-2kappa exceed one.

Consequently the accumulated Haar velocities are exactly the initial Haar
velocities plus the FINITE sum of preceding slow increments. On the current
reserved hull there is a polynomial a_n(S) with nonnegative coefficients
such that

    |bar(v)|_n, |bar(gamma)|_n <= epsilon a_n(S),
    |bar(beta)|_n <= epsilon^2 b_n(S).                      (1)

For example an initial majorant a0,n(S) and past slow majorants l_j,n(S)
give a_n=a0,n+sum_j l_j,n, since epsilon^(C_j-1)<=1. This is a finite
sum, with no limiting summability or new threshold assertion. On adding a
new slow step the ledger updates by adding its polynomial. It does not
replace the full mean rate mu=1-kappa with one: only the Haar part has (1).

## 2. Explicit preslow scalar budget through current updates

Let I contain the complete supports used in the angular moment, including
the current full wave, and let M2=(d-a)max(a^2,d^2), I=[a,d]. A separate
reserved-hull constant M2s can be used for the slow-only output. No wider
wave support is assigned to I4.

Suppose the current scalar has the bound epsilon^C j_n(S). At each of the
sourced and signed wave steps (s and d), assume the retained current map
bounds for the COMPLETE velocities

    |w_l|_n <= epsilon^(1/2) w_l,n(S),
    |q_l|_n <= epsilon^(C-1/2) q_l,n(S), l=s,d.             (2)

The old wave at the second operation is the UPDATED wave. Its polynomial
must bound that current field. Each covariance change is exactly
bar(w tensor q+q tensor w+q tensor q), so its theta-z contribution changes
J by at most

    M2[2 epsilon^C w_l,n q_l,n
                       +epsilon^(2C-1) q_l,n^2].          (3)

Other covariance entries still enter the other moment and pressure equations;
they have not been set to zero. No signed-polarization cancellation of the
theta-z entry is assumed. Formula (3) uses complete velocities, including
both quadratures, curl/localizer terms, old/new crosses and squares.

For the fast step let old full tangential means have norm at most
epsilon^mu o_n(S), and each actual fast tangential component have norm at
most epsilon^H f_n(S), mu=1-kappa, H=C-kappa. Exact current subtraction gives

    |delta J_fast|_n <= M2[2 epsilon^(mu+H) o_n f_n
                              +epsilon^(2H) f_n^2].       (4)

For C>=C0 the extra exponents relative to C are

    C-1>0, mu+H-C=1-2kappa>0, 2H-C=C-2kappa>0.

Thus on the SAME domain epsilon<=1, without selecting a later first band,

    |J_postfast|_n <= epsilon^C Pi_n(S),
    Pi_n=j_n+M2 sum_(l=s,d)(2 w_l,n q_l,n+q_l,n^2)
                +M2(2 o_n f_n+f_n^2).                    (5)

Every input polynomial in (5) refers to the current state at its operation.
This inequality uses no assumed outgoing scalar order. It supplies both
pure-time jets and, by M=epsilon partial_Z J,

    <J_postfast>_n <= epsilon^C Pi_n(S),
    |M_postfast|_(n-1) <= n epsilon^(C+1) Pi_n(S), n>=1.   (6)

It therefore supplies both scalar channels required by the parent pressure
join, relative to (2), the fast-map inputs and the incoming scalar bound.

## 3. Fresh slow cancellation and the next scalar input

The fixed five-map reads the RECOMPUTED postfast P,J,Jz and satisfies the
current angular linear row int R^2(G u+V g)=-J_postfast. With W held fixed,
exact recomputation gives

    J_new=int R^2[bar(gamma)u+bar(v)g+u g] dR.              (7)

On the reserved hull use (1) and the actual tangential inverse bounds

    |u|_n <= epsilon^C U_n(S), |g|_n <= epsilon^C Z_n(S).

Here U depends on CURRENT P,Jz and the angular block inverse; Z depends on
CURRENT J and the axial inverse, including its fixed 1/lambda cost.
Their tangential estimates do not lose an input slow derivative. The radial
curl component does, and remains necessary in the full residual. It follows

    |J_new|_n <= M2s[epsilon^(C+1)a_n(U_n+Z_n)
                                     +epsilon^(2C)U_n Z_n]
      <= epsilon^(C+1) Theta_n(S),
    Theta_n=M2s[a_n(U_n+Z_n)+U_n Z_n], C>=1.              (8)

Using only the weaker full mean rate mu would give C+1-kappa and would
not justify (8). The sharper Haar ledger (1) is essential. The desired
next cycle order Cnext=C+1/2-kappa is below C+1 by 1/2+kappa, so (8)
supplies the next incoming angular scalar polynomial j_next,n=Theta_n.
The new Haar polynomial is a_next,n=a_n+U_n+Z_n. These are finite-stage
updates, not uniform-in-stage or infinite-series estimates.

Other defects P,Jz must still have their current S_C bounds for the five-map.
Their radial/time/viscous budgets require higher input orders. In particular
the existing m+2 radial-increment budget and its curl construction can
require angular target jets through m+3. The finite smooth cycle is allowed
those higher jets; no C^m-to-C^m contraction or infinite derivative ladder
has been proved here. The initial canonical seed remains an input.

## 4. Substitute into the pressure estimates without losing logarithms

Let d_m and k_m be the retained fixed-profile mismatch and inverse-coefficient
norms. From (5), the sharper direct bounds for the two mismatch contributions
in the axial equation are

    |epsilon^5 partial_Z[D_rho partial_Z^3(KJ_postfast)]|_m
       <= q_4(m)d_(m+1) k_(m+4) epsilon^(C+5) Pi_(m+4)(S),
    |epsilon^3 partial_Z[D_rho partial_T partial_Z(KJ_postfast)]|_m
       <= q_3(m)d_(m+1) k_(m+3) epsilon^(C+3) Pi_(m+3)(S). (9)

Alternatively inserting (6) into the parent two-channel formulas adds at
most factors m+5 and m+4 respectively. Every moving-density commutator
is retained by substituting |KJ|_n<=k_n epsilon^C Pi_n into the full
six-term budget of the fixed-profile checkpoint. Postslow versions use
C+1 and Theta. These conclusions concern the identified pressure terms;
all other viscosity, cutoffs and the complete residual remain to be bounded.

The polynomials Pi and Theta must NOT be replaced by constants independent
of ell. A term A epsilon^eta S^p with positive margin eta can be tested on
an ALREADY accepted tail ell>=ell_*>=1 by a rational sufficient condition.
For integer p>=0 set N=2p+2 and a=h eta/2. Since log 2>=1/2 and
exp(a ell)>=(a ell)^N/N!,

    A epsilon^eta S^p <= A N!/(a^N ell^2)
                        <= A N!/(a^N ell_*^2).            (10)

For noninteger nonnegative p, replace it by its ceiling. Sum the finite
term bounds and compare with the tolerance. Formula (10) checks a given
tail; it does not replace or certify the actual ell_* or field constants.
At eta=0,p>0 the loss is unbounded, and at eta=p=0 one needs A<=tolerance.
No all-orders or all-stage common threshold follows from finitely many checks.

Next independently discharge the retained current wave, fast/slow inverse,
other-defect and canonical seed inputs, including their full supports and
edge weights, and join all pressure, viscosity and nonlinear terms. The
angular scalar input is now routed through those SAME hypotheses instead
of an extra independent high-jet assumption. No unconditional pulse,
complete residual improvement or blow-up theorem is established.

## 5. Evidence

Rational positive-polynomial and exponent ledgers check the norm recurrence
and its whole-half-line margins. Independent current-state evaluations
verify fresh angular-row cancellation with nonzero old Haar means and a
fixed complete covariance. Dropping the old cross terms or using a stale
linear row is detected. Exact rational exponential-series bounds check
(10) without constructing enormous physical q. Historical evidence remains
unchanged. Diagnostics support the algebra; sections 1--4 prove the stated
conditional estimates and explicitly retain their premises.

```sh
python3 research/navier_stokes_cascade/results/core_scalar_cycle_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_scalar_cycle_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_scalar_cycle_budget/evidence.json \
  --output /tmp/core_scalar_cycle_budget_fresh.json
```
