# Q34: correct the raw/common-chart to initial Haar norm interface

Parent: `b49bb07f3c36d7d8ecced0b940bc6403fac2e6f1` (Q33).
Q33's covariance/selected-flux inequalities are raw common-chart bounds.
Its statement that their returned polynomials supply the initial Haar-row
inputs needs a qualification: matching Python polynomial types is not
enough to match derivative norms or radial charts. The historical
`core_symbolic_initial_haar_join` requires through-order factorial R-chart
norms. This checkpoint supplies that conversion and the resulting initial
scalar input. Historical files and frozen evidence remain unchanged.

## 1. A precise obstruction to the unconverted interface

For five continuous coordinates x, distinguish

    M_n(f)=max_(|I|<=n) |partial^I f|,
    N_n(f)=sum_(|I|<=n) |partial^I f|/I!.

Even at a single point M_1(exp(x1+...+x5))=1 but N_1=6.
Thus N_n<=M_n is false; testing that a raw bound has the right polynomial
type does not establish the norm inequality. A valid conversion is

    N_n(f)<=C_n M_n(f), C_n=sum_(j=0)^n 5^j/j!.        (1)

This follows from the multinomial theorem and is sharp for the stated
pointwise exponential fixture. It is a finite-jet control, not an actual
field example. Q33's existing raw estimates and exact physical assembly
remain valid; their direct promotion to this normalized interface is not
used here. No new physical parameter, support or accepted band is chosen.

## 2. Covariance: reverse the radial chart before the Haar recipes

The native root/compact-map chart uses U=R^d with the retained
1/3<d<5/2. The four other coordinates are unchanged. The common-band
atlas transfer in Q33 is linear in U and does not perform this nonlinear
return to R. Use the existing `core_weighted_reverse_chart` bound on the
same closed positive hull Rmin<=R<=Rmax:

    B=max(1,Rmax^2,1/Rmin),
    a_j=4*1_(j=1)+3B/Rmin^(j-1), j>=1,
    e_0=1, e_n=sum_(j=1)^n (j/n)a_j e_(n-j),
    L_n=n! sum_(j=0)^n e_j.                            (2)

The a_j bound the sums of normalized coordinate increments. Finite Taylor
composition gives (2); its additional n! makes it a bound even for the
sum of raw output derivatives, hence also for N_n. This transfers a raw
through-n input bound, not a quotient differentiated after removing zeta.
All weights and poles stay evaluated at the same physical point. Q33
proved X,eta,delta,zeta invariant under native-band changes; R-to-U is
also a reparameterization of that same point. No inverse edge floor or
epsilon loss is introduced. Derivatives of delta are never taken.

Consequently, for Q33's raw W_n,Dcurl_n in U,

    W_R,n=L_n W_n, Dcurl_R,n=L_n Dcurl_n.              (3)

L_n ALREADY dominates the normalized derivative sum; multiplying it by
C_n again is unnecessary. The polynomial conversion preserves every
coefficient exactly. It does not replace the coordinate cost by one.

## 3. Background: evaluate in R and apply the missing norm count

Q32 can evaluate the higher stress directly in the common R,Z,T chart,
by omitting its optional R(U) composition. Let H_R,n=4 Cnative_n S_n be
that RAW R-chart constant after removing epsilon^3 zeta delta^(-3n-3).
It has no auxiliary derivatives, so (1) is conservative but valid. Set

    Hthrough_n=C_n H_R,n,
    D_R,n=Dcurl_R,n+Hthrough_n.                        (4)

Do not apply the U-to-R factor L_n to a background already bounded in R.
The new `r_inputs` deliberately accepts no optional U geometry for this
background route. It returns W_R,n,D_R,n in the historical polynomial
type. It retains Q33's epsilon orders 1 and 2-2kappa respectively, with
higher stress order3, and common pole p_n=12n+18. The exact cancellation
of the principal selected fluxes is unchanged. The remaining four tensor
entries and complete curl squares remain in W_R,n.

This is a norm correction, not a new physical cancellation or improved
epsilon rate. It does not assert that the unconverted Q33 constants happen
to fail numerically on the selected field; it establishes that the prior
argument did not justify that interface and supplies a sufficient one.

## 4. Feed the original initial scalar identities on the whole hull

For the INITIAL real primary with zero old means, retain the existing
compact-boundary and Haar identities in `core_sharp_initial_scalar_join`:

    P0=int(Wtt-Wrr)/R - epsilon partial_Z int Wrz,
    Jtheta0=int R^2 Wtheta,z,
    Jz0=int R[Wzz-(Wrr+Wtt)/2]
                         +(epsilon/2)partial_Z int R^2 Wrz. (5)

Every tensor here is Haar averaged. Let [a,b] contain the WHOLE primary
and fast potential support, including A1; it need not be the reserved I4.
With M_j=(b-a)max(a^j,b^j), choose a single E dominating
sqrt(zeta)delta^-p for the COMMON finite pole p>=p_(n+1), at the actual
fixed width. Put Omega_j=E W_R,j. Since zeta<=sqrt(zeta), this removes
the weight once. Then the safe normalized through-order estimates are

    IP_n=2M_-1 Omega_n+(n+1)M_0 Omega_(n+1),
    Itheta_n=M_2 Omega_n,
    Iz_n=2M_1 Omega_n+((n+1)/2)M_2 Omega_(n+1).        (6)

The factor n+1 pays N_n(partial_Z f)<=(n+1)N_(n+1)(f).
There is no graph derivative cost after the EXACT full Haar and radial
identities; taking absolute values earlier would not justify that claim.
`initial_scalars` implements (6) and retains the off-diagonal radial/axial
pressure term. It does not compute the current postslow P or reset old
means. Fresh fast quadratic changes, the slow five-map remainder, the
moving c_rho P product and the two distinct projected moments remain
the next composition steps. An n+1 request is before expanding the
complete-curl/root dependency and Q14's primitive ceiling.

## 5. Validation and scope

Independent tests exhibit the norm counterexample, compare mixed
coordinate derivatives with (2), preserve exact covariance polynomials,
and exercise (6) and its extra pressure jet. Parent and reverse-chart
controls are replayed. Tests are finite manufactured diagnostics, not
continuum estimates of the enormous physical construction.
At inspection on 30 October 2026 Q33 focused CI 41102631384 passed;
broad CI 41102631538 was still running.

Next insert these correctly typed R-chart inputs into the complete sharp
fast/slow scalar and projected signed-target chain, keeping the actual
hulls, symbolic edge constants and full pressure memory. Actual canonical
and fixed-map premises remain conditional. No unconditional source bound,
complete residual improvement, all-scale pulse theorem or blow-up follows.
Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json` here.
