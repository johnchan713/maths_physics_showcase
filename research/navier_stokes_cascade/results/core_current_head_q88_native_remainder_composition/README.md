# Q88: compose the nonlinear radial bound on one native hull

Parent: `58bbffdde9087aef859429849d955942ef867802` (Q87).
Q84 supplies a positive remainder polynomial; Q85--Q87 supply its wave,
mean and angular coefficients in the inherited R chart. This checkpoint
composes them in native U=R^d on the SAME first-real pressure state. Each
weighted input is unweighted at its own proved pole before composition.
It adds no coefficient-free assertion about the actual physical norms.

## 1. Derive a common native hull rather than rename the chart

The profile and original root chart share 0<r<=R<=H. Retain the original
1/3<d<5/2. A safe enclosure of U=R^d is

    u0=min(1,r³), u1=max(1,H³).                              (1)

For R<=1, R^d>=R³; for R>=1, R^d>=1. The upper bound is analogous.
Thus (1) applies to the same hull for every admitted fixed exponent; it
is not an independent assertion about a new integration interval.
Use Q87's finite positive normalized-R to factorial-through-U factor L_N
with Rmax=H, Umin=u0. The four other variables are unchanged. This preserves
the epsilon order and loses no derivative, though the fixed constants may
be enormous. Do not use the raw-to-factorial factor 243 as a radial chart
conversion: it performs a different operation.

The native reciprocal rfun=1/R and slow graph coefficient g=dU/dR need
not be imported from an R-chart norm with an unstated conversion. Put
Bgeom=max(1,H²,1/r), Qgeom=max(1,1/u0). At native order N sufficient
factorial-through bounds are

    r_N=r^-1 (N+1) binom(N+2,2) Qgeom^N,
    g_N=3 Bgeom (N+1)² Qgeom^N.                             (2)

Indeed rfun=U^(-1/d), with |1/d|<3, and g=d U^(1-1/d), with
|1-1/d|<2. Their normalized degree-j derivatives are bounded by their
value bounds times respectively binom(j+2,2)u0^-j and (j+1)u0^-j.
The value bound for g is 3Bgeom. Summing through N proves (2).
These are pure geometric multiplier estimates; a graph derivative of a
mean still carries the original M<=epsilon^-kappa loss in Q84.

## 2. Construct every field input on that hull

At radial-remainder output m put N=m+6 and j=m+4. Require the Q47
Postslow factory reserve m+7, as derived in Q87. The composition checks
Q86's same-state/root/full-pressure identities, then obtains:

* x,u,z from Q86.fields(N), followed by L_N;
* W,Z from Q85.defect(state,N), first unweighted using their returned
  covariance pole and original width, then multiplied by L_N;
* e from Q87.native_coefficient(j,Umin=u0), retaining its exact angular
  cancellation, same signed covariance and native conversion;
* b0,V,G,c from the original initial-budget fixed(N,p) coefficients,
  interpreted as their retained unweighted R-chart base bounds and
  converted by L_N; r_N,g_N from (2).

The base b coefficient is the original epsilon-stripped radial base
bound, not a deletion of the base radial velocity. The existing fixed
callback's actual complete-background/unweighted contract remains required.
The same full base and original inverse callbacks remain attached to the
factory; no new field, covariance or outgoing angular callback is accepted.
The adapter rechecks pressure-state identity at each composition.

This route uses an unweighted native estimate AFTER the input edge bounds.
It does not embed all fields into one arbitrarily renamed weighted pole.
The initial Haar, signed covariance, old fast and new fast fields retain
their distinct proved edge factors before unweighting. Their closed source
hulls, whole reserved potential support and accepted fixed-q domain must
still be the same. The later compact output has its own existing transition
weight floor and raw/native conversion; none is inferred from input norm
finiteness. Both exceptional pressure orders, moving compensator, axial
A1 velocity, signed squares and all scalar refreshes remain upstream.

## 3. A direct unweighted temporal bound in the native chart

The temporal inverse does not change U, so no radial chart factor belongs
to its output. We give a deliberately coarse unweighted proof instead of
silently discarding weighted conversion constants from an older recipe.
For vt=(1,1+sqrt(2)), nonzero integer k obeys

    |vt.k|^-1 <= 3 |k|_1,

by multiplying its quadratic conjugate and using the nonzero integer norm.
On a shell K<=|k|_1<2K there are at most 25K² modes; the inverse symbol
is at most 6K after dropping favorable 2pi constants. Four derivatives in
the larger auxiliary coordinate bound each coefficient by 16K^-4 times
a raw fourth-derivative norm. The shell sum is at most 2400/K and the
dyadic sum at most 4800. An additional factor two safely retains the Haar
projection. Arbitrary fixed coefficient derivatives commute with this
frozen inverse. Raw input maxima are bounded by (m+4)! times the factorial
through norm, and the output factorial sum costs at most 6^m. Therefore

    |c_i^-1 Nt^-1 Q f|_m
      <= 9600 6^m (m+4)! c_i^-1 |f|_(m+4).                (3)

Use the ORIGINAL compact(0,p)['I'] callback, whose contract is I0>=c_i^-1/S.
Then the required Q84 temporal polynomial is

    T_m(S)=9600 6^m (m+4)! I0(S) S.                       (4)

This is an unweighted same-torus inequality, not a new physical inverse or
a weighted theorem at an unknown pole. Periodicity, smooth input, the frozen
label, nonzero-mode projection and the original c_i floor remain premises.
It removes an independent temporal-norm callback from this particular
unweighted composition. No Gaussian reciprocal, new band or pressure reset
is introduced. Constants in (2)--(4) are fixed for each m, not uniform as
m tends to infinity.

## 4. Evaluate the same positive polynomial over the coefficient ring

The executable adapter evaluates Q84's formula over the inherited positive
EP coefficient ring. The old Q84.envelope function is a rational control
evaluator and rejects EP objects; it is NOT bypassed by claiming EP objects
are measured field values. Independent exact tests compare the new polynomial
with that rational evaluator on manufactured inputs.

The resulting EP polynomial R_m(S) bounds Q Rem/epsilon^(2-3kappa)
on the retained physical input hypotheses. At m=n+199996, Q82's shell
count 99998 gives its conditional compact contribution at 3-5kappa.
The principal part retains Q83's separate leading table and normal/column/
propagator hypotheses. Factory reserve n+200003 and conservative profile
reserve n+200030 are retained. Those indices do not expand the entire
nested historical source DAG or evaluate actual giant constants.

This completes an explicit native coefficient composition for the nonlinear
radial remainder, relative to the stated original contracts. It does NOT
certify their actual continuum inequalities or the complete physical residual.
The next obligation is to assemble the first-real regular and additive-tail
ledgers on the same exact-parameter domain, keeping the stage-specific axial
and radial pressure sources. An improved later sourced regular residual cannot
be paired with these first-real compact tails. Uniform actual-pulse bounds,
infinite realization, full residual improvement and blow-up remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. Low-order typed factory fixtures use
manufactured coefficient leaves; exact algebra and geometry checks do not
prove the actual physical norm contracts.
