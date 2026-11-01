# Q85: same-root covariance defect for the nonlinear radial envelope

Parent: `def77cfab82a2d941fbaa8d8f3847abef66558e8` (Q84).
Q84 requires both the full FINAL wave covariance at order one and its
difference from the ORIGINAL principal tensor at order 2−2kappa. Here those
two coefficient inputs come from the existing Q37 same-root covariance
object and Q51's stronger signed target, rather than independent callbacks.
This remains conditional on the actual Q51, root/curl/atlas and common-domain
contracts. No older field, root, target, support or frozen evidence changes.

## 1. Distinguish a stronger norm estimate from a new signed wave

Q37's ledger retains alpha=2−3kappa for its signed target. Q51 proves
order TWO for that same target on its canonical hypotheses, with the same
positive coefficient bounds. The definition delta a=H^-1 Xi/(2 epsilon a)
is linear in the FIXED target Xi; curl is linear in the resulting potential.
Consequently a target bound epsilon^2 times its coefficient polynomial
propagates to signed principal and curl coefficients at orders 3/2 and
2−kappa. It is not legitimate to infer this merely from the old alpha
entry in the ledger. The new adapter records the rate-2 premise explicitly.
It does not multiply the actual target by epsilon^(3kappa) or rebuild it.

Let A,B be Q37's transferred coefficient bounds for the original principal
and full curl vectors; L,S are the bounds from that SAME root for the real
phase-locked signed principal and curl vectors. Their sufficient orders are

    original principal: 1/2, original curl: 1−kappa,
    signed principal: 3/2, signed curl: 2−kappa.                (1)

These are coefficients of the complete potential-curl construction, with
actual temporal and spatial cutoffs inside the potential. The normal,
pressure, fixed chart conversions, partition, old root and exact signed
formula remain those of Q37. No sourced complex wave has intervened.

## 2. Expand the REAL covariance with all entries and orientations

Real cosine/sine quadrature of a phase-locked wave cancels the principal-
curl cross terms in the averaged covariance, as in Q57. Auxiliary torus
Haar is NOT taken at this step. Denote tensor symmetrization by
sym(u,v)=u v^T+v u^T. In schematic transferred units, including the fixed
quadrature factor 1/2, the difference from the original principal tensor is

    Delta W = (1/2)[epsilon² sym(Avec,Lvec)+epsilon³ Lvec Lvec^T
                +epsilon^(2−2kappa) Bvec Bvec^T
                +epsilon^(3−2kappa) sym(Bvec,Svec)
                +epsilon^(4−2kappa) Svec Svec^T].             (2)

Equation (2) describes physical vectors whose stripped norms are bounded
by A,B,L,S; it does not identify these scalar norm bounds with vectors.
Label sums and real quadrature are taken with the original atlas. Closed
support separation removes distinct labels as in Q37/Q57, not arbitrary
crosses within a label. Both signed squares, all six tensor entries and
both orientations in each symmetrized product remain.

For through-n transferred native raw jets let Tn be Q37's original chart
transfer. The SAME positive coefficient programs already give

    oldcurl_n = Tn 2^n B_n²,
    signed_n = Q37.covariance(n)['complete']
             = Tn [2^(n+1) A_n L_n
                    +2^n(2B_n S_n+L_n²+S_n²)].               (3)

The oldcurl factor is conservative relative to Q57's 1/2 quadrature bound.
Every exponent in (2) is at least 2−2kappa. For 0<epsilon<=1 the sum

    Z_n=oldcurl_n+signed_n                                  (4)

therefore bounds Delta W/epsilon^(2−2kappa) in the RAW through-order norm. With

    principal_n=Tn 2^n A_n², W_n=principal_n+Z_n,             (5)

W_n bounds the full current wave covariance divided by epsilon in that raw
norm. Q84 instead uses the factorial through-order sum in FIVE variables.
If every raw mixed derivative through n is bounded by L, that sum is at most
L sum_(j<=n) 5^j/j! <= exp(5)L <243L, using exp(1)<3. The adapter returns
243 times BOTH (4) and (5), with the same common weight and pole, as Q84's
Z and W slots. Even at n=0 this padding is harmless. It converts norms,
not physical velocities or covariance normalization. Equations
(4)--(5) are coefficient bounds, NOT a triangle estimate between two
order-one physical covariances. The subtraction was expanded first. The
factory reuses Q37's programs rather than square-rooting an updated column
or deleting L². It checks the original and signed root object identity and
the finite declared request order, and returns the original covariance pole.
The signed pole dominates the original curl/product pole on this family.

The final fast and sharp slow changes alter the MEAN velocities, not this
wave covariance. Thus (4)--(5) apply at Q84's FINAL mean state. Correlations
of the fast and slow mean velocities belong to Q84's explicit mean products,
not to a second insertion into W. The signed covariance is added once.

## 3. Common domains and the Q84 slot

At Q84 radial output n take the covariance request N=n+200002. The Q37
object must reserve target derivatives through N+3 for its full pressure
ledger, even though coefficients(n) itself may use fewer. We retain this
conservative bound in the adapter, with all shape, root, potential and
original-amplitude requests of its constructors. A test at low N is not an
evaluation at N=200002, nor is that index the total depth of the source DAG.
The coefficient programs are finite positive polynomials at each fixed N;
the actual uniform inequalities are still premises of the canonical maps.

W_n and Z_n carry the same declared NEW covariance pole. Q84's T and all
mean/base inputs must be evaluated on a common archived weighted domain
that dominates it, or be explicitly unweighted on the same positive hull
with the inherited fixed-profile costs. The adapter does not relabel a pole
or declare monotonicity between unrelated weighted norms. The full T0,
Pi1/Pi2 memory, moving pressure moment, actual axial gd−A1 gd velocity and
all corrections remain in the mean-update programs. This velocity-covariance
identity does not assert that the pressure change is zero anywhere.

## 4. What is supplied and what remains

Equations (4)--(5) remove the two independent W and Z coefficient inputs in
Q84 relative to the Q37/Q51 canonical physical hypotheses. A caller providing
only the historical alpha target bound has not supplied the rate-2 premise:
without it, the principal signed change is only O_(2−3kappa), which is
insufficient to infer Q84's stated Z normalization. The adapter refuses that
rate. This is a corrected interface, not evidence of failure of the actual
stronger target. The old conservative ledger remains valid and frozen.

Next supply Q84's FINAL mean-field and angular-row coefficients from the
corrected full-pressure Q47/Q39/Q41 update graph, including the original
initial axial tail and the whole reserved-potential support. Keep the
first-real domain distinct from any later sourced state. Unconditional
pulse bounds, complete residual improvement, infinite realization and
blow-up remain unproved.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. Exact tests expand the real tensor,
detect deleted squares/orientations, preserve the actual root and pole,
and reject the old-rate shortcut. Rational fixtures are not continuum
suprema or proof of the selected enormous field inequalities.
