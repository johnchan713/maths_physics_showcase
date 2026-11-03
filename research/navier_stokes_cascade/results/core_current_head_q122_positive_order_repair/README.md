# Q122: all-order positive-coefficient repair bound and current tuple join

Parent: `d8455586d7de4461534b6b5989c7f6dd87e506bb` (Q121).
On 3 November 2026 the freshly fetched remote agreed with that parent.
Q121 research CI 41614701355 and field CI 41614701499 were in progress;
Q119 field CI 41601360763 had passed. Current notes, frozen Q121 evidence,
the background jet foundation, the whole-support first-coefficient repair
and the all-order modulation/leading-root map were inspected first.
The core_growth_activation result remains local leading core/first collar;
its width C_setup^-20 and all historical evidence remain unchanged.

Q121 needs finalized F,U jets at arbitrary fixed order. This checkpoint
provides the positive-order LINEAR repair part of those bounds, using the
actual reserved-patch normalizer. The compact leading-profile modulation
and its nonlinear accepted root already have an all-order MAP theorem in
modulation_repair_all_order; they are not repeated or replaced here.
Actual incoming moments, cut-seed jets and bump jets remain inputs.

## 1. Original five rows, not a pressure reset

For each positive coefficient n>=1, the background foundation's (7) gives
five linear variations in ordinary moment order:

    (int delta U, int sqrt(2X)delta E, int E0 delta E/X,
     int sqrt(2X)E0 delta U, -int E0 delta E).

Keep their full incoming discrepancies d1,...,d5 as functions of eta.
They include all lower-order products, full pressure-source integrals through
Xv, radial memory and viscosity, including the axis boundary contribution.
A discrepancy computed only to Xcut, I3 or I4 is not the input in this
checkpoint. Zero local velocity does not set any incoming moment to zero.
No current-order quadratic term occurs: it has formal order 2n>n.

The original patch has E0=e(eta) exp((-.5-lambda_patch)y), U0=0,
e(eta)=e_star/(1+eta^2), e_star>0, 0<lambda_patch<=.01. This precise e
identity is stated in global_first_coefficient_c3 (9). For a generic smooth
positive e the polynomial calculation below is unavailable; the API rejects
that different identity. An identity label itself does not prove that a
physical patch satisfies it.

Put rho=rho3, K_g=(sqrt(2)rho^(3/2))^-1, P=1/e=(1+eta^2)/e_star.
In the background foundation's normalized/reordered coordinates the target is

    z=(d1/rho, (d1/rho-K_g P d4)/lambda_patch,
       K_g d2, -P d5/rho, P d3).                    (1)

The same matrix B is independent of eta and n; c=-B^-1 z. The foundation
proves ||B^-1||_infinity<=1000 for the actual ordered unit-mass bump supports
and the entire positive exponent range stated above. This conservative
constant is used here; no improvement over the earlier first-order 800
bound is claimed. There is no small-root or contraction premise for this
positive-order solve. In contrast the n=0 leading repair is quadratic, so
this API rejects n=0 rather than applying (1) to that root.

## 2. Exact normalization has only three derivative rows

On |eta|<=1, raw Leibniz differentiation gives exactly

    (P f)^(k)=[(1+eta^2)f^(k)+2k eta f^(k-1)
                     +k(k-1)f^(k-2)]/e_star,        (2)

with negative-order terms absent. Thus, for supplied raw bounds a_j,

    ||(P f)^(k)|| <=[2a_k+2k a_(k-1)
                          +k(k-1)a_(k-2)]/e_star.  (3)

There are no unknown derivatives of P above two. If every discrepancy has
raw derivatives through m bounded by its supplied D_i, put D=max_i D_i and

    p_m=(m^2+m+2)/e_star,
    G=max(1,rho^-2).

Equation (3) bounds every needed P d_i jet by p_m D. All geometric factors
rho^-1 and K_g are <=G: for rho>=1 both are <=1; for rho<=1,
sqrt(2)rho^(3/2)>=rho^2. Consequently (1) and lambda_patch<=1 imply

    max_(0<=k<=m) ||z^(k)||_infinity
       <= G D (1+p_m)/lambda_patch,
    max_(0<=k<=m) ||c^(k)||_infinity
       <= 1000 G D (1+p_m)/lambda_patch =: C_m.    (4)

This is an analytical implication for every fixed m. Its normalizer cost
is quadratic in raw order, rather than an artificial factorial loss from
an unspecified reciprocal estimate. In normalized factorial jets (2) has
weights (2,2,1)/e_star. Their sum is 5/e_star at every order >=2, agreeing
with the historical complete C_eta^3 norm W3. A factorial C_eta bound on
incoming discrepancies must still be converted before insertion as a RAW
D_i; that conversion is not erased by the polynomial identity.

## 3. Finalized profile and Q121 interface

The unchanged corrections are sums of two U bumps and three E bumps with
angular coefficients c_i. Their fixed supports are mutually disjoint and
flat at the endpoints. Let B_U,r bound radial derivatives through r of
each U bump, and B_F,r those of each E bump DIVIDED BY sqrt(2X), on the
required positive annulus. These bounds concern the actual fixed bumps;
no Borel, variance or compact-X cutoff bound is substituted. In particular
B_F,r retains every inverse-radius derivative.

If Fcut_m is the raw mixed norm of the same tentatively cut F through m,
and Ucut_(m+1) is the corresponding U norm on [0,Xhi], then disjointness
and separation of X/eta variables give

    Fraw_m <= Fcut_m + B_F,m C_m,
    Uraw_(m+1) <= Ucut_(m+1) + B_U,(m+1) C_(m+1).  (5)

No binomial coefficient enters a mixed derivative of c(eta) beta(X).
Disjointness makes the sum bound a maximum at each radius. A common
conservative input D through m+1 permits C_(m+1) on both right sides.
The function finalized_tuple uses this latter convention and passes (5)
to the unchanged Q121 closed tuple recipe. Thus its explicit remaining
inputs are the cut-seed norms, ALL five full incoming discrepancy norms
through m+1, the same patch parameters and the actual radial bump norms.
No high angular derivative of e remains an independent input on this patch.

Q121 then supplies its m! tau^-m tuple recipe and a binary upper budget.
Q120's existing-scale inequality may use that budget if the actual cutoff
identity and retained log2(a_n) lower bounds are certified. The large-order
100005 diagnostic requests discrepancy order 100006, retains factorials as
recipes, and evaluates only binary budgets. Its constants are manufactured;
it does not certify an actual source, profile, scale, or cone bound.

## 4. Validation and outstanding work

Independent symbolic differentiation checks (2) beyond the polynomial
profile's degree, including eta endpoints. A manufactured five-row constant
linear system checks target ordering, signs and raw derivative domination;
it is not a sampled determinant or a substitute for the inherited continuum
inverse proof. The tests also check factorial/raw conventions, zero-order
scope and missing-row rejection, and the Q121/high-order recipe join.
The audit recursively replays Q121 and all its declared source/evidence
hashes. Fresh reproduction must match the frozen record byte-for-byte.
Passing these checks is not a continuum proof.

Next derive the actual whole-interval cut-seed, fixed bump and full incoming
moment-discrepancy bounds at the required high orders, with complete
pressure memory and existing domain coverage. Bind the patch identity and
support conditions to the actual coefficients. Then verify the retained
scale inequality and carry actual constants through the current high-pole
source and compact-modulation/exact-repair cone join. Global pressure,
moment and reserved-mean support identities, complete-wave covariance mean
updates including axial velocity and pressure, uniform actual-pulse bounds
at every physical scale, and the full viscous/cutoff/interaction residual
remain open. No field, root, scale, support, pressure datum or mean is
reselected. No unconditional pulse bound, residual improvement or blow-up
conclusion follows.
