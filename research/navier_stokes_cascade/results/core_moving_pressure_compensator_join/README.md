# Moving pressure compensator and new-pole moment substitution

Parent: `ab0c64127a9aef62cab16f7ba7d88c97e8b31e33` (24 October 2026).
The remote branch matched this parent before editing. Its Actions run
40084453931 was still running with no reported failed step; the parent's
new signed-covariance step had passed. This checkpoint
replaces THREE independent inputs to the preceding complete signed-pressure
join: the moving rho norm, whole c_rho norm, and radial moment norm Kmom.
The compensator identity and coordinate theorem already exist in
`core_reserved_geometry_recheck`; the new result is their finite mixed
R,Z,T profile composition and substitution at the enlarged covariance pole.
Historical evidence, actual profiles, width, band and pressure gauge are
unchanged. The actual sigma_e family is a DIFFERENT profile: none of its
missing derivative bounds is inferred from this rho calculation.

## 1. Fixed-profile data, not a new moving-family hypothesis

The RETAINED pressure compensator is

    rho(R,Z,T)=u(Z,T) rhohat(R u(Z,T)), u=s^(-1/2),
    integral rhohat(x) dx=1, supp rhohat subset [l,u0].       (1)

Here 0<l<u0 and rhohat is the same fixed smooth compact profile, with
smooth zero extension. It need not be nonnegative. Its finite fixed
seminorms H_k=sup |rhohat^(k)|/k! remain input data, uniform simply because
the fixed profile is independent of the physical band. The API requires
them through every requested order. It does not certify their values,
the identification of that profile, or its support/mass by metadata.

Use the original implicit coordinate s-Z^2 s^(2h)=T, 1/4<=s<=4,
|eta|<=1, 0<h<1/100. The previous coordinate proof bounds the TOTAL
degree-j absolute normalized Taylor coefficients of log s by

    L_j=(1/3) sum_(i=0)^j 32^i 64^(j-i), j>=1.

Define U_0=2 and

    U_n=(1/(2n)) sum_(j=1)^n j L_j U_(n-j).                 (2)

This is 2 exp((1/2) sum L_j t^j). It bounds the dilation u's degree-n
jets. There is NO additional x-log coordinate here: adding the log-chart
radial increment would bound a different map and obscure the actual
R-derivative. On the original whole radial norm hull [a,b], the inner
coordinate R u has positive increment series

    C(t)=sum_(j>=1) (b U_j+U_(j-1)) t^j.

The second term is precisely the mixed product of the R increment with
u; it includes mixed radial/slow derivatives. Finite Taylor composition
and the absolute-coefficient product inequality give

    ||rho||_(n,unweighted)
       <= sum_(i+j<=n) U_i sum_(k=0)^j H_k [t^j] C(t)^k.  (3)

The norm sums |partial_R^r partial_Z^z partial_T^t rho|/(r!z!t!)
through total order n. Smooth finite jets suffice for (3); rhohat need
not be analytic and no infinite Taylor expansion is asserted. All
suprema are uniform on the retained real chart, including one-sided
physical endpoint limits supplied by the coordinate theorem.

## 2. Weight and whole pressure moment on their proper hulls

The ACTUAL moving support is [sqrt(s)l,sqrt(s)u0]. Require this entire
support to fit the same norm hull [a,b]. The established reserved
geometry computes E from the unchanged annulus and width and gives
zeta>=exp(-E) on that moving support. Thus exp(E) times (3) bounds
the weighted rho jets there, with smooth zero extension elsewhere.
This also bounds rho in every larger nonnegative edge-pole class,
since delta<=1. No rectangular-product weight floor is asserted for
the integration union [l/2,2u0]. The original annulus and support
identifications remain premises, rather than a new choice of rho.

The exact change of variable in (1) gives c_rho=s*c_hat, where
c_hat=(1/2) integral x^2 rhohat(x) dx. A sufficient fixed coefficient is

    cstar=(u0-l) u0^2 H_0/2,
    ||c_rho||_n <= cstar (4+sum_(1<=i+j<=n)32^i64^j).      (4)

The estimate uses |rhohat|, not the false implication that signed unit
mass gives L1 norm one. It supplies the WHOLE c_rho norm, including its
nonzero Z derivatives. Inserting (4) at n+1 preserves
epsilon partial_Z(Jz+c_rho P), including c_rho,Z P and both old and new
scalar defects; it does not replace total current pressure by its change.

## 3. A direct moment bound at every requested new pole

For a complete smooth zero-extended input supported in the same fixed
radial integration hull [a,b], suppose its through-order coefficient
bound is zeta delta^(-p) H. The existing edge absorption inequality is

    sqrt(zeta) delta^(-p) <= E_p,
    E_p=max(1,N! (2/w^2)^N), N=ceil(p/2).                 (5)

Consequently its UNWEIGHTED through norm is <=E_p H. The fixed-limit
integral and normalized Haar commute with Z,T differentiation. For
each radial weight R^j, j=-1,0,1,2, a single sufficient bound is

    M_p=(b-a)(a^-1+1+b+b^2) E_p.                         (6)

Indeed each weight is bounded by the corresponding summand on [a,b].
The scalar output has no radial derivatives. Summing the normalized
slow jets cannot exceed the incoming full R,Z,T through norm, so (6)
has no derivative loss or omitted factorial. Haar averaging is a
contraction. This covers the pressure mass and all radial weights
in the historical wave_change functional after its already-recorded
integration by parts; the smooth endpoint terms vanish. It is NOT a
bound for a weighted pressure primitive or for its radial derivatives.

Lower-pole covariance inputs embed in the largest requested pole since
delta<=1. Thus (6) can be used as Kmom at the NEW pole requested by
RealCovariance.pressure, rather than reusing an old-domain callback.
The full wave support must fit [a,b]; the smaller reserved rho support
does not bound wave tails. E_p is an exact finite factorial recipe;
the exp(E) floor stays symbolic. No practical high-order materialization
or uniformity as p tends to infinity is claimed.

## 4. Substitute without changing the current state

MovingCompensator.operators requests only order,pole,Kp,primitive,bump
from the remaining map supplier. It checks their coverage and inserts
(3)--(6) into the preceding all-six-entry pressure/scalar recipe.
The new covariance pole, signed squares, axial covariance, old scalar
additions and whole moving pressure memory are retained unchanged.
This is a conditional analytical reduction of the lower input list;
it is not a fine residual gain and does not change any epsilon exponent.

The precise remaining map obligations at the requested order and pole
are the ACTUAL weighted pressure primitive inequality

    ||I0(f-rho integral Haar f)||_(n,weighted,p)
       <= Kp(n,p) ||f-rho integral Haar f||_(n,weighted,p),

with its same zero-mass support and oriented edges, plus the original
target primitive and actual sigma_e projector seminorms. The metadata
checks are necessary domain checks, not proofs of these inequalities.
The existing two-edge symbolic kernel route can be used without changing
physical profiles; (6) alone cannot supply that primitive. No independent
uniform moving rho norm or radial moment callback is now needed beyond
the stated fixed-profile/coordinate/support data.

Next justify those remaining pressure/projector inputs and recompute
the post-signed fast/slow means from the CURRENT rows and defects.
Actual uniform generated-field/profile enclosures, coefficient/support
compatibility, complete viscosity/cutoff/interactions and the full
residual remain open. Unconditional pulse bounds and blow-up are unproved.

## 5. Evidence

Exact controls independently expand the positive dilation and composition
series, check local mixed R,Z,T jets of the actual implicit coordinate,
and integrate endpoint-flat signed polynomial diagnostic profiles under
three dilations. The polynomial profile tests are marked manufactured
diagnostics, not substitutions for rhohat. Negative controls reject
missing jets, impossible unit-mass height/support data, changed width,
out-of-hull support, and stale remaining operator domains. Fresh replays
bind the historical frozen sources. Passing tests is not a continuum proof.

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
