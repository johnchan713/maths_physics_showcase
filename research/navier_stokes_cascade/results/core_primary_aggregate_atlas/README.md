# Common-atlas primary aggregate bounds and matching band normalization

Parent: `1e5bc69040c0aaa0d740819dc3eaccb56aff5833`.
The parent retained aggregate primary coefficient norms as inputs. This
checkpoint derives them from the retained UNIFORM PER-LABEL canonical
jet estimates, using the actual frozen atlas, explicit edge constants and
a consistent common band. It does not re-prove the primary ODE/root/core
lemmas. Qualitative support separation and covering bounds were proved in
`auxiliary_geometry_audit`; the new result is their quantitative substitution
into the initial scalar chain and the normalization caveat below.

All original profiles, radius r0=1/6057072, tilt, width C_setup^-20,
carriers, absolute fields and common tail remain unchanged. Selecting a
common coordinate representation does not alter any of these data.

## 1. Edge losses have explicit finite unweighted constants

For the same t1=C_setup^-20 in (0,1], put

    zeta=exp(-t1^2/ya^2-4/yb^2), delta=min(1,ya,yb),
    alpha=t1^2/2, N=ceil(b/2),
    E_b=max(1,N! alpha^-N), b>=0.                        (1)

Then sqrt(zeta)delta^-b<=E_b on the entire active shell. If delta<1,
the exponent is at least alpha/delta^2. With x=delta^-2>=1,
x^(b/2)<=x^N and exp(alpha x)>=(alpha x)^N/N! prove the bound.
For delta=1 simply use sqrt(zeta)<=1. The profile extends smoothly by
zero; no derivative of the quotient by zeta is being bounded here.

Thus a per-label native estimate C_n S_i^p sqrt(zeta)delta^-b P,
with P<=1 and every cutoff already included, gives the unweighted bound
C_n E_b S_i^max(0,ceil(p)). Constants can be enormous. For integer b
and N>=1, log E_b=sum_(j=1)^N log j+N log 2-2N log t1;
the actual width need not be materialized numerically. There is no
uniformity in b or t1->0 and no new small-q threshold assertion.

## 2. Fix the coarsest band on the whole common neighborhood

Take the SAME band set used by the constructive common-torus atlas on
a small (z,t) neighborhood, independent of the radial integration point.
Its band span is at most four. Set ell*=min ell_i, Q*=2^-ell*,
epsilon*=Q*^h and S*=ell*^2. The covering index is nondecreasing on
ell>=16, so i*=i(ell*) is the atlas's minimum index. In H=J^i*Y,

    Y_i=J^Delta_i H, 0<=Delta_i<=3,
    J=((3,1),(1,5)), ||J^Delta||_1=||J^Delta||_infinity<=188.

The coordinate neighborhood is chosen inside the existing enlarged chart
1/4<=q/Q*<=4. At its center only finitely many closed dyadic supports
meet; bands outside that closed set have a positive local gap and can
be excluded by shrinking the neighborhood. The closed supports have
1/2<=q/Q_i<=2, leaving the stated enlargement. This is a local coordinate
choice, not a new accepted small-q tail or a change of endpoint labels.

No index is differentiated. Relative native coordinates are linear:
U_i=(Q*/Qi)^(d_r/2)U*, Z_i=(Q*/Qi)^D Z*, T_i=(Q*/Qi)T*.
Because Q*/Qi<=16, 0<rho_g<1 and 0<h<1/100,
d_r/2<101/100<5/4. Their derivative multipliers are bounded by
32,4,16, respectively. The auxiliary chain rule therefore gives the
common raw coefficient jet factor 188^n at order n. Higher derivatives
of these coordinate changes vanish. The constants include the actual
fixed conversion between R and U=R^d_r already used by the parent.

Keep the complete velocity in physical units before converting:

    Qi^-A sqrt(epsilon_i)[Ai cos-ki^-1 Bi sin]
      =Q*^-A sqrt(epsilon*)[Atilde_i cos-ki^-1 Btilde_i sin].

Both coefficients acquire r^(A-h/2), r=Q*/Qi. For the normalized curl
coefficient epsilon*^kappa Btilde_i the factor is r^(A-h/2+h kappa).
Both are bounded by 8 because their exponents lie in (0,3/4).
This also agrees with the reference normalized curl: n*=sqrt(r)ni
and curl_* =sqrt(r)curl_i, so the normal reciprocal and curl factors
cancel. Every cylindrical connection remains in Bi.

Moreover S_i<=(25/16)S*. Most importantly, coarsest normalization gives

    ki^-2<=epsilon_i<=epsilon*,
    c*=Tg^i* Q*^(1+h), (Tg S*)^-1<c*<=S*^-1.           (2)

Hence the parent's carrier and fast-time floor bounds hold together
in this same representation. The band set and common index stay fixed
along the entire radial integral; absolute-variable gluing handles chart
changes without new labels or derivative-dependent band choices.

## 3. Actual closed separation removes the aggregate counting assumption

The atlas has at most 162 slow labels and 2744 deck copies, but its
1,250-color CLOSED enlarged support separation permits at most ONE
actual wave label at a point of the extended common domain. Smooth
zero extension keeps every derivative in that support, including at
boundaries. Different copies of one band rectangle are also disjoint.
Thus the label-l2 derivative norm is bounded by the largest per-label
norm, without a square root of either label or lift count.

Let piA_n(S),piB_n(S) be the uniform positive polynomial native bounds
for Ai and epsilon_i^kappa Bi obtained from the retained per-label
canonical estimates and (1). All ODE, inverse/root, radius and cutoff
constants remain in these polynomials. They are bounds for complete
coefficients, not positive lower bounds for oscillating velocities.
With K=25/16, the aggregate polynomials required by the parent are

    a_n(S)=8*188^n piA_n(KS),
    b_n(S)=8*188^n piB_n(KS),
    Omega_n(S)=2^(n-1)[a_n(S)^2+b_n(S)^2].               (3)

The aggregate inequality is now derived from the actual atlas and the
native per-label hypotheses, rather than imposed independently. It is
pointwise in the common torus BEFORE Haar averaging. Normalized Haar
measure is preserved by the covering, so no factor 2744 is added to
moments. Whole-hull radial integration can likewise use the parent M_j
constants after this pointwise bound. There is no need to sum per-label
suprema along that integral. Closed supports and locally finite sums
justify this reasoning; different angular frequencies alone do not.

Substituting (3) in the parent's SeedBudget supplies its explicit initial
S_2 scalar and epsilon^3 whole-pressure-moment budgets. The parent's
n+10/n+11 complete-coefficient jet counts remain; the derivative and
width constants in (1),(3) are retained at each of those orders.

## 4. Normalization caveat and exact failure controls

Native bounds must not be silently reused with a finer reference Q and
the minimum covering index. Even if ell_ref-ell_i=1, then
epsilon_ref=2^(-h)epsilon_i. Both dyadic supports can be active together.
The literal bound ki^-2<=epsilon_ref
can fail: ki<=epsilon_i^-1/2+1 gives

    ki^-2>=epsilon_i/(1+sqrt(epsilon_i))^2>epsilon_ref

for all sufficiently small epsilon_i. This is a rounding/floor
normalization issue, not an obstruction to the physical construction.
For example h=1/200, ell_i=40000 has ki=2^100 exactly, while the
adjacent finer reference has strictly smaller epsilon. This example is not an
evaluation of the selected physical parameters or accepted first band.

An arbitrary reference among these bands is still usable with explicit
constants: epsilon_i/epsilon_ref<=16^h<2 and
c_common^-1<Tg^4 S_ref. The aggregate covariance bound then uses
2^(n-1)(a_n^2+2 b_n^2), and the fast-map floor constant uses Tg^4,
not Tg. Since Tg<6, rational constants 2 and 1296 suffice. Equation (2)
avoids both extra costs by the matched coarsest choice. Earlier native
floor bounds remain correct in their matched charts; their unrestricted
mixed-reference use is what these controls reject. No historical evidence
is rewritten and no epsilon exponent or support width is changed.

Remaining obligations are the retained uniform per-label canonical
ODE/root/edge estimates relative to the exact core/background premises,
the actual signed/sourced seed and source supports, and the complete
pressure, viscosity, cutoff and nonlinear pointwise residual. The aggregate
join does not remove those dependencies or prove blow-up.

## 5. Evidence

The exact atlas separation, covering and counting certificates are bound
and replayed. Rational finite coefficient expansions check the common
derivative multiplier, edge-series constants and normalization exponents.
Deck-copy measure controls distinguish normalized Haar average from an
unnormalized copy sum. Positive polynomial substitutions feed (3) through
the parent's initial scalar chain in both reference conventions. These
finite controls are not continuum proofs; sections 1--4 give the scoped
analytical argument under the retained native pulse hypotheses.

```sh
python3 research/navier_stokes_cascade/results/core_primary_aggregate_atlas/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_primary_aggregate_atlas/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_primary_aggregate_atlas/evidence.json \
  --output /tmp/core_primary_aggregate_atlas_fresh.json
```
