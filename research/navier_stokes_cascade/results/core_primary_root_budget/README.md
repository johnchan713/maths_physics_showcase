# Quantitative canonical root, edge and signed-linearization budgets

Parent: `3a1ac5299d888033d079ca3a2fd71173691c1c7a`.
The parent derives aggregate bounds from uniform per-label coefficient
polynomials. This checkpoint replaces their implicit root and edge constants
by finite positive-polynomial recurrences in the normalized canonical
column, weighted target and homogeneous-pulse jets. It supplies the SAME
canonical roots and original-root signed map, rather than rebuilding a
different pulse. The root regularity theorem itself is historical; its
quantitative substitution into the aggregate/initial-scalar chain is new.

This remains conditional on the canonical column/ODE estimates, exact
core/background and strict value/cone premises. Those inputs are identified
below, not independently certified here. Fixed radius r0=1/6057072,
width t1=C_setup^-20, carriers, supports and the accepted common tail are
unchanged. No selected enormous physical profile is numerically evaluated.

## 1. Use factorial-normalized jets before doing root arithmetic

In one frozen native coefficient chart let

    J_n(F) = sum_{|I|=n} ||partial^I F||/I!,

using absolute value, vector maximum or induced matrix infinity norm.
At a point, products obey J_n(FG)<=sum_{j=0}^n J_j(F)J_{n-j}(G).
The scalar envelope below is therefore valid in any fixed finite chart
dimension; it does not identify ordinary derivatives with Taylor
coefficients. Every raw derivative through n is bounded by
n! sum_{j=0}^n J_j. The fixed R-to-U conversion and every chart coefficient
are included in the input jets. No integer index is differentiated.

Write H=s0 Hbar for the CANONICAL principal column matrix, with
s0=r0^2/sqrt(L) and L=2r0/ci. L and s0 are constant under continuous
coefficient differentiation in this frozen band. Let

    Hbar b = That,  That=T/zeta,  y=(zeta/s0)b,
    zeta=exp(-t1^2/ya^2-4/yb^2), delta=min(1,ya,yb).

Here T is the original leading stress target, not the full corrected stress.
The canonical cone bounds give, componentwise,

    b_sigma >= beta := gamma*vstar/(8*jmax) > 0,
    ||Hbar^-1||_infinity <= K0 := 6/(jmin*min(Amin,u)).       (1)

Pad K0 upward if necessary. The actual weighted-stress estimates and smooth
coordinate composition give positive polynomials H_n(S),T_n(S) such that

    J_n(Hbar)<=H_n(S),
    J_n(That)<=T_n(S) delta^(-3n-3).                       (2)

The second bound follows by differentiating T/zeta, including derivatives
of zeta^-1: derivative order j of that reciprocal costs at most delta^-3j,
and the weighted T bound costs delta^(-3(n-j)-3). This is a polynomial
edge loss, not boundedness of That at the outer edge. Indeed That can
grow like yb^-3. These estimates retain the exact inner and outer flat
constants t1^2 and 4 and all fixed coordinate/width constants.

Define positive polynomials, with the same K0 as (1),

    K_n = K0 sum_{j=1}^n H_j K_{n-j},             n>=1,
    B_n = sum_{j=0}^n K_j T_{n-j}.                          (3)

Differentiating Hbar^-1 Hbar=I, then Hbar b=That, proves
J_n(Hbar^-1)<=K_n and J_n(b)<=B_n delta^(-3n-3).
No root is used to prove these matrix/column estimates.

## 2. Explicit reciprocal-root and flat-factor constants

Choose a fixed R>=beta^-1/2; this is a VALUE floor, not a new higher-jet
acceptance condition. For a formal variable z put Bplus(z)=sum_{j>=1}B_j z^j
and c_k=|binom(-1/2,k)|. At each requested finite order define

    R_0=R,
    R_n=R sum_{k=1}^n c_k beta^-k [z^n]Bplus(z)^k, n>=1,
    Q_n=sum_{j=0}^n B_j R_{n-j}.                            (4)

The finite Taylor formula for x^-1/2 at b_sigma>=beta gives
J_n(b_sigma^-1/2)<=R_n delta^-6n: each product with k factors
has edge exponent 3n+3k<=6n. Multiplying by b_sigma gives
J_n(sqrt(b_sigma))<=Q_n delta^(-6n-3). The deliberately padded
zero-order bound uses B0 R, not an unjustified bounded normalized root.
The formulas use Taylor jets with remainder; no infinite power series
convergence at an edge is assumed.

Let ell_n(S) bound delta^(3n) J_n(log sqrt(zeta)) for n>=1, including
the smooth ya,yb chart maps. Such bounds follow from differentiating
-t1^2/(2ya^2)-2/yb^2: each term costs at most delta^(-n-2)<=delta^-3n.
Define E0=1 and

    E_n=(1/n) sum_{j=1}^n j*ell_j E_{n-j}.                  (5)

Then J_n(sqrt(zeta))<=sqrt(zeta) E_n delta^-3n. The opposite
reciprocal exponential has the same differentiated polynomial losses,
but is never asserted bounded or smoothly extendible at its own edge.

Let W_n(S) bound J_n(chi*psi*t_sigma)/P for the COMPLETE cutoff canonical
homogeneous vector before multiplication by its root. These are the
retained differentiated ODE, Gaussian and fixed-cutoff estimates, and
include all parameter, auxiliary and pulse derivatives. Equations (3)--(5)
give the coefficient of the sqrt(epsilon) primary velocity:

    A_n = Cmass*S sum_{i+j+k=n} E_i Q_j W_k,
    J_n(A)<=sqrt(zeta) P A_n delta^(-6n-3).                (6)

Since L<2r0*Tg*S, sqrt(s0)^-1 < (2Tg/r0^3)^(1/4) S^(1/4).
The rational choice Cmass=2/r0 suffices for r0<=1,Tg<6, and
S^(1/4)<=S for S>=1. Thus (6) keeps an explicit fixed physical Haar/radius
cost. It does not declare the estimate uniform as r0 or t1 varies.

## 3. Complete curl and the aggregate/pressure-moment chain

Express the SAME normalized physical curl as a first-order coefficient
operator on A:

    epsilon^kappa B = D0 A + sum_mu Dmu partial_mu A.

The matrices include n cross /|n|^2, cylindrical connections, radial
graph derivatives and all slow/auxiliary terms. The normal floor is |n|,
not an angular component. Let D_j(S) bound the SUM of J_j norms of these
matrices. This is a retained normal/background/graph input, not permission
to omit any operator term. The fixed-scale floor M epsilon^kappa<=1
and epsilon<=1 allow positive-polynomial bounds with no remaining
epsilon loss. Product differentiation yields

    C_n=sum_{j=0}^n D_j[A_{n-j}+(n-j+1)A_{n-j+1}],
    J_n(epsilon^kappa B)
       <=sqrt(zeta) P C_n delta^(-6n-9).                 (7)

Indeed sum_mu J_m(partial_mu A)=(m+1)J_{m+1}(A) for scalar
coefficient absolute sums, and the corresponding inequality suffices for
vector norms. Derivatives of the Dmu cannot be dropped.

For the parent's raw through-order-n input use

    piA_n=2 n! Eedge(6n+3) sum_{j=0}^n A_j,
    piB_n=2 n! Eedge(6n+9) sum_{j=0}^n C_j,                (8)

where the fixed factor 2>=sqrt(3) also covers the conversion from the
three-component maximum norm to its Euclidean norm, and Eedge is the
parent's fixed-width series constant. This unweights
sqrt(zeta) and P<=1 only AFTER all roots and curls are formed. Substitute
(8) into the actual-atlas formulas a_n=8*188^n*piA_n((25/16)S),
b_n=8*188^n*piB_n((25/16)S), with no second edge-unweighting factor,
and Omega_n=2^(n-1)(a_n^2+b_n^2). The matched coarsest carrier/fast-time
convention, support separation and normalized Haar measure remain.

This removes independent root/reciprocal-root and complete coefficient
polynomials as inputs RELATIVE TO H_j,T_j,W_j,D_j and the value floors.
It does not remove these lower-level canonical ODE/column/normal premises.
Scalar output n still uses complete coefficients through n+10; the WHOLE
moving axial-pressure moment uses n+11. Resolving B into its amplitude
therefore needs H,T,W through n+12 and D through n+11 for that moment.
Pressure reconstruction, its moving cutoff, Pi1/Pi2 and the scalar moment
are inherited in full from that chain, not inferred from the principal
root equation. No complete pointwise residual follows from this join.

## 4. Same-root signed update and the next source gate

For a prescribed real, auxiliary-independent signed target
Xi=epsilon^alpha zeta xi with J_n(xi)<=X_n delta^(-3n-d), d>=0,
use EXACTLY

    dvec=H^-1(Xi/epsilon), delta a=dvec/(2a),
    a=sqrt(zeta/s0)*sqrt(b).

Let V_n=sum K_j X_{n-j}. Then

    Z_n=(Cmass*S/2) sum_{i+j+k+l=n} E_i V_j R_k W_l       (9)

bounds the signed principal velocity coefficient by
epsilon^(alpha-1/2) sqrt(zeta) P Z_n delta^(-6n-d).
The same curl rule adds one input jet, edge power at most 6n+d+6,
and gives the remainder order alpha-kappa. This is a mapping for an
INPUT xi satisfying the full weighted bound and support, not verification
of the actual next signed target or its complete nonlinear sources.

Exactly 2 epsilon H diag(a) delta a=Xi. The signed square
epsilon H (delta a)^2 generally survives. Replacing a by an updated
positive root, omitting its differentiated reciprocal, or treating the
whole covariance as the principal cross identity changes the map.
No positivity or smallness of the updated signed amplitudes is used.

Next discharge the remaining actual H,T,W,D/core premises and verify
that the CURRENT signed target and then the sourced residual satisfy
these displayed weighted norms and support conditions. Full axial means,
pressure, A0/A1, viscosity, cutoffs and all interactions remain in the
complete residual obligation. Infinite realization and blow-up are unproved.

## 5. Evidence and reproduction

Exact noncommuting matrix Taylor controls, multivariate reciprocal-root
and original-root signed controls test the formulas independently of
the positive-polynomial engine. Negative controls retain inverse/root
derivatives, complete curl coefficients and the signed square. Degree
ledgers and substitution into the whole initial scalar/pressure chain
retain its n+12 underlying-input cost. These are manufactured finite jets,
not continuum suprema or actual giant-parameter evaluations. Sections
1--4 are the conditional analytical derivation.

```sh
python3 research/navier_stokes_cascade/results/core_primary_root_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_primary_root_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_primary_root_budget/evidence.json \
  --output /tmp/core_primary_root_budget_fresh.json
```
