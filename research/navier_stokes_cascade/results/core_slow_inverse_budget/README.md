# Explicit fixed-profile slow inverse in the coupled scalar budget

Parent: `b8da70e057ca500f508ba8f6f19ead1499a0e253`.
Its coupled scalar estimate retained abstract slow inverse majorants
L_u,L_g,L_B. This checkpoint replaces those by explicit finite-jet
expressions for the SAME five-map, including its distinct physical row
normalizations and the radial potential. No profile, positive tilt,
width C_setup^-20, reserved hull or accepted tail is selected again.

The five-map's existence, sharp fixed 1/lambda axial conditioning and
qualitative all-fixed-order mapping were already proved conditionally in
`compact_mean_foundation`. They are not new results here. The advancement
is the explicit response-profile representation and quantitative substitution
into the current coupled recurrence, with all normalization and jet costs.
The retained background floor and smoothness remain hypotheses. Canonical
seed, sourced/fast inverse and full residual obligations remain separate.

## 1. Freeze three response profiles, not updated target-dependent roots

Let x=r/sqrt(q), g=exp(d)>1 and v=exp(-2d lambda), lambda>0 fixed.
For the originally fixed geometric bumps eta_j(x)=g^-j eta_0(x/g^j),
put mu_p=int x^p eta_0(x) dx>0. Let

    Atheta(a)=diag(mu_2,2a mu_(-2-2lambda),-a mu_(-2lambda))
                   Vandermonde(g^2,v/g^2,v),
    Az(a)=diag(mu_1,a mu_(1-2lambda)) Vandermonde(g,gv).

Rows evaluate powers 0,...,n-1 at the displayed nodes. Define ONCE

    bP=Atheta(1)^-1 (0,-1,0),
    bz=Atheta(1)^-1 (0,0,-1), bg=Az(1)^-1 (0,-1),
    FP=sum_(j=0)^2 bP_j eta_j,
    Fz=sum_(j=0)^2 bz_j eta_j, Fg=sum_(j=0)^1 bg_j eta_j,
    Psi(x)=x^-1 int_left^x y Fg(y)dy.                     (1)

The matrices are invertible for the retained lambda>0. Their constants
keep the actual moments, bump separation and tilt. In particular bg has
the sharp factor (1-v)^-1; it is not uniformly bounded as lambda->0.
The angular block has no such singularity at zero. The construction is
linear for arbitrary finite signed targets: no target smallness test.

The zero first rows give int x^2 FP=int x^2 Fz=int x Fg=0.
Consequently Psi is smooth with compact support in the entire convex
hull of the axial bumps. It can fill gaps between them. The velocity and
potential hulls remain within the SAME reserved I4. Also

    Psi'+Psi/x=Fg.                                      (2)

## 2. The three different scale factors survive in the actual response

Write s=q/Q in [1/4,4], epsilon=Q^h, A=1/2+h, and use Q-normalized
CURRENT scalar defects P,Jtheta,Jz. The physical-to-q target factors are

    Pq=s^(2A)P, Jtheta_q=s^(2A-3/2)Jtheta,
    Jz_q=s^(2A-1)Jz.

The q-to-Q velocity factor is s^-A. Dividing the nonzero moment rows by
the SAME a=a(eta), with |a|>=a0>0, gives the exact current slow update

    u=a^-1[s^A P FP(R/sqrt(s))
                      +s^(A-1)Jz Fz(R/sqrt(s))],
    gamma=a^-1 s^(A-3/2)Jtheta Fg(R/sqrt(s)),
    psi=a^-1 s^(A-1)Jtheta Psi(R/sqrt(s)),
    beta=-epsilon partial_Z psi.                        (3)

In particular the two angular-response columns do not have a common
scale factor. The complete reserved swirl in this chart is
V=a s^(lambda-h)R^(-1-2lambda), not a R^(-1-2lambda).
Substitution and R=sqrt(s)x give EXACTLY

    int R^2 u=0, int R gamma=0,
    int (2V/R)u=-P, int R^2 V gamma=-Jtheta,
    -int RVu=-Jz.                                      (4)

The induced radial beta is retained. Equation (2) proves
D1 beta+epsilon partial_Z gamma=0. Since the response is independent of
the auxiliary variables, this ordinary compact primitive is the prescribed
T1 response, with A1=0 exactly for this slow update. This does not set A1
to zero for the fast map. The actual base radial velocity b is unchanged
and remains in the nonlinear remainder. No pressure tail has been removed.

## 3. Explicit all-fixed-order composition and reciprocal bounds

Use the inherited unweighted normalized slow norm |.|_m and the same
coordinate logarithm u_s=log(s/s0) at a real center. Its nonconstant
Taylor coefficients are bounded by the positive formal series

    L(z,t)=(1/3)[1/((1-32z)(1-64t))-1].

For any fixed real p, radial order r and smooth compact profile F define

    c=r-2p, D_(p,r,k)=||(x partial_x+c)^k partial_x^r F||_infinity,
    ell_(m,k)=sum_(i+j<=m)[z^i t^j] L(z,t)^k,
    G_(p,r,m)(F)=2^ceil(|c|)
               sum_(k=0)^m D_(p,r,k) ell_(m,k)/(2^k k!). (5)

At each center finite REAL Taylor expansion in u_s proves

    |partial_R^r[s^p F(R/sqrt(s))]|_m<=G_(p,r,m)(F).      (6)

Indeed the generator is x partial_x+r-2p and the center prefactor is
s0^(p-r/2)<=2^|c|. The profile need not be analytic; only the coordinate
has a complex extension. This extends the earlier density-weight formula
to the actual fractional response exponents in (3).

The reciprocal also needs no small high-derivative hypothesis. Let Aplus_m
bound the sum of positive-order normalized jets of a. At any center the
finite reciprocal series, with zero-constant increment, gives

    |a^-1|_m <= H_m=a0^-1 sum_(k=0)^m (Aplus_m/a0)^k.     (7)

For completeness Aplus_m has a finite fixed-profile bound. Put
D=1/2-h, eta=Z s^-D, a_k=sup_|eta|<=1 |a^(k)(eta)|/k!.
The coordinate power bound implies

    E_m=(2+1_(m>=1))(16/3) sum_(i+j<=m)32^i64^j,
    Aplus_m <= sum_(k=1)^m a_k E_m^k.                    (8)

Here |Z|<=2 on the real chart, and |Z|_m<=2+1_(m>=1).
Only positive-order eta jets enter Taylor composition; E_m is a safe
overestimate for their sum. At m=0 take Aplus_0=0. Finite Taylor's theorem
on the real coefficient range proves (8), including one-sided physical
endpoints, under the retained fixed smooth profile hypothesis. These are
finite constants uniform in the band, not uniform in order or a0->0.
Alternatively supplied positive polynomial S jet bounds give a positive
polynomial H_m(S) by (7); no hidden negative epsilon power is introduced.

## 4. Substitute explicit inverse norms into the coupled recurrence

Define Kp_m=G_(A,0,m)(FP), Kz_m=G_(A-1,0,m)(Fz),
Kg_m=G_(A-3/2,0,m)(Fg), Kpsi_m=G_(A-1,0,m)(Psi).
Leibniz and (3),(6),(7) give

    |u|_m <= H_m(Kp_m |P|_m+Kz_m |Jz|_m),
    |gamma|_m <= H_m Kg_m |Jtheta|_m,
    |beta|_m <= epsilon(m+1)H_(m+1)Kpsi_(m+1)
                                                |Jtheta|_(m+1). (9)

Mixed radial derivative bounds use (5) at the requested r with the same
one extra slow jet for beta. On the actual moving hull the earlier
reserved-support comparison adds its fixed exp(E_geometry) factor to
weighted norms. The union hull is used only for integration.

For incoming scalars bounded by epsilon^C times Pi_P,m,Pi_z,m,Pi_theta,m,
the parent's concrete polynomials are now

    U_m=H_m(Kp_m Pi_P,m+Kz_m Pi_z,m),
    Z_m=H_m Kg_m Pi_theta,m,
    X_m=(m+1)H_(m+1)Kpsi_(m+1)Pi_theta,m+1.              (10)

Substitute these into its Theta_P,m and Theta_z,m at field order m+2.
This closes the SLOW-MAP part of that conditional quantitative recurrence
without an abstract L_u,L_g,L_B assumption. Its output still has order
C+1, and the whole pressure moment epsilon partial_Z(Jz+c_rho P)
still uses incoming angular scalar jets through m+4. Both full pressure
memories and all six current covariances remain in the scalar inputs;
(10) does not replace those inputs by projected targets.

No actual first-band constants or canonical seed are evaluated. The
retained background/core floor, actual sourced/fast input norms, full
pressure/cutoff/source support bounds and the pointwise residual remain
to be discharged. There is no infinite-cycle, unconditional pulse or
blow-up conclusion. The next obligation is to join the actual canonical
seed and sourced/fast maps to these explicit scalar and slow-map budgets.

## 5. Evidence

Exact rational matrix solves check all five signed rows and the distinct
scale exponents. Separate compact piecewise-polynomial bumps realize the
response profiles at a manufactured lambda=1/2, not the selected physical
tilt. Their ordinary compact potential detects the filled support gap.
Independent bivariate generalized powers check (5), including its shifted
Euler generator. Reciprocal jets test (7) even when the positive-order
jet norm is large. Positive polynomial substitution checks (10) against
the prior coupled budget. Finite fixtures are diagnostics, not proofs of
the continuum estimates or evaluations of actual canonical fields.

```sh
python3 research/navier_stokes_cascade/results/core_slow_inverse_budget/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_slow_inverse_budget/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_slow_inverse_budget/evidence.json \
  --output /tmp/core_slow_inverse_budget_fresh.json
```
