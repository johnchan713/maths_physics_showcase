# Same-weight primary covariance and initial mean-row substitution

Parent: `b919fe2154652bfc4df6586ef22918a61c147c07`.
The parent CI run 40042800072 was still in progress without failed steps
on 24 October 2026; its target/cutoff step passed. This checkpoint joins
the explicit primary norm factory to the EXISTING weighted covariance and
initial mean-row formulas. Their conditional rate improvement is historical,
not a new iteration or a proof of the full physical residual. All historical
records, physical roots, supports, widths and accepted bands are unchanged.

## 1. Weight equality precedes the atlas estimate

For each active band Q_i, the actual physical evaluation gives

    R_i=r/sqrt(Q_i), s_i=q/Q_i, Z_i=z/Q_i^D,
    X_i=R_i^2/(2s_i)=r^2/(2q), eta_i=Z_i/s_i^D=z/q^D.

Thus x=log(X/Xa), yb=log(Xb/X), zeta and delta are EXACTLY the same
physical functions for every overlapping band. This is equality of the
actual pulled-back weight, not comparison of distinct flat exponentials.
It retains both t1^2 and 4. The full common-chart chain rule still acts on
the coefficients; derivatives of the weight have already been included
in the root bounds. No quotient is differentiated as if it were constant.

Use the retained matched COARSEST band, fixed common index along radial
integrals, S_i<=(25/16)S and common derivative multiplier 188^n. All native
root inputs are first consistently expressed in U=R^d_r by the preceding
factory. Closed enlarged support separation implies at most one wave label
at a point, including derivative supports. Normalized finite-cover Haar
averages have no deck-count multiplier. Without those original geometry
premises, sums/cross-label terms would have to be estimated separately.

## 2. Explicit weighted complete covariance coefficients

Let A_j and C_j be the SAME primary amplitude and normalized full curl
polynomials returned by that U-chart root recipe. They bound normalized
exact-degree jets with weights sqrt(zeta) P delta^(-6j-3) and
sqrt(zeta) P delta^(-6j-9), respectively. Define

    a_n=16*188^n*n! sum_(j=0)^n A_j((25/16)S),
    b_n=16*188^n*n! sum_(j=0)^n C_j((25/16)S),
    W_n=2^(n-1)(a_n^2+b_n^2),
    Dcurl_n=2^(n-1)b_n^2, p_n=12n+18.                 (1)

The constant16 is the existing atlas factor8 times a safe component-to-
Euclidean factor2. These are raw through-order bounds BEFORE edge
unweighting. No Eedge, inverse zeta, label count or deck count is inserted.

The exact angular tensor of each real complete curl is

    W_ab=(epsilon/2)[A_a A_b+k^-2 B_a B_b].             (2)

All six entries remain. Mixed sine/cosine terms vanish by the real angular
identity, not pointwise. Since k^-2<=epsilon and the full B has the retained
epsilon^-kappa cost, Leibniz and sum_i binom(n,i)=2^n give

    |W|_n <= epsilon*zeta*delta^-p_n W_n(S),
    |W-Wprincipal|_n
       <= epsilon^(2-2kappa)*zeta*delta^-p_n Dcurl_n(S). (3)

P^2<=1 is safely dropped; zeta is retained. The pole in (1) covers the
worst square of the curl jet. The exact canonical principal assembly uses
the original squared slow partitions once, both signs already inside H,
and the physical target scale; it has not been rebuilt. Hence only the two
matched Haar entries acquire the fine difference estimate after subtracting
Sigma0. The other four entries are still order epsilon. Later complex waves
have nonzero principal-sine/curl cross terms and cannot use this square-only
reduction. This bridge explicitly concerns the initial real primary.

## 3. Substitute into the complete initial mean rows

For output n, the existing initial-row recipe reads fields through N=n+3
and full pressure inputs through n+1. Let S_N be the retained weighted
higher-stress bound for Sigma-Sigma0 at epsilon order3 and pole <=p_N;
the original weighted background proof supplies that input with pole
3N+3<=p_N. Then its original covariance slots can be replaced by

    W=Tr_N W_N, D=Tr_N Dcurl_N+S_N.                   (4)

Tr_N is a REQUIRED same-weight U-to-R tensor-jet transfer, not a constant
silently set to one. The precise unresolved contract is

    |F|_(R,through N) <= Tr_N(S) zeta delta^-p_N K_N(S)
    whenever |F|_(U,through N)<=zeta delta^-p_N K_N(S), (5)

for the actual common hull and coordinate composition. It must cover all
higher U(R) chain terms, coefficient-normalization and hull conventions.
Weight equality alone does not prove (5). The existing forward R-to-U
recipe cannot be reused as its reverse without a new derivation. The
current implementation therefore requires an explicit `row_transfer` norm
input. The diagnostic Tr=1 below concerns identical manufactured charts
only and does not certify this actual transfer.

The higher-stress order3 embeds into 2-2kappa without an epsilon loss.
`initial_rows` substitutes (4) into the ORIGINAL `row_budget`, retaining
every other base, weighted fast/slow field, axial, multiplier and full
pressure input. It rejects missing inputs. Independent diagnostic polynomial
classes are converted by exact coefficient copies; this is not a change of
physical coordinates or fields. The original row algebra, weighted map
constants, pole enlargements for products and pressure, and its sharper
INITIAL Haar rate remain precisely their historical conditional statements.
The bridge does not infer a complete row from covariance alone.

In particular the full pressure is reconstructed from its radial source
and rho times the WHOLE moving pressure moment. Pi1/Pi2, the axial velocity,
A0/A1 products and reserved supports remain in those lower inputs. Deleting
the pressure moment changes the axial row in the exact control. The older
unweighted aggregate/seed bounds cannot replace the remaining weighted
field/map inputs; this checkpoint makes no such inference.

## 4. Evidence and the next exact obligation

The code reproduces the historical weighted covariance coefficients and
complete row substitution, checks physical weight invariance, all six exact
real angular entries, the later-complex cross obstruction, factorial/Leibniz
counts and normalized Haar measure. Tests also compose the explicit original
primary factory into this bridge. These are manufactured finite controls,
not evaluation of the actual giant field or a continuum proof. Five frozen
dependency records and all their declared source hashes are rebound.

The next exact requirement is to derive (5) from the original positive
radial hull and fixed exponent, then supply the SAME-HULL weighted fast/slow
fields and complete pressure/profile inputs to `initial_rows`, including
their fixed-map/pole contracts. Existing conditional initial Haar and later
source/mean closures should be composed using those contracts, not reproved
as new gains. Full background/phase and fixed-profile validity, original
value/cone gates, all higher pressure/moment identities, supports, viscosity,
cutoffs, interactions and the complete physical residual remain open in the
overall ledger. No unconditional pulse, full residual improvement, infinite
realization or blow-up is certified.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
