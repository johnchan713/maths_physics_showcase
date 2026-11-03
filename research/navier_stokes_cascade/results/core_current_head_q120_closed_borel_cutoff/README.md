# Q120: closed implemented Borel cutoff norms and retained-scale slack

Parent: `c5f1d8e54d187ee17ab850b693c760baf184d0fb` (Q119).
On 3 November 2026 the freshly fetched remote and local tracked tree agreed.
Q119 research CI 41601360864 and field CI 41601360763 were running; Q118
research CI 41594347916 and field CI 41594347819 were also running.
The current notes and frozen evidence were read before this continuation.
Historical records, the leading-core/first-collar scope of
core_growth_activation, and the width C_setup^-20 remain unchanged.

This checkpoint closes the factorial growth of the already implemented
Borel cutoff certificate and derives an explicit sufficient inequality for
using that bound at an EXISTING scale. It does not certify that the actual
infinite sequence satisfies this inequality. No scale or cutoff is chosen,
replaced or enlarged. This follows Q119's closed geometry factors rather
than repeating the original compact-modulation/exact-repair audit.

## 1. The fixed function and its separate roles

The function in core_borel_cutoff_certificate is

    g(x) = exp(-1/x^2) for x>0, and 0 for x<=0,
    sigma(x) = g(x)/(g(x)+g(1-x)),
    chi(s) = 1-sigma(2s-1).

It has chi=1 for s<=1/2 and chi=0 for s>=1. Write E=s partial_s and
K_m=max(1, max_{0<=j<=m} ||E^j chi||_infinity). Derivatives outside the
transition vanish. The old Gaussian polynomial proof already establishes
smooth flat extension; this checkpoint does not assert analyticity of chi.
It also already proves K_0<=1, K_1<=64, K_2<=4928. Those sharper low-order
certificates are retained in the provider; they are not new results.

The theorem is specific to this squared-flat BOREL cutoff. Compact-X and
variance cutoffs have different roles and are not assigned this bound.
The source identity and AST bindings remain in the inherited evidence.
An unspecified smooth cutoff in a historical analytical ansatz does not
become this implemented function merely because a norm provider accepts
an identity string. Actual-field identity remains an explicit premise.

## 2. Closed all-order bound

All bounds in this section are analytical. Integer checks below only test
their implementation and finite consequences.

For x>0, g^(j)(x)=P_j(1/x) exp(-1/x^2), with

    P_0=1,   P_{j+1}=-u^2 P'_j+2u^3 P_j.

The degree is at most 3j. If L_j is the coefficient l1 norm, then
L_{j+1}<=(3j+2)L_j, hence L_j<=3^j j!. The existing monomial bound gives
A_j=sum_d |[u^d]P_j| ceil(d/2)! as a raw derivative bound. Since
ceil(3j/2)<=2j and (2j)!<=4^j(j!)^2, for every integer j>=0,

    A_j <= 12^j (j!)^3.                              (1)

On [0,1], H=g(x)+g(1-x)>=exp(-4)>1/81. The old reciprocal recurrence is

    B_0=81,
    B_j=162 sum_{k=1}^j binom(j,k) A_k B_{j-k}.

Induction using (1) proves

    B_j <= 81 * 2048^j * (j!)^3.                     (2)

Indeed after dividing by the right side, the sum is at most
162 sum_{k>=1}(12/2048)^k=486/509<1. The omitted factorial ratio is exactly
1/binom(j,k)^2<=1, so no unproved submultiplicativity is used.
Leibniz for sigma=g/H similarly gives raw bounds S_j satisfying

    S_j <= 81 * 2048^j * (j!)^3 /(1-12/2048)
         < 82 * 2048^j * (j!)^3.                    (3)

The chain factor 2^j gives ||partial_s^j chi||<=82*4096^j*(j!)^3
for j>=1. Expand Euler derivatives with Stirling numbers of the second kind:

    E^m = sum_{j=1}^m S2(m,j) s^j partial_s^j.

For m>=1 and 1<=j<=m, ordered set partitions imply

    j! S2(m,j) <= m! binom(m-1,j-1).

For example, arrange a permutation and choose j-1 cuts; this surjects onto
ordered partitions, so the inequality has the stated direction. Since
j!<=m! and |s|<=1 on the transition, (3) gives

    ||E^m chi|| <= 82 (m!)^3
                   sum_{j=1}^m 4096^j binom(m-1,j-1)
                 = 82 (m!)^3 *4096*4097^(m-1)
                 <= 82 *4097^m*(m!)^3.

This upper bound increases with m and also dominates 1. Therefore

    K_m <= 82 *4097^m*(m!)^3,   m>=1.                (4)

This conservative Gevrey-3 estimate is not asserted sharp. The code uses
1,64,4928 through orders 0,1,2 and (4) thereafter. Moderate-order comparison
with the independent old polynomial/reciprocal recurrences checks (1)--(4).

## 3. Existing-scale obligation, without reselection

Retain the COMPLETE full-background tuple norm M_(n,m) from
core_background_finite_head_budget, including pressure, its moments and
the radial cutoff-curl terms. Its normalized increment constant is

    C_(n,m)(K)=max(1, A_(n,m) K_(m+1)),
    A_(n,m)=2(1+2nh)^m M_(n,m),   0<=m<=n.         (5)

The extra cutoff derivative m+1 cannot be dropped. For the supplied
existing scale a_n, let ell_n<=log2(a_n) be a certified lower bound. A
sufficient original diagonal row is

    nh ell_n >= n + ceil(log2 C_(n,m)).             (6)

This implies the conservative bound C_(n,m) a_n^(-nh)<=2^(-n),
and applies on the smaller-q support used by the original prescription.
Input provenance and the other separation/domain/selector conditions
remain necessary.

Suppose an old certified bound Kold was used in the same row and the new
bound is Knew. With rho=max(1,Knew/Kold), positivity in (5) proves
Cnew<=rho Cold. Thus the existing row transfers if its old binary slack

    nh ell_n - n - ceil(log2 Cold)
        >= ceil(log2(max(1,Knew/Kold))).             (7)

The implementation also checks (6) directly with Cnew. Neither test chooses
an a_n. Failure of these sufficient tests leaves that row unresolved; it
does not prove the physical row false. Passing manufactured rows does not
certify the actual sequence or even an actual finite prefix.

For high requested orders, no factorial has to be materialized. From
82<=2^7, 4097<=2^13 and m!<=m^m, (4) yields the integer binary budget

    b_m=7+13m+3m ceil(log2 m), m>=3,
    b_0=0, b_1=6, b_2=13.

If the ACTUAL full tuple obeys M_(n,m)<=2^B_(n,m), B_(n,m)>=0, then

    ceil(log2 C_(n,m))
      <= 1+b_(m+1)+m ceil(log2(1+2nh))+B_(n,m).     (8)

Consequently the concrete deferred obligation at every retained row is

    nh ell_n >= n+1+b_(m+1)
                    +m ceil(log2(1+2nh))+B_(n,m).  (9)

The large-order diagnostic evaluates (8) at n=m=100010 using integer bit
budgets only. Its profile budget and scale lower bound are manufactured;
it is not evidence for the actual requested native source norm. A missing
actual B_(n,m) or ell_n is not silently replaced by these fixtures.

## 4. Low-order claims and the current interface

At s=3/4 the implemented function has chi_s=-16 and Echi=E^2chi=-12.
Thus K_m>=12 for m>=1. A proposed K_2=3 is impossible for this function;
the old manufactured low-order interface fixtures cannot establish a
physical comparison/selector bound. This lower bound was already known.
Values between 12 and a certified upper bound are undecided by these
bounds, not automatically false. The claim-bracket diagnostic states this
explicitly, and never asserts actual-field identity.

Q119's ClosedBackground accepts the identified provider and retains its
same full tuple, finite head, diagonal tail and closed real slow-chart
factors. The integration test compares its m=0,1,2 results with the OLD
BackgroundBudget using the SAME certified callback, not a different old
fixture. It checks the m=4 finite head and the requested cutoff order 5.
This establishes interface compatibility, not that old physical Cphi,
normal-floor or selector numerical values survive an increased input.
A conservative replacement upper bound must satisfy (6) or (7) and all
other dependent gates before an earlier scale argument can use it.

## 5. Validation and remaining obligation

The audit recursively verifies Q119 and its immutable dependency closure,
checks all inherited source/evidence hashes and freezes the new proof and
code hashes. Tests check exact reciprocal constants, ordered partitions,
finite recurrence domination, low-order impossibility/undecided cases,
identity guards, retained-row transfer, deferred large-order budgets and
current-background integration. A subprocess independently reproduces the
frozen record; passing these checks is not a continuum proof.

Next: bound the ACTUAL full profile tuple and verify the supplied retained
scale intervals against (9) at every required order (or prove (6) with a
sharper bound). Bind the physical Borel cutoff identity explicitly. Then
carry those actual constants through Q118/Q119 into the high-pole native
source budget and the compact-modulation/exact-repair cone join, retaining
complete pressure memory. Global pressure/moment and reserved-mean support
identities, uniform actual-pulse estimates over every physical scale,
complete-wave covariance means including axial velocity and pressure,
and viscosity/cutoff/all-interaction residual estimates remain open.
No unconditional actual-pulse bound, residual improvement or blow-up
conclusion follows from this checkpoint.
