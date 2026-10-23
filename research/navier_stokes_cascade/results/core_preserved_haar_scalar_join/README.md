# Preserved Haar means and the three complete current scalar remainders

Parent: `7e3acac4a78c3064bf2227cc35f9268892fae423`.
The latest prefix has regular nonzero-angular order D=3/2-2kappa,
complete regular mean-row order C=3/2-kappa, and a conservative postslow
scalar bound B+C, B=1-kappa. This checkpoint refines the THREE scalar
remainders to C+1, conditional on the original weighted canonical input
recipes. It supplies a finite positive-polynomial route for the old Haar
means needed by that refinement. No incoming fine scalar/mean callback
is substituted for that route.

The qualitative preservation of Haar orders 1 and 2 and the angular-only
cycle estimate were already known (`core_scalar_cycle_budget`). The
complete slow moment identities were already known
(`core_coupled_defect_budget`). The new result is their quantitative join
to the recent first-signed/sourced prefix, all three scalar channels,
the moving axial pressure product and current radial compensator. The
older B+C estimate remains valid; it was conservative, not an error.

The construction still depends on original weighted seed, complete cut
signed coefficient, chart-transfer, reserved-hull and operator estimates.
This checkpoint explicitly lists their missing actual inequalities. Its
diagnostics evaluate manufactured lower norm data, not those actual
fields. Historical frozen evidence and the actual roots, pulses, width
C_setup^-20, carriers, cutoffs, accepted band and pressure are unchanged.

## 1. Average only the factor that is multiplied by a slow increment

Write bar for normalized auxiliary Haar average. An actual slow
increment (X,u,z) is auxiliary-independent, so

    bar(beta X)=bar(beta) X,
    bar(v u)=bar(v) u, bar(gamma z)=bar(gamma) z.       (1)

This does NOT replace bar(v gamma) by bar(v)bar(gamma) in an old scalar.
For an actual fast increment, including its compact A1 velocity, each
component has zero Haar, but bar(v u_fast) generally does not vanish.
Its complete correlations keep the full old mean estimates at B. No
improved fast size, full pressure-field rate or nonzero-angular source
order follows from (1).

For slow changes the retained fixed base has b=epsilon*b0 and G=0 on
the WHOLE reserved potential hull. V is not constant or zero. After
subtracting the exact current linear five-map rows put

    S=2(b+bar(beta))X+X^2,
    H=(b+bar(beta))z+bar(gamma)X+Xz,
    A=2bar(v)u+u^2,
    K=bar(gamma)u+bar(v)z+uz,
    F=2bar(gamma)z+z^2.

For Lr=partial_R^2+R^-1 partial_R-R^-2 the radial source remainder is

    Rg=epsilon partial_T X + epsilon^3 partial_Z^2 X
          -D1 S-epsilon partial_Z H+A/R+epsilon Lr X. (2)

The time sign in a different t* convention changes the corresponding
signed time term, never its norm bound. Every mean is still the same
actual accumulated mean and W is the same complete six-entry covariance
during this slow operation. Incoming scalar cancellation uses FRESH
postfast P,Jtheta,Jz; stale rows leave their exact uncancelled differences.

Compact whole-line extension gives int Lr X=int R^2 Lr X=0. It does
not delete pointwise radial viscosity; axial viscosity remains. Thus,
with tau X=epsilon partial_T X+epsilon^3 partial_Z^2 X,

    P+=tau int X+int (A-S)/R-epsilon partial_Z int H,
    Q+=tau int R^2 X+int R(S+A)-epsilon partial_Z int R^2 H,
    Jtheta+=int R^2 K,
    Jz+=int R F-Q+/2.                                (3)

These are the postslow remainders AFTER the exact linear solve, not just
uncorrected scalar differences. They retain the axial velocity covariance
and the weighted pressure moment. With the SAME moving unit-mass rho,
c_rho=(int R^2 rho)/2 and p+=primitive(Rg-rho P+),

    epsilon partial_Z int R(F+p+)
       =epsilon partial_Z(Jz+ + c_rho P+).             (4)

In particular (partial_Z c_rho)P+ is retained. The total current pressure
contains previous wave/mean pressure and has not been assigned the rate
of this remainder. Pi1/Pi2 and all compact pressure memory stay in the
existing mean equations. The additive A0/A1 residual tails remain in the
parent's separate ledger; actual A1 velocities stay in regular products.

## 2. A complete positive-polynomial C+1 budget

At n use through-(n+2) normalized coefficient bounds, on the SAME
weighted chart and complete reserved hull:

    bar(beta) : epsilon^2 ox, bar(v),bar(gamma) : epsilon(ou,oz),
    X : epsilon^(C+1) x, u,z : epsilon^C(u0,z0), C>1.

Let M dominate the required plain/weighted radial moments with the
original edge absorption and moving chart derivatives. Put

    s=2(b0+ox)x+x^2, h=(b0+ox)z0+oz*x+x*z0,
    a=2ou*u0+u0^2, k=oz*u0+ou*z0+u0*z0,
    f=2oz*z0+z0^2, t_n=(n+1)+(n+1)(n+2),
    ThetaP=ThetaQ=M[t_n*x+s+a+(n+1)h],
    ThetaTheta=M*k, ThetaZ=M*f+ThetaQ/2.               (5)

Positive products and epsilon<=1 give (3) bounded by epsilon^(C+1)
times (5). Limiting terms are old Haar tangential times new tangential,
of order C+1; the new squares have order 2C>C+1. Radial time/base/axial
terms have order at least C+2. `haar_slow_reset` implements exactly this
envelope by the inherited complete scalar formula. An independent raw
physical calculation keeps all epsilon factors and checks domination.

Consequently the current C=3/2-kappa prefix has

    (P+,Jtheta+,Jz+) : S_(5/2-kappa),
    rho P+ : M_(5/2-kappa),
    epsilon partial_Z(Jz+ + c_rho P+) : S_(7/2-kappa).  (6)

The last coefficient is
(n+1)[ThetaZ_(n+1)+C_rho,n+1 ThetaP_(n+1)], where C_rho,n+1 bounds
ALL jets of the moving multiplier. Compared with the former B+C=5/2-2kappa
bound this gains kappa, not a new full iteration. Full regular order is
still D, and the parent's negative current physical exponent remains.

## 3. Derive old Haar bounds from the same finite prefix

The original initial fast step has zero Haar in all components, including
A1. Hence its initial Haar field is exactly the original initial slow
increment. Use `SeedBudget.fields(n)` for its tangent order 1/radial 2
and `SeedBudget.output(n)` for its complete initial S2 defects. A SAME
chart/hull transfer T_n must include all normalization, chart and weight
derivatives; an unweighted bound alone is not silently weighted.

Let alpha=2-3kappa. Bound the original complete CUT first real signed
velocity at alpha-1/2 by V_n, obtained by the original target/root/curl
recipe, with cutoff curl terms. Use a SEPARATE original complete cut
primary coefficient bound U_n at order 1/2. A bound for the updated
total wave alone cannot bound that original summand: cancellation can
make the total smaller. Both U_n and V_n must come from their original
lower coefficient recipes. All six covariance changes are bounded by

    Delta_n=2 U_n V_n+V_n^2, at order alpha.            (7)

This conservative l1 estimate retains every mixed term and signed square;
it does not apply the real-primary quadrature shortcut to later complex
waves. Feed (7) and its next jet into `wave_change` for complete P,Jtheta,Jz.
The first actual fast change is `fast_change` on the FULL initial means
and `before.fast`, with order 2B. Therefore a derived postfast envelope is

    A_n=T_n SeedBudget.output(n)+wave_change_n
                                    +first_fast_change_n,
    first_postfast order=min(2,alpha,2B)=alpha.         (8)

The SAME fixed slow response applied to A_n,A_(n+1) gives a finite
`first_slow(n)` bound at alpha/alpha+1. It uses current postfast data,
including pressure and axial velocity terms, not a relabelled coarse
slow norm. Since alpha>1 for 0<kappa<1/4, its Haar field can be embedded
in tangential order 1 and radial order 2. Subsequent first sourcing changes no mean velocity, and
the latest actual fast step has zero Haar. Just before the latest slow
step the exact Haar majorant is the FINITE sum

    old_Haar_n=T_n SeedBudget.fields(n)+first_slow(n).  (9)

No infinite sum is used. `PreservedHaarScalarJoin` composes (5),(7)--(9)
with the existing scalar-derived `ScalarGainBudget.slow_increment` at C.
All lower-pole fields are conservatively embedded into the current 12p
moment domain by the existing L^ceil(12p) constant. Products, integration
and slow maps keep the same original weight/hull contracts; no arbitrary
unweighted inverse constant replaces them. For every fixed request choose
the original finite pole to cover all nested jets, then evaluate the
existing estimate shifts at that pole. The physical fields do not change.

The initial-primary scalar channel of (9) at output n reads Omega through
n+13 (n+14 for the whole axial residual moment). These counts concern
the already-curled aggregate primary coefficients in the initial recipe.
Other inputs, especially the latest sourced/fast recipes and their huge
compact-integration jets, retain their earlier larger requirements. This
is not a total original-amplitude derivative count or an evaluation of
those jets. The deferred A0 raw-source order 100005 and all nested A1
counts from the preceding checkpoints remain unchanged.

## 4. What the recipe still requires from actual fields

The generic `PrefixBudget.initial` stores only B/B+1 size bounds. It does
NOT certify the sharper Haar norms in (9). The precise missing inequalities
are the same-hull WEIGHTED initial fields/defects bounded by T_n times the
original seed expressions, and the separate original complete cut primary
and original-root signed bounds U_n,V_n in (7), for every requested jet
on the whole accepted domain.
Their original weighted coefficient/support/response premises must be
validated; no test or newly supplied fine output callback establishes them.

Indeed an epsilon^B norm bound alone allows a slow Haar field of that
order, whose epsilon^-1 rescaling grows like epsilon^-kappa. At
epsilon_j=2^(-j/kappa) this factor is 2^j, so no fixed S-polynomial
coefficient can justify replacing that coarse bound by a fine order-one
Haar bound on an entire tail. This is an obstruction to that inference,
not an actual canonical counterexample. The preserved-sum route (7)--(9)
uses the extra structural identities and lower weighted estimates.

The executable certificate lists those inequalities and refuses to certify
actual input validity. This checkpoint consequently supplies a conditional
quantitative refinement, not independent actual-input validation or the
next complete signed/fast/slow cycle. Next compose the same original-root
GENERIC signed update with the current target/whole moments and prove its
complete covariance cancellation at 2-2kappa, then fresh mean rows and
compact tails with their original high-jet/support domains. Physical
absolute decay, a controlled infinite diagonal, force regularity and
velocity growth remain open obligations before any blow-up claim.
