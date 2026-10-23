# Derive the current scalar gain from the coarse actual mean bounds

Parent: `a68fe336add908e03b0790688235c48bd6469a05`.
The parent's sufficient gain route left both post-wave fast rows at
D=3/2-2kappa and fresh postfast scalars at C=3/2-kappa unproved.
Here the **scalar part follows from the existing coarse mean constructor**:
no new incoming scalar estimate, sharper initial Haar theorem, or fast-row
D estimate is needed. The fast-row part remains open. The complete
nonzero source therefore still has only the parent's B=1-kappa guarantee.

The exact moment subtraction identities are historical results of
`core_dynamic_mean_recheck` and `core_coupled_defect_budget`. This checkpoint
applies their complete three-scalar version at the SUBUNIT order B<1 and
composes it with the actual first-prefix and generated-wave bound recipes.
It neither reproves those identities from numerical samples nor advertises
a new general mean-cycle theorem. Historical frozen evidence is unchanged.

## 1. Retained premises and the precise conditional statement

Use the SAME fields, common chart, finite requested jet order, actual A1
fast velocities, fixed five-map and unmodified supports as the parent.
The retained coarse mean theorem gives tangential means and increments
at epsilon^B and radial ones at epsilon^(B+1), with B=1-kappa,
0<kappa<1/4, epsilon<=1. The base b has an extra epsilon; V and G are
bounded. The generated complete sourced velocity has order B, its full
covariance CHANGE has order C=B+1/2, and its pressure is retained.

The actual five-map's linear identities use the original reserved support
contract: **G=0 on the ENTIRE slow potential hull and V is the same
power-law background there**, with the actual fixed inverse floor.
This is already a premise in `compact_mean_foundation` section 2 and
`complex_mean_relaxation` section 2; it is not newly certified here.
The hull includes the gaps between the axial bumps. Smoothness or support
of the bump values alone does not prove this property for the summed base.
Original physical input, frame/root identification and support enclosures
remain conditional exactly as in the parent.

Relative to those retained premises, all three ACTUAL postfast scalars
P,Jtheta,Jz have order C, with explicit positive polynomial majorants
derived below, and the next actual slow increment has tangential order C
and radial order C+1. This is not a claim that the pointwise mean pressure
has order C, or that the full tangential mean rows have order D.

## 2. Reset the old scalars using the last prescribed slow map

Let (beta,v,gamma) be the actual postfast old mean just BEFORE the last
slow solve in `means.before`, and (X,u,z) that solve's slow increment.
Its wave is held fixed. Put bars for auxiliary Haar averages and define

    S=2(b+bar(beta))X+X^2,
    H=(b+bar(beta))z+(G+bar(gamma))X+Xz,
    U=2bar(v)u+u^2, K=bar(gamma)u+bar(v)z+uz,
    Zf=2bar(gamma)z+z^2,
    tau=epsilon partial_T+epsilon^3 partial_Z^2.

Signs in tau are immaterial for norms. Exact current subtraction, followed
by the five-map's cancellation of the three CURRENT linear scalar rows,
leaves

    Pnew=tau int X + int (U-S)/R - epsilon partial_Z int H,
    Qrem=tau int R^2 X + int R(S+U)
                                  -epsilon partial_Z int R^2 H,
    Jtheta,new=int R^2 K,
    Jz,new=int R Zf - Qrem/2.                           (1)

The old incoming scalars have disappeared by their signed linear solve,
not by replacing the old mean state with zero. Radial viscosity has zero
plain and R^2 moments by compact integration by parts; axial viscosity
is the epsilon^3 term in tau and remains. The exact angular-zero pressure
is still T0(gr-rho P), and the original full background Pi1/Pi2 is retained.
Both old/new mean correlations and the nonlinear new squares occur in (1).

At tangential order B and radial B+1, the potentially weakest terms U,K,Zf
have order 2B. Time, base-radial and G*X flux terms have order B+2; the
other radial terms have order 2B+2. Thus **0<B<=2** implies order 2B
for every scalar in (1). This uses only the parent's coarse old fields
and coarse actual slow increment, not their outgoing scalar bound.

For through-n factorial slow norms, use fields through n+2. Let M bound
all radial moment operators with weights -1,0,1,2 and their chart jets;
the retained Kmom recipes provide such an enlarged bound. Replace the
six fields by their nonnegative norm polynomials x0,u0,z0,x,u,z, form
the positive versions S,H,U,K,Zf above, and put

    t_n=(n+1)+(n+1)(n+2),
    A_n=M[t_n x+S+U+(n+1)H].

Then `slow_reset` returns (A_n, M K, M Zf+A_n/2) as sufficient P,Jtheta,Jz
polynomials at order 2B. M separately dominates each radial weight, so
using the same A_n for the plain and weighted remainder only enlarges it.
G*X is conservatively retained even though the reserved-hull premise
annuls it on the actual slow support.

**Subunit-order correction.** It would be wrong to reuse the historical
postslow C+1 conclusion with C=B<1. Its square has order 2B<B+1.
For a compact swirl with int R^2 u=0 and zero old means, int u^2/R>0;
its remainder is exactly of order 2B. The ratio to epsilon^(B+1) is
epsilon^(B-1), unbounded as epsilon tends to zero. We use 2B throughout.
No historical record asserting the C>=1 theorem is changed or refuted.

## 3. Generated wave and fresh actual fast correlations

Starting from that old POSTSLOW state, the sourced wave changes complete
covariance by bar(w tensor q+q tensor w+q tensor q). Its order is C;
all six entries and both quadratures remain. The full pressure-moment
formula already established by `core_sourced_covariance_moment_join`
gives a three-scalar change of order C. In particular Dzz enters Jz,
Dtheta,z enters Jtheta, and the axial derivative of Drz enters P and Jz.
The wave support is the whole actual wave hull, not the reserved slow hull.

Next let (Xf,uf,zf) be the actual fresh fast increment, INCLUDING gd-A1 gd.
Every component is exactly zero-Haar. Base-linear, time and viscous changes
therefore vanish only after Haar averaging. With old full means, define

    Sf=bar(2beta Xf+Xf^2),
    Hf=bar(beta zf+gamma Xf+Xf zf),
    Uf=bar(2v uf+uf^2), Zf=bar(2gamma zf+zf^2),
    Kf=bar(gamma uf+v zf+uf zf).

The scalar changes are (1) without tau, using these fast products. In
particular Jz uses -Qrem/2, not just int R Zf. Each has order 2B using
the coarse fast bound ALREADY proved by the parent. No fast-row gain or
vanishing old Haar means is used. `fast_change` keeps separate component
norms so that axial and angular correlations cannot be silently equated.

Since

    2B-C=1/2-kappa>0,

the old reset at 2B, complete wave change at C, and fast change at 2B
sum to an ACTUAL postfast scalar bound at C on the SAME accepted band.
`ScalarGainBudget.postfast` computes exactly this sum. Its constructor
takes only the existing `FreshMeanBudget`; it accepts no new scalar
callback. `slow_increment` applies the SAME slow profiles to these fresh
derived scalars, including the extra Jtheta derivative for radial velocity.

## 4. Pressure memory, finite poles and scope of the join

For c_rho=int R^2 rho/2 the full axial moment is Jz+c_rho P. Its derivative
is the WHOLE product derivative, including (partial_Z c_rho)P. Calling
`whole_axial_moment` with a through-order c_rho norm keeps all these jets;
one extra derivative supplies the residual moment at C+1. This is an
integrated pressure identity; pointwise pressure still has its old B bound.
Actual A0 pressure errors and exact additive temporal tails stay in their
original ledger and are never supplied to an inverse.

Use the parent's fixed x=S+K, higher-pole recipes and unchanged physical
width C^-20. Old means have at most 2p and fresh fast means at most 6p;
put all field norms at 6p before products at 12p. The existing moment
recipe covers 12p. Lower powers are embedded with the explicit harmless
constant E=max(1,L)^ceil(12p), also used when applying that moment bound
to a lower-pole expression. Indeed delta^q<=max(1,L)^(q-r)delta^r for
0<=r<=q<=12p on delta<=L. This avoids an unjustified assumption that
moment constants are monotone in the pole on a long shell. Scalar outputs
are integrated; no new physical support or frequency scale is introduced.
All polynomials use finite original jets; fixtures report their actual
highest dependency, rather than certifying original field enclosures.

The fresh slow mean/current-wave feedback now has order C>D. The fresh
FAST feedback still has only B<D. Thus the remaining sufficient gain
obligation is the actual post-wave zero-Haar tangential row estimate at D,
or the parent's exact combined feedback inequality (8). This checkpoint
does not change its source/inverse interface to advertise D. Uniform
actual inputs/supports, all pointwise mean equations, viscosity/cutoffs,
complete residual, infinite construction and blow-up are still open.

If G is nonzero on any part of the slow potential support, the block solve
leaves linear leaks int R^2 G u and 2 int R G z in the angular/axial
scalars. They can have order B. Our compact exact control makes both
nonzero even with the two zero velocity moments. Therefore this scalar
gain must NOT be presented as an independent validation of that original
support contract, or extended to an arbitrary background with the same
coarse size bounds.

## 5. Obligation ledger and reproduction

Removed relative to the retained actual-map/support premises: the extra
fresh postfast scalar-C hypothesis in the parent's sufficient route,
and its fresh slow-C size consequence. Added analytical assumptions: none.
No fast-row estimate, actual support enclosure, primary-frame identity,
new accepted-band threshold or improved full residual is marked proved.

Run `python test_audit.py -v`, then
`python audit.py --verify-record evidence.json --output /fresh/path.json`
from this folder. Exact controls compare the reduced polynomial with the
independent unreduced dynamic budget at rational subunit powers, test
nonzero axial/support leaks and whole pressure products, and compose the
same generated inverse and original fast/slow recipes. They are finite
algebra/provenance checks. The continuum estimates are the derivation
above, conditional on its explicitly retained analytical premises.
