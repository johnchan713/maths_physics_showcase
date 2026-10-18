# Canonical value domain: exact all-later-band certificate and first-band limit

Parent: `edf21536d1268a5d8061d663427128c361ec4970`.
Status: `conditional-dyadic-value-domain-certificate-not-first-band-verification`.
On 18 October 2026 the parent's hosted CI was running without a failed step.
Historical defining maps, sources, evidence and schedules are unchanged.

This checkpoint resolves the QUANTIFIER in the next gate: every finite set
of the reviewed actual value budgets admits one explicitly bounded later
half-line, common to all fixed derivative orders. It does NOT certify that
the already named first band is on that half-line. It gives the precise
remaining inequality for that claim, without changing qstar, the field,
pulse radius, width C_setup^-20, tilt, carriers, supports or either schedule.
Full pressure, including normal-source terms, axial means and Pi1/Pi2
memory, remains in the construction. No unconditional residual or blow-up
theorem follows. Exact finite controls below are not continuum proofs.

## 1. Actual inputs and two logically different gates

Use the SAME closed enlarged charts and frozen leading profile as
[leading_frame_foundation](../leading_frame_foundation/README.md), Sections
1--4. Its alpha, nu, f_star, g_star, lambda_star and c_star are analytic
positive floors, not the small rational inputs() example in frame_budget.py.
Multiply F,g,lambda floors by 4^(-1-h), and ceilings by 4^(1+h), on the
1/4<=q/Q<=4 charts. Carrier calibration gives, for chart lambda bounds,

    sqrt(lambda_minus/[4(1+u^2)^(3/2)]) <= Bs
    <= sqrt(lambda_plus/(1+u^2)^(3/2)).

A certificate call rejects missing primitive keys; it never fills them
from the parent's manufactured defaults.

Choose rational OUTWARD bounds for these finitely many positive values,
Rmin,Rmax,Fmax,gmin,gmax,cmin,cmax,lambda_min,lambda_max,bmin,bmax,u,
the actual background norms Db,Mz, fixed box bound Dbox and Haar bounds
0<Jmin<=Jmax. The fixed actual u may instead be enclosed monotonically
term by term; the implementation accepts rational u fixtures, not an
assertion that the actual u is rational. A fully numerical actual certificate
requires supplying and verifying those enclosures, including cone gamma.
All expressions below are real-variable inequalities before enclosure.
Finite C2 leading/C1 background-comparison bounds supply Db; the current
background jet bounds supply Mz. Their finiteness is relative to the same
reviewed core, compact modulation, exact repair and exterior defining maps.
No new independent spectral-floor assumption replaces that derivation.

From frame_column_foundation and canonical_primary_foundation, retain
CEa, CEr, Cd, mu>0, Ds>0, B0, D0a,D0r,D1 and normal threshold Ln.
The conservative mu=lambda_min/(1+3u/2) is valid in both frames. Let

    c=lambda_min*u^2/2 * [(1+9u^2/4)^(-2)+2(1+u^2)^(-2)],
    C=3lambda_max*u^2/2 * [(1+u^2/4)^(-1)+2(1+u^2)^(-1)],
    Dj=Cd+9CEj/8, m=min(cmin,u), j=a,r,
    j_j^- = Jmin exp(-2/3)(7/8)^2 exp(-2Dj-C),
    j_j^+ = Jmax (9/8)^2 exp(2Dj) sqrt(pi/c),
    Cdir_j=D0j+D1 (9/7)^2 exp(4Dj+C)/(2c).

The ordinary canonical CONSTRUCTION gate is

    Ls >= Ldirect := max(25,Ln,8CEr/mu,
                         [4 Cdir_r max(1/cmin,1/u)/gamma]^2).       (1)

Together with the stated target/edge hypotheses, (1) gives its OWN positive
roots yr=Hr^-1 T. It does not need closeness to historical roots. Historical
and canonical roots are not identified. The additional COMPARISON gate is

    Ls >= Lcompare := [12 CH/(j_a^- m)](2 kappa+1),
    CH=2Jmax B0^2 Ds exp(Dstar+Da_star) sqrt(pi/c),
    Dstar=Cd+max(CEa,CEr), Da_star=Cd+CEa,
    kappa=24 j_a^+/(gamma j_a^-).                                  (2)

It implies theta<=1/(2kappa+1) and the existing relative squared-root
error <=1/2. Failing to verify (2) is NOT proof that canonical roots fail
to exist. Neither is checking a small manufactured matrix a proof of (1).
For a common comparison domain we include both inherited and canonical
normal/Riccati/cone bounds and (2) in Lreq. These are sufficient, not
necessary or sharp, thresholds. Manuscript Lemma 7.1/Proposition 7.5 allow
one value-domain selection before derivatives; that is not permission to
retroactively assert that an unspecified existing band passed it.

## 2. Rational dyadic majorants: no giant exponential evaluation

Write b(x)=ceil(log2 x) for positive rational x, and ceil for rational
ceiling. Compute b exactly by integer bit lengths and a single comparison;
no float, underflow, transcendental rounding or huge 2^exponent is used.
The inequalities e^x<=2^(2x) for x>=0 and
sqrt(pi/c)<=G:=2(1+1/c) give

    Ja = b(Jmin^-1 (8/7)^2) + ceil(2(2/3+2Da+C)),
    Pa = b(Jmax (9/8)^2 G) + ceil(4Da),
    H  = b(2Jmax B0^2 Ds G) + ceil(2(Dstar+Da_star)),
    K  = b(24/gamma)+Pa+Ja,
    T  = b(12/m)+H+Ja+max(K+1,0)+1.                              (3)

Then 1/j_a^-<=2^Ja, j_a^+<=2^Pa, CH<=2^H, kappa<=2^K and
Lcompare<=2^T. The last addition bounds 2kappa+1 by
2^(max(K+1,0)+1), retaining the +1 rather than dropping it.
For both j=a,r put

    Vj=b(D0j+D1(9/7)^2/(2c))+ceil(8Dj+2C),
    Ej=2[b(4 max(1/cmin,1/u)/gamma)+Vj].                          (4)

Since 4Dj+C>=0, distributing its exponential over the positive D0j
term gives Cdir_j<=2^Vj. Thus the cone length is <=2^Ej.
Set

    M=max(b(25),b(Ln),b(8CEa/mu),b(8CEr/mu),Ea,Er,T).            (5)

Lreq<=2^M, and M depends on finitely many VALUE budgets, never a desired
jet order. The large constants are retained, not replaced by convenient
diagnostic values. Rational enclosures exist for each fixed finite real
budget, but this existence argument is not an implemented actual enclosure.

## 3. One explicit half-line, with the ORIGINAL physical scales

Q=2^-ell, S=ell^2, epsilon=Q^h, k=ceil(epsilon^-1/2), and

    Ls=2r0/ci >= 2r0 ell^2,  0<r0<=1/6057072.                   (6)

Ls is the pulse length, NOT log(1/Q). There is no radius-free replacement
of (6). Let ell_old>=16 be the already accepted band floor and define

    B=max(4, b(ell_old), b(512/h^2),
          ceil([M+b(1/(2r0))]/2)), E=2^B.                        (7)

For every integer ell>=E, (6) implies Ls>=2^M. Also k^-1<=epsilon^1/2,
so S^2(epsilon+epsilon^2+k^-1)<=3ell^4 exp(-h ell/4), using ln2>=1/2.
The eighth exponential-series term gives the ALL-LATER estimate

    3ell^4 exp(-h ell/4)
    <=3*8!*4^8/(h^8 ell^4)
    <=3*8!*4^8/512^4 < 1,   ell>=512/h^2.                       (8)

Thus (7) simultaneously ensures the whole-pulse rounding/axial condition,
both normal/Riccati/cone domains and relative-root comparison for ALL later
bands. It does not differentiate rounded carriers or select a new start
for each fixed jet. Every fixed jet uses these same positive denominators.
The source-support, weighted edge and complete covariance estimates still
need their own retained hypotheses; value acceptance alone is not them.

## 4. Unchanged domain versus a certified observation tail

The exact sufficient inequality for the old analytic floor is

    2r0 ell_old^2 >= Lreq,
    ell^4(2^(-h ell)+2^(-2h ell)+k^-1)<=1 for all ell>=ell_old.    (9)

For the particular previously fixed observation cylinder use its actual
minimum contributing band instead of ell_old. None of the frozen records
supplies an actual-enclosure evaluation of (9); they explicitly disclaim
it. Consequently this checkpoint does NOT promote the whole old domain or
the already localized force to verified acceptance. (9) is sufficient and
its missing evaluation is an obstruction to that certification, not a
necessary nonexistence inequality or a counterexample to the construction.

A rigorously certified tail can be DESCRIBED without mutating the old field:
intersect its domain with q<qtail=min(qstar/64,2^(-E-2)). Any intersecting
enlarged chart has Q<=4q<2^-E, hence ell>E, and 4Q<qstar/4. This restricts
observation only. If the canonical fields have only been defined under an
acceptance premise, it certifies where that same defining formula satisfies
the value gates; it does not assume a field exists outside its domain.
No earlier cutoff, localized force or schedule is silently replaced.
The former endpoint cylinder, with q_e=min(qstar/16,2^(-ell0-2)), has
Q<=3*2^-ell0/8 and hence ell>=ell0+2. Its acceptance follows, for example,
if ell0+2>=E, but that inequality also must not be presumed.

Next obligation: verify actual rational outward enclosures and the existing
minimum-band inequality (9), OR explicitly formulate a new theorem on the
certified observation tail and separately propagate its localization/force
choices without identifying them with the historical ones. Then discharge
the remaining actual core/background and finite inverse/edge/contraction
premises. No percentage, completion time or unconditional proof is inferred.

## 5. Meaningful checks

`acceptance.py` implements (3)--(8) with exact Fraction/integer arithmetic,
keeps direct and comparison gates separate, and records no giant band or
qtail number. Rational test fixtures verify binary-log rounding, majorant
composition, whole-pulse scaling, Haar-radius costs, the retained +1, and
the difference between direct root existence and relative root comparison.
An independent mpmath logarithmic check tests the sufficient majorants on
manufactured budgets only. Exact domain witnesses test chart margins and
show why finite old-band membership cannot be deduced from eventual
acceptance. They are parameter-space controls, NOT fluid counterexamples.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output fresh.json` here.
The output must be fresh. Historical source bindings are checked recursively.
