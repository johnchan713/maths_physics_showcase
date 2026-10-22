# The retained Borel prescription implies the infinite diagonal

Parent: `18faf0b9ee43efe277d0b57c3e662222d80c6bd7`.
This closes an independently assumed diagonal budget for the construction
DEFINED BY the existing real scale prescription. No numerical infinite
list is needed. It also removes identification with the diagnostic Borel
cutoff as a premise of the preceding tangential-floor recipe.

This is a closure audit of an existing definition, not a new Borel theorem
or a newly selected schedule. Its conclusion is relative to the retained
full coefficient regularity certificate. It does not certify a separately
fixed sequence whose adherence to that definition is unknown, actual
leading gaps/chart coverage, the original first band's normal floor, a
complete residual improvement, unconditional pulses, or blow-up.

## 1. Separate three statements previously conflated

1. The exact-real prescription defines an admissible infinite sequence and
   implies its own diagonal for every n. This is proved below.
2. A diagnostic program evaluates a finite list with manufactured norms.
   This is not a certified evaluation of actual Navier–Stokes scales.
3. A separately supplied earlier sequence satisfies a later, more
   conservative bound. This requires its own inequality or identification;
   it is not a consequence of statements 1 or 2.

The repeated phrase “actual infinite schedule unverified” was too broad
when used as an additional obstacle to statement 1. The frozen
`background_jet_foundation`, equation (11), already supplies its defining
choice rule, including the correct Euler norms. The latest normal-floor
recipe only needs the diagonal that follows from that rule, not an
evaluation of every scale. We preserve every historical record and retain
the narrower unresolved issue in statement 3.

## 2. Constants exist before selecting any Borel scale

Fix the SAME h in (0,1/100), leading profile, positive compact native norm
annulus, and any SAME fixed smooth cutoff chi with the retained plateaus.
Use the all-order coefficient construction in `background_jet_foundation`,
sections 2–4, relative to its leading regularity, analytic collar, reserved
patch and exterior certificate. It gives every globally finalized
coefficient at finite order, with the full pressure, lower products,
viscosity, primitive and five-moment repair. It is an uncut formal hierarchy;
its coefficient equations do not depend on the later Borel scales a_n.

Thus for each finite n,m the seminorm M_(n,m) of the COMPLETE tuple
(E_n,U_n,W_n^r,J_n^r) on that compact norm domain is finite. For this SAME
chi, its own padded Euler seminorm K_j is finite for every finite j: all
positive cutoff derivatives have compact transition support. Neither a
uniform bound in n nor analyticity of chi is needed. These seminorms can
be defined as exact suprema; numerical upper enclosures are not necessary
for an existence theorem. No sampled norm is substituted for a supremum.

At each n,m, the product rule gives the already retained constant

    C_(n,m)=max(1,2 K_(m+1)(1+2nh)^m M_(n,m)).                 (1)

It controls the full normalized increment through total order m by
C_(n,m) q^(2nh). The radial term is chi W_n^r-(E chi)J_n^r,
so K_(m+1), not K_m, is required. This is in Euler q derivatives at
fixed X,eta. If frozen derivatives were used to form a numerical bound,
the old Bell/Stirling conversion remains compulsory. The theorem uses
the correct seminorms directly; it does not identify those operators.

Full pressure/moment endpoints and reserved supports are not changed by
this summation argument. No pressure is reset on a local zero-velocity
plateau. A missing leading/coefficient regularity certificate would still
block finiteness of M and hence this theorem's application.

## 3. Prove the ORIGINAL rule at arbitrary n, then use induction

Write b_n=log a_n, ell=log 2, B=log(q0^-1), q0<=1, and g_n=nh>0.
Equation (11) in the frozen foundation, combined as instructed there with
all inherited finite physical/weighted requirements, is

    b_n >= 1+max(b_(n-1)+ell, B+ell,
                max_(0<=m<=n) (log C_(n,m)+n ell)/g_n,
                all inherited finite threshold rows).        (2)

The inherited rows are retained, not replaced by the diagonal rows.
For example the frozen `background_summation` log-power requirement
C(1+|log q|)^P q^gamma<=2^-n, C>=1, P>=0, gamma>0, has threshold

    H=2(log C+n ell)/gamma+P^2/gamma^2.                        (3)

There are finitely many rows at each stage by the retained construction.
All constants are finite before the current choice. Given finite b_(n-1),
their maximum is finite. Recursion on the natural numbers therefore
defines nonempty admissible choices at EVERY finite n. Every sequence
admitted by (2), not merely a specially selected example, satisfies

    b_n-b_(n-1)-ell >=1,
    b_n-B-ell >=1,
    g_n b_n-log C_(n,m)-n ell >=g_n>0, EVERY m<=n.              (4)

The final identity follows by multiplying the relevant max-row slack by
g_n. It holds at arbitrary n,m; finite checks are not used to infer its
universal quantifier. Consequently a_n>=2e a_(n-1), and on the ENTIRE
support 0<q<=a_n^-1,

    C_(n,m) q^(nh) <=2^-n exp(-nh) <=2^-n, EVERY m<=n.         (5)

Outside that support the increment and all its cutoff jets vanish.
This is exactly the infinite diagonal needed by the preceding checkpoint,
and actually contains a harmless strict margin. It costs no extra Q
power and does not require a sharper cutoff constant or a changed band.

The inherited row (3) is compatible with the same proof: for s=-log q,
log(1+s)<=sqrt(s) and Young's inequality give

    log C+P log(1+s)-gamma s
      <=log C+P^2/(2gamma)-gamma s/2.

If b_n>=1+H, every s>=b_n makes this <=-n ell-gamma/2.
Thus retaining such rows in the same maximum preserves their all-smaller-q
guarantees too. This does not establish all physical residual hypotheses:
it preserves requirements whose finite constants/gains were justified.

`proof.py` checks the polynomial identities behind (4) and this last
implication, with symbolic n,gamma,log C and log 2, after clearing positive
denominators. It also binds the historical helper's all-row enumeration,
fixed gain, boundary, separation and padding expression using its AST.
These are algebra/source checks. The helper uses mpmath, not interval
arithmetic; no universal numerical-rounding guarantee is claimed for it.
The continuum theorem concerns the exact REAL rule documented in (2).

## 4. Feed the closed diagonal into the existing floor recipe

At total Euler/X/eta orders m<=1, every n>=1 obeys (5). Hence on the same
q interval, t=q^h<=1,

    |full normalized correction tuple|_1
      <=sum_(n>=1) 2^-n t^n=t/(2-t)<=t.                       (6)

This uses the complete tuple and its radial curl commutator. The preceding
`core_tangential_normal_closure` can therefore use Jfull1=(CF1+CG1)
(Lead1+1) and its Ctan for ALL choices admitted by (2), without a separate
infinite-diagonal assumption. Any fixed smooth chi works with ITS OWN K_j
in (1)–(2). Equality to the diagnostic squared-flat cutoff, and the
numerical upper bound K2=4928 for that cutoff, are not needed for this
particular floor. They remain relevant if a different estimate explicitly
substitutes that numerical K2 into its constants.

What is still needed for the actual unchanged physical band is

    S=ell_band^2,
    S^2(epsilon+epsilon^2+k^-1)<=1,
    4 Ctan^2<=beta S^2, beta=lambda_min/[4(1+U^2)^2],           (7)

with the actual physical chart/representatives, leading first-jet bounds,
shear/growth gaps and existing scale. The original length compatibility
beta Lreq^2>=16 r0^2 Ctan^2 is a sufficient alternative. None is inferred
from existence of (2). Full radial phase, sourced inverse, complete
covariance/mean and physical residual obligations remain separate.

Removed independent premises for this rule-defined floor construction:
the infinite diagonal and identity with the diagnostic Borel chi.
Added analytical assumptions: none. Retained premises: the original
leading/coefficient certificate, SAME real rule (2), all original finite
physical constraints and actual chart/gap/selector bounds. No a_n, chi,
qstar, physical support, carrier, width C^-20 or historical evidence is
chosen again or relabeled. If an independently fixed sequence is meant
instead of the rule-defined construction, its exact remaining condition
is b_n satisfying the SAME rows (2), or the direct inequalities (5);
unknown adherence is not cured by rerunning the choice rule.

Next obligation: verify physical representative/neighborhood coverage and
substitute the actual leading first norms and gaps into (7). Do not ask
for an infinite numerical scale list as a separate premise of this rule.

## 5. Reproduction

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output fresh.json`. Universal polynomial certificates,
source mutations, constant/order conversion controls and supplied-scale
counterchecks are exact. Finite fixtures are not actual NS coefficients
or continuum proofs. Historical evidence remains frozen.
