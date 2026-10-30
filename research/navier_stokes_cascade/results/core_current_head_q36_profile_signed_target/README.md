# Q36: one fixed profile through Q35 and the original signed target

Parent: `e4aceb61c9c54d2f35174a402edf1a6e3ba9380a` (Q35).
This checkpoint binds Q4's chosen rho-squared projectors and the historical
logarithmic/two-edge estimates to the corrected Q35 initial prefix. It
removes independent density, moving-moment and projector norm callbacks.
The underlying formulas are inherited analytical results, not new
cancellations. The new executable composition does not require a signed
covariance object to define a pressure profile that is needed BEFORE that
signed target exists. It recomputes the chosen target from the current
prefix. Historical frozen evidence is unchanged.

## 1. Fixed profile before the signed target

Retain Q4's one fixed smooth profile rhohat compactly supported in (l,u),
0<l<u, with signed integral one. Set H_j=sup|rhohat^(j)|/j!. These are finite
for each requested order; H_j are inputs for this SAME profile. Neither
positivity of rhohat nor uniformity over different profiles is assumed.
The existing scale satisfies 1/4<=s<=4, with no R dependence, and

    rho=s^(-1/2) rhohat(R/sqrt(s)),
    c_rho=(s/2) integral x^2 rhohat(x) dx,
    sigma_e=R^(-e) sqrt(s) rho^2/A2, A2=integral rhohat^2, e=1,2.   (1)

Cauchy--Schwarz gives A2^(-1)<=u-l. Changing variables R=sqrt(s)x
proves integral R^e sigma_e dR=1 exactly. The actual support is
sqrt(s) supp(rhohat), strictly within its moving enclosing interval;
[l/2,2u] is only the common integration hull. The moving support, not a
rectangular product hull, supplies the existing reserved weight floor.
No identity with an unspecified pre-Q4 historical sigma is asserted.

Here is an explicit normalized through-order majorant. Let L_j be the
existing exact-degree log(s/s0) bound, and define

    d_0=2,
    d_n=(1/(2n)) sum_{j=1}^n j L_j d_{n-j},
    a_j=Rmax*d_j+d_{j-1}  (j>=1), a_0=0,
    P_{n,k}=[t^n](sum_{j>=1} a_j t^j)^k,
    B_n=sum_{i+j<=n} sum_{k<=j} d_i H_k P_{j,k}.                 (2)

The finite exponential chain rule bounds both sqrt(s) and its inverse by
d_j. Multiplication by R gives a_j for the increment of R/sqrt(s).
Finite Taylor composition and multiplication prove |rho|_through n<=B_n
in R,Z,T. In particular this norm is a sum of normalized derivatives, not
Q33's raw maximum. Products of through-order norms are submultiplicative.
Consequently

    Sigma_e,n=(u-l) [sum_{j=0}^n binom(e+j-1,j) Rmin^(-e-j)]
                     [sum_{j=0}^n d_j] B_n^2                 (3)

bounds sigma_e in that same norm. No current covariance or target is
needed for (2)--(3). The helper's finite ledger checks all requested orders.

For the pressure moment put cstar=(u-l)u^2 H_0/2. The inherited scale jets
then give

    C_rho,n=cstar [4+sum_{j=1}^n sum_{i=0}^j 32^i 64^(j-i)].    (4)

This is a bound for the WHOLE moving c_rho, including derivatives. The
implicit equation gives s_Z=2Z s^(2h)/(1-2h Z^2 s^(2h-1)), generally
nonzero. The axial term remains epsilon partial_Z(Jz+c_rho P), with both
c_rho,Z P and c_rho P_Z. It is not replaced by Jz alone or by a constant
second moment. Q35's full pressure inputs and tails remain required.

## 2. Same profile in every initial and projector slot

The factory accepts Q35's remaining fixed data with the rho key removed,
and supplies rho from (2) times exp(E), the reserved support weight factor.
It supplies c_rho from (4). It rejects an independent c_rho callback or an
additional rho key. All other original background, pressure-map, compact,
slow-map, support and canonical-root inputs remain conditional and must
refer to this exact profile/choice. This interface does NOT prove those
operators were constructed from it. The slow potential hull and profile
support hull are separately supplied; neither is silently enlarged.

For maximum target order m, supply profile jets through m+13, the same
safe ceiling as Q35's fixed-data/covariance requests, and retain

    p=12*(m+13)+18.                                             (5)

This reserve covers intermediate fast and fresh scalar calculations, not
just the projector's visible m jets. Root/canonical derivative expansion
has its own retained higher-order ledger; (5) is not its primitive ceiling.

The existing three-variable logarithmic pullback PB_e,n acts on the SAME
R hull and physical weight. With (3), its weighted projector coefficient is

    Bump_e,n=exp(E) PB_e,n Sigma_e,n,
    H_e,n=PB_e,n Row_e,n+Bump_e,n Moment_e,n.                     (6)

The support floor is used once in Bump; (3) is unweighted. The separate
floor used for rho belongs to the density's own slot. Both use the same
profile support. Moment_1 includes (4) and the complete pressure defect;
Moment_2 is the angular identity. Q35 returns (6) with epsilon order
alpha=2-3/100000, weight zeta and pole 2p. Exact zero moment follows from
(1) and the original row moment identities, not from these upper bounds.

## 3. Two edges and the unchanged signed root

The new source adapter transfers exact rational exponential coefficients
between the historical helper rings, including every exponent; it does
not evaluate or round giant exponentials. It passes Q35's actual norm
recipe outputs to the existing SignedTarget. No manufactured covariance
subclass, independent current row or independently assumed signed norm is
used in the factory.

The provider recovers the same shell ratio Xb/Xa, integrates forward at
the left edge and backward at the right using zero total moment, and
retains both terms in the full weight. Its native U-chart target estimate is

    |J_n xi_e| <= V_e,n delta^(-2p-3n),                          (7)

where V_e,n is the explicit finite convolution of (6), the two-edge
primitive, R^(-e), reciprocal weight and nonlinear native coordinate jets.
All corresponding weight and geometry hypotheses remain those of the
cited providers; no shell, width or accepted band is changed. The target
route requests source jets only through n. Complete signed curl requests
n+1 and is refused if that exceeds m.

The signed amplitude recipe uses the SAME original positive primary root,
inverse, reciprocal and complete curl data. It does not solve a new
positive root at the signed target. Its principal epsilon order is
3/2-3/100000 and curl-remainder order is 2-4/100000. As before, these
conditional bounds do not cancel the signed square or later covariance.
The historical target ledger conservatively continues to list actual
projector/root premises. Here its projector norm is supplied by (3);
actual fixed profile mass/support/jets and all original field premises
remain required. No historical ledger or frozen assertion is overwritten.

## 4. Validation and remaining obligation

Controls check exact weighted normalization and the moving second moment
for a finite-regularity polynomial fixture, independently differentiate its
R,T density and projector jets, retain nonzero higher scale jets, and
check exact exponential-coefficient transfer. The end-to-end norm fixture
uses Q35's real recipe chain (with manufactured lower constants), verifies
the profile slot and distinct moments, the original root identity, and the
complete curl's next target jet. It is NOT the selected continuum profile
or an evaluation of the actual giant field. Historical mixed-coordinate,
reciprocal-weight and two-edge tests are also rerun.

Q35 and its frozen closure are replayed. On 30 October 2026 parent focused
CI 41118800252 passed; broad CI 41118800258 was still running at inspection.
Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.

Next compose this freshly chosen signed target with all six complete-wave
covariance entries and fresh fast/slow means, retaining old means and
both pressure archives, then bind the first full unaveraged source. The
remaining precise promotion requirement is an actual same-field weighted
source estimate, uniformly over the certified physical bands, with the
recorded derivative and pole ledger. These positive norm recipes alone do
not verify its fixed-map/root premises. Uniform actual pulses, complete
viscosity/cutoff/interaction residual improvement and blow-up remain open.
