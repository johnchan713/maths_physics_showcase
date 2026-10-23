# Canonical Haar-column jets and the joint column/wave root input

Parent: `f749704c9ee88e0758640e189a2c024dbe1f579b`.
Status: `conditional-canonical-Haar-column-jet-closure`.
On 23 October 2026 the parent's CI run 40024231582 was still in progress,
without failed steps; its canonical entrance/mixed-jet step had passed.
All historical evidence, actual roots, widths, bands and profiles remain
unchanged. This closes the next quantitative INPUT identified by the parent.

## 1. Which covariance is being differentiated

Use the ORIGINAL real canonical homogeneous pulse t_sigma, its reference
entrance Br(0)eplus and all full moving-normal pressure terms. Write
x_sigma=(t_sigma)_r. The two columns defining its original positive roots
are the angular and auxiliary Haar averages of the PRINCIPAL wave flux:

    H_sigma=(J/2)(integral chi^2 dxi) ci
             integral_0^L psi(v/L)^2 x_sigma t_sigma,tan dv,
    ci=2r0/L, s0=r0^2/sqrt(L), Hbar=H/s0.                (1)

J=|det(vr,vt)|=1+(sqrt(2)-1)^2<2. The angular factor is EXACTLY 1/2:
the original kp is a nonzero integer, so the angular average of the shifted
cosine square is 1/2 before any slow derivative is taken. There is no large
phase derivative left over from that exact identity. The original fixed
transverse cutoff has integral chi^2<=r0; 0<=psi<=1 is the same fixed
relative temporal cutoff. No actual integral is replaced by its envelope
as an identity. Equation (1) is the original full principal-column integral.

The same formula with x_sigma t_sigma,tan replaced by
(t_sigma)_a(t_sigma)_b defines each of the SIX entries of the unrooted
principal tensor. Those entries are not all zero or determined by the
two selected columns. The full physical curl-wave covariance also contains
the potential/curl remainder, its squares and correlations. It is NOT
identified with this principal tensor. Those complete-wave terms remain
in the existing initial/generic mean and residual budgets.

## 2. Moving maps cancel under the complete Haar integral

In the original auxiliary rectangle use Y_i=center+xi vr+ci(v-L/2)vt.
Any parameter-dependent torus translation can be included in center. The
full phase formulas in `core_phase_coefficient_budget` depend on the slow
point and v, not on xi. Thus the canonical entrance, coefficients and
homogeneous ODE are xi-independent, and uniqueness makes t_sigma so too.
The rectangle Jacobian is the constant J ci. The original cutoff is flat
before the boundary; its zero extension is smooth on the torus.

For a smooth periodic f(x,Y), any smooth center shift a(x) satisfies

    integral_T2 f(x,Y-a(x)) dY = integral_T2 f(x,Y) dY.   (2)

Differentiate this equality at each finite order. Terms from a'(x), a''(x)
and their products are total torus derivatives and have zero integral.
They do not vanish pointwise. Equation (2), not an estimate that discards
chain terms, justifies differentiating (1) at fixed v on the fixed interval
[0,L]. Labels, L, ci, frozen c0, actual tilt and reference spectral data
remain fixed under continuous differentiation, as in the parent. The slow
partition/root factors are outside the per-label unrooted column (1).

The common-torus map Y_i=Jmatrix^Delta H has degree 14^Delta, but preserves
NORMALIZED Haar average: a nonzero integer character stays nonzero under
(Jmatrix^Delta)^T. Thus no deck-copy count multiplies (1). Changing the
physical torus representation neither changes H nor adds a center-jet
input. These conclusions rely on the original full torus integral and
smooth flat extension. A truncated nonperiodic integration domain would
have boundary terms. A varying length or a nontranslation chart would
require its own endpoint/Jacobian terms and is outside this fixed-label
statement. The retained common-chart-to-path jet conversions are still
required for the parent U_j; no new conversion is silently inferred here.

## 3. Explicit all-order column estimate, retaining the physical Haar scale

The parent derives positive polynomials U_j(S), uniform over the original
accepted labels and both signs, such that every raw mixed derivative of
total order j obeys |D^I t_sigma|<=P(v)U_j(S) on 0<=v<=L. These are
derivatives of t, not t/P: all Gaussian derivative terms are already in U_j.
The same fixed positive actual tilt and spectral floor give the explicit
constant c>0 in `gaussian_mass_foundation`, with

    P(v)^2 <= exp[-c(v-L/2)^2/L],
    integral_0^L P^2 dv <= sqrt(pi L/c).                (3)

For raw convolution (U*U)_j=sum_i binom(j,i) U_i U_(j-i), Leibniz and
(1)--(3) prove, for every fixed slow multiindex |I|=j,

    |D^I(H_sigma/s0)|_2
       <= [(J/2)r0(2r0/L)/s0] sqrt(pi L/c)(U*U)_j
       = J sqrt(pi/c)(U*U)_j.                          (4)

The same right-hand side bounds any one normalized unrooted principal
tensor entry. In particular the full radius r0^2/sqrt(L) cancels ONLY
after correct normalization. Dropping it earlier would change the physical
root scale. The temporal cutoff remains in the exact column; only its
value bound psi^2<=1 is used in the inequality. Its slow derivatives are
zero in the fixed-v integral, because L and its fixed profile are frozen.
Its auxiliary derivatives vanish from H only AFTER averaging as in (2).
No claim that temporal-cutoff derivatives vanish from the actual wave is made.

Let G=2 max(1,1/c). Since pi<4 and G^2 c>=4, G>=sqrt(pi/c). The matrix
infinity norm sums TWO sign columns. For the parent's factorial exact-degree
jet sum J_j in the d=5 coefficient chart, pad the slow count by d and obtain

    H_j(S)=2*(d^j/j!)*2*G*(U*U)_j(S),
    J_j(Hbar)<=H_j(S).                                 (5)

The first factor 2 is the row sum of the two columns; the second is the
upper bound for the original Jacobian, not the angular factor. The count
sum_(|I|=j)1/I!=d^j/j! retains Taylor normalization. Auxiliary derivatives
of Hbar are zero. Every fixed H_j is a finite positive polynomial with
NO new epsilon, radius, edge or S loss from Haar normalization. The
already derived homogeneous U_j still carries all its legitimate losses.

The production factory derives c from the original ACTUAL fixed positive
tilt and original lambda_min/lambda_max using the prior curvature formula.
A mere upper bound on tilt does not imply a positive uniform c: c tends
to zero as the actual tilt tends to zero. `tilt` must therefore identify
the already fixed positive tilt, at most the retained upper bound; it is
not a new tilt choice, extra accepted-band gate or sampled curvature fit.

## 4. Substitute H and W into the same canonical root recipe

`root_join` now supplies BOTH (5) and the parent's

    W_j=(d^j/j!)(cutoff*U)_j,                           (6)

from the SAME canonical homogeneous recipe. Equation (6) keeps every
actual-wave cutoff derivative, even though those terms integrated out in
(2). All target, full curl-operator and loghalf callbacks and original
K0, beta, R, mass, width and trace are preserved. This changes only
positive NORM recipes. It does not solve for different physical roots,
reuse inherited-family roots with canonical columns, alter the value/cone
domain, or apply edge unweighting twice. The existing differentiated inverse,
root, signed-linearization and complete curl recurrences can use (5)--(6).

Thus independent H_j and W_j norm inputs have been removed RELATIVE TO
the full original phase/background/canonical ODE and fixed Gaussian/geometry
premises. This is a quantitative version of the previously qualitative
finite-jet statement, not a new canonical-primary existence theorem.

## 5. Evidence and the next lower obligation

Exact manufactured polynomial integrals retain all six principal tensor
entries, both sign columns, nonconstant Gaussian-factor derivatives,
the temporal cutoff and the physical radius. They are enclosed through
order3 on the whole parameter interval by rational coefficient bounds.
Separate periodic translation integrals test all moving-center terms through
order3 and contrast a domain with boundary. Original-factory controls bind
the unchanged accepted tail, derive the fixed-tilt curvature and feed the
joint H/W inputs through the original inverse, root, signed and full curl
norm recipes. These controls are manufactured and do not certify actual
background suprema or a complete physical covariance. Eight historical
records and all their declared source hashes are preserved and rebound.

Next make the remaining full curl-operator coefficient input D_j explicit
from the already derived normal reciprocal and actual cylindrical/graph
operators, including every coefficient derivative and normal/potential
factor; or identify its exact blocking inequality. Actual full background
jet validity, weighted target/freezing/value/cone premises and original
physical domains still apply. Full-wave covariance means must retain
axial velocity, pressure primitives/moments and reserved supports. Viscosity,
cutoffs, every interaction, the full physical residual and infinite
realization remain separate obligations. No unconditional pulse bound,
complete residual improvement or blow-up theorem is asserted.

Reproduce with `python test_audit.py -v` and then
`python audit.py --verify-record evidence.json --output /fresh/path.json`.
Install the declared SymPy requirement if necessary. The analytical
all-order derivation is (1)--(6); finite checks do not replace that proof.
