# Q74: the terminal modulation defect and its five-row repair

At the entrance to the FIRST reserved repair patch the five incoming
moment defects are O_k(N^-2), with an explicit quadratic phase average
and O_k(N^-3) error. This sharpens the previous O_k(N^-1) absolute
integral bound at that endpoint. It does not give the same expansion
at an arbitrary cut through the active loop. All estimates are relative
to the retained profile, low-order acceptance and floors at the original
N=1+floor(H^32); the constants need not be small on a physical band.

## Use the phase in which the primitives have zero mean

Write y=log X and let [a,b] run from a preserved entrance to the repair
entrance b=log(rhop), containing the whole modulation. The actual fields
are E exp(A/N), U+B/N at phase phi=Ny, where A,B denote the unique
period-one, zero-mean primitives in [the realization](../stress_realization_audit/README.md)
section 4. They vanish, with their slow derivatives, on the end collars.
All averages below are in phi AFTER the moving inverse; averaging at
uniform original theta would use the wrong measure. Full incoming Pi
and Pi_eta enter through p2 in that loop. Only the differences of moments
start at zero; the moments themselves do not.

For any fixed k use the angular algebra norm sum_(j=0)^k
sup|partial_eta^j f|/j!, taking suprema in y and phase as needed.
The Q73 regional bounds and Q66,Q65,Q64 in that order bound all required
fixed-phase slow jets of E,U,A,B by finite operations at the retained
parameters. In particular Q64 (4) at fast order zero supplies their
slow phase-primitive bounds. This is pre-substitution differentiation:
d/dy of A(y,eta,Ny) is not bounded by the slow derivative.
Q64 (5) and the ordinary product rule then give the modulated field
mixed bounds (radial order r costs up to N^(r-1) in the differences).
The estimate below needs only two slow radial derivatives, regardless
of angular order k; the phase has no eta dependence.

## Endpoint expansion, including pressure

Set tau=1/N, w=X sqrt(2X). At fixed y,eta,phi expand the five EXACT
y-densities with E(tau)=E exp(tau A), U(tau)=U+tau B. In the order
M,I,J,S,Cp their changes have linear and quadratic coefficients:

| row | L (coefficient of tau) | Q (coefficient of tau²) |
|---|---|---|
| M | X B | 0 |
| I | w E A | w E A²/2 |
| J | w E (UA+B) | w E (UA²/2+AB) |
| S | X(2UB-E² A) | X(B²-E² A²) |
| Cp | E² A | E² A² |

These follow from the separate densities XU,wE,wEU,X(U²-E²/2),E²/2.
In particular the J cross term survives, and the last row has no X.
For a row let H_k dominate the angular algebra norm of its entire
density on the auxiliary circle |tau|=2, uniformly in y,phi.
This is explicit: use the original density formula with positive sums,
E replaced by ||E|| exp(2||A||), U by ||U||+2||B|| and the supremum
of its geometric weights. Cauchy's coefficient bound H_k/2^j gives
a tail after degree two bounded by H_k/N³ for real N>=1.
No holomorphic slow cutoff extension is involved.

Every L has mean zero. For zero-mean g let Pg be its periodic
zero-mean phase primitive. Integration from 0 to phi and subtraction
of that primitive's mean give ||Pg||<=2||g||. It commutes with slow
derivatives and preserves vanishing on the collars. Along phi=Ny,
g=(d_y(Pg)(y,eta,Ny)-(P partial_y g)(y,eta,Ny))/N, where d_y is
the TOTAL derivative. Integrate twice for L and once for Q-<Q>.
The boundary terms vanish because the interval contains the complete
support and its endpoint jets vanish. Thus for each row, in the same
angular norm, the terminal modulated-minus-original difference satisfies

    Delta = N^(-2) V + error,     V=int_a^b <Q> dy,
    ||error|| <= (b-a) [4||L_yy||+4||Q_y||+H_k]/N³.        (1)

Here the factor 4 on Q_y includes subtracting its phase average.
The moments retain their genuine inputs throughout; (1) is only an
increment. At interior cuts the primitive boundary value must be kept.
For real fields the pressure component V_Cp=int E²<A²> dy is nonnegative,
strictly positive at an eta if the swirl primitive is nonzero on an open
set there. A blanket o(N^-2) claim would therefore require additional
cancellation or degeneracy. This sign criterion does not infer strict
positivity for every eta or estimate the tail of physical wave pressure.

## Normalize before comparing the repair root

Let T act on ordinary moment differences by first forming
m=M/(rhop e), j=J/(sqrt(2) rhop^(3/2) e²),
i=I/(sqrt(2) rhop^(3/2) e), s=S/(rhop e²), p=Cp/e²,
then returning (m,(m-j)/lambda,i,-s,p). Here e is the ORIGINAL
eta-dependent reserved swirl amplitude, e=E(rhop,eta) on the unchanged
power patch; lowercase p here denotes only the normalized fifth row.
The actual target is z=-T Delta. Thus (1) yields

    ||z + T V/N²||_k <= ||T||_k max_row error_constant / N³. (2)

For an explicit operator bound, replace every reciprocal scale in the
five normalized rows by its angular algebra norm, using the sum of the
m and j bounds divided by the positive lambda for the second row.
Use Q73's E and the retained positive floor at rhop in Q67 (4) to
compute e^-1 jets, and convolve for e^-2. This includes derivatives
of the normalizers. Neither the exact zero V_M nor disjoint repair
bumps remove the 1/lambda cost of general incoming J. Equations (1),(2)
are fixed-parameter estimates with constants independent of N when the
native data and loop are held fixed. They do not choose a new frequency.

For the same accepted solution of Bmat c+Qmat(c,c)=z, the constant
matrix has ||Bmat^-1||<=1000 and the pointwise inverse Jacobian <=2000.
Use (2) to bound each normalized z jet by Z_j/N², taking the RHS
constant plus ||TV||_j since N>=1. The low-order acceptance already
ensures beta q ||c||_infty<1/2; hence ||c||<=2 beta ||z|| pointwise.
For j>=1 remove the two linear-in-highest-derivative terms as in
[the fixed-frequency repair](../modulation_repair_all_order/README.md).
Explicit bounds ||c^(j)||/j!<=D_j/N² follow with D_0=2000 Z_0 and
D_j=2000(Z_j+10^5 sum_(l=1)^(j-1) D_l D_(j-l)), using N^-4<=N^-2.
No contraction in a high-order norm is used. Moreover the exact identity
c=Bmat^-1 z-Bmat^-1 Qmat(c,c) shows

    c= -Bmat^-1 T V/N² + O_k(N^-3).

The quadratic term here is O_k(N^-4) by the preceding bounds and
Leibniz. Q63's fixed bump products now give repair field derivatives
at any fixed mixed order with N^-2 in place of N^-1. The oscillating
modulation itself still has the Q64 radial powers. At the patch exit all
five moments match exactly by the accepted equation, not just to the
order displayed here. No axis datum, width, support or gauge has changed.

Next insert the resulting finalized-profile bounds in the actual
physical-band pulse/propagator and axial moment inequalities. Neither
their smallness nor the separate radial pressure tail, complete covariance
mean corrections, full residual improvement or blow-up is established.

Run test_audit.py and audit.py --verify-record evidence.json --output a
fresh path. Manufactured controls are not a continuum proof.
