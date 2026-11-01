# Q70: arbitrary radial order through the shape and earlier annulus

For the actual PREMOD field from Xi=110 through Xsep and the earlier
annular restoration/join to R=X_R, there are explicit bounds at every
finite radial order with angular order <=6. The full state including p1
and p2 consequently has bounds through angular order five, independent of
how many radial derivatives are taken. This uses the inherited exact root
and quantitative sixth-order estimates, not an all-order small-root claim.
It complements Q69's outer prefix; it does not extend Q69's unrestricted
angular order to this earlier interval.

Use the real angular algebra norm ||v||_6=sum_(j=0)^6 sup|v^(j)|/j!.
It is submultiplicative by Leibniz. Throughout let
B_j=162·3260^j(j!)³ be the Q63 raw derivative bound for the fixed
squared-flat step sigma. The same sigma occurs here, but the affine
arguments and the bump normalizations differ.

## Xi to Xsep: exponentiate the mixed log jets at their actual value

Use [the sixth native angular enclosure](../native_sixth_angular/README.md),
section 3. In z=log(X/Xi), T=2000 B0, 0<=z<=T,
s=sigma(z/T), ell=log(C E(Xi)), G=U(Xi),

    log E=-log C+z/10+(1-s)ell+s log f,   U=G,   f=1/(1+eta²).

The inherited bounds are ||ell||_6<=J6<B0/100, ||log f||_6<6,
||G||_6<9, ||E||_6<=Eb=exp(204B0)/C. Consequently for r>=1

    ||partial_z^r log E||_6/r! <= H_r
         = 1_{r=1}/10 + (J6+6) B_r/(r! T^r).
    D_0=1;   D_r=(1/r) sum_(k=1)^r k H_k D_(r-k).
    ||partial_z^r E||_6/r! <= Eb D_r.                       (1)

This is the differentiated identity E_z=E(log E)_z in the angular
algebra; the zeroth coefficient is the actual angular norm of E.
It does not exponentiate an enormous high-derivative bound. The U bound
is 9 at r=0 and zero thereafter. A real lower bound for E on this segment
is exp(-J6-6)/C by the value formula and z>=0. All width losses remain
in T^(-r); there is no assumption of radial analyticity.

## Xsep to R: relative swirl edits and additive axial edits

Write y=log(X/R). Before the bumps E0=P f exp(y/10); up to y=-8
U=4eta+g, with ||g||_6<=gb=3e-18. On [-8,-7] replace g by
g(1-sigma(y+8)), and continue the latter as zero. On the annulus
[-6,-5] use [the exact original columns](../axis_matching_audit/README.md):
E=E0(1+sum e_l beta_l), U=4eta+sum u_l beta_l, with three e's and two
u's. Outside their supports these edits are zero. The accepted earlier
root has ||e_l||_6,||u_l||_6<=b=5.72e-14. Its defining quadratic map
depends on eta; this bound already charges those derivatives.

Crucially beta_l=sigma'((y+6-center_l)/w+1/2), w=.08, with NO 1/w
prefactor. Let L_0=64 and L_r=B_(r+1) w^(-r) for r>=1. These bound
raw beta derivatives, including the flat extensions. Let V_0=1 and
V_r=B_r for r>=1, bounding the width-one restoration. For every r>=0,

    ||partial_y^r E||_6 <= 28P [10^(-r) +
             3b sum_(k=0)^r binom(r,k)10^(-(r-k)) L_k],
    ||partial_y^r U||_6 <= 8·1_{r=0}+gb V_r+2b L_r.           (2)

Indeed ||f||_6<=28 and exp(y/10)<=1. Leibniz acts on the relative
swirl factor E0; omitting that product would miss cross derivatives.
The term 8 is the norm of 4eta. Positive sums cover gaps and the ideal
segment after the join. The real floor is
E>=(P/2)(Xsep/R)^(1/10)(1-192b)>0 since |sum e_l beta_l|<=192b.
No factor of the later compact modulation frequency is needed here.

## Physical derivatives and full, partial pressure memory

Call the raw bounds above F_r for either field (multiply (1) by r!).
On either segment with left edge X0, for j<=6 the normalized (i,j)
field coefficient is at most

    X0^(-i)/i! sum_(r=0)^i [i,r] F_r,                         (3)

where [i,r] are unsigned Stirling numbers of the first kind. This is
the falling-factorial identity for differentiating a function of log X;
both log coordinates differ from log X only by a constant. Take the
maximum of the two segment bounds to cover the hull. The original fields
match smoothly: the shape exit is flat onto the ideal power, the
entry into restoration retains g, and all bumps end flat.

For every X in [Xi,R] use the genuine matched native trace at R from
[Q68](../core_current_head_q68_outer_prefix_traces/README.md), eqs. (1),(2)
with B=R. For each of the five primitives write
P_row(X)=P_row(R)-integral_X^R density(s)ds using its OWN pre-modulation
density U, sqrt(2X)E, sqrt(2X)EU, U²-E²/2 or E²/(2X).
This use of the right endpoint is justified by the earlier exact five-row
match and propagation on [-5,0], not by assuming the finished heat trace.
If t_j bounds its normalized eta derivative at R, the pure eta row is
bounded by t_j+(R-Xi) sup g_0j; a positive radial row is bounded by
sup g_(i-1,j)/i in normalized coefficients. Form g by ordinary
convolution from (3) and the geometric factors on X>=Xi as in Q67.
Thus partial moment and pressure increments inside a bump survive even
when the final five increments cancel. Differentiation under this finite
integral uses the smooth matched pieces and fixed endpoints.

Use Q67's reciprocal recurrence with the smaller of the two displayed
positive E floors. Its full p2 expression still includes
X(4Aphys eta Pi-(1-eta²)Pi_eta); p1 retains I_eta and J_eta as in Q69.
For j<=5 the shifted primitive row j+1 is <=6. Local E_X,U_X use one
more radial derivative, available at every order in (1)--(3). This gives
a quantitative state enclosure throughout [Xi,R], not merely existence
of a different finite constant for each derivative. In particular no
independent unknown pressure trace at Xi or Xsep is required here.

Remaining: quantify the prescribed continuation before Xi, and bound
angular orders beyond six for the shape data and earlier joining target
when required by the output order. Their smoothness alone gives no
quantitative norm. Then close the normalized repair defect and check the
physical smallness inequalities. Uniform actual pulses, the pressure tails,
complete residual improvement and blow-up remain unproved.

Run test_audit.py and audit.py --verify-record evidence.json --output a fresh
path. The small-order algebra checks and provenance replay are distinct
from the analytical estimate.
