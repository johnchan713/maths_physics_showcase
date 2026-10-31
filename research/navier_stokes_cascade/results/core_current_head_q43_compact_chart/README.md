# Q43: same-weight compact operator chart conversion

Parent `c612d8415503d4f2fc327e307be2eee39adfe985` (Q42).
On 31 October 2026 its focused CI 41157417054 passed; broad CI
41157417019 was still running. Historical evidence is preserved.

Q42 supplies RAW native-U compact constants. Its statement that these can
be passed to Q40 is conditional on a matching input norm; it does not itself
establish the normalized R-chart operator inequality. This checkpoint
supplies that missing sufficient conversion. It leaves physical fields,
width, shifts, pressure memory and supports unchanged.

## Finite coordinate-jet bounds

Write U=R^d, 1/3<d<5/2. On the same fixed positive hull choose positive
rational enclosures Rmin<=R<=Rmax and Umin<=U<=Umax. They must enclose
corresponding points of this SAME power chart; mere positivity of four
numbers is not certification of that geometric premise. Four other
coefficient coordinates are unchanged. Cylindrical scalar components are
re-expressed, not transformed to coordinate-basis tensor components.
For a chart x=y^a,

    |partial_y^j x|/j! <= H |binom(a,j)| l^(-j), j>=1,

where x<=H and y>=l. Define A(t)=4t+sum_(j>=1) a_j t^j using this bound.
Only finite jets are used; no analyticity of the physical fields is assumed.

Let N_n be the sum of absolute derivatives divided by their multiindex
factorials through n, and M_n the raw derivative maximum through n.
If N_n(F)<=K w at the input point, each normalized Taylor coefficient is
at most K w. At outer degree k the sum of coordinate monomials is bounded
by A(t)^k (multinomial coefficients only increase that positive sum).
Thus a sufficient normalized-input to raw-output factor is

    Lforward_n = n! sum_(j=0)^n [t^j](1-A(t))^(-1).       (1)

For RAW input M_n(F)<=K w, the outer coefficient includes 1/alpha!.
The multinomial theorem instead gives exp(A(t)); hence a sufficient factor
for raw-input to normalized-output is

    Lreverse_n = n! sum_(j=0)^n [t^j]exp(A(t)).           (2)

The factorial in (2) is conservative and also bounds raw outputs. Both
formulas are positive finite rational recurrences. Using exp(A) in (1)
without first controlling raw input derivatives would be unjustified.

Use exponent 1/d, H=Rmax,l=Umin in (1), and exponent d,H=Umax,l=Rmin
in (2). All derivatives are evaluated at the same physical point. The
weight w=zeta delta^(-v) is carried as a pointwise right-hand side and is
never differentiated. Thus its pole, both edges and epsilon order remain
unchanged. These constants are uniform over the fixed hull, not as its
lower endpoint tends to zero or as n tends to infinity.

## Compose with the actual shifted operator

For either Q42 operator O=T or A,

    N_R,n(O F) <= Lreverse_n Kraw_(n,v) Lforward_n N_R,n(F), (3)

in the weighted operator sense (input/output weights at their respective
points). First convert the input bound at EVERY shifted input point, apply
Q29's integral theorem on the same whole support, then convert the output
jets. This is not a claim that nonlinear chart changes commute with the
shift. The shift remains in U; neither large M nor shifted endpoints are
differentiated in this argument. Formula (3) adds no derivative loss.

`ChartedCompact` wraps the exact Q42 provider and multiplies KT and KA
by the two factors. It preserves order, exact pole and the original Q40
time-inverse/S input. The latter has its own retained norm convention and
is not an operator output of this radial chart. Q40's existing factorial
padding is left intact: additional overestimation is safe. Q42's cutoff
ceilings m+2 for fast output m, and m+7 for Q41 output m, remain sufficient.

Exact tests compare low-order Taylor coefficients and actual mixed jets,
including a nonlinear power chart, and check the wrapped interface. These
are finite controls of the derivation, not selected-field certification.
Frozen Q42 evidence and its transitive source bindings are replayed.

The actual profile jets, plateau gap, raw multiplier bounds and geometric
enclosures still require verification for the original fixed choices.
Next bind those inputs and the SAME initial A1 tail to the all-order source
graph; operator boundedness alone supplies no M^(-p) tail improvement.
Full current source bounds, uniform pulse estimates, complete residual
improvement and blow-up remain unproved.
