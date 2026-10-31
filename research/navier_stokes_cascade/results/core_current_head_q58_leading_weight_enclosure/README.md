# Q58: enclose the four axial moment weights from the leading background

Q57 left the four assembled slow-weight norms as inputs. For the actual
fixed cutoff and positive annular hull these norms in the oscillatory
estimate can be replaced by one leading-profile seminorm and explicit
geometric constants. The full-background error uses only low jets here.
The replacement is **only in the multipliers of Q56's moment**; the wave,
roots, pressure and earlier remainder estimates still use their true state.

## Pressure shape and exact multiplier algebra

Let χ₀=χ̂(R²/(2s)) with the original plateaux [c₀,c₁] in X. On the
whole fixed integration hull of Q56, change variable in its integrated
cutoff derivative:

    cχ = ∫ (R²/2) ∂Rχ₀ dR = s μ,
    μ = ∫ X χ̂′(X)dX = c₁ − ∫_(c₀)^(c₁) χ̂(X)dX.

Thus |μ|≤c₁+(c₁−c₀)sup|χ̂|; the moment uses no derivative norm of χ̂,
even if its transition is very steep. Its Z derivative is μ s_Z,
not zero. This argument uses the full transition, not a truncated
integration interval. The axial output cutoff B₁ is a different
multiplier and still requires Q46's derivative bounds.

Write F=V/R and use U=R^d, D=R_U/R=1/(dU), C=sμ. Direct substitution
of w=C−R²/2 and a=2Fw gives Q57's four dU weights in order:

    D R²,   D(C−R²/2),   2R G_U,   (2C−R²)F_U − 4C D F.       (1)

In particular the terms proportional to R R_U F in the last weight
cancel. Treating V as F, or dropping the derivative of w, spoils (1).

## Charge the Borel correction at low order before oscillatory estimation

Use the [finite-head background estimate, §3](../core_background_finite_head_budget/README.md):
the full chart F,G minus their leading terms have raw R,Z,T derivatives
through j bounded by ε² Jerr_j, on the same accepted hull and schedule.
This bound includes all coefficient and cutoff corrections. Transfer
these finite jets to U by the positive chart of Q43. Products in (1)
then give O₂ difference for either shear weight, through order n+1
using background errors through n+2. The diagonal weights do not change.

In Q57 the shear input is ε² times a shifted integral of a Z derivative
of weight times I QP, where P is a principal-amplitude product of order
zero. Use the direct whole-line bound here: at most n+1 derivatives hit
the weight, and n+5 raw amplitude jets suffice for the direct temporal
inverse. The output cutoff and transition floor cost no ε power.
Hence replacing the shear weights by their leading versions costs O₄,
smaller than Q56's O_(3−2κ) remainder. This argument does not put an
ε² correction through 99996 integrations by parts or assume that all
high-order Borel coefficients have already been enclosed.

## A uniform geometric factor at the requested order

Here is a coarse alternative to evaluating a multi-index composition
table at order 99998. Fix the actual hull r≤R≤H, r>0, u=min(1,r³)
and δ=u/1024. Write L_j for the maximum of the raw (X,η) derivatives
through j of the two finalized leading profiles E₀,U₀ on the required
image of this hull. The leading fields in (1) are

    F⁰=s^(-A) E₀(X,η)/R,   G⁰=s^(-A) U₀(X,η),
    X=R²/(2s),  η=Z s^(-D₀),  A=1/2+h, D₀=1/2−h.

For a finite Taylor jet use the sum of absolute normalized coefficients
times δ^(total degree). It is submultiplicative with truncation.
At an arbitrary chart center U≥u, the binomial series bounds the norms
of R and 1/R by 2H and 2/r: the absolute coefficient series of
(1+z)^a is dominated by (1−|z|)^(-|a|), with |1/d|<3.
The norm of 1/(dU) is ≤6/u. From the [coordinate log bound](../core_primary_target_cutoff_budget/README.md)
the increment of log s has norm at most
((1−32δ)^(-1)(1−64δ)^(-1)−1)/3 < 1/8.
Thus the norms of s and s^(-A) are ≤8, and 1/(2s) has norm ≤4.
For R² directly use exponent 2/d<6 and (1−1/1024)^(-6)<2.
The increments of X and η have norms ≤8H² and 3 respectively
(the latter uses |η_center|≤1, s_center^(-D₀)≤4 and |ΔZ|≤δ).
Truncating the ordinary Taylor formula for each smooth profile gives
norm ≤L_j exp(8H²+3); this uses only real jets at the center and
does not analytically continue the profiles. Consequently at order j+1
both F⁰ and G⁰ have norm bounded by

    B = 8 max(1,2/r) L_(j+1) exp(8H²+3).

After taking one U derivative, the through-j norm at radius δ/2
is ≤2B/δ, since k/2^(k−1)≤2 for k≥1. Products in (1) at this smaller
radius are bounded by substituting |R|≤2H, |R²|≤4H², |D|≤6/u,
|C|≤8|μ| and these field bounds, with plus signs. Multiply each
result by j! (2/δ)^j to bound every raw derivative through j.
This is a band-independent upper bound from the original geometric
and leading-profile data, finite for every fixed j. It can be left in
factored form; no unrepresentably large floating number is required.

Take j=n+99997 for Q57. The diagonal constants need no field data;
the shear constants request L through n+99998, with no high Borel
head or pressure-cutoff shape derivatives. The isotropic profile
seminorm used here is sufficient, not minimal. Actual numerical
enclosures of that leading seminorm and the principal-amplitude
envelopes at the Q57 indices have not been evaluated. Next bind those
amplitude inputs to the actual root and pulse ODE on every band. The
radial A₀ term, full residual and unconditional pulse remain open.
