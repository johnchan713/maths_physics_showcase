# Q51: the radial component of the primary curl and the Q36 target

Q50's necessary value test for an isolated leading viscous term did not
use all the inequalities already available on the branch. The frozen
[initial target rate gain](../core_initial_target_rate_gain/README.md)
retains the Haar rate beta=2-2kappa through the linear two-edge inverse,
rather than downgrading it to alpha=2-3kappa. Thus under the existing
physical norm hypotheses Q50's test was already met for eta=kappa.
The value test alone is not a strict-gain proof. Here a component
estimate improves the **full finite-jet initial signed target** to order
two on the canonical curl route, with the unchanged Q36 polynomials.

## One curl factor has no graph loss

Use the [complete cylindrical curl](../core_complete_curl_operator_budget/README.md)
of C=h cross A, h=normal/|normal|^2. A includes the cutoff and root.
In the frozen native band its components are exactly

    B_r=-epsilon partial_Z C_theta,
    B_theta=epsilon partial_Z C_r-D_r C_z,
    B_z=D_r C_theta+C_theta/R.

No radial graph or cylindrical connection enters the radial component.
The Z derivative acts on both h and A. To propagate this fact through
the actual finite norm chain, bound the operator for **B_r/epsilon**
in its own right. Its only first-order coefficient applies to partial_Z A
and its zeroth-order coefficient contains partial_Z h. For raw mixed
jets of degree j the sum of the two induced matrix norms is bounded by
2(h_j+h_(j+1)), using the complete-normal raw lift seminorms. Convert
to a factorial exact-degree bound by 5^j/j!, then apply the original
positive R(U) composition. Z is unchanged by this radial chart: there
is no extra g on the Z derivative. No estimate for bare graph derivatives
of the other two components has been used.

This program is coefficientwise dominated by the historical canonical
normalized-U curl operator: its first convolution contains the order-zero
unit direction in V_U, and its shifted lift convolution contains that in
V_R. All remaining terms and the chart composition are positive. Therefore
the original root recurrence for complete_curl(j), including the next
amplitude jet and half-weight pole 6j+9, bounds BOTH epsilon^kappa B
and B_r/epsilon, separately. This is not an inference from an arbitrary
norm bound on epsilon^kappa B. The adapter checks domination against the
actual root polynomial through n+3 and refuses a smaller archive. The
lift must enclose the same full normal; matching polynomial classes or
R endpoints is not evidence of that physical hypothesis.

In the [coarsest common atlas](../core_primary_aggregate_atlas/README.md),
the stripped radial normalization adds epsilon_i/epsilon_* <=1 to the
usual velocity transfer, while the other factor keeps the old normalized
curl transfer. Linear jet conversions, finite support separation and
normalized Haar averaging are unchanged. Thus the SAME aggregate b_N
bounds the two separately normalized factors; no growing label count is
introduced. The exact real quadrature and matched principal Haar stresses
give, for e=theta,z, the complete selected curl defect term

    (epsilon_*/2) sum_i k_i^-2 Haar(Btilde_i,r Btilde_i,e).

Since k_i^-2<=epsilon_*, its raw through-N coefficient bound is
2^(N-1)b_N^2 times epsilon_*^(3-kappa) zeta delta^(-12N-18).
This improves only entries with a radial factor. The full W for pressure
and axial products still has order one. The higher background stress in
the selected defect has order three; the Q34 R transfer carries the same
weight, and 3-kappa <= 3. Hence the original Dcurl+S polynomial suffices
at 3-kappa. No assertion of physical cancellation in the curl square is
required.

## Consequence for both retained projector channels

Reinsert the improved D in [the original Haar equations (6)-(10)](../core_initial_haar_row_budget/README.md).
In the two radial fluxes b times the slow tangential velocity and V or G
times the slow radial field start at order two. The fast radial/tangential
correlations start at 3-2kappa. The only previously worse-than-two term
was the selected flux defect. Slow time and viscosity are order two; the
axial derivative costs an explicit epsilon multiplying the full order-one
pressure/axial flux. Haar removes graph derivatives only after all
products are formed. These estimates yield order two in BOTH complete
Haar rows with their original positive coefficients and doubled pole.
The order-three projector moments retain both terms of partial_Z(c_rho P)
and the distinct angular moment. Pullback, rho-squared projector,
two-edge primitive, reciprocal-weight and chart return preserve order.
For the actual Q36 target, on precisely those lower-input premises,

    J_n(Xi_e / (epsilon^2 zeta)) <= V_e,n(S) delta^(-2p-3n),

where V_e,n is the unaltered native target polynomial (exact jet degree n).
This is a uniform-in-band analytic implication at each fixed finite order,
not a numerical evaluation or proof of the underlying field enclosures.

Relative to Q50, this supports eta=3kappa for the target, and, using the
original linear signed principal-profile estimate with the needed jets,
order 3-3kappa for its isolated leading viscous term (epsilon times three
graph costs and the polynomial temporal floor). This is NOT the rate of
the complete fast difference. In particular the original fast difference
in the input to the new fast map still has only order alpha; a better
signed target alone does not improve that minimum. Q49's original axial
tail at 2-kappa already exceeds alpha. The later shifted axial remainder
and the physical radial pressure tail have different sources.

Next bind the sharpened target to the complete first-real final rows,
particularly the later A1 and A0 source and support estimates. No uniform
actual pulse enclosure, complete residual improvement or blow-up follows.
