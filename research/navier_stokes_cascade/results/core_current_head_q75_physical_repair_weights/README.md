# Q75: the compact repair in the physical axial-moment weights

Q74 bounds the five repair coefficients. Here we transfer that bound to
the two actual leading shear multipliers of Q57 on every accepted band,
and isolate its cost in the axial estimate. This is a decomposition of
the WEIGHTS only. The tensor being multiplied remains that of the original
wave (or the same conditional Q61 comparison); it is never recomputed
after deleting the repair. No new root, partition, cutoff or band is used.

## Two exact perturbations, not four independent errors

Use [Q58](../core_current_head_q58_leading_weight_enclosure/README.md)'s
notation: q=R^d for the radial integration coordinate, D=1/(d q),
s=s(Z,T), Cchi=s mu the complete moving pressure-cutoff moment,
X=R²/(2s), eta=Z s^(-D0), A=1/2+h. On the first reserved patch write
y=log(X/rhop). Let P_E=e(eta) sum_E c_l(eta) beta_l(y) and
P_U=e(eta) sum_U c_l(eta) beta_l(y) be the original additive changes.
Subtract only these smooth, compact fields from the final leading fields
to define the comparison weights. All five bumps are the width-.1,
unit-mass bumps of Q63; no hard window at either patch end is inserted.

At fixed Z,T, R_q=D R and y_q=2D, eta_q=0. Thus the changes of
F=s^(-A)E/R and G=s^(-A)U are explicit. In Q58's order the two
diagonal weights have ZERO change; the two shear changes are

    delta w3 = 4 R D s^(-A) partial_y P_U,
    delta w4 = D s^(-A)/R [
       (2Cchi-R²)(2 partial_y P_E-P_E) - 4 Cchi P_E ].        (1)

In particular even a locally constant swirl edit has a derivative
contribution through 1/R. The moving pressure coefficient is not
frozen when the physical Z derivative later hits (1). At the level of
leading weights, w3 is supported on the two axial bumps on this patch
(the native axial field there is zero); delta w4 uses only the three
swirl bumps. It is delta w4, not the full w4, that vanishes in their
gaps. Pullback moves those supports with s. The estimates below
differentiate the actual smooth translated functions, not indicators.

## Uniform finite constants in the accepted chart

Fix output n. For each Q57 integration count let j=n+1+count and
l=n+1. Its weight derivatives are a subset of total order <=j with
at most l combined Z,T derivatives. Denote by E_k a uniform bound on
|e^(k)|/k!, and by D_k the Q74 bounds on N²|c_l^(k)|/k! (maximum over
bump labels). Define L_k=sum_(i=0)^k E_i D_(k-i).
These are fixed profile quantities supplied by Q73--74 including the
earlier pressure memory, scale multipliers and lambda loss. The dependence
on the fixed modulation frequency is displayed separately from physical
epsilon; e and the Q74 bounds can be chosen independent of N in that
frozen-family comparison.

Use r<=R<=H on the original band hull, u=min(1,r³), delta=u/1024
and the absolute normalized jet norm of Q58 at radius delta. Work in the
downward-closed truncation just described. From Q58,
||R||<=2H, ||1/R||<=2/r, ||D||<=6/u, ||s^(-A)||,||s||<=8,
||eta-eta0||<=3. Also
||y-y0|| <= (2/d) sum_(k>=1) (delta/u)^k/k + ||log s-log s0||
< 6/1023 + 1/8 < 1.
Only the logarithmic coordinate is expanded here. Taylor substitution
of the smooth compact bump is finite, with derivatives at the actual
real center. No analytic continuation of a flat bump is made.

For m=0,1 let
    V_m=(sum_(k=0)^l 3^k L_k) (sum_(i=0)^j B_(i+m)/i!),
where B_i=162·32600^(i+1)((i+1)!)³ is Q63's raw unit-bump bound.
Because eta has no q increment, no order above l occurs in the first
sum. Thus each pulled-back partial_y^m P has jet norm <=V_m/N².
At a given center at most one bump contributes even to differentiated
jets; at a support endpoint its jets vanish. The same bound applies to
the sum without multiplying by the number of bands. Submultiplicativity
in (1) now gives upper bounds for raw derivatives through this stencil:

    D3= j! delta^(-j) (4·2H·6/u·8) V_1 /N²,
    D4= j! delta^(-j) (6/u)·8·(2/r) /N² times
        [(16|mu|+4H²)(2V_1+V_0)+32|mu| V_0].                (2)

All geometric bounds use the fixed positive hull, not its behavior as
r goes to zero. The full absolute jet norm includes derivatives of Cchi,
including mu s_Z. By linearity Q74's leading coefficient may be pulled
through (1): substitute c=-Bmat^-1 T V/N² there, with mixed weight
remainder O_(j,l)(N^-3) on the same bumps. The physical shear weights
in (1) concern the slow background, not the principal amplitudes.

## Which frequency gain can be used in the axial test?

In [Q57](../core_current_head_q57_primary_tail_budget/README.md)'s
shear estimate replace its weight-derivative input D by (2), keeping
its original label-square amplitude input H_amp, the complete Z product
derivative, temporal divisor and cutoff B1. In particular the repair-only
part of the sufficient D H_amp² factor is bounded by
max_(counts)(max(D3,D4) H_amp²), with **no unknown repair-weight jet**.
This proves that this part of the weight meets the fixed-order polynomial
weight requirement uniformly on the accepted physical charts (a constant
is a degree-zero polynomial). Whether the accompanying H_amp, propagator
and inverse satisfy their required uniform bounds is a separate test.
No N^-2 improvement of the principal tensor itself is inferred.

There are two distinct small parameters. A direct whole-line bound with
low-order amplitude inputs would instead charge this term as
epsilon² K_n(S)/N², with all low-order constants, including the divisor
and the weight cutoff, in K_n. For such an upper bound alone to fit a
proposed allowance epsilon^(3-5kappa) P_n(S), one must check

    K_n(S) <= N² epsilon^(1-5kappa) P_n(S).                   (3)

A fixed N is not a power of the shrinking physical epsilon. Equation (3)
is the exact comparison of the two positive budgets, not a counterexample
to the signed integral. If (3) is unavailable retain the Q57 oscillatory
estimate for this compact term; do not erase it as a higher physical
order. Unlike Q61's Borel error it has no automatic epsilon² factor
inside the weight. The exact moment matching after the patch does not
set the interior shear multipliers to zero.

Next enclose the FULL leading weights and generator on the rest of the
physical hull, in particular the exterior beyond Xrep up to the terminal
stress edge Xb, and verify the simultaneous low-order physical inverse
and propagator gates. The axial rate remains conditional; the separate
radial pressure tail, uniform actual pulse, complete covariance/mean and
viscous residual improvement, and blow-up are unproved.

Run test_audit.py and audit.py --verify-record evidence.json --output a
fresh path. The exact tests check the physical pullback and frequency
comparison, not the continuum hypotheses.
