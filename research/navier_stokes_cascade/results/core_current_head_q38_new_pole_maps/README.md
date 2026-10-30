# Q38: same-profile radial maps on Q37's enlarged pressure domain

Parent: `18fe97ad4956c3e25133013d178d90cb86d9e164` (Q37).
The historical same-choice map formula is now justified directly in the
normalized R-chart derivative norm and attached to Q37's actual requested
pole. This removes its independent Kp, primitive and bump callbacks.
It does not establish an actual continuum source bound or fresh means.
Parent focused CI 41134296419 and broad CI 41134296463 were both running
when inspected on 30 October 2026. Historical frozen records are preserved.

## 1. The two-edge kernel at the new pole

Use the unchanged shell L=log(Xb/Xa), width w and
zeta=exp(-w^2/x^2-4/(L-x)^2), delta=min(1,x,L-x).
Q36 supplies rational Llo<=L<=Lhi. For any nonnegative integer v put
r=min(1,Llo/4), a=w^2, b=4, k=max(0,ceil((v-3)/2)), and

    D(c)=(1/(2c)) sum_{j=0}^k binom(k,j)j! (r^2/c)^j,
    C(v)=D(a) exp(b/(Llo-r)^2)+D(b) exp(a/(Llo-r)^2)
                  +Lhi Eedge(v,w) exp((a+b)/r^2).             (1)

This is the existing sufficient sum of the left-edge, right-edge and
middle-interval primitive bounds. The right-edge branch requires exact
zero total integral; a forward-only estimate is insufficient there.
Constants remain exact rational exponential sums. No huge exponential is
numerically evaluated. The proof collar r is not a new physical support.

For a zero-integral source g in the R chart, dR=(R/2)dx implies

    |I0 g| <= (Rmax/2) C(v) K zeta delta^-v                  (2)

whenever |g|<=K zeta delta^-v. Both oriented integrals of g are available,
so the same Jacobian bound applies on either edge. The actual moving
annulus stays in the original common positive radial norm hull.

## 2. Through-order derivative bound without differentiating a bound

Assume the original smooth flat zero extensions and exact integral zero
hold as functions of the slow variables. Then for every slow multi-index b,

    partial_slow^b I0 g = I0(partial_slow^b g),
    integral partial_slow^b g dR = 0.                         (3)

Moving endpoint terms vanish by flatness. Equivalently integrate the zero
extension from a fixed exterior radius; the identity is then immediate.
Apply (2) to the already differentiated sources in (3) at each fixed slow
point. This does NOT differentiate a pointwise upper bound, its weight or
the coordinate Jacobian. The nonlinear scale derivatives are already in
the assumed R-chart source jets. The change to x is used only to estimate
the radial integral at that fixed point.

For a>=1 radial derivatives, the fundamental theorem gives exactly

    partial_R^a partial_slow^b(I0 g)/(a! b!)
      =(1/a) partial_R^(a-1) partial_slow^b g/((a-1)! b!).    (4)

Sum (3) over pure slow derivatives and (4) over the remaining indices.
Each radial term maps injectively to a source derivative, and 1/a<=1.
Thus in the normalized THROUGH-order R,Z,T norm, for every finite n,

    ||I0 g||_n <= Kp(v) K zeta delta^-v,
    Kp(v)=1+(Rmax/2) C(v),                                  (5)

provided the source through-n norm is at most K zeta delta^-v.
No order-independent bound for an arbitrary smooth source is asserted:
K is its through-n source bound. The operator factor itself need not
increase with n in this normalized norm. Auxiliary parameters, when
present, can be included among the commuting derivatives under the same
hypotheses; Q37's Haar coefficient application uses R,Z,T.

For pressure apply this to g-rho integral g. Its exact integral is zero
because the retained rho has unit mass. Q37 separately bounds that full
projected source using its derived rho and whole-hull moment constants.
Kp is not a bound on an unprojected nonzero-mass primitive: such a source
has a nonzero exterior tail, incompatible with a vanishing flat weight.

## 3. Weighted targets and fixed projectors

For e=1,2, the oriented target primitive is
T_e f=R^-e I0(R^e f), with integral R^e f=0. Set

    P_e,n=sum_{j=0}^{min(e,n)} binom(e,j) Rmax^(e-j),
    N_e,n=sum_{j=0}^n binom(e+j-1,j) Rmin^(-e-j).

Normalized derivative norms are submultiplicative. Applying (5) between
the two radial multipliers yields T_e norm <=N_e,n P_e,n Kp(v).
Their sum safely supplies Q37's two-channel primitive input. R powers
have no slow derivatives in this chart; no missing scale differentiation
is hidden in these factors.

Q36's unweighted fixed-profile sigma bounds and actual moving-support
floor exp(E) supply

    Bump_n=exp(E)(Sigma_1,n+Sigma_2,n).                       (6)

This is a physical R-chart projector bound, so no logarithmic pullback
factor belongs in (6). The floor is applied once. The normalized mass,
support and fixed-profile jet hypotheses remain Q36's, with no new bump.
Since delta<=1, a pole-zero bump bound also suffices at nonnegative v.

## 4. Attach the actual requested pole and retain pressure memory

`NewPoleMaps` accepts the actual Q37 state and derives all three remaining
map entries from its same target geometry and profile. Q37 calls it at
order n+1 and

    v=2 max(6(n+2)+9,6(n+2)+d+6), d=Q36 target pole.         (7)

The kernel is recomputed at v, not reused at d. The output keeps Q37's
old scalar additions and both whole-pressure moments, including the moving
c_rho product. All six covariance entries and signed squares remain.
The factory needs no independent numerical map constants. Its analytical
validity still requires (3), flat zero extensions and the actual weighted
source/input contracts. Merely instantiating it is not their verification.

Tests compare (1) to the historical formula at several poles, independently
check normalized radial multiplier derivatives and (4), detect the
nonzero-mass obstruction, and run Q36/Q37 through the new maps with
manufactured lower bounds. The polynomial fixture tests identities only;
it is not a smooth flat weighted continuum source. Exact coefficient
comparison avoids decimal conversion of giant constants. Frozen Q37 and
its dependency closure are replayed.

Next construct fresh fast/slow corrections on the CURRENT post-signed
state, retaining existing means and pressure archives, and derive the first
complete unaveraged source estimate. Background/canonical/input support
premises, uniform actual pulses, full viscosity/cutoff interactions,
complete residual improvement and blow-up remain unresolved.

Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.
