# Q132: closed Gevrey-3 jets of the SAME leading exact five-row repair

Parent: `604d4566f6b9153bd1dcbab694a82d5a18d2bb28` (Q131).
Fresh fetch on 3 November 2026 agreed with this head. Q131 research CI
41709384070 and field CI 41709384251 were in progress; Q130 field CI
41699897210 had passed. Current notes, Q131 frozen evidence, the original
quadratic repair, Q63, Q122 and the repair-state neighborhood were read
first. Historical evidence, bump supports, pressure gauge, accepted roots,
phase origin, N, scales and C_setup^-20 width remain unchanged. The setup
activation conclusion remains local leading core/first collar.

Q131 bounds the modulated fields and all five incoming moment differences.
Here we close the subsequent NORMALIZED leading quadratic repair's angular
Gevrey recurrence. The result concerns the SAME already accepted small C1
root. Actual target binding, positive-field and low-order cone gates remain
premises. There is no new root, frequency choice or contraction in every
high derivative norm. This is distinct from Q122's positive-order LINEAR
repair and Q127's earlier eta-dependent quadratic coefficients.

## 1. Full physical deficit and the actual leading normalization

On the original reserved patch X=rho*exp(y), 0<y<5,

    U0=0, E0=e(eta) exp((-1/2-lambda)y),
    e(eta)=e_star/(1+eta^2), e_star>0, 0<lambda<=.01.

The five field corrections are e times the original coefficient/bump sum.
Their physical moment changes have EXACT scale factors in ordinary order
(M,I,J,S,Pi)

    rho*e, sqrt(2)*rho^(3/2)*e,
    sqrt(2)*rho^(3/2)*e^2, rho*e^2, e^2.        (1)

Let d be the full desired physical correction. It includes both the
PREMOD target-minus-profile defect and minus the Q131 modulated-minus-
PREMOD cumulative difference. Require bounds for ALL FIVE baseline defects;
zero baseline entries require an actual preceding target-match identity.
Do not silently infer them from the fact that the modulation has ended.
The modulation's incoming traces, including Pi, remain in Q131's memory.

If Q131 gives prefactors Dmem_i and base Lmem, and the five baseline
functions have raw bounds T_i L_i^m(m!)^3, set

    D_i=Dmem_i+T_i, L=max(Lmem, all L_i)>=1.

Then |d_i^(m)|<=D_i L^m(m!)^3. Signs are retained in the exact map; only
the envelope uses absolute values. Normalize (1) to ordinary m,i,j,s,cp.
The actual constant B,Q map is

    Bc+Q(c,c)=z,
    z=(m,(m-j)/lambda,i,-s,cp).                  (2)

Here j=J/(sqrt(2)*rho^(3/2)*e^2), NOT the e^-1 factor appropriate in a
positive-order linear system. The signs and the E/U square and pressure
square terms of Q must remain. Disjoint U/E supports kill only the J
cross term. B,Q are independent of eta only after this SAME normalization.
Equation (2) is stated for exact unit-mass integrals, not quadrature matrices.

The normalized eta Taylor coefficient bounds for 1/e on [-1,1] are
(2,2,1)/e_star and zero above degree two. For 1/e^2 they are
(4,8,8,4,1)/e_star^2 and zero above degree four. Multiplication by a
polynomial with coefficient bounds p_j costs at most sum p_j in a
Gevrey-3 envelope with base L>=1: (m-j)!^2<=m!^2 and L^(m-j)<=L^m.
Thus the exact polynomial costs are 5/e_star and 25/e_star^2.

Set g=max(1,rho^-2), which bounds rho^-1 and
(sqrt(2)*rho^(3/2))^-1. Conservative ordinary prefactors are

    A_M=g*5D_M/e_star, A_I=g*5D_I/e_star,
    A_J=g*25D_J/e_star^2, A_S=g*25D_S/e_star^2,
    A_Pi=g*25D_Pi/e_star^2.

Consequently

    Z=max(A_M,(A_M+A_J)/lambda,A_I,A_S,A_Pi),
    ||z^(m)||infinity<= Z L^m(m!)^3.             (3)

The inverse-lambda cost, both amplitude powers and every full pressure
row survive. A value at a sampled angle is not one of these function norms.
This is a conditional binding from Q131 PLUS the actual baseline target
identity/data, not a proof that those data vanish in the physical profile.

## 2. Same accepted root: no higher-norm smallness premise

Use the retained exact beta=||B^-1||<=1000 and q=||Q||<=10^5,
or valid upper bounds for the same map. Let cstar bound the ALREADY
accepted C1 root and retain the original pointwise gate

    2 beta q cstar <1/2.                        (4)

It gives ||(B+2Q(c,.))^-1||<=gamma=2beta. Real IFT and uniqueness in
the accepted ball identify the same smooth root, including the same
pressure-square correction. The separate original positive-field/state
neighborhood and finite-frequency acceptance criteria are not replaced
by (4). Neither entering cstar nor checking this rational inequality
proves the existence/acceptance of the physical root.

For m>=1, write d_m=||c^(m)||/m! and use constant B,Q to obtain

    d_m<=gamma[ Z L^m(m!)^2
                          +q sum_(j=1)^(m-1) d_j d_(m-j)]. (5)

The highest terms are in the Jacobian, not the lower sum. Dividing by
(m!)^2 and using (j! (m-j)!)^2<=(m!)^2 shows that the transformed real
jets are dominated by the nonnegative analytic majorant

    V=a Lx/(1-Lx)+b V^2, V(0)=0,
    a=gamma Z, b=gamma q.                       (6)

For a>0 let T=min(a,1/(8b)) (T=a if b=0), and

    R=[8L max(1,8ab)]^-1=T/(8aL).

On |x|<=R, |V|<=T, the forcing is <=T/6, the quadratic image <=T/8
and its derivative <=1/4. Thus the map takes the ball strictly inside
T/2 and is a contraction. Its positive formal coefficients agree with
Picard iteration and dominate the actual FINITE real inverse recurrence.
Cauchy gives V_m<=T R^-m<=a R^-m. Therefore

    ||c^(m)|| <= (gamma Z)
            [8L max(1,8 gamma^2 q Z)]^m (m!)^3, m>=1. (7)

For Z=0 all positive jets vanish by (5); no division by Z is needed.
At m=0 keep cstar, not gamma Z or T. Only the majorant is analytic in x;
no analytic eta extension through the original flat cutoff is asserted.
The bound permits arbitrarily large higher jets without a new acceptance
test in each Ck. Its base depends on Z, so it is not an unconditional
uniform-in-N smallness claim. If Z<=Zconst/N with fixed L and N>=1,
one may keep gamma Zconst/N in front and replace Z by Zconst in the base.
Q131's actual input/frequency dependence must be checked before using that
corollary. Low-order acceptance must continue to use its sharper original
C1 estimates rather than (7)'s coarse envelope.

## 3. Repair fields, supports and conditional pressure exit identity

On each eta disk of radius 1/4 about [-1,1], Re(1+eta^2)>=15/16,
so |e|<=2e_star and its normalized derivatives are <=2e_star*4^m.
Let K be the base in (7), C=max(cstar,gamma Z). A finite product gives

    |partial_eta^m(e*c_l)|
           <=2e_star*C*(4+K)^m(m!)^3.           (8)

No analyticity of c is used. Q63's original bumps have fixed disjoint
supports, so at any y only one term contributes to each field. Using its
actual bound B_r=162*32600^(r+1)*((r+1)!)^3 gives

    |partial_y^r partial_eta^m Delta(E,U)|
       <=2e_star*C*(4+K)^m(m!)^3 * B_r.        (9)

Ordinary X derivatives additionally pay x_min^-r*r! on these same
supports by the signed Stirling identity. Inverse widths and the full
radial factorial remain. No physical support or original N is altered.
The actual positive-field gate is still the old small-coefficient gate;
(8),(9) alone do not prove repaired positivity or a cone comparison.

Conditional on binding d to the true target, equation (2) restores ALL
five terminal moments as functions of eta. Fields already coincide after
the last fixed bump. Their later densities are therefore identical, and
the terminal moment equality plus the same pressure gauge propagates full
pressure equality to every later radius. All eta derivatives match by
the same functional identities. Restoring four rows would leave c(eta)
in Pi and 4Aphys eta c-(1-eta^2)c_eta in the forward state. No pressure
reset, missing pressure-square term or sample-angle match supplies this
exit identity. This does not establish the complete wave-covariance
mean correction with axial velocity and pressure in physical pulses.

## 4. Checks, binding and remaining obligation

Exact polynomial controls check both reciprocal powers, five physical
normalizations, row ordering, signs and lambda loss. A scalar square-root
model and independent five-component quadratic jets check (5)--(7),
including non-small high derivatives of a small C1 root and countercontrols
for dropping either energy or pressure squares. Rational controls verify
both contraction gates, zero target, linear limit and binary budgets at
order 100005 without factorial expansion. Field product controls check
(8); no numerical bump quadrature replaces exact support/integral premises.
The audit replays Q131 and checks the frozen source/evidence hash closure;
fresh output must match byte for byte. Manufactured diagnostics and passing
tests are not a continuum proof.

Next bind the actual baseline target, PREMOD/loop/value caps and the SAME
low-order accepted-root and positive-field/state gates. Then assemble the
final native profile norms and join them to the inner and outer cones.
Global coefficient/reserved-mean supports, uniform actual physical pulses,
complete wave-covariance mean corrections and the viscous/cutoff/all-
interaction residual remain open. No unconditional pulse bound, complete
residual improvement or blow-up is established. Historical evidence is
preserved in place; this new checkpoint records the leading/positive-order
normalization distinction rather than rewriting the older records.
