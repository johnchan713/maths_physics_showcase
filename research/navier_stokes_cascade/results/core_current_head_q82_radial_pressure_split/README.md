# Q82: separate the leading covariance in the final radial pressure error

Parent: `d92bb5392a4ab5fe11762a1aac92edc2b81c28b9` (Q81).
`core_initial_a0_physical_tail_closure` already gives a conditional
3−5κ bound by requesting the entire final radial source through n+399995.
Here only two original principal products need the largest oscillatory
estimates. The other, nonlinear final-state source has a lower sufficient
jet order. This is an analytical reduction of inputs, not a numerical
enclosure of the selected fields or a new complete residual rate.

## Take the quotient on the correct state

Use the canonical first-real prefix AFTER the new fast and sharp slow
updates. The projected compact term is −A0 Q gr; its sign in the final
row is minus. Here Q removes Haar, I=c_i^-1 Nt^-1 on its range, and
B=1−κ. Use the final-state identity (5)--(7) of
`core_angular_pressure_normal_form` with a=2V/R, a1=2a/R−a_R. In its
notation, modulo one full characteristic derivative D_r K1,

    Q gr ≡ −a1 I Q J − a I Q Fsrc + Q(v²/R + Rrest),
    J=bv+V beta+beta v+W_rθ,
    Rrest=−tstar beta−(2b beta+beta²+Wrr)/R
          −Dz(b gamma+G beta+beta gamma+Wrz)+Wθθ/R+ε Dz² beta.

The omitted derivative here means equality AFTER A0 only: A0 D_r QK1=0
on the smooth compact zero extension. The pressure itself retains K1 and
the primitive of r1−rho Pmean. Fsrc includes the final angular residual,
slow time, axial flux, radial AND axial angular viscosity. It is not a
background F field. Its existing O_(2−3κ) bound comes from the angular
equation, without using the radial error being estimated. The scalar
compensator −rho Pmean in the final radial row is a separate term.

Write Pij for Q57's ORIGINAL dimensionless principal tensor (so that
its covariance contribution is ε Pij/2). On Q51/Q56's signed-target and
Q57's complete-curl hypotheses the difference W−εP/2 is O_(2−2κ):
the curl square has that order, and the complete signed change is O₂.
The subsequent mean updates do not change this wave tensor. Thus exactly,
with a remainder defined by subtraction on the right,

    A0 Q gr = (ε/2) A0 [ Q(Pθθ−Prr)/R − a1 I Q Prθ ] + A0 Q Rem.   (1)

The two terms in (1) have different temporal symbols. The radial weight
in the shear is the FULL a1=6V/R²−2V_R/R, not just 6V/R².
The original dU kernel multiplier f0=dR/dU applies to both terms.

Under the inherited all-fixed-jet bounds Rem is O_(2−3κ). Indeed
J−εPrθ/2 is O_(2−2κ) since beta=O_(2−κ), b=O₁ and v=O_B;
v²=O_(2−2κ). After removal of the indicated diagonal, the remaining
Rrest terms are O_(2−2κ) or better: tstar beta is O_(2−κ), Dz adds
one power to the order-one mixed flux, and ε Dz² beta is retained.
The factor I costs a fixed four extra jets and polynomial scale factors,
not epsilon powers. Finally I Q Fsrc has order 2−3κ. These statements
concern derivatives on the fixed native chart before moving to a physical
graph. All support, edge-weight conversions, slow multipliers and finite
primitive requests of the historical estimates remain necessary.

## A different integration count for each source

Take κ=1/100000 and target t=3−5κ. On a dyadic auxiliary Fourier shell
of size K, the joint divisor bound of Q57, s radial integrations and q
pure auxiliary derivatives give M^-s K^(s−r−q) times c_i^-r, where r=1
for I and r=0 otherwise. Haar is removed before inversion. For p an
integer >=3 take s=p+1 below sqrt(M) and s=p−1 above it; q=2p−r−s.
Relative to M^-p the two shell exponents are +2 and −2. The geometric
sums cost a constant; M=c_i^ρ ε^-κ and all fixed inverse c_i powers
are paid from the retained polynomial floor. Endpoint terms vanish by
smooth compact support. Output derivatives at fixed shift add n slots;
only n+s can fall on a slow kernel multiplier. Unlike Q57's axial tail,
there is NO exterior Dz factor or extra Z slot in (1).

| Input inside A0 | Starting order | p | Extra raw source jets (2p−r) | Extra slow multiplier jets (p+1) |
|---|---:|---:|---:|---:|
| Pθθ−Prr, with ε/2 | 1 | 199995 | 399990 | 199996 |
| Prθ, with ε/2 and I | 1 | 199995 | 399989 | 199996 |
| Q Rem, AFTER its internal temporal inverses | 2−3κ | 99998 | 199996 | 99999 |

Each row reaches starting order + pκ = t. These are sufficient counts for
this split, not sharp lower bounds on physical waves. In particular the
last row requests the reconstructed remainder through n+199996; computing
it requires additional lower-field jets, including the four slots for I.
It does not request the WHOLE final pressure source through n+399990.
This statement is not a bound on the maximum depth of the entire DAG.

## Full-background corrections cost order three on the comparison branch

On Q61's successful common-domain branch its normalized principal products
differ from their fixed-carrier leading-background comparisons by O₂ at
low jets. Here the prefactor is ε rather than the axial ε². Their direct
contribution to (1) is consequently O₃, still sufficient for t. Multiplying
by a1−a1⁰=O₂ gives the same order on the true order-one covariance.
This last bound uses the inherited full-minus-leading first derivative
and the positive-radius chart; keep its derivative and product costs at
higher fixed output orders. Q81's fixed-slow-fiber argument localizes both
input differences to X<=Xvel. The reference uses the SAME target, partition,
Haar maps and complete moving-normal equation; only the tensor coefficient
is compared. No pressure or physical wave is replaced. The large principal
rows in the table can therefore use the leading comparison and leading
multipliers on that branch. Their integration hull is still the FULL hull.

On Q61's other branch ε_*² Bcomp(S*)>1/2. The true leading-order terms
in (1) have a direct O₁ bound (four extra jets for I); they can be charged
at t since for 0<ε_*<=1,

    ε_*^(-(t−1)) <= ε_*^-2 < 2 Bcomp(S*).                  (2)

This uses t−1=2−5κ, not the smaller axial gap 1−5κ. Bcomp is the same
fixed-neighborhood polynomial constructed from VALUE comparison tests,
not differentiated or used as a physical cutoff. No reference inverse
is assumed on that branch. This argument bounds the principal terms;
it neither localizes the fallback nor removes Rem. It cannot pay an
arbitrary target t>3 by this same inequality.

## What remains

Conditional on the table's same-domain weighted polynomial envelopes,
(1)--(2) give the compact radial error at 3−5κ, retaining both pressure
traces, the full Haar compensator and the complete axial velocity. Q81
localizes errors of comparison only. A shifted pressure output can extend
outside their support. The axial Q57 high-order moment, the radial table's
leading envelopes and nonlinear remainder archive must now be joined to
actual common-domain bounds; neither finite checks nor the conditional
historical full-source theorem independently proves those inputs. No
unconditional pulse estimate, complete residual improvement or blow-up.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output FRESH.json`. Controls check exact radial integration
with a nonconstant base, the joint-divisor shell exponents and branch cost.
