# Q67: pre-modulation state jets with incoming pressure memory

Parent: `a009049a36d8a10f159ab370cd7cbc6ff1c6812b` (Q66).
This supplies a finite forward enclosure for Q66's a,b_s,p2 inputs on a
compact pre-modulation interval. It identifies the required input orders
and an obstruction to dropping the incoming trace. It does not evaluate
all regional field seminorms. Historical evidence, supports, pressure gauge,
frequency and accepted repair root are unchanged.

## Forward memory, before modulation or repair

Fix a genuine pre-modulation interval [x0,x1], x0>0, with E>=e_*>0,
|eta|<=1 and 0<h<1/2. Both endpoints are fixed in eta. Use the same
profile and its actual incoming traces at x0; do not set these to zero.
Let r=sqrt(2X), Aphys=1/2+h, D=1/2-h, d=1-eta², L=1-2h eta².
The five primitives (M,I,J,S,Pi) have respective X densities

    U, rE, rEU, U²-E²/2, E²/(2X).                       (1)

Pi includes its incoming pressure datum. For each density let g_ij bound
|partial_X^i partial_eta^j g|/(i!j!) uniformly on the whole interval,
and let t_j bound the corresponding normalized trace at x0. Then

    P_0j <= t_j+(x1-x0)g_0j,
    P_ij <= g_(i-1,j)/i, i>=1.                         (2)

The first formula is differentiation under the finite integral, the second
is the fundamental theorem of calculus with factorial normalization.
The bounds require actual uniform density bounds, not values at a sampled
point. Intervening supports are included in the integral. One-sided jets
at eta endpoints follow from the same smooth extensions as the input.

Define w=X-2D eta M-d M_eta. The full original forward map gives

    a=1-2X E_X/E, b_s=2X U_X/E,
    p2=(-w U+D(M-eta M_eta)+4h eta S-d S_eta
                         +X(4Aphys eta Pi-d Pi_eta))/(L E). (3)

These identities follow directly from the
[complete forward map](../core_activation_volterra_recheck/README.md),
with XW=w. No terminal pressure reconstruction or repair cancellation is
used here. I,J are retained in (1) for the full five-row map, although
(3)'s three inputs use only M,S,Pi. This distinction does not authorize
dropping either of the other rows in the eventual repair.

## An explicit finite table algorithm and derivative count

Multiply nonnegative normalized jet tables by ordinary convolution.
Differentiation replaces f_ij by (i+1)f_(i+1,j) or (j+1)f_(i,j+1).
For E inverse use

    (1/E)_00 <= 1/e_*,
    (1/E)_alpha <= (1/e_*) sum_(0<beta<=alpha)
                                E_beta (1/E)_(alpha-beta). (4)

This follows by differentiating E(1/E)=1 and removing the zero-beta term.
For X and eta use their exact affine increment polynomials with value
bounds x1 and 1. For 1/X use coefficient x0^(-i-1), with no eta entries.
For r use sqrt(2x1)|binom(1/2,i)|x0^(-i). For L inverse use its positive
floor 1-2h and the finite geometric series with incremental majorant
(4h z+2h z²)/(1-2h), z denoting the eta increment. Take absolute signs
in (1),(3), use (2),(4), and truncate after every operation to the desired
downward-closed set. This is an explicit positive finite-sum bound for
all three state inputs; their positive-degree entries supply Q66 and Q65.
It assumes no slow analyticity, and does not extrapolate a low-order norm.

For output i+j<=m, j<=l (0<=l<=m), the local a,b_s terms request
E,U through (i+1,j). A differentiated primitive in (3) requests trace
order j+1 when i=0; its integral then needs density order (0,j+1).
When i>=1 it instead requests density order (i-1,j+1), with factor
(j+1)/i. Thus sufficient common inputs are

    E,U: total order <=m+1, eta order <=l+1;
    incoming M,S,Pi traces: eta order <=l+1.             (5)

The table algorithm keeps the finer stencil just described. In particular
high radial orders never demand high radial derivatives of an incoming
trace. Products, inverses and multiplication by X,eta only distribute
subindices. This proves (5) for the finite state map, not for the complete
pulse ODE or the full finalized profile. Q63's final-profile index set
cannot simply be relabeled as a pre-modulation state bound without this
additional bookkeeping.

## Why the incoming angular derivative is indispensable

Keep E,U and all their jets fixed and replace only incoming Pi by
Pi+c(eta). The forward map changes by the exact amount

    Delta p2=X(4Aphys eta c-(1-eta²)c')/(L E).             (6)

For c(eta)=eta at eta=0 this is -X/E, despite identical local velocity.
More sharply, c(eta)=sin(k eta) has uniformly bounded value but at eta=0
produces |Delta p2|=kX/E. Consequently no state bound depending only on
local velocity jets and the VALUE of incoming pressure can hold uniformly.
This is an obstruction for the forward map, not a claim that arbitrary c
is admissible in the full global construction. That construction fixes c;
its actual angular trace still has to be bounded, not suppressed.

Unlike Q62's finalized target reconstruction, (2) requires no future exact
repair or terminal heat identity. Using Q62's terminal trace for a different
pre-modulation profile would require proving equality of that trace first.
The next obligation is to evaluate (2)--(5) on the original annular and
exterior inputs and close the normalized repair defect on the same regions.
No uniform physical pulse, axial-tail rate, radial pressure-tail estimate,
complete residual improvement or blow-up follows from this conditional map.

Run test_audit.py and audit.py --verify-record evidence.json --output a fresh
path. Exact manufactured checks and provenance replay are not a continuum proof.
