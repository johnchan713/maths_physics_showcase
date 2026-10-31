# Q60: a finite path stencil for the canonical primary pulse

[Q59](../core_current_head_q59_squared_root_budget/README.md) separates the
slow order d from the total order J in the axial covariance test. Here
the positive path orders of the uncut homogeneous vector require **no
mixed generator-jet input**. The price is a recurrence with four
numerator coefficients beyond the constant and two denominator
coefficients beyond the constant. Its slow input is the full background
F,G through d+1. This is a conditional estimate on the original bands;
the order count does not enclose the actual background seminorms.

## Clear the denominator only on the physical solution

In the fixed lifted coefficient chart let x be the slow variables, v the
path time and prime differentiation in v. Discrete data, including L,
are frozen. The [phase formulas](../core_phase_coefficient_budget/README.md)
give n=(x0,pθ/R,pz)+v m with m=(-H_R,0,-ε H_Z), H=pθ F+pz G.
K has rows (0,-2F,0), (2F+R F_R,0,0), (G_R,0,0); m,K are independent
of v. Put N=nᵀn and γ=εk². For the original |harmonic|=1 vector t,

    N t' = B t,
    B = -N K + n(nᵀK-mᵀ) - γ N² I.                         (1)

Indeed the full ambient generator of the [entrance
lemma](../core_canonical_entrance_jet_closure/README.md) adds
(γN+p)n nᵀ/N to B/N, where p=P'/P. Its equation propagates
(nᵀt)'=p(nᵀt); the actual canonical entrance is tangent. Thus (1)
holds identically in x,v for this family. Differentiate this identity,
including derivatives of its coefficients: differentiated vectors are
generally **not** tangent. The term -n mᵀ is the moving-normal pressure;
dropping it gives (nᵀt)'=mᵀt at tangency instead of zero. No
pressure response in the background or forced inverse is removed.

Keep the ambient propagator, actual entrance and rate in Q59's
pure-slow Duhamel line. We have not proved a Gaussian bound for the
unextended ambient propagator or set derivatives of p to zero. In this
fixed chart p is a reference function of frozen constants and v.
In particular the positive-path recurrence does not need high reference
rate jets; this is a property of the homogeneous physical solution,
not a change of rate for general forcing.

## Positive recurrence on a closed path

At each real base point use w for the increment of v. Write brackets
for coefficients of w. Exactly,

    n(w)=n+w m,    N(w)=N+2w nᵀm+w² mᵀm,    B(w) as in (1).

All slow derivatives preserve degrees 2 and 4. For each i through d let
D_ij, E_ij bound respectively |∂x^i [w^j]N(w)|/i! and
||∂x^i [w^j]B(w)||/i!, uniformly in the base point and in choices
of i slow coordinate slots. Use Euclidean/operator norms. A concrete
positive recipe is to use phase (4) at fixed v for n,m,K, divided by
i!, insert n and m as the two coefficients of n(w), and perform the
displayed products by ordinary convolution in (i,j), replacing signs
by plus and contractions by products of norms. Grouping i labeled
derivatives gives binomial multiplicities canceled by i!. Use the
original whole-path bound Q≥1/N; no individual component is divided by.

Let T_ab bound |∂x^a ∂v^b t|/(a! b! P(v)). Initialize T_a0 with
Q59's pure-slow U_(a,0)/a!, using its original ambient value and
entrance bounds. For b≥1 use

    T_ab = (Q/b) [
      Σ_(i≤a, j≤min(4,b−1)) E_ij T_(a−i,b−1−j)
      + Σ_(i≤a, j≤min(2,b−1), (i,j)≠(0,0))
                  (b−j) D_ij T_(a−i,b−j) ].               (2)

Indices are nonnegative. This is coefficient comparison in (1) at
path degree b−1: the isolated term is b N times the unknown. P(v)
in the norm is a common value at the base point; (1) differentiates
t, **not** t/P. The second sum includes j=0,i>0. Fill increasing b,
then increasing a, on a≤d and a+b≤J. The same-b terms use strictly
smaller a; all other terms have smaller b. Smooth finite jets suffice,
including one-sided jets at the path ends. There is no analyticity
premise on the slow data or smooth cutoffs. The implementation evaluates
finite orders; at the stated large orders (2) defines constants without
expanding them. Each fixed output is polynomial in the input bounds.

## What enters the axial test

The affine phase uses first slow derivatives of the full F,G; hence
(1)–(2) and the pure-slow generator through d use at most d+1 such
derivatives (the low-order propagator and frame gates are retained).
These are the finalized fields, with modulation, exact repair and all
pressure memory. The [finite-head estimate](../core_background_finite_head_budget/README.md)
is a recipe for bounding them from the corresponding profile, chart
and coefficient data; its premises are not implied by the path stencil.
For Q59's largest d=n+99997 the generator side thus requests at most
n+99998 slow derivatives, even for total order n+199991.

For the explicit coefficient lift of [Q21](../core_current_head_q21_lifted_path_join/README.md),
Y_i=Jcov^Δ Haux and (ξ,c_i v-r0) are the linear transverse/longitudinal
coordinates of Y_i minus its frozen center and copy. On the native
band let lξ,lv be the maximum absolute first coordinate derivatives of
ξ,v in its auxiliary chart (lv includes 1/c_i); there are no higher
ones. For a slow derivative of order a and b auxiliary slots, the
localized vector χ(ξ)ψ(v/L)t has raw norm at most

    P(v) Σ_(i+j+k=b) b!/(i! j! k!)
        lξ^i lv^(j+k) ||χ^(i)||∞ L^(-j) ||ψ^(j)||∞ a! k! T_ak.

This bounds W in Q59 using P≤1, on the downward-closed index set.
χ includes the physical radius of the transverse cutoff. Use the usual
slow R/U composition if changing radial variable and apply the common
atlas conversion once; alternatively use the common auxiliary slopes
directly, without counting that covering twice. A moving physical graph
requires its own chain rule. The Haar translation cancellation alone
would not justify ignoring a varying center in a pointwise estimate.
The unchanged support and flat extension apply to the localized vector.

This discharges the high path-generator and frozen auxiliary-map slots
of the principal moment *relative to* the stated slow inputs. The next
quantitative step is to enclose the unstripped slow target, partition,
column/entrance and full-background inputs together on the original
integration hull, and check the resulting axial bound in the tail norm.
The distinct radial pressure remainder, uniform physical-scale pulse
bounds, complete mean corrections and viscous/cutoff/interaction residual
are not certified. No blow-up conclusion follows.

Run `python test_audit.py -v` and `python audit.py --verify-record
evidence.json --output /path/to/fresh.json` in this directory.
