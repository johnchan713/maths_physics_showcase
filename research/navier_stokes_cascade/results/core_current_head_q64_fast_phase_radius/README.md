# Q64: an explicit analytic radius in the moving loop's fast phase

The moving-loop input left by Q63 need not include independent high
fast-phase seminorms. The original real slow profile can be merely smooth,
yet its loop is analytic in the phase variable on a uniform strip.
Here is an explicit, conservative radius from the unchanged loop cap
and phase-density floor. It is independent of the physical band and
does not change the selected modulation frequency N.

## A complex phase disk with no denominator zeros

Use exactly the [loop bounds](../loop_modulation_bounds/README.md):
G≥2 bounds |μ|,|p|,|t_s|,a,v,|E| and |μp| on the frozen compact
rectangle; v≥2, a>0 and real w≥m=G⁻³. Put z=μp, and retain the
removable representation

    t=t_s+d₀ μ ∫₀¹ q_z(uz,θ)du,
    q(z,θ)=exp(z sin θ)/I₀(z),       0<d₀≤1.

The parameters are REAL here, even when θ is complex. For |Im θ|≤1,
|sin θ|,|cos θ|<2, I₀(uz)≥1 and |(log I₀)'(uz)|≤1.
Thus |q_z|≤3exp(2G) and |q_(zθ)|≤(2+6G)exp(2G). Define

    T=4G exp(2G), D=8G² exp(2G), M=G T D,
    ρ=min(1/2,m/(4M)), δ=mρ/4.                           (1)

Then |t|≤T, |tθ|≤D, and |wθ|≤M for
w=a(1+t²)/(2πv). At every real θ₀, throughout |θ−θ₀|≤ρ,
|w−w(θ₀)|≤m/4 and |w|≥3m/4. In particular neither w nor
1+t² vanishes on this disk.

Let φ(x,θ)=∫₀^θ w(x,r)dr with the actual fixed origin.
For |ζ−φ(x,θ₀)|<δ solve θ=θ₀+u by the contraction

    u=(ζ−φ(x,θ₀))/w(x,θ₀)
       −[φ(x,θ₀+u)−φ(x,θ₀)−w(x,θ₀)u]/w(x,θ₀).          (2)

It maps |u|≤ρ into |u|≤ρ/2 and has Lipschitz constant ≤1/4.
Its analytic solution agrees on the real axis with the original
moving inverse Θ. Local branches glue by the identity theorem;
the real orientation-preserving circle map and uniqueness identify
overlaps. The resulting pulled-back periodic loop is holomorphic on
|Im ζ|<δ. The lift Θ itself is not a bounded periodic function;
only periodic compositions and its positive slow derivatives are used.

Since 1/(1+t²)=a/(2πv w), its modulus is at most G/m on
these inverse disks. For a_L=v/(1+t²), b_L=−vt/(1+t²), both
periodic sources of the [original primitives](../stress_realization_audit/README.md),

    Acal_ζ=−(a_L−a)/2,   Bcal_ζ=E(b_L−b_s)/2,

are bounded by F=2G³(T+1)/m. Here |b_s|=|a t_s|≤G².
The real zero-mean primitive has sup norm ≤2F; integrate vertically
over distance ≤δ<1 to extend it with bound 3F. The zero real mean
of each source gives zero integral around every horizontal period
by holomorphy, so the primitive is single-valued and periodic.
No additive constant or source mean is changed.

For either primitive f, Cauchy's inequality gives

    sup_real |∂ζ^j f| ≤ 3F j! (2/δ)^j.                   (3)

One can also shift a Fourier contour by ±δ/2:
|fhat(k)|≤3F exp(−πδ|k|). This is decay in the auxiliary loop
phase, not a claim of exponential decay in the physical torus
variables of an entire localized wave or in slow radius.

## Mixed slow derivatives use finite real jets, not slow analyticity

The [signed-root and moving-phase theorem](../modulation_repair_all_order/README.md)
supplies μ(x) as a smooth function even at zero variance. Keep its
real finite slow jets, those of t_s,a,v,E, and the original positive
slow denominators; do not divide by p or μ. For any slow multiindex
α=(r,k), finite composition of the removable integral bounds
∂x^α t and ∂x^α w on |Im θ|≤1. I₀ derivatives occur only
at real arguments; their fixed-order bounds can be obtained by
differentiating its real integral, then its positive reciprocal.
This requires parameter jets only at subindices of α.

Here is a precise finite procedure for mixed inverse bounds. At a
complex base point Θ₀ over |ζ−ζ₀|<δ, expand only in the REAL
slow increment h. Once lower-degree coefficients of Θ are bounded,
extract the coefficient of h^α in
φ(x+h,Θ_<α(h))−ζ, omitting the unknown degree-|α| inverse
coefficients. Divide its absolute majorant by m/2. The actual
denominator w(x,Θ₀) has modulus ≥3m/4. For the required θ
Taylor derivatives use Cauchy radius 1/2 around Θ₀ (which has
|Im Θ₀|≤ρ/2≤1/4); a uniform slow-derivative bound on the
larger strip gives j! 2^j. For pure slow derivatives of φ use
a real integral of length at most 2π and a complex segment of
length ≤1. Its θ derivatives recover w. This fixes all constants
using finite operations and the stated slow inputs. Higher-degree
unknowns occur only linearly, as in the original IFT recurrence.

Compose the finite jets into the two sources and integrate their
zero-mean primitives. If H_α bounds the resulting normalized
slow jet sup |∂x^α f|/α! on the phase strip, then

    |∂y^r ∂η^k ∂ζ^j f| ≤ r! k! H_(r,k) j! (2/δ)^j.     (4)

These H require slow parameter jets only at subindices of (r,k).
No high j-th slow derivative is requested by (4). Differentiation
is at fixed ζ, including φ_x and the motion of Θ; differentiating
at fixed θ instead would be a different function. The construction
is local in slow coordinates and patches by the original real
solution. It assumes only finite real differentiability of the
slow fields, not a complex extension through a flat cutoff.
One-sided slow endpoint jets follow by their smooth limits.

## Substitute before the fixed-frequency pullback

For f=Acal or Bcal the exact chain rule at ζ=Ny yields

    |∂y^r ∂η^k [f(y,η,Ny)/N]|
      ≤ (r! k!/N) Σ_(j=0)^r (2N/δ)^j H_(r−j,k).        (5)

The factor binom(r,j) cancels the factorials from (4).
Thus the high radial modulation bound of Q63 can be formed from
slow input jets and the explicit phase radius, without tabulating
independent mixed fast-phase derivatives. Apply the original
exponential finite product rule for E exp(Acal/N) and retain its
actual value bound; U+Bcal/N is linear. Ordinary X derivatives
retain the full logarithmic Stirling conversion. This reproduces
the legitimate N^(r−1) scale for the difference, with fixed-order
constants, rather than a false O(1/N) radial error.

The next remaining inputs are the pre-modulation slow profile jets,
the signed variance-root right-hand side and its derivatives, and
the low-order normalized five-row defect used by Q63. Complete
pressure memory remains in p=p₂ and in that defect. This checkpoint
does not bound them from the earlier regions or verify the physical
pulse/inverse gates. The axial estimate remains conditional; the
radial pressure tail, full residual and blow-up are not proved.

Run test_audit.py and audit.py --verify-record evidence.json
--output /path/to/fresh.json. The controls test rational inverse
disks, the moving inverse chain rule and frequency factorials;
they do not numerically evaluate the enormous actual G.
