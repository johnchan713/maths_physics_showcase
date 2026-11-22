# Q220 shifted two-edge weighted-kernel source audit

This checkpoint audits the compact modulation/exact-repair kernel used by
Q216.  It keeps the physical transverse shift and the complete shifted
moment.  The result is a typed operator theorem: for every **fixed** output
jet order and every **fixed finite** endpoint pole, the correctly oriented
two-edge operator preserves the full flat weight and the same pole.  It does
not prove that every generated source belongs to that input class.

## Directed shifted primitives

For the smooth zero extension of an input `F`, write

```text
I_-(U,y) = int_{a<0} F(U+a, y+M a v_r) da,
I_+(U,y) = int_{a>0} F(U+a, y+M a v_r) da.
```

The compact inverse has the form

```text
T_e F = g_e ((1-chi) I_- - chi I_+),
A_e F = k_e chi_U int_R F(U+a, y+M a v_r) da.
```

The signs and half-lines are essential.  Near the left edge, `1-chi`
selects the backward primitive, whose integration ray stays away from the
right edge; near the right edge, `chi` selects the forward primitive, whose
ray stays away from the left edge.  Reversing either orientation leaves a
nonzero whole primitive at the opposite flat edge and cannot preserve the
flat weight.

Differentiation is taken at fixed integration variable `a`.  The translated
point `y+M a v_r` does not depend on the differentiated coefficient
variables, so this step creates no additional power of `M`.  Smooth zero
extension removes endpoint terms.  The `A_e` term retains the whole shifted
moment; zero Haar mean by itself does not delete it.  Its cutoff derivative
is supported in an interior transition set where the weight has a positive
lower bound.  Output jets through order `m` require cutoff jets through
`m+1`.

## Every fixed finite pole

Use

```text
zeta(x) = exp(-w^2/x^2 - 4/(L-x)^2),
delta(x) = min(1,x,L-x),
r = min(1,L/4,g),
```

where `w,L,g` are fixed positive data and `g` is the plateau gap.  If each
requested input derivative is bounded by

```text
epsilon^b zeta delta^-p,
```

then the directed primitives and the transition term are bounded by the
same weight and the same `delta^-p`.  The source proof splits the constant
into three exact branches:

```text
inner:    D(w^2) exp(4/(L-r)^2),
outer:    D(4)   exp(w^2/(L-r)^2),
interior: L E_p  exp((w^2+4)/r^2).
```

Here `D` is the Gaussian-tail bound and
`E_p = max(1, ceil(p/2)! (2/w^2)^ceil(p/2))`.  Thus every fixed finite `p`
has a finite constant: the Gaussian absorbs every fixed polynomial pole on
the two collars, and the interior branch uses the positive weight floor.
The map returns `p`, rather than recursively requesting a larger pole.  The
constant is not uniform in `p`; for the Q215 request
`p=2,626,513,221,816` it must remain an exact symbolic finite expression,
not a floating-point number.

This is an analytical estimate inherited from the cited kernel and weighted
map derivations.  The executable record checks its typing, orientations,
quantifiers and frozen dependency closure; those checks are not a continuum
proof.

## Corrected source-class boundary

The theorem's input type cannot be weakened silently.  A half-weight source
does not enter the full-weight class for any finite pole, since at the left
edge

```text
zeta^-1 x^p sqrt(zeta) -> infinity  as x downarrow 0.
```

Consequently, the blanket claim that every compact source has full `zeta`
weight remains false.  Harmonic rows in the inherited construction carry
only `sqrt(zeta)`.  The physical covariance targets audited in Q200 carry a
full `zeta` factor directly from their product structure and exact state
differences; that is a source-specific route.  Unweighted jets, compact
support and moment cancellation alone still do not imply the weighted
target estimate.

## Established scope and remaining obligation

Q216's two-edge kernel quantifier is supported for typed full-weight inputs:
all fixed jet orders and all fixed finite poles, with the physical shift,
correct edge orientation, complete pressure memory and unchanged pole.  The
actual membership of every generated source must be checked separately by
source family.  This checkpoint therefore establishes no blanket
full-weight pulse estimate, hierarchy-uniform constant, complete residual
improvement, force exclusion or blow-up result.

The next obligation is to audit the remaining compact/slow source theorems:
the zero-Haar Fourier inverse and fixed slow response, retaining axial
velocity, the moving pressure product, and the finite `lambda^-1` loss.
