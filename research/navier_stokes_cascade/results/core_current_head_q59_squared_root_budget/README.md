# Q59: bound the covariance product before differentiating the root

For the first-real axial tail, Q57 used the product of two independent
amplitude bounds. The two factors actually have the **same root**.
Here is a sufficient direct estimate on the canonical labels in Q57,
using Q58's leading weights. The high-order linear column solve and
the entrance data need only n+99997 slow derivatives, instead of the
n+199991 total-order amplitude archive. The large path/auxiliary orders
are still needed. This does not change the root used in the wave.

## Exact per-label coefficient and its slow solve

In the native band let hσ=χψ tσ be the localized canonical homogeneous
vector, with the slow partition ηβ outside it. For a single label the
principal tensor coefficient, before the common ε/2 normalization, is

    P_ab = ηβ² yσ hσ,a hσ,b,       y = Hbar⁻¹ Tβ / s₀,
    H=s₀ Hbar,   Tβ=s^(-1-h) T₀(X,η),   s₀=r₀²/√L.

This follows by squaring the positive physical root in the [actual
column solve](../core_localized_stress_assembly/README.md). H is the
whole normalized Haar column with the real angular factor 1/2 inside.
No partition, cutoff or new Haar factor multiplies its equation.
There is one ηβ² in P. Distinct-label products vanish by the existing
support separation. The linear formula applies even to diagonal entries
which are *not* matched by the two-column solve. It computes their root
factor only. The complete curl term retains the low-order Q57 bound.

Hbar, Tβ and ηβ are independent of the auxiliary coordinate, including
on the common covering. Use the unstripped smooth Tβ itself, not Tβ/ζ:
the original target is flat at both edges by the [target factorizations](../core_primary_target_cutoff_budget/README.md).
No inverse edge weight or reciprocal square-root derivatives are needed
in this unweighted moment calculation. This does not assert smoothness
of arbitrary divided targets. On the active shell the original inverse
value bound K₀ from the [root budget](../core_primary_root_budget/README.md)
is still needed. Smooth flat extension of the matched product and the
original positive cone domain are still premises of the physical wave.

Let H_j and T_j bound exact-degree factorial sums of the normalized
column and unstripped target in **slow** native variables. In a frozen
band s₀ is constant. A direct majorant for the exact-degree norm of
s₀ y is

    Y₀=K₀ T₀*,   Y_j=K₀ (T_j + Σ_(i=1)^j H_i Y_(j−i)),

where T₀* means the target norm at degree zero, not the profile.
This is Leibniz applied to Hbar(s₀ y)=Tβ, followed once at each
degree by its value inverse. In particular no power of the inverse
floor for an individual squared root is spent here. Raw through-order
d costs d! Σ_(j≤d) Y_j / s₀; 1/s₀=√L/r₀² costs a polynomial
in S by the fixed length bound. The target norms can be taken by
finite Taylor composition of the original smooth stress with the
Q58 scale/annular chart (and the factor s^(-1-h)), including all
incoming modulation, pressure and repair contributions.

## Count path and slow derivatives separately

For each Q57 shell choice let s be the integration count (unrelated
to the scale above), q the pure auxiliary count, and put
d=n+1+s, J=d+q. In a derivative of P of the prescribed Q57 type,
at most d slots are slow. None of the q pure auxiliary derivatives
can land on the root or slow partition. The pullback from the common
torus and slow band chart is block diagonal as in the [aggregate
atlas](../core_primary_aggregate_atlas/README.md); linear combinations
of native auxiliary derivatives are allowed. If M_d bounds raw slow
partition jets, and W_(d,J) bounds hσ jets with at most d slow slots
and total order ≤J, a sufficient native product bound is

    5^J M_d² (d! Σ_(j≤d)Y_j / s₀) W_(d,J)².

There are five factors including the two partitions. Each differentiated
term stays in the actual support, so summing on the common torus costs
no increasing number of labels. For example the original chart bounds
give a factor 64·188^J and the scale substitution S_i≤(25/16)S for
the product of transferred amplitudes (two value transfers ≤8 each).
Use this product norm **instead of** the amplitude-square Leibniz
bound in Q57, with a further factor two for a diagonal difference.
No pointwise estimate of a Fourier coefficient is made from its Haar
average.

For clarity the remaining W can be obtained by the following triangular
version of the [canonical entrance/ODE argument](../core_canonical_entrance_jet_closure/README.md).
Use fixed path coordinates x,v; x denotes the slow parameters. The
solution is independent of transverse ξ. Let A_(a,b) bound generator
derivatives with a parameter slots and b path slots (the full moving-normal
generator Ahat). Let I_a bound the corresponding actual entrance derivative
divided by P(0), K the ambient propagator constant, and L the path length.
For |∂x^a ∂v^b t| ≤ P(v) U_(a,b), sufficient raw bounds are

    U_(a,0)=K (I_a + L Σ_(i=1)^a binom(a,i) A_(i,0) U_(a−i,0)),
    U_(a,b)=Σ_(i≤a,j<b) binom(a,i)binom(b−1,j)
                      A_(i,j) U_(a−i,b−1−j)      (b≥1).

The first line is differentiated Duhamel with P(w) canceling the
propagator denominator; the second is differentiation of t′=Ahat t.
These sums index exact degree maxima; the parameter binomials sum the
multiindex Leibniz coefficients. They close on a≤d, a+b≤J. In
particular I is required only through d, and the high-total-order
generator requests have at most d parameter derivatives. Conversion
to the actual chart must include the complete path map and cutoffs:
compose the finite Taylor jet of t with x=x(slow), v=v(slow,y)
and multiply by χψ before taking W. Each parameter derivative needs
at least one slow slot, even if the center moves; derivatives of v
need not. This count does not bound those chart derivatives by one.
For the Haar-column H_j the fixed-v identity and Gaussian mass
in [the column lemma](../core_canonical_haar_column_jet_closure/README.md)
use only U_(a,0), a≤j. The pressure contribution to Ahat and every
cutoff derivative in the actual pulse remain.

At p=99995 the largest d is n+99997 in both channels; J is n+199991
(diagonal) or n+199990 (shear). These indices and K₀, target,
partition, path and generator bounds suffice on the original accepted
bands. They do not constitute numerical enclosures of those data.
Next trace the restricted generator and chart jets to the actual
phase/normal and background; the radial pressure remainder, uniform
physical pulse enclosure and full residual estimate are still open.
