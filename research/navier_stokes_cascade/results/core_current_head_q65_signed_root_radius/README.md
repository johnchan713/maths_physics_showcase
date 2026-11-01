# Q65: explicit signed variance-root derivatives from its scalar inputs

Q64 retains the slow jets of μ. They can instead be bounded by finite
composition of the two actual scalar inputs p=p₂ and
s=ζ(v_s)√((2+δ_L/2−v_s)/a), extended smoothly by zero outside the
original cutoff support. This checkpoint supplies a uniform local analytic
radius for the signed inverse μ=Ψ(p,s). It does not assert that p,s
are analytic as functions of the slow variables. No cutoff, root, phase,
frequency, pressure datum or support is changed.

## An explicit zero-free disk for the removable variance factor

Use the retained cap |μ|,|p|,|μp|≤G with G≥2, 0<d₀≤1, and choose
a symbolic upper bound E≥exp(G). Set

    A=4E, C=3EA, b=1/(2E⁴).
    q(z,θ)=exp(z sin θ)/I₀(z),
    R(z,θ)=∫₀¹ q_z(tz,θ)dt, V(z)=average_θ R(z,θ)².      (1)

The average is normalized real angular integration, with no complex
conjugation in V. At real z this agrees with the positive variance.
For real |z₀|≤G, I₀(z₀)≥1 and |I₀^(j)(z₀)|≤E. On the complex
disk |z−z₀|≤1/A,
|I₀(z)−I₀(z₀)|≤E(exp(1/A)−1)<1/2, so q has no denominator
zero. The reciprocal normalized Taylor coefficients are bounded by A^j:
induction uses E Σ_(k≥1) A^(-k)/k!<1. Product with the numerator
then gives |q^(j)(z₀,θ)|/j!≤3E A^j, uniformly for real θ.

Differentiating the removable integral, the factor (j+1) from q_z
cancels ∫₀¹t^jdt. Hence |R^(j)|/j!≤C A^j and
|V^(j)|/j!≤(j+1)C² A^j. These series converge in the smaller disks
used below; differentiation under the t and θ integrals is uniform.

The positive real floor also has an explicit value. With g=log I₀,
the tilted density is at least exp(−2|z|) times uniform measure, so
g''(z)≥exp(−2|z|)/2 by minimization over the subtracted constant.
For either sign of real z,

    g(2z)−2g(z)=z²∫₀¹∫₀¹g''((u+v)z)du dv
                  ≥z² exp(−4G)/2.

Since V=(exp(g(2z)−2g(z))−1)/z², extended at zero,
V(z)≥b for |z|≤G, including V(0)=1/2. This argument never divides
a numerical residual by a small z.

Let Rv=b/(16AC²). For |u|≤Rv the nonconstant Taylor part of V
is at most C²[(1−A|u|)^(-2)−1]≤3C² A|u|<b/2.
Thus the positive square root of V at each real center extends
holomorphically to that disk, with modulus ≤2C. This covers both
zero pressure and zero variance uniformly.

## A quantitative inverse for the original signed equation

The exact defining map is F(μ,p)=d₀ μ√V(μp).
Its [retained real derivative floor](../loop_modulation_bounds/README.md)
is F_μ≥d₀ exp(−4G)/32. Use m=d₀/(32E⁴).
At each real accepted base (μ₀,p₀), put

    τ=min(1,Rv/(4G+2)), M=2d₀(G+1)C.
    a=min(τ/4, mτ²/(64M)), P=4M/τ,
    σ=min(a, ma/[4(1+P)]).                               (2)

The μ,p polydisk of radius τ maps into the z disk of radius Rv/2,
because |Δ(μp)|≤2Gτ+τ². Therefore |F|≤M there.
On its half-radius polydisk Cauchy's formula gives
|F_μμ|≤8M/τ², |F_μp|≤4M/τ² and |F_p|≤P.
For |u|,|v|≤a, F_μ(μ₀+u,p₀+v) differs from the REAL
F_μ(μ₀,p₀) by at most 12Ma/τ²≤m/4.

For |v|,|w|<σ solve F(μ₀+u,p₀+v)=s₀+w, s₀=F(μ₀,p₀),
by dividing its linear u term by F_μ(μ₀,p₀).
The remaining u-dependent map has Lipschitz constant ≤1/4.
Its constant term is at most (|w|+P|v|)/m≤a/4.
Thus it maps |u|≤a into |u|≤a/2 and has one holomorphic fixed point.
On the real domain it agrees with the accepted signed root by the
original local uniqueness; local branches patch to that same root.
No inverse of V_μ, p, μ or s is introduced.

For normalized input Taylor coefficients of positive total degree,
Cauchy's formula consequently gives the convenient safe bound

    |∂p^i ∂s^j Ψ|/(i!j!) ≤ a (2/σ)^(i+j),  i+j≥1.       (3)

The value is bounded by G on the actual family. The radius σ is
only an estimating device in scalar input space: it is not a new
physical band selector or a change to the accepted low-order root.

## Finite slow composition and the remaining input

At any real slow center let P₊(h),S₊(h) be nonnegative polynomials
bounding the normalized, zero-constant Taylor increments of the ACTUAL
p,s through the desired downward-closed index set. For each nonzero α,

    |∂x^α μ|/α! ≤
       a Σ_(k=1)^|α| (2/σ)^k [h^α](P₊+S₊)^k.            (4)

This follows from (3) and
Σ_(i+j=k)P₊^i S₊^j≤(P₊+S₊)^k coefficientwise.
It is a finite real chain-rule identity: no infinite slow Taylor
series, analyticity across a flat cutoff or lower bound on s is needed.
Subindices remain within the requested set, preserving Q63's separation
of radial and η orders. Formula (4) replaces the independent μ-jet
input in Q64's mixed phase estimate.

The cutoff in s is the [explicit selected variance cutoff](../explicit_variance_cutoff/README.md),
not an arbitrary member of a smooth cutoff class. Where it is nonzero,
the radicand has the original positive gap δ_L/4 and a has its
retained floor; its jets follow by finite square-root, reciprocal and
cutoff composition. Flat extension covers the two transition boundaries.
The pressure input p₂ must include the full incoming five moments
and both 4AηΠ and −(1−η²)Πη terms of the forward map.
The estimate does not supply those input seminorms or permit a
pressure reset. Their pre-modulation regional enclosure and the
low-order normalized repair defect are the next quantitative obligation.
Uniform physical-pulse bounds, the radial pressure tail, complete
residual improvement and blow-up remain unproved.

Run test_audit.py and audit.py --verify-record evidence.json
--output /path/to/fresh.json. Manufactured exact tests check the disk
inequalities, zero-pressure Taylor factor and finite composition.
