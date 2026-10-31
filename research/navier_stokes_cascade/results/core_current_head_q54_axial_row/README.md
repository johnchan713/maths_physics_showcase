# Q54: separate the two leading channels of the final axial row

On the same canonical first-real prefix as Q51–Q53, write W for the
**original** complete normalized covariance entry W_rz, Q for identity
minus auxiliary Haar, I=c_i^-1 Nt^-1 on zero-Haar functions,
L=v_r·grad_y, j=d R^(d−1). All orders O_a below are conditional fixed
coefficient-jet estimates with the original weights, finite poles and
polynomial band losses. Use a common positive hull and the lower estimates
at the extra finite orders used by differentiation and inversion.

Let f₁=Q Ez_postsigned, before the latest fast/slow means, and let
T=A₁ f₁. With the unchanged signs and normalization, the result is

    Ez_final = T + ε² M⁵ j⁵ I² L⁵ QW + O_(3−4κ).                 (1)

The coefficient of T is +1. It is not part of the two-viscosity symbol,
and (1) does not assign it a small norm. On any open set off the axial
cutoff derivative support it vanishes exactly. No axial stress floor
or lower bound on physical phase evaluations follows from Q53.

## Keep both occurrences of the axial correction

At either fast step the exact axial update for input f is
z=−I(f−A₁ f), not −If. The temporal inverse commutes with the shifted
operator on this frozen chart (auxiliary translations and a slow cutoff).
Take first the projected original axial source
f₀=Q[D₁(W−Σz)+Dz(Wzz+pm₀)]. The target is slow.
The graph piece of D₁ is MjL; the slow and connection pieces on W
are O₁. With the [full original pressure](../core_current_head_q49_initial_flux_tail/README.md),
Dz(Wzz+pm₀) is O_(2−κ), and A₁ f₀ is O_(2−κ) by compact graph-flux
annihilation. Thus z₀=−Mj I L QW+O₁.

After subtracting c_i Nt z, the exact fast difference is
−ε∂T z + D₁ ΔSrz + Dz ΔAz − ε Δ₀ z, including changed pressure in
ΔAz. For a tangential increment at order H (radial companion H+1)
the only term below H+1−κ is −ε M² j² L² z: every other radial viscosity
term has at most one graph loss, including ∂R j; axial viscosity has
none. Linear radial flux costs at worst H+1−κ, because either the
base radial velocity or the radial companion supplies ε. Old-mean
products and the new square also have a radial factor there. ΔAz
retains 2Gz, both old/new axial correlations, the new square and the
full pressure change: the latter has order H before Dz, from the
[fast difference](../core_fast_row_gain_join/README.md) with the
[corrected T₀ pressure](../core_current_head_q47_full_pressure_correction/README.md).
This uses a difference of pressures on the prescribed state, not
smallness of the total pressure at H.

First set H=1−κ. The old slow difference is O_(2−2κ); the complete
signed covariance response is O_(2−κ) using Q51's target order two.
The initial leftover +A₁ f₀ is O_(2−κ), not zero. Hence

    f₁ = +ε M³ j³ I L³ QW + O_(2−2κ).

Crucially the [same-state axial potential estimate](../core_axial_potential_pressure_gain/README.md)
already gives T=O₂ by the **direct** low-order operator bound, at each
fixed required jet order. It suffices inside z₁=−I(f₁−T):
z₁=−ε M³ j³ I² L³ QW+O_(2−2κ). It does not suffice for the unmultiplied
+T in the final row. At H=2−3κ the other fast-difference terms have
order ≥3−4κ. Applying the two-graph viscosity to the velocity error costs
2κ and gains one power, also giving 3−4κ. The post-signed Haar row and
latest slow difference are at least 3−4κ (Q51 target two and the
[complete first-real rate groups](../core_initial_target_rate_gain/README.md);
the slow difference uses old means PLUS the latest fast increment).
The [exact final decomposition](../core_initial_final_row_obstruction/README.md)
now yields (1).

## Exact remaining inequality

The source of T is the post-signed state including its initial means;
T=A₁ Dz QΦ, Φ=−2G_R Ψ+γ²+Wzz_postsigned+pm_postsigned. Here γ=D₁Ψ
includes the actual axial correction; pm is the cumulative shifted
pressure with its local terms and moving compensator. Dz acts on the
whole product and pressure, before A₁. No commutation with the moving
cutoff is implied. The angular row has no such additive term. The
final-state radial pressure remainder still has a different source.

For ν=3−6κ and any fixed 0<δ≤κ, the regular axial row has order
ν+δ **if and only if** T does, on the same weighted domain, allowing
finite polynomial losses. This follows by subtracting (1): the
viscous term is O_(3−5κ) already. For an order above 3−5κ it is
the **sum** of T and the viscous term that must be tested. This
equivalence isolates a physical-input obligation; it is not a proof
of the large-jet bounds. A sufficient route at δ=κ uses the
[shifted Fourier estimate](../core_current_head_q48_original_axial_tail/README.md)
on H=Dz QΦ, of order two: p=99995 integrations, raw native H through
n+2p+1=n+199991 and Φ through n+199992. Use the actual full compact
support, zero extension, native coordinate conversions, transition
floor, and high input pole on this post-signed state. Only the
nonzero auxiliary modes gain M^(-p); M=c_i^ρ ε^(-κ) contributes
a fixed polynomial floor cost. These counts refer to this tail only.

Next bound this post-signed tail directly or validate those same-state
inputs, and then treat the distinct radial A₀ remainder. The test
fields check signs, perturbations and rates; they do not enclose an
actual pulse, prove a complete residual estimate, or establish blow-up.
