# Three-Body Problem: Pure Mathematical Formulation

## Plain English Explanation

**How the equations work together:** Start with initial positions and velocities for all three bodies. At each time step: (1) Calculate gravitational forces between each pair using F = Gm₁m₂/r², (2) Sum the forces on each body, (3) Use F = ma to get accelerations, (4) Update velocities (v_new = v_old + a·Δt), (5) Update positions (r_new = r_old + v·Δt). Repeat these steps millions of times. The Hamiltonian H = T + V provides energy conservation checks. Symplectic integrators ensure the geometry of phase space is preserved, preventing numerical errors from accumulating. The result: complete trajectories showing how the three bodies orbit, scatter, or collide.

**Why Lyapunov exponents matter:** The Lyapunov exponent λ measures how fast nearby trajectories diverge. We use it because even though we can calculate trajectories numerically, we need to know if our solution is predictable or chaotic. It works by comparing two nearly identical starting conditions: if separation grows exponentially as e^(λt), then λ is the Lyapunov exponent. When λ > 0, the system is chaotic—tiny measurement errors double every 1/λ time units, making long-term prediction impossible. When λ = 0, motion is regular and predictable. For the three-body problem, typical λ ≈ 0.1-1.0 per orbit, meaning prediction breaks down after ~10-20 orbits regardless of computational precision. The crucial property: λ connects directly to the Hamiltonian through phase space geometry. Symplectic systems have paired exponents (λᵢ, -λᵢ) that sum to zero, reflecting time-reversibility. This tells us the three-body problem is fundamentally unpredictable, not due to computational limits, but due to the mathematics itself.

---

## System Definition

**Three masses**: m₁, m₂, m₃
**Position vectors**: **r**₁, **r**₂, **r**₃ ∈ ℝ³
**Momentum vectors**: **p**₁, **p**₂, **p**₃ ∈ ℝ³
**Phase space**: Γ = {(**r**₁, **r**₂, **r**₃, **p**₁, **p**₂, **p**₃)} ∈ ℝ¹⁸

---

## Hamiltonian

```
H: ℝ¹⁸ → ℝ

H(**r**, **p**) = T(**p**) + V(**r**)
```

**Kinetic Energy**:
```
T = Σᵢ₌₁³ |**pᵢ**|²/(2mᵢ) = |**p₁**|²/(2m₁) + |**p₂**|²/(2m₂) + |**p₃**|²/(2m₃)
```

**Potential Energy**:
```
V = - Σᵢ<ⱼ (Gmᵢmⱼ)/|**rᵢⱼ**|

  = -Gm₁m₂/|**r₁** - **r₂**| - Gm₂m₃/|**r₂** - **r₃**| - Gm₁m₃/|**r₁** - **r₃**|
```

---

## Hamilton's Equations

```
d**rᵢ**/dt = ∂H/∂**pᵢ** = **pᵢ**/mᵢ                    (i = 1,2,3)

d**pᵢ**/dt = -∂H/∂**rᵢ** = -∇ᵣᵢV                      (i = 1,2,3)
```

**Gravitational Force**:
```
∇ᵣᵢV = Σⱼ≠ᵢ Gmᵢmⱼ(**rᵢ** - **rⱼ**)/|**rᵢ** - **rⱼ**|³

**Fᵢ** = -∇ᵣᵢV
```

**Explicit form**:
```
mᵢ d²**rᵢ**/dt² = Σⱼ≠ᵢ Gmᵢmⱼ(**rⱼ** - **rᵢ**)/|**rᵢⱼ**|³
```

---

## Conservation Laws

**Total Energy**:
```
E = H(**r**, **p**) = const.

dE/dt = {H, H} = 0
```

**Linear Momentum**:
```
**P** = Σᵢ₌₁³ **pᵢ** = const.

d**P**/dt = 0
```

**Angular Momentum**:
```
**L** = Σᵢ₌₁³ **rᵢ** × **pᵢ** = const.

d**L**/dt = 0
```

**Center of Mass**:
```
**R**_cm = (Σᵢ mᵢ**rᵢ**)/(Σᵢ mᵢ)

**P** = M**V**_cm,  where M = m₁ + m₂ + m₃
```

---

## Poisson Bracket Formulation

For observables A, B:
```
{A, B} = Σᵢ₌₁³ (∂A/∂**rᵢ** · ∂B/∂**pᵢ** - ∂A/∂**pᵢ** · ∂B/∂**rᵢ**)
```

**Canonical commutation relations**:
```
{rᵢᵅ, pⱼᵝ} = δᵢⱼδᵅᵝ
{rᵢᵅ, rⱼᵝ} = 0
{pᵢᵅ, pⱼᵝ} = 0
```

**Time evolution**:
```
dA/dt = {A, H} + ∂A/∂t
```

---

## Symplectic Structure

**Phase space form**:
```
ω = Σᵢ₌₁³ d**pᵢ** ∧ d**rᵢ**
```

**Liouville's theorem**:
```
dω/dt = 0

∇ · (**ṙ**, **ṗ**) = Σᵢ (∂ṙᵢ/∂rᵢ + ∂ṗᵢ/∂pᵢ) = 0
```

**Poincaré invariant**:
```
∫_Ω ω^n = ∫_φₜ(Ω) ω^n,  where n = 9
```

---

## Integrability vs Chaos

**Degrees of freedom**: N = 9

**Known integrals**:
```
I₁ = E              (energy)
I₂, I₃, I₄ = Pₓ, Pᵧ, Pᵤ   (linear momentum)
I₅, I₆, I₇ = Lₓ, Lᵧ, Lᵤ   (angular momentum)
I₈, I₉, I₁₀ = R_cm components
```

**Liouville-Arnold theorem**: System is integrable iff ∃ N independent integrals in involution.

**Three-body problem**:
```
# of integrals: 10
# needed: 9
BUT: Not all in involution ⇒ Non-integrable
```

**Result**: Chaotic dynamics for generic initial conditions.

---

## Jacobi Coordinates

**Reduce to relative motion**:

```
**ρ**₁ = **r₁** - **r₂**
**ρ**₂ = **r₃** - (**r₁** + **r₂**)/2
**R** = **R**_cm

Conjugate momenta: **π**₁, **π**₂, **P**
```

**Reduced Hamiltonian** (center-of-mass frame, **P** = 0):
```
H_red = |**π**₁|²/(2μ₁₂) + |**π**₂|²/(2μ₃,₁₂) + V(**ρ**₁, **ρ**₂)

μ₁₂ = m₁m₂/(m₁ + m₂)
μ₃,₁₂ = m₃(m₁ + m₂)/(m₁ + m₂ + m₃)
```

**Phase space dimension**: 18 → 12 (reduced by 6)

---

## Numerical Integration

**Symplectic Euler** (1st order):
```
**pᵢ**^(n+1) = **pᵢ**^(n) + Δt · **Fᵢ**(**r**^(n))
**rᵢ**^(n+1) = **rᵢ**^(n) + Δt · **pᵢ**^(n+1)/mᵢ
```

**Störmer-Verlet** (2nd order):
```
**rᵢ**^(n+1/2) = **rᵢ**^(n) + (Δt/2) · **pᵢ**^(n)/mᵢ
**pᵢ**^(n+1) = **pᵢ**^(n) + Δt · **Fᵢ**(**r**^(n+1/2))
**rᵢ**^(n+1) = **rᵢ**^(n+1/2) + (Δt/2) · **pᵢ**^(n+1)/mᵢ
```

**Symplectic property**:
```
ω(**ż**, **ż**') = ω(**z**, **z**'),  where **z** = (**r**, **p**)

det(∂**z**^(n+1)/∂**z**^(n)) = 1
```

---

## Lyapunov Exponent

**Tangent space evolution**:
```
d/dt(δ**z**) = J(**z**) · δ**z**

where J = (∂**ż**/∂**z**) is Jacobian
```

**Maximal Lyapunov exponent**:
```
λ_max = lim_(t→∞) (1/t) ln(|δ**z**(t)|/|δ**z**(0)|)
```

**Chaos criterion**:
```
λ_max > 0 ⇒ Chaotic
λ_max = 0 ⇒ Quasi-periodic
λ_max < 0 ⇒ Stable
```

**Spectrum** (ordered): λ₁ ≥ λ₂ ≥ ... ≥ λ₁₈

**Symplectic constraint**: λᵢ + λ₁₉₋ᵢ = 0

---

## Special Periodic Solutions

### Lagrange Equilateral (1772)

**Configuration**: |**r₁₂**| = |**r₂₃**| = |**r₁₃**| = R

**Equations**:
```
ω² = G(m₁ + m₂ + m₃)/R³

**rᵢ**(t) = R · (cos(ωt + φᵢ), sin(ωt + φᵢ))
φ₁ = 0, φ₂ = 2π/3, φ₃ = 4π/3
```

**Stability**: Stable iff
```
27(m₁m₂ + m₂m₃ + m₃m₁) > (m₁ + m₂ + m₃)²
```

### Euler Collinear (1767)

**Configuration**: **r₁**, **r₂**, **r₃** collinear, rotating

**Always unstable**: λ_max > 0

### Figure-8 (Chenciner-Montgomery, 2000)

**Equal masses**: m₁ = m₂ = m₃ = m

**Symmetry**: ℤ₃ × ℤ₂

**Period**: T ≈ 6.32591398

**Topological constraint**: Braiding in configuration space

---

## Restricted Three-Body Problem

**Limit**: m₃ → 0 (test particle)

**Circular case**: m₁, m₂ in circular orbits

**Effective potential** (rotating frame):
```
U_eff = -Gm₁/r₁ - Gm₂/r₂ - ½ω²(x² + y²)

where ω² = G(m₁ + m₂)/a³, a = |**r₁** - **r₂**|
```

**Lagrange points**: ∇U_eff = 0

```
L₁, L₂, L₃: Collinear (unstable)
L₄, L₅: Equilateral triangle (stable if m₂/m₁ < 0.0385)
```

**Jacobi integral**:
```
C_J = 2U_eff - (ẋ² + ẏ² + ż²) = const.
```

**Zero-velocity curves**: U_eff(x, y, z) = C_J/2

---

## Poincaré Section

**Surface of section**: Σ ⊂ Γ

**Poincaré map**:
```
P: Σ → Σ
**z**_n+1 = P(**z**_n)
```

**Example**: Σ = {(**r**, **p**) : z₃ = 0, ṗ_z₃ > 0}

**Structure**:
- Integrable: Smooth curves (KAM tori)
- Chaotic: Scattered points (homoclinic tangle)

---

## Library Component Mapping

```
Hamiltonian H                  → classical_hamiltonian.hpp::HamiltonianSystem
Hamilton's equations           → hamiltonEquations()
∂H/∂q, ∂H/∂p                  → dH_dq(), dH_dp()
Gravitational force Gmᵢmⱼ/r²  → gravitation.hpp::universalGravitationForce()
Phase point (**r**, **p**)     → PhasePoint struct
Symplectic integrator          → stepSymplecticEuler(), stepVerlet()
RK4 integrator                 → ode_dynamical_systems.hpp::rk4()
Lyapunov exponent λ            → ode_dynamical_systems.hpp (Chapter 9)
Poisson bracket {·,·}          → classical_hamiltonian.hpp
Conservation: dH/dt = 0        → Energy monitoring
Phase space volume             → classical_phase_space.hpp
```

---

## Existence & Uniqueness

**Theorem** (Picard-Lindelöf):

For Lipschitz **F**:
```
d**z**/dt = **F**(**z**, t)
**z**(t₀) = **z**₀
```

∃! local solution **z**(t) for t ∈ [t₀ - ε, t₀ + ε]

**Three-body problem**:
```
**F** = (**p**/m, -∇V)

Singularities: |**rᵢⱼ**| = 0 (collisions)
```

**Generic solutions**: Exist globally except at collisions

**Sundman's theorem** (1913): No triple collisions in finite time for non-zero angular momentum.

---

## Summary

**System**: 9 DOF, 18-dimensional phase space
**Integrals**: 10 known (E, **P**, **L**, **R**)
**Result**: Non-integrable → Chaos
**Solutions**: Numerical (symplectic preferred)
**Library**: Complete coverage of all mathematical operations
