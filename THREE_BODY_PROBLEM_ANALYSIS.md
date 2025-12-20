# Three-Body Problem: Mathematical Analysis & Library Mapping

## Overview

The gravitational three-body problem asks: Given three masses with initial positions and velocities, how do they move under mutual gravitational attraction?

**Key Result**: No general analytical solution exists (Poincaré, 1890s), but numerical integration reveals rich chaotic dynamics.

## Mathematical Formulation

### Hamiltonian Structure

For three bodies with masses m₁, m₂, m₃ at positions **r₁**, **r₂**, **r₃**:

**Hamiltonian**:
```
H(q, p) = T + V

Kinetic Energy:  T = |p₁|²/(2m₁) + |p₂|²/(2m₂) + |p₃|²/(2m₃)

Potential Energy: V = -Gm₁m₂/|r₁₂| - Gm₂m₃/|r₂₃| - Gm₁m₃/|r₁₃|

where rᵢⱼ = rᵢ - rⱼ
```

**Degrees of Freedom**: 9 coordinates (3 bodies × 3D) → 18-dimensional phase space (q, p)

**Conservation Laws**:
- Total energy: H = constant
- Linear momentum: P = p₁ + p₂ + p₃ = constant
- Angular momentum: L = r₁×p₁ + r₂×p₂ + r₃×p₃ = constant

These reduce the effective dimensionality but still leave chaos possible.

---

## Proof: Existing Library Components Are Sufficient

### **Theorem**: The library contains all necessary mathematical machinery to simulate and analyze the three-body problem.

### **Proof**:

**1. Hamiltonian Formulation** ✓
   - **Required**: Phase space (q, p), Hamilton's equations
   - **Library provides**: `classical_hamiltonian.hpp`
     - `PhasePoint` structure for (q, p)
     - `hamiltonEquations()` computes (dq/dt, dp/dt)
     - Numerical differentiation ∂H/∂q, ∂H/∂p

**2. Gravitational Potential** ✓
   - **Required**: V = -Σᵢ<ⱼ Gmᵢmⱼ/|rᵢⱼ|
   - **Library provides**: `gravitation.hpp`, `orbital.hpp`
     - `universalGravitationForce(m1, m2, r)` → F = Gm₁m₂/r²
     - `calculateOrbitalPotentialEnergy(m, r)` → U = -GMm/r
     - Can be composed for three pairwise interactions

**3. Numerical Integration** ✓
   - **Required**: Solve Hamilton's equations numerically
   - **Library provides**: `ode_dynamical_systems.hpp`
     - **Runge-Kutta 4th order** (RK4): O(h⁴) accuracy
     - **Heun's method**: 2nd order predictor-corrector
     - **Euler method**: 1st order baseline
   - **Symplectic integrators**: `classical_hamiltonian.hpp`
     - `stepSymplecticEuler()`: Preserves symplectic structure exactly
     - `stepVerlet()`: 2nd order symplectic, conserves H better

**4. Chaos Analysis** ✓
   - **Required**: Identify chaotic behavior
   - **Library provides**: `ode_dynamical_systems.hpp`
     - Chapter 9: Higher dimensional chaos
     - Lyapunov exponent calculations
     - Poincaré sections
     - Phase space trajectory analysis

**5. Phase Space Tools** ✓
   - **Required**: Analyze 18D phase space
   - **Library provides**: `classical_phase_space.hpp`
     - Phase space density evolution
     - Trajectory tracking
     - Liouville theorem verification

∎ **Q.E.D.**: All mathematical components exist.

---

## Why The Three-Body Problem Is Hard

### **Poincaré's Insight** (1887-1889)

Even with 10 conserved quantities (from symmetries), the system has:
- 18 phase space dimensions
- 10 constraints
- **Effective dimension: 8**

**Result**: Still high-dimensional enough for chaos!

### **Demonstration of Chaos**

**Initial Condition Sensitivity**:
```
Let Δr₀ = 10⁻⁹ meters (atomic scale perturbation)

After time t, separation grows as:
Δr(t) ≈ Δr₀ · exp(λt)

where λ ≈ 0.1 to 1.0 per orbital period (typical Lyapunov exponent)
```

**Implication**: Prediction horizon is ~10-20 orbital periods, even with perfect arithmetic.

The library's Lyapunov exponent calculators can measure λ numerically.

---

## Special Solutions

Despite chaos, some periodic solutions exist:

### 1. **Lagrange's Equilateral Triangle Solution** (1772)
   - Three masses at vertices of rotating equilateral triangle
   - **Stability**: Stable if mass ratios satisfy specific inequalities
   - **Library check**: Hamiltonian is constant, system is periodic

### 2. **Euler's Collinear Solution** (1767)
   - Three masses on a rotating line
   - **Stability**: Always unstable (!)
   - **Library verification**: Lyapunov exponent > 0

### 3. **Figure-8 Orbit** (Moore 1993, Chenciner-Montgomery 2000)
   - Three equal masses chase each other in figure-8 pattern
   - **Period**: ~6.3 time units
   - **Library test**: H conserved to machine precision with symplectic integrator

---

## Numerical Validation Test

### Energy Conservation Check

For any three-body trajectory, the library can verify:

```
ΔH/H₀ = |H(t) - H(0)|/|H(0)| < tolerance

Theoretical requirement: ΔH = 0 (exact conservation)

RK4 (non-symplectic):      ΔH/H₀ ~ O(h⁴)     (accumulates)
Symplectic Verlet:         ΔH/H₀ ~ O(h²)     (oscillates, bounded)
```

**Verdict**: Symplectic integrators from `classical_hamiltonian.hpp` are superior for long-time integration.

**Mathematical proof**: Symplectic maps preserve phase space volume (Liouville's theorem), preventing artificial energy drift.

---

## Restricted Three-Body Problem

Special case: m₃ << m₁, m₂ (e.g., Sun-Earth-Satellite)

### Simplification

The small mass doesn't affect the two large masses:
- **Reduced to**: Test particle in time-dependent potential
- **Degrees of freedom**: 3 (instead of 9)
- **Still chaotic**: Due to time-dependent forcing

### Lagrange Points

Five equilibrium points L₁, L₂, L₃, L₄, L₅ exist where:

```
∇U_eff = 0    (in rotating frame)

where U_eff = -Gm₁/r₁ - Gm₂/r₂ - ½ω²|r_perp|²
```

**Library calculation**:
- Use `gravitationalFieldStrength()` for each mass
- Find zeros of gradient numerically
- Verify stability via linearization (Hamiltonian eigenvalues)

**Physical examples**:
- L₄, L₅: Trojan asteroids (stable)
- L₂: James Webb Space Telescope (unstable, requires station-keeping)

---

## Conclusion

### **Mathematical Completeness Statement**

The library provides:

1. ✓ **Hamiltonian framework** for formulating the problem
2. ✓ **Gravitational potentials** for physical forces
3. ✓ **Numerical integrators** (both standard and symplectic)
4. ✓ **Chaos diagnostics** (Lyapunov exponents, Poincaré sections)
5. ✓ **Conservation law verification** (energy, momentum)

### **What's Missing**

Only high-level orchestration:
- Combining 3 bodies' mutual interactions
- Setting up initial conditions
- Visualization of trajectories

But all **mathematical operations** exist in the library.

### **Historical Note**

The three-body problem has driven major mathematical developments:
- **Poincaré** (1890s): Discovered chaos, invented topology
- **KAM theory** (1954-1963): Persistence of quasi-periodic orbits
- **Numerical symplectic integration** (1980s): Structure-preserving algorithms

All these topics appear in the library:
- Chaos: `ode_dynamical_systems.hpp` (Chapter 9)
- KAM theorem: `ode_dynamical_systems.hpp` (mentioned in Hamiltonian systems)
- Symplectic methods: `classical_hamiltonian.hpp`

---

## References

**Library Components**:
- `include/physics/classical_hamiltonian.hpp` - Hamiltonian mechanics
- `include/maths/ode_dynamical_systems.hpp` - ODE solvers and chaos
- `include/physics/gravitation.hpp` - Universal gravitation
- `include/physics/orbital.hpp` - Two-body orbital mechanics
- `include/physics/classical_phase_space.hpp` - Phase space analysis

**Theoretical Background**:
- Poincaré, H. (1890). "Sur le problème des trois corps"
- Kolmogorov-Arnold-Moser (KAM) Theory (1954-1963)
- Chenciner, A., Montgomery, R. (2000). "A remarkable periodic solution of the three-body problem"

**Key Insight**: The two-body problem (`orbital.hpp`) is **integrable** (exact solutions). Adding one more body makes it **chaotic** (no general solution, requires numerical methods from `ode_dynamical_systems.hpp`).
