#ifndef PHYSICS_ADVANCED_FLUID_DYNAMICS_GOVERNING_EQUATIONS_HPP
#define PHYSICS_ADVANCED_FLUID_DYNAMICS_GOVERNING_EQUATIONS_HPP

#include "maths/vectors.hpp"
#include "maths/matrices.hpp"
#include <cmath>
#include <stdexcept>
#include <functional>

/**
 * @file governing_equations.hpp
 * @brief Fundamental equations of fluid dynamics
 *
 * Implements:
 * - Continuity Equation (mass conservation)
 * - Navier-Stokes Equations (momentum conservation)
 * - Euler Equations (inviscid flow)
 * - Bernoulli's Equation (energy conservation for inviscid flow)
 * - Energy Equation (first law of thermodynamics)
 */

namespace physics::advanced::fluid_dynamics {

/**
 * @struct FluidState
 * @brief Complete state of a fluid element
 */
struct FluidState {
    double density;              // ρ (kg/m³)
    maths::linear_algebra::Vector velocity;    // u (m/s)
    double pressure;             // p (Pa)
    double temperature;          // T (K)
    double internal_energy;      // e (J/kg)

    FluidState() : density(1.0), velocity(3),
                   pressure(101325.0), temperature(300.0),
                   internal_energy(0.0) {}
};

/**
 * @class ContinuityEquation
 * @brief Mass conservation equation
 *
 * ∂ρ/∂t + ∇·(ρu) = 0
 *
 * Incompressible form: ∇·u = 0
 */
class ContinuityEquation {
public:
    /**
     * @brief Check incompressibility constraint
     *
     * For incompressible flow: ∇·u = 0
     *
     * @param divergence Velocity divergence ∇·u
     * @param tolerance Numerical tolerance
     * @return true if flow is incompressible
     */
    static bool isIncompressible(double divergence, double tolerance = 1e-6) {
        return std::abs(divergence) < tolerance;
    }

    /**
     * @brief Calculate velocity divergence (∇·u) numerically
     *
     * ∇·u = ∂u/∂x + ∂v/∂y + ∂w/∂z
     *
     * @param velocity_field Function u(x,y,z)
     * @param position Point (x,y,z)
     * @param dx Finite difference step size
     * @return Divergence value
     */
    static double velocityDivergence(
        const std::function<maths::linear_algebra::Vector(const maths::linear_algebra::Vector&)>& velocity_field,
        const maths::linear_algebra::Vector& position,
        double dx = 1e-6) {

        double div = 0.0;

        for (int i = 0; i < 3; ++i) {
            maths::linear_algebra::Vector pos_plus = position;
            maths::linear_algebra::Vector pos_minus = position;
            pos_plus[i] += dx;
            pos_minus[i] -= dx;

            maths::linear_algebra::Vector u_plus = velocity_field(pos_plus);
            maths::linear_algebra::Vector u_minus = velocity_field(pos_minus);

            div += (u_plus[i] - u_minus[i]) / (2.0 * dx);
        }

        return div;
    }

    /**
     * @brief Time derivative of density (compressible continuity)
     *
     * ∂ρ/∂t = -∇·(ρu) = -ρ(∇·u) - u·∇ρ
     *
     * @param state Current fluid state
     * @param velocity_divergence ∇·u
     * @param density_gradient ∇ρ
     * @return dρ/dt
     */
    static double densityTimeDerivative(
        const FluidState& state,
        double velocity_divergence,
        const maths::linear_algebra::Vector& density_gradient) {

        return -state.density * velocity_divergence -
               state.velocity.dot(density_gradient);
    }

    /**
     * @brief Mass flux through a surface
     *
     * Flux = ρu·n̂ dA
     *
     * @param state Fluid state
     * @param normal Surface normal vector
     * @param area Surface area
     * @return Mass flux (kg/s)
     */
    static double massFlux(const FluidState& state,
                          const maths::linear_algebra::Vector& normal,
                          double area) {
        return state.density * state.velocity.dot(normal) * area;
    }
};

/**
 * @class NavierStokesEquations
 * @brief Momentum conservation for viscous fluids
 *
 * ρ(∂u/∂t + u·∇u) = -∇p + μ∇²u + ρg
 *
 * Incompressible form:
 * ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + g
 */
class NavierStokesEquations {
public:
    /**
     * @brief Calculate convective acceleration (u·∇)u
     *
     * (u·∇)u = u∂u/∂x + v∂u/∂y + w∂u/∂z (for u-component)
     *
     * @param velocity u vector
     * @param velocity_gradient ∇u (3x3 Jacobian)
     * @return Convective acceleration
     */
    static maths::linear_algebra::Vector convectiveAcceleration(
        const maths::linear_algebra::Vector& velocity,
        const maths::linear_algebra::Matrix& velocity_gradient) {

        return velocity_gradient * velocity;
    }

    /**
     * @brief Calculate viscous term μ∇²u
     *
     * For incompressible: ∇²u = ∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²
     *
     * @param laplacian Laplacian ∇²u
     * @param dynamic_viscosity μ (Pa·s)
     * @return Viscous force per unit mass
     */
    static maths::linear_algebra::Vector viscousTerm(
        const maths::linear_algebra::Vector& laplacian,
        double dynamic_viscosity,
        double density) {

        return laplacian * (dynamic_viscosity / density);
    }

    /**
     * @brief Calculate kinematic viscosity
     *
     * ν = μ/ρ
     *
     * @param dynamic_viscosity μ (Pa·s)
     * @param density ρ (kg/m³)
     * @return Kinematic viscosity ν (m²/s)
     */
    static double kinematicViscosity(double dynamic_viscosity, double density) {
        if (density <= 0.0) {
            throw std::invalid_argument("Density must be positive");
        }
        return dynamic_viscosity / density;
    }

    /**
     * @brief Complete Navier-Stokes acceleration
     *
     * Du/Dt = -∇p/ρ + ν∇²u + g
     *
     * @param pressure_gradient ∇p
     * @param laplacian ∇²u
     * @param density ρ
     * @param kinematic_viscosity ν
     * @param gravity g
     * @return Total acceleration
     */
    static maths::linear_algebra::Vector totalAcceleration(
        const maths::linear_algebra::Vector& pressure_gradient,
        const maths::linear_algebra::Vector& laplacian,
        double density,
        double kinematic_viscosity,
        const maths::linear_algebra::Vector& gravity = maths::linear_algebra::Vector({0, 0, -9.81})) {

        return pressure_gradient * (-1.0 / density) +
               laplacian * kinematic_viscosity +
               gravity;
    }

    /**
     * @brief Vorticity form of Navier-Stokes
     *
     * ∂ω/∂t + (u·∇)ω = (ω·∇)u + ν∇²ω
     *
     * where ω = ∇×u (vorticity)
     */
    static maths::linear_algebra::Vector vorticityEquation(
        const maths::linear_algebra::Vector& vorticity,
        const maths::linear_algebra::Vector& velocity,
        const maths::linear_algebra::Matrix& velocity_gradient,
        const maths::linear_algebra::Matrix& vorticity_gradient,
        const maths::linear_algebra::Vector& vorticity_laplacian,
        double kinematic_viscosity) {

        maths::linear_algebra::Vector convection = velocity_gradient * vorticity;
        maths::linear_algebra::Vector stretching = vorticity_gradient * velocity;
        maths::linear_algebra::Vector diffusion = vorticity_laplacian * kinematic_viscosity;

        return convection * (-1.0) + stretching + diffusion;
    }
};

/**
 * @class EulerEquations
 * @brief Inviscid flow equations (zero viscosity)
 *
 * ρ(∂u/∂t + u·∇u) = -∇p + ρg
 *
 * Euler equations are Navier-Stokes with μ = 0
 */
class EulerEquations {
public:
    /**
     * @brief Euler momentum equation (inviscid)
     *
     * Du/Dt = -∇p/ρ + g
     *
     * @param pressure_gradient ∇p
     * @param density ρ
     * @param gravity g
     * @return Acceleration
     */
    static maths::linear_algebra::Vector eulerAcceleration(
        const maths::linear_algebra::Vector& pressure_gradient,
        double density,
        const maths::linear_algebra::Vector& gravity = maths::linear_algebra::Vector({0, 0, -9.81})) {

        return pressure_gradient * (-1.0 / density) + gravity;
    }

    /**
     * @brief Check if flow can be treated as inviscid
     *
     * Inviscid approximation valid when Re >> 1 (except boundary layers)
     *
     * @param reynolds_number Re
     * @return true if inviscid approximation is valid
     */
    static bool isInviscidValid(double reynolds_number) {
        return reynolds_number > 1000.0;  // Rule of thumb
    }

    /**
     * @brief Euler equations in conservation form
     *
     * ∂U/∂t + ∂F/∂x + ∂G/∂y + ∂H/∂z = S
     *
     * where U = [ρ, ρu, ρv, ρw, ρE]^T
     */
    struct ConservationForm {
        maths::linear_algebra::Vector U;  // Conservative variables [ρ, ρu, ρv, ρw, ρE]
        maths::linear_algebra::Vector F;  // Flux in x-direction
        maths::linear_algebra::Vector G;  // Flux in y-direction
        maths::linear_algebra::Vector H;  // Flux in z-direction
        maths::linear_algebra::Vector S;  // Source terms

        ConservationForm() : U(5), F(5), G(5), H(5), S(5) {
            // All vectors are already zero-initialized by default
        }
    };

    /**
     * @brief Convert primitive to conservative variables
     */
    static maths::linear_algebra::Vector primitiveToConservative(const FluidState& state) {
        maths::linear_algebra::Vector U(5);
        double rho = state.density;
        const maths::linear_algebra::Vector& u = state.velocity;
        double u_norm = u.norm();
        double E = state.internal_energy + 0.5 * u_norm * u_norm;

        U[0] = rho;
        U[1] = rho * u[0];
        U[2] = rho * u[1];
        U[3] = rho * u[2];
        U[4] = rho * E;

        return U;
    }
};

/**
 * @class BernoulliEquation
 * @brief Energy conservation for steady, inviscid, incompressible flow
 *
 * p + ½ρu² + ρgh = constant along streamline
 *
 * Or in head form: p/(ρg) + u²/(2g) + h = constant
 */
class BernoulliEquation {
public:
    /**
     * @brief Calculate total pressure (dynamic + static)
     *
     * p₀ = p + ½ρu²
     *
     * @param static_pressure p
     * @param density ρ
     * @param velocity u
     * @return Total pressure p₀
     */
    static double totalPressure(double static_pressure,
                                double density,
                                const maths::linear_algebra::Vector& velocity) {
        double v_norm = velocity.norm();
        double dynamic_pressure = 0.5 * density * v_norm * v_norm;
        return static_pressure + dynamic_pressure;
    }

    /**
     * @brief Calculate dynamic pressure
     *
     * q = ½ρu²
     *
     * @param density ρ
     * @param velocity u
     * @return Dynamic pressure q
     */
    static double dynamicPressure(double density,
                                  const maths::linear_algebra::Vector& velocity) {
        double v_norm = velocity.norm();
        return 0.5 * density * v_norm * v_norm;
    }

    /**
     * @brief Calculate velocity from pressure difference
     *
     * u = √(2(p₁-p₂)/ρ)
     *
     * @param pressure_diff p₁ - p₂
     * @param density ρ
     * @return Velocity magnitude
     */
    static double velocityFromPressure(double pressure_diff, double density) {
        if (pressure_diff < 0.0) {
            throw std::invalid_argument("Pressure difference must be positive");
        }
        return std::sqrt(2.0 * pressure_diff / density);
    }

    /**
     * @brief Calculate pressure from velocity (Bernoulli)
     *
     * p₂ = p₁ + ½ρ(u₁² - u₂²) + ρg(h₁ - h₂)
     *
     * @param p1 Pressure at point 1
     * @param rho Density
     * @param v1 Velocity at point 1
     * @param v2 Velocity at point 2
     * @param h1 Height at point 1
     * @param h2 Height at point 2
     * @param g Gravity
     * @return Pressure at point 2
     */
    static double pressureFromVelocity(
        double p1, double rho,
        const maths::linear_algebra::Vector& v1,
        const maths::linear_algebra::Vector& v2,
        double h1, double h2,
        double g = 9.81) {

        double v1_norm = v1.norm();
        double v2_norm = v2.norm();
        double v1_sq = v1_norm * v1_norm;
        double v2_sq = v2_norm * v2_norm;

        return p1 + 0.5 * rho * (v1_sq - v2_sq) + rho * g * (h1 - h2);
    }

    /**
     * @brief Calculate Bernoulli constant (total head)
     *
     * H = p/(ρg) + u²/(2g) + z
     *
     * @param pressure p
     * @param density ρ
     * @param velocity u
     * @param height z
     * @param g Gravity
     * @return Total head H (m)
     */
    static double totalHead(double pressure, double density,
                           const maths::linear_algebra::Vector& velocity,
                           double height, double g = 9.81) {
        double v_norm = velocity.norm();
        return pressure / (density * g) +
               (v_norm * v_norm) / (2.0 * g) +
               height;
    }

    /**
     * @brief Check if Bernoulli equation is applicable
     *
     * Requirements:
     * - Steady flow
     * - Inviscid (Re >> 1)
     * - Incompressible (Ma << 1)
     * - Along streamline
     */
    struct ApplicabilityCheck {
        bool is_steady;
        bool is_inviscid;
        bool is_incompressible;
        bool along_streamline;

        bool isValid() const {
            return is_steady && is_inviscid &&
                   is_incompressible && along_streamline;
        }
    };
};

/**
 * @class EnergyEquation
 * @brief First law of thermodynamics for fluid flow
 *
 * ρ(∂e/∂t + u·∇e) = -p∇·u + ∇·(k∇T) + Φ
 *
 * where:
 * - e is internal energy
 * - k is thermal conductivity
 * - Φ is viscous dissipation
 */
class EnergyEquation {
public:
    /**
     * @brief Calculate enthalpy
     *
     * h = e + p/ρ
     *
     * @param internal_energy e (J/kg)
     * @param pressure p (Pa)
     * @param density ρ (kg/m³)
     * @return Enthalpy h (J/kg)
     */
    static double enthalpy(double internal_energy,
                          double pressure,
                          double density) {
        return internal_energy + pressure / density;
    }

    /**
     * @brief Calculate total energy
     *
     * E = e + ½u²
     *
     * @param internal_energy e
     * @param velocity u
     * @return Total energy E
     */
    static double totalEnergy(double internal_energy,
                             const maths::linear_algebra::Vector& velocity) {
        double v_norm = velocity.norm();
        return internal_energy + 0.5 * v_norm * v_norm;
    }

    /**
     * @brief Calculate total enthalpy (stagnation enthalpy)
     *
     * h₀ = h + ½u²
     *
     * @param enthalpy h
     * @param velocity u
     * @return Total enthalpy h₀
     */
    static double totalEnthalpy(double enthalpy,
                                const maths::linear_algebra::Vector& velocity) {
        double v_norm = velocity.norm();
        return enthalpy + 0.5 * v_norm * v_norm;
    }

    /**
     * @brief Thermal conduction term
     *
     * Q̇ = ∇·(k∇T)
     *
     * @param temperature_laplacian ∇²T
     * @param thermal_conductivity k (W/(m·K))
     * @return Heat addition rate per unit volume (W/m³)
     */
    static double thermalConduction(double temperature_laplacian,
                                   double thermal_conductivity) {
        return thermal_conductivity * temperature_laplacian;
    }

    /**
     * @brief Viscous dissipation function
     *
     * Φ = μ[(∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ)²/2 - (2/3)(∇·u)²]
     *
     * Converts mechanical energy to heat
     *
     * @param velocity_gradient ∇u
     * @param dynamic_viscosity μ
     * @return Dissipation rate Φ (W/m³)
     */
    static double viscousDissipation(
        const maths::linear_algebra::Matrix& velocity_gradient,
        double dynamic_viscosity) {

        // Strain rate tensor: Sᵢⱼ = (∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ)/2
        maths::linear_algebra::Matrix strain_rate =
            (velocity_gradient + velocity_gradient.transpose()) * 0.5;

        // Trace (divergence) - manual calculation
        double div = 0.0;
        for (int i = 0; i < 3; ++i) {
            div += velocity_gradient(i, i);
        }

        // Dissipation
        double dissipation = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dissipation += strain_rate(i,j) * strain_rate(i,j);
            }
        }
        dissipation *= 2.0 * dynamic_viscosity;

        // Subtract bulk viscosity term
        dissipation -= (2.0/3.0) * dynamic_viscosity * div * div;

        return dissipation;
    }

    /**
     * @brief Energy equation time derivative
     *
     * ∂(ρE)/∂t = -∇·(ρEu) - ∇·(pu) + ∇·(k∇T) + ∇·(τ·u)
     */
    static double energyTimeDerivative(
        const FluidState& state,
        const maths::linear_algebra::Vector& heat_flux,
        double viscous_dissipation,
        double work_by_pressure) {

        // Simplified form
        return -work_by_pressure + heat_flux.norm() + viscous_dissipation;
    }
};

} // namespace physics::advanced::fluid_dynamics

#endif // PHYSICS_ADVANCED_FLUID_DYNAMICS_GOVERNING_EQUATIONS_HPP
