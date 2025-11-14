#ifndef PHYSICS_ROTATIONAL_DYNAMICS_HPP
#define PHYSICS_ROTATIONAL_DYNAMICS_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace rotational_dynamics {

/**
 * @brief Rotational Dynamics
 *
 * Functions for analyzing rotational motion, analogous to linear dynamics.
 *
 * Analogies:
 * - Linear: F = ma  →  Rotational: τ = Iα
 * - Linear: p = mv  →  Rotational: L = Iω
 * - Linear: KE = (1/2)mv²  →  Rotational: KE = (1/2)Iω²
 *
 * Standard notation:
 * - θ (theta): angular position (radians)
 * - ω (omega): angular velocity (rad/s)
 * - α (alpha): angular acceleration (rad/s²)
 * - τ (tau): torque (N⋅m)
 * - I: moment of inertia (kg⋅m²)
 * - L: angular momentum (kg⋅m²/s)
 */

// ============================================================================
// Angular Acceleration
// ============================================================================

/**
 * @brief Calculate angular acceleration from change in angular velocity
 *
 * α = Δω/Δt
 *
 * @param initialAngularVel Initial angular velocity (in rad/s)
 * @param finalAngularVel Final angular velocity (in rad/s)
 * @param time Time interval (in seconds, must be > 0)
 * @return Angular acceleration (in rad/s²)
 * @throws std::invalid_argument if time <= 0
 */
inline double calculateAngularAcceleration(double initialAngularVel,
                                           double finalAngularVel, double time) {
    if (time <= 0) {
        throw std::invalid_argument("Time must be positive");
    }
    return (finalAngularVel - initialAngularVel) / time;
}

/**
 * @brief Calculate final angular velocity with constant angular acceleration
 *
 * ω_f = ω_i + αt (rotational analog of v = v₀ + at)
 *
 * @param initialAngularVel Initial angular velocity (in rad/s)
 * @param angularAccel Angular acceleration (in rad/s²)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Final angular velocity (in rad/s)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateFinalAngularVelocity(double initialAngularVel,
                                            double angularAccel, double time) {
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return initialAngularVel + angularAccel * time;
}

/**
 * @brief Calculate angular displacement with constant angular acceleration
 *
 * θ = ω₀t + (1/2)αt² (rotational analog of s = v₀t + (1/2)at²)
 *
 * @param initialAngularVel Initial angular velocity (in rad/s)
 * @param angularAccel Angular acceleration (in rad/s²)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Angular displacement (in radians)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateAngularDisplacement(double initialAngularVel,
                                           double angularAccel, double time) {
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return initialAngularVel * time + 0.5 * angularAccel * time * time;
}

/**
 * @brief Calculate tangential acceleration from angular acceleration
 *
 * a_tangential = rα
 *
 * @param radius Distance from axis of rotation (in meters, must be > 0)
 * @param angularAccel Angular acceleration (in rad/s²)
 * @return Tangential acceleration (in m/s²)
 * @throws std::invalid_argument if radius <= 0
 */
inline double calculateTangentialAcceleration(double radius, double angularAccel) {
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return radius * angularAccel;
}

// ============================================================================
// Change in Direction of Angular Velocity
// ============================================================================

/**
 * @brief Calculate angle turned through in time t
 *
 * For constant angular velocity: θ = ωt
 *
 * @param angularVelocity Angular velocity (in rad/s)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Angle turned (in radians)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateAngleTurned(double angularVelocity, double time) {
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return angularVelocity * time;
}

/**
 * @brief Calculate change in angular velocity vector magnitude
 *
 * When angular velocity changes direction (e.g., precession),
 * Δω_magnitude = 2ω sin(Δθ/2) for small angle changes
 *
 * @param angularVelocity Magnitude of angular velocity (in rad/s, must be >= 0)
 * @param angleChange Angle through which ω vector rotates (in radians)
 * @return Magnitude of change in angular velocity (in rad/s)
 * @throws std::invalid_argument if angularVelocity < 0
 */
inline double calculateAngularVelocityChange(double angularVelocity, double angleChange) {
    if (angularVelocity < 0) {
        throw std::invalid_argument("Angular velocity must be non-negative");
    }
    return 2.0 * angularVelocity * std::sin(angleChange / 2.0);
}

// ============================================================================
// Torque
// ============================================================================

/**
 * @brief Calculate torque from force and lever arm
 *
 * τ = r × F = rF sin(θ)
 * For perpendicular force: τ = rF
 *
 * @param radius Lever arm (perpendicular distance, in meters, must be > 0)
 * @param force Applied force (in Newtons, must be >= 0)
 * @return Torque (in N⋅m)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTorque(double radius, double force) {
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    if (force < 0) {
        throw std::invalid_argument("Force must be non-negative");
    }
    return radius * force;
}

/**
 * @brief Calculate torque with angle between r and F
 *
 * τ = rF sin(θ)
 *
 * @param radius Distance from axis (in meters, must be > 0)
 * @param force Applied force (in Newtons, must be >= 0)
 * @param angleRadians Angle between r and F (in radians)
 * @return Torque (in N⋅m)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTorqueWithAngle(double radius, double force, double angleRadians) {
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    if (force < 0) {
        throw std::invalid_argument("Force must be non-negative");
    }
    return radius * force * std::sin(angleRadians);
}

// ============================================================================
// Angular Acceleration Caused by Torque
// ============================================================================

/**
 * @brief Calculate angular acceleration from torque
 *
 * From τ = Iα, we get α = τ/I
 *
 * @param torque Net torque (in N⋅m)
 * @param momentOfInertia Moment of inertia (in kg⋅m², must be > 0)
 * @return Angular acceleration (in rad/s²)
 * @throws std::invalid_argument if momentOfInertia <= 0
 */
inline double calculateAngularAccelFromTorque(double torque, double momentOfInertia) {
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    return torque / momentOfInertia;
}

/**
 * @brief Calculate torque required for desired angular acceleration
 *
 * τ = Iα
 *
 * @param momentOfInertia Moment of inertia (in kg⋅m², must be > 0)
 * @param angularAccel Desired angular acceleration (in rad/s²)
 * @return Required torque (in N⋅m)
 * @throws std::invalid_argument if momentOfInertia <= 0
 */
inline double calculateRequiredTorque(double momentOfInertia, double angularAccel) {
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    return momentOfInertia * angularAccel;
}

// ============================================================================
// Angular Momentum
// ============================================================================

/**
 * @brief Calculate angular momentum
 *
 * L = Iω
 *
 * @param momentOfInertia Moment of inertia (in kg⋅m², must be > 0)
 * @param angularVelocity Angular velocity (in rad/s)
 * @return Angular momentum (in kg⋅m²/s)
 * @throws std::invalid_argument if momentOfInertia <= 0
 */
inline double calculateAngularMomentum(double momentOfInertia, double angularVelocity) {
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    return momentOfInertia * angularVelocity;
}

/**
 * @brief Calculate angular momentum for point mass
 *
 * L = mvr (for circular motion)
 *
 * @param mass Mass (in kilograms, must be > 0)
 * @param velocity Tangential velocity (in m/s)
 * @param radius Distance from axis (in meters, must be > 0)
 * @return Angular momentum (in kg⋅m²/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculatePointMassAngularMomentum(double mass, double velocity, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return mass * velocity * radius;
}

/**
 * @brief Calculate change in angular momentum
 *
 * ΔL = τ⋅Δt (impulse-momentum theorem for rotation)
 *
 * @param torque Net torque (in N⋅m)
 * @param time Time interval (in seconds, must be >= 0)
 * @return Change in angular momentum (in kg⋅m²/s)
 * @throws std::invalid_argument if time < 0
 */
inline double calculateAngularMomentumChange(double torque, double time) {
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return torque * time;
}

// ============================================================================
// Kinetic Energy of Rotating Body
// ============================================================================

/**
 * @brief Calculate rotational kinetic energy
 *
 * KE_rot = (1/2)Iω²
 *
 * @param momentOfInertia Moment of inertia (in kg⋅m², must be > 0)
 * @param angularVelocity Angular velocity (in rad/s)
 * @return Rotational kinetic energy (in Joules)
 * @throws std::invalid_argument if momentOfInertia <= 0
 */
inline double calculateRotationalKE(double momentOfInertia, double angularVelocity) {
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    return 0.5 * momentOfInertia * angularVelocity * angularVelocity;
}

/**
 * @brief Calculate angular velocity from rotational kinetic energy
 *
 * ω = √(2KE/I)
 *
 * @param kineticEnergy Rotational kinetic energy (in Joules, must be >= 0)
 * @param momentOfInertia Moment of inertia (in kg⋅m², must be > 0)
 * @return Angular velocity (in rad/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateAngularVelocityFromKE(double kineticEnergy, double momentOfInertia) {
    if (kineticEnergy < 0) {
        throw std::invalid_argument("Kinetic energy must be non-negative");
    }
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    return std::sqrt(2.0 * kineticEnergy / momentOfInertia);
}

// ============================================================================
// Moment of Inertia Formulas
// ============================================================================

/**
 * @brief Moment of inertia of point mass
 *
 * I = mr²
 *
 * @param mass Mass (in kilograms, must be > 0)
 * @param radius Distance from axis (in meters, must be >= 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if mass <= 0 or radius < 0
 */
inline double momentOfInertiaPointMass(double mass, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius < 0) {
        throw std::invalid_argument("Radius must be non-negative");
    }
    return mass * radius * radius;
}

/**
 * @brief Moment of inertia of thin rod about center
 *
 * I = (1/12)ML² (axis through center, perpendicular to rod)
 *
 * @param mass Mass of rod (in kilograms, must be > 0)
 * @param length Length of rod (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaRodCenter(double mass, double length) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    return (mass * length * length) / 12.0;
}

/**
 * @brief Moment of inertia of thin rod about end
 *
 * I = (1/3)ML² (axis through one end, perpendicular to rod)
 *
 * @param mass Mass of rod (in kilograms, must be > 0)
 * @param length Length of rod (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaRodEnd(double mass, double length) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    return (mass * length * length) / 3.0;
}

/**
 * @brief Moment of inertia of solid cylinder about central axis
 *
 * I = (1/2)MR²
 *
 * @param mass Mass of cylinder (in kilograms, must be > 0)
 * @param radius Radius of cylinder (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaCylinder(double mass, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return 0.5 * mass * radius * radius;
}

/**
 * @brief Moment of inertia of hollow cylinder (thin-walled)
 *
 * I = MR² (all mass at radius R)
 *
 * @param mass Mass of cylinder (in kilograms, must be > 0)
 * @param radius Radius of cylinder (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaHollowCylinder(double mass, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return mass * radius * radius;
}

/**
 * @brief Moment of inertia of solid sphere about diameter
 *
 * I = (2/5)MR²
 *
 * @param mass Mass of sphere (in kilograms, must be > 0)
 * @param radius Radius of sphere (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaSolidSphere(double mass, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return 0.4 * mass * radius * radius;
}

/**
 * @brief Moment of inertia of hollow sphere (thin-walled)
 *
 * I = (2/3)MR²
 *
 * @param mass Mass of sphere (in kilograms, must be > 0)
 * @param radius Radius of sphere (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaHollowSphere(double mass, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return (2.0 / 3.0) * mass * radius * radius;
}

/**
 * @brief Moment of inertia of thin disk about central axis
 *
 * I = (1/2)MR²
 *
 * @param mass Mass of disk (in kilograms, must be > 0)
 * @param radius Radius of disk (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaDisk(double mass, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return 0.5 * mass * radius * radius;
}

/**
 * @brief Moment of inertia of thin hoop about central axis
 *
 * I = MR²
 *
 * @param mass Mass of hoop (in kilograms, must be > 0)
 * @param radius Radius of hoop (in meters, must be > 0)
 * @return Moment of inertia (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double momentOfInertiaHoop(double mass, double radius) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return mass * radius * radius;
}

// ============================================================================
// Parallel Axis Theorem
// ============================================================================

/**
 * @brief Calculate moment of inertia about parallel axis
 *
 * Parallel Axis Theorem: I = I_cm + Md²
 * where I_cm is moment about center of mass,
 * d is distance between parallel axes
 *
 * @param momentAboutCM Moment of inertia about center of mass (in kg⋅m², must be >= 0)
 * @param mass Total mass (in kilograms, must be > 0)
 * @param distance Distance between parallel axes (in meters, must be >= 0)
 * @return Moment of inertia about parallel axis (in kg⋅m²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double parallelAxisTheorem(double momentAboutCM, double mass, double distance) {
    if (momentAboutCM < 0) {
        throw std::invalid_argument("Moment of inertia must be non-negative");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (distance < 0) {
        throw std::invalid_argument("Distance must be non-negative");
    }
    return momentAboutCM + mass * distance * distance;
}

// ============================================================================
// Compound Pendulum
// ============================================================================

/**
 * @brief Calculate period of compound pendulum (physical pendulum)
 *
 * T = 2π√(I/(mgd))
 * where I is moment of inertia about pivot,
 * d is distance from pivot to center of mass
 *
 * @param momentOfInertia Moment of inertia about pivot (in kg⋅m², must be > 0)
 * @param mass Total mass (in kilograms, must be > 0)
 * @param distanceToCM Distance from pivot to CM (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Period (in seconds)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateCompoundPendulumPeriod(double momentOfInertia, double mass,
                                              double distanceToCM, double gravity = 9.81) {
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (distanceToCM <= 0) {
        throw std::invalid_argument("Distance to CM must be positive");
    }
    return 2.0 * M_PI * std::sqrt(momentOfInertia / (mass * gravity * distanceToCM));
}

/**
 * @brief Calculate equivalent length of simple pendulum (compound pendulum)
 *
 * L_eq = I/(md) - length of simple pendulum with same period
 *
 * @param momentOfInertia Moment of inertia about pivot (in kg⋅m², must be > 0)
 * @param mass Total mass (in kilograms, must be > 0)
 * @param distanceToCM Distance from pivot to CM (in meters, must be > 0)
 * @return Equivalent length (in meters)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateEquivalentLength(double momentOfInertia, double mass, double distanceToCM) {
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (distanceToCM <= 0) {
        throw std::invalid_argument("Distance to CM must be positive");
    }
    return momentOfInertia / (mass * distanceToCM);
}

// ============================================================================
// Center of Percussion
// ============================================================================

/**
 * @brief Calculate center of percussion distance from pivot
 *
 * Center of percussion is the point where an impulse produces no reaction
 * at the pivot (sweet spot).
 *
 * d_percussion = I/(m⋅d_CM)
 * where d_CM is distance from pivot to center of mass
 *
 * @param momentOfInertia Moment of inertia about pivot (in kg⋅m², must be > 0)
 * @param mass Total mass (in kilograms, must be > 0)
 * @param distanceToCM Distance from pivot to CM (in meters, must be > 0)
 * @return Distance to center of percussion from pivot (in meters)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateCenterOfPercussion(double momentOfInertia, double mass,
                                         double distanceToCM) {
    if (momentOfInertia <= 0) {
        throw std::invalid_argument("Moment of inertia must be positive");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (distanceToCM <= 0) {
        throw std::invalid_argument("Distance to CM must be positive");
    }
    return momentOfInertia / (mass * distanceToCM);
}

/**
 * @brief Calculate center of percussion for uniform rod
 *
 * For uniform rod pivoted at one end:
 * Center of percussion is at distance (2/3)L from pivot
 *
 * @param length Length of rod (in meters, must be > 0)
 * @return Distance to center of percussion (in meters)
 * @throws std::invalid_argument if length <= 0
 */
inline double calculateRodCenterOfPercussion(double length) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    return (2.0 / 3.0) * length;
}

} // namespace rotational_dynamics
} // namespace physics

#endif // PHYSICS_ROTATIONAL_DYNAMICS_HPP
