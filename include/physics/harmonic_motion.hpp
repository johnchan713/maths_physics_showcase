#ifndef PHYSICS_HARMONIC_MOTION_HPP
#define PHYSICS_HARMONIC_MOTION_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace harmonic_motion {

/**
 * @brief Simple Harmonic Motion (SHM)
 *
 * Functions for analyzing oscillatory motion where the restoring force
 * is proportional to displacement: F = -kx
 *
 * Standard equations:
 * - Displacement: x(t) = A cos(ωt + φ)
 * - Velocity: v(t) = -Aω sin(ωt + φ)
 * - Acceleration: a(t) = -Aω² cos(ωt + φ) = -ω²x
 * - Force: F = -kx = -mω²x
 * - Angular frequency: ω = √(k/m) = 2πf = 2π/T
 *
 * Standard notation:
 * - A: amplitude (maximum displacement)
 * - ω (omega): angular frequency (rad/s)
 * - f: frequency (Hz)
 * - T: period (seconds)
 * - φ (phi): phase constant
 * - k: spring constant or restoring force constant
 * - m: mass
 */

// ============================================================================
// Acceleration in Simple Harmonic Motion
// ============================================================================

/**
 * @brief Calculate acceleration in SHM from displacement
 *
 * Acceleration in SHM: a = -ω²x
 * The negative sign indicates acceleration is always toward equilibrium
 *
 * @param angularFrequency Angular frequency ω (in rad/s, must be > 0)
 * @param displacement Displacement from equilibrium (in meters)
 * @return Acceleration (in m/s², negative when displaced in positive direction)
 * @throws std::invalid_argument if angularFrequency <= 0
 */
inline double calculateAcceleration(double angularFrequency, double displacement) {
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    return -angularFrequency * angularFrequency * displacement;
}

/**
 * @brief Calculate maximum acceleration in SHM
 *
 * Maximum acceleration occurs at maximum displacement:
 * a_max = ω²A
 *
 * @param angularFrequency Angular frequency ω (in rad/s, must be > 0)
 * @param amplitude Maximum displacement A (in meters, must be > 0)
 * @return Maximum acceleration magnitude (in m/s²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateMaxAcceleration(double angularFrequency, double amplitude) {
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    return angularFrequency * angularFrequency * amplitude;
}

/**
 * @brief Calculate acceleration at time t
 *
 * a(t) = -Aω² cos(ωt + φ)
 *
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @param angularFrequency Angular frequency (in rad/s, must be > 0)
 * @param time Time (in seconds, must be >= 0)
 * @param phaseConstant Phase constant (in radians, default: 0)
 * @return Acceleration at time t (in m/s²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateAccelerationAtTime(double amplitude, double angularFrequency,
                                          double time, double phaseConstant = 0.0) {
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return -amplitude * angularFrequency * angularFrequency *
           std::cos(angularFrequency * time + phaseConstant);
}

// ============================================================================
// Force in Simple Harmonic Motion
// ============================================================================

/**
 * @brief Calculate restoring force in SHM (Hooke's Law)
 *
 * F = -kx
 * The negative sign indicates force is always directed toward equilibrium
 *
 * @param springConstant Spring constant k (in N/m, must be > 0)
 * @param displacement Displacement from equilibrium (in meters)
 * @return Restoring force (in Newtons, negative when displaced positively)
 * @throws std::invalid_argument if springConstant <= 0
 */
inline double calculateRestoringForce(double springConstant, double displacement) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    return -springConstant * displacement;
}

/**
 * @brief Calculate maximum force in SHM
 *
 * F_max = kA
 * Maximum force occurs at maximum displacement
 *
 * @param springConstant Spring constant k (in N/m, must be > 0)
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @return Maximum force magnitude (in Newtons)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateMaxForce(double springConstant, double amplitude) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    return springConstant * amplitude;
}

/**
 * @brief Calculate force using mass and angular frequency
 *
 * Since k = mω², we have F = -mω²x
 *
 * @param mass Mass (in kilograms, must be > 0)
 * @param angularFrequency Angular frequency (in rad/s, must be > 0)
 * @param displacement Displacement from equilibrium (in meters)
 * @return Restoring force (in Newtons)
 * @throws std::invalid_argument if mass <= 0 or angularFrequency <= 0
 */
inline double calculateForceFromMass(double mass, double angularFrequency, double displacement) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    return -mass * angularFrequency * angularFrequency * displacement;
}

// ============================================================================
// Angular Frequency and Period
// ============================================================================

/**
 * @brief Calculate angular frequency from spring constant and mass
 *
 * ω = √(k/m)
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param mass Mass (in kilograms, must be > 0)
 * @return Angular frequency (in rad/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateAngularFrequency(double springConstant, double mass) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    return std::sqrt(springConstant / mass);
}

/**
 * @brief Calculate period of oscillation
 *
 * T = 2π/ω = 2π√(m/k)
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param mass Mass (in kilograms, must be > 0)
 * @return Period (in seconds)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculatePeriod(double springConstant, double mass) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    return 2.0 * M_PI * std::sqrt(mass / springConstant);
}

/**
 * @brief Calculate frequency of oscillation
 *
 * f = 1/T = ω/(2π) = (1/2π)√(k/m)
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param mass Mass (in kilograms, must be > 0)
 * @return Frequency (in Hz)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateFrequency(double springConstant, double mass) {
    double omega = calculateAngularFrequency(springConstant, mass);
    return omega / (2.0 * M_PI);
}

// ============================================================================
// Energy of a Vibrating Mass
// ============================================================================

/**
 * @brief Calculate kinetic energy at displacement x
 *
 * For SHM: KE = (1/2)mv² = (1/2)k(A² - x²)
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @param displacement Current displacement (in meters, must be |x| <= A)
 * @return Kinetic energy (in Joules)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateKineticEnergy(double springConstant, double amplitude, double displacement) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    if (std::abs(displacement) > amplitude) {
        throw std::invalid_argument("Displacement cannot exceed amplitude");
    }
    return 0.5 * springConstant * (amplitude * amplitude - displacement * displacement);
}

/**
 * @brief Calculate potential energy at displacement x
 *
 * PE = (1/2)kx²
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param displacement Displacement from equilibrium (in meters)
 * @return Potential energy (in Joules)
 * @throws std::invalid_argument if springConstant <= 0
 */
inline double calculatePotentialEnergy(double springConstant, double displacement) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    return 0.5 * springConstant * displacement * displacement;
}

/**
 * @brief Calculate total mechanical energy of oscillator
 *
 * E_total = (1/2)kA² = constant
 * Total energy equals maximum PE or maximum KE
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @return Total mechanical energy (in Joules)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateTotalEnergy(double springConstant, double amplitude) {
    if (springConstant <= 0) {
        throw std::invalid_argument("Spring constant must be positive");
    }
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    return 0.5 * springConstant * amplitude * amplitude;
}

/**
 * @brief Verify energy conservation in SHM
 *
 * At any point: E_total = KE + PE = (1/2)kA²
 *
 * @param springConstant Spring constant (in N/m, must be > 0)
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @param displacement Current displacement (in meters)
 * @return True if energy is conserved within tolerance
 */
inline bool verifyEnergyConservation(double springConstant, double amplitude,
                                     double displacement, double tolerance = 1e-6) {
    double totalEnergy = calculateTotalEnergy(springConstant, amplitude);
    double ke = calculateKineticEnergy(springConstant, amplitude, displacement);
    double pe = calculatePotentialEnergy(springConstant, displacement);
    return std::abs((ke + pe) - totalEnergy) < tolerance;
}

// ============================================================================
// Displacement and Velocity
// ============================================================================

/**
 * @brief Calculate displacement at time t
 *
 * x(t) = A cos(ωt + φ)
 *
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @param angularFrequency Angular frequency (in rad/s, must be > 0)
 * @param time Time (in seconds, must be >= 0)
 * @param phaseConstant Phase constant (in radians, default: 0)
 * @return Displacement at time t (in meters)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateDisplacement(double amplitude, double angularFrequency,
                                    double time, double phaseConstant = 0.0) {
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return amplitude * std::cos(angularFrequency * time + phaseConstant);
}

/**
 * @brief Calculate velocity at time t
 *
 * v(t) = -Aω sin(ωt + φ)
 *
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @param angularFrequency Angular frequency (in rad/s, must be > 0)
 * @param time Time (in seconds, must be >= 0)
 * @param phaseConstant Phase constant (in radians, default: 0)
 * @return Velocity at time t (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateVelocity(double amplitude, double angularFrequency,
                                double time, double phaseConstant = 0.0) {
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    if (time < 0) {
        throw std::invalid_argument("Time must be non-negative");
    }
    return -amplitude * angularFrequency *
           std::sin(angularFrequency * time + phaseConstant);
}

/**
 * @brief Calculate maximum velocity
 *
 * v_max = Aω
 * Maximum velocity occurs at equilibrium position
 *
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @param angularFrequency Angular frequency (in rad/s, must be > 0)
 * @return Maximum velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateMaxVelocity(double amplitude, double angularFrequency) {
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    return amplitude * angularFrequency;
}

/**
 * @brief Calculate velocity at displacement x
 *
 * From energy: v = ±ω√(A² - x²)
 *
 * @param angularFrequency Angular frequency (in rad/s, must be > 0)
 * @param amplitude Maximum displacement (in meters, must be > 0)
 * @param displacement Current displacement (in meters, |x| <= A)
 * @return Velocity magnitude at displacement x (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateVelocityAtDisplacement(double angularFrequency,
                                              double amplitude, double displacement) {
    if (angularFrequency <= 0) {
        throw std::invalid_argument("Angular frequency must be positive");
    }
    if (amplitude <= 0) {
        throw std::invalid_argument("Amplitude must be positive");
    }
    if (std::abs(displacement) > amplitude) {
        throw std::invalid_argument("Displacement cannot exceed amplitude");
    }
    return angularFrequency * std::sqrt(amplitude * amplitude - displacement * displacement);
}

// ============================================================================
// Simple Pendulum
// ============================================================================

/**
 * @brief Calculate period of simple pendulum
 *
 * For small angles: T = 2π√(L/g)
 * Period is independent of mass and amplitude (for small angles)
 *
 * @param length Length of pendulum (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Period (in seconds)
 * @throws std::invalid_argument if length <= 0
 */
inline double calculatePendulumPeriod(double length, double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    return 2.0 * M_PI * std::sqrt(length / gravity);
}

/**
 * @brief Calculate frequency of simple pendulum
 *
 * f = 1/T = (1/2π)√(g/L)
 *
 * @param length Length of pendulum (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Frequency (in Hz)
 * @throws std::invalid_argument if length <= 0
 */
inline double calculatePendulumFrequency(double length, double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    return (1.0 / (2.0 * M_PI)) * std::sqrt(gravity / length);
}

/**
 * @brief Calculate angular frequency of simple pendulum
 *
 * ω = √(g/L)
 *
 * @param length Length of pendulum (in meters, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Angular frequency (in rad/s)
 * @throws std::invalid_argument if length <= 0
 */
inline double calculatePendulumAngularFrequency(double length, double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    return std::sqrt(gravity / length);
}

/**
 * @brief Calculate length of pendulum for desired period
 *
 * From T = 2π√(L/g), we get L = gT²/(4π²)
 *
 * @param period Desired period (in seconds, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Required length (in meters)
 * @throws std::invalid_argument if period <= 0
 */
inline double calculatePendulumLength(double period, double gravity = 9.81) {
    if (period <= 0) {
        throw std::invalid_argument("Period must be positive");
    }
    return (gravity * period * period) / (4.0 * M_PI * M_PI);
}

/**
 * @brief Calculate maximum velocity of pendulum
 *
 * For small angle θ₀: v_max = θ₀√(gL)
 * where θ₀ is initial angular displacement in radians
 *
 * @param length Length of pendulum (in meters, must be > 0)
 * @param angularAmplitude Initial angular displacement (in radians, must be > 0, typically < 0.2)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Maximum velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculatePendulumMaxVelocity(double length, double angularAmplitude,
                                           double gravity = 9.81) {
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (angularAmplitude <= 0) {
        throw std::invalid_argument("Angular amplitude must be positive");
    }
    return angularAmplitude * std::sqrt(gravity * length);
}

/**
 * @brief Calculate total energy of pendulum
 *
 * For small angles: E = (1/2)mgLθ₀²
 *
 * @param mass Mass of pendulum bob (in kilograms, must be > 0)
 * @param length Length of pendulum (in meters, must be > 0)
 * @param angularAmplitude Maximum angular displacement (in radians, must be > 0)
 * @param gravity Gravitational acceleration (in m/s², default: 9.81)
 * @return Total energy (in Joules)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculatePendulumEnergy(double mass, double length, double angularAmplitude,
                                      double gravity = 9.81) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (length <= 0) {
        throw std::invalid_argument("Length must be positive");
    }
    if (angularAmplitude <= 0) {
        throw std::invalid_argument("Angular amplitude must be positive");
    }
    return 0.5 * mass * gravity * length * angularAmplitude * angularAmplitude;
}

} // namespace harmonic_motion
} // namespace physics

#endif // PHYSICS_HARMONIC_MOTION_HPP
