#ifndef PHYSICS_GRAVITATION_HPP
#define PHYSICS_GRAVITATION_HPP

#include <cmath>
#include <stdexcept>

namespace physics {
namespace gravitation {

/**
 * @brief Universal Gravitation and Kepler's Laws
 *
 * Functions for analyzing gravitational interactions and orbital mechanics,
 * including Newton's synthesis of terrestrial and celestial mechanics.
 *
 * Historical note: Newton connected the fall of an apple with the Moon's
 * orbit, showing that the same gravitational force governs both.
 */

// ============================================================================
// Physical Constants
// ============================================================================

namespace constants {
    constexpr double G = 6.674e-11;              // Gravitational constant (N⋅m²/kg²)
    constexpr double EARTH_MASS = 5.972e24;      // kg
    constexpr double EARTH_RADIUS = 6.371e6;     // meters
    constexpr double MOON_MASS = 7.342e22;       // kg
    constexpr double MOON_ORBITAL_RADIUS = 3.844e8; // meters
    constexpr double MOON_ORBITAL_PERIOD = 2.36e6;  // seconds (27.3 days)
    constexpr double SUN_MASS = 1.989e30;        // kg
}

// ============================================================================
// Universal Gravitation
// ============================================================================

/**
 * @brief Newton's Law of Universal Gravitation
 *
 * F = G(m₁m₂)/r²
 *
 * Every particle in the universe attracts every other particle with a force
 * proportional to the product of their masses and inversely proportional
 * to the square of the distance between them.
 *
 * @param mass1 First mass (in kilograms, must be > 0)
 * @param mass2 Second mass (in kilograms, must be > 0)
 * @param distance Distance between centers (in meters, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Gravitational force (in Newtons)
 * @throws std::invalid_argument if parameters out of range
 */
inline double universalGravitationForce(double mass1, double mass2, double distance,
                                       double G = constants::G) {
    if (mass1 <= 0 || mass2 <= 0) {
        throw std::invalid_argument("Masses must be positive");
    }
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (G * mass1 * mass2) / (distance * distance);
}

/**
 * @brief Calculate gravitational field strength at distance r
 *
 * g = GM/r² (acceleration due to gravity)
 *
 * @param mass Source mass creating the field (in kilograms, must be > 0)
 * @param distance Distance from center of mass (in meters, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Gravitational field strength (in m/s²)
 * @throws std::invalid_argument if parameters out of range
 */
inline double gravitationalFieldStrength(double mass, double distance, double G = constants::G) {
    if (mass <= 0) {
        throw std::invalid_argument("Mass must be positive");
    }
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return (G * mass) / (distance * distance);
}

// ============================================================================
// Moon's Motions Connected with Fall of Apple
// ============================================================================

/**
 * @brief Calculate acceleration of falling apple at Earth's surface
 *
 * Newton's insight: g_surface = GM/R²
 * This is the acceleration we measure for falling objects on Earth
 *
 * @param earthMass Mass of Earth (default: 5.972e24 kg)
 * @param earthRadius Radius of Earth (default: 6.371e6 m)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Surface gravity (in m/s², should be ≈ 9.81)
 */
inline double calculateSurfaceGravity(double earthMass = constants::EARTH_MASS,
                                      double earthRadius = constants::EARTH_RADIUS,
                                      double G = constants::G) {
    return gravitationalFieldStrength(earthMass, earthRadius, G);
}

/**
 * @brief Calculate centripetal acceleration of Moon
 *
 * Moon's centripetal acceleration: a_c = v²/r = (2πr/T)²/r = 4π²r/T²
 *
 * Newton showed this equals GM/r², proving same force governs both
 * apple and Moon!
 *
 * @param orbitalRadius Moon's orbital radius (default: 3.844e8 m)
 * @param orbitalPeriod Moon's period (default: 2.36e6 s = 27.3 days)
 * @return Moon's centripetal acceleration (in m/s², ≈ 0.00272)
 */
inline double calculateMoonCentripetalAccel(double orbitalRadius = constants::MOON_ORBITAL_RADIUS,
                                            double orbitalPeriod = constants::MOON_ORBITAL_PERIOD) {
    return (4.0 * M_PI * M_PI * orbitalRadius) / (orbitalPeriod * orbitalPeriod);
}

/**
 * @brief Verify inverse square law using Moon and apple
 *
 * Newton's test: (g_surface / g_moon) should equal (r_moon / R_earth)²
 *
 * This ratio should be approximately 3600 (60²), showing the inverse
 * square relationship.
 *
 * @param earthMass Mass of Earth (default: 5.972e24 kg)
 * @param earthRadius Radius of Earth (default: 6.371e6 m)
 * @param moonDistance Moon's orbital radius (default: 3.844e8 m)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Ratio of surface to orbital gravity (should be ≈ 3600)
 */
inline double verifyInverseSquareLaw(double earthMass = constants::EARTH_MASS,
                                     double earthRadius = constants::EARTH_RADIUS,
                                     double moonDistance = constants::MOON_ORBITAL_RADIUS,
                                     double G = constants::G) {
    double surfaceG = calculateSurfaceGravity(earthMass, earthRadius, G);
    double moonG = gravitationalFieldStrength(earthMass, moonDistance, G);
    return surfaceG / moonG; // Should equal (moonDistance/earthRadius)²
}

/**
 * @brief Calculate ratio of distances (Moon orbit / Earth radius)
 *
 * r_moon / R_earth ≈ 60
 * Therefore (r_moon / R_earth)² ≈ 3600
 *
 * @param moonDistance Moon's orbital radius (default: 3.844e8 m)
 * @param earthRadius Earth's radius (default: 6.371e6 m)
 * @return Distance ratio
 */
inline double calculateDistanceRatio(double moonDistance = constants::MOON_ORBITAL_RADIUS,
                                     double earthRadius = constants::EARTH_RADIUS) {
    return moonDistance / earthRadius;
}

// ============================================================================
// Mass of the Earth
// ============================================================================

/**
 * @brief Calculate Earth's mass from surface gravity
 *
 * From g = GM/R², we get M = gR²/G
 *
 * This was first calculated by Cavendish in 1798 using a torsion balance
 * to measure G, thus "weighing the Earth".
 *
 * @param surfaceGravity Surface gravity (default: 9.81 m/s²)
 * @param earthRadius Earth's radius (in meters, default: 6.371e6 m)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Earth's mass (in kilograms, should be ≈ 5.972e24 kg)
 */
inline double calculateEarthMass(double surfaceGravity = 9.81,
                                 double earthRadius = constants::EARTH_RADIUS,
                                 double G = constants::G) {
    return (surfaceGravity * earthRadius * earthRadius) / G;
}

/**
 * @brief Calculate Earth's mass from Moon's orbit
 *
 * From Kepler's Third Law and gravitation:
 * T² = (4π²/GM)r³
 * Therefore: M = 4π²r³/(GT²)
 *
 * @param moonOrbitalRadius Moon's orbital radius (default: 3.844e8 m)
 * @param moonOrbitalPeriod Moon's period (default: 2.36e6 s)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Earth's mass (in kilograms)
 */
inline double calculateEarthMassFromMoon(double moonOrbitalRadius = constants::MOON_ORBITAL_RADIUS,
                                         double moonOrbitalPeriod = constants::MOON_ORBITAL_PERIOD,
                                         double G = constants::G) {
    double r3 = moonOrbitalRadius * moonOrbitalRadius * moonOrbitalRadius;
    double T2 = moonOrbitalPeriod * moonOrbitalPeriod;
    return (4.0 * M_PI * M_PI * r3) / (G * T2);
}

/**
 * @brief Calculate planet/star mass from satellite orbit
 *
 * General formula: M = 4π²r³/(GT²)
 *
 * @param satelliteRadius Orbital radius (in meters, must be > 0)
 * @param satellitePeriod Orbital period (in seconds, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Central body mass (in kilograms)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateCentralMassFromOrbit(double satelliteRadius, double satellitePeriod,
                                            double G = constants::G) {
    if (satelliteRadius <= 0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }
    if (satellitePeriod <= 0) {
        throw std::invalid_argument("Orbital period must be positive");
    }
    double r3 = satelliteRadius * satelliteRadius * satelliteRadius;
    double T2 = satellitePeriod * satellitePeriod;
    return (4.0 * M_PI * M_PI * r3) / (G * T2);
}

// ============================================================================
// Kepler's Third Law - Significance
// ============================================================================

/**
 * @brief Verify Kepler's Third Law
 *
 * Kepler's Third Law (empirical): T² ∝ r³
 * Newton's derivation: T² = (4π²/GM)r³
 *
 * Significance: The proportionality constant (4π²/GM) depends ONLY on
 * the central mass M, not on the orbiting mass. This means:
 * 1. All satellites around same body obey same T²/r³ ratio
 * 2. We can determine central mass from ANY satellite's orbit
 * 3. Law works for planets around Sun, moons around planets, etc.
 *
 * @param orbitalRadius Orbital radius (in meters, must be > 0)
 * @param centralMass Mass of central body (in kg, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Orbital period (in seconds)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculatePeriodKeplerThird(double orbitalRadius, double centralMass,
                                         double G = constants::G) {
    if (orbitalRadius <= 0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }
    if (centralMass <= 0) {
        throw std::invalid_argument("Central mass must be positive");
    }
    double r3 = orbitalRadius * orbitalRadius * orbitalRadius;
    return 2.0 * M_PI * std::sqrt(r3 / (G * centralMass));
}

/**
 * @brief Calculate T²/r³ ratio (Kepler's constant)
 *
 * K = T²/r³ = 4π²/(GM)
 *
 * This constant is the same for ALL satellites of the same central body.
 * Significance: Proves all satellites respond to same gravitational field.
 *
 * @param centralMass Mass of central body (in kg, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Kepler's constant (in s²/m³)
 * @throws std::invalid_argument if centralMass <= 0
 */
inline double calculateKeplerConstant(double centralMass, double G = constants::G) {
    if (centralMass <= 0) {
        throw std::invalid_argument("Central mass must be positive");
    }
    return (4.0 * M_PI * M_PI) / (G * centralMass);
}

/**
 * @brief Verify Kepler's Third Law for specific orbit
 *
 * Checks if T²/r³ equals expected value 4π²/(GM)
 *
 * @param orbitalRadius Orbital radius (in meters, must be > 0)
 * @param orbitalPeriod Orbital period (in seconds, must be > 0)
 * @param centralMass Mass of central body (in kg, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @param tolerance Acceptable error (default: 1e-6)
 * @return True if law is satisfied within tolerance
 */
inline bool verifyKeplerThirdLaw(double orbitalRadius, double orbitalPeriod,
                                 double centralMass, double G = constants::G,
                                 double tolerance = 1e-6) {
    double observedRatio = (orbitalPeriod * orbitalPeriod) /
                          (orbitalRadius * orbitalRadius * orbitalRadius);
    double expectedRatio = calculateKeplerConstant(centralMass, G);
    return std::abs(observedRatio - expectedRatio) / expectedRatio < tolerance;
}

/**
 * @brief Compare orbits of two satellites around same body
 *
 * For two satellites of same central body:
 * T₁²/T₂² = r₁³/r₂³
 *
 * This ratio is independent of central mass and G!
 *
 * @param radius1 First satellite's orbital radius (in meters, must be > 0)
 * @param radius2 Second satellite's orbital radius (in meters, must be > 0)
 * @return Ratio T₁/T₂ (periods ratio)
 * @throws std::invalid_argument if parameters out of range
 */
inline double compareSatellitePeriods(double radius1, double radius2) {
    if (radius1 <= 0 || radius2 <= 0) {
        throw std::invalid_argument("Radii must be positive");
    }
    // T₁/T₂ = (r₁/r₂)^(3/2)
    double radiusRatio = radius1 / radius2;
    return std::pow(radiusRatio, 1.5);
}

// ============================================================================
// Gravitational Potential Energy
// ============================================================================

/**
 * @brief Calculate gravitational potential energy
 *
 * U = -GMm/r (negative because we define U = 0 at infinity)
 *
 * @param mass1 First mass (in kg, must be > 0)
 * @param mass2 Second mass (in kg, must be > 0)
 * @param distance Separation distance (in meters, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Gravitational potential energy (in Joules, negative)
 * @throws std::invalid_argument if parameters out of range
 */
inline double gravitationalPotentialEnergy(double mass1, double mass2, double distance,
                                          double G = constants::G) {
    if (mass1 <= 0 || mass2 <= 0) {
        throw std::invalid_argument("Masses must be positive");
    }
    if (distance <= 0) {
        throw std::invalid_argument("Distance must be positive");
    }
    return -(G * mass1 * mass2) / distance;
}

/**
 * @brief Calculate escape velocity from celestial body
 *
 * v_escape = √(2GM/R)
 *
 * Significance: Speed needed to escape to infinity (where U = 0)
 *
 * @param centralMass Mass of body (in kg, must be > 0)
 * @param radius Distance from center (in meters, must be > 0)
 * @param G Gravitational constant (default: 6.674e-11 N⋅m²/kg²)
 * @return Escape velocity (in m/s)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateEscapeVelocity(double centralMass, double radius, double G = constants::G) {
    if (centralMass <= 0) {
        throw std::invalid_argument("Central mass must be positive");
    }
    if (radius <= 0) {
        throw std::invalid_argument("Radius must be positive");
    }
    return std::sqrt(2.0 * G * centralMass / radius);
}

// ============================================================================
// Binary Systems
// ============================================================================

/**
 * @brief Calculate reduced mass for binary system
 *
 * μ = m₁m₂/(m₁ + m₂)
 *
 * Used in two-body problem where both objects orbit their common center of mass
 *
 * @param mass1 First mass (in kg, must be > 0)
 * @param mass2 Second mass (in kg, must be > 0)
 * @return Reduced mass (in kg)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateReducedMass(double mass1, double mass2) {
    if (mass1 <= 0 || mass2 <= 0) {
        throw std::invalid_argument("Masses must be positive");
    }
    return (mass1 * mass2) / (mass1 + mass2);
}

/**
 * @brief Calculate center of mass distance for binary system
 *
 * For mass m₁ at distance r₁ and mass m₂ at distance r₂:
 * r₁ = (m₂/(m₁ + m₂)) × d
 * where d is the separation
 *
 * @param mass1 First mass (in kg, must be > 0)
 * @param mass2 Second mass (in kg, must be > 0)
 * @param separation Total separation (in meters, must be > 0)
 * @return Distance of mass1 from center of mass (in meters)
 * @throws std::invalid_argument if parameters out of range
 */
inline double calculateCMDistance(double mass1, double mass2, double separation) {
    if (mass1 <= 0 || mass2 <= 0) {
        throw std::invalid_argument("Masses must be positive");
    }
    if (separation <= 0) {
        throw std::invalid_argument("Separation must be positive");
    }
    return (mass2 / (mass1 + mass2)) * separation;
}

} // namespace gravitation
} // namespace physics

#endif // PHYSICS_GRAVITATION_HPP
