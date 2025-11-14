# Physics Showcase - C++ Implementation

A comprehensive C++ library implementing fundamental physics concepts with well-documented standalone functions. This showcase demonstrates classical mechanics principles through clean, reusable code organized in header-only libraries.

## Overview

This project implements core physics functionality as standalone C++ functions, categorized into three main modules:

1. **Newton's Laws of Motion** (`newton_laws.hpp`)
2. **Kinematics** - Motion with constant acceleration (`kinematics.hpp`)
3. **Dynamics** - Force causing rectilinear motion (`dynamics.hpp`)

All functions are thoroughly documented with parameter descriptions, units, return values, and exception handling.

## Project Structure

```
physics_showcase/
├── include/
│   └── physics/
│       ├── newton_laws.hpp      # Newton's three laws of motion
│       ├── kinematics.hpp       # Motion equations
│       └── dynamics.hpp         # Force and motion integration
├── examples/
│   └── main.cpp                 # Comprehensive usage examples
├── LICENSE
└── README.md
```

## Features

### Newton's Laws of Motion (`newton_laws.hpp`)

#### First Law - Law of Inertia
- `isInEquilibrium()` - Check if object is in equilibrium
- `calculateNetForce()` - Calculate net force from multiple forces

#### Second Law - F = ma
- `calculateForce()` - Calculate force from mass and acceleration
- `calculateAcceleration()` - Calculate acceleration from force and mass
- `calculateMass()` - Calculate mass from force and acceleration
- `calculateWeight()` - Calculate gravitational force

#### Third Law - Action-Reaction
- `calculateReactionForce()` - Get reaction force for given action
- `verifyActionReactionPair()` - Verify if two forces form valid pair

### Kinematics (`kinematics.hpp`)

Motion in a straight line with constant acceleration using all kinematic equations:

#### First Equation: v = v₀ + at
- `calculateFinalVelocity()` - Calculate final velocity
- `calculateAccelerationFromVelocities()` - Calculate acceleration
- `calculateTimeFromVelocities()` - Calculate time required

#### Second Equation: s = v₀t + ½at²
- `calculateDisplacement()` - Calculate displacement
- `calculateAccelerationFromDisplacement()` - Calculate acceleration

#### Third Equation: v² = v₀² + 2as
- `calculateFinalVelocityFromDisplacement()` - Calculate final velocity
- `calculateAccelerationFromVelocitySquared()` - Calculate acceleration
- `calculateDisplacementFromVelocities()` - Calculate displacement

#### Fourth Equation: s = ((v + v₀) / 2) × t
- `calculateDisplacementFromAverageVelocity()` - Using average velocity
- `calculateAverageVelocity()` - Calculate average velocity
- `calculateTimeFromAverageVelocity()` - Calculate time

#### Utility Functions
- `calculateDistance()` - Calculate distance traveled (always positive)
- `calculateStoppingDistance()` - Distance required to stop
- `calculateStoppingTime()` - Time required to stop

### Dynamics (`dynamics.hpp`)

Integration of forces with motion - how forces cause acceleration and resulting motion:

#### Force-Acceleration Relationships
- `calculateNetForce()` - Sum multiple forces
- `calculateAccelerationFromForce()` - F = ma → a = F/m
- `calculateRequiredForce()` - Force needed for desired acceleration

#### Force-Motion Integration: Velocity
- `calculateFinalVelocityFromForce()` - Final velocity after force acts
- `calculateVelocityChange()` - Change in velocity (Δv)
- `calculateTimeForVelocityChange()` - Time for velocity change

#### Force-Motion Integration: Displacement
- `calculateDisplacementFromForce()` - Displacement when force acts
- `calculateFinalVelocityFromForceAndDisplacement()` - Final velocity from force and displacement

#### Stopping Problems
- `calculateBrakingForce()` - Force needed to stop in given distance
- `calculateStoppingDistanceFromForce()` - Distance to stop with given force
- `calculateStoppingTimeFromForce()` - Time to stop with given force

#### Friction
- `calculateFrictionForce()` - Frictional force (f = μN)
- `calculateAccelerationWithFriction()` - Net acceleration with friction
- `calculateMinimumForceToOvercomeFriction()` - Minimum force to start motion

#### Work and Power
- `calculateWork()` - Work done by force (W = F·s)
- `calculatePower()` - Power (rate of work, P = F·v)

## Units

All functions use SI units:
- **Mass**: kilograms (kg)
- **Distance/Displacement**: meters (m)
- **Time**: seconds (s)
- **Velocity**: meters per second (m/s)
- **Acceleration**: meters per second squared (m/s²)
- **Force**: Newtons (N)
- **Work/Energy**: Joules (J)
- **Power**: Watts (W)

## Building and Running

### Prerequisites
- C++ compiler with C++11 support or later (g++, clang++, MSVC)
- Make (optional, for using Makefile)
- CMake 3.10+ (optional, for using CMake)

### Option 1: Direct Compilation

```bash
g++ -std=c++11 -I./include examples/main.cpp -o physics_demo
./physics_demo
```

### Option 2: Using Make

```bash
make
./physics_demo
```

### Option 3: Using CMake

```bash
mkdir build
cd build
cmake ..
make
./physics_demo
```

## Usage Examples

### Example 1: Newton's Second Law

```cpp
#include "physics/newton_laws.hpp"

// Calculate force required to accelerate 10kg object at 5 m/s²
double mass = 10.0;           // kg
double acceleration = 5.0;    // m/s²
double force = physics::newton::calculateForce(mass, acceleration);
// Result: 50.0 N
```

### Example 2: Kinematics - Free Fall

```cpp
#include "physics/kinematics.hpp"

// Object dropped from rest, falling for 3 seconds
double initialVelocity = 0.0;     // m/s
double gravity = 9.81;             // m/s²
double time = 3.0;                 // s

double finalVelocity = physics::kinematics::calculateFinalVelocity(
    initialVelocity, gravity, time);
// Result: 29.43 m/s

double distance = physics::kinematics::calculateDisplacement(
    initialVelocity, gravity, time);
// Result: 44.145 m
```

### Example 3: Dynamics - Braking Force

```cpp
#include "physics/dynamics.hpp"

// Calculate braking force to stop 1500kg car in 40m
double mass = 1500.0;           // kg
double speed = 25.0;            // m/s (90 km/h)
double distance = 40.0;         // m

double brakingForce = physics::dynamics::calculateBrakingForce(
    mass, speed, distance);
// Result: -11718.75 N (negative = opposite to motion)
```

### Example 4: Motion with Friction

```cpp
#include "physics/dynamics.hpp"

// Box pushed across floor with friction
double mass = 20.0;                  // kg
double appliedForce = 100.0;         // N
double frictionCoefficient = 0.25;   // dimensionless

double acceleration = physics::dynamics::calculateAccelerationWithFriction(
    mass, appliedForce, frictionCoefficient);
// Result: 2.5475 m/s²
```

## Error Handling

All functions validate input parameters and throw `std::invalid_argument` exceptions for invalid inputs:
- Mass must be greater than zero
- Time must be non-negative
- Certain denominators must be non-zero
- Physical constraints (e.g., v² cannot be negative in real solutions)

Example:
```cpp
try {
    double force = physics::newton::calculateForce(-10.0, 5.0);
} catch (const std::invalid_argument& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    // Output: "Error: Mass must be greater than zero"
}
```

## Function Documentation

All functions include detailed Doxygen-style comments:
- **@brief**: Brief description of what the function does
- **@param**: Description of each parameter with units
- **@return**: Description of return value with units
- **@throws**: Exceptions that may be thrown

Example:
```cpp
/**
 * @brief Calculate force using Newton's Second Law: F = ma
 *
 * @param mass Mass of the object (in kilograms, must be > 0)
 * @param acceleration Acceleration of the object (in m/s²)
 * @return Force acting on the object (in Newtons)
 * @throws std::invalid_argument if mass <= 0
 */
inline double calculateForce(double mass, double acceleration);
```

## Design Principles

1. **Header-Only Library**: All functions are inline in headers for easy integration
2. **Namespace Organization**: Functions organized in `physics::newton`, `physics::kinematics`, `physics::dynamics`
3. **Comprehensive Documentation**: Every function fully documented with units and constraints
4. **Error Handling**: Input validation with meaningful exception messages
5. **SI Units**: Consistent use of SI units throughout
6. **Pure Functions**: No side effects, easy to test and reason about

## Future Extensions

This showcase can be extended with additional physics modules:
- 2D and 3D motion (vectors)
- Projectile motion
- Circular motion and rotation
- Energy and momentum
- Oscillations and waves
- Gravitation
- Fluid mechanics
- Thermodynamics

## License

This project is licensed under the terms specified in the LICENSE file.

## Contributing

When adding new physics functions:
1. Follow the existing documentation style
2. Include parameter descriptions with units
3. Add input validation with appropriate exceptions
4. Use SI units consistently
5. Add usage examples to the main.cpp
6. Update this README with new functionality

## References

The physics equations implemented are based on standard classical mechanics:
- Newton's Laws of Motion
- Kinematic equations for constant acceleration
- Dynamics combining forces and motion
- SI unit system (International System of Units)
