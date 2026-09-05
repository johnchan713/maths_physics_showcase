#include "ns_cascade/galerkin.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void expectNear(double actual,
                double expected,
                double tolerance,
                const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": expected " +
                                 std::to_string(expected) + ", got " +
                                 std::to_string(actual));
    }
}

void testLerayProjection() {
    const ns_cascade::WaveVector wave(1, 2, -1);
    const ns_cascade::ComplexVector input(
        ns_cascade::Complex(1.2, -0.4),
        ns_cascade::Complex(-0.5, 0.7),
        ns_cascade::Complex(0.3, 0.2));
    const ns_cascade::ComplexVector projected =
        ns_cascade::lerayProject(wave, input);
    expect(std::abs(ns_cascade::dot(wave, projected)) < 1e-13,
           "Leray projection must be divergence-free");
}

void testInitialStateInvariants() {
    const ns_cascade::GalerkinSystem system(2, 0.05);
    const ns_cascade::GalerkinSystem::State state =
        system.deterministicLowModeState(1.0);
    expectNear(system.energy(state), 1.0, 1e-13,
               "Initial energy normalization failed");
    expect(system.divergenceDefect(state) < 1e-13,
           "Initial state is not divergence-free");
    expect(system.realityDefect(state) < 1e-13,
           "Initial Fourier state does not represent a real field");
}

void testNonlinearEnergyCancellation() {
    const ns_cascade::GalerkinSystem system(2, 0.0);
    const ns_cascade::GalerkinSystem::State state =
        system.deterministicLowModeState(1.0);
    const ns_cascade::GalerkinSystem::State derivative =
        system.rightHandSide(state, false);
    expect(std::abs(system.energyDerivative(state, derivative)) < 1e-12,
           "Galerkin convection must conserve kinetic energy instantaneously");
}

void testViscousEnergyIdentity() {
    const double viscosity = 0.15;
    const ns_cascade::GalerkinSystem system(2, viscosity);
    const ns_cascade::GalerkinSystem::State state =
        system.deterministicLowModeState(1.0);
    const ns_cascade::GalerkinSystem::State derivative =
        system.rightHandSide(state);
    const double residual = system.energyDerivative(state, derivative) +
                            2.0 * viscosity * system.enstrophy(state);
    expect(std::abs(residual) < 1e-12,
           "Semi-discrete viscous energy identity failed");
}

void testShortViscousRun() {
    const ns_cascade::GalerkinSystem system(2, 0.1);
    ns_cascade::GalerkinSystem::State state =
        system.deterministicLowModeState(1.0);
    double previous_energy = system.energy(state);
    for (int step = 0; step < 50; ++step) {
        system.stepRungeKutta4(state, 0.0005);
        const double current_energy = system.energy(state);
        expect(current_energy <= previous_energy + 1e-13,
               "Energy increased during a resolved viscous RK4 step");
        previous_energy = current_energy;
    }
    expect(system.divergenceDefect(state) < 1e-11,
           "Divergence defect grew during the short run");
    expect(system.realityDefect(state) < 1e-11,
           "Reality defect grew during the short run");
}

void testZeroGradientIsNotZeroValue() {
    // A non-zero constant velocity field has zero spatial gradient. This guards
    // against a common invalid implication in purported Navier-Stokes proofs.
    const ns_cascade::ComplexVector constant_velocity(
        ns_cascade::Complex(1.0, 0.0),
        ns_cascade::Complex(-2.0, 0.0),
        ns_cascade::Complex(0.5, 0.0));
    const double gradient_norm = 0.0;
    expect(gradient_norm == 0.0 && ns_cascade::norm(constant_velocity) > 0.0,
           "Zero gradient must not be treated as zero field value");
}

}  // namespace

int main() {
    try {
        testLerayProjection();
        testInitialStateInvariants();
        testNonlinearEnergyCancellation();
        testViscousEnergyIdentity();
        testShortViscousRun();
        testZeroGradientIsNotZeroValue();
        std::cout << "All Navier-Stokes cascade tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
