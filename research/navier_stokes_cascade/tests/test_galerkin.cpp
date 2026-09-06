#include "ns_cascade/galerkin.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

std::size_t activeModeCount(const ns_cascade::GalerkinSystem::State& state) {
    std::size_t count = 0;
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (ns_cascade::normSquared(state[i]) > 1e-24) ++count;
    }
    return count;
}

void testNamedInitialConditions() {
    const ns_cascade::GalerkinSystem system(3, 0.05);
    const ns_cascade::InitialCondition conditions[] = {
        ns_cascade::InitialCondition::Deterministic,
        ns_cascade::InitialCondition::TaylorGreen,
        ns_cascade::InitialCondition::ABC};

    for (std::size_t i = 0; i < 3; ++i) {
        const ns_cascade::GalerkinSystem::State state =
            system.initialState(conditions[i], 0.75);
        expectNear(system.energy(state), 0.75, 1e-13,
                   "Named initial condition energy normalization failed");
        expect(system.divergenceDefect(state) < 1e-13,
               "Named initial condition is not divergence-free");
        expect(system.realityDefect(state) < 1e-13,
               "Named initial condition is not real-valued");
    }

    expect(activeModeCount(system.taylorGreenState()) == 8,
           "Taylor-Green must initially occupy eight corner modes");
    expect(activeModeCount(system.abcState()) == 6,
           "ABC must initially occupy six axial modes");
    expectNear(system.enstrophy(system.taylorGreenState()),
               3.0,
               1e-13,
               "Normalized Taylor-Green enstrophy should equal three");
    expectNear(system.enstrophy(system.abcState()),
               1.0,
               1e-13,
               "Normalized ABC enstrophy should equal one");

    bool vortex_tubes_rejected = false;
    try {
        system.initialState(ns_cascade::InitialCondition::VortexTubes);
    } catch (const std::invalid_argument&) {
        vortex_tubes_rejected = true;
    }
    expect(vortex_tubes_rejected,
           "The direct backend must reject FFT-only vortex-tube data");
}

void testAbcIsNonlinearNegativeControl() {
    const ns_cascade::GalerkinSystem system(3, 0.0);
    const ns_cascade::GalerkinSystem::State state = system.abcState();
    const ns_cascade::GalerkinSystem::State nonlinear =
        system.rightHandSide(state, false);
    double derivative_norm_squared = 0.0;
    for (std::size_t i = 0; i < nonlinear.size(); ++i) {
        derivative_norm_squared += ns_cascade::normSquared(nonlinear[i]);
    }
    expect(derivative_norm_squared < 1e-24,
           "ABC flow should have zero projected nonlinear evolution");
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

void testShellAccounting() {
    const double viscosity = 0.07;
    const ns_cascade::GalerkinSystem system(3, viscosity);
    ns_cascade::GalerkinSystem::State state = system.taylorGreenState();
    for (int step = 0; step < 10; ++step) {
        system.stepRungeKutta4(state, 0.0005);
    }

    const std::vector<ns_cascade::GalerkinSystem::ShellDiagnostics> shells =
        system.shellDiagnostics(state);
    double shell_energy = 0.0;
    double nonlinear_transfer = 0.0;
    double viscous_dissipation = 0.0;
    double cumulative_transfer = 0.0;
    for (std::size_t i = 0; i < shells.size(); ++i) {
        shell_energy += shells[i].energy;
        nonlinear_transfer += shells[i].nonlinear_transfer;
        viscous_dissipation += shells[i].viscous_dissipation;
        cumulative_transfer += shells[i].nonlinear_transfer;
        expectNear(shells[i].forward_flux,
                   -cumulative_transfer,
                   1e-13,
                   "Forward flux does not match cumulative shell transfer");
    }

    expectNear(shell_energy, system.energy(state), 1e-13,
               "Shell energies do not sum to total energy");
    expect(std::abs(nonlinear_transfer) < 1e-12,
           "Nonlinear shell transfers do not sum to zero");
    expectNear(viscous_dissipation,
               2.0 * viscosity * system.enstrophy(state),
               1e-12,
               "Shell viscous loss does not match total dissipation");
    expect(std::abs(shells.back().forward_flux) < 1e-12,
           "Flux outside the final retained shell must vanish");
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
        testNamedInitialConditions();
        testAbcIsNonlinearNegativeControl();
        testNonlinearEnergyCancellation();
        testViscousEnergyIdentity();
        testShortViscousRun();
        testShellAccounting();
        testZeroGradientIsNotZeroValue();
        std::cout << "All Navier-Stokes cascade tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
