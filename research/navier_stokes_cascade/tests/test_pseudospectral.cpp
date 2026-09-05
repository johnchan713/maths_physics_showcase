#include "ns_cascade/galerkin.hpp"
#include "ns_cascade/pseudospectral.hpp"

#include <cmath>
#include <complex>
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

double compactDifference(
    const ns_cascade::GalerkinSystem& compact_system,
    const ns_cascade::GalerkinSystem::State& compact_state,
    const ns_cascade::PseudospectralSystem& fft_system,
    const ns_cascade::PseudospectralSystem::State& fft_state) {
    double difference = 0.0;
    for (std::size_t i = 0; i < compact_system.modes().size(); ++i) {
        const std::size_t fft_index = fft_system.indexOf(compact_system.modes()[i]);
        difference = std::max(
            difference,
            ns_cascade::norm(compact_state[i] - fft_state[fft_index]));
    }
    return difference;
}

ns_cascade::GalerkinSystem::State fullSpectrumState(
    const ns_cascade::GalerkinSystem& system) {
    ns_cascade::GalerkinSystem::State state = system.zeroState();
    for (std::size_t i = 0; i < system.modes().size(); ++i) {
        const ns_cascade::WaveVector& wave = system.modes()[i];
        const double decay = std::exp(-0.12 * wave.normSquared());
        const double phase_x = 0.19 * wave.x + 0.31 * wave.y - 0.23 * wave.z;
        const double phase_y = -0.29 * wave.x + 0.17 * wave.y + 0.37 * wave.z;
        const double phase_z = 0.41 * wave.x - 0.13 * wave.y + 0.11 * wave.z;
        const ns_cascade::ComplexVector raw(
            decay * ns_cascade::Complex(std::cos(phase_x), std::sin(phase_x)),
            decay * ns_cascade::Complex(std::cos(phase_y), std::sin(phase_y)),
            decay * ns_cascade::Complex(std::cos(phase_z), std::sin(phase_z)));
        state[i] = ns_cascade::lerayProject(wave, raw);
    }
    const double scale = std::sqrt(1.0 / system.energy(state));
    for (std::size_t i = 0; i < state.size(); ++i) state[i] = state[i] * scale;
    return state;
}

ns_cascade::PseudospectralSystem::State expandCompactState(
    const ns_cascade::GalerkinSystem& compact_system,
    const ns_cascade::GalerkinSystem::State& compact_state,
    const ns_cascade::PseudospectralSystem& fft_system) {
    ns_cascade::PseudospectralSystem::State fft_state = fft_system.zeroState();
    for (std::size_t i = 0; i < compact_system.modes().size(); ++i) {
        fft_state[fft_system.indexOf(compact_system.modes()[i])] = compact_state[i];
    }
    return fft_state;
}

void testGridAndDealiasingRules() {
    const ns_cascade::PseudospectralSystem system(8, 0.05);
    expect(system.gridSize() == 8, "Unexpected FFT grid size");
    expect(system.safeDealiasCutoff() == 2, "Incorrect strict 2/3 cutoff");
    expect(system.cutoff() == 2, "Default cutoff must use the safe limit");
    expect(system.modeCount() == 124, "K=2 must retain 124 non-zero modes");

    bool rejected = false;
    try {
        const ns_cascade::PseudospectralSystem invalid(8, 0.05, 3);
        (void)invalid;
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    expect(rejected, "An aliased cutoff must be rejected");
}

void testScalarFftRoundTrip() {
    const ns_cascade::PseudospectralSystem system(8, 0.05);
    std::vector<ns_cascade::Complex> physical(system.gridPointCount());
    for (std::size_t i = 0; i < physical.size(); ++i) {
        const double coordinate = static_cast<double>(i);
        physical[i] = ns_cascade::Complex(
            std::sin(0.17 * coordinate) + 0.3 * std::cos(0.07 * coordinate),
            0.2 * std::sin(0.11 * coordinate));
    }
    const std::vector<ns_cascade::Complex> recovered =
        system.inverseTransform(system.forwardTransform(physical));
    double maximum_error = 0.0;
    for (std::size_t i = 0; i < physical.size(); ++i) {
        maximum_error = std::max(maximum_error,
                                 std::abs(physical[i] - recovered[i]));
    }
    expect(maximum_error < 2e-15, "3D FFT round trip lost accuracy");
}

void testInitialDataMatchesCompactOracle() {
    const ns_cascade::GalerkinSystem compact(2, 0.05);
    const ns_cascade::PseudospectralSystem fft(8, 0.05, 2);
    const ns_cascade::InitialCondition conditions[] = {
        ns_cascade::InitialCondition::Deterministic,
        ns_cascade::InitialCondition::TaylorGreen,
        ns_cascade::InitialCondition::ABC};
    for (std::size_t i = 0; i < 3; ++i) {
        const ns_cascade::GalerkinSystem::State compact_state =
            compact.initialState(conditions[i]);
        const ns_cascade::PseudospectralSystem::State fft_state =
            fft.initialState(conditions[i]);
        expect(compactDifference(compact, compact_state, fft, fft_state) < 1e-15,
               "FFT and compact initial coefficients differ");
    }
}

void testNonlinearTermMatchesCompactOracle() {
    const double viscosity = 0.05;
    const ns_cascade::GalerkinSystem compact(2, viscosity);
    const ns_cascade::PseudospectralSystem fft(8, viscosity, 2);
    const ns_cascade::GalerkinSystem::State compact_state =
        compact.initialState(ns_cascade::InitialCondition::Deterministic);
    const ns_cascade::PseudospectralSystem::State fft_state =
        fft.initialState(ns_cascade::InitialCondition::Deterministic);

    const ns_cascade::GalerkinSystem::State compact_nonlinear =
        compact.rightHandSide(compact_state, false);
    const ns_cascade::PseudospectralSystem::State fft_nonlinear =
        fft.rightHandSide(fft_state, false);
    expect(compactDifference(compact,
                             compact_nonlinear,
                             fft,
                             fft_nonlinear) < 2e-14,
           "Dealiased FFT nonlinear term differs from direct convolution");

    const ns_cascade::GalerkinSystem::State compact_full =
        compact.rightHandSide(compact_state);
    const ns_cascade::PseudospectralSystem::State fft_full =
        fft.rightHandSide(fft_state);
    expect(compactDifference(compact, compact_full, fft, fft_full) < 2e-14,
           "FFT and direct full right-hand sides differ");
}

void testMaximumCutoffMatchesOracleForFullSpectrum() {
    const ns_cascade::GalerkinSystem compact(5, 0.0);
    const ns_cascade::PseudospectralSystem fft(16, 0.0, 5);
    const ns_cascade::GalerkinSystem::State compact_state =
        fullSpectrumState(compact);
    const ns_cascade::PseudospectralSystem::State fft_state =
        expandCompactState(compact, compact_state, fft);
    expect(compact.divergenceDefect(compact_state) < 2e-15,
           "Full-spectrum oracle state is not divergence-free");
    expect(compact.realityDefect(compact_state) < 2e-15,
           "Full-spectrum oracle state is not real-valued");

    const ns_cascade::GalerkinSystem::State compact_nonlinear =
        compact.rightHandSide(compact_state, false);
    const ns_cascade::PseudospectralSystem::State fft_nonlinear =
        fft.rightHandSide(fft_state, false);
    expect(compactDifference(compact,
                             compact_nonlinear,
                             fft,
                             fft_nonlinear) < 5e-13,
           "Maximum-cutoff FFT disagrees with full-spectrum direct convolution");
}

void testEnergyAndShellIdentities() {
    const double viscosity = 0.07;
    const ns_cascade::PseudospectralSystem system(16, viscosity, 5);
    const ns_cascade::PseudospectralSystem::State state =
        system.initialState(ns_cascade::InitialCondition::TaylorGreen);
    const ns_cascade::PseudospectralSystem::State nonlinear =
        system.rightHandSide(state, false);
    expect(std::abs(system.energyDerivative(state, nonlinear)) < 2e-13,
           "Dealiased FFT nonlinearity must conserve energy");

    const std::vector<ns_cascade::PseudospectralSystem::ShellDiagnostics> shells =
        system.shellDiagnostics(state);
    double energy_sum = 0.0;
    double transfer_sum = 0.0;
    double dissipation_sum = 0.0;
    for (std::size_t i = 0; i < shells.size(); ++i) {
        energy_sum += shells[i].energy;
        transfer_sum += shells[i].nonlinear_transfer;
        dissipation_sum += shells[i].viscous_dissipation;
    }
    expectNear(energy_sum, system.energy(state), 2e-13,
               "FFT shell energies do not sum to total energy");
    expect(std::abs(transfer_sum) < 2e-13,
           "FFT nonlinear shell transfers do not sum to zero");
    expectNear(dissipation_sum,
               2.0 * viscosity * system.enstrophy(state),
               2e-13,
               "FFT shell viscous losses do not sum correctly");
    expect(std::abs(shells.back().forward_flux) < 2e-13,
           "FFT flux outside the retained set must vanish");
}

void testShortTrajectoryMatchesCompactOracle() {
    const ns_cascade::GalerkinSystem compact(2, 0.05);
    const ns_cascade::PseudospectralSystem fft(8, 0.05, 2);
    ns_cascade::GalerkinSystem::State compact_state =
        compact.initialState(ns_cascade::InitialCondition::TaylorGreen);
    ns_cascade::PseudospectralSystem::State fft_state =
        fft.initialState(ns_cascade::InitialCondition::TaylorGreen);
    for (int step = 0; step < 10; ++step) {
        compact.stepRungeKutta4(compact_state, 0.0005);
        fft.stepRungeKutta4(fft_state, 0.0005);
    }
    expect(compactDifference(compact, compact_state, fft, fft_state) < 2e-13,
           "FFT trajectory departed from the direct-convolution oracle");
    expect(fft.divergenceDefect(fft_state) < 2e-13,
           "FFT trajectory developed a divergence defect");
    expect(fft.realityDefect(fft_state) < 2e-13,
           "FFT trajectory lost Fourier reality symmetry");
}

void testAbcNegativeControl() {
    const ns_cascade::PseudospectralSystem system(16, 0.0, 5);
    const ns_cascade::PseudospectralSystem::State state =
        system.initialState(ns_cascade::InitialCondition::ABC);
    const ns_cascade::PseudospectralSystem::State nonlinear =
        system.rightHandSide(state, false);
    double derivative_norm_squared = 0.0;
    for (std::size_t i = 0; i < nonlinear.size(); ++i) {
        derivative_norm_squared += ns_cascade::normSquared(nonlinear[i]);
    }
    expect(derivative_norm_squared < 1e-24,
           "FFT backend reported a false ABC cascade");
}

}  // namespace

int main() {
    try {
        testGridAndDealiasingRules();
        testScalarFftRoundTrip();
        testInitialDataMatchesCompactOracle();
        testNonlinearTermMatchesCompactOracle();
        testMaximumCutoffMatchesOracleForFullSpectrum();
        testEnergyAndShellIdentities();
        testShortTrajectoryMatchesCompactOracle();
        testAbcNegativeControl();
        std::cout << "All pseudospectral Navier-Stokes tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
