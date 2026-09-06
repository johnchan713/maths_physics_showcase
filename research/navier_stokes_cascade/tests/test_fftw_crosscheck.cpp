#include "ns_cascade/fftw_reference.hpp"
#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void expectRelativeNear(double actual,
                        double expected,
                        double tolerance,
                        const std::string& message) {
    const double scale = std::max(1.0, std::abs(expected));
    if (std::abs(actual - expected) > tolerance * scale) {
        throw std::runtime_error(message);
    }
}

double maximumDifference(const std::vector<ns_cascade::Complex>& left,
                         const std::vector<ns_cascade::Complex>& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("Cannot compare FFT arrays of different sizes");
    }
    double difference = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        difference = std::max(difference, std::abs(left[i] - right[i]));
    }
    return difference;
}

double relativeStateDifference(
    const ns_cascade::PseudospectralSystem::State& left,
    const ns_cascade::PseudospectralSystem::State& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("Cannot compare states of different sizes");
    }
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        difference_squared += ns_cascade::normSquared(left[i] - right[i]);
        reference_squared += ns_cascade::normSquared(right[i]);
    }
    return std::sqrt(difference_squared /
                     std::max(reference_squared, 1e-300));
}

ns_cascade::PseudospectralSystem::State fullSpectrumState(
    const ns_cascade::PseudospectralSystem& system) {
    ns_cascade::PseudospectralSystem::State state = system.zeroState();
    const std::vector<ns_cascade::WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < modes.size(); ++i) {
        const ns_cascade::WaveVector& wave = modes[i];
        const int maximum_component = std::max(
            std::abs(wave.x), std::max(std::abs(wave.y), std::abs(wave.z)));
        if (wave.normSquared() == 0 || maximum_component > system.cutoff()) {
            continue;
        }
        const double decay = std::exp(-0.12 * wave.normSquared());
        const double phase_x =
            0.19 * wave.x + 0.31 * wave.y - 0.23 * wave.z;
        const double phase_y =
            -0.29 * wave.x + 0.17 * wave.y + 0.37 * wave.z;
        const double phase_z =
            0.41 * wave.x - 0.13 * wave.y + 0.11 * wave.z;
        const ns_cascade::ComplexVector raw(
            decay * ns_cascade::Complex(std::cos(phase_x), std::sin(phase_x)),
            decay * ns_cascade::Complex(std::cos(phase_y), std::sin(phase_y)),
            decay * ns_cascade::Complex(std::cos(phase_z), std::sin(phase_z)));
        state[i] = ns_cascade::lerayProject(wave, raw);
    }
    const double scale = std::sqrt(10.0 / system.energy(state));
    for (std::size_t i = 0; i < state.size(); ++i) {
        state[i] = state[i] * scale;
    }
    return state;
}

void compareDiagnostics(
    const ns_cascade::PseudospectralSystem::Diagnostics& internal,
    const ns_cascade::FftwReferenceDiagnostics& reference,
    double tolerance) {
    expectRelativeNear(internal.energy,
                       reference.energy,
                       tolerance,
                       "Independent energy diagnostic differs");
    expectRelativeNear(internal.enstrophy,
                       reference.enstrophy,
                       tolerance,
                       "Independent enstrophy diagnostic differs");
    expectRelativeNear(internal.critical_h_half,
                       reference.critical_h_half,
                       tolerance,
                       "Independent H1/2 diagnostic differs");
    expectRelativeNear(internal.critical_l3_sample,
                       reference.critical_l3_sample,
                       tolerance,
                       "Independent sampled L3 diagnostic differs");
    expectRelativeNear(internal.sampled_vorticity_max,
                       reference.sampled_vorticity_max,
                       tolerance,
                       "Independent sampled vorticity diagnostic differs");
    expectRelativeNear(internal.high_shell_energy_fraction,
                       reference.cutoff_shell_energy_fraction,
                       tolerance,
                       "Independent cutoff-energy diagnostic differs");
    expectRelativeNear(internal.divergence_defect,
                       reference.divergence_defect,
                       tolerance,
                       "Independent divergence diagnostic differs");
    expectRelativeNear(internal.reality_defect,
                       reference.reality_defect,
                       tolerance,
                       "Independent reality diagnostic differs");
}

void testReferenceGridAndValidation() {
    const ns_cascade::PseudospectralSystem internal(16, 0.02, 5);
    const ns_cascade::FftwReferenceSystem reference(16, 0.02, 5);
    expect(reference.gridSize() == internal.gridSize() &&
               reference.cutoff() == internal.cutoff() &&
               reference.gridPointCount() == internal.gridPointCount(),
           "Independent FFTW grid metadata differs");
    for (std::size_t i = 0; i < internal.gridPointCount(); ++i) {
        const ns_cascade::WaveVector& left = internal.gridModes()[i];
        const ns_cascade::WaveVector& right = reference.gridModes()[i];
        expect(left.x == right.x && left.y == right.y && left.z == right.z,
               "Independent FFTW mode ordering differs");
    }

    bool invalid_cutoff_rejected = false;
    try {
        const ns_cascade::FftwReferenceSystem invalid(16, 0.02, 6);
        (void)invalid;
    } catch (const std::invalid_argument&) {
        invalid_cutoff_rejected = true;
    }
    expect(invalid_cutoff_rejected,
           "Independent FFTW path accepted an aliased cutoff");

    bool invalid_state_rejected = false;
    try {
        reference.rightHandSide(ns_cascade::FftwReferenceSystem::State(3));
    } catch (const std::invalid_argument&) {
        invalid_state_rejected = true;
    }
    expect(invalid_state_rejected,
           "Independent FFTW path accepted a mismatched state");
}

void testTransformsAgainstFftw() {
    const ns_cascade::PseudospectralSystem system(16, 0.05, 5);
    const ns_cascade::FftwReferenceSystem reference(16, 0.05, 5);
    std::vector<ns_cascade::Complex> physical(system.gridPointCount());
    for (std::size_t i = 0; i < physical.size(); ++i) {
        const double coordinate = static_cast<double>(i);
        physical[i] = ns_cascade::Complex(
            std::sin(0.013 * coordinate) + 0.2 * std::cos(0.071 * coordinate),
            0.3 * std::sin(0.037 * coordinate));
    }
    const std::vector<ns_cascade::Complex> internal_forward =
        system.forwardTransform(physical);
    const std::vector<ns_cascade::Complex> fftw_forward =
        reference.forwardTransform(physical);
    expect(maximumDifference(internal_forward, fftw_forward) < 2e-14,
           "In-repository forward FFT differs from FFTW");

    const std::vector<ns_cascade::Complex> internal_inverse =
        system.inverseTransform(internal_forward);
    const std::vector<ns_cascade::Complex> fftw_inverse =
        reference.inverseTransform(internal_forward);
    expect(maximumDifference(internal_inverse, fftw_inverse) < 2e-13,
           "In-repository inverse FFT differs from FFTW");
}

void testNavierStokesRightHandSideAgainstFftw() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const ns_cascade::FftwReferenceSystem reference_system(16, 0.02, 5);
    const ns_cascade::PseudospectralSystem::State state =
        system.vortexTubePairState(
            ns_cascade::VortexTubeParameters(0.7, 1.2, 0.3, 2), 10.0);
    const ns_cascade::PseudospectralSystem::State internal =
        system.rightHandSide(state);
    const ns_cascade::PseudospectralSystem::State reference =
        reference_system.rightHandSide(state);
    double difference = 0.0;
    for (std::size_t i = 0; i < internal.size(); ++i) {
        difference = std::max(
            difference, ns_cascade::norm(internal[i] - reference[i]));
    }
    expect(difference < 2e-12,
           "Navier-Stokes right-hand side differs from independent FFTW path");
}

void testFullSpectrumRightHandSideAgainstFftw() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const ns_cascade::FftwReferenceSystem reference_system(16, 0.02, 5);
    const ns_cascade::PseudospectralSystem::State state =
        fullSpectrumState(system);
    expect(system.divergenceDefect(state) < 2e-14,
           "Full-spectrum comparison state is not divergence-free");
    expect(system.realityDefect(state) < 2e-14,
           "Full-spectrum comparison state is not real-valued");
    const ns_cascade::PseudospectralSystem::State internal =
        system.rightHandSide(state);
    const ns_cascade::PseudospectralSystem::State reference =
        reference_system.rightHandSide(state);
    double difference = 0.0;
    for (std::size_t i = 0; i < internal.size(); ++i) {
        difference = std::max(
            difference, ns_cascade::norm(internal[i] - reference[i]));
    }
    expect(difference < 5e-12,
           "Full-spectrum right-hand side differs from the FFTW path");
}

void testIndependentDiagnostics() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const ns_cascade::FftwReferenceSystem reference_system(16, 0.02, 5);
    const ns_cascade::PseudospectralSystem::State state =
        system.vortexTubePairState(
            ns_cascade::VortexTubeParameters(0.7, 1.2, 0.3, 2), 10.0);
    compareDiagnostics(system.diagnostics(state),
                       reference_system.diagnostics(state),
                       3e-13);

    const ns_cascade::AdaptiveStepInfo internal_step =
        system.chooseAdaptiveTimeStep(state, 0.01, 0.35, 2.0);
    const ns_cascade::FftwReferenceStepInfo reference_step =
        reference_system.chooseAdaptiveTimeStep(state, 0.01, 0.35, 2.0);
    expectRelativeNear(internal_step.time_step,
                       reference_step.time_step,
                       2e-15,
                       "Independent adaptive time step differs");
    expectRelativeNear(internal_step.advective_cfl_upper_bound,
                       reference_step.advective_cfl_upper_bound,
                       2e-15,
                       "Independent CFL diagnostic differs");
}

void testFullIndependentTrajectory() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const ns_cascade::FftwReferenceSystem reference_system(16, 0.02, 5);
    ns_cascade::PseudospectralSystem::State internal_state =
        system.vortexTubePairState(
            ns_cascade::VortexTubeParameters(0.7, 1.2, 0.3, 2), 10.0);
    ns_cascade::FftwReferenceSystem::State reference_state = internal_state;

    double peak_relative_state_difference = 0.0;
    const double time_step = 0.0005;
    for (int step = 0; step < 40; ++step) {
        system.stepRungeKutta4(internal_state, time_step);
        reference_system.stepRungeKutta4(reference_state, time_step);
        peak_relative_state_difference = std::max(
            peak_relative_state_difference,
            relativeStateDifference(internal_state, reference_state));
        if ((step + 1) % 10 == 0) {
            compareDiagnostics(system.diagnostics(internal_state),
                               reference_system.diagnostics(reference_state),
                               2e-10);
        }
    }

    expect(peak_relative_state_difference < 2e-11,
           "Independent FFTW RK4 trajectory departed from the main solver");
    const ns_cascade::FftwReferenceDiagnostics final_diagnostics =
        reference_system.diagnostics(reference_state);
    expect(final_diagnostics.divergence_defect < 2e-12,
           "Independent FFTW trajectory developed a divergence defect");
    expect(final_diagnostics.reality_defect < 2e-12,
           "Independent FFTW trajectory lost Fourier reality symmetry");
}

}  // namespace

int main() {
    try {
        testReferenceGridAndValidation();
        testTransformsAgainstFftw();
        testNavierStokesRightHandSideAgainstFftw();
        testFullSpectrumRightHandSideAgainstFftw();
        testIndependentDiagnostics();
        testFullIndependentTrajectory();
        std::cout << "Independent FFTW trajectory cross-check passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
