#include "ns_cascade/galerkin.hpp"
#include "ns_cascade/checkpoint.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
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

void testCriticalAndEnstrophyDiagnosticsMatchCompactOracle() {
    const double viscosity = 0.07;
    const ns_cascade::GalerkinSystem compact(2, viscosity);
    const ns_cascade::PseudospectralSystem fft(8, viscosity, 2);
    const ns_cascade::GalerkinSystem::State compact_state =
        compact.deterministicLowModeState(1.0);
    const ns_cascade::PseudospectralSystem::State fft_state =
        expandCompactState(compact, compact_state, fft);
    const ns_cascade::GalerkinSystem::Diagnostics compact_diagnostics =
        compact.diagnostics(compact_state);
    const ns_cascade::PseudospectralSystem::Diagnostics fft_diagnostics =
        fft.diagnostics(fft_state);

    expectNear(fft_diagnostics.critical_h_half,
               compact_diagnostics.critical_h_half,
               2e-13,
               "FFT critical H1/2 norm differs from direct convolution");
    expectNear(fft_diagnostics.nonlinear_enstrophy_production,
               compact_diagnostics.nonlinear_enstrophy_production,
               2e-12,
               "FFT nonlinear enstrophy production differs from direct convolution");
    expectNear(fft_diagnostics.viscous_enstrophy_destruction,
               compact_diagnostics.viscous_enstrophy_destruction,
               2e-13,
               "FFT viscous enstrophy destruction differs from direct convolution");

    const ns_cascade::PseudospectralSystem::State full_derivative =
        fft.rightHandSide(fft_state);
    expectNear(fft.enstrophyDerivative(fft_state, full_derivative),
               fft_diagnostics.net_enstrophy_rate,
               2e-12,
               "FFT enstrophy budget does not match the full evolution");
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

double nonzeroAxialEnergyFraction(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::PseudospectralSystem::State& state) {
    double axial_energy = 0.0;
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (system.gridModes()[i].z != 0) {
            axial_energy += 0.5 * ns_cascade::normSquared(state[i]);
        }
    }
    return axial_energy / system.energy(state);
}

void testVortexTubeInitialData() {
    const ns_cascade::PseudospectralSystem system(32, 0.05, 10);
    const ns_cascade::VortexTubeParameters straight_parameters(
        0.55, 1.6, 0.0, 1);
    const ns_cascade::VortexTubeParameters bent_parameters(
        0.55, 1.6, 0.25, 2);
    const ns_cascade::PseudospectralSystem::State straight =
        system.vortexTubePairState(straight_parameters, 1.0);
    const ns_cascade::PseudospectralSystem::State bent =
        system.vortexTubePairState(bent_parameters, 1.0);

    expectNear(system.energy(straight), 1.0, 2e-13,
               "Straight vortex tubes were not energy-normalized");
    expectNear(system.energy(bent), 1.0, 2e-13,
               "Bent vortex tubes were not energy-normalized");
    expect(system.divergenceDefect(bent) < 2e-13,
           "Bent vortex tubes are not divergence-free");
    expect(system.realityDefect(bent) < 2e-13,
           "Bent vortex tubes do not represent a real velocity field");
    expect(nonzeroAxialEnergyFraction(system, straight) < 1e-24,
           "Straight tubes unexpectedly contain three-dimensional modes");
    expect(nonzeroAxialEnergyFraction(system, bent) > 1e-5,
           "Bent tubes failed to populate three-dimensional modes");
    expect(system.diagnostics(bent).high_shell_energy_fraction < 0.01,
           "Default bent tubes begin under-resolved");

    bool invalid_core_rejected = false;
    try {
        system.vortexTubePairState(
            ns_cascade::VortexTubeParameters(0.0, 1.6, 0.25, 1));
    } catch (const std::invalid_argument&) {
        invalid_core_rejected = true;
    }
    expect(invalid_core_rejected, "A zero tube core radius must be rejected");

    bool invalid_axial_mode_rejected = false;
    try {
        system.vortexTubePairState(
            ns_cascade::VortexTubeParameters(0.55, 1.6, 0.25, 11));
    } catch (const std::invalid_argument&) {
        invalid_axial_mode_rejected = true;
    }
    expect(invalid_axial_mode_rejected,
           "An unresolved tube axial mode must be rejected");
}

void testAdaptiveTimeStepControl() {
    const ns_cascade::PseudospectralSystem system(16, 0.2, 5);
    ns_cascade::PseudospectralSystem::State state =
        system.vortexTubePairState(ns_cascade::VortexTubeParameters());
    const double target_cfl = 0.08;
    const double diffusion_safety = 1.5;
    const ns_cascade::AdaptiveStepInfo information =
        system.chooseAdaptiveTimeStep(
            state, 0.1, target_cfl, diffusion_safety);
    expect(information.time_step > 0.0 && information.time_step < 0.1,
           "Adaptive control did not reduce an unsafe trial step");
    expect(information.advective_cfl_upper_bound <=
               target_cfl * (1.0 + 2e-15),
           "Adaptive step exceeds its conservative CFL target");
    expect(information.viscous_stability_number <=
               diffusion_safety * (1.0 + 2e-15),
           "Adaptive step exceeds its viscous stability bound");

    const ns_cascade::PseudospectralSystem diffusion_limited_system(16, 10.0, 5);
    const ns_cascade::PseudospectralSystem::State diffusion_state =
        diffusion_limited_system.initialState(
            ns_cascade::InitialCondition::TaylorGreen);
    const ns_cascade::AdaptiveStepInfo diffusion_limited =
        diffusion_limited_system.chooseAdaptiveTimeStep(
            diffusion_state, 0.1, 100.0, diffusion_safety);
    expectNear(diffusion_limited.viscous_stability_number,
               diffusion_safety,
               2e-15,
               "Viscous control did not become active when required");

    const ns_cascade::AdaptiveStepInfo accepted =
        system.stepAdaptiveRungeKutta4(
            state, 0.1, target_cfl, diffusion_safety);
    expectNear(accepted.time_step, information.time_step, 1e-16,
               "Adaptive stepping used a different proposed step");
    const ns_cascade::PseudospectralSystem::Diagnostics diagnostics =
        system.diagnostics(state);
    expect(std::isfinite(diagnostics.energy),
           "Adaptive step produced non-finite energy");
    expect(diagnostics.divergence_defect < 2e-13,
           "Adaptive step developed a divergence defect");
    expect(diagnostics.reality_defect < 2e-13,
           "Adaptive step lost Fourier reality symmetry");

    bool invalid_cfl_rejected = false;
    try {
        system.chooseAdaptiveTimeStep(state, 0.01, 0.0, 1.0);
    } catch (const std::invalid_argument&) {
        invalid_cfl_rejected = true;
    }
    expect(invalid_cfl_rejected, "A non-positive CFL target must be rejected");
}

void testRescaledSpectrumProfile() {
    const ns_cascade::PseudospectralSystem system(16, 0.05, 5);
    const ns_cascade::PseudospectralSystem::State unit_state =
        system.initialState(ns_cascade::InitialCondition::TaylorGreen, 1.0);
    const ns_cascade::PseudospectralSystem::State scaled_state =
        system.initialState(ns_cascade::InitialCondition::TaylorGreen, 9.0);
    const ns_cascade::SpectrumProfile unit_profile =
        ns_cascade::rescaledSpectrumProfile(system, unit_state, 16, 4.0);
    const ns_cascade::SpectrumProfile scaled_profile =
        ns_cascade::rescaledSpectrumProfile(system, scaled_state, 16, 4.0);
    expectNear(system.cutoffShellEnergyFraction(unit_state),
               system.diagnostics(unit_state).high_shell_energy_fraction,
               0.0,
               "Public cutoff-shell fraction disagrees with diagnostics");

    expectNear(unit_profile.characteristic_wavenumber,
               std::sqrt(3.0),
               2e-15,
               "Taylor-Green RMS wavenumber is incorrect");
    double fraction_sum = 0.0;
    for (std::size_t i = 0; i < unit_profile.energy_fractions.size(); ++i) {
        fraction_sum += unit_profile.energy_fractions[i];
    }
    expectNear(fraction_sum, 1.0, 2e-15,
               "Rescaled spectral energy fractions do not sum to one");
    expect(ns_cascade::spectrumProfileL1Distance(unit_profile, scaled_profile) <
               2e-15,
           "Rescaled spectrum changed under amplitude-only energy scaling");
    expectNear(ns_cascade::spectrumProfileOverlap(unit_profile, scaled_profile),
               1.0,
               2e-15,
               "Identical rescaled spectra do not have unit overlap");
    const ns_cascade::SpectrumProfileChange amplitude_change =
        ns_cascade::spectrumProfileChange(unit_profile, scaled_profile);
    expect(!amplitude_change.scale_normalized_drift_valid,
           "An unchanged spectral scale should not report a normalized drift");

    ns_cascade::SpectrumProfile exactly_rescaled_profile = unit_profile;
    exactly_rescaled_profile.characteristic_wavenumber *= std::exp(0.125);
    const ns_cascade::SpectrumProfileChange exactly_rescaled_change =
        ns_cascade::spectrumProfileChange(
            exactly_rescaled_profile, unit_profile);
    expect(exactly_rescaled_change.scale_normalized_drift_valid,
           "A finite spectral-scale movement did not produce a drift metric");
    expectNear(exactly_rescaled_change.absolute_log_scale_change,
               0.125,
               2e-15,
               "Logarithmic spectral-scale movement is incorrect");
    expectNear(exactly_rescaled_change.l1_per_log_scale_change,
               0.0,
               2e-15,
               "An exactly preserved rescaled shape has nonzero drift");

    ns_cascade::PseudospectralSystem::State low_shell = system.zeroState();
    ns_cascade::PseudospectralSystem::State high_shell = system.zeroState();
    for (std::size_t i = 0; i < system.gridModes().size(); ++i) {
        const ns_cascade::WaveVector& wave = system.gridModes()[i];
        if (wave.y == 0 && wave.z == 0 && std::abs(wave.x) == 1) {
            low_shell[i] = ns_cascade::ComplexVector(0.0, 1.0, 0.0);
        }
        if (wave.y == 0 && wave.z == 0 && std::abs(wave.x) == 2) {
            high_shell[i] = ns_cascade::ComplexVector(0.0, 1.0, 0.0);
        }
    }
    const ns_cascade::SpectrumProfile low_shell_profile =
        ns_cascade::rescaledSpectrumProfile(system, low_shell, 16, 4.0);
    const ns_cascade::SpectrumProfile high_shell_profile =
        ns_cascade::rescaledSpectrumProfile(system, high_shell, 16, 4.0);
    const ns_cascade::SpectrumProfileChange shell_rescaling =
        ns_cascade::spectrumProfileChange(
            high_shell_profile, low_shell_profile);
    expectNear(high_shell_profile.characteristic_wavenumber /
                   low_shell_profile.characteristic_wavenumber,
               2.0,
               0.0,
               "Synthetic spectrum did not move by its exact scale factor");
    expectNear(shell_rescaling.l1_per_log_scale_change,
               0.0,
               0.0,
               "An exact discrete spectral rescaling has nonzero shape drift");

    ns_cascade::PseudospectralSystem::State evolved = unit_state;
    for (int step = 0; step < 20; ++step) {
        system.stepRungeKutta4(evolved, 0.0005);
    }
    const ns_cascade::SpectrumProfile evolved_profile =
        ns_cascade::rescaledSpectrumProfile(system, evolved, 16, 4.0);
    const double distance =
        ns_cascade::spectrumProfileL1Distance(unit_profile, evolved_profile);
    expect(distance >= 0.0 && distance <= 2.0 + 2e-15,
           "Rescaled-profile L1 distance is outside its mathematical range");
}

ns_cascade::SimulationCheckpoint makeTestCheckpoint(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::PseudospectralSystem::State& state,
    std::uint64_t step,
    double time) {
    const ns_cascade::PseudospectralSystem::Diagnostics diagnostics =
        system.diagnostics(state);
    const ns_cascade::SpectrumProfile profile =
        ns_cascade::rescaledSpectrumProfile(system, state, 16, 4.0);
    ns_cascade::SimulationCheckpoint checkpoint;
    checkpoint.configuration.grid_size = system.gridSize();
    checkpoint.configuration.cutoff = system.cutoff();
    checkpoint.configuration.viscosity = system.viscosity();
    checkpoint.configuration.initial_condition =
        ns_cascade::InitialCondition::TaylorGreen;
    checkpoint.configuration.initial_energy = 1.0;
    checkpoint.configuration.adaptive = false;
    checkpoint.configuration.maximum_time_step = 0.0005;
    checkpoint.configuration.target_cfl = 0.4;
    checkpoint.configuration.diffusion_safety = 2.0;
    checkpoint.configuration.profile_bin_count = 16;
    checkpoint.configuration.profile_maximum_coordinate = 4.0;
    checkpoint.progress.step = step;
    checkpoint.progress.time = time;
    checkpoint.progress.bkm_sampled_integral =
        time * diagnostics.sampled_vorticity_max;
    checkpoint.progress.previous_vorticity_max =
        diagnostics.sampled_vorticity_max;
    checkpoint.progress.previous_diagnostic_time = time;
    checkpoint.progress.initial_critical_l3 = diagnostics.critical_l3_sample;
    checkpoint.progress.initial_critical_h_half = diagnostics.critical_h_half;
    checkpoint.progress.initial_sampled_vorticity =
        diagnostics.sampled_vorticity_max;
    checkpoint.progress.peak_high_shell_fraction =
        diagnostics.high_shell_energy_fraction;
    checkpoint.progress.peak_critical_l3 = diagnostics.critical_l3_sample;
    checkpoint.progress.peak_critical_h_half = diagnostics.critical_h_half;
    checkpoint.progress.peak_sampled_vorticity =
        diagnostics.sampled_vorticity_max;
    checkpoint.progress.maximum_production_to_dissipation =
        diagnostics.enstrophy_production_to_dissipation;
    checkpoint.progress.peak_forward_flux = 0.0;
    checkpoint.progress.minimum_time_step = 0.0005;
    checkpoint.progress.maximum_time_step = 0.0005;
    checkpoint.progress.maximum_cfl_bound = 0.1;
    checkpoint.progress.maximum_viscous_number = 0.1;
    checkpoint.initial_profile = profile;
    checkpoint.previous_profile = profile;
    checkpoint.state = state;
    return checkpoint;
}

void testCheckpointRoundTripAndRestartTrajectory() {
    const std::string checkpoint_path = "ns_cascade_checkpoint_test.bin";
    const ns_cascade::PseudospectralSystem system(16, 0.05, 5);
    ns_cascade::PseudospectralSystem::State uninterrupted =
        system.initialState(ns_cascade::InitialCondition::TaylorGreen, 1.0);
    ns_cascade::PseudospectralSystem::State split = uninterrupted;
    const double time_step = 0.0005;
    for (int step = 0; step < 8; ++step) {
        system.stepRungeKutta4(uninterrupted, time_step);
    }
    for (int step = 0; step < 3; ++step) {
        system.stepRungeKutta4(split, time_step);
    }

    ns_cascade::saveSimulationCheckpoint(
        checkpoint_path, makeTestCheckpoint(system, split, 3U, 3.0 * time_step));
    const ns_cascade::SimulationCheckpoint recovered =
        ns_cascade::loadSimulationCheckpoint(checkpoint_path);
    expect(recovered.configuration.grid_size == 16 &&
               recovered.configuration.cutoff == 5,
           "Checkpoint configuration was not preserved");
    expect(recovered.progress.step == 3U,
           "Checkpoint step count was not preserved");
    expectNear(recovered.progress.time, 3.0 * time_step, 0.0,
               "Checkpoint physical time was not preserved bit-for-bit");
    double round_trip_difference = 0.0;
    for (std::size_t i = 0; i < split.size(); ++i) {
        round_trip_difference = std::max(
            round_trip_difference, ns_cascade::norm(split[i] - recovered.state[i]));
    }
    expectNear(round_trip_difference, 0.0, 0.0,
               "Checkpoint changed a Fourier coefficient");

    ns_cascade::PseudospectralSystem::State restarted = recovered.state;
    for (int step = 3; step < 8; ++step) {
        system.stepRungeKutta4(restarted, time_step);
    }
    double trajectory_difference = 0.0;
    for (std::size_t i = 0; i < restarted.size(); ++i) {
        trajectory_difference = std::max(
            trajectory_difference,
            ns_cascade::norm(uninterrupted[i] - restarted[i]));
    }
    expectNear(trajectory_difference, 0.0, 0.0,
               "Restarted and uninterrupted trajectories differ");

    {
        std::fstream corrupt(checkpoint_path.c_str(),
                             std::ios::binary | std::ios::in | std::ios::out);
        expect(static_cast<bool>(corrupt), "Could not reopen checkpoint for test");
        corrupt.seekg(40);
        char byte = 0;
        corrupt.read(&byte, 1);
        byte ^= 0x01;
        corrupt.seekp(40);
        corrupt.write(&byte, 1);
    }
    bool corruption_rejected = false;
    try {
        ns_cascade::loadSimulationCheckpoint(checkpoint_path);
    } catch (const std::runtime_error&) {
        corruption_rejected = true;
    }
    std::remove(checkpoint_path.c_str());
    expect(corruption_rejected,
           "Checkpoint checksum failed to reject a changed byte");
}

}  // namespace

int main() {
    try {
        testGridAndDealiasingRules();
        testScalarFftRoundTrip();
        testInitialDataMatchesCompactOracle();
        testNonlinearTermMatchesCompactOracle();
        testCriticalAndEnstrophyDiagnosticsMatchCompactOracle();
        testMaximumCutoffMatchesOracleForFullSpectrum();
        testEnergyAndShellIdentities();
        testShortTrajectoryMatchesCompactOracle();
        testAbcNegativeControl();
        testVortexTubeInitialData();
        testAdaptiveTimeStepControl();
        testRescaledSpectrumProfile();
        testCheckpointRoundTripAndRestartTrajectory();
        std::cout << "All pseudospectral Navier-Stokes tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
