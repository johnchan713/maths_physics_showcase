#include "ns_cascade/candidate_evidence.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void expect(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

void testShiftedDecayingShear() {
    const double pi = 3.1415926535897932384626433832795;
    const ns_cascade::PseudospectralSystem coarse(8, 0.05);
    const ns_cascade::PseudospectralSystem fine(32, 0.05);
    ns_cascade::OptimizationState state = coarse.zeroState();
    // u=(cos(2y+pi/8),0,0): the exact shear solution has no stretching.
    const ns_cascade::Complex amplitude =
        0.5 * ns_cascade::Complex(std::cos(pi / 8), std::sin(pi / 8));
    state[coarse.indexOf(ns_cascade::WaveVector(0, 2, 0))].x = amplitude;
    state[coarse.indexOf(ns_cascade::WaveVector(0, -2, 0))].x = std::conj(amplitude);
    const auto native = coarse.diagnostics(state);
    const auto common = ns_cascade::candidateDiagnostics(coarse, state, fine);
    const auto padded = ns_cascade::liftOptimizationState(coarse, state, fine);
    const auto fine_values = ns_cascade::candidateDiagnostics(fine, padded, fine);
    expect(native.sampled_vorticity_max < 1.9,
           "Shifted shear should reveal coarse maximum undersampling");
    expect(std::abs(common.sampled_vorticity_max - 2.0) < 1e-13,
           "Padded shear samples should attain its known maximum");
    expect(std::abs(common.critical_l3_sample - fine_values.critical_l3_sample) < 1e-13 &&
           std::abs(common.sampled_vorticity_max - fine_values.sampled_vorticity_max) < 1e-13,
           "Identical Fourier fields must agree on common spatial samples");
    expect(common.energy == native.energy && common.enstrophy == native.enstrophy &&
           common.net_enstrophy_rate == native.net_enstrophy_rate &&
           common.high_shell_energy_fraction == native.high_shell_energy_fraction,
           "Padding must not alter the evolving truncation's spectral budget");
    expect(std::abs(common.energy - 0.25) < 1e-14 &&
           std::abs(common.enstrophy - 1.0) < 1e-14 &&
           std::abs(common.viscous_enstrophy_destruction - 0.4) < 1e-13 &&
           std::abs(common.nonlinear_enstrophy_production) < 1e-13,
           "Shear enstrophy budget must match analytic diffusion");
    for (int i = 0; i < 100; ++i) coarse.stepRungeKutta4(state, 0.001);
    const auto evolved = ns_cascade::candidateDiagnostics(coarse, state, fine);
    expect(std::abs(evolved.enstrophy - std::exp(-0.04)) < 1e-12 &&
           std::abs(evolved.sampled_vorticity_max - 2.0 * std::exp(-0.02)) < 1e-12,
           "Evidence must reproduce the exact decaying PDE shear trajectory");
    bool rejected = false;
    try {
        ns_cascade::candidateDiagnostics(fine, padded, coarse);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    expect(rejected, "Evidence must reject a sampling grid that drops modes");
}
}  // namespace

int main() {
    try {
        testShiftedDecayingShear();
        std::cout << "Candidate evidence analytic shear and sampling tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
