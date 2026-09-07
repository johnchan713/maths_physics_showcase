#include "ns_cascade/continuation.hpp"

#include <iostream>

namespace {
using namespace ns_cascade;
void expect(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
template <typename F> void rejects(F operation, const char* message) {
    bool failed = false;
    try { operation(); } catch (const std::exception&) { failed = true; }
    expect(failed, message);
}
void testShearAndRestart() {
    ContinuationConfiguration c;
    c.grid = c.sampling_grid = 16; c.cutoff = 5; c.dense_sampling_grid = 32;
    c.viscosity = 0.05; c.maximum_time_step = 0.001; c.observation_interval = 0.01;
    ContinuationSystem engine(c);
    auto state = engine.spectralSystem().zeroState();
    const double phase = std::acos(-1.0) / 8.0;
    const Complex amplitude = 0.5 * std::exp(Complex(0.0, phase));
    state[engine.spectralSystem().indexOf(WaveVector(0, 2, 0))].x = amplitude;
    state[engine.spectralSystem().indexOf(WaveVector(0, -2, 0))].x = std::conj(amplitude);
    auto full = engine.initialize(state);
    expect(full.initial.vorticity_sample < 1.97 && std::abs(full.initial.vorticity_dense - 2.0) < 1e-13,
           "Dense sampling must resolve the shifted shear maximum");
    engine.advance(full, 0.03, [](const ContinuationCheckpoint&) {});
    expect(std::abs(full.latest.energy - 0.25 * std::exp(-0.012)) < 1e-13 &&
           std::abs(full.latest.h_half - std::exp(-0.006)) < 1e-13 &&
           std::abs(full.latest.vorticity_dense - 2.0 * std::exp(-0.006)) < 1e-13 &&
           std::abs(full.latest.stretching) < 1e-13,
           "Continuation must reproduce the analytic Navier-Stokes shear solution");
    auto split = engine.initialize(state);
    engine.advance(split, 0.01, [](const ContinuationCheckpoint&) {});
    const std::string path = "ns-continuation-unit-test.chk";
    saveContinuationCheckpoint(path, split);
    split = loadContinuationCheckpoint(path);
    engine.advance(split, 0.03, [](const ContinuationCheckpoint&) {});
    expect(encodeContinuationCheckpoint(full) == encodeContinuationCheckpoint(split),
           "Split and uninterrupted continuations must be identical, including clock/baselines");
    c.fftw = false;
    ContinuationSystem radix(c);
    rejects([&] { radix.advance(split, 0.04, [](const ContinuationCheckpoint&) {}); },
            "Engine must reject a checkpoint from another backend");
    auto independent = radix.initialize(state);
    radix.advance(independent, 0.03, [](const ContinuationCheckpoint&) {});
    expect(optimizationStateNorm(addOptimizationStates(full.state, independent.state, -1.0)) < 1e-12,
           "Radix2 continuation must agree with FFTW");
    rejects([&] { engine.advance(split, 0.035, [](const ContinuationCheckpoint&) {}); },
            "Off-clock endpoint must be rejected");
    auto invalid = full;
    invalid.state[0].x = Complex(1.0, 0.0);
    rejects([&] { encodeContinuationCheckpoint(invalid); }, "Unretained mean must be rejected");
    invalid = full;
    invalid.state[1].x = Complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
    rejects([&] { encodeContinuationCheckpoint(invalid); }, "Nonfinite coefficient must be rejected");
    invalid = full;
    invalid.completed_observations = 2;
    rejects([&] { encodeContinuationCheckpoint(invalid); }, "Inconsistent checkpoint clock must be rejected");
    std::fstream corrupt(path.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    corrupt.seekg(40);
    char byte = 0; corrupt.read(&byte, 1); byte ^= 1;
    corrupt.seekp(40); corrupt.write(&byte, 1); corrupt.close();
    rejects([&] { loadContinuationCheckpoint(path); }, "Corrupt checkpoint must fail its checksum");
    saveContinuationCheckpoint(path, full);
    std::ofstream trailing(path.c_str(), std::ios::binary | std::ios::app);
    trailing.put('x'); trailing.close();
    rejects([&] { loadContinuationCheckpoint(path); }, "Trailing checkpoint bytes must be rejected");
    std::remove(path.c_str());
}

void testInterObservationCutoffStop() {
    ContinuationConfiguration c;
    c.grid = c.sampling_grid = c.dense_sampling_grid = 8; c.cutoff = 2;
    c.maximum_time_step = 0.001; c.cutoff_limit = 1e-16;
    ContinuationSystem engine(c);
    auto p = engine.initialize(engine.spectralSystem().initialState(InitialCondition::TaylorGreen, 1.0));
    int observations = 0;
    engine.advance(p, 0.02, [&](const ContinuationCheckpoint&) { ++observations; });
    expect(p.cutoff_stopped && p.step == 1 && p.latest.time < c.observation_interval && observations == 1,
           "Cutoff must stop at the first violating step, between scheduled observations");
    const std::string path = "ns-continuation-stop-test.chk";
    saveContinuationCheckpoint(path, p);
    const auto restored = loadContinuationCheckpoint(path);
    expect(encodeContinuationCheckpoint(p) == encodeContinuationCheckpoint(restored),
           "Failed cutoff state and reason must survive a checkpoint round trip");
    rejects([&] { engine.advance(p, 0.02, [](const ContinuationCheckpoint&) {}); },
            "Cutoff-stopped checkpoint must not restart");
    std::remove(path.c_str());
}

void testInteractingBackends() {
    ContinuationConfiguration c;
    c.grid = c.sampling_grid = c.dense_sampling_grid = 16; c.cutoff = 5;
    c.maximum_time_step = 0.0005; c.observation_interval = 0.005;
    ContinuationSystem fftw(c);
    const auto state = fftw.spectralSystem().initialState(InitialCondition::Deterministic, 1.0);
    c.cutoff_limit = 0.9;
    ContinuationSystem relaxed(c);
    auto a = relaxed.initialize(state);
    c.fftw = false;
    ContinuationSystem radix(c);
    auto b = radix.initialize(state);
    relaxed.advance(a, 0.01, [](const ContinuationCheckpoint&) {});
    radix.advance(b, 0.01, [](const ContinuationCheckpoint&) {});
    expect(optimizationStateNorm(addOptimizationStates(a.state, b.state, -1.0)) < 1e-12 &&
           std::abs(a.latest.stretching - b.latest.stretching) < 1e-10,
           "Interacting trajectory and stretching must agree across independent backends");
}

void testCompactRungeKuttaIdentity() {
    for (int grid : {8, 16, 32}) {
        const PseudospectralSystem seed(grid, 0.02);
        const FftwReferenceSystem fftw(grid, 0.02);
        for (const auto initial : {InitialCondition::TaylorGreen, InitialCondition::Deterministic}) {
            auto full = seed.initialState(initial, 10.0);
            auto compact = full;
            for (int step = 0; step < 12; ++step) {
                const double dt = fftw.chooseAdaptiveTimeStep(full, 0.0003).time_step;
                fftw.stepRungeKutta4(full, dt);
                fftw.stepRungeKutta4Compact(compact, dt);
                expect(optimizationStateNorm(addOptimizationStates(full, compact, -1.0)) == 0.0,
                       "Compact RK4 must preserve every Fourier coefficient bit for bit");
            }
        }
    }
}
}  // namespace

int main() {
    try {
        testShearAndRestart();
        testInterObservationCutoffStop();
        testInteractingBackends();
        testCompactRungeKuttaIdentity();
        std::cout << "Continuation analytic, independent evolution, restart, corruption and stop gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
