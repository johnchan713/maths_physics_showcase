// Bounded cost probe of the existing FFTW RK4 path, not a .14 trajectory.
#include "ns_cascade/fftw_reference.hpp"
#include "ns_cascade/optimization_state_csv.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>

using namespace ns_cascade;
using Clock = std::chrono::steady_clock;

struct Probe {
    int grid;
    int cutoff;
    OptimizationState state;
    std::vector<double> seconds;
    double initial_energy;
    double initial_h;
    double final_energy;
    double final_h;
    double divergence;
    double reality;
    double cutoff_fraction;
};

Probe run(const std::string& initial, int grid) {
    PseudospectralSystem spectral(grid, .02);
    FftwReferenceSystem fftw(grid, .02);
    Probe p;
    p.grid = grid;
    p.cutoff = spectral.cutoff();
    p.state = readOptimizationStateCsv(initial, spectral).state;
    p.initial_energy = spectral.energy(p.state);
    p.initial_h = spectral.criticalHOneHalf(p.state);
    for (int step = 0; step < 2; ++step) {
        const double dt = 1e-5;
        if (fftw.chooseAdaptiveTimeStep(p.state, dt, .4, 2.).time_step != dt) {
            throw std::runtime_error("Cost probe requires two matching safe steps");
        }
        const auto before = Clock::now();
        fftw.stepRungeKutta4Compact(p.state, dt);
        p.seconds.push_back(std::chrono::duration<double>(Clock::now() - before).count());
        std::cerr << "grid=" << grid << " step=" << step + 1
                  << " RK4_seconds=" << p.seconds.back() << std::endl;
    }
    p.final_energy = spectral.energy(p.state);
    p.final_h = spectral.criticalHOneHalf(p.state);
    p.divergence = spectral.divergenceDefect(p.state);
    p.reality = spectral.realityDefect(p.state);
    p.cutoff_fraction = spectral.cutoffShellEnergyFraction(p.state);
    if (!std::isfinite(p.final_energy) || p.final_energy > p.initial_energy * (1 + 1e-10) ||
        p.divergence > 1e-9 || p.reality > 1e-9 || p.cutoff_fraction > .008) {
        throw std::runtime_error("Cost-probe invariant failed");
    }
    return p;
}

void print(const Probe& p) {
    std::cout << "{\"grid\":" << p.grid << ",\"cutoff\":" << p.cutoff
              << ",\"steps\":2,\"final_time\":0.00002,\"rk4_seconds\":["
              << p.seconds[0] << ',' << p.seconds[1]
              << "],\"initial_energy\":" << p.initial_energy
              << ",\"initial_h_half\":" << p.initial_h
              << ",\"final_energy\":" << p.final_energy
              << ",\"final_h_half\":" << p.final_h
              << ",\"divergence_defect\":" << p.divergence
              << ",\"reality_defect\":" << p.reality
              << ",\"cutoff_fraction\":" << p.cutoff_fraction << '}';
}

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 4) throw std::runtime_error("Use INITIAL_CSV [COARSE_GRID FINE_GRID]");
        const int coarse = argc == 4 ? parseOptimizationStateCsvNumber<int>(argv[2], "coarse grid") : 128;
        const int fine = argc == 4 ? parseOptimizationStateCsvNumber<int>(argv[3], "fine grid") : 256;
        if (coarse < 16 || coarse > 128 || fine != 2 * coarse) {
            throw std::runtime_error("Probe needs a nested pair from 16/32 through 128/256");
        }
        const Probe a = run(argv[1], coarse);
        const Probe b = run(argv[1], fine);
        double squared_error = 0., squared_reference = 0., maximum_error = 0.;
        for (int x = 0; x < fine; ++x) {
            const int kx = x <= fine / 2 ? x : x - fine;
            for (int y = 0; y < fine; ++y) {
                const int ky = y <= fine / 2 ? y : y - fine;
                for (int z = 0; z < fine; ++z) {
                    const int kz = z <= fine / 2 ? z : z - fine;
                    const std::size_t index = (static_cast<std::size_t>(x) * fine + y) * fine + z;
                    ComplexVector reference;
                    if (std::max(std::abs(kx), std::max(std::abs(ky), std::abs(kz))) <= a.cutoff) {
                        const int cx = (kx + coarse) % coarse, cy = (ky + coarse) % coarse;
                        const int cz = (kz + coarse) % coarse;
                        reference = a.state[(static_cast<std::size_t>(cx) * coarse + cy) * coarse + cz];
                    }
                    const double error = normSquared(b.state[index] - reference);
                    squared_error += error;
                    squared_reference += normSquared(b.state[index]);
                    maximum_error = std::max(maximum_error, std::sqrt(error));
                }
            }
        }
        const double relative_error = std::sqrt(squared_error / squared_reference);
        std::cout << std::setprecision(17) << "{\"probes\":[";
        print(a); std::cout << ','; print(b);
        std::cout << "],\"relative_state_l2_gap\":" << relative_error
                  << ",\"maximum_coefficient_gap\":" << maximum_error << "}\n";
        if (!std::isfinite(relative_error) || relative_error > 1e-10) {
            throw std::runtime_error("Short-probe state agreement failed");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
