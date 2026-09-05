#include "ns_cascade/galerkin.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

struct Options {
    int cutoff = 2;
    double viscosity = 0.05;
    double time_step = 0.0005;
    int steps = 400;
    int diagnostic_every = 20;
    int sample_points = 0;
    double initial_energy = 1.0;
    std::string output = "navier_stokes_cascade.csv";
};

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Small, auditable Fourier-Galerkin Navier-Stokes experiment.\n\n"
        << "Options:\n"
        << "  --cutoff N              Fourier cube cutoff (default: 2)\n"
        << "  --viscosity NU          Non-negative viscosity (default: 0.05)\n"
        << "  --dt DT                 RK4 step size (default: 0.0005)\n"
        << "  --steps N               Number of RK4 steps (default: 400)\n"
        << "  --diagnostic-every N    CSV sampling interval (default: 20)\n"
        << "  --sample-points N       Grid points per axis; 0 selects a safe default\n"
        << "  --energy E              Initial normalized kinetic energy (default: 1)\n"
        << "  --output PATH           CSV path (default: navier_stokes_cascade.csv)\n"
        << "  --help                  Show this message\n";
}

template <typename T>
T parseNumber(const std::string& text, const std::string& flag) {
    std::istringstream stream(text);
    T value;
    char trailing;
    if (!(stream >> value) || (stream >> trailing)) {
        throw std::invalid_argument("Invalid value for " + flag + ": " + text);
    }
    return value;
}

std::string requireValue(int& index, int argc, char** argv) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(std::string("Missing value after ") + argv[index]);
    }
    ++index;
    return argv[index];
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (flag == "--cutoff") {
            options.cutoff = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--viscosity") {
            options.viscosity = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--dt") {
            options.time_step = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--steps") {
            options.steps = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diagnostic-every") {
            options.diagnostic_every =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--sample-points") {
            options.sample_points = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    if (options.cutoff < 1) throw std::invalid_argument("--cutoff must be >= 1");
    if (options.viscosity < 0.0) {
        throw std::invalid_argument("--viscosity must be non-negative");
    }
    if (options.time_step <= 0.0) throw std::invalid_argument("--dt must be positive");
    if (options.steps < 1) throw std::invalid_argument("--steps must be >= 1");
    if (options.diagnostic_every < 1) {
        throw std::invalid_argument("--diagnostic-every must be >= 1");
    }
    if (options.sample_points != 0 &&
        options.sample_points < 2 * options.cutoff + 1) {
        throw std::invalid_argument("--sample-points must be 0 or >= 2*cutoff+1");
    }
    if (options.initial_energy <= 0.0) {
        throw std::invalid_argument("--energy must be positive");
    }
    if (options.output.empty()) throw std::invalid_argument("--output cannot be empty");
    return options;
}

void writeHeader(std::ostream& output) {
    output
        << "step,time,energy,enstrophy,palinstrophy,critical_l3_sample,"
        << "sampled_vorticity_max,vorticity_sup_upper_bound,spectral_centroid,"
        << "high_shell_energy_fraction,divergence_defect,reality_defect,"
        << "energy_balance_residual,bkm_sampled_integral\n";
}

void writeRow(std::ostream& output,
              int step,
              double time,
              const ns_cascade::GalerkinSystem::Diagnostics& values,
              double bkm_integral) {
    output << step << ',' << time << ',' << values.energy << ',' << values.enstrophy
           << ',' << values.palinstrophy << ',' << values.critical_l3_sample << ','
           << values.sampled_vorticity_max << ','
           << values.vorticity_sup_upper_bound << ',' << values.spectral_centroid
           << ',' << values.high_shell_energy_fraction << ','
           << values.divergence_defect << ',' << values.reality_defect << ','
           << values.energy_balance_residual << ',' << bkm_integral << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const ns_cascade::GalerkinSystem system(options.cutoff, options.viscosity);
        ns_cascade::GalerkinSystem::State state =
            system.deterministicLowModeState(options.initial_energy);

        std::ofstream csv(options.output.c_str());
        if (!csv) throw std::runtime_error("Cannot open output file: " + options.output);
        csv << std::setprecision(17);
        writeHeader(csv);

        ns_cascade::GalerkinSystem::Diagnostics diagnostics =
            system.diagnostics(state, options.sample_points);
        double bkm_sampled_integral = 0.0;
        double previous_vorticity_max = diagnostics.sampled_vorticity_max;
        double peak_high_shell_fraction = diagnostics.high_shell_energy_fraction;
        double peak_critical_l3 = diagnostics.critical_l3_sample;
        int previous_diagnostic_step = 0;
        writeRow(csv, 0, 0.0, diagnostics, bkm_sampled_integral);

        for (int step = 1; step <= options.steps; ++step) {
            system.stepRungeKutta4(state, options.time_step);
            if (step % options.diagnostic_every != 0 && step != options.steps) continue;

            diagnostics = system.diagnostics(state, options.sample_points);
            const double interval =
                (step - previous_diagnostic_step) * options.time_step;
            bkm_sampled_integral +=
                0.5 * interval *
                (previous_vorticity_max + diagnostics.sampled_vorticity_max);
            previous_vorticity_max = diagnostics.sampled_vorticity_max;
            previous_diagnostic_step = step;
            peak_high_shell_fraction =
                std::max(peak_high_shell_fraction,
                         diagnostics.high_shell_energy_fraction);
            peak_critical_l3 =
                std::max(peak_critical_l3, diagnostics.critical_l3_sample);
            writeRow(csv,
                     step,
                     step * options.time_step,
                     diagnostics,
                     bkm_sampled_integral);
        }

        std::cout << std::setprecision(8)
                  << "Completed " << options.steps << " steps with "
                  << system.modeCount() << " non-zero Fourier modes.\n"
                  << "Final normalized energy: " << diagnostics.energy << '\n'
                  << "Peak sampled L3 norm: " << peak_critical_l3 << '\n'
                  << "Sampled BKM integral: " << bkm_sampled_integral << '\n'
                  << "Peak cutoff-shell energy fraction: "
                  << peak_high_shell_fraction << '\n'
                  << "Diagnostics written to " << options.output << '\n';

        if (peak_high_shell_fraction > 0.01) {
            std::cout
                << "Resolution warning: more than 1% of energy reached the cutoff shell.\n";
        }
        std::cout
            << "This finite-dimensional run cannot prove regularity or blow-up.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
