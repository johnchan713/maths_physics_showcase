#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    int grid_size = 16;
    int cutoff = 0;
    double viscosity = 0.05;
    double time_step = 0.0005;
    int steps = 200;
    int diagnostic_every = 20;
    double initial_energy = 1.0;
    ns_cascade::InitialCondition initial_condition =
        ns_cascade::InitialCondition::TaylorGreen;
    std::string output = "navier_stokes_fft.csv";
    std::string shell_output = "navier_stokes_fft_shells.csv";
};

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

ns_cascade::InitialCondition parseInitialCondition(const std::string& value) {
    if (value == "deterministic") return ns_cascade::InitialCondition::Deterministic;
    if (value == "taylor-green") return ns_cascade::InitialCondition::TaylorGreen;
    if (value == "abc") return ns_cascade::InitialCondition::ABC;
    throw std::invalid_argument(
        "Unknown initial condition: " + value +
        " (expected deterministic, taylor-green, or abc)");
}

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Strictly dealiased, radix-2 FFT Navier-Stokes experiment.\n\n"
        << "Options:\n"
        << "  --grid N                Power-of-two grid size (default: 16)\n"
        << "  --cutoff K              Retained cube cutoff; 0 uses floor((N-1)/3)\n"
        << "  --viscosity NU          Non-negative viscosity (default: 0.05)\n"
        << "  --dt DT                 RK4 step size (default: 0.0005)\n"
        << "  --steps N               Number of RK4 steps (default: 200)\n"
        << "  --diagnostic-every N    CSV sampling interval (default: 20)\n"
        << "  --energy E              Initial normalized energy (default: 1)\n"
        << "  --initial-condition C   deterministic, taylor-green, or abc\n"
        << "                            (default: taylor-green)\n"
        << "  --output PATH           Time-series CSV path\n"
        << "  --shell-output PATH     Shell-transfer CSV path\n"
        << "  --help                  Show this message\n";
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (flag == "--grid") {
            options.grid_size = parseNumber<int>(requireValue(i, argc, argv), flag);
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
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--initial-condition") {
            options.initial_condition =
                parseInitialCondition(requireValue(i, argc, argv));
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else if (flag == "--shell-output") {
            options.shell_output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    if (options.cutoff < 0) throw std::invalid_argument("--cutoff cannot be negative");
    if (options.viscosity < 0.0) {
        throw std::invalid_argument("--viscosity must be non-negative");
    }
    if (options.time_step <= 0.0) throw std::invalid_argument("--dt must be positive");
    if (options.steps < 1) throw std::invalid_argument("--steps must be >= 1");
    if (options.diagnostic_every < 1) {
        throw std::invalid_argument("--diagnostic-every must be >= 1");
    }
    if (options.initial_energy <= 0.0) {
        throw std::invalid_argument("--energy must be positive");
    }
    if (options.output.empty() || options.shell_output.empty()) {
        throw std::invalid_argument("Output paths cannot be empty");
    }
    if (options.output == options.shell_output) {
        throw std::invalid_argument("--output and --shell-output must be different");
    }
    return options;
}

void writeTimeHeader(std::ostream& output) {
    output
        << "step,time,grid_size,cutoff,energy,enstrophy,palinstrophy,"
        << "critical_l3_sample,sampled_vorticity_max,"
        << "vorticity_sup_upper_bound,spectral_centroid,"
        << "high_shell_energy_fraction,divergence_defect,reality_defect,"
        << "energy_balance_residual,bkm_sampled_integral\n";
}

void writeTimeRow(std::ostream& output,
                  int step,
                  double time,
                  const ns_cascade::PseudospectralSystem& system,
                  const ns_cascade::PseudospectralSystem::Diagnostics& values,
                  double bkm_integral) {
    output << step << ',' << time << ',' << system.gridSize() << ','
           << system.cutoff() << ',' << values.energy << ',' << values.enstrophy
           << ',' << values.palinstrophy << ',' << values.critical_l3_sample << ','
           << values.sampled_vorticity_max << ','
           << values.vorticity_sup_upper_bound << ',' << values.spectral_centroid
           << ',' << values.high_shell_energy_fraction << ','
           << values.divergence_defect << ',' << values.reality_defect << ','
           << values.energy_balance_residual << ',' << bkm_integral << '\n';
}

void writeShellHeader(std::ostream& output) {
    output
        << "initial_condition,step,time,grid_size,cutoff,shell,lower_radius,"
        << "upper_radius,shell_energy,nonlinear_transfer,viscous_dissipation,"
        << "forward_flux\n";
}

double writeShellRows(
    std::ostream& output,
    ns_cascade::InitialCondition initial_condition,
    int step,
    double time,
    const ns_cascade::PseudospectralSystem& system,
    const std::vector<ns_cascade::PseudospectralSystem::ShellDiagnostics>& shells) {
    double largest_forward_flux = 0.0;
    for (std::size_t i = 0; i < shells.size(); ++i) {
        const ns_cascade::PseudospectralSystem::ShellDiagnostics& values = shells[i];
        output << ns_cascade::initialConditionName(initial_condition) << ',' << step
               << ',' << time << ',' << system.gridSize() << ',' << system.cutoff()
               << ',' << values.shell << ',' << values.lower_radius << ','
               << values.upper_radius << ',' << values.energy << ','
               << values.nonlinear_transfer << ',' << values.viscous_dissipation
               << ',' << values.forward_flux << '\n';
        largest_forward_flux =
            std::max(largest_forward_flux, values.forward_flux);
    }
    return largest_forward_flux;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const ns_cascade::PseudospectralSystem system(
            options.grid_size, options.viscosity, options.cutoff);
        ns_cascade::PseudospectralSystem::State state =
            system.initialState(options.initial_condition, options.initial_energy);

        std::ofstream csv(options.output.c_str());
        std::ofstream shell_csv(options.shell_output.c_str());
        if (!csv) throw std::runtime_error("Cannot open output file: " + options.output);
        if (!shell_csv) {
            throw std::runtime_error("Cannot open shell output file: " +
                                     options.shell_output);
        }
        csv << std::setprecision(17);
        shell_csv << std::setprecision(17);
        writeTimeHeader(csv);
        writeShellHeader(shell_csv);

        ns_cascade::PseudospectralSystem::Diagnostics diagnostics =
            system.diagnostics(state);
        double bkm_sampled_integral = 0.0;
        double previous_vorticity_max = diagnostics.sampled_vorticity_max;
        double peak_high_shell_fraction = diagnostics.high_shell_energy_fraction;
        double peak_critical_l3 = diagnostics.critical_l3_sample;
        int previous_diagnostic_step = 0;
        writeTimeRow(csv, 0, 0.0, system, diagnostics, bkm_sampled_integral);
        double peak_forward_flux = writeShellRows(
            shell_csv,
            options.initial_condition,
            0,
            0.0,
            system,
            system.shellDiagnostics(state));

        const std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
        for (int step = 1; step <= options.steps; ++step) {
            system.stepRungeKutta4(state, options.time_step);
            if (step % options.diagnostic_every != 0 && step != options.steps) continue;

            diagnostics = system.diagnostics(state);
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
            writeTimeRow(csv,
                         step,
                         step * options.time_step,
                         system,
                         diagnostics,
                         bkm_sampled_integral);
            peak_forward_flux = std::max(
                peak_forward_flux,
                writeShellRows(shell_csv,
                               options.initial_condition,
                               step,
                               step * options.time_step,
                               system,
                               system.shellDiagnostics(state)));
        }
        const std::chrono::duration<double> elapsed =
            std::chrono::steady_clock::now() - start;

        std::cout << std::setprecision(8)
                  << "Initial condition: "
                  << ns_cascade::initialConditionName(options.initial_condition) << '\n'
                  << "Grid: " << system.gridSize() << "^3, retained cutoff K="
                  << system.cutoff() << " (safe maximum "
                  << system.safeDealiasCutoff() << ")\n"
                  << "Retained non-zero modes: " << system.modeCount() << '\n'
                  << "Final normalized energy: " << diagnostics.energy << '\n'
                  << "Peak sampled L3 norm: " << peak_critical_l3 << '\n'
                  << "Sampled BKM integral: " << bkm_sampled_integral << '\n'
                  << "Peak cutoff-shell energy fraction: "
                  << peak_high_shell_fraction << '\n'
                  << "Peak forward shell flux: " << peak_forward_flux << '\n'
                  << "Evolution time: " << elapsed.count() << " seconds\n"
                  << "Diagnostics written to " << options.output << '\n'
                  << "Shell transfers written to " << options.shell_output << '\n';
        if (peak_high_shell_fraction > 0.01) {
            std::cout
                << "Resolution warning: more than 1% of energy reached the cutoff shell.\n";
        }
        std::cout
            << "This dealiased finite grid still cannot prove regularity or blow-up.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
