#ifndef NS_CASCADE_OPTIMIZATION_STATE_CSV_HPP
#define NS_CASCADE_OPTIMIZATION_STATE_CSV_HPP

#include "ns_cascade/state_optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ns_cascade {

struct OptimizationStateCsvMetadata {
    std::string family;
    int source_grid = 0;
    int simulation_cutoff = 0;
    int seed_bandwidth = 0;
    double target_energy = 0.0;
};

struct LoadedOptimizationState {
    OptimizationStateCsvMetadata metadata;
    OptimizationState state;
};

inline const char* optimizationStateCsvHeader() {
    return "family,source_grid,simulation_cutoff,seed_bandwidth,"
           "target_energy,kx,ky,kz,ux_real,ux_imag,uy_real,uy_imag,"
           "uz_real,uz_imag";
}

inline std::vector<std::string> splitOptimizationStateCsvRow(
    const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    if (!line.empty() && line[line.size() - 1U] == ',') fields.push_back("");
    return fields;
}

template <typename T>
inline T parseOptimizationStateCsvNumber(const std::string& text,
                                         const std::string& field) {
    std::istringstream stream(text);
    T value = T();
    char trailing = '\0';
    if (!(stream >> value) || (stream >> trailing)) {
        throw std::runtime_error(
            "Invalid optimization-state CSV " + field + ": " + text);
    }
    return value;
}

inline double optimizationStateCsvRelativeDifference(double left,
                                                      double right) {
    return std::abs(left - right) /
           std::max(1e-300, std::max(std::abs(left), std::abs(right)));
}

inline LoadedOptimizationState readOptimizationStateCsv(
    const std::string& path,
    const PseudospectralSystem& system) {
    std::ifstream input(path.c_str());
    if (!input) {
        throw std::runtime_error("Could not open optimization-state CSV: " +
                                 path);
    }
    std::string line;
    if (!std::getline(input, line) || line != optimizationStateCsvHeader()) {
        throw std::runtime_error("Optimization-state CSV header is invalid");
    }

    LoadedOptimizationState loaded;
    loaded.state = system.zeroState();
    std::vector<bool> seen(system.gridPointCount(), false);
    std::size_t rows = 0U;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields =
            splitOptimizationStateCsvRow(line);
        if (fields.size() != 14U) {
            throw std::runtime_error(
                "Optimization-state CSV row has the wrong width");
        }
        const int row_grid = parseOptimizationStateCsvNumber<int>(
            fields[1], "source grid");
        const int row_cutoff = parseOptimizationStateCsvNumber<int>(
            fields[2], "simulation cutoff");
        const int row_bandwidth = parseOptimizationStateCsvNumber<int>(
            fields[3], "seed bandwidth");
        const double row_energy = parseOptimizationStateCsvNumber<double>(
            fields[4], "target energy");
        if (rows == 0U) {
            loaded.metadata.family = fields[0];
            loaded.metadata.source_grid = row_grid;
            loaded.metadata.simulation_cutoff = row_cutoff;
            loaded.metadata.seed_bandwidth = row_bandwidth;
            loaded.metadata.target_energy = row_energy;
        }
        if (fields[0].empty() || fields[0] != loaded.metadata.family ||
            row_grid != loaded.metadata.source_grid ||
            row_cutoff != loaded.metadata.simulation_cutoff ||
            row_bandwidth != loaded.metadata.seed_bandwidth ||
            !std::isfinite(row_energy) || row_energy <= 0.0 ||
            optimizationStateCsvRelativeDifference(
                row_energy, loaded.metadata.target_energy) > 1e-13) {
            throw std::runtime_error(
                "Optimization-state CSV metadata is inconsistent");
        }
        const bool source_grid_is_power_of_two =
            row_grid >= 8 && (row_grid & (row_grid - 1)) == 0;
        if (!source_grid_is_power_of_two || row_bandwidth < 1 ||
            row_bandwidth > system.cutoff() ||
            row_cutoff < row_bandwidth ||
            row_cutoff > (row_grid - 1) / 3) {
            throw std::runtime_error(
                "Optimization-state CSV spectral metadata is invalid");
        }

        const WaveVector wave(
            parseOptimizationStateCsvNumber<int>(fields[5], "kx"),
            parseOptimizationStateCsvNumber<int>(fields[6], "ky"),
            parseOptimizationStateCsvNumber<int>(fields[7], "kz"));
        if (wave.normSquared() == 0 ||
            stateMaximumComponent(wave) > row_bandwidth) {
            throw std::runtime_error(
                "Optimization-state CSV contains a forbidden mode");
        }
        const std::size_t index = system.indexOf(wave);
        if (seen[index]) {
            throw std::runtime_error(
                "Optimization-state CSV repeats a Fourier mode");
        }
        const double ux_real = parseOptimizationStateCsvNumber<double>(
            fields[8], "ux real");
        const double ux_imag = parseOptimizationStateCsvNumber<double>(
            fields[9], "ux imaginary");
        const double uy_real = parseOptimizationStateCsvNumber<double>(
            fields[10], "uy real");
        const double uy_imag = parseOptimizationStateCsvNumber<double>(
            fields[11], "uy imaginary");
        const double uz_real = parseOptimizationStateCsvNumber<double>(
            fields[12], "uz real");
        const double uz_imag = parseOptimizationStateCsvNumber<double>(
            fields[13], "uz imaginary");
        const double components[] = {
            ux_real, ux_imag, uy_real, uy_imag, uz_real, uz_imag};
        for (std::size_t component = 0;
             component < sizeof(components) / sizeof(components[0]);
             ++component) {
            if (!std::isfinite(components[component])) {
                throw std::runtime_error(
                    "Optimization-state CSV contains a non-finite coefficient");
            }
        }
        loaded.state[index] = ComplexVector(
            Complex(ux_real, ux_imag),
            Complex(uy_real, uy_imag),
            Complex(uz_real, uz_imag));
        seen[index] = true;
        ++rows;
    }
    if (!input.eof()) {
        throw std::runtime_error("Failed while reading optimization-state CSV");
    }
    if (rows == 0U) {
        throw std::runtime_error(
            "Optimization-state CSV contains no coefficients");
    }
    const std::size_t expected_rows =
        stateOptimizationDegreesOfFreedom(
            loaded.metadata.seed_bandwidth) /
        2U;
    if (rows != expected_rows ||
        optimizationStateCsvRelativeDifference(
            system.energy(loaded.state),
            loaded.metadata.target_energy) > 2e-13 ||
        system.divergenceDefect(loaded.state) > 2e-12 ||
        system.realityDefect(loaded.state) > 2e-12) {
        throw std::runtime_error(
            "Optimization-state CSV fails completeness or invariant checks");
    }
    return loaded;
}

inline void writeOptimizationStateCsv(
    const std::string& path,
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const std::string& family,
    int seed_bandwidth,
    double target_energy) {
    requireOptimizationState(system, state);
    if (path.empty() || family.empty() || seed_bandwidth < 1 ||
        seed_bandwidth > system.cutoff() || !std::isfinite(target_energy) ||
        target_energy <= 0.0 ||
        optimizationStateCsvRelativeDifference(
            system.energy(state), target_energy) > 2e-13 ||
        system.divergenceDefect(state) > 2e-12 ||
        system.realityDefect(state) > 2e-12) {
        throw std::invalid_argument(
            "Optimization state or CSV metadata is invalid");
    }
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        if ((modes[i].normSquared() == 0 ||
             stateMaximumComponent(modes[i]) > seed_bandwidth) &&
            normSquared(state[i]) != 0.0) {
            throw std::invalid_argument(
                "Optimization state contains an out-of-band coefficient");
        }
    }
    std::ofstream output(path.c_str());
    if (!output) {
        throw std::runtime_error("Could not open optimization-state CSV: " +
                                 path);
    }
    output << std::setprecision(17) << optimizationStateCsvHeader() << '\n';
    for (std::size_t i = 0; i < state.size(); ++i) {
        const WaveVector& wave = modes[i];
        if (wave.normSquared() == 0 ||
            stateMaximumComponent(wave) > seed_bandwidth) {
            continue;
        }
        output << family << ',' << system.gridSize() << ',' << system.cutoff()
               << ',' << seed_bandwidth << ',' << target_energy << ','
               << wave.x << ',' << wave.y << ',' << wave.z << ','
               << std::real(state[i].x) << ',' << std::imag(state[i].x) << ','
               << std::real(state[i].y) << ',' << std::imag(state[i].y) << ','
               << std::real(state[i].z) << ',' << std::imag(state[i].z) << '\n';
    }
    output.close();
    if (!output) {
        throw std::runtime_error("Failed while writing optimization-state CSV");
    }
}

}  // namespace ns_cascade

#endif
