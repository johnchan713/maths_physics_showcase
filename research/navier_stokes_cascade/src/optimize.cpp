#include "ns_cascade/parameter_optimizer.hpp"
#include "ns_cascade/pseudospectral.hpp"
#include "ns_cascade/spectral_profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const double kPi = 3.1415926535897932384626433832795;

struct Options {
    int grid_size = 16;
    int cutoff = 0;
    double viscosity = 0.02;
    double maximum_time_step = 0.005;
    double final_time = 0.08;
    double target_cfl = 0.35;
    double diffusion_safety = 2.0;
    double initial_energy = 10.0;
    int carrier_wavenumber = 1;
    double initial_width = 1.1;
    double initial_weight = 0.75;
    double initial_phase = kPi / 3.0;
    int iterations = 2;
    double gradient_step = 0.03;
    double gradient_tolerance = 0.15;
    double initial_line_step = 0.15;
    double minimum_line_step = 0.002;
    double cutoff_fraction_threshold = 0.01;
    int profile_bin_count = 32;
    double profile_maximum_coordinate = 4.0;
    std::string output = "navier_stokes_packet_optimization.csv";
};

struct Evaluation {
    double objective = -std::numeric_limits<double>::infinity();
    bool valid = false;
    int steps = 0;
    double width = 0.0;
    double weight = 0.0;
    double phase = 0.0;
    double peak_h_ratio = 0.0;
    double peak_l3_ratio = 0.0;
    double final_h_ratio = 0.0;
    double final_l3_ratio = 0.0;
    double peak_vorticity_ratio = 0.0;
    double characteristic_ratio = 0.0;
    double profile_l1 = 0.0;
    double profile_drift = 0.0;
    double peak_cutoff_fraction = 0.0;
    double maximum_divergence_defect = 0.0;
    double maximum_reality_defect = 0.0;
};

template <typename T>
T parseNumber(const std::string& text, const std::string& flag) {
    std::istringstream stream(text);
    T value = T();
    char trailing = '\0';
    if (!(stream >> value) || (stream >> trailing)) {
        throw std::invalid_argument("Invalid value for " + flag + ": " + text);
    }
    return value;
}

std::string requireValue(int& index, int argc, char** argv) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(std::string("Missing value after ") +
                                    argv[index]);
    }
    return argv[++index];
}

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Optimize three continuous wave-packet coordinates with checked "
           "finite differences and backtracking. This is a numerical search, "
           "not a proof.\n\n"
        << "  --grid N                 FFT grid (default: 16)\n"
        << "  --cutoff K               Retained cutoff; zero selects safe maximum\n"
        << "  --viscosity NU           Positive viscosity (default: 0.02)\n"
        << "  --dt DT                  Maximum adaptive step (default: 0.005)\n"
        << "  --final-time T           Objective horizon (default: 0.08)\n"
        << "  --cfl C                  Conservative CFL target (default: 0.35)\n"
        << "  --diffusion-safety S     Viscous RK4 bound (default: 2)\n"
        << "  --energy E               Initial normalized energy (default: 10)\n"
        << "  --carrier-mode M         Fixed integer carrier (default: 1)\n"
        << "  --initial-width W        Initial envelope width (default: 1.1)\n"
        << "  --initial-weight W       Initial secondary weight (default: 0.75)\n"
        << "  --initial-phase P        Initial phase (default: pi/3)\n"
        << "  --iterations N           Maximum accepted steps (default: 2)\n"
        << "  --gradient-step H        Normalized difference step (default: 0.03)\n"
        << "  --gradient-tolerance R   Coarse/refined disagreement gate\n"
        << "  --line-step A            Initial normalized line step (default: 0.15)\n"
        << "  --minimum-line-step A    Backtracking stop (default: 0.002)\n"
        << "  --cutoff-threshold F     Hard cutoff-energy gate (default: 0.01)\n"
        << "  --output PATH            Optimization trace CSV\n"
        << "  --help                   Show this message\n";
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
            options.maximum_time_step =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--final-time") {
            options.final_time = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cfl") {
            options.target_cfl = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diffusion-safety") {
            options.diffusion_safety =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--carrier-mode") {
            options.carrier_wavenumber =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--initial-width") {
            options.initial_width =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--initial-weight") {
            options.initial_weight =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--initial-phase") {
            options.initial_phase =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--iterations") {
            options.iterations = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--gradient-step") {
            options.gradient_step =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--gradient-tolerance") {
            options.gradient_tolerance =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--line-step") {
            options.initial_line_step =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--minimum-line-step") {
            options.minimum_line_step =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff-threshold") {
            options.cutoff_fraction_threshold =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    if (options.cutoff < 0 || options.iterations < 0 ||
        options.carrier_wavenumber < 1 || options.output.empty()) {
        throw std::invalid_argument("Integer or output option is invalid");
    }
    const double positive_values[] = {
        options.viscosity,
        options.maximum_time_step,
        options.final_time,
        options.target_cfl,
        options.diffusion_safety,
        options.initial_energy,
        options.gradient_step,
        options.gradient_tolerance,
        options.initial_line_step,
        options.minimum_line_step,
        options.cutoff_fraction_threshold};
    for (std::size_t i = 0;
         i < sizeof(positive_values) / sizeof(positive_values[0]); ++i) {
        if (!std::isfinite(positive_values[i]) || positive_values[i] <= 0.0) {
            throw std::invalid_argument("Positive optimizer option is invalid");
        }
    }
    if (options.gradient_step >= 0.25 ||
        options.cutoff_fraction_threshold >= 1.0 ||
        options.minimum_line_step > options.initial_line_step) {
        throw std::invalid_argument("Optimizer scale or threshold is invalid");
    }
    return options;
}

double mapCoordinate(double normalized, double lower, double upper) {
    return lower + normalized * (upper - lower);
}

double normalizeCoordinate(double value, double lower, double upper) {
    return (value - lower) / (upper - lower);
}

std::string pointKey(const std::vector<double>& point) {
    std::ostringstream stream;
    stream << std::setprecision(17);
    for (std::size_t i = 0; i < point.size(); ++i) stream << point[i] << ';';
    return stream.str();
}

class TrajectoryObjective {
public:
    explicit TrajectoryObjective(const Options& options)
        : options_(options),
          system_(options.grid_size, options.viscosity, options.cutoff) {}

    double operator()(const std::vector<double>& point) {
        return get(point).objective;
    }

    const Evaluation& get(const std::vector<double>& point) {
        const std::string key = pointKey(point);
        std::map<std::string, Evaluation>::const_iterator found = cache_.find(key);
        if (found != cache_.end()) return found->second;
        return cache_.insert(std::make_pair(key, evaluate(point))).first->second;
    }

    int evaluationCount() const { return static_cast<int>(cache_.size()); }

private:
    const Options& options_;
    ns_cascade::PseudospectralSystem system_;
    std::map<std::string, Evaluation> cache_;

    Evaluation evaluate(const std::vector<double>& point) const {
        if (point.size() != 3) {
            throw std::invalid_argument("Wave-packet optimizer expects three coordinates");
        }
        Evaluation result;
        result.width = mapCoordinate(point[0], 0.7, 1.5);
        result.weight = mapCoordinate(point[1], 0.5, 1.5);
        result.phase = mapCoordinate(point[2], -kPi, kPi);
        ns_cascade::PseudospectralSystem::State state =
            system_.interactingWavePacketState(
                ns_cascade::WavePacketParameters(
                    result.width,
                    options_.carrier_wavenumber,
                    result.weight,
                    result.phase),
                options_.initial_energy);
        const ns_cascade::PseudospectralSystem::Diagnostics initial =
            system_.diagnostics(state);
        const ns_cascade::SpectrumProfile initial_profile =
            ns_cascade::rescaledSpectrumProfile(
                system_, state, options_.profile_bin_count,
                options_.profile_maximum_coordinate);
        double peak_h = initial.critical_h_half;
        double peak_l3 = initial.critical_l3_sample;
        double peak_vorticity = initial.sampled_vorticity_max;
        result.peak_cutoff_fraction = initial.high_shell_energy_fraction;
        result.maximum_divergence_defect = initial.divergence_defect;
        result.maximum_reality_defect = initial.reality_defect;

        double time = 0.0;
        while (time < options_.final_time) {
            const ns_cascade::AdaptiveStepInfo step =
                system_.chooseAdaptiveTimeStep(
                    state,
                    options_.maximum_time_step,
                    options_.target_cfl,
                    options_.diffusion_safety);
            const double dt = std::min(step.time_step, options_.final_time - time);
            if (!std::isfinite(dt) || dt <= 0.0) {
                throw std::runtime_error("Optimizer selected an invalid timestep");
            }
            system_.stepRungeKutta4(state, dt);
            time += dt;
            ++result.steps;
            const ns_cascade::PseudospectralSystem::Diagnostics current =
                system_.diagnostics(state);
            peak_h = std::max(peak_h, current.critical_h_half);
            peak_l3 = std::max(peak_l3, current.critical_l3_sample);
            peak_vorticity =
                std::max(peak_vorticity, current.sampled_vorticity_max);
            result.peak_cutoff_fraction = std::max(
                result.peak_cutoff_fraction,
                current.high_shell_energy_fraction);
            result.maximum_divergence_defect = std::max(
                result.maximum_divergence_defect, current.divergence_defect);
            result.maximum_reality_defect = std::max(
                result.maximum_reality_defect, current.reality_defect);
        }

        const ns_cascade::PseudospectralSystem::Diagnostics final =
            system_.diagnostics(state);
        const ns_cascade::SpectrumProfile final_profile =
            ns_cascade::rescaledSpectrumProfile(
                system_, state, options_.profile_bin_count,
                options_.profile_maximum_coordinate);
        const ns_cascade::SpectrumProfileChange profile_change =
            ns_cascade::spectrumProfileChange(final_profile, initial_profile);
        result.peak_h_ratio = peak_h / initial.critical_h_half;
        result.peak_l3_ratio = peak_l3 / initial.critical_l3_sample;
        result.final_h_ratio = final.critical_h_half / initial.critical_h_half;
        result.final_l3_ratio =
            final.critical_l3_sample / initial.critical_l3_sample;
        result.peak_vorticity_ratio =
            peak_vorticity / initial.sampled_vorticity_max;
        result.characteristic_ratio =
            final_profile.characteristic_wavenumber /
            initial_profile.characteristic_wavenumber;
        result.profile_l1 = profile_change.l1_distance;
        double profile_squared_distance = 0.0;
        for (std::size_t i = 0;
             i < initial_profile.energy_fractions.size(); ++i) {
            const double difference =
                final_profile.energy_fractions[i] -
                initial_profile.energy_fractions[i];
            profile_squared_distance += difference * difference;
        }
        const double scale_log = std::log(result.characteristic_ratio);
        // A smooth endpoint analogue of scale-normalized profile drift.  Tiny
        // floors avoid a singular quotient without changing resolved runs.
        result.profile_drift =
            std::sqrt(profile_squared_distance + 1e-16) /
            std::sqrt(scale_log * scale_log + 1e-8);

        const double critical_reward =
            std::log(result.final_h_ratio) +
            0.10 * std::log(result.final_l3_ratio);
        const double scale_reward = 0.15 * scale_log;
        const double profile_cost = 0.05 * std::log1p(result.profile_drift);
        const double cutoff_cost = 0.04 * std::log1p(
            final.high_shell_energy_fraction /
            options_.cutoff_fraction_threshold);
        result.objective = critical_reward + scale_reward -
                           profile_cost - cutoff_cost;
        result.valid =
            std::isfinite(result.objective) &&
            result.peak_cutoff_fraction <= options_.cutoff_fraction_threshold &&
            result.maximum_divergence_defect <= 1e-10 &&
            result.maximum_reality_defect <= 1e-10;
        return result;
    }
};

double gradientDotStep(const std::vector<double>& gradient,
                       const std::vector<double>& from,
                       const std::vector<double>& to) {
    double product = 0.0;
    for (std::size_t i = 0; i < gradient.size(); ++i) {
        product += gradient[i] * (to[i] - from[i]);
    }
    return product;
}

void writeRow(std::ostream& output,
              int iteration,
              const char* stage,
              const Evaluation& value,
              const ns_cascade::FiniteDifferenceGradient* gradient,
              double line_step,
              int evaluations) {
    output << iteration << ',' << stage << ',' << value.objective << ','
           << (value.valid ? "true" : "false") << ',' << value.width << ','
           << value.weight << ',' << value.phase << ',' << value.steps << ','
           << value.peak_h_ratio << ',' << value.peak_l3_ratio << ','
           << value.final_h_ratio << ',' << value.final_l3_ratio << ','
           << value.peak_vorticity_ratio << ',' << value.characteristic_ratio
           << ',' << value.profile_l1 << ',' << value.profile_drift << ','
           << value.peak_cutoff_fraction << ','
           << value.maximum_divergence_defect << ','
           << value.maximum_reality_defect << ',';
    if (gradient == NULL) {
        output << ",,,,";
    } else {
        output << gradient->extrapolated[0] << ','
               << gradient->extrapolated[1] << ','
               << gradient->extrapolated[2] << ','
               << gradient->maximum_relative_disagreement << ',';
    }
    output << line_step << ',' << evaluations << '\n';
}

int run(const Options& options) {
    TrajectoryObjective objective(options);
    std::vector<double> point = {
        normalizeCoordinate(options.initial_width, 0.7, 1.5),
        normalizeCoordinate(options.initial_weight, 0.5, 1.5),
        normalizeCoordinate(options.initial_phase, -kPi, kPi)};
    const std::vector<double> difference_steps(
        3, options.gradient_step);
    ns_cascade::validateOptimizerPoint(point, difference_steps);

    std::ofstream csv(options.output.c_str());
    if (!csv) throw std::runtime_error("Could not open optimization CSV");
    csv << std::setprecision(17)
        << "iteration,stage,objective,valid,width,secondary_weight,phase,steps,"
           "peak_h_half_ratio,peak_l3_ratio,final_h_half_ratio,final_l3_ratio,"
           "peak_vorticity_ratio,"
           "characteristic_wavenumber_ratio,profile_l1,profile_drift,"
           "peak_cutoff_fraction,max_divergence_defect,max_reality_defect,"
           "gradient_width_coordinate,gradient_weight_coordinate,"
           "gradient_phase_coordinate,gradient_relative_disagreement,"
           "line_step,evaluation_count\n";

    const Evaluation* current = &objective.get(point);
    if (!current->valid) {
        throw std::runtime_error("Initial optimizer point fails the hard gates");
    }
    writeRow(csv, 0, "initial", *current, NULL, 0.0, objective.evaluationCount());

    int accepted_steps = 0;
    for (int iteration = 1; iteration <= options.iterations; ++iteration) {
        const ns_cascade::FiniteDifferenceGradient gradient =
            ns_cascade::checkedCentralDifferenceGradient(
                objective, point, difference_steps);
        if (gradient.maximum_relative_disagreement >
            options.gradient_tolerance) {
            writeRow(csv,
                     iteration,
                     "gradient-rejected",
                     *current,
                     &gradient,
                     0.0,
                     objective.evaluationCount());
            std::cout << "iteration " << iteration
                      << " gradient rejected: refined disagreement="
                      << gradient.maximum_relative_disagreement << '\n';
            break;
        }
        const std::vector<double> direction =
            ns_cascade::normalizedAscentDirection(gradient.extrapolated);
        double line_step = options.initial_line_step;
        bool accepted = false;
        std::vector<double> accepted_point;
        const Evaluation* accepted_value = NULL;
        while (line_step >= options.minimum_line_step) {
            std::vector<double> trial =
                ns_cascade::boundedStep(point, direction, line_step);
            for (std::size_t i = 0; i < trial.size(); ++i) {
                trial[i] = std::max(
                    options.gradient_step,
                    std::min(1.0 - options.gradient_step, trial[i]));
            }
            const Evaluation& trial_value = objective.get(trial);
            const double predicted_ascent =
                gradientDotStep(gradient.extrapolated, point, trial);
            if (trial_value.valid && predicted_ascent > 0.0 &&
                trial_value.objective >
                    current->objective + 1e-4 * predicted_ascent) {
                accepted = true;
                accepted_point = trial;
                accepted_value = &trial_value;
                break;
            }
            line_step *= 0.5;
        }
        if (!accepted) {
            writeRow(csv,
                     iteration,
                     "line-search-rejected",
                     *current,
                     &gradient,
                     0.0,
                     objective.evaluationCount());
            std::cout << "iteration " << iteration
                      << " stopped: no valid improving backtracking step\n";
            break;
        }
        point = accepted_point;
        current = accepted_value;
        ++accepted_steps;
        writeRow(csv,
                 iteration,
                 "accepted",
                 *current,
                 &gradient,
                 line_step,
                 objective.evaluationCount());
        std::cout << std::setprecision(10)
                  << "iteration " << iteration << " objective="
                  << current->objective << " width=" << current->width
                  << " weight=" << current->weight
                  << " phase=" << current->phase
                  << " H1/2=" << current->peak_h_ratio
                  << " k_rms=" << current->characteristic_ratio
                  << " profile-drift=" << current->profile_drift
                  << " cutoff=" << current->peak_cutoff_fraction
                  << " gradient-check="
                  << gradient.maximum_relative_disagreement << '\n';
    }
    if (!csv) throw std::runtime_error("Failed while writing optimization CSV");

    std::cout << std::setprecision(12)
              << "Optimization complete: accepted " << accepted_steps << '/'
              << options.iterations << " steps using "
              << objective.evaluationCount() << " distinct trajectories.\n"
              << "Final objective/H1/2/k_rms/profile-drift/cutoff: "
              << current->objective << '/' << current->peak_h_ratio << '/'
              << current->characteristic_ratio << '/' << current->profile_drift
              << '/' << current->peak_cutoff_fraction << '\n'
              << "This finite-dimensional floating-point optimum is not a "
                 "PDE proof or an adjoint certificate.\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(parseOptions(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
