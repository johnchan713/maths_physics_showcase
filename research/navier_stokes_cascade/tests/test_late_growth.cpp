// Compile the actual trajectory/adjoint/option code in this test translation
// unit. The test must not pass by reproducing the same wiring in a second loop.
#define main ns_state_optimizer_program_main
#include "../src/state_optimize.cpp"
#undef main

#include <functional>

namespace {
using namespace ns_cascade;

void requireLate(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void requireLateThrows(const std::function<void()>& action) {
    bool rejected = false;
    try { action(); } catch (const std::exception&) { rejected = true; }
    requireLate(rejected, "Invalid late-growth data were accepted");
}

Options lateTestOptions(const std::vector<std::string>& extra = {}) {
    std::vector<std::string> text = {"test", "--search-track", "amplification",
        "--growth-objective", "late-rate", "--grid", "8", "--fine-grid", "16",
        "--seed-bandwidth", "2", "--final-time", "0.0041", "--dt", "0.0007",
        "--late-window-start", "0.37", "--late-rate-samples", "4",
        "--profile-path-samples", "3", "--energy", "2"};
    text.insert(text.end(), extra.begin(), extra.end());
    std::vector<char*> arguments;
    for (std::string& argument : text) arguments.push_back(&argument[0]);
    return parseOptions(static_cast<int>(arguments.size()), arguments.data());
}

void testLateAggregation() {
    const LateGrowthSummary result = summarizeLateGrowthRates({0.1, 0.8, -0.2, 0.3}, 0.1);
    double total = 0.0;
    for (double weight : result.gradient_weights) total += weight;
    requireLate(std::abs(total - 1.0) < 1e-14, "Late weights do not sum to one");
    requireLate(result.soft_minimum >= result.minimum &&
                result.soft_minimum <= result.minimum + 0.1 * std::log(4.0),
                "Normalized soft minimum has the wrong bound");
    requireLate(result.gradient_weights[2] > result.gradient_weights[0],
                "Weakest growth did not get the strongest weight");
    const LateGrowthSummary equal = summarizeLateGrowthRates({-2.0, -2.0, -2.0}, 0.1);
    requireLate(equal.soft_minimum == -2.0, "Equal rates are not preserved");
    requireLate(std::abs(equal.gradient_weights[0] - 1.0 / 3.0) < 1e-15,
                "Equal rates have unequal weights");
    const LateGrowthSummary optimistic = summarizeLateGrowthRates({-0.01, 1.0}, 1.0);
    requireLate(optimistic.soft_minimum > 0 && optimistic.minimum < 0,
                "Test no longer demonstrates why the true minimum needs its own gate");
    const LateGrowthSummary extremes = summarizeLateGrowthRates({-1e200, 1e200}, 1e-200);
    requireLate(extremes.gradient_weights[0] == 1.0 &&
                extremes.gradient_weights[1] == 0.0, "Extreme soft minimum is unstable");
    std::vector<double> shifted = result.rates;
    for (double& rate : shifted) rate += 2.0;
    requireLate(std::abs(summarizeLateGrowthRates(shifted, 0.1).soft_minimum -
                         result.soft_minimum - 2.0) < 1e-14,
                "Soft minimum is not translation-equivariant");
    std::vector<double> duplicated = result.rates;
    duplicated.insert(duplicated.end(), result.rates.begin(), result.rates.end());
    requireLate(std::abs(summarizeLateGrowthRates(duplicated, 0.1).soft_minimum -
                         result.soft_minimum) < 1e-14, "Duplicate sampling changed normalization");
    for (std::size_t i = 0; i < result.rates.size(); ++i) {
        std::vector<double> plus = result.rates, minus = result.rates;
        plus[i] += 1e-6;
        minus[i] -= 1e-6;
        const double fd = (summarizeLateGrowthRates(plus, 0.1).soft_minimum -
                           summarizeLateGrowthRates(minus, 0.1).soft_minimum) / 2e-6;
        requireLate(std::abs(fd - result.gradient_weights[i]) < 1e-9,
                    "Soft-minimum scalar gradient failed finite differences");
    }
    requireLateThrows([] { summarizeLateGrowthRates({}, 0.1); });
    requireLateThrows([] { summarizeLateGrowthRates({1.0}, 0.0); });
    requireLateThrows([] { summarizeLateGrowthRates({std::numeric_limits<double>::infinity()}, 0.1); });
}

void testLateOptions() {
    requireLate(lateTestOptions().objective_weights.endpoint_critical_weight == 0.0,
                "Late mode double-counts the endpoint critical reward");
    requireLate(lateTestOptions({"--growth-objective", "endpoint"})
                    .objective_weights.endpoint_critical_weight == 1.0,
                "Historical endpoint reward changed");
    for (const auto& argument : std::vector<std::vector<std::string>>{
             {"--late-window-start", "0"}, {"--late-window-start", "1"},
             {"--late-rate-samples", "1"}, {"--late-rate-samples", "65"},
             {"--late-rate-temperature", "0"}, {"--late-rate-temperature", "-1"},
             {"--growth-objective", "typo"}, {"--search-track", "profile"},
             {"--final-time", "1e-16"},
             {"--late-rate-output", "navier_stokes_optimized_state.csv"},
             {"--growth-objective", "endpoint", "--late-rate-output", "rates.csv"}}) {
        requireLateThrows([&] { lateTestOptions(argument); });
    }
}

void testLateExactShear() {
    Options options = lateTestOptions();
    const PseudospectralSystem system(8, 0.02);
    OptimizationState shear = system.zeroState();
    shear[system.indexOf(WaveVector(2, 0, 0))].y = Complex(1.0, 0.0);
    shear[system.indexOf(WaveVector(-2, 0, 0))].y = Complex(1.0, 0.0);
    const Evaluation evaluation = evaluateTrajectory(system, shear, options, true, true);
    for (double rate : evaluation.late_growth.rates) {
        requireLate(std::abs(rate + 0.08) < 1e-14,
                    "Shear critical rate should be exactly -nu*k^2");
    }
    requireLate(evaluation.late_growth.soft_minimum < 0.0,
                "Decaying shear passed the positive-growth objective");
    requireLateThrows([&] { stateCriticalLogRate(system, system.zeroState()); });
    requireLateThrows([&] { stateCriticalLogRateGradient(system, system.zeroState()); });
}

void testLateFullGradient(int grid, bool aligned) {
    Options options = lateTestOptions(aligned
        ? std::vector<std::string>{"--late-window-start", "0.5", "--late-rate-samples", "3",
                                   "--profile-path-samples", "4"}
        : std::vector<std::string>{});
    const PseudospectralSystem system(grid, options.viscosity);
    const OptimizationState initial = initialState(options, system);
    const Evaluation evaluation = evaluateTrajectory(system, initial, options, true, true);
    requireLate(evaluation.completed && evaluation.fixed_step_safe, "Test trajectory incomplete");
    const OptimizationState gradient = fullInitialGradient(system, initial, evaluation, options);
    const std::vector<double> clock = lateGrowthObservationTimes(options.late_growth, options.final_time);
    for (std::size_t i = 0; i < clock.size(); ++i) {
        double time = 0.0;
        for (int step = 0; step < evaluation.late_growth_state_indices[i]; ++step) {
            time += evaluation.time_steps[static_cast<std::size_t>(step)];
        }
        requireLate(std::abs(time - clock[i]) < 1e-14, "Source inserted at the wrong physical time");
    }
    requireLate(evaluation.late_growth_state_indices.back() == evaluation.steps,
                "Terminal late source was omitted");
    for (int direction_index : {1, 2, 3}) {
        const OptimizationState direction = deterministicOptimizationDirection(
            system, initial, options.seed_bandwidth, direction_index);
        const double predicted = stateRealInnerProduct(gradient, direction);
        double best_error = 1.0;
        for (double epsilon : {1e-3, 1e-4, 1e-5}) {
            const Evaluation plus = evaluateTrajectory(system, stateEnergySphereStep(
                system, initial, direction, epsilon), options, true, false);
            const Evaluation minus = evaluateTrajectory(system, stateEnergySphereStep(
                system, initial, direction, -epsilon), options, true, false);
            requireLate(plus.fixed_step_safe && minus.fixed_step_safe,
                        "Finite differences changed the safe-step branch");
            const double measured = (plus.objective.total - minus.objective.total) / (2 * epsilon);
            const double error = std::abs(predicted - measured) /
                std::max(1e-9, std::max(std::abs(predicted), std::abs(measured)));
            best_error = std::min(best_error, error);
            std::cout << grid << ',' << (aligned ? "aligned" : "split") << ','
                      << direction_index << ',' << epsilon << ',' << predicted << ','
                      << measured << ',' << error << '\n';
        }
        requireLate(best_error < 2e-6, "Full late-window discrete adjoint failed finite differences");
    }
    Options cadence = options;
    cadence.diagnostic_every = 101;
    const Evaluation repeat = evaluateTrajectory(system, initial, cadence, true, true);
    requireLate(repeat.time_steps == evaluation.time_steps &&
                repeat.late_growth.rates == evaluation.late_growth.rates &&
                repeat.objective.total == evaluation.objective.total,
                "Diagnostic cadence changed the late objective or its clock");
    Evaluation incomplete = evaluation;
    incomplete.late_growth_state_indices.pop_back();
    requireLateThrows([&] { fullInitialGradient(system, initial, incomplete, options); });
    incomplete = evaluation;
    incomplete.trajectory.clear();
    incomplete.time_steps.clear();
    requireLateThrows([&] { fullInitialGradient(system, initial, incomplete, options); });
}
}  // namespace

int main() {
    try {
        testLateAggregation();
        testLateOptions();
        testLateExactShear();
        std::cout << std::setprecision(17)
                  << "grid,clock,direction,epsilon,adjoint_slope,finite_difference_slope,relative_error\n";
        testLateFullGradient(8, false);
        testLateFullGradient(16, true);
        std::cout << "Late-growth tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Late-growth test failure: " << error.what() << '\n';
        return 1;
    }
}
