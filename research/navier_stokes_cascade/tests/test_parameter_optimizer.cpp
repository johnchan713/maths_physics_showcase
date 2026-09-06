#include "ns_cascade/parameter_optimizer.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void expect(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double quadratic(const std::vector<double>& point) {
    const double dx = point[0] - 0.7;
    const double dy = point[1] - 0.2;
    return 2.0 - 3.0 * dx * dx - 0.5 * dy * dy;
}

void testCheckedGradient() {
    const std::vector<double> point = {0.4, 0.6};
    const std::vector<double> steps = {0.04, 0.03};
    const ns_cascade::FiniteDifferenceGradient gradient =
        ns_cascade::checkedCentralDifferenceGradient(
            quadratic, point, steps);
    expect(std::abs(gradient.extrapolated[0] - 1.8) < 1e-12,
           "First checked gradient component is inaccurate");
    expect(std::abs(gradient.extrapolated[1] + 0.4) < 1e-12,
           "Second checked gradient component is inaccurate");
    expect(gradient.maximum_relative_disagreement < 1e-12,
           "Gradient refinement falsely disagreed on a quadratic");
}

void testAscentAndBounds() {
    const std::vector<double> point = {0.4, 0.6};
    const std::vector<double> gradient = {1.8, -0.4};
    const std::vector<double> direction =
        ns_cascade::normalizedAscentDirection(gradient);
    const std::vector<double> trial =
        ns_cascade::boundedStep(point, direction, 0.2);
    expect(quadratic(trial) > quadratic(point),
           "Normalized gradient direction did not ascend");

    const std::vector<double> boundary_trial =
        ns_cascade::boundedStep(
            std::vector<double>{0.99, 0.01},
            std::vector<double>{1.0, -1.0},
            0.2);
    expect(boundary_trial[0] == 1.0 && boundary_trial[1] == 0.0,
           "Bounded step left the normalized parameter box");
}

void testInvalidDifferencePointRejected() {
    bool rejected = false;
    try {
        ns_cascade::checkedCentralDifferenceGradient(
            quadratic,
            std::vector<double>{0.01, 0.5},
            std::vector<double>{0.02, 0.02});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    expect(rejected, "Out-of-bounds central difference was accepted");
}

}  // namespace

int main() {
    try {
        testCheckedGradient();
        testAscentAndBounds();
        testInvalidDifferencePointRejected();
        std::cout << "All parameter-optimizer tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
