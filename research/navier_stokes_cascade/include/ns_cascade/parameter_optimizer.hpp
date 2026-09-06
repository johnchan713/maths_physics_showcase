#ifndef NS_CASCADE_PARAMETER_OPTIMIZER_HPP
#define NS_CASCADE_PARAMETER_OPTIMIZER_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ns_cascade {

struct FiniteDifferenceGradient {
    std::vector<double> coarse;
    std::vector<double> refined;
    std::vector<double> extrapolated;
    double maximum_relative_disagreement = 0.0;
};

inline void validateOptimizerPoint(const std::vector<double>& point,
                                   const std::vector<double>& steps) {
    if (point.empty() || point.size() != steps.size()) {
        throw std::invalid_argument(
            "Finite-difference point and steps must have equal non-zero size");
    }
    for (std::size_t i = 0; i < point.size(); ++i) {
        if (!std::isfinite(point[i]) || !std::isfinite(steps[i]) ||
            steps[i] <= 0.0 || point[i] - steps[i] < 0.0 ||
            point[i] + steps[i] > 1.0) {
            throw std::invalid_argument(
                "Finite differences require finite interior normalized coordinates");
        }
    }
}

template <typename Objective>
FiniteDifferenceGradient checkedCentralDifferenceGradient(
    Objective& objective,
    const std::vector<double>& point,
    const std::vector<double>& steps) {
    validateOptimizerPoint(point, steps);
    FiniteDifferenceGradient result;
    result.coarse.resize(point.size());
    result.refined.resize(point.size());
    result.extrapolated.resize(point.size());

    for (std::size_t coordinate = 0; coordinate < point.size(); ++coordinate) {
        std::vector<double> plus = point;
        std::vector<double> minus = point;
        plus[coordinate] += steps[coordinate];
        minus[coordinate] -= steps[coordinate];
        const double coarse_plus = objective(plus);
        const double coarse_minus = objective(minus);
        result.coarse[coordinate] =
            (coarse_plus - coarse_minus) / (2.0 * steps[coordinate]);

        plus[coordinate] = point[coordinate] + 0.5 * steps[coordinate];
        minus[coordinate] = point[coordinate] - 0.5 * steps[coordinate];
        const double refined_plus = objective(plus);
        const double refined_minus = objective(minus);
        result.refined[coordinate] =
            (refined_plus - refined_minus) / steps[coordinate];
        // Cancel the leading O(h^2) central-difference error.
        result.extrapolated[coordinate] =
            (4.0 * result.refined[coordinate] -
             result.coarse[coordinate]) /
            3.0;

        if (!std::isfinite(result.coarse[coordinate]) ||
            !std::isfinite(result.refined[coordinate]) ||
            !std::isfinite(result.extrapolated[coordinate])) {
            throw std::runtime_error(
                "Finite-difference gradient produced a non-finite value");
        }
        const double scale = std::max(
            1e-10,
            std::max(std::abs(result.coarse[coordinate]),
                     std::abs(result.refined[coordinate])));
        result.maximum_relative_disagreement = std::max(
            result.maximum_relative_disagreement,
            std::abs(result.coarse[coordinate] -
                     result.refined[coordinate]) /
                scale);
    }
    return result;
}

inline std::vector<double> normalizedAscentDirection(
    const std::vector<double>& gradient) {
    if (gradient.empty()) {
        throw std::invalid_argument("Cannot normalize an empty gradient");
    }
    double norm_squared = 0.0;
    for (std::size_t i = 0; i < gradient.size(); ++i) {
        if (!std::isfinite(gradient[i])) {
            throw std::invalid_argument("Gradient contains a non-finite value");
        }
        norm_squared += gradient[i] * gradient[i];
    }
    if (norm_squared <= std::numeric_limits<double>::min()) {
        throw std::runtime_error("Gradient is numerically zero");
    }
    const double inverse_norm = 1.0 / std::sqrt(norm_squared);
    std::vector<double> direction(gradient.size());
    for (std::size_t i = 0; i < gradient.size(); ++i) {
        direction[i] = gradient[i] * inverse_norm;
    }
    return direction;
}

inline std::vector<double> boundedStep(const std::vector<double>& point,
                                       const std::vector<double>& direction,
                                       double step_size) {
    if (point.empty() || point.size() != direction.size() ||
        !std::isfinite(step_size) || step_size <= 0.0) {
        throw std::invalid_argument("Bounded step arguments are invalid");
    }
    std::vector<double> trial(point.size());
    for (std::size_t i = 0; i < point.size(); ++i) {
        trial[i] = std::max(
            0.0, std::min(1.0, point[i] + step_size * direction[i]));
    }
    return trial;
}

}  // namespace ns_cascade

#endif
