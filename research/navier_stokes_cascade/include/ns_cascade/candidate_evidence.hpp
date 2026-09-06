#ifndef NS_CASCADE_CANDIDATE_EVIDENCE_HPP
#define NS_CASCADE_CANDIDATE_EVIDENCE_HPP

#include "ns_cascade/state_optimizer.hpp"

namespace ns_cascade {

// Keep the evolving system's spectral budget. Only physical-space quadrature
// is evaluated on a common, zero-padded grid. Padding supplies no missing modes.
inline PseudospectralSystem::Diagnostics candidateDiagnostics(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const PseudospectralSystem& sampling_system) {
    if (sampling_system.gridSize() < system.gridSize() ||
        sampling_system.cutoff() < system.cutoff()) {
        throw std::invalid_argument("Candidate sampling grid cannot drop modes");
    }
    PseudospectralSystem::Diagnostics result = system.diagnostics(state);
    if (sampling_system.gridSize() != system.gridSize()) {
        const OptimizationState padded =
            liftOptimizationState(system, state, sampling_system);
        const std::pair<double, double> norms =
            sampling_system.sampledPhysicalNorms(padded);
        result.critical_l3_sample = norms.first;
        result.sampled_vorticity_max = norms.second;
    }
    return result;
}

}  // namespace ns_cascade
#endif
