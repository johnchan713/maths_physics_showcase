"""Q213: join the reconstructed primitive bundle and isolate exact core/root premises."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q213-RECONSTRUCTED-PRIMITIVE-BUNDLE-EXACT-ROOT-FRONTIER"
PROVIDERS = (
    "core_current_head_q212_existing_target_degree_raw_map_frontier",
    "core_leading_background_seminorm_reconstruction",
    "core_leading_profile_certificate",
    "core_same_band_gap_reconstruction",
    "core_full_path_ode_reconstruction",
    "core_fixed_profile_projector_reconstruction",
    "core_compact_slow_inverse_reconstruction",
    "core_current_head_q189_same_profile_envelope_gap_audit",
    "core_current_head_q204_global_ledger_uniform_primary",
    "core_all_order_foundation",
    "annular_exterior_foundation",
)
EXPECTED = {
    PROVIDERS[0]: "Q212-existing-target-degree-contract-joined-conditionally-raw-map-frontier-restored",
    PROVIDERS[1]: "fixed-order-background-seminorm-and-real-borel-selector-reconstructed-relative-to-leading-certificate",
    PROVIDERS[2]: "finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses",
    PROVIDERS[3]: "same-band-continuum-gap-contract-reconstructed-on-finalized-leading-profile",
    PROVIDERS[4]: "canonical-full-path-ode-reconstructed-on-current-same-band-tail",
    PROVIDERS[5]: "fixed-profile-same-choice-projector-contract-reconstructed-on-current-full-path-tail",
    PROVIDERS[6]: "compact-slow-weighted-inverse-contract-reconstructed-on-current-same-choice-tail",
    PROVIDERS[7]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
    PROVIDERS[8]: "Q204-global-pressure-coefficient-uniform-primary-and-covariance-ledgers-joined-conditionally",
    PROVIDERS[9]: "local-core-all-order-joint-analytic-reciprocal-and-Cartesian-scale-foundation",
    PROVIDERS[10]: "annular-exterior-all-order-leading-regularity-and-global-moment-foundation-relative-to-exact-roots",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q213 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def finite_bundle(majorants, gaps, inverse_constants, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q213 identity required")

    def values(rows, positive=False):
        out = []
        for row in rows:
            if isinstance(row, bool):
                raise ValueError("finite rational bundle entries required")
            try:
                row = Q(row)
            except (TypeError, ValueError, ZeroDivisionError):
                raise ValueError("finite rational bundle entries required") from None
            if row < 0 or (positive and row <= 0):
                raise ValueError("bundle majorants are nonnegative and gaps are positive")
            out.append(row)
        if not out:
            raise ValueError("nonempty finite bundle rows required")
        return out

    majorants = values(majorants)
    gaps = values(gaps, positive=True)
    inverse_constants = values(inverse_constants)
    return {
        "A_R": max(majorants),
        "g_R": min(gaps),
        "K_R": max(inverse_constants),
        "finite": True,
    }


def reconstructed_bundle(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q213 identity required")
    return {
        "leading_majorant": (
            "the finalized leading-profile certificate and fixed-order background reconstruction "
            "give every requested selected-profile/first-coefficient weighted seminorm before "
            "the real Borel scale is selected"
        ),
        "borel_selector": (
            "after the finite seminorms are known, the real selector chooses a_n so every "
            "fixed row satisfies C_(n,m) q^(nh)<2^(-n)"
        ),
        "gap_floor": (
            "the same-band reconstruction selects one accepted tail with positive value, "
            "spectral, cone, determinant, support and phase-normal floors"
        ),
        "full_path": (
            "the homogeneous/sourced ODE has fixed-jet Duhamel bounds and retains the "
            "unprojected n dot f pressure and every cutoff derivative"
        ),
        "projector_inverse": (
            "one explicit same-choice projector family, both exact moments, shifted compact "
            "maps, zero-Haar temporal inverse and fixed slow response share that tail"
        ),
        "complete_pressure_memory": [
            "Pi1/Pi2", "all six covariances", "radial and axial velocity",
            "whole moving pressure product", "both moment identities",
        ],
        "Q212_raw_majorant_gap_inverse_obligation_closed_conditionally": True,
        "one_common_tail_and_same_profile": True,
        "numerical_high_jet_bundle_evaluated": False,
        "unconditional_foundation_claimed": False,
    }


def exact_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q213 identity required")
    return {
        "core_fixed_point": (
            "for the exact pressure-family datum, prove the infinite Banach-space map is a "
            "self-map and contraction and that its fixed point satisfies the original core equations"
        ),
        "exact_root_rows": (
            "for each selected modulation, annular, axial-strength and heat map prove "
            "F(c_*)=0 exactly and a positive Jacobian/inverse floor on the accepted root neighborhood"
        ),
        "pressure_identification": (
            "identify the pressure-family datum used by the local core with the finalized "
            "forward/backward full-pressure profile and unchanged gauge"
        ),
        "precise_joint_gate": (
            "q_core<1, F_root(c_*)=0, inf_N ||D F_root(c)^(-1)||^(-1)>0, "
            "and Pi_core=Pi_final on the same selected profile"
        ),
        "named_analytical_providers_exist": True,
        "independently_rederived_in_Q213": False,
        "tests_or_hashes_are_continuum_proof": False,
        "unconditional_residual_or_blowup": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q213 identity required")
    status = statuses()
    bundle = reconstructed_bundle()
    frontier = exact_frontier()
    sample = finite_bundle([1, 9, Q(7, 2)], [Q(1, 5), Q(1, 9)], [3, 11])
    checks = {
        "provider_statuses_match": all(status.values()),
        "finite_bundle_quantifier_closes": (
            sample["A_R"] == 9 and sample["g_R"] == Q(1, 9)
            and sample["K_R"] == 11 and sample["finite"]
        ),
        "Q212_bundle_is_already_reconstructed_conditionally": (
            bundle["Q212_raw_majorant_gap_inverse_obligation_closed_conditionally"]
            and bundle["one_common_tail_and_same_profile"]
            and not bundle["numerical_high_jet_bundle_evaluated"]
        ),
        "exact_core_root_frontier_not_promoted": (
            frontier["named_analytical_providers_exist"]
            and not frontier["independently_rederived_in_Q213"]
            and not frontier["tests_or_hashes_are_continuum_proof"]
        ),
    }
    return {
        "status": "Q213-reconstructed-primitive-bundle-joined-conditionally-exact-core-root-frontier-isolated",
        "checks": checks,
        "scope": {
            "raw_majorant_gap_inverse_bundle_joined_conditionally": True,
            "actual_high_jet_constants_numerically_evaluated": False,
            "exact_core_and_root_hypotheses_independently_reproved": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
        "next": (
            "audit the exact core pressure-family fixed point and the selected implicit roots "
            "against the finalized same-pressure profile, beginning with the core pressure datum"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q213 identity required")
    sample = finite_bundle([1, 9, Q(7, 2)], [Q(1, 5), Q(1, 9)], [3, 11])
    return {
        "status": review()["status"],
        "providers": statuses(),
        "finite_bundle_example": {
            "A_R": str(sample["A_R"]), "g_R": str(sample["g_R"]),
            "K_R": str(sample["K_R"]), "diagnostic_only": True,
        },
        "reconstructed_bundle": reconstructed_bundle(),
        "frontier": exact_frontier(),
        "review": review(),
    }
