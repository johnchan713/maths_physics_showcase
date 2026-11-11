"""Q5 operator/source audit on the Q173 fixed-projector object."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q174-Q5-operator-source-audit"
PROVIDERS = (
    "core_current_head_q173_q4_source_audit",
    "core_current_head_q5_weighted_inverse_source_join",
    "core_direct_compact_slow_provider_reconstruction",
    "core_compact_slow_inverse_reconstruction",
    "core_current_head_q18_edge_weight_obstruction",
    "core_current_head_q19_typed_weight_cycle",
    "core_current_head_q20_gaussian_quotient_join",
    "core_current_head_q21_lifted_path_join",
    "core_current_head_q26_source_domain_join",
    "core_current_head_q94_full_regular_source",
    "core_current_head_q117_current_covariant_response",
    "core_current_head_q149_primary_pulse_join",
    "core_current_head_q160_physical_target_induction",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_canonical_full_path_ode_join",
)
EXPECTED = {
    "core_current_head_q173_q4_source_audit": "Q173-Q4-fixed-profile-source-audited-with-generated-source-boundary",
    "core_current_head_q5_weighted_inverse_source_join": "current-head-Q5-weighted-inverse-and-log-source-domain-bound",
    "core_direct_compact_slow_provider_reconstruction": "direct-compact-slow-operators-proved-with-actual-source-supremum-frontier",
    "core_compact_slow_inverse_reconstruction": "compact-slow-weighted-inverse-contract-reconstructed-on-current-same-choice-tail",
    "core_current_head_q18_edge_weight_obstruction": "blanket-full-edge-source-implication-obstructed",
    "core_current_head_q19_typed_weight_cycle": "conditional-typed-radial-cycle-invariant",
    "core_current_head_q20_gaussian_quotient_join": "conditional-whole-path-Gaussian-quotient-jet-join",
    "core_current_head_q21_lifted_path_join": "conditional-affine-whole-path-source-transfer",
    "core_current_head_q26_source_domain_join": "fixed-q-local-source-domain-closure",
    "core_current_head_q94_full_regular_source": "same-state-complete-first-real-regular-source-join",
    "core_current_head_q117_current_covariant_response": "current-source-covariant-response-and-full-restoration",
    "core_current_head_q149_primary_pulse_join": "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
    "core_current_head_q160_physical_target_induction": "Q160-physical-all-finite-cycle-target-moment-induction-joined",
    "core_current_head_q161_source_tail_simultaneous_induction": "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    "core_canonical_full_path_ode_join": "same-tail-canonical-full-path-ode-pressure-and-cutoff-contract-joined",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q174 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def half_weight_obstruction(pole, index):
    """Exact lower exponent for sqrt(zeta)/zeta at x=1/index, w=1."""
    if type(pole) is not int or pole < 0 or type(index) is not int or index < 1:
        raise ValueError("nonnegative pole and positive integer index required")
    lower_log_ratio = Q(index * index, 2) - pole * index
    return {
        "pole": pole,
        "index": index,
        "lower_log_ratio": lower_log_ratio,
        "positive_after_four_poles": index >= max(1, 4 * pole) and lower_log_ratio > 0,
        "blanket_full_weight_follows": False,
    }


def operator_gate(source_membership, exact_moment, pressure_rows):
    """The direct compact/slow estimates are one-way mapping statements."""
    if type(source_membership) is not bool or type(exact_moment) is not bool:
        raise ValueError("boolean mapping premises required")
    pressure_rows = tuple(pressure_rows)
    required = ("Jz_Z", "crho_Z_P", "crho_P_Z", "Pi1", "Pi2", "A1_axial")
    return {
        "source_membership": source_membership,
        "exact_weighted_moment": exact_moment,
        "complete_pressure_memory": all(row in pressure_rows for row in required),
        "operator_output_bound": source_membership and exact_moment and all(row in pressure_rows for row in required),
        "operator_bound_implies_source_membership": False,
    }


def family_join(canonical_family, physical_family, identification_proof=False):
    if not canonical_family or not physical_family or type(identification_proof) is not bool:
        raise ValueError("two named families and a boolean proof flag required")
    same_name = canonical_family == physical_family
    return {
        "canonical_family": canonical_family,
        "physical_family": physical_family,
        "same_name": same_name,
        "identification_proof_supplied": identification_proof,
        "may_transfer_actual_source_bounds": same_name or identification_proof,
    }


def typed_source_gate(harmonic_quotient_bound=None, mean_target_bound=None,
                      family_identified=False):
    if type(family_identified) is not bool:
        raise ValueError("boolean family-identification flag required")
    supplied = harmonic_quotient_bound is not None and mean_target_bound is not None
    if harmonic_quotient_bound is not None and Q(harmonic_quotient_bound) < 0:
        raise ValueError("nonnegative harmonic bound required")
    if mean_target_bound is not None and Q(mean_target_bound) < 0:
        raise ValueError("nonnegative target bound required")
    return {
        "harmonic_quotient_bound_supplied": harmonic_quotient_bound is not None,
        "mean_target_bound_supplied": mean_target_bound is not None,
        "family_identified": family_identified,
        "gate_closed": supplied and family_identified,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q174 identity required")
    status = statuses()
    pressure = ("Jz_Z", "crho_Z_P", "crho_P_Z", "Pi1", "Pi2", "A1_axial")
    conditional_map = operator_gate(True, True, pressure)
    missing_source = operator_gate(False, True, pressure)
    families = family_join("canonical-reference-mixing", "inherited-actual-mixing")
    open_gate = typed_source_gate()
    physical_fixture = typed_source_gate(Q(7, 3), Q(11, 5), True)
    obstruction = half_weight_obstruction(7, 32)
    canonical = evidence("core_canonical_full_path_ode_join")
    q94 = evidence("core_current_head_q94_full_regular_source")
    q117 = evidence("core_current_head_q117_current_covariant_response")
    q149 = evidence("core_current_head_q149_primary_pulse_join")
    canonical_flag = canonical["flags"]["canonical_and_inherited_primary_families_identified"]
    q149_scope = q149["symbolic"]["scope"]
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q173_fixed_projector_object_is_direct": PROVIDERS[0] == "core_current_head_q173_q4_source_audit",
        "compact_slow_maps_close_when_typed_inputs_are_supplied": conditional_map["operator_output_bound"],
        "operator_estimates_do_not_supply_their_input": not missing_source["operator_output_bound"] and not missing_source["operator_bound_implies_source_membership"],
        "half_weight_does_not_promote_to_full_weight": obstruction["positive_after_four_poles"] and not obstruction["blanket_full_weight_follows"],
        "canonical_provider_preserves_nonidentification": canonical_flag is False,
        "physical_join_preserves_nonidentification": "not identified" in q149_scope,
        "Q94_regular_source_polynomial_keeps_native_inequalities_conditional": q94["symbolic"]["native_provider_inequalities_certified"] is False,
        "Q117_covariant_response_keeps_native_entrance_contract_open": q117["symbolic"]["actual_native_entrance_or_high_pole_contracts_certified"] is False,
        "Q161_physical_bounds_cannot_be_substituted_into_canonical_chain": not families["may_transfer_actual_source_bounds"],
        "typed_canonical_gate_remains_open": not open_gate["gate_closed"],
        "finite_same_family_fixture_closes_only_its_fixture": physical_fixture["gate_closed"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "operator_contract": conditional_map,
        "half_weight_control": {k: str(v) if isinstance(v, Q) else v for k, v in obstruction.items()},
        "family_interface": families,
        "typed_open_gate": {
            "harmonic": "sup_{ell>=ell_*} epsilon_ell^{-B_j}(1+|log epsilon_ell|)^{-K_{j,R}} max_{|I|<=R} sup_{lifted Omega_{j,ell}} zeta^{-1/2} delta^{p_{j,I}} |D^I(f_{j,ell}/P_{j,ell})| < infinity",
            "mean_target": "sup_{ell>=ell_*} epsilon_ell^{-C_j}(1+|log epsilon_ell|)^{-L_{j,R}} max_{|I|<=R} sup_{Omega_{j,ell}} zeta^{-1} delta^{q_{j,I}} |D^I Xi_{j,ell}| < infinity",
            "moment_memory": "both exact tangential moments, all six covariances, gd-A1 gd, dZ(Jz+c_rho P), old pressure and Pi1/Pi2",
            "fixture": open_gate,
        },
        "correction": {
            "Q173_pure_power_gate_is_sufficient_not_exact": True,
            "Q5_blanket_full_zeta_generated_source_domain_is_invalid": True,
            "correct_domain_is_typed_and_log_renormalized": True,
            "Q160_Q161_theorem_is_for_inherited_physical_family": True,
            "canonical_Q3_Q4_family_is_not_identified_with_that_family": True,
        },
        "scope": {
            "Q5_operator_mapping_audited_on_Q173_fixed_coefficients": True,
            "actual_canonical_typed_source_membership_proved": False,
            "physical_source_theorem_promoted_to_canonical_source_theorem": False,
            "complete_residual_replayed_here": False,
            "unconditional_pulse_force_exclusion_or_blowup_claimed": False,
        },
        "next": "either prove the two typed bounds for the canonical reference-mixing family, or rebuild the current Q3-Q5 spine explicitly on the inherited physical family before replaying Q6",
    }


def record():
    return {
        "status": "Q174-Q5-operator-source-audited-with-family-and-typed-domain-boundary",
        "established": "the shifted compact maps, zero-Haar temporal inverse and fixed slow response preserve their declared typed weighted domains at every fixed order, conditional on matching source membership and exact moments",
        "correction": "the exact source frontier is typed and log-renormalized; Q5's blanket full-zeta generated-source statement and Q173's pure-power display are not the exact current domain",
        "obstruction": "Q160-Q161 prove source/target bounds for the inherited physical primary family, while the canonical Q3-Q4 provider explicitly leaves that family unidentified; those bounds cannot be substituted into the canonical same-object Q5 chain",
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
