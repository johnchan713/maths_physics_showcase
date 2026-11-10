"""Audit the Q147--Q157 physical dependency frontier and close its first loose binding."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q158-Q147-Q157-physical-dependency-and-quantifier-frontier"
PROVIDERS = (
    "core_current_head_q157_physical_core_exclusion",
    "core_current_head_q145_hierarchy_gate",
    "compact_jet_envelope",
    "loop_modulation_bounds",
    "core_current_head_q149_primary_pulse_join",
    "core_current_head_q150_initial_covariance_mean_join",
    "core_current_head_q152_physical_stress_feedback_join",
    "core_current_head_q153_physical_source_tail_diagonal_join",
    "core_current_head_q154_physical_endpoint_join",
    "core_current_head_q155_physical_force_energy_join",
    "core_current_head_q156_physical_pressure_comparison",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def physical_route(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q158 identity required")
    return (
        "Q147 exact repair, value cone and five pressure/moment rows",
        "Q148 finite coefficient and reserved-pressure supports",
        "Q149 uniform inherited physical primary family",
        "Q150 complete initial covariance and means",
        "Q152 same-physical-family stress and feedback",
        "Q153 source, tails and complete-potential diagonal",
        "Q154 punctured closed-time endpoint",
        "Q155 localization, compact force and preterminal energy",
        "Q156 fixed-gauge comparison",
        "Q157 protected physical core and conditional exclusion",
    )


def primitive_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q158 identity required")
    return (
        {
            "id": "P1",
            "contract": "selected PREMOD compact value and signed-root envelope",
            "quantifiers": "one fixed finalized profile; supremum on its whole compact modulation hull",
            "status": "discharged-by-Q158-provider-binding",
        },
        {
            "id": "P2",
            "contract": "rule-defined Borel coefficients and physical primary fixed-jet, support and full-path bounds",
            "quantifiers": "for every fixed coefficient order and jet before each Borel scale is selected",
            "status": "conditional-import",
        },
        {
            "id": "P3",
            "contract": "complete covariance and mean cycle with all six entries, axial velocity and full pressure memory",
            "quantifiers": "every fixed physical label and finite cycle state on the same family",
            "status": "conditional-import",
        },
        {
            "id": "P4",
            "contract": "generated-source Gaussian quotient, endpoint zero extension and complete pressure response",
            "quantifiers": "every nonzero generated harmonic and every fixed output jet requested by A0/A1",
            "status": "conditional-import",
        },
        {
            "id": "P5",
            "contract": "compact/slow inverse weighted domains and all temporal-Gaussian, A0/A1 and exterior tails",
            "quantifiers": "every fixed stage j, physical jet m and requested power N",
            "status": "conditional-import",
        },
        {
            "id": "P6",
            "contract": "complete-potential diagonal, punctured endpoint and localized force jets",
            "quantifiers": "one cutoff sequence chosen after finite stage constants; every fixed punctured compact and jet",
            "status": "conditional-import",
        },
        {
            "id": "P7",
            "contract": "preterminal pressure/comparison and protected-core implication",
            "quantifiers": "T<1 fixed, then R to infinity, then tau down to zero",
            "status": "conditional-consequence-of-P2-through-P6",
        },
    )


def first_binding(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q158 identity required")
    return {
        "missing_edge_at_Q145": (
            "Q145 invoked the original whole-profile compact envelope but depended only "
            "on Q144, so its selected-profile value premise was not explicitly frozen"
        ),
        "compact_provider": "compact_jet_envelope: |p2|<=A on the selected PREMOD compact hull",
        "root_provider": "loop_modulation_bounds: |mu|<G and |mu*p2|<G for the same state",
        "hierarchy": "A=C^4096<G=exp(A^32); exp(G)<H=exp(G^4)<N=1+floor(H^32)",
        "choices": "Eexp=V0=H",
        "conclusion": "G>=max(2,|mu|,|p2|,|mu*p2|), Eexp>=exp(G), V0<=N",
        "same_profile_required": True,
        "numerical_materialization_required": False,
        "discharged": True,
    }


def quantifier_order(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q158 identity required")
    return (
        "fix the exact core, finalized profile, physical family, roots, gauges and supports",
        "for each fixed coefficient order and jet prove finite seminorms, then choose that Borel scale",
        "choose one derivative-independent admitted later tail",
        "for each fixed label and output jet prove full-path pulse/source estimates",
        "for each fixed (j,m,N) choose finite compact integration and source-jet requests",
        "choose each complete-potential cutoff only after its finite stage constants are known",
        "take one-sided endpoint limits on each fixed punctured compact",
        "fix T<1, send comparison radius R to infinity, then send tau to zero on the core path",
    )


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q158 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q157-same-physical-core-path-and-competitor-exclusion-conditionally",
        PROVIDERS[1]: "original-hierarchy-discharges-Q144-scalar-gate",
        PROVIDERS[2]: "actual-compact-input-envelope-bounded",
        PROVIDERS[3]: "fixed-phase-loop-and-modulation-errors-bounded",
        PROVIDERS[4]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
        PROVIDERS[5]: "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined",
        PROVIDERS[6]: "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
        PROVIDERS[7]: "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
        PROVIDERS[8]: "Q154-physical-diagonal-punctured-endpoint-bound-conditionally",
        PROVIDERS[9]: "Q155-physical-localization-compact-force-and-energy-bound-conditionally",
        PROVIDERS[10]: "Q156-physical-fixed-gauge-pressure-and-expanding-ball-comparison-conditionally",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError(f"Q158 provider status mismatch: {name}")
    envelope = data[PROVIDERS[2]]["flags"]
    loop = data[PROVIDERS[3]]["flags"]
    if not envelope["actual_compact_jet_envelope_verified"]:
        raise ValueError("Selected compact envelope is not certified")
    if not loop["actual_fixed_phase_loop_first_derivatives_bounded"]:
        raise ValueError("Signed-root provider is not certified")
    return {
        "all_statuses_match": True,
        "selected_compact_envelope_bound": True,
        "same_state_signed_root_bound": True,
        "Q145_missing_provider_edge_repaired": True,
        "later_Q149_Q157_conclusions_remain_conditional": True,
    }


def record():
    return {
        "status": "Q158-physical-dependency-frontier-audited-and-first-loose-provider-bound",
        "correction": (
            "Q147--Q157 form one same-physical-object implication, but Q145's scalar gate "
            "did not freeze the compact-envelope and signed-root providers that justify its "
            "whole-profile value premise"
        ),
        "established": (
            "the selected PREMOD compact envelope and same-state signed-root estimates give "
            "the Q145 scalar gates exactly; the physical route and its quantifier order are "
            "acyclic and contain no request-dependent replacement field"
        ),
        "scope": (
            "one missing provider binding and a dependency/quantifier audit; P2--P6 remain "
            "conditional continuum imports, and tests prove provenance and exact implications only"
        ),
        "next": (
            "independently close P2 on the exact Q149 physical primary family, beginning with "
            "the fixed-order Borel/whole-path pulse seminorms before any scale selection, then "
            "propagate the resulting constants through P3--P6"
        ),
        "first_loose_provider_binding_discharged": True,
        "unconditional_admissible_force_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
