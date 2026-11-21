"""Q193: differentiated Duhamel and phase Faà-di-Bruno recurrence audit."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q193-DUHAMEL-PHASE-RECURRENCE"

PROVIDERS = (
    "core_current_head_q192_frame_graph_jet_transfer",
    "core_current_head_q191_physical_pulse_seed_audit",
    "physical_pulse_foundation",
    "actual_band_pulse_bounds",
)

EXPECTED = {
    PROVIDERS[0]: "Q192-auxiliary-frame-jets-closed-and-physical-graph-loss-corrected",
    PROVIDERS[1]: "Q191-physical-pulse-chart-frame-pressure-and-scaling-seed-rederived",
    PROVIDERS[2]: "actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas",
    PROVIDERS[3]: "conditional-actual-band-pulse-primary-covariance-and-linear-bounds",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _order(n):
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("nonnegative integer order required")
    return n


def duhamel_degrees(order, coefficient_powers=None, identity=IDENTITY):
    """Fixed-v homogeneous differentiated-Duhamel polynomial degrees."""
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    n = _order(order)
    powers = list(range(n + 1)) if coefficient_powers is None else list(coefficient_powers)
    if len(powers) != n + 1 or any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in powers):
        raise ValueError("one nonnegative integer coefficient degree per order required")
    answer = [0]
    for j in range(1, n + 1):
        answer.append(max(powers[i] + 1 + answer[j - i] for i in range(1, j + 1)))
    return answer


def sourced_degrees(order, coefficient_powers, source_powers, identity=IDENTITY):
    """Zero-data sourced recurrence; no assertion that an actual source satisfies it."""
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    n = _order(order)
    coefficient_powers, source_powers = list(coefficient_powers), list(source_powers)
    arrays = (coefficient_powers, source_powers)
    if any(len(a) != n + 1 for a in arrays) or any(
        isinstance(x, bool) or not isinstance(x, int) or x < 0 for a in arrays for x in a
    ):
        raise ValueError("complete nonnegative integer degree arrays required")
    forcing = [max(coefficient_powers[i] + source_powers[j - i] for i in range(j + 1))
               for j in range(n + 1)]
    answer = [forcing[0] + 1]
    for j in range(1, n + 1):
        answer.append(max([forcing[j] + 1] + [
            coefficient_powers[i] + 1 + answer[j - i] for i in range(1, j + 1)
        ]))
    return {"forcing": forcing, "solution": answer}


def conservative_amplitude(order=6, identity=IDENTITY):
    """Carry Q192's nonsharp frame powers rather than importing sharp 2j."""
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    n = _order(order)
    sharp = duhamel_degrees(n)
    conservative_coefficients = [0] + [2 * i + 4 for i in range(1, n + 1)]
    conservative = duhamel_degrees(n, conservative_coefficients)
    return {
        "order": n,
        "provider_sharp_degrees": sharp,
        "Q192_conservative_coefficient_degrees": conservative_coefficients,
        "independent_conservative_degrees": conservative,
        "sharp_formula": [2 * j for j in range(n + 1)],
        "conservative_formula": [7 * j for j in range(n + 1)],
        "finite_polynomial_loss_for_every_fixed_order": conservative == [7 * j for j in range(n + 1)],
        "sharp_two_j_claimed_as_independently_needed_here": False,
    }


def phase_partition(spatial_order, time_order, factors, h=Q(1, 200), identity=IDENTITY):
    """Exact Q exponent for one Faà-di-Bruno partition."""
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    a, b, j = _order(spatial_order), _order(time_order), _order(factors)
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    total = a + b
    if total == 0:
        if j != 0:
            raise ValueError("zero derivative has no nonempty phase partition")
    elif not 1 <= j <= total:
        raise ValueError("phase factor count must lie between one and total order")
    raw = Q(a, 2) + b * (1 + h) + j * h / 2
    target = a * (Q(1, 2) + h / 2) + b * (1 + 3 * h / 2)
    return {
        "spatial_order": a,
        "time_order": b,
        "phase_factors": j,
        "raw_loss": raw,
        "target_loss": target,
        "margin": target - raw,
        "expected_margin": h * (total - j) / 2,
        "carrier_bound": "k^j<=2^j*Q^(-h*j/2) for fixed harmonic m",
        "absorbed": raw <= target,
    }


def all_phase_partitions(spatial_order, time_order, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    a, b = _order(spatial_order), _order(time_order)
    if a + b == 0:
        return [phase_partition(0, 0, 0)]
    return [phase_partition(a, b, j) for j in range(1, a + b + 1)]


def uncut_primary_request(spatial_order=2, time_order=1, identity=IDENTITY):
    """Join conservative amplitude and phase recurrences before field assembly."""
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    a, b = _order(spatial_order), _order(time_order)
    amplitude = conservative_amplitude(a + b)
    partitions = all_phase_partitions(a, b)
    return {
        "spatial_order": a,
        "time_order": b,
        "amplitude_polynomial_degree": amplitude["independent_conservative_degrees"][a + b],
        "phase_Q_loss": partitions[-1]["target_loss"],
        "every_phase_partition_absorbed": all(x["absorbed"] for x in partitions),
        "uncut_frame_amplitude_bound": (
            "|D^(a,b)t_m|<=C_ab*S^B_ab*P for each fixed nonzero harmonic"
        ),
        "pressure": (
            "pi=i*(n.Kbg*t-nprime.t)/(k*m*|n|^2) for the homogeneous primary; "
            "the n.f channel remains mandatory for sourced inverses"
        ),
        "no_new_Q_loss_from_fixed_v_Duhamel": True,
        "positive_lower_bound_only_for_uncut_radial_amplitude": True,
    }


def correction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    return {
        "sharp_provider_degree": "2j when coefficient order i costs S^i",
        "independent_degree_used_here": "7j from Q192's conservative S^(2i+4) ledger",
        "reason": (
            "the all-fixed conclusion needs a finite polynomial S loss, so importing the sharper "
            "2j without rederiving its sharper coefficient ledger is unnecessary"
        ),
        "Q_exponent_unchanged": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    return {
        "closed_here": (
            "fixed-order uncut homogeneous amplitude Duhamel recurrence and every phase "
            "Faadi-Bruno Q-exponent partition, relative to Q192's frame inputs"
        ),
        "blocking_inequality": (
            "assemble vector potentials and pressures with the actual transverse, pulse and shell "
            "cutoffs, take the complete cylindrical/Cartesian curl, and prove the stated velocity, "
            "pressure and potential suprema plus smooth closed-support extension on one domain"
        ),
        "complete_cutoff_curl_pressure_potential_assembly_rederived_here": False,
        "generated_source_endpoint_extension_rederived_here": False,
        "provider_supplies_primary_assembly_relative_to_inputs": True,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    status = statuses()
    amplitude = conservative_amplitude(10)
    phase = [phase_partition(a, b, j) for a in range(6) for b in range(5)
             for j in (range(1, a + b + 1) if a + b else (0,))]
    requests = [uncut_primary_request(a, b) for a in range(4) for b in range(3)]
    fix = correction()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "conservative_Duhamel_recurrence_closes_every_fixed_order": amplitude["finite_polynomial_loss_for_every_fixed_order"],
        "every_checked_phase_partition_fits_target_Q_loss": all(x["absorbed"] and x["margin"] == x["expected_margin"] for x in phase),
        "all_checked_uncut_requests_have_no_new_Q_loss": all(x["every_phase_partition_absorbed"] and x["no_new_Q_loss_from_fixed_v_Duhamel"] for x in requests),
        "sharper_degree_not_silently_imported": not amplitude["sharp_two_j_claimed_as_independently_needed_here"],
        "complete_field_assembly_not_promoted": not boundary["complete_cutoff_curl_pressure_potential_assembly_rederived_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "amplitude_recurrence": amplitude,
        "phase_partitions": phase,
        "sample_requests": requests,
        "correction": fix,
        "frontier": boundary,
        "scope": {
            "uncut_homogeneous_all_fixed_amplitude_recurrence_closed_relative_to_inputs": True,
            "all_fixed_phase_partition_Q_ledger_closed": True,
            "complete_primary_field_suprema_rederived": False,
            "complete_residual_improvement_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q193 identity required")
    return {
        "status": "Q193-uncut-amplitude-Duhamel-and-phase-partitions-closed-conditionally",
        "established": (
            "Q192's conservative all-fixed frame ledger gives a finite 7j Duhamel degree and "
            "every Faadi-Bruno phase partition fits the unchanged sx,st Q exponents"
        ),
        "correction": correction()["reason"],
        "next": frontier()["blocking_inequality"],
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
