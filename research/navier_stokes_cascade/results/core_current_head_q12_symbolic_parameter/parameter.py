"""Exact parameter propagation; rational fixtures never substitute for the root."""
from fractions import Fraction as Q

KAPPA = Q(1, 100000)
DELTA = Q(1, 2) - 3 * KAPPA
B0 = Q(3, 2) - 2 * KAPPA
MU = 1 - KAPPA
HMAX = Q(1, 100)
ROOT_H = "exp(-8*(exp(64)+10))"
GROUPS = ("time", "viscosity", "pressure", "old_new", "new_old", "new_square")


def order(n):
    if type(n) is not int or n < 0:
        raise ValueError("nonnegative integer required")
    return n


def parameter(h):
    if isinstance(h, (float, bool)):
        raise ValueError("exact rational fixture required, not a floating approximation")
    h = Q(h)
    if not 0 < h < HMAX:
        raise ValueError("0<h<1/100 required")
    return h


def ceil(x):
    x = Q(x)
    return -(-x.numerator // x.denominator)


def stage_coefficient(m, target):
    return 2 * (order(target) + 17 + 4 * order(m)) / DELTA


def state_ceiling(stage, n=0):
    j, n = order(stage), order(n)
    return n + 500038 + 99994*j*j + 400035*j


def compact_request(h, stage, m, target):
    """No stage list: propagate both extra integrations to the original primitive."""
    h = parameter(h)
    j, m, target = order(stage), order(m), order(target)
    threshold = (target + Q(3, 2) + 2*m) / h + 2
    rates = {"A0": MU, "A1": B0 + DELTA*j}
    counts = {name: max(1, (threshold-rate)//KAPPA + 1)
              for name, rate in rates.items()}
    raw = {name: m + 2*p + 3 for name, p in counts.items()}
    # Generic trace maximum after p0,p1 are each enlarged by one.
    previous_sourced = m + max(2*counts["A0"] + 19, 2*counts["A1"] + 9)
    if j:
        physical_ceiling = state_ceiling(
            j-1, previous_sourced + 200018 + 99994*(j-1))
    else:
        # First-sourced normalized_request(n) = n + 300021.
        physical_ceiling = previous_sourced + 300021
    return {
        "integrations_before_extra": counts,
        "raw_before_extra": raw,
        "raw_after_extra": {name: n+2 for name, n in raw.items()},
        "base_strict_margins": {
            name: h*(rate + counts[name]*KAPPA - threshold)
            for name, rate in rates.items()},
        "state_primitive_ceiling": state_ceiling(j, m),
        "enlarged_physical_primitive_ceiling": physical_ceiling,
        "combined_primitive_ceiling": max(state_ceiling(j, m), physical_ceiling),
    }


def rational_plan(h, m=0, target=0):
    """Exact diagnostics for explicit rational h; no default selected parameter."""
    h = parameter(h)
    K = stage_coefficient(m, target)
    J0 = ceil(K/h)
    J = J0 + 1
    regular = h*(B0 + DELTA*J) - (Q(5, 2)+2*h+2*m) - h*DELTA/2
    gain = h*DELTA*(J+1)
    tail = gain/2 - (7+2*h+2*m)
    cross = tail - (8+2*h+2*m)
    groups = dict(zip(GROUPS, (tail, tail, tail, cross, cross, 2*tail)))
    return {
        "h": h, "jet": m, "target": target, "K": K, "base_stage": J0, "stage": J,
        "rounding_remainder": h*J0-K,
        "request": compact_request(h, J, m, target),
        "fixed_state_exponent": regular,
        "groups": groups, "complete_minimum": min(regular, *groups.values()),
        "reserves": {"A0": h*KAPPA/2, "A1": h*KAPPA/2,
                     "regular": h*DELTA/2, "cutoff": gain/2},
        "jacobian_floor": 1-2*h,
    }


def universal_bounds(m, target):
    """Lower bounds using r=h*J0-K>=0 and 0<h<1/100; no evaluation of h."""
    K = stage_coefficient(m, target)
    cross = DELTA*K/2 - 15 - 4*m + min(Q(0), DELTA-4)*HMAX
    regular = DELTA*K - Q(5, 2)-2*m + min(Q(0), B0+DELTA/2-2)*HMAX
    a1_dominance = DELTA*K-target-Q(3, 2)-2*m + min(Q(0), B0+DELTA-2)*HMAX
    return {"cross_lower": cross, "regular_lower": regular,
            "A1_rate_minus_threshold_times_h_lower": a1_dominance,
            "complete_lower": min(cross, regular),
            "jacobian_lower": 1-2*HMAX}


def cutoff_requirement(stage, log_constant, degree, h):
    """A necessary input to this sufficient majorant, not a new cutoff selection."""
    j, P = order(stage), order(degree)
    h, c = parameter(h), Q(log_constant)
    if j == 0 or c < 0:
        raise ValueError("positive stage and nonnegative log constant required")
    gamma = h*DELTA*j/2
    return 2*(c+j)/gamma + Q(P*P)/(gamma*gamma)


def symbolic_root_plan(m=0, target=0):
    bounds = universal_bounds(m, target)
    assert bounds["complete_lower"] > target+1
    assert bounds["A1_rate_minus_threshold_times_h_lower"] > 0
    return {
        "h": ROOT_H, "h_domain": "0<h<1/513<1/100", "kappa": str(KAPPA),
        "K": str(stage_coefficient(m, target)), "jet": m, "target": target,
        "J0": "ceil(K/h)", "J": "J0+1", "rounding": "0<=h*J0-K<h",
        "A0_integrations": "floor(((N+3/2+2*m)/h+2-(1-kappa))/kappa)+1",
        "A1_integrations": "1",
        "extra_integration_raw_orders": {"A0": "m+2*p0+5", "A1": "m+7"},
        "state_primitive_ceiling": "m+500038+99994*J^2+400035*J",
        "enlarged_physical_primitive_ceiling":
            "state_ceiling(J-1,m+2*p0+19+200018+99994*(J-1))",
        "reserve_owners": {"A0": "h*kappa/2", "A1": "h*kappa/2",
                           "regular": "h*delta/2", "cutoff": "h*delta*(J+1)/2"},
        "bounds": {key: str(value) for key, value in bounds.items()},
        "cutoff_condition":
            "s_j>max(s_(j-1)+1,log(qstar^-1)+1,1,4*(c_j+ j)/(h*delta*j)+4*P_j^2/(h^2*delta^2*j^2))",
        "actual_cutoff_condition_verified": False,
        "same_field_join_restored": False,
        "endpoint_bound_proved": False,
    }
