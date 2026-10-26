"""Finite-stage compact request, pole and physical-order recurrence."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


seed = load("finite_compact_seed", HERE.parent/"core_post_signed_seed_sync/ledger.py")
generic = load("finite_compact_generic", HERE.parent/"core_post_signed_mean_closure/budgets.py")
sourced = load("finite_compact_sourced", HERE.parent/"core_next_sourced_mean_tail_join/budgets.py")


def order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def parameters(kappa=Q(1, 100000)):
    k = Q(kappa)
    if not Q(0) < k < Q(1, 6):
        raise ValueError("use the retained 0<kappa<1/6 induction regime")
    return dict(
        kappa=k, mu=Q(1)-k, B0=Q(3, 2)-2*k,
        gain=Q(1, 2)-3*k)


def strict_integrations(source, target, kappa):
    source, target, kappa = Q(source), Q(target), Q(kappa)
    if kappa <= 0:
        raise ValueError("positive kappa required")
    return max(1, (target-source)//kappa+1)


def rates(j, kappa=Q(1, 100000)):
    j = order(j, "j")
    p = parameters(kappa)
    B = p["B0"]+j*p["gain"]
    C = B+Q(1, 2)
    return dict(stage=j, regular=B, mean=C, **p)


def layer(j, kappa=Q(1, 100000)):
    """One sourced then generic-signed cycle, traced backwards.

    The derivative offsets come from the actual request_trace and
    dependency_trace functions.  Constants/harmonics may grow with j, but
    these integer derivative offsets remain finite at every fixed stage.
    """
    j = order(j, "j")
    current, nxt = rates(j, kappa), rates(j+1, kappa)
    k, mu = current["kappa"], current["mu"]
    p0s = strict_integrations(mu, current["mean"], k)
    p1s = strict_integrations(nxt["regular"], current["mean"], k)
    p0g = strict_integrations(mu, nxt["mean"], k)
    p1g = strict_integrations(nxt["regular"], nxt["mean"], k)
    st = sourced.request_trace(0, p1s, p0s)
    gt = generic.dependency_trace(0, p1g, p0g)
    so = st["maxima"]["original_inverse_normal"]
    go = gt["maxima"]["next_signed_target_input"]
    return dict(
        from_stage=j, to_stage=j+1,
        sourced_integrations=dict(A0=p0s, A1=p1s),
        generic_integrations=dict(A0=p0g, A1=p1g),
        sourced_backward_offset=so,
        generic_backward_offset=go,
        full_cycle_backward_offset=so+go,
        sourced_trace=st, generic_trace=gt,
    )


def state_ceiling(j, n=0, kappa=Q(1, 100000)):
    """Safe first-real primitive order for stage-j state output order n."""
    j, n = order(j, "j"), order(n, "n")
    base = seed.normalized_seed_request(n, kappa)["first_real_ceiling"]
    layers = [layer(i, kappa) for i in range(j)]
    ceiling = base+sum(x["full_cycle_backward_offset"] for x in layers)
    return dict(stage=j, output_order=n, seed_ceiling=base,
                layers=layers, first_real_ceiling=ceiling,
                stage_decreases_on_every_recursive_edge=True)


def physical_request(j, m=0, N=0, h=Q(1, 200),
                     kappa=Q(1, 100000)):
    """Same-state compact physical tails at corrected post-signed stage j."""
    j, m, N = order(j, "j"), order(m, "m"), order(N, "N")
    h, p = Q(h), parameters(kappa)
    if not Q(0) < h < Q(1, 100):
        raise ValueError("use the retained 0<h<1/100 regime")
    r = rates(j, p["kappa"])
    target = (Q(N)+Q(3, 2)+2*h+2*m)/h
    p0 = strict_integrations(p["mu"], target, p["kappa"])
    p1 = strict_integrations(r["regular"], target, p["kappa"])
    trace = generic.dependency_trace(m, p1, p0)
    previous_sourced_order = trace["maxima"]["next_signed_target_input"]
    if j == 0:
        ceiling = seed.first.normalized_request(
            previous_sourced_order, p["kappa"])["largest_first_real_ceiling"]
        path = dict(first_sourced_target=previous_sourced_order)
    else:
        previous_layer = layer(j-1, p["kappa"])
        previous_state_order = (previous_sourced_order
                                + previous_layer["sourced_backward_offset"])
        prior = state_ceiling(j-1, previous_state_order, p["kappa"])
        ceiling = prior["first_real_ceiling"]
        path = dict(previous_sourced_order=previous_sourced_order,
                    previous_state_order=previous_state_order,
                    previous_state_request=prior)
    loss = Q(3, 2)+2*h+2*m
    return dict(
        stage=j, physical_jet=m, requested_q_power=N,
        source_orders=dict(A0=p["mu"], A1=r["regular"]),
        integrations=dict(A0=p0, A1=p1),
        raw_sources=dict(A0=m+2*p0+3, A1=m+2*p1+3),
        strict_margins=dict(A0=h*(p["mu"]+p0*p["kappa"])-loss-N,
                            A1=h*(r["regular"]+p1*p["kappa"])-loss-N),
        first_real_ceiling=ceiling, backward_path=path,
    )


def pole_support_ledger(j, m=0, N=0, h=Q(1, 200),
                        kappa=Q(1, 100000), original_pole=1,
                        retained_wave_pole=30):
    """Finite distinct-pole archive; no monotonicity of operator constants."""
    j, original_pole = order(j, "j"), Q(original_pole)
    if original_pole < 0:
        raise ValueError("nonnegative original pole required")
    physical = physical_request(j, m, N, h, kappa)
    requests = [physical["first_real_ceiling"]]
    requests += [state_ceiling(i, 0, kappa)["first_real_ceiling"]
                 for i in range(j+1)]
    q_values = []
    retained = Q(retained_wave_pole)
    for L in requests:
        q = max(retained, 6*original_pole, 6*L+9,
                6*L+12*original_pole+6)
        q_values.append(q)
        retained = q
    return dict(
        stage=j, derivative_requests=requests,
        distinct_common_poles=q_values,
        operator_poles=[dict(KT1=6*q, KA1=6*q, Kp=(6*q, 12*q),
                             Kmom=(6*q, 12*q), oriented=6*q)
                        for q in q_values],
        common_argument=(
            "S + sum of every archived K_6q_12q and K_Ie_6q; "
            "no distinct-pole monotonicity used"),
        harmonics_at_most=4*2**j,
        same_closed_supports_preserved_if_inherited_support_contracts_hold=True,
        finite_for_fixed_stage=all(isinstance(x, int) for x in requests),
        actual_constants_or_high_jets_materialized=False,
    )


def first_stage_for_power(N=0, m=0, h=Q(1, 200),
                          kappa=Q(1, 100000), reserve=1):
    N, m = order(N, "N"), order(m, "m")
    h, reserve = Q(h), Q(reserve)
    if not Q(0) < h < Q(1, 100) or reserve <= 0:
        raise ValueError("retained h and positive polynomial reserve required")
    p = parameters(kappa)
    loss = Q(3, 2)+2*h+2*m+reserve
    needed = (Q(N)+loss-h*p["B0"])/(h*p["gain"])
    j = max(0, -(-needed.numerator//needed.denominator))
    return dict(stage=j, preceding_stage=j-1 if j else None,
                achieved_power=h*rates(j, p["kappa"])["regular"]-loss,
                preceding_power=(h*rates(j-1, p["kappa"])["regular"]-loss
                                 if j else None),
                loss=loss, reserve=reserve)


def closure_statement():
    return dict(
        finite_stage_compact_graph_synchronized_conditionally=True,
        all_fixed_order_and_pole_primitive_contracts_required=True,
        stage_constants_may_grow_without_uniform_bound=True,
        actual_large_derivatives_factorials_and_operator_constants_materialized=False,
        fixed_stage_complete_residual_powers_follow_conditionally=True,
        one_field_all_targets_diagonal_proved=False,
        localization_forcing_energy_or_blowup_proved=False,
    )
