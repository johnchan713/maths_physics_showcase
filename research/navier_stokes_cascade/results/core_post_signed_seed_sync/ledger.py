"""Post-signed induction-seed compact graph and remaining stage interface."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


first = load("seed_sync_first", HERE.parent/"core_first_sourced_physical_tail_closure/closure.py")
post = load("seed_sync_post", HERE.parent/"core_post_signed_residual_join/bounds.py")


def _order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def normalized_seed_request(n=0, kappa=Q(1, 100000)):
    n = _order(n, "n")
    r = post.mean.rates(kappa)
    trace = post.mean.dependency_trace(
        n, r["A1_integrations"], r["A0_integrations"])
    target = trace["maxima"]["next_signed_target_input"]
    inherited = first.normalized_request(target, kappa)
    return dict(
        output_order=n,
        post_signed_integrations=dict(A0=r["A0_integrations"], A1=r["A1_integrations"]),
        post_signed_raw_sources=dict(
            A0=trace["A0_source_order"], A1=trace["A1_source_order"]),
        first_sourced_target_order=target,
        first_real_ceiling=inherited["largest_first_real_ceiling"],
        inherited_first_sourced_request=inherited,
        trace=trace,
    )


def physical_seed_request(m=0, N=0, h=Q(1, 200),
                          kappa=Q(1, 100000), illustrative_original_pole=1):
    m, N = _order(m, "m"), _order(N, "N")
    request = post.physical_request(m, N, h, kappa)
    target = request["trace"]["maxima"]["next_signed_target_input"]
    inherited = first.normalized_request(target, kappa)
    domain = post.whole_domain(
        request["trace"], previous_wave_pole=30,
        original_pole=Q(illustrative_original_pole))
    return dict(
        output_order=m, requested_q_power=N, h=Q(h),
        current_A0=request["A0"], current_A1=request["A1"],
        first_sourced_target_order=target,
        first_real_ceiling=inherited["largest_first_real_ceiling"],
        inherited_first_sourced_request=inherited,
        post_signed_trace=request["trace"],
        illustrative_domain=domain,
    )


def corrected_seed(kappa=Q(1, 100000)):
    k = Q(kappa)
    if not Q(0) < k < Q(1, 4):
        raise ValueError("use the retained 0<kappa<1/4 regime")
    D = Q(3, 2)-2*k
    C0 = Q(2)-2*k
    old_B0 = Q(3, 2)-6*k
    gain = Q(1, 2)-3*k
    return dict(
        regular_seed=D, mean_seed=C0, old_half_line_minimum=old_B0,
        seed_margin_over_old_minimum=D-old_B0,
        mean_relation=C0-D, finite_stage_gain=gain,
        half_line_rate_recurrence_applicable=D >= old_B0 and C0 == D+Q(1, 2),
    )


def stage_obligation(j=1, m=0, N=0, h=Q(1, 200),
                     kappa=Q(1, 100000)):
    j, m, N = _order(j, "j"), _order(m, "m"), _order(N, "N")
    if j < 1:
        raise ValueError("later-stage obligation starts at j=1")
    h, k = Q(h), Q(kappa)
    if not Q(0) < h < Q(1, 100):
        raise ValueError("use the retained 0<h<1/100 regime")
    seed = corrected_seed(k)
    Bj = seed["regular_seed"]+j*seed["finite_stage_gain"]
    mu = Q(1)-k
    target = (Q(N)+Q(3, 2)+2*h+2*m)/h
    p0 = max(1, (target-mu)//k+1)
    p1 = max(1, (target-Bj)//k+1)
    J0, J1 = m+2*p0+3, m+2*p1+3
    return dict(
        stage=j, regular_order=Bj, radial_source_order=mu,
        target_order=target, integrations=dict(A0=p0, A1=p1),
        raw_sources=dict(A0=J0, A1=J1),
        required_inequalities=(
            f"N_{J0}(Qcirc gr_{j}) <= epsilon^({mu}) P_{{{j},{J0}}}(S+K_{j})",
            f"N_{J1}(Qcirc Ez_postwave_{j}) <= epsilon^({Bj}) Q_{{{j},{J1}}}(S+K_{j})",
        ),
        strict_margins=dict(
            A0=h*(mu+p0*k)-Q(3, 2)-2*h-2*m-N,
            A1=h*(Bj+p1*k)-Q(3, 2)-2*h-2*m-N,
        ),
        stage_indexed_backward_request_recurrence_supplied=False,
        stage_indexed_common_pole_and_support_recurrence_supplied=False,
    )


def interface():
    return dict(
        seed_compact_graph_synchronized_conditionally=True,
        old_half_line_rate_recurrence_compatible_with_corrected_seed=True,
        old_half_line_rate_recurrence_alone_proves_later_compact_inputs=False,
        exact_missing_objects=(
            "stage-indexed backward derivative request recurrence",
            "stage-indexed common pole/operator/support recurrence",
            "the two displayed actual-source inequalities at every j>=1",
        ),
        complete_finite_stage_residual_decay_proved=False,
        infinite_realization_or_blowup_proved=False,
    )
