"""Exact finite request graph for the corrected first-real A1 final-row tail."""
from fractions import Fraction as Q


def _order(n):
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError("A nonnegative integer output order is required")
    return n


def _kappa(kappa):
    k = Q(kappa)
    if not Q(0) < k <= Q(1, 100):
        raise ValueError("Use the retained positive small-kappa regime")
    return k


def request_graph(n=0, kappa=Q(1, 100000)):
    """A sufficient same-state graph, without materializing any giant norm.

    The +8 state padding is the conservative pressure-field dependency of the
    already proved axial-potential normal form.  The +3/+4 paddings are the
    existing same-choice signed target and fixed-profile interfaces.  The
    InitialHaar projected moment reads initial.output(n+1), so its lower
    SharpInitial ledger must have one additional slot.
    """
    n = _order(n)
    k = _kappa(kappa)
    source = Q(2)
    nu = Q(3) - 6 * k
    p = (nu - source) // k + 1
    dz_phi = n + 2 * p + 3
    phi = dz_phi + 1
    state = phi + 8
    signed_target = state + 3
    fixed_profile = state + 4
    initial_haar = signed_target
    sharp_initial = initial_haar + 1
    initial_fields = initial_haar + 3
    canonical_covariance = initial_haar + 13
    achieved = source + p * k
    return dict(
        output_order=n,
        kappa=k,
        source_order=source,
        required_order=nu,
        integrations=p,
        achieved_order=achieved,
        strict_margin=achieved - nu,
        raw_Dz_Phi_order=dz_phi,
        raw_Phi_order=phi,
        conservative_same_state_order=state,
        signed_target_order=signed_target,
        fixed_profile_order=fixed_profile,
        initial_Haar_order=initial_haar,
        sharp_initial_order=sharp_initial,
        initial_field_order=initial_fields,
        canonical_covariance_safe_order=canonical_covariance,
        largest_immediate_order=canonical_covariance,
        old_combined_A0_A1_target_order_at_output_zero=400009,
        A0_graph_deferred_to_complete_residual=True,
        support_trace_enumeration_required=False,
        actual_band_uniform_field_inequalities_proved=False,
    )


def pole_domains(n=0, old_pole=1, target_allowance=12, retained_wave_pole=30):
    """Illustrative sufficient common pole using the corrected A1-only state."""
    for value in (old_pole, target_allowance, retained_wave_pole):
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError("Nonnegative integer pole data are required")
    g = request_graph(n)
    N = g["conservative_same_state_order"]
    q = max(retained_wave_pole, 6 * old_pole, 6 * N + 9,
            6 * N + target_allowance + 6)
    return dict(q=q, KT1=2*q, KA1=2*q, Kp=(2*q, 4*q), Kmom=(2*q, 4*q),
                state_order=N, illustrative_only=True)


def unresolved_inequalities(n=0):
    g = request_graph(n)
    return dict(
        Phi=("N_%d(Qcirc Phi_postsigned) <= epsilon P_%d(S)"
             % (g["raw_Phi_order"], g["raw_Phi_order"])),
        DzPhi=("N_%d(Dz Qcirc Phi_postsigned) <= epsilon^2 Q_%d(S)"
               % (g["raw_Dz_Phi_order"], g["raw_Dz_Phi_order"])),
        weighted_prefix=("all actual weighted initial/root/mean inputs through order %d "
                         "obey their recorded epsilon, zeta and edge-pole bounds"
                         % g["largest_immediate_order"]),
        same_root_support=("the identical positive root, whole compact hull, endpoint traces, "
                           "moving rho-squared projector and pressure memory are retained"),
        passing_tests_are_continuum_proof=False,
    )


def qualitative_closure(n=0):
    g = request_graph(n)
    return dict(
        every_requested_index_is_finite=all(isinstance(v, int) for key, v in g.items()
            if key.endswith("order") and key not in ("source_order", "required_order", "achieved_order")),
        fixed_C_infinity_profiles_have_finite_requested_seminorms=True,
        rho_squared_projector_has_band_uniform_fixed_profile_recipe=True,
        finite_support_grammar_gives_smooth_zero_extension_and_matching_traces=True,
        numerical_enumeration_of_all_traces_or_factorials_required=False,
        qualitative_smoothness_implies_weighted_field_rate=False,
    )
