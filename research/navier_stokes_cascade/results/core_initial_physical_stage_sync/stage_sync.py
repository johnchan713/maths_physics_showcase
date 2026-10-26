"""Exact stage ledger for the initial compact and first sourced residuals."""
from fractions import Fraction as Q


def _params(h=Q(1, 200), kappa=Q(1, 100000), m=0, N=0, eta=1):
    h, kappa, eta, N = Q(h), Q(kappa), Q(eta), Q(N)
    if not Q(0) < h < Q(1, 100):
        raise ValueError("the retained regime requires 0<h<1/100")
    if not Q(0) < kappa <= Q(1, 100):
        raise ValueError("the retained regime requires small positive kappa")
    if not isinstance(m, int) or isinstance(m, bool) or m < 0:
        raise ValueError("m must be a nonnegative integer")
    if N < 0 or eta <= 0:
        raise ValueError("N must be nonnegative and eta positive")
    return h, kappa, m, N, eta


def stages(h=Q(1, 200), kappa=Q(1, 100000), m=0, N=0, eta=1):
    h, k, m, N, eta = _params(h, kappa, m, N, eta)
    B, D = Q(1)-k, Q(3, 2)-2*k
    Lambda = Q(3, 2)+2*h
    eR = h*B-Lambda-2*m
    eS = h*D-Lambda-2*m
    threshold = (N+Lambda+2*m+eta)/h
    return dict(
        h=h, kappa=k, physical_jet=m, requested_power=N, reserve=eta,
        normalization=Lambda,
        first_real=dict(name="R", regular_order=B,
                        physical_exponent=eR,
                        compact_state="gr_R,Phi_R"),
        first_sourced=dict(name="S", regular_order=D,
                           physical_exponent=eS,
                           compact_state="gr_S,Phi_S"),
        sourced_gain=D-B, physical_gain=eS-eR,
        sufficient_later_order_threshold=threshold,
        first_real_meets_requested_bound=eR-eta >= N,
        first_sourced_meets_requested_bound=eS-eta >= N,
    )


def compact_difference_identity():
    return dict(
        old="(-A0 Qcirc gr_R,0,+A1 Dz Qcirc Phi_R)",
        new="(-A0 Qcirc gr_S,0,+A1 Dz Qcirc Phi_S)",
        difference="(-A0 Qcirc Delta_gr,0,+A1 Dz Qcirc Delta_Phi)",
        covariance_change="Bil(w_R,v)+Cov(v)",
        pressure_change=(
            "Delta_gp=-D1 Delta_Wrr-Dz Delta_Wrz+Delta_Wtt/R; "
            "Delta_P=int Haar Delta_gp; "
            "Delta_pm=T0(Delta_gp-rho Delta_P)"),
        initial_graph_contains_first_sourced_wave=False,
        initial_compact_bound_transfers_to_sourced_state=False,
    )


def manufactured_change():
    """A finite auxiliary fixture showing both changed compact inputs survive.

    Qcirc is represented by subtracting the average of two auxiliary samples.
    The entries are exact rational values at one fixed radius, with D1 and Dz
    represented by exact manufactured derivative data.
    """
    # Two auxiliary samples of a six-entry covariance change and its derivatives.
    rows = (
        dict(D1Wrr=Q(3), DzWrz=Q(1), Wtt_over_R=Q(5),
             Wzz=Q(7), pm=Q(2), gamma=Q(1), dgamma=Q(2)),
        dict(D1Wrr=Q(1), DzWrz=Q(2), Wtt_over_R=Q(1),
             Wzz=Q(2), pm=Q(1), gamma=Q(1), dgamma=Q(0)),
    )
    gr = tuple(-r["D1Wrr"]-r["DzWrz"]+r["Wtt_over_R"] for r in rows)
    # Delta(gamma^2)=2 gamma dgamma+dgamma^2; the remaining -2 G_R Psi
    # term is set to zero only in this algebraic noncancellation fixture.
    phi = tuple(r["Wzz"]+r["pm"]+2*r["gamma"]*r["dgamma"]+r["dgamma"]**2
                for r in rows)
    def qcirc(pair):
        avg = sum(pair, Q(0))/2
        return tuple(x-avg for x in pair)
    return dict(gr=gr, phi=phi, Qcirc_gr=qcirc(gr), Qcirc_phi=qcirc(phi),
                radial_changed=any(qcirc(gr)), axial_changed=any(qcirc(phi)))


def scope():
    return dict(
        initial_compact_tail_closed_conditionally=True,
        post_sourced_regular_gain_closed_conditionally=True,
        synchronized_post_sourced_compact_DAG_closed_here=False,
        primitive_continuum_contracts_certified_here=False,
        complete_physical_residual_decay_proved=False,
        blowup_proved=False,
    )
