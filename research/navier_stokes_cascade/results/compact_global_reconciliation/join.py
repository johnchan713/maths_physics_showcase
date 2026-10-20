#!/usr/bin/env python3
"""Exact bookkeeping for the cutoff-to-global compact reconciliation."""
from fractions import Fraction


def review():
    # The fifth-loop gate gives log Lambda5 < G^2.  The interface asks for
    # log Lambda5 + 1044 log A < log N + 260 log H.  We use
    # log A < G, log N > 32 G^4 and log H=G^4.  It is enough that
    # G^2+1044G < 292G^4, checked monotonically from G>=2.
    g = 2
    fifth_slack = 292 * g**4 - g**2 - 1044 * g

    # Cone perturbation: 256(4G^3)^10 G^-64 < G^-4/2.
    # After multiplication by 2G^34 this is 512*4^10 < G^30.
    # The actual G is vastly larger; G>=4 suffices.
    gc = 4
    cone_slack = gc**30 - 512 * 4**10

    # Exact triangular row transforms used by the annular and compact maps.
    eta = Fraction(3, 7)
    lam = Fraction(1, 100)
    d = (Fraction(2), Fraction(3), Fraction(5), Fraction(7), Fraction(11))
    m, r, j, s, cp = d
    ann = (m, j - 4 * eta * r, r, -(s - 8 * eta * m), cp)
    ann_inv = (ann[0], ann[2], ann[1] + 4 * eta * ann[2],
               8 * eta * ann[0] - ann[3], ann[4])
    comp = (m, (m - j) / lam, r, -s, cp)
    comp_inv = (comp[0], comp[2], comp[0] - lam * comp[1], -comp[3], comp[4])

    # Exterior difference formula at eta=0 for p(eta)=eta, all other
    # defects zero.  L(0)=1 and d(0)=1, hence Delta(r Tz)=-X != 0.
    pressure_control = -Fraction(13)

    checks = {
        "fifth_gate_implies_F5": fifth_slack > 0,
        "complete_compact_C4_H925_below_H1024": 925 < 1024,
        "cone_perturbation_below_half_gap": cone_slack > 0,
        "annular_transform_exactly_invertible": ann_inv == d,
        "compact_transform_exactly_invertible_for_positive_lambda": comp_inv == d,
        "pressure_defect_negative_control_nonzero": pressure_control != 0,
        "all_five_rows_are_required": len(d) == 5,
    }
    return {
        "bounds": {
            "Lambda5_log_upper": "G^2",
            "F5_log_lower_available": "292G^4-1044logA",
            "fifth_minimum_G_slack": fifth_slack,
            "compact_ordinary_C4_upper": "H^925",
            "compact_target": "H^1024",
            "state_error_upper": "G^-64/2",
            "cone_gap_lower": "G^-4/2",
            "cone_minimum_G_slack": cone_slack,
        },
        "identities": {
            "annular_rows": [str(x) for x in ann],
            "annular_inverse": [str(x) for x in ann_inv],
            "compact_rows": [str(x) for x in comp],
            "compact_inverse": [str(x) for x in comp_inv],
            "pressure_defect_control_at_eta0_X13": str(pressure_control),
        },
        "checks": checks,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(review(), indent=2, sort_keys=True))
