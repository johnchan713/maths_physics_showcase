"""Exact bounds for the common-support and row comparison; no sampled field."""
import json
from fractions import Fraction as F
from pathlib import Path


def union_area(labels, radius):
    if not isinstance(labels, int) or labels < 0 or radius < 0:
        raise ValueError("invalid support bounds")
    return min(1, labels * (8 * F(radius))**2)


def normalized_floor(mean_floor, curl_loss, inverse_cost, row_loss):
    """Losses already divided by epsilon and epsilon^(3-5*kappa), respectively."""
    if min(mean_floor, curl_loss, row_loss) < 0 or inverse_cost <= 0:
        raise ValueError("invalid floor, loss or inverse cost")
    return max(0, F(mean_floor-curl_loss)/inverse_cost-row_loss)


def record():
    src = Path(__file__).resolve().parent.parent/'auxiliary_geometry_audit/protocol.json'
    geom = json.loads(src.read_text())
    area = union_area(geom['original_pointwise_slow_overlap_bound'], F(geom['radius_upper_bound']))
    assert area < 1
    from angular import ledger
    rates = ledger()
    k = F(1, 100000)
    gap = rates['final_error'] - rates['final_leading']
    assert gap == k and rates['final_leading'] == 3-5*k
    return dict(union_area_upper=str(area), open_complement_area_lower=str(1-area),
                curl_relative_power=str(2-k), row_relative_power=str(gap),
                angular_exponent=str(rates['final_leading']), auxiliary_derivative_order=8,
                applies_to_physical_phase_graph=False)
