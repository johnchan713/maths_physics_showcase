"""Original projected axial tail; finite norm programs, not field certificates."""
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'q48_full', HERE.parent / 'core_current_head_q47_full_pressure_correction/full_rows.py')
full = importlib.util.module_from_spec(spec)
spec.loader.exec_module(full)
moving = full.moving
radial = moving.radial
order, EP = full.order, full.EP


def ledger(n, p=99999):
    order(n)
    order(p)
    if p < 1:
        raise ValueError('At least one radial integration required')
    J = n + 2*p + 1
    return dict(output=n, integrations=p, raw_source=J, covariance=J+3,
                pressure=J+1, pressure_shape=J+2, axial_shape=n+1,
                radial_multiplier=n+p, source_pole=12*(J+3)+18, output_pole=0)


def kt(j, pole, jets, shell):
    """KT0 only: the authentic Q42 constants and Q43 conjugation, no time map."""
    g = jets.geo
    data = dict(shell.remaining(jets)(j, pole)['bounds'], **g.multipliers(j, 0))
    raw = radial.provider.constants(j, data)['KT']
    return (radial.chart.transfer(j, 1/g.d, g.Rmax, g.Umin, True) * raw
            * radial.chart.transfer(j, g.d, g.Umax, g.Rmin, False))


def initial_shell(h):
    """Invert the rational witnesses stored by the ORIGINAL budget's geometry()."""
    left, right = h['left_log_gap_lower'], h['right_log_gap_lower']
    l, u = 2*h['Rminus'], h['Rplus']/2
    return (l*l*(2-left)/(2*(2+left)), u*u*(2+right)/(2*(2-right)))


def coefficient(current, n, p, *, d, pressure_cutoff, axial_cutoff,
                covariance_coverage, pole_coverage):
    """Each cutoff = (left witness, right witness, normalized physical shape jets).
    Coverage is a declaration for the SAME original field; it is not its proof.
    Keep the enormous floor polynomial factored. Only small orders are replayed.
    """
    if not isinstance(current, full.refresh.maps.NewPoleMaps):
        raise TypeError('Original Q38 state required')
    L = ledger(n, p)
    j, N, v = L['raw_source'], L['covariance'], L['source_pole']
    order(covariance_coverage)
    order(pole_coverage)
    if covariance_coverage < N or pole_coverage < v:
        raise ValueError('Original covariance archive does not cover the tail')
    b = current.target.source.haar.initial.budget
    chart = b.primary.chart
    width = b.primary.primary.root.width
    Xa, Xb = initial_shell(b.hull)

    def cutoff(args, needed):
        c0, c1, shape = args
        if len(shape) <= needed:
            raise ValueError('Original fixed cutoff jet coverage too short')
        return (moving.Moving(d, chart.Rmin, chart.Rmax, Xb, width, shape),
                moving.directed.Shell(Xa, Xb, c0, c1, width))

    pj, ps = cutoff(pressure_cutoff, L['pressure_shape'])
    aj, ax = cutoff(axial_cutoff, L['axial_shape'])
    geo = aj.geo
    W = full.convert(b.primary.covariance_R(N)['W'])
    r = sum(chart.Rmin**(-k-1) for k in range(N+1))
    g = 2 + 6*chart.B*sum((k+1)*chart.Rmin**(-k) for k in range(N+1))
    D = j+4
    gr = (D*g+D+2*r)*W
    pm = kt(j+1, v, pj, ps)*gr
    ez = (D*g+r)*W + D*(W+pm)

    forward = radial.chart.transfer(j, 1/geo.d, geo.Rmax, geo.Umin, True)
    reverse = radial.chart.transfer(n, geo.d, geo.Umax, geo.Rmin, False)
    B = 2**n * geo.multipliers(n, 1)['k'] * aj(n, v)['bounds']['chi_U']
    A = geo.power(n+p, 2-geo.d, 1/geo.d)
    # Only the n+p non-torus derivatives can land on the scalar radial power.
    kernel = (8 * 2**(p+1) * 12**p * 2**(2*n+p) * (geo.Umax-geo.Umin) * B*A)
    base = (2 * kernel * forward * reverse * ax.bounds(v)['edge']
            * ax.bounds(0)['interior'] * ez)
    floor = 1 + full.convert(b.compact(0, b.pole)['I'])*EP({1: 1})
    return dict(ledger=L, projected_source=2*ez, full_pressure=pm,
                coefficient=base, floor_base=floor, floor_power=p)


def record():
    return dict(status='conditional-original-axial-tail', row_zero=ledger(0),
                q41_output_zero=ledger(10), materialized_physical_high_jets=False,
                uniform_actual_pulse_or_complete_residual=False)
