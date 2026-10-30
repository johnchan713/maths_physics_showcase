"""Finite exact stress majorants; input seminorms are analytical premises."""
from fractions import Fraction as Q
from math import comb, factorial


def order(m):
    if type(m) is not int or m < 0:
        raise ValueError('Nonnegative integer jet order required')


def nonnegative(x):
    if type(x) not in (int, Q):
        raise TypeError('Exact integer/rational diagnostic bound required')
    x = Q(x)
    if x < 0:
        raise ValueError('Nonnegative seminorm required')
    return x


def stirling_abs(n):
    order(n)
    row = [1]
    for j in range(1, n + 1):
        nxt = [0] * (j + 1)
        for k, v in enumerate(row):
            nxt[k] += (j - 1) * v
            nxt[k + 1] += v
        row = nxt
    return row


def weight_jets(m, a, xlo):
    """Raw derivative constants for BOTH zeta and its reciprocal, factored."""
    order(m)
    a, xlo = nonnegative(a), nonnegative(xlo)
    if a <= 0 or xlo <= 0:
        raise ValueError('Original positive flat constant and annular floor required')
    log = [Q(0)]
    edge = [Q(1)]
    for j in range(1, m + 1):
        log.append((a + 4) * xlo**(-j) * sum(
            c * factorial(k + 1) for k, c in enumerate(stirling_abs(j))))
        edge.append(sum(Q(comb(j-1, i-1)) * log[i] * edge[j-i]
                        for i in range(1, j+1)))
    return log, edge


def compact_tail(m, head):
    """n>=2 only; the same diagonal bounds both quotient families."""
    order(m)
    start = max(2, m)
    if set(head) != set(range(2, start)):
        raise ValueError('Exactly the compact finite head 2<=n<max(2,m) required')
    return sum((nonnegative(c) for c in head.values()), Q(0)) + Q(1, 2**(start-1))


def stress_budget(m, h, a, xlo, leading, first_weighted, first_relative,
                  cutoff, compact_constants, derivative_basis='euler'):
    """All callbacks bound raw through-order profile jets, never sampled values."""
    order(m)
    h = nonnegative(h)
    if not 0 < h < Q(1, 100) or derivative_basis != 'euler':
        raise ValueError('Fixed positive construction h and Euler derivative basis required')
    _, edge = weight_jets(m, a, xlo)
    tail, target, correction, relative = [], [], [], []
    for j in range(m+1):
        K = nonnegative(cutoff(j))
        if K < 1:
            raise ValueError('Cutoff Euler bounds must include its value one')
        tail.append(compact_tail(j, {n: compact_constants(n, j)
                                    for n in range(2, max(2, j))}))
        target.append(sum(Q(comb(j, i))*edge[i]*nonnegative(leading(j-i))
                          for i in range(j+1)))
        correction.append(K*(1+2*h)**j*nonnegative(first_weighted(j)) +
                          sum(Q(comb(j, i))*edge[i]*tail[j-i] for i in range(j+1)))
        relative.append(K*(1+2*h)**j*nonnegative(first_relative(j)) + tail[j])
    return dict(target=max(target), correction=max(correction), relative=max(relative),
                pole=3*m+3, compact_tail=tail, weight=edge)


def fixture(m):
    return stress_budget(m, Q(1, 1000), Q(1, 16), 1,
                         leading=lambda j: j+2, first_weighted=lambda j: 2*j+3,
                         first_relative=lambda j: j+4, cutoff=lambda j: j+1,
                         compact_constants=lambda n, j: n+j+1)


def record():
    f = fixture(4)
    return dict(result='finite-head-weighted-stress-input-budget',
                diagnostic={k: str(v) for k, v in f.items() if k not in ('compact_tail','weight')},
                compact_tail=[str(v) for v in f['compact_tail']],
                original_first_coefficient_kept_out_of_compact_quotient_norm=True,
                both_weight_constants_retained=True, same_schedule_required=True,
                numeric_root_h_evaluated=False, pressure_replaced_by_stress=False,
                full_source_or_residual_proved=False)
