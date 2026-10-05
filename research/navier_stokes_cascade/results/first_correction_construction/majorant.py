"""Convergent Picard majorant from the sparse angular-derivative block."""
import mpmath as mp


def allowed_word(word):
    """One denotes A1 d_eta, zero denotes A0; adjacent ones vanish."""
    if any(v not in (0, 1) for v in word):
        raise ValueError('Only multiplication or derivative factors allowed')
    return all(a+b < 2 for a, b in zip(word, word[1:]))


def log_term(k, operator_bound, radius, strip_loss):
    C, a, delta = map(mp.mpf, (operator_bound, radius, strip_loss))
    if type(k) is not int or k < 0 or C < 1 or a <= 0 or not 0 < delta <= 1:
        raise ValueError('Invalid Picard bound or analytic strip loss')
    count = (k+1)//2
    return ((k+1)*mp.log(2*C*a)-mp.loggamma(k+2)
            +(count*mp.log(max(1, mp.mpf(count)/delta)) if count else 0))


def log_sum_bound(operator_bound, radius, strip_loss):
    """log[(z+e*z^2/delta) exp(e*z^2/delta)], z=2*C*a."""
    C, a, delta = map(mp.mpf, (operator_bound, radius, strip_loss))
    if C < 1 or a <= 0 or not 0 < delta <= 1:
        raise ValueError('A proved coefficient bound and positive geometry required')
    z = 2*C*a
    b = mp.e*z*z/delta
    return mp.log(z+b)+b


def proof_ledger():
    return dict(derivative_factors_in_length_k='ceil(k/2)',
                diagonal_exponents=[0, 0, 2, 0, 3, 1],
                construction='sum_{k>=0} [G(A0+A1*d_eta)]^k G f1',
                coefficient_bound_must_dominate_actual_matrices=True,
                actual_core_operator_envelope='C^16',
                numerical_actual_coefficient_bound_materialized=False,
                radius_may_exceed_inverse_operator_bound=True,
                parity=['even', 'even', 'even', 'even', 'odd', 'odd'],
                sum_bound='(z+e*z^2/Delta)*exp(e*z^2/Delta), z=2*C_op*a',
                conditional_on_inherited_holomorphic_inner_profile=True)
