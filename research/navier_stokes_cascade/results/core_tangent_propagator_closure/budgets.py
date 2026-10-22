"""Close the tangent/length inputs from ORIGINAL retained frame budgets."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


transfer=load('tangent_closure_transfer',HERE.parent/'core_sourced_normal_transfer/budgets.py')
normal=transfer.normal
P,exact=transfer.P,transfer.exact


def certificate(p,selector_inputs,ell):
    """Analytical implication from valid ORIGINAL bounds, not field validation.

    The exponential is retained as a finite constant recipe. It is not a
    scale-dependent Gronwall loss or an inverse Gaussian bound.
    """
    p=normal.frame_inputs(p);old=normal.verify_existing(p,selector_inputs,ell)
    if not old['conditional_normal_floor_implication_passes']:
        raise ValueError('Original accepted-band and retained length premises must hold')
    b=old['original_budget']['original_budgets']
    H=b['CE']+b['Cd'];N=-(-H.numerator//H.denominator)
    kappa=(1+b['S0'])*max(Q(1),p['cmax'])/min(Q(1),p['cmin'])
    return dict(original_normal_certificate=old,frame_conversion=kappa,
        integrated_error=H,integer_exponential_upper_power=N,
        tangent_constant_recipe='frame_conversion * 3^integer_exponential_upper_power',
        length_polynomial=P({1:1}).record(),damping_floor=p['bmin']**2/4,
        gates=dict(original_normal_floor=old['conditional_normal_floor_implication_passes'],
                   original_length_below_S=old['gates']['original_geometric_length_upper_is_below_S'],
                   nonnegative_finite_integrated_error=H>=0,
                   rational_frame_conversion_at_least_one=kappa>=1,
                   integer_upper_exponent=N>=H),
        tangent_and_length_inputs_derived_from_retained_frame_bounds=True,
        new_acceptance_or_smallness_condition=False,
        actual_primitive_enclosures_certified=False,
        homogeneous_primary_frames_or_roots_identified=False,
        complete_covariance_means_or_residual_certified=False)


def tangent_polynomial(c):
    """Materialize a finite rational upper bound; never approximate exp(H).

    Intended for supplied finite inputs, which may be manufactured. The
    mathematical certificate can retain its exact recipe without allocating
    an enormous selected integer. This helper does not certify the inputs.
    """
    H=exact(c['integrated_error']);kappa=exact(c['frame_conversion']);N=c['integer_exponential_upper_power']
    if not all(c['gates'].values()) or H<0 or kappa<1 or type(N) is not int or N!=-(-H.numerator//H.denominator):
        raise ValueError('Failed tangent certificate or exponential recipe')
    return P(kappa*3**N)


def forced_from_frame(p,selector_inputs,ell,*,F,G,btilde,path,source,harmonic=1):
    """Same full-path inverse: derive T and length rather than assume them."""
    c=certificate(p,selector_inputs,ell)
    phase,forced,old=transfer.forced_from_original(p,selector_inputs,ell,
        F=F,G=G,btilde=btilde,path=path,source=source,
        tangent=tangent_polynomial(c),length=P({1:1}),harmonic=harmonic)
    return phase,forced,c


def residence_upper(p,length,harmonic):
    """Tangent zero-order residence only; NOT the ambient inverse length."""
    p=normal.frame_inputs(p);length=exact(length)
    if type(harmonic) is not int or not harmonic or length<=0:
        raise ValueError('Nonzero integer harmonic and positive length required')
    delta=harmonic**2-1
    return length if not delta else min(length,4/(delta*p['bmin']**2))
