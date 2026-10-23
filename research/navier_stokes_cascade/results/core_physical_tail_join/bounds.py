"""Physical tail transfer and the remaining absolute regular-order inequality."""
from fractions import Fraction as Q
from math import ceil,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


final=load('physical_tail_final_mean',HERE.parent/'core_final_mean_row_closure/budgets.py')
exact,order=final.exact,final.order


def scales(h,kappa=Q(1,100000)):
    h=exact(h)
    if not 0<h<Q(1,100):raise ValueError('Original fixed 0<h<1/100 required')
    r=final.rates(kappa);A=Q(1,2)+h;normalization=2*A+Q(1,2)
    return dict(h=h,kappa=exact(kappa),A=A,normalization=normalization,**r)


def physical_loss(m,h):
    order(m);return scales(h)['normalization']+2*m


def compact_orders(e,m,N,h,kappa=Q(1,100000)):
    """Minimal integer for this strict-margin raw-kernel transfer rule.

    Not a lower bound on what every possible proof must use. N is a
    physical q power; source_rate and the operator gain are epsilon powers.
    """
    order(m);N=exact(N)
    if e not in (0,1) or N<0:raise ValueError('e=0/1 and nonnegative physical power required')
    s=scales(h,kappa);a=s['B'] if e==0 else s['D'];L=s['normalization']+2*m
    target=(N+L)/s['h'];p=max(1,(target-a)//s['kappa']+1)
    E=s['h']*(a+p*s['kappa'])-L
    return dict(e=e,physical_jet=m,requested_q_power=N,source_rate=a,
                normalized_target=target,integration_order=p,raw_source_order=m+2*p+3,
                physical_exponent=E,strict_margin=E-N,
                Q_to_q_factor_power=ceil(abs(E)),
                previous_integer_margin=s['h']*(a+(p-1)*s['kappa'])-L-N)


def absorption_constant(coefficients,K,margin):
    """Moderate exact diagnostic for sum c_j(16+K)^j(2j)!/margin^(2j).

    x=ln(1/q)>=1, S<=16 x^2. The actual huge polynomial remains its
    finite formula; evaluation here is not needed for the analytical bound.
    """
    K,margin=exact(K),exact(margin)
    if K<0 or margin<=0:raise ValueError('Nonnegative estimate shift and positive margin required')
    if any(type(j) is not int or j<0 or exact(c)<0 for j,c in coefficients.items()):
        raise ValueError('Positive polynomial coefficients required')
    return sum((exact(c)*(16+K)**j*factorial(2*j)/margin**(2*j) for j,c in coefficients.items()),Q(0))


def gaussian_recipe(c_band,D,b,N):
    c_band,D,b,N=map(exact,(c_band,D,b,N))
    if c_band<=0 or min(D,b,N)<0:raise ValueError('Positive Gaussian separation and nonnegative orders required')
    # ell>=4, q/Q in [1/4,4], x=ln(1/q): x<3ell/2 and
    # ell<4x. Thus c_band ell^2 >= (c_band/4)x^2.
    c=c_band/4;n=ceil(b/2)
    return dict(physical_Gaussian_constant=str(c),polynomial_index=n,
                coefficient=str(Q(factorial(n))*(2/c)**n),
                exponential_log=str((D+N)**2/(2*c)),
                bound='coefficient * exp(exponential_log) * q^N',
                same_domain='original ell>=16; no jet-dependent starting band')


def regular_transfer(m,h,kappa=Q(1,100000),reserve=Q(1)):
    order(m);s=scales(h,kappa);reserve=exact(reserve)
    if reserve<=0:raise ValueError('Positive polynomial absorption reserve required')
    L=s['normalization']+2*m;E=s['h']*s['D']-L
    return dict(physical_jet=m,normalized_regular_order=str(s['D']),
                native_physical_exponent=str(E),absorbed_physical_exponent=str(E-reserve),
                required_regular_order_for_qN=f'alpha > (N + {L} + {reserve}) / h',
                polynomial_reserve=str(reserve),absolute_decay_certified=E-reserve>0,
                actual_residual_lower_bound=False)


def hypothetical_steps(m,N,h,kappa=Q(1,100000),reserve=Q(1)):
    """Only arithmetic IF a later full same-family induction supplies orders.

    No future stage or infinite diagonal is certified by this formula.
    """
    order(m);N=exact(N);reserve=exact(reserve);s=scales(h,kappa)
    if N<0 or reserve<=0:raise ValueError('Nonnegative target and positive reserve required')
    gamma=Q(1,2)-s['kappa'];L=s['normalization']+2*m+reserve
    j=max(0,(N+L-s['h']*s['D'])//(s['h']*gamma)+1)
    return dict(additional_steps=j,rate_if_supplied=str(s['D']+j*gamma),
                physical_margin=str(s['h']*(s['D']+j*gamma)-L-N),
                induction_supplied=False,arbitrary_diagonal_supplied=False)


class PhysicalTailJoin:
    """Same final prefix and original operator/profile geometry.

    No additive residual tail is used as a sourced inverse input, nor
    inserted as a velocity in a nonlinear interaction.
    """
    def __init__(self,closure,h):
        self.closure=closure;self.s=scales(h,closure.kappa)

    def compact_recipe(self,e,m,N):
        a=compact_orders(e,m,N,self.s['h'],self.s['kappa']);c=self.closure
        p0=Q(c.means.reweighting['original_pole'])
        source='Q1 current final gr' if e==0 else 'Q1 postwave Ez'
        r=final.oscillatory_recipe(m,a['integration_order'],e=e,source=source,
            source_rate=a['source_rate'],target_rate=a['normalized_target'],
            kappa=self.s['kappa'],geometry=c.means.reweighting,pole=12*p0)
        if e==0:
            src=dict(method='FinalMeanRowClosure.current_radial_source_polynomial',order=a['raw_source_order'],
                     field_input_order=a['raw_source_order']+3,embedding=str(c.embedding))
        else:
            old=final.fast.flat_recipe(0,self.s['kappa'])['raw_source_order']
            src=dict(method='FastRowGain.postwave_recipe',order=a['raw_source_order'],
                     nested_old_source_order=a['raw_source_order']+old,
                     nested_original_initial_mean_order=a['raw_source_order']+old+3)
        return dict(orders={k:str(v) for k,v in a.items()},original_kernel_recipe=r,
            source_recipe=src,
            chart_factor='4^'+str(a['Q_to_q_factor_power']),
            polynomial_absorption='sum_j c_j (16+K)^j (2j)! / strict_margin^(2j)',
            exact_tail='-A0 Q1 gr' if e==0 else '+A1 Q1 Ez_postwave',
            materialized_actual_high_jet_bound=False)

    def certificate(self,m,N):
        return dict(A0=self.compact_recipe(0,m,N),A1=self.compact_recipe(1,m,N),
            regular=regular_transfer(m,self.s['h'],self.s['kappa']),
            base_exterior='independently retained exact flat remainder AFTER its stress divergence is assigned to mean equations',
            Gaussian='sum of original homogeneous temporal commutators and sourced (1-psi)CURRENT_f + Tpsi',
            total_identity='Rphysical = Rregular + Tcompact + TGaussian + Tbase_exterior',
            spatial_slow_cutoffs='remain inside complete regular curls/pressure/viscosity/interactions',
            later_diagonal_and_localization='not Gaussian; require their own product/jet estimates',
            full_actual_inputs_or_absolute_residual_improvement_certified=False)
