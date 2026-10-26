"""Exact rate join for the sharper initial Haar target."""
from fractions import Fraction as Q


def exponents(kappa=Q(1,100000)):
    k=Q(kappa)
    if not 0<k<Q(1,4):
        raise ValueError('need 0<kappa<1/4')
    alpha=2-3*k
    D=alpha-k
    nu=D+1-2*k
    beta=2-2*k
    source=beta-k
    final=source+1-2*k
    return dict(kappa=k,alpha=alpha,D=D,nu=nu,beta=beta,source=source,final=final)


def divisor_data(k1,k2):
    if not isinstance(k1,int) or not isinstance(k2,int) or (k1,k2)==(0,0):
        raise ValueError('a nonzero integer Fourier mode is required')
    r=abs(k1)+abs(k2);pell=(k1+k2)**2-2*k2*k2
    if pell==0:raise AssertionError('irrational radial direction acquired an integer kernel')
    return dict(mode=(k1,k2),radius=r,pell_norm=pell,inverse_power_bound=27*r**4)


def shell_inverse_bound(jet_loss=7):
    if not isinstance(jet_loss,int) or jet_loss<7:
        raise ValueError('this absolute two-dimensional proof needs J>=7')
    shell_exponent=5-jet_loss
    return dict(jet_loss=jet_loss,shell_exponent=shell_exponent,
                summable=shell_exponent<=-2,coefficient=108*2**jet_loss,
                series_bound=2 if jet_loss==7 else 1+Q(1,jet_loss-6),
                estimate='||Qcirc q||_C0 <= C_J ||L^2 Nt^-1 L q||_CJ')


def target_transfer(jet_loss=7):
    return dict(inverse=shell_inverse_bound(jet_loss),
        support_identity=('q_beta(y_open,beta)=0 implies eta_beta^2 Xi='
                          '-Qcirc(q_beta)(y_open,beta); summing beta uses sum eta_beta^2=1'),
        operator_cannot_create_epsilon_gain=True,
        inherited_signed_maps_preserve_epsilon_order=True)


def target_gain(kappa=Q(1,100000)):
    o=exponents(kappa);k=o['kappa']
    return dict(orders=o,initial_Haar_rows=o['beta'],compatibility_moments=Q(3),
        projected_source=min(o['beta'],Q(3)),compact_target=o['beta'],
        target_gain=o['beta']-o['alpha'],signed_velocity=o['beta']-Q(1,2),
        signed_curl=o['beta']-k,real_old_new_sine=o['beta']+1-2*k,
        real_signed_square=2*o['beta']-1,real_signed_curl_square=2*o['beta']-2*k,
        real_pair_remainder=min(o['beta']+1-2*k,2*o['beta']-1,2*o['beta']-2*k),
        pointwise_covariance_response=o['source'])


def final_row(kappa=Q(1,100000)):
    o=exponents(kappa);k=o['kappa'];g=target_gain(k)
    source_groups=dict(original_fast_difference=o['alpha'],
        original_slow_difference=o['beta'],signed_covariance_response=o['source'],
        corrected_original_A1_complement=Q(2))
    fast_groups=dict(slow_time=o['source']+1,
        base_radial_flux=o['source']+1-k,radial_viscosity=o['final'],
        old_mean_transport=(1-k)+o['source']+1-k,
        changed_axial_pressure=o['source']+1,
        new_mean_square=2*o['source']+1-k)
    final_groups=dict(postsigned_Haar_row=g['real_pair_remainder'],
        fresh_fast_difference=min(fast_groups.values()),
        retained_A1_raw_bound=3-5*k,sharp_slow_difference=3-3*k)
    return dict(orders=o,target=g,source_groups=source_groups,
        complete_fast_source=min(source_groups.values()),fast_groups=fast_groups,
        final_groups=final_groups,conditional_final_order=min(final_groups.values()),
        strict_margin=min(final_groups.values())-o['nu'],
        inherited_required_order=o['nu'],conditional_strict_absorption=True)
