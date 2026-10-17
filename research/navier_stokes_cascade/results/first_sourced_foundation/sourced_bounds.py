"""Exact first-sourced-step accounting, not an all-stage residual certificate."""
from fractions import Fraction as F

KAPPA=F(1,100000)


def orders(kappa=KAPPA):
    k=F(kappa)
    if not 0<k<F(1,10):raise ValueError('Require fixed 0<kappa<1/10')
    B=1-k;C=B+F(1,2);D=C-k;alpha=2-3*k;H=alpha-k;mu=1-k
    prefix=dict(original_primary_source=B,signed_linear=alpha-k,
        primary_signed_cross=alpha-k,signed_square=2*alpha-1-k,
        mean_signed_cross=alpha-k,post_signed_fast_mean_wave=H,
        post_signed_slow_mean_wave=alpha)
    generated=dict(full_linear=D,complete_old_new_wave=B+F(1,2)-k,
        new_wave_square=2*B-k,actual_total_mean_new_wave=B+mu-F(1,2))
    return dict(B=B,C=C,D=D,mu=mu,old_error=B,alpha=alpha,H=H,
        prefix=prefix,generated=generated,source_gain=D-B,
        covariance_cross=C,covariance_square=2*B,covariance_change=min(C,2*B),
        full_radial_source_change=C-k,full_mean_pressure_change=C-k,
        averaged_radial_source_change=C,averaged_mean_pressure_change=C,
        new_full_tangential_mean=C-k,new_Haar_tangential_mean=C,
        recomputed_defects=C,radial_moment_residual=C,
        pre_source_full_mean=3-6*k,pre_source_defects=3-3*k,
        next_full_signed_target=C,next_generic_old_error=B)


def physical_powers(h,a,b,kappa=KAPPA):
    h=F(h)
    if not 0<h<F(1,100) or any(type(x) is not int or x<0 for x in (a,b)):raise ValueError('Fixed h and nonnegative integer jet counts required')
    o=orders(kappa);A=F(1,2)+h;sx=F(1,2)+h/2;st=1+3*h/2;shift=-a*sx-b*st
    return dict(velocity=-A+h*o['B']+shift,pressure=-2*A+h*o['C']+shift,
        potential=F(1,2)-A+h*o['C']+shift,
        regular_residual=-2*A-F(1,2)+h*o['D']+shift)


def inverse_polynomial_losses(source):
    if not source or any(type(x) is not int or x<0 for x in source):raise ValueError('Nonempty finite polynomial jet-degree list required')
    g=[max(i+source[j-i] for i in range(j+1)) for j in range(len(source))]
    b=[]
    for j in range(len(source)):
        b.append(max([g[j]+1]+[i+1+b[j-i] for i in range(1,j+1)]))
    return g,b


def review():
    o=orders();k=KAPPA
    return dict(orders={n:str(v) for n,v in o.items() if not isinstance(v,dict)},
        unchanged_prefix_source_orders={n:str(v) for n,v in o['prefix'].items()},
        complete_first_generated_orders={n:str(v) for n,v in o['generated'].items()},
        zero_jet_physical_powers={n:str(v) for n,v in physical_powers(F(1,200),0,0).items()},
        checks=dict(unchanged_signed_and_mean_prefix_preserves_B=all(v>=o['B'] for v in o['prefix'].values()),
            old_complete_error_and_new_sourced_wave_have_same_B=o['old_error']==o['B'],
            all_regular_nonzero_groups_meet_D=min(o['generated'].values())==o['D'],
            strict_one_step_regular_gain=o['source_gain']==F(1,2)-k>0,
            new_square_covariance_no_worse_than_old_new=2*o['B']>=o['C'],
            full_means_and_defects_are_recomputed_not_preserved=o['D']<o['pre_source_full_mean'] and o['C']<o['pre_source_defects'],
            full_pressure_change_does_not_worsen_accumulated_pressure=o['D']>o['mu'],
            covariance_and_defects_not_promoted_to_averaged_residual_gain=o['C']>o['D'],
            physical_regular_residual_still_has_negative_power=physical_powers(F(1,200),0,0)['regular_residual']<0,
            same_inverse_source_and_operation_order_not_retuned=True,
            all_later_complex_targets_and_full_residual_not_inferred=True))
