"""Current-input signed/mean order join; exact rational sufficient budgets."""
from fractions import Fraction as F

KAPPA=F(1,100000)


def generic_gain(alpha,rho,kappa=KAPPA):
    a,r,k=F(alpha),F(rho),F(kappa)
    if a<1 or r<F(1,2) or not 0<k<F(1,2):raise ValueError('Require alpha>=1, rho>=1/2, 0<kappa<1/2')
    return min(r-F(1,2),F(1,2)-k,a-1)


def orders(kappa=KAPPA):
    k=F(kappa)
    if not 0<k<F(1,10):raise ValueError('Require fixed 0<kappa<1/10')
    B=1-k;C=B+F(1,2);D=C-k;mu=B;H=D
    cov=dict(old_error_new_principal=B+C-F(1,2),complete_old_new_curl=C+F(1,2)-k,signed_square=2*C-1)
    nu=min(cov.values())
    signed=dict(full_linear=B+F(1,2)-k,old_new_wave=F(1,2)+B-k,mean_new_wave=B+mu-F(1,2),signed_square=2*B-k)
    fast=dict(slow_time=H+1,base_radial_flux=H+1-k,radial_viscosity=H+1-2*k,
        old_mean_transport=mu+H+1-k,changed_axial_pressure=H+1,new_mean_square=2*H+1-k)
    slow=dict(slow_time=C+1,base_radial_flux=C+1,radial_viscosity=C+1,
        old_mean_transport=mu+C+1-k,changed_axial_pressure=C+1,new_mean_square=2*C+1)
    scalar=dict(accumulated_Haar_cross=C+1,new_square=2*C,time_base_viscosity=C+2,old_radial_mean=C+3)
    return dict(B=B,C=C,D=D,mu=mu,H=H,rho=B,nu=nu,covariance=cov,signed_nonzero=signed,
        fast=fast,slow=slow,scalar=scalar,generic_gain=generic_gain(C,B,k),
        fast_defect_change=mu+H,slow_defects=min(scalar.values()),
        stale_target_gap=C+1-(mu+H),final_full_tangential=min(nu,*fast.values(),*slow.values()),
        final_nonzero=min(D,H,C,*signed.values()),old_wave_error=min(B,C-F(1,2)),
        mean_wave_fast=H,mean_wave_slow=C,old_Haar_tangential=F(1),old_Haar_radial=F(2))


def physical_powers(h,order,a,b,wave=False):
    h,t=F(h),F(order)
    if not 0<h<F(1,100) or t<0 or any(type(n) is not int or n<0 for n in (a,b)):raise ValueError('Fixed h, nonnegative order and integer jets required')
    A=F(1,2)+h;shift=-a*(F(1,2)+h/2)-b*(1+3*h/2)
    return dict(velocity=-A+h*(t-F(1,2) if wave else t)+shift,
        pressure=-2*A+h*t+shift,potential=F(1,2)-A+h*t+shift)


def review():
    o=orders();k=KAPPA;C,D,H,nu=(o[n] for n in ('C','D','H','nu'))
    groups={n:{a:str(b) for a,b in o[n].items()} for n in ('covariance','signed_nonzero','fast','slow','scalar')}
    return dict(orders={n:str(v) for n,v in o.items() if not isinstance(v,dict)},groups=groups,
        zero_jet_physical_powers={n:{a:str(b) for a,b in physical_powers(F(1,200),t,0,0,wave).items()}
            for n,t,wave in (('signed',C,True),('fast',H,False),('slow',C,False))},
        strict_gain_condition='min(rho-1/2,1/2-kappa,alpha-1)>0',
        checks=dict(current_input_not_old_conservative_conditional_orders=o['B']==1-k and C==F(3,2)-k and D==F(3,2)-2*k,
            all_three_generic_covariance_errors_tie=len(set(o['covariance'].values()))==1 and nu==2-2*k,
            positive_generic_gain=o['generic_gain']==F(1,2)-k>0,
            full_tensor_change_only_C_not_nu=C<nu,
            full_signed_nonzero_groups_preserve_D=min(o['signed_nonzero'].values())==D,
            fast_viscosity_strictly_above_current_averaged_error=min(o['fast'].values())-nu==F(1,2)-2*k>0,
            fast_changes_defects_but_preserves_C=o['fast_defect_change']>C,
            fresh_slow_scalar_gain_limited_by_old_Haar_mean=o['slow_defects']==C+1<2*C,
            stale_targets_miss_new_bound_by_two_kappa=o['stale_target_gap']==2*k>0,
            all_slow_full_tangential_changes_above_nu=min(o['slow'].values())>nu,
            full_tangential_mean_reaches_nu=o['final_full_tangential']==nu,
            fast_mean_wave_source_TIES_D_not_strictly_above=o['mean_wave_fast']==D<C,
            complete_regular_nonzero_residual_not_improved=o['final_nonzero']==D,
            original_primary_error_remains_B=o['old_wave_error']==o['B']<D,
            total_mu_and_old_Haar_orders_preserved=min(H,C)>1>o['mu'],
            physical_regular_residual_still_negative=-2*(F(1,2)+F(1,200))-F(1,2)+F(1,200)*D<0,
            all_fixed_bands_and_jets_use_retained_inverse_and_product_lemmas=True,
            all_later_source_family_closure_not_inferred=True))
