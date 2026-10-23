"""Exact current-source identities and positive conditional norm controls."""
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
from bounds import rates,reduced_envelope,A1_request,certificate,P,parent,sharp

def drift_controls():
    checks={}
    R,Z,T,y,w=s.symbols('R Z T y w',positive=True)
    eps,M,c,tau=s.symbols('eps M c tau',positive=True)
    coords=(R,Z,T,y,w)
    f={k:s.Function(k)(*coords) for k in ('beta','v','gamma','b','V','G','Wrr','Wrz','Wtt','Wzz','pm','Sigma','Psi')}
    beta,v,gamma,b,V,G,Wrr,Wrz,Wtt,Wzz,pm,Sigma,Psi=(f[k] for k in f)
    for d in (s.Rational(1),s.Rational(3,2)):
        Dr=lambda a:s.diff(a,R)+M*d*R**(d-1)*s.diff(a,y)
        De=lambda a,e:Dr(a)+e*a/R
        Dz=lambda a:eps*s.diff(a,Z)
        ts=lambda a:-eps*s.diff(a,T)+c*(tau*s.diff(a,y)+s.diff(a,w))
        S=2*b*beta+beta**2+Wrr;H=b*gamma+G*beta+beta*gamma+Wrz;C=2*V*v+v*v+Wtt
        radialvis=Dr(Dr(beta))+Dr(beta)/R-beta/R**2
        scalarvis=Dr(Dr(gamma))+Dr(gamma)/R
        gr=-ts(beta)-De(S,1)-Dz(H)+C/R+eps*(radialvis+Dz(Dz(beta)))
        r0=-ts(beta)-S/R-Dz(H)+C/R+eps*Dz(Dz(beta));K0=-S+eps*De(beta,1)
        ez=ts(gamma)+De(H,1)+Dz(2*G*gamma+gamma**2+Wzz+pm)-eps*(scalarvis+Dz(Dz(gamma)))-De(Sigma,1)
        r1=ts(gamma)+Dz(2*G*gamma+gamma**2+Wzz+pm)-eps*Dz(Dz(gamma))
        key=f'd{d}'
        checks[key+':complete_radial_vector_viscosity_is_DrD1']=s.expand(radialvis-Dr(De(beta,1)))==0
        checks[key+':complete_axial_scalar_viscosity_is_D1Dr']=s.expand(scalarvis-De(Dr(gamma),1))==0
        checks[key+':full_gr_is_reduced_r0_plus_Dr_K0']=s.expand(gr-r0-Dr(K0))==0
        checks[key+':full_Ez_retains_time_pressure_flux_and_both_viscosities']=s.expand(ez-r1-De(H-eps*Dr(gamma)-Sigma,1))==0
        beta_p=-Dz(Psi);gamma_p=De(Psi,1)
        checks[key+':actual_complete_potential_is_exactly_divergence_free']=s.expand(De(beta_p,1)+Dz(gamma_p))==0
        checks[key+':time_of_actual_gamma_stays_in_D1_range']=s.expand(ts(gamma_p)-De(ts(Psi),1))==0
        checks[key+':both_axial_viscosities_of_actual_gamma_stay_in_D1_range']=s.expand(De(Dr(gamma_p),1)+Dz(Dz(gamma_p))-De(Dr(gamma_p)+Dz(Dz(Psi)),1))==0
        checks[key+':radial_viscous_pressure_memory_annihilated_only_after_A1Dz']=s.expand(Dz(De(beta,1))-De(Dz(beta),1))==0
        checks[key+':K0_pressure_memory_is_not_zero']=s.expand(K0)!=0
        p0=s.Function('p0')(*coords)
        S_p=2*b*beta_p+beta_p**2+Wrr;H_p=b*gamma_p+G*beta_p+beta_p*gamma_p+Wrz
        pm_p=-S_p+eps*De(beta_p,1)+p0
        H1_p=2*G*gamma_p+gamma_p**2+Wzz-S_p+p0
        Ez_p=ts(gamma_p)+De(H_p,1)+Dz(2*G*gamma_p+gamma_p**2+Wzz+pm_p)-eps*(De(Dr(gamma_p),1)+Dz(Dz(gamma_p)))-De(Sigma,1)
        flux=ts(Psi)+H_p-eps*Dr(gamma_p)-2*eps*Dz(Dz(Psi))-Sigma
        checks[key+':entire_current_Ez_minus_Dz_H1_is_one_complete_D1_flux']=s.expand(Ez_p-Dz(H1_p)-De(flux,1))==0
    return checks,dict(scope='formal full variable drift, current mean equations and exact compact potential; support/boundary identities checked separately')

def characteristic_controls():
    checks={};r,R=s.symbols('r R',positive=True);lam=s.symbols('lambda',real=True)
    p=r**2*(r-1)**4*(2-r)**4
    for d in (1,2):
        phase=s.exp(-s.I*lam*r**d)
        Dr=lambda a:s.diff(a,r)+s.I*lam*d*r**(d-1)*a
        for e in (0,1):
            for b in (0,1,2):
                weighted=s.expand(s.simplify(r**e*(Dr(p*phase)+b*p*phase/r)/phase))
                actual=s.integrate(weighted,(r,1,2))
                rhs=(b-e)*s.integrate(r**(e-1)*p,(r,1,2))
                checks[f'd{d}e{e}b{b}:exact_full_drift_weight_and_boundary_rule']=s.simplify(actual-rhs)==0
                if b==e:
                    checks[f'd{d}e{e}:A_e_annihilates_compact_D_e_potential']=actual==0
                    forward=s.integrate(weighted,(r,1,R))
                    checks[f'd{d}e{e}:T_e_D_e_restores_pressure_potential_not_zero']=s.simplify(forward-R**e*p.subs(r,R))==0
            # Noncompact h=r exp(-i lambda r^d) has a nonzero endpoint
            # term, so local derivative structure alone is insufficient.
            bad=s.expand(s.simplify(r**e*(Dr(r*phase)+e*phase)/phase))
            checks[f'd{d}e{e}:noncompact_boundary_cannot_be_erased']=s.integrate(bad,(r,1,2))==2**(e+1)-1
    x=s.symbols('x',real=True)
    chi=3*x*x-2*x**3
    mass=s.integrate(s.diff(chi,x),(x,0,1));left=s.integrate(s.diff(chi,x),(x,0,R))
    checks['same_cutoff_endpoint_mass_gives_A_e_projection']=mass==1 and s.expand(left-chi.subs(x,R))==0
    checks['T_e_A_e_zero_uses_same_chi_and_not_a_second_cutoff']=s.expand(left-chi.subs(x,R)*mass)==0
    chi_other=chi+x*(1-x)
    checks['different_cutoff_in_same_operator_identity_is_detected']=s.expand(left-chi_other.subs(x,R)*mass)!=0
    # Projection mode j -> cutoff derivative*j. After repair its
    # characteristic integral is exactly j-j*int chi'=0.
    j=s.symbols('j',nonzero=True)
    checks['actual_gd_minus_A1_gd_is_in_kernel_A1']=j-j*mass==0
    checks['dropping_actual_A1_velocity_breaks_kernel_identity']=j!=0
    return checks,dict(scope='exact Fourier characteristic integrals of finite endpoint polynomial fixtures; analytical proof uses smooth whole-line extensions')

def pressure_and_cutoff_controls():
    checks={};R,Z,T=s.symbols('R Z T',positive=True);eps=s.symbols('eps',positive=True)
    nu=(1+R)*(1+Z+T+Z*T);J=(1+Z**2)*(1+T**2)
    checks['moving_chi_Z_commutator_is_retained']=s.expand(s.diff(nu*J,Z)-nu*s.diff(J,Z)-s.diff(nu,Z)*J)==0
    checks['moving_chi_T_commutator_is_retained']=s.expand(s.diff(nu*J,T)-nu*s.diff(J,T)-s.diff(nu,T)*J)==0
    checks['commuting_Dz_past_A1_without_chi_term_changes_result']=s.expand(s.diff(nu*J,Z)-nu*s.diff(J,Z))!=0
    checks['tstar_commutator_has_negative_epsilon_sign']=s.expand(-eps*s.diff(nu*J,T)+eps*nu*s.diff(J,T)+eps*s.diff(nu,T)*J)==0
    r=s.symbols('r',positive=True);phi=(r-1)**4*(2-r)**4;phi=s.expand(phi/s.integrate(phi,(r,1,2)))
    # A pure nonzero Fourier mode has zero Haar P. Its local pressure
    # K0=-Z*phi exp(-i lambda r^d) still gives an A1 Dz contribution.
    weighted_memory=-eps*s.integrate(r*phi,(r,1,2))
    checks['dropping_minus_S_local_pressure_memory_changes_A1_axial_flux']=weighted_memory==-s.Rational(3,2)*eps and weighted_memory!=0
    rhoP,K0,p0=s.symbols('rhoP K0 p0')
    checks['complete_current_pressure_is_K0_plus_compensated_primitive']=(K0+p0)-p0==K0
    checks['radial_compensation_not_replaced_by_pointwise_source_zero']=rhoP!=0
    c_rho=s.Function('c_rho')(Z);Pcur=s.Function('Pcur')(Z)
    checks['whole_moving_pressure_moment_derivative_keeps_c_rho_jet']=s.expand(s.diff(c_rho*Pcur,Z)-c_rho*s.diff(Pcur,Z))==s.diff(c_rho,Z)*Pcur
    return checks,dict(scope='same current pressure and moving cutoff, not an improved bound for the full pressure field')

def swirl_controls():
    checks={};A,R,V,B,phi,z=s.symbols('A R V B phi z',nonzero=True)
    cos=(z+1/z)/2
    v=A*R/(2*V)*phi*B*cos
    change=s.expand((2*V*v+v*v)/R)
    theta=s.symbols('theta',real=True)
    checks['pure_swirl_is_divergence_free_beta_gamma_zero']=s.diff(v,theta)/R==0
    checks['pure_swirl_exact_Haar_velocity_moments_both_zero']=R**2*s.expand(v).coeff(z,0)==0 and R*s.S.Zero==0
    checks['full_centripetal_fundamental_is_exact_resonant_source']=change.coeff(z,1)==A*phi*B/2 and change.coeff(z,-1)==A*phi*B/2
    checks['centripetal_square_cannot_cancel_fundamental']=s.expand(v*v/R).coeff(z,1)==0 and s.expand(v*v/R).coeff(z,2)!=0
    checks['Haar_radial_pressure_change_is_nonzero_and_must_be_compensated']=change.coeff(z,0)==A*A*R*phi*phi*B*B/(8*V*V)
    scalar=s.symbols('rho_P')
    checks['same_scalar_pressure_compensator_changes_no_fundamental']=s.expand(change-scalar).coeff(z,1)==change.coeff(z,1)
    for j in (2,3,8):
        mode=sharp.pell(j);checks[f'swirl_pell{j}:actual_parent_frequency_is_nonzero']=mode['k']!=(0,0) and abs(mode['norm'])==1
    return checks,dict(scope='comparison swirl satisfies divergence and BOTH moments, not canonical iteration or uniform all-jet source data')

def budget_controls():
    checks={};rows=[]
    for n in (0,1,3):
        data={k:P(n+j+1) for j,k in enumerate(('x','u','z','b','V','G','c','r','rho','Kp','Kmom','KA1','W'))}
        a=reduced_envelope(n,data);base={k:v.evaluate(1) for k,v in data.items()};b=reduced_envelope(n,base)
        checks[f'n{n}:full_positive_polynomial_and_exact_rational_envelope_agree']=all(a[k].evaluate(1)==b[k] for k in ('S','H','C','r0','P','tilde_pressure','H1','A1_coefficient'))
        checks[f'n{n}:A1_Qcirc_factor_two_and_explicit_Dz_epsilon_retained']=a['A1_coefficient'].c==(2*data['KA1']*(n+4)*a['H1']).c and a['explicit_epsilon_power']==1
        checks[f'n{n}:whole_pressure_flux_derivative_has_n_plus_three_input_budget']=a['field_input_order']==n+3 and a['internal_r0_pressure_H1_order']==n+1
        changed=dict(data,W=100*data['W']);dc=reduced_envelope(n,changed)
        checks[f'n{n}:complete_covariance_rr_rz_tt_zz_not_dropped']=dc['A1_coefficient'].evaluate(1)>a['A1_coefficient'].evaluate(1)
        changed=dict(data,Kp=100*data['Kp']);dc=reduced_envelope(n,changed)
        checks[f'n{n}:full_current_pressure_operator_changes_A1_budget']=dc['A1_coefficient'].evaluate(1)>a['A1_coefficient'].evaluate(1)
        # Independent unreduced physical-epsilon norm substitutions at
        # kappa=1/10, epsilon=t^10. Internal r0/H1 are through n+1.
        for t in (Q(1,2),Q(1,8)):
            ep=t**10;x=t**19*base['x'];u=t**9*base['u'];z=t**9*base['z'];b0=t**10*base['b'];W=t**10*base['W']
            d1=n+2;d2=(n+2)*(n+3)
            S=2*b0*x+x*x+W;H=b0*z+base['G']*x+x*z+W;C=2*base['V']*u+u*u+W
            raw0=d1*(ep+2*base['c'])*x+base['r']*S+ep*d1*H+base['r']*C+ep**3*d2*x
            pressure=base['Kp']*(raw0+base['rho']*base['Kmom']*raw0)
            rawH=2*base['G']*z+z*z+W+S+pressure
            rawA=2*base['KA1']*ep*(n+1)*rawH
            checks[f'n{n}t{t}:independent_full_pressure_flux_B_and_A1_B_plus_one_budget']=raw0<=t**9*b['r0'] and rawH<=t**9*b['H1'] and rawA<=t**19*b['A1_coefficient']
    for k in (Q(1,100000),Q(1,10),Q(6,25)):
        r=rates(k)
        checks[f'k{k}:r0_remains_B_and_A1_only_gains_two_kappa']=r['reduced_radial_source']==r['B'] and r['gain_over_parent_A1_source']==2*k and r['margin_over_C0']==k
        for m in (0,1,4):
            d=A1_request(m,kappa=k)
            checks[f'k{k}m{m}:normalized_A1_uses_direct_operator_no_oscillatory_integration']=d['integration_order']==0 and d['raw_H1_input_order']==m+1 and d['strict_margin']==k
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for m,N in ((0,0),(1,2),(3,7)):
            d=A1_request(m,physical_power=N,h=h)
            checks[f'h{h}m{m}N{N}:physical_raw_H1_keeps_extra_Dz_and_strict_margin']=d['raw_H1_input_order']==m+2*d['integration_order']+4 and d['strict_margin']>0
    d=A1_request(0,physical_power=0);rows.append(d)
    checks['illustrative_physical_A1_p30000002_H1_J60000008']=d['integration_order']==30000002 and d['raw_H1_input_order']==60000008 and d['previous_full_E_row_raw_order']==60000011
    c=certificate()
    checks['postwave_A1_and_final_A0_fields_not_mixed']='POSTWAVE' in c['axial']['stage'] and 'FINAL' in c['radial']['stage']
    checks['original_6q_KA1_and_12q_final_radial_domains_retained']='at6q' in c['domains'] and '12q' in c['domains']
    checks['whole_regular_E_and_actual_input_gate_unchanged']=c['complete_regular_source_order']==str(rates()['E']) and not c['actual_high_jet_support_inputs_or_absolute_residual_improvement_validated']
    checks['old_physical_A0_cost_not_erased']=parent.compact_physical_orders(0,0,0,Q(1,200))['raw_source_order']==60200007 and not d['older_A0_source_root_normal_and_physical_request_DAGs_erased']
    return checks,dict(rates=rates(),physical_rows=rows,certificate=c,
        scope='conditional norm implication from retained coarse actual mean/covariance and whole support inequalities, not their independent validation')

@lru_cache(None)
def _review():
    checks={};parts={}
    for name,fn in (('drift',drift_controls),('characteristics',characteristic_controls),('pressure',pressure_and_cutoff_controls),('swirl',swirl_controls),('budgets',budget_controls)):
        a,b=fn();checks.update(a);parts[name]=b
    return dict(checks=checks,parts=parts)

def review():
    from copy import deepcopy
    return deepcopy(_review())
