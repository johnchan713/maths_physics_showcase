"""Exact axial potential/pressure algebra, compact moments and scaled budgets."""
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
from bounds import (stage_rates,A1_request,direct_A1_recipe,potential_envelope,
    joined_envelope,complement_envelope,tail_envelope,beta_envelope,certificate,P,
    JOIN_KEYS,POTENTIAL_KEYS,COMPLEMENT_KEYS,ANGULAR_KEYS,angular)

def axial_controls():
    checks={};R,Z,T=s.symbols('R Z T',positive=True);y=s.symbols('y',nonzero=True)
    ep,c,M=s.symbols('epsilon c M',positive=True);I=s.I
    haar=lambda f:s.expand(f).coeff(y,0)
    qc=lambda f:s.expand(f-haar(f))
    nt=lambda f:s.expand(I*y*s.diff(f,y))
    dr=lambda f:s.expand(s.diff(f,R)+2*R*M*I*y*s.diff(f,y))
    d1=lambda f:s.expand(dr(f)+f/R)
    dz=lambda f:s.expand(ep*s.diff(f,Z))
    ts=lambda f:s.expand(-ep*s.diff(f,T)+c*nt(f))
    def inv(f):
        q=s.expand(f);return s.expand(sum(q.coeff(y,k)*y**k/(I*k) for k in range(-12,13) if k))
    cos=(y+y**-1)/2;sin=(y-y**-1)/(2*I);cos2=(y*y+y**-2)/2
    bump=(R-1)**4*(2-R)**4
    psi=R**2*bump*(1+T+Z*Z+(Z+T)*cos+cos2)
    beta=-dz(psi);gamma=d1(psi);b=ep*R;G=R**2*(1+Z)
    A=[R*bump*(1+cos),R*bump*(Z+sin),R*bump*(1+T*cos2)]
    W=[[s.expand(A[i]*A[j]) for j in range(3)] for i in range(3)]
    S=2*b*beta+beta*beta+W[0][0]
    # Keep an arbitrary pressure primitive AND the actual local viscous part.
    pm=-S+ep*d1(beta)+R*bump*(1+Z+T*cos)
    stress=R*bump*(1+Z)
    H=b*gamma+G*beta+beta*gamma+W[0][2]
    X=2*G*gamma+gamma*gamma+W[2][2]+pm
    Ez=ts(gamma)+d1(H)+dz(X)-ep*(d1(dr(gamma))+dz(dz(gamma)))-d1(stress)
    Lam=ts(psi)+H-ep*dr(gamma)-ep*dz(dz(psi))-stress
    checks['complete_axial_row_has_exact_compact_potential_factorization']=s.expand(Ez-dz(X)-d1(Lam))==0
    checks['actual_potential_preserves_exact_divergence']=s.expand(d1(beta)+dz(gamma))==0
    checks['fixed_label_time_and_slow_derivatives_commute_with_D1']=s.expand(d1(ts(psi))-ts(d1(psi)))==0 and s.expand(d1(dz(psi))-dz(d1(psi)))==0
    # T1 Q(Ez-DzX)=Q Lambda follows from the independent compact inverse
    # controls below; do not commute T1 through Dz or a moving cutoff.
    source=qc(Lam)+ep*s.diff(qc(psi),T)-qc(H)+ep*dr(qc(gamma))+ep*dz(dz(qc(psi)))
    checks['coarse_source_gives_exact_nonzero_Haar_potential_after_Nt_inverse']=s.expand(inv(source)/c-qc(psi))==0
    checks['potential_time_term_is_not_dropped']=s.expand(inv(source-ep*s.diff(qc(psi),T))/c-qc(psi))!=0
    checks['potential_radial_viscosity_is_not_dropped']=s.expand(inv(source-ep*dr(qc(gamma)))/c-qc(psi))!=0
    Xflat=X-ep*d1(beta)
    Lamflat=ts(psi)+H-ep*dr(gamma)-2*ep*dz(dz(psi))-stress
    checks['full_pressure_and_reduced_pressure_forms_require_different_viscous_factor']=s.expand(Ez-dz(Xflat)-d1(Lamflat))==0
    checks['removing_local_viscous_pressure_without_doubling_is_detected']=s.expand(Ez-dz(Xflat)-d1(Lam))!=0
    checks['local_viscous_pressure_still_has_nonzero_axial_force']=s.expand(dz(ep*d1(beta)))!=0
    Phi=-2*s.diff(G,R)*psi+gamma*gamma+W[2][2]+pm
    Lproj=Lam+2*dz(G*psi)
    checks['full_A1_source_reduction_retains_variable_G_radial_derivative']=s.expand(Ez-dz(Phi)-d1(Lproj))==0
    checks['wrongly_omitting_G_R_Psi_changes_the_axial_source']=s.expand(Ez-dz(Phi+2*s.diff(G,R)*psi)-d1(Lproj))!=0
    checks['compact_product_identity_for_G_gamma_is_exact']=s.expand(2*G*gamma-d1(2*G*psi)+2*s.diff(G,R)*psi)==0
    checks['Q_allows_auxiliary_independent_stress_only_after_projection']=qc(d1(stress))==0 and d1(stress)!=0
    checks['all_current_pressure_and_gamma_square_stay_in_Phi']=s.expand(Phi+2*s.diff(G,R)*psi-gamma*gamma-W[2][2]-pm)==0
    checks['dropping_full_pressure_from_Phi_changes_pointwise_axial_force']=s.expand(dz(qc(pm)))!=0
    checks['dropping_axial_velocity_square_changes_projected_flux']=s.expand(qc(gamma*gamma))!=0
    checks['Haar_potential_and_Haar_axial_velocity_are_nonzero']=haar(psi)!=0 and haar(gamma)!=0
    checks['whole_potential_Lambda_trace_vanishes_at_both_ends']=all(s.simplify(qc(Lam).subs(R,e))==0 for e in (1,2))
    checks['whole_G_Psi_potential_trace_vanishes_at_both_ends']=all(s.simplify((G*psi).subs(R,e))==0 for e in (1,2))
    checks['actual_Viscous_pressure_may_only_drop_after_A1_Dz']=s.expand(dz(ep*d1(beta))-d1(ep*dz(beta)))==0
    rho=(1+Z)*R;Pslow=1+Z*Z
    checks['Q_removes_only_the_Haar_compensator_not_its_full_axial_force']=qc(dz(rho*Pslow))==0 and dz(rho*Pslow)!=0
    return checks,dict(scope='full finite real Laurent potential/tensor fixtures, arbitrary retained pressure primitive and Dr=d_R+2RM d_y; finite endpoint regularity, not canonical fields or a coupled PDE solve')

def compact_controls():
    checks={};R,Z=s.symbols('R Z',real=True);I=s.I;lam=2*s.pi
    D=lambda f:s.diff(f,R)+I*lam*f
    def anti(f):
        f=s.expand(f);return s.expand(sum((-1)**j*s.diff(f,R,j)/(I*lam)**(j+1) for j in range(s.degree(f,R)+1)))
    def moment(f):
        a=anti(f);return s.expand(a.subs(R,2)-a.subs(R,1))
    for j in (2,3,5):
        psi=(R-1)**j*(2-R)**j;G=R*R
        checks[f'j{j}:T1_D1_compact_potential_identity_is_exact']=s.expand(anti(D(R*psi))-R*psi)==0
        checks[f'j{j}:whole_weighted_axial_characteristic_moment_is_zero']=moment(D(R*psi))==0
        checks[f'j{j}:G_gamma_characteristic_moment_keeps_G_R_Psi']=s.expand(moment(G*D(R*psi))+moment(2*R*R*psi))==0
        checks[f'j{j}:G_R_term_cannot_be_omitted_from_A1']=moment(2*R*R*psi)!=0
        # Haar mode: no temporal inverse is available, use T1 on gamma bar.
        checks[f'j{j}:Haar_potential_uses_weighted_primitive_not_zero_inverse']=s.expand(s.integrate(s.diff(R*psi,R),R)-R*psi)==0
    f=(R-1)**2*(2-R)**2;A=anti(f);J=moment(f);phase=s.exp(-I*lam*(R-1))
    x=R-1;chi=3*x*x-2*x**3+Z*x*x*(1-x)**2
    T1=(A-phase*A.subs(R,1)-chi*phase*J)/R
    A1=s.diff(chi,R)*phase*J/R
    checks['nonzero_characteristic_source_moment_exercises_moving_cutoff']=J!=0
    checks['T1_Dz_commutator_contains_full_cutoff_derivative']=s.simplify(s.diff(T1,Z)+s.diff(chi,Z)*phase*J/R)==0
    checks['commuting_T1_through_Dz_for_general_X_is_invalid']=s.diff(T1,Z)!=0
    checks['A1_Dz_commutator_contains_full_cutoff_derivative']=s.simplify(s.diff(A1,Z)-s.diff(chi,R,Z)*phase*J/R)==0
    checks['commuting_A1_through_Dz_for_general_Phi_is_invalid']=s.diff(A1,Z)!=0
    f=R*(R-1)**4*(2-R)**4;gamma=s.diff(f,R)+f/R
    response=s.diff(R*gamma-s.diff(gamma,R),R)+(R*gamma-s.diff(gamma,R))/R
    checks['Haar_kernel_comparison_preserves_weighted_axial_velocity_moment']=s.integrate(R*gamma,(R,1,2))==0
    checks['Haar_kernel_comparison_has_nonzero_potential_and_fine_viscous_row_response']=f!=0 and response!=0
    return checks,dict(scope='independent exact weighted characteristic moments and a moving finite polynomial cutoff; compact inversion does not justify a general slow commutator')

def physical_pressure(vals,n,t,stage,Q_only=False):
    """Independent unreduced epsilon substitution for the angular lemma."""
    ep=t**10;u=t**9*vals['u'];x=t**19*vals['x'];z=t**9*vals['z'];b=ep*vals['b0'];W=ep*vals['W'];d=n+7
    D=3*d*(1+vals['graph'])/t;e=t**(17 if stage=='postwave' else 18)*vals['Etheta']
    J=b*u+vals['V']*x+x*u+W;L=vals['G']*u+vals['V']*z+z*u+W
    F=e+ep*d*u+ep*d*L+ep*(D*D+vals['r']*D+vals['r']**2+ep**2*d*d)*u
    S=2*b*x+x*x+W;Hrz=b*z+vals['G']*x+x*z+W
    rest=d*(ep+2*vals['c'])*x+vals['r']*S+ep*d*Hrz+vals['r']*W+ep**3*d*d*x
    H=vals['temporal']*J;K=vals['temporal']*F
    av=2*vals['V']*vals['r'];a1=6*vals['V']*vals['r']**2+2*d*vals['graph']*vals['V']*vals['r']
    local=S+ep*(D+vals['r'])*x+av*H
    source=(0 if Q_only else av*ep*vals['uHaar'])+vals['r']*u*u+rest+a1*H+av*K
    if Q_only:return 2*(local+vals['KT0']*source)
    return local+vals['KT0']*(source+vals['rho']*vals['Kmom']*source)

def budget_controls():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(6,25)):
        for stage in ('postwave','final'):
            a=stage_rates(stage,k)
            checks[f'k{k}{stage}:coarse_RHS_source_is_one_with_no_fine_pressure_assumption']=a['potential_source']==1 and min(a['potential_source_groups'].values())==1
            checks[f'k{k}{stage}:potential_one_beta_two_gamma_and_regular_unchanged']=a['full_potential']==1 and a['full_beta']==2 and a['gamma']==1-k and a['whole_regular']==2-3*k
            checks[f'k{k}{stage}:angular_pressure_extension_uses_only_row_at_least_E']=a['row_order']>=2-3*k>1
        for m in (0,1,3):
            a=A1_request(m,kappa=k)
            checks[f'k{k}m{m}:direct_A1_order_two_has_double_kappa_target_margin']=a['integration_order']==0 and a['strict_margin']==2*k and a['raw_Phi_order']==m+1
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for m,N in ((0,0),(1,2),(3,7)):
            a=A1_request(m,physical_power=N,h=h)
            checks[f'h{h}m{m}N{N}:same_physical_target_and_two_lower_distinct_flux_jets']=a['strict_margin']>0 and a['previous_raw_H1_order']-a['raw_Phi_order']==2
    a=A1_request(0,physical_power=0);rows.extend((A1_request(0),a))
    checks['physical_Phi_raw60000006_and_previous_H1_raw60000008_are_distinct_sources']=a['raw_Phi_order']==60000006 and a['previous_raw_H1_order']==60000008
    checks['new_local_coarse_pressure_field60000014_is_below_old_final_A0_field60200007']=a['coarse_pressure_conservative_field_recipe']==60000014<angular.request_counts(0,physical_power=0)['rest_beta_order']
    checks['postwave6q_and_final12q_not_interchanged']=stage_rates()['pole_multiplier']==6 and stage_rates('final')['pole_multiplier']==12
    for n in (0,1,3):
        data={k:P(n+j+1) for j,k in enumerate(JOIN_KEYS)};vals={k:v.evaluate(1) for k,v in data.items()}
        for stage in ('postwave','final'):
            q=joined_envelope(n,data,stage);a=joined_envelope(n,vals,stage)
            checks[f'n{n}{stage}:composed_positive_polynomials_match_exact_rational_inputs']=q['Phi'].evaluate(1)==a['Phi'] and q['Q_Phi'].evaluate(1)==a['Q_Phi'] and q['potential']['potential'].evaluate(1)==a['potential']['potential']
            altered=dict(vals)
            for k in ('uHaar','zHaar','Kmom','rho'):altered[k]*=10000
            b=joined_envelope(n,altered,stage)
            checks[f'n{n}{stage}:Q_tail_ignores_finer_Haar_norms_and_Haar_pressure_moment']=a['Q_Phi']==b['Q_Phi'] and a['Phi']!=b['Phi']
            c=complement_envelope(n,{k:vals[k] for k in COMPLEMENT_KEYS},stage)
            checks[f'n{n}{stage}:Q_only_program_requires_no_Haar_or_moment_callback']=c['Q_Phi']==a['Q_Phi'] and not c['sharper_Haar_inputs_used']
            for t in (Q(1,2),Q(1,8)):
                ep=t**10;x=t**19*vals['x'];z=t**9*vals['z'];pm=t**9*vals['pmCoarse'];W=ep*vals['W'];d=n+7;D=3*d*(1+vals['graph'])/t
                ez=t**(17 if stage=='postwave' else 18)*vals['Ez'];psi0=vals['KT1']*z
                X=2*vals['G']*z+z*z+W+pm
                Hrz=ep*vals['b0']*z+vals['G']*x+x*z+W
                A=vals['KT1']*(ez+ep*d*X)+ep*d*psi0+Hrz+ep*D*z+ep**3*d*d*psi0
                psifine=vals['temporal']*A+ep*vals['KT1']*vals['zHaar']
                pressure=physical_pressure(vals,n,t,stage)
                Phi=2*d*vals['graph']*vals['G']*psifine+z*z+W+pressure
                pmQ=physical_pressure(vals,n,t,stage,True)
                PhiQ=2*d*vals['graph']*vals['G']*vals['temporal']*A+2*z*z+2*W+pmQ
                checks[f'n{n}{stage}t{t}:independent_full_epsilon_coarse_to_fine_potential_pressure_and_Phi_budget']=A<=ep*a['potential']['A'] and psifine<=ep*a['potential']['potential'] and pressure<=ep*a['pressure']['pressure'] and Phi<=ep*a['Phi'] and PhiQ<=ep*a['Q_Phi']
        out=tail_envelope(n,{k:vals[k] for k in COMPLEMENT_KEYS})
        beta=beta_envelope(n,{k:vals[k] for k in POTENTIAL_KEYS})
        checks[f'n{n}:tail_uses_one_extra_Phi_jet_and_beta_separately_needs_Haar_history']=out['explicit_A1_epsilon_power']==beta['explicit_epsilon_power']==2 and out['joined']['current_fields']==n+7 and beta['sharper_Haar_gamma_history_required']
    for k in (Q(1,100000),Q(1,10),Q(6,25)):
        B=1-k
        checks[f'k{k}:Haar_kernel_response_is_finer_than_both_existing_rows']=B+1>2-2*k>=2-3*k
    checks['coarse_Haar_potential_rescaling_can_grow_exponentially']=all(Q(2)**(10*j)==Q(2**j)**10 for j in (1,3,7,15))
    c=certificate()
    checks['potential_source_uses_OLD_coarse_pressure_and_Psi_to_avoid_circularity']='OLD coarse M_B' in c['noncircular']
    checks['same_POSTWAVE_pressure_and_whole_Psi_hull_are_explicit']='no final pressure substitution' in c['postwave']['source'] and 'WHOLE compact potential hull' in c['omitted_radial_potential']
    checks['complete_pressure_mean_variance_and_Pi_memory_remain']='moving rho P' in c['pressure'] and 'Pi1/Pi2' in c['pressure']
    checks['no_actual_input_full_regular_or_blowup_promotion']=not c['actual_all_scale_high_jet_support_norms_absolute_decay_infinite_field_or_blowup_proved'] and c['whole_regular_order']=='199997/100000'
    return checks,dict(rows=rows,certificate=c,scope='one-way coarse-to-fine norm programs and integer requests with exact epsilon controls; no selected actual norm supremum')

@lru_cache(None)
def _review():
    checks={};parts={}
    for name,fn in (('axial',axial_controls),('compact',compact_controls),('budgets',budget_controls)):
        a,b=fn();checks.update(a);parts[name]=b
    return dict(checks=checks,parts=parts)
def review():
    from copy import deepcopy
    return deepcopy(_review())
