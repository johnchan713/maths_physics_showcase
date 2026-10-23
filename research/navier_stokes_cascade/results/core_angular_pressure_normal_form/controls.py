"""Exact full-row normal form, independent compact moments and norm budgets."""
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
from bounds import rates,request_counts,pressure_envelope,P,mixed,certificate

def row_controls():
    checks={};R,Z,T=s.symbols('R Z T',positive=True);y=s.symbols('y',nonzero=True)
    ep,c,M=s.symbols('epsilon c M',positive=True);I=s.I
    haar=lambda f:s.expand(f).coeff(y,0)
    qc=lambda f:s.expand(f-haar(f))
    nt=lambda f:s.expand(I*y*s.diff(f,y))
    dr=lambda f:s.expand(s.diff(f,R)+2*R*M*I*y*s.diff(f,y))
    de=lambda f,e:s.expand(dr(f)+e*f/R)
    dz=lambda f:s.expand(ep*s.diff(f,Z))
    ts=lambda f:s.expand(-ep*s.diff(f,T)+c*nt(f))
    def inv(f):
        q=s.expand(f);return s.expand(sum(q.coeff(y,k)*y**k/(I*k) for k in range(-12,13) if k))
    cos=(y+y**-1)/2;cos2=(y*y+y**-2)/2;sin=(y-y**-1)/(2*I)
    # Compact polynomial fields with finite endpoint regularity only.
    bump=(R-1)**2*(2-R)**2
    # Small radial/slow degrees keep the independent exact algebra readable.
    V=R*R*(1+Z);b=ep*R;G=1+Z
    v=bump*((1+T)+(Z+T)*cos+cos2)
    psi=bump*Z*cos;beta=-dz(psi);gamma=de(psi,1)
    # Six entries of a positive rank-two real tensor, not a canonical wave.
    A=[bump*(1+cos),bump*(Z+sin),bump*(1+T*cos)]
    B=[bump*sin,bump*cos2,bump*(1+Z)]
    W=[[s.expand(A[i]*A[j]+B[i]*B[j]) for j in range(3)] for i in range(3)]
    st=bump*(1+Z)
    J=b*v+V*beta+beta*v+W[0][1]
    L=G*v+V*gamma+gamma*v+W[1][2]
    visc=dr(de(v,1))+dz(dz(v))
    E=ts(v)+de(J,2)+dz(L)-ep*visc-de(st,2)
    F=-E-ep*s.diff(v,T)+dz(L)-ep*visc
    H=inv(qc(J))/c;K=inv(qc(F))/c
    checks['full_current_angular_row_yields_exact_swirl_normal_form']=s.expand(qc(v)+de(H,2)+K)==0
    checks['no_current_pressure_or_A0_enters_angular_equation']='pressure' not in str(E.free_symbols) and len(E.free_symbols)==7
    checks['full_vector_angular_viscosity_has_correct_cylindrical_connection']=s.expand(dr(de(v,1))-(dr(dr(v))+dr(v)/R-v/R**2))==0
    checks['wrong_viscosity_sign_in_F_destroys_exact_swirl_identity']=s.expand(qc(v)+de(H,2)+inv(qc(F+2*ep*visc))/c)!=0
    checks['missing_slow_time_in_F_destroys_exact_swirl_identity']=s.expand(qc(v)+de(H,2)+inv(qc(F+ep*s.diff(v,T)))/c)!=0
    checks['auxiliary_independent_selected_stress_drops_only_after_Q']=qc(de(st,2))==0 and de(st,2)!=0
    a=2*V/R;a1=2*a/R-s.diff(a,R)
    checks['varying_V_requires_the_full_radial_multiplier_derivative']=s.expand(a1-(6*V/R**2-2*s.diff(V,R)/R))==0
    checks['weighted_radial_product_identity_includes_full_graph_drift']=s.expand(a*de(H,2)-dr(a*H)-a1*H)==0
    checks['omitting_V_derivative_is_detected']=s.expand(a*de(H,2)-dr(a*H)-2*a*H/R)!=0
    S=2*b*beta+beta*beta+W[0][0]
    flux=b*gamma+G*beta+beta*gamma+W[0][2]
    rest=-ts(beta)-S/R-dz(flux)+W[1][1]/R+ep*dz(dz(beta))
    K0=-S+ep*de(beta,1)
    gr=-ts(beta)-de(S,1)-dz(flux)+(2*V*v+v*v+W[1][1])/R+ep*(dr(de(beta,1))+dz(dz(beta)))
    r0=a*v+v*v/R+rest
    K1=K0-a*H
    r1=a*haar(v)+v*v/R+rest-a1*H-a*K
    checks['old_radial_source_reconstructed_with_both_viscosities']=s.expand(gr-r0-dr(K0))==0
    checks['full_new_radial_source_and_local_pressure_are_exact']=s.expand(gr-r1-dr(K1))==0
    checks['pressure_source_rewrite_is_exact_before_projection']=s.expand(r0-r1+dr(a*H))==0
    checks['H_K_have_zero_Haar_but_are_not_zero']=haar(H)==haar(K)==0 and H!=0 and K!=0
    checks['all_Haar_pressure_density_and_moving_P_are_unchanged']=s.expand(haar(r0)-haar(r1))==0
    checks['full_fast_variance_is_not_Haar_v_squared']=s.expand(haar(v*v)-haar(v)**2)!=0
    checks['mean_pressure_keeps_full_variance']=s.expand(haar(r1)-a*haar(v)-haar(v*v)/R-haar(rest))==0
    checks['omitting_new_local_pressure_changes_axial_pressure_force']=s.expand(dz(a*H))!=0
    checks['axial_velocity_and_pressure_flux_stay_in_radial_remainder']=gamma!=0 and s.expand(dz(flux))!=0
    checks['entire_K0_compact_trace_vanishes_at_both_hull_ends']=all(s.simplify(K0.subs(R,e))==0 for e in (1,2))
    checks['new_local_potential_has_zero_both_endpoint_traces']=all(s.simplify((a*H).subs(R,e))==0 for e in (1,2))
    return checks,dict(scope='exact full finite Laurent tensor/mean data with Dr=d_R+2RM d_y, Nt=d_y; not canonical C-infinity fields or numerical suprema')

def compact_controls():
    checks={};r=s.symbols('r',real=True);lam=2*s.pi;I=s.I
    D=lambda f:s.diff(f,r)+I*lam*f
    def anti(f):
        f=s.expand(f);return s.expand(sum((-1)**j*s.diff(f,r,j)/(I*lam)**(j+1) for j in range(s.degree(f,r)+1)))
    def moment(f):
        h=anti(f);return s.expand(h.subs(r,2)-h.subs(r,1))
    for j in (2,3,5):
        h=(r-1)**j*(2-r)**j;a=2*r;a1=2
        q=s.expand(a*h)
        checks[f'compact{j}:same_full_characteristic_moment_kills_exact_radial_derivative']=s.expand(moment(D(q)))==0
        checks[f'compact{j}:weighted_angular_flux_reduces_without_graph_derivative']=s.expand(moment(a*(D(h)+2*h/r))-moment(a1*h))==0
        checks[f'compact{j}:pressure_primitive_of_derivative_is_the_potential_not_zero']=s.expand(anti(D(q))-q)==0 and q!=0
        checks[f'compact{j}:dropping_local_pressure_would_leave_nonzero_radial_force']=D(q)!=0
    V=s.Function('V');U,d=s.symbols('U d',positive=True);R=s.symbols('R',positive=True)
    a1=6*V(R)/R**2-2*s.diff(V(R),R)/R
    # Substitute V_R=d R^(d-1) V_U before changing to U=R^d.
    V0,VU=s.symbols('V0 VU')
    kernel=s.expand(R**(1-d)/d*a1.subs({s.diff(V(R),R):d*R**(d-1)*VU,V(R):V0}))
    checks['actual_new_J_kernel_has_the_extra_background_radial_jet']=s.simplify(kernel-(6*V0/(d*R**d*R)-2*VU/R))==0
    checks['actual_F_kernel_is_the_retained_reciprocal_U_times_two_V']=s.simplify(R**(1-d)/d*(2*V0/R)-2*V0/(d*R**d))==0
    return checks,dict(scope='independent exact polynomial Fourier moments with nonzero drift; actual profile identity holds for arbitrary positive d')

def budget_controls():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(6,25)):
        a=rates(k)
        checks[f'k{k}:J_is_one_F_is_E_and_both_strictly_controlled']=a['J']==1 and a['F']==a['E'] and a['E']>1
        checks[f'k{k}:mean_pressure_improves_conditionally_but_regular_E_A1_do_not']=a['current_mean_pressure']==1 and a['whole_regular']==2-3*k and a['parent_A1']==2-k
        for m in (0,1,3):
            q=request_counts(m,kappa=k)
            checks[f'k{k}m{m}:all_four_same_target_strict_margins']=all(q[x]['strict_margin']>0 for x in ('J','F','quadratic','rest'))
            checks[f'k{k}m{m}:normalized_F_requires_two_integrations_quadratic_one']=q['F']['integrations']==2 and q['quadratic']['integrations']==1
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for m,N in ((0,0),(1,2),(3,7)):
            q=request_counts(m,physical_power=N,h=h)
            checks[f'h{h}m{m}N{N}:physical_request_keeps_every_distinct_source_and_profile']=all(q[x]['strict_margin']>0 for x in ('J','F','quadratic','rest')) and q['J_multiplier_V_order']==q['J']['raw_order']+1
    n=request_counts(0);p=request_counts(0,physical_power=0);rows.extend((n,p))
    checks['normalized_raw_orders_are_J200000_F6_square5_rest200001']=[n[x]['raw_order'] for x in ('J','F','quadratic','rest')]==[200000,6,5,200001]
    checks['physical_raw_orders_are_J60200004_F60000010_square60000009_rest60200005']=[p[x]['raw_order'] for x in ('J','F','quadratic','rest')]==[60200004,60000010,60000009,60200005]
    for q in (n,p):
        checks[f"target{q['target']}:retained_older_local_field_max_is_not_declared_reduced"]=q['rest_beta_order']==mixed.radial_counts(0,physical_power=q['physical_power'])['local_current_field_order'] and q['older_entire_DAG_retained']
    keys=('x','u','z','uHaar','b0','V','G','r','graph','c','Etheta','W','temporal','KT0','Kmom','rho')
    for n in (0,1,3):
        data={k:P(n+j+1) for j,k in enumerate(keys)};q=pressure_envelope(n,data)
        vals={k:v.evaluate(1) for k,v in data.items()};a=pressure_envelope(n,vals)
        checks[f'n{n}:positive_polynomial_pressure_program_agrees_with_exact_rational']=all(q[k].evaluate(1)==a[k] for k in ('J','F','H','K','rest','local','source','moment','pressure'))
        checks[f'n{n}:four_temporal_jets_plus_two_viscous_jets_retained']=a['current_field_input_order']==n+6 and a['current_angular_row_input_order']==n+4
        # Independently substitute all epsilon powers, including M<=eps^-k.
        for t in (Q(1,2),Q(1,8)):
            ep=t**10;u=t**9*vals['u'];x=t**19*vals['x'];z=t**9*vals['z'];bh=ep*vals['b0'];W=ep*vals['W'];e=t**18*vals['Etheta'];d=n+7;D=3*d*(1+vals['graph'])/t
            J=bh*u+vals['V']*x+x*u+W
            L=vals['G']*u+vals['V']*z+z*u+W
            F=e+ep*d*u+ep*d*L+ep*(D*D+vals['r']*D+vals['r']**2+ep**2*d*d)*u
            S=2*bh*x+x*x+W;rz=bh*z+vals['G']*x+x*z+W
            rest=d*(ep+2*vals['c'])*x+vals['r']*S+ep*d*rz+vals['r']*W+ep**3*d*d*x
            H=vals['temporal']*J;K=vals['temporal']*F
            av=2*vals['V']*vals['r'];aq=6*vals['V']*vals['r']**2+2*d*vals['graph']*vals['V']*vals['r']
            local=S+ep*(D+vals['r'])*x+av*H
            src=av*ep*vals['uHaar']+vals['r']*u*u+rest+aq*H+av*K
            pressure=local+vals['KT0']*(src+vals['rho']*vals['Kmom']*src)
            checks[f'n{n}t{t}:independent_epsilon_J1_FE_and_full_pressure1_budget']=J<=ep*a['J'] and F<=t**17*a['F'] and pressure<=ep*a['pressure']
    c=certificate()
    checks['explicit_non_circular_input_is_the_existing_angular_component_only']='no current pm or current A0' in c['noncircular']
    checks['complete_Haar_variance_local_K1_and_old_pressure_memory_kept']='Haar(v^2)' in c['Haar'] and 'K1' in c['pressure'] and 'Pi1/Pi2' in c['pressure']
    checks['actual_input_and_complete_residual_frontier_not_promoted']=not c['total_original_minimal_jet_actual_source_validity_absolute_decay_or_blowup_proved']
    return checks,dict(rows=rows,certificate=c,scope='derived norm programs and exact epsilon substitutions, not actual canonical high-jet enclosures')

@lru_cache(None)
def _review():
    checks={};parts={}
    for name,fn in (('rows',row_controls),('compact',compact_controls),('budgets',budget_controls)):
        a,b=fn();checks.update(a);parts[name]=b
    return dict(checks=checks,parts=parts)
def review():
    from copy import deepcopy
    return deepcopy(_review())
