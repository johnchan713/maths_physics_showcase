"""Exact whole-vector energy and frame controls; no continuum data sampling."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from budgets import HERE,load,normal,transfer,P,certificate,tangent_polynomial,forced_from_frame,residence_upper

identity=load('tangent_closure_identity',HERE.parent/'frame_column_foundation/frame_identity.py')


def fixture():
    p=normal.frame.inputs(bmin=1,bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),
        cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator))
    ell=normal.selector.select_tail(sp)['ell_tail'] # manufactured fixture ONLY
    kw=dict(F=lambda n:P(n+2),G=lambda n:P(n+3),btilde=lambda n:P(n+1),
        path=lambda n:P({1:1}) if n==0 else P(int(n==1)),source=lambda n:P(n+1))
    return p,sp,ell,kw


def psd2(M):
    return M[0][1]==M[1][0] and min(M[0][0],M[1][1])>=0 and M[0][0]*M[1][1]-M[0][1]**2>=0


def energy_review():
    checks={};count=0
    for lam,eta,m,E in product((Q(1,7),Q(3)),(Q(-2,3),Q(0),Q(4,5)),(-3,-1,1,2),
            ([[Q(0),Q(2,5)],[Q(-1,3),Q(0)]],[[Q(1,5),Q(2,7)],[Q(3,8),Q(-1,4)]])):
        d=Q(7,6);delta=m*m-1;bound=sum(abs(x) for row in E for x in row)
        R=[[E[i][j]+(Q(0) if i!=j else (-2*lam if i else 0)-eta-delta*d) for j in range(2)] for i in range(2)]
        upper=bound+abs(eta)-delta*d
        gap=[[(upper if i==j else 0)-(R[i][j]+R[j][i])/2 for j in range(2)] for i in range(2)]
        checks[f'case{count}_whole_complex_vector_logarithmic_energy_bound']=psd2(gap)
        count+=1
    # Negative damping deviation produces positive normalized growth even
    # with E=0. Dropping its absolute cost gives a false estimate.
    checks['omitting_damping_deviation_breaks_energy_bound']=not psd2([[Q(-1,2),Q(0)],[Q(0),Q(3,2)]])
    return dict(checks=checks,cases=count,scope='exact positive-semidefinite 2x2 matrices certify all real/complex vectors for each manufactured coefficient case')


def frame_review():
    checks={};rows=[]
    for i,(a,(s,w),c) in enumerate(product(((Q(3,5),Q(-4,5)),(Q(-4,5),Q(-3,5))),
            ((Q(0),Q(1)),(Q(3,4),Q(5,4)),(Q(-3,4),Q(5,4))), (Q(-1,2),Q(2)))):
        r=identity.matrices(a,s,w,c,Q(1),(Q(-8,5),Q(-6,5)),Q(1,50),Q(1,100))
        gram=identity.mul(identity.transpose(r['B']),r['B'])
        expected=[[w*w*(1+c*c),w*w*(1-c*c)],[w*w*(1-c*c),w*w*(1+c*c)]]
        checks[f'case{i}_exact_Gram_and_singular_value_squares']=gram==expected and gram[0][0]+gram[0][1]==2*w*w and gram[0][0]-gram[0][1]==2*c*c*w*w
        checks[f'case{i}_specified_left_inverse_restricts_to_identity']=identity.mul(r['left'],r['B'])==[[Q(1),Q(0)],[Q(0),Q(1)]]
        rows.append(dict(s=str(s),c=str(c),singular_value_squares=[str(2*w*w),str(2*c*c*w*w)]))
    # A coordinate contraction estimate alone cannot omit endpoint frames.
    f0=identity.matrices((Q(3,5),Q(-4,5)),0,1,Q(-1,2),1,(Q(-8,5),Q(-6,5)),0,0)
    f1=identity.matrices((Q(3,5),Q(-4,5)),Q(3,4),Q(5,4),Q(-1,2),1,(Q(-8,5),Q(-6,5)),0,0)
    v=[[Q(1)],[Q(1)]];t0=identity.mul(f0['B'],v);t1=identity.mul(f1['B'],v)
    norm2=lambda t:sum(row[0]**2 for row in t)
    checks['coordinate_identity_can_expand_physical_tangent_norm']=norm2(t1)/norm2(t0)==Q(25,16)>1
    return dict(checks=checks,rows=rows)


def gauge_review():
    checks={};args=((Q(3,5),Q(-4,5)),Q(3,4),Q(5,4),Q(-1,2),Q(1),(Q(-8,5),Q(-6,5)),Q(1,50),Q(1,100))
    a0=identity.matrices(*args);r0=identity.matrices(*args,mixing_w=Q(5,3),mixing_wprime=Q(1,70))
    args1=list(args);args1[1:3]=[Q(-3,4),Q(5,4)]
    a1=identity.matrices(*args1);r1=identity.matrices(*args1,mixing_w=Q(13,12),mixing_wprime=Q(1,80))
    T0=identity.mul(a0['Minv'],r0['M']);Ti1=identity.mul(r1['Minv'],a1['M'])
    Phi=[[Q(2),Q(3,5)],[Q(-1,7),Q(1,2)]]
    Phi_r=identity.mul(identity.mul(Ti1,Phi),T0)
    physical_a=identity.mul(identity.mul(a1['B'],Phi),a0['left'])
    physical_r=identity.mul(identity.mul(r1['B'],Phi_r),r0['left'])
    checks['two_time_physical_tangent_maps_are_identical_under_full_gauge']=physical_a==physical_r
    checks['zero_entrance_stays_zero_in_either_frame']=identity.mul(T0,[[Q(0)],[Q(0)]])==[[Q(0)],[Q(0)]]
    checks.update({'inherited:'+k:v for k,v in identity.gauge_review()['checks'].items()})
    return dict(checks=checks,scope='exact gauge algebra including both endpoints and moving generator; no primary or root identification')


def harmonic_review():
    I=[[Q(1),Q(0)],[Q(0),Q(1)]];A=[[Q(1),Q(1)],[Q(0),Q(1)]];B=[[Q(1),Q(0)],[Q(2),Q(1)]]
    scale=lambda M,a:[[a*x for x in row] for row in M]
    checks={'coefficient_evolutions_need_not_commute':identity.mul(A,B)!=identity.mul(B,A)}
    for m in (-3,-2,-1,1,2,3):
        delta=m*m-1;r,s=Q(1,2),Q(4,5)
        checks[f'm{m}_scalar_damping_factors_through_noncommuting_evolution']=identity.mul(scale(B,s**delta),scale(A,r**delta))==scale(identity.mul(B,A),(r*s)**delta)
        checks[f'm{m}_extra_actual_damping_never_amplifies']=(r*s)**delta<=1
    p,_,_,_=fixture()
    checks['normal_floor_supplies_positive_damping_floor']=p['bmin']**2/4>0
    checks['residence_gain_for_high_modes_not_fundamental_modes']=residence_upper(p,100,1)==100 and residence_upper(p,100,2)==Q(4,3)
    return dict(checks=checks,scope='exact scalar-factor algebra and residence bounds; higher jet damping terms are still charged by the inherited inverse')


def join_review():
    p,sp,ell,kw=fixture();c=certificate(p,sp,ell);b,f,cc=forced_from_frame(p,sp,ell,**kw)
    T=f.tangent;N=c['integer_exponential_upper_power']
    checks=dict(original_CE_plus_Cd_is_the_entire_integrated_cost=c['integrated_error']==normal.frame.budgets(p)['CE']+normal.frame.budgets(p)['Cd'],
        exact_exponential_upper_recipe_without_float=N>=c['integrated_error'] and T.c[0]==c['frame_conversion']*3**N,
        tangent_bound_is_constant_in_S=set(T.c)=={0},
        original_L_le_S_closes_whole_path_length_input=f.length.record()==P({1:1}).record(),
        original_normal_denominator_unchanged=f.inverse_normal_square==4/p['bmin']**2,
        all_original_gates_and_no_new_acceptance=all(c['gates'].values()) and not c['new_acceptance_or_smallness_condition'],
        source_callback_not_estimated_or_changed=f.source is kw['source'],
        no_actual_enclosure_primary_identity_or_residual_promotion=not c['actual_primitive_enclosures_certified'] and not c['homogeneous_primary_frames_or_roots_identified'] and not c['complete_covariance_means_or_residual_certified'])
    for j in range(3):
        # No need to stringify the large, but exact, rational coefficient.
        checks[f'order{j}_joined_velocity_and_pressure_are_finite_positive_polynomials']=f.velocity(j).evaluate(1)>0 and f.pressure(j).evaluate(1)>0
        checks[f'order{j}_full_normal_and_moving_normal_callbacks_retained']=f.normal(j).record()==b.normal(j).record() and f.moving_normal(j).record()==b.moving_normal(j).record()
    return dict(checks=checks,frame_conversion=str(c['frame_conversion']),integrated_error=str(c['integrated_error']),integer_exponent=N,
        tangent_numerator_bit_length=T.c[0].numerator.bit_length(),scope='manufactured original-rule inputs; actual primitive validity remains retained')


@lru_cache(None)
def review():
    parts=dict(energy=energy_review(),Euclidean_frame=frame_review(),canonical_gauge=gauge_review(),harmonics=harmonic_review(),inverse_join=join_review())
    return dict(checks={name+':'+k:v for name,r in parts.items() for k,v in r['checks'].items()},parts=parts)
