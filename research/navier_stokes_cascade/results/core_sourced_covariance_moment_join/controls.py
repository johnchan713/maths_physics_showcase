"""Exact complete-tensor, curl, moving-moment and edge-domain controls."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial
from budgets import HERE,load,tangent,prefix,P,GeneratedWave,curl_budget,moment_budget,generated_moments,generated_residual_moments,next_pole_ledger,operator_coverage,orders

means=load('generated_moment_exact',HERE.parent/'complete_covariance_initial_means/means.py')
R,Z,T,U=means.R,means.Z,means.T,means.U


def fixture(trace=None):
    p=tangent.normal.frame.inputs(bmin=1,bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=tangent.normal.selector.fixture();L=tangent.normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator));ell=tangent.normal.selector.select_tail(sp)['ell_tail'] # manufactured ONLY
    pp=prefix.P
    b=prefix.PrefixBudget(fixed=lambda n:{k:pp(i+n+1) for i,k in enumerate(prefix.FIXED)},
        initial=lambda n:{k:pp(i+n+1) for i,k in enumerate(prefix.STATE)},
        wave=lambda n:{k:pp(i+n+1) for i,k in enumerate(('T','R','Pi','U'))},
        linear=lambda n:{k:pp(i+n+1) for i,k in enumerate(prefix.LINEAR)},
        KT1=lambda n:pp(n+2),KA1=lambda n:pp(n+1),time_inverse=lambda n:6*pp({1:1}),
        slow=lambda n:{k:pp(i+n+1) for i,k in enumerate(('H','Kp','Kz','Kg','Kpsi'))},trace=trace)
    w=GeneratedWave(p,sp,ell,b,F=lambda n:P(n+2),G=lambda n:P(n+3),btilde=lambda n:P(n+1),
                    path=lambda n:P({1:1}) if not n else P(int(n==1)))
    return w


def covariance_review():
    old={1:([Q(1),Q(2),Q(3)],[Q(2),Q(-1),Q(4)])}
    inc={1:([Q(2),Q(1),Q(-1)],[Q(1),Q(3),Q(2)]),2:([Q(3),Q(-2),Q(1)],[Q(-1),Q(2),Q(4)])}
    def cov(data):return [[sum((a[i]*a[j]+b[i]*b[j] for a,b in data.values()),Q(0))/2 for j in range(3)] for i in range(3)]
    new={m:tuple([[v+z for v,z in zip(x,y)] for x,y in zip(old.get(m,([0]*3,[0]*3)),inc[m])]) for m in inc}
    W0,W1=cov(old),cov(new)
    delta=[[W1[i][j]-W0[i][j] for j in range(3)] for i in range(3)]
    exact=[[sum((old.get(m,([0]*3,[0]*3))[0][i]*a[j]+a[i]*old.get(m,([0]*3,[0]*3))[0][j]
                   +old.get(m,([0]*3,[0]*3))[1][i]*b[j]+b[i]*old.get(m,([0]*3,[0]*3))[1][j]+a[i]*a[j]+b[i]*b[j]
                  for m,(a,b) in inc.items()),Q(0))/2 for j in range(3)] for i in range(3)]
    checks=dict(complete_tensor_matches_difference_of_two_covariances=delta==exact,
        all_six_entries_change=all(delta[i][j]!=0 for i in range(3) for j in range(i,3)),
        absent_old_second_mode_still_has_nonzero_square=all(v!=0 for row in cov({2:inc[2]}) for v in row))
    old_l1=sum(abs(x) for pair in old.values() for v in pair for x in v)
    new_l1=sum(abs(x) for pair in inc.values() for v in pair for x in v)
    checks['full_tensor_each_entry_has_l1_cross_plus_square_bound']=all(abs(x)<=2*old_l1*new_l1+new_l1**2 for row in delta for x in row)
    return dict(checks=checks,delta=[[str(x) for x in row] for row in delta],scope='exact two-frequency real quadratures, not actual sourced ODE samples')


def poly_norm(f,n):
    total=Q(0)
    for alpha in product(range(n+1),repeat=4):
        if sum(alpha)>n:continue
        df=f
        for axis,j in enumerate(alpha):
            for _ in range(j):df=df.diff(axis)
        denominator=1
        for j in alpha:denominator*=factorial(j)
        total+=sum((abs(c)*max(Q(1),Q(2)**r) for (r,z,t,u),c in df.terms.items()),Q(0))/denominator
    return total


def curl_review():
    A=[R*Z+U*R,R**2*T+U*R,R**3+Z*T+U*R**2];psi=(R-1)*(2-R)*(1+Z)
    cut=[psi*a for a in A];eps=Q(1,200);dr=lambda f:f.diff(0)+R**2*f.diff(3)
    rem=[-eps*cut[1].diff(1),eps*cut[0].diff(1)-dr(cut[2]),dr(cut[1])+cut[1].shift(-1)]
    checks={};rows=[]
    for n in range(3):
        pot=sum(poly_norm(f,n) for f in cut);pn=sum(poly_norm(f,n+1) for f in cut)
        graph=1+poly_norm(R**2,n);recip=sum(Q(1) for _ in range(n+1))
        bound=curl_budget(n,Q(0),pot,pn,graph,recip)['R'];actual=sum(poly_norm(f,n) for f in rem)
        checks[f'order{n}_full_graph_cylindrical_curl_norm_is_bounded']=actual<=bound
        rows.append(dict(order=n,actual_majorant=str(actual),upper=str(bound)))
    checks['cutoff_derivatives_cannot_be_erased']=any(dr(psi*a)!=psi*dr(a) for a in A)
    checks['graph_and_cylindrical_connection_are_nonzero']=dr(cut[2])!=cut[2].diff(0) and bool(cut[1].shift(-1))
    return dict(checks=checks,rows=rows,scope='independent full-box polynomial coefficient norms; finite fixture, not a continuum pulse enclosure')


def moment_review():
    f=(R-1)**4*(2-R)**4
    Drr=f*(1+Z+T);Dtt=2*f*(1+Z**2);Drz=f*(1+2*Z+T);Dtz=3*f*(1+Z*T);Dzz=4*f*(1+Z+Z**2)
    eps=Q(1,200);dg=-means.radial(Drr)-eps*Drz.diff(1)+Dtt.shift(-1)
    # R^-1 is outside the fixture integrator: use an additional R factor.
    Drr,Dtt,Drz,Dtz,Dzz=(R*x for x in (Drr,Dtt,Drz,Dtz,Dzz))
    dg=-means.radial(Drr)-eps*Drz.diff(1)+Dtt.shift(-1)
    Pd=dg.integral();Jt=(R**2*Dtz).integral();Jz=(R*Dzz).integral()-(R**2*dg).integral()/2
    rho=f/f.integral().scalar()+Z*(R-Q(3,2))*f/f.integral().scalar();cr=(R**2*rho).integral()/2
    pm=(dg-rho*Pd).primitive();pm=pm-pm.at(0,1)
    checks=dict(full_axis_pressure_defect_identity=Pd==((Dtt-Drr).shift(-1)).integral()-eps*Drz.integral().diff(1),
        full_axial_covariance_and_pressure_moment=Jz==(R*(Dzz-(Drr+Dtt)/2)).integral()+eps*(R**2*Drz).integral().diff(1)/2,
        compact_pressure_both_endpoints=pm.at(0,1)==pm.at(0,2)==0,
        full_pressure_first_moment=(R*pm).integral()==-(R**2*dg).integral()/2+cr*Pd,
        whole_moving_pressure_derivative=((R*(Dzz+pm)).integral()).diff(1)==(Jz+cr*Pd).diff(1),
        moving_compensator_derivative_survives=bool(cr.diff(1)*Pd),
        axial_pressure_not_supplied_by_two_selected_fluxes=bool(Dzz) and bool(Drz) and bool(dg))
    for n in range(3):
        V=max(poly_norm(x,n) for x in (Drr,Dtt,Drz,Dtz,Dzz));Vn=max(poly_norm(x,n+1) for x in (Drr,Dtt,Drz,Dtz,Dzz))
        b=moment_budget(n,V,Vn,1,2,poly_norm(cr,n))
        for name,field in dict(P=Pd,Jtheta=Jt,Jz=Jz,pressure_moment=(R*pm).integral(),whole_axial_moment=Jz+cr*Pd).items():
            checks[f'order{n}_{name}_independent_polynomial_majorant_is_bounded']=poly_norm(field,n)<=b[name]
    return dict(checks=checks,scope='exact compact finite-regularity polynomial moment/pressure controls, with nonconstant compensator')


def pole_review():
    checks={};rows=[]
    for p in (Q(0),Q(1),Q(3,2)):
        b=next_pole_ledger(p);old=dict(KT1=[2*p],KA1=[2*p],Kp=[2*p,4*p],Kmom=[2*p,4*p]);new=b['required_operator_poles']
        checks[f'p{p}_sourced_self_square_and_fresh_map_poles']=b['covariance']==6*p and b['next_mean_ledger']['fresh_post_fast_rows_and_moments']==12*p
        checks[f'p{p}_correctly_recomputed_domains_cover_next_maps']=operator_coverage(p,new)['pass_all']
        checks[f'p{p}_old_domains_are_sufficient_only_at_zero_pole']=operator_coverage(p,old)['pass_all']==(p==0)
        rows.append(dict(pole=str(p),required={k:list(map(str,v)) for k,v in new.items()}))
    ratios=[]
    for y in (Q(1,2),Q(1,4),Q(1,8)):
        ratios.append(y**-2) # p=1: (zeta*y^-6)/(zeta*y^-4)
    checks['old_4p_input_norm_cannot_bound_generic_6p_input']=ratios==[4,16,64]
    recipes=[prefix.haar.edge.primitive_recipe(p,4,Q(1,2)) for p in (6,12)]
    checks['same_fixed_width_primitive_recipe_exists_at_both_required_poles']=all(r['width']=='1/2' for r in recipes) and recipes[0]['kernel_degree']<recipes[1]['kernel_degree']
    return dict(checks=checks,rows=rows,scope='proof-domain correction, not an obstruction to existence of the fixed operator at a higher finite pole')


def generated_review():
    trace=set();w=fixture(trace);checks={};rows=[]
    for n in (0,1):
        a=w.addition(n);c=w.covariance(n);old=P(w.before.wave(n)['U'].c)
        checks[f'order{n}_all_signed_modes_and_complete_cut_curl_retained']=a['U'].c==sum((w.mode(n,m)['U'] for m in w.harmonics),P()).c and a['U'].evaluate(1)>a['T'].evaluate(1)>0
        checks[f'order{n}_generated_covariance_cross_and_square_retained']=c.c==(2*old*a['U']+a['U']*a['U']).c
        checks[f'order{n}_full_cut_pressure_keeps_unprojected_source']=a['Pi'].evaluate(1)>0
        rows.append(dict(order=n,complete_degree=max(a['U'].c),covariance_degree=max(c.c)))
    c_rho=lambda n:P(n+2)
    m=generated_moments(w,0,1,Q(1,2),c_rho);r=generated_residual_moments(w,0,1,Q(1,2),c_rho)
    checks['all_three_generated_moments_and_whole_pressure_response_are_positive']=all(v.evaluate(1)>0 for v in m.values()) and all(v.evaluate(1)>0 for v in r.values())
    highest=max(n for name,n in trace if name=='initial')
    checks['residual_order_zero_uses_finite_before_prefix_jets_through_18']=highest==18
    checks['normal_denominator_and_derived_tangent_still_used']=all(f.inverse_normal_square==4 for f in w.modes.values()) and w.certificate['tangent_and_length_inputs_derived_from_retained_frame_bounds']
    checks['no_actual_enclosure_or_complete_residual_promotion']=not w.certificate['actual_primitive_enclosures_certified'] and not w.certificate['complete_covariance_means_or_residual_certified']
    for k in (Q(1,100000),Q(1,1000),Q(1,100)):
        o=orders(k);checks[f'k{k}_all_power_margins_and_extra_axial_epsilon_retained']=o['self_square_slack']>0 and o['curl_slack']>0 and o['residual_moments']==o['covariance']+1
    return dict(checks=checks,rows=rows,highest_initial_jet=highest,scope='exact positive-polynomial composition on manufactured retained inputs; no actual generated-fluid bounds certified')


@lru_cache(None)
def review():
    parts=dict(tensor=covariance_review(),complete_curl=curl_review(),pressure_moments=moment_review(),weighted_domains=pole_review(),generated_chain=generated_review())
    return dict(checks={name+':'+k:v for name,r in parts.items() for k,v in r['checks'].items()},parts=parts)
