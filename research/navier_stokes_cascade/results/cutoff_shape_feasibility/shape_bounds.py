"""Sharp smooth-cutoff budget feasibility; auxiliary existence is not the retained cutoff."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import importlib.util
import sys

HERE=Path(__file__).resolve().parent
def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def gate_module():
    p=HERE.parent/'cutoff_cone_scale_gate/gate_bounds.py'
    s=importlib.util.spec_from_file_location('shape_retained_affine_gate',p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m

def capacity(M,N):
    """Integral of min(M,N*s,N*(1-s)), on the unchanged normalized interval."""
    M,N=map(rat,(M,N))
    if min(M,N)<=0:raise ValueError('Positive finite derivative enclosures required')
    if 2*M<=N:return M-M*M/N
    return N/4
def shape_budget(M,N):
    """Iff existence of SOME C-infinity cutoff with these upper enclosures, not the old K."""
    M,N=map(rat,(M,N));A=capacity(M,N)
    return dict(capacity=A,smooth_shape_exists=A>1,clipped_slope_saturates_height=2*M<=N)
def smoothing_witness(M,N):
    """Bounds for the analytical shrink/mollify/normalize construction, not profile samples.

    V_eps=(1-2eps)*W((s-eps)/(1-2eps)), then convolve by any nonnegative
    unit-mass C-infinity kernel of support [-eps/2,eps/2], and normalize mass.
    """
    M,N=map(rat,(M,N));A=capacity(M,N)
    if A<=1:raise ValueError('Strict slope capacity above one required for smooth existence')
    eps=(A-1)/(16*A);mass=(1-2*eps)**2*A
    return dict(epsilon=eps,pre_normalization_mass=mass,
        first_norm_upper=(1-2*eps)*M/mass,second_norm_upper=N/mass,
        slope_support=(eps/2,1-eps/2),checks=dict(
            positive_padding_inside_original_interval=0<eps<F(1,16),
            mass_stays_strictly_above_one=mass>1,
            first_norm_strictly_below_supplied_enclosure=(1-2*eps)*M/mass<M,
            second_norm_strictly_below_supplied_enclosure=N/mass<N,
            normalized_enclosures_have_strict_capacity=capacity((1-2*eps)*M/mass,N/mass)>1))

def normalized_rows(rows):
    if set(rows)!=set(('unit','value','derivative','stress')):raise ValueError('The four retained affine rows required')
    rows={k:tuple(map(rat,v)) for k,v in rows.items()}
    if any(len(v)!=3 or min(v)<0 for v in rows.values()):raise ValueError('Three nonnegative affine coefficients required')
    if any(rows[k][1:]!=(0,0) for k in ('unit','value')) or rows['derivative'][1]<=0 or rows['derivative'][2]!=0 or rows['stress'][2]<=0:
        raise ValueError('Constant value rows, positive first slope and positive viscous second slope required')
    return rows
def affine_value(row,M,N):return row[0]+row[1]*M+row[2]*N
def shape_gate(rows,B):
    """Eliminate derivative budgets for existence in THIS sufficient bound only.

    Every accepted result contains rational upper enclosures for an auxiliary
    smooth shape. They never identify the unchanged cutoff or certify other gates.
    """
    rows=normalized_rows(rows);B=rat(B)
    if B<=0:raise ValueError('Positive fixed scale error budget required')
    if max(rows['unit'][0],rows['value'][0])>B:
        return dict(possible=False,reason='constant_row_exceeds_budget')
    D0,D1,_=rows['derivative'];A0,A1,A2=rows['stress']
    cap=min(F(1),(B-D0)/D1-1);beta=B-A0-A1-2*A2;P=A1+A2
    if cap<=0:return dict(possible=False,reason='first_shape_budget_at_most_one',t_cap=cap,beta=beta)
    if beta<=0:return dict(possible=False,reason='nonpositive_stress_surplus',t_cap=cap,beta=beta)
    t=min(cap,beta/(2*P));q=P*t*t-beta*t+A2
    base=dict(t_cap=cap,beta=beta,quadratic_minimizer=t,quadratic_minimum=q)
    if q>=0:return dict(possible=False,reason='strict_joint_shape_quadratic_not_negative',**base)
    M=1+t;Nmin=M*M/t;Nmax=(B-A0-A1*M)/A2;N=(Nmin+Nmax)/2
    vals={k:affine_value(v,M,N) for k,v in rows.items()}
    return dict(possible=True,reason='auxiliary_smooth_shape_budget_exists_not_retained_shape_verified',**base,
        M=M,N=N,N_frontier=Nmin,N_stress_cap=Nmax,row_values=vals,
        smoothing=smoothing_witness(M,N),checks=dict(
            first_budget_between_one_and_two=1<M<=2,
            strict_frontier_slack=Nmin<N<Nmax,
            capacity_strictly_above_one=capacity(M,N)>1,
            all_four_original_rows_within_budget=all(v<=B for v in vals.values())))

def stress_floor_test(A0,A1,A2,B):
    """Necessary strict B>A0+A1+2A2+2sqrt(A2*(A1+A2)), by rational squaring."""
    A0,A1,A2,B=map(rat,(A0,A1,A2,B))
    if min(A0,A1)<0 or A2<=0 or B<=0:raise ValueError('Nonnegative constant/first slope and positive viscosity/budget required')
    beta=B-A0-A1-2*A2;product=A2*(A1+A2)
    return dict(shift=beta,square_threshold=4*product,
        passes_necessary_strict_stress_floor=beta>0 and beta*beta>4*product)

def symbolic_review():
    import sympy as s
    z,M,N,P,Q,beta,t,t0,A1,A2=s.symbols('z M N P Q beta t t0 A1 A2',positive=True)
    area=2*s.integrate(N*z,(z,0,M/N))+M*(1-2*M/N)
    triangle=2*s.integrate(N*z,(z,0,s.Rational(1,2)))
    Fcap=M-M*M/N
    h=(A1+A2)*t+A2/t
    q=P*t*t-beta*t+Q;vertex=beta/(2*P)
    return dict(checks=dict(
        exact_clipped_trapezoid_integral=s.simplify(area-Fcap)==0,
        exact_triangular_integral=s.simplify(triangle-N/4)==0,
        branches_agree_at_M_equals_N_over_two=s.simplify(Fcap.subs(M,N/2)-N/4)==0,
        strict_frontier_from_area_one=s.simplify(Fcap.subs(N,M*M/(M-1))-1)==0,
        t_parameter_frontier=s.simplify(((M*M/(M-1)).subs(M,1+t))-(2+t+1/t))==0,
        square_complete_for_quadratic_minimum=s.simplify(q-(P*(t-vertex)**2+Q-beta**2/(4*P)))==0,
        capped_minimum_difference=s.simplify(q-q.subs(t,t0)-(t0-t)*(beta-P*(t+t0)))==0,
        weighted_infimum_at_stationary_t=s.simplify(h.subs(t,s.sqrt(A2/(A1+A2)))-2*s.sqrt(A2*(A1+A2)))==0,
        higher_M_comparison_square=s.expand((A1+2*A2)**2-4*A2*(A1+A2)-A1**2)==0,
        rectangular_lower_corner_capacity_is_three_quarters=s.simplify(Fcap.subs({M:1,N:4})-s.Rational(3,4))==0))

def capacity_review():
    # Independent exact trapezoid integration, no sampling/finite differencing.
    rows=[]
    for M,N in product((F(1),F(3,2),F(2),F(8),F(32)),(F(3),F(4),F(9,2),F(8),F(1216))):
        d=min(M/N,F(1,2));integral=N*d*d+M*max(F(0),1-2*d)
        # When the height cap is inactive, the middle interval has zero length.
        A=capacity(M,N);r=shape_budget(M,N)
        rows.append(dict(M=str(M),N=str(N),capacity=str(A),exists=r['smooth_shape_exists'],checks=dict(
            independent_clipped_linear_integration=integral==A,
            no_M_at_most_one_can_work=M>1 or A<=1,
            no_N_at_most_four_can_work=N>4 or A<=1)))
    near=[]
    for k in (1,2,8,30):
        e=F(1,2**k);M=F(3,2);N0=M*M/(M-1)
        lo,at,hi=(shape_budget(M,N0+d) for d in (-e,F(0),e))
        w=smoothing_witness(M,N0+e)
        near.append(dict(k=k,padding=serial(w),checks=dict(
            both_sides_and_exact_frontier_strict=not lo['smooth_shape_exists'] and not at['smooth_shape_exists'] and hi['smooth_shape_exists'],
            analytic_smoothing_padding_keeps_both_norms=all(w['checks'].values()))))
    return dict(exact_rows=rows,near_frontier=near,checks=dict(
        all_twenty_five_exact_capacity_controls=all(all(r['checks'].values()) for r in rows),
        all_four_strict_frontier_and_smoothing_controls=all(all(r['checks'].values()) for r in near),
        rectangular_lower_corner_unrealizable=capacity(1,4)==F(3,4),
        loose_first_upper_bound_does_not_force_second_above_its_square=capacity(1000,1216)>1 and 1216<1000**2))

def gate_review():
    # Toy slopes deliberately make the old (1,4) corner insufficient.
    rows=dict(unit=(F(1),F(0),F(0)),value=(F(1),F(0),F(0)),
        derivative=(F(0),F(1),F(0)),stress=(F(0),F(3),F(1)))
    g=gate_module();controls=[]
    for B in (F(7),F(8),F(9),F(10),F(20)):
        r=shape_gate(rows,B)
        controls.append(dict(B=str(B),output=serial(r),old_rectangular_corner_passes=g.shape_region(rows,B,1,4)['accepted']))
    capped=dict(rows,derivative=(F(0),F(3),F(0)),stress=(F(0),F(0),F(1)))
    boundary=shape_gate(capped,F(9,2));after=shape_gate(capped,F(9,2)+F(1,100))
    p=g.fixture();full=g.affine_gate(**p)['rows'];B=2*g.cstar(full,32,1216)
    joined=shape_gate(full,B);erased=dict(p,leading=(*p['leading'][:5],F(0),F(0),F(0)))
    no_memory=g.affine_gate(**erased)['rows']
    floor=stress_floor_test(0,3,1,9);above=stress_floor_test(0,3,1,F(9)+F(1,1000))
    return dict(toy_rows=serial(rows),exact_budget_controls=controls,capped_boundary=serial(boundary),capped_above=serial(after),
        full_inherited_rows=serial(full),full_inherited_budget=str(B),full_inherited_join=serial(joined),
        checks=dict(old_corner_pass_can_coexist_with_no_smooth_shape_certificate=controls[1]['old_rectangular_corner_passes'] and not controls[1]['output']['possible'],
            exact_unattainable_stress_infimum_rejected=not controls[2]['output']['possible'] and not floor['passes_necessary_strict_stress_floor'],
            strict_above_infimum_can_work=controls[3]['output']['possible'] and above['passes_necessary_strict_stress_floor'],
            derivative_cap_equality_has_no_smooth_shape=not boundary['possible'] and boundary['quadratic_minimum']==0,
            derivative_cap_strict_surplus_has_smooth_shape=after['possible'] and all(after['checks'].values()),
            complete_corrected_width_pressure_and_viscosity_rows_composed=joined['possible'] and all(joined['checks'].values()) and full['stress'][0]>no_memory['stress'][0],
            no_source_term_or_old_field_replaced_by_elimination=full==g.affine_gate(**p)['rows']))

def weighted_review():
    rows=[]
    for A1,A2,t in ((F(0),F(1),F(1)),(F(3),F(1),F(1,2)),(F(8),F(1),F(1,3)),(F(5),F(4),F(2,3)),(F(15),F(1),F(1,4))):
        M=1+t;N0=M*M/t;inf=A1+2*A2+2*A2/t
        for eps in (F(1,2),F(1,100),F(1,2**20)):
            w=smoothing_witness(M,N0+eps);upper=A1*w['first_norm_upper']+A2*w['second_norm_upper']
            rows.append(dict(A1=str(A1),A2=str(A2),t=str(t),infimum=str(inf),epsilon=str(eps),
                smoothed_cost_upper=str(upper),checks=dict(
                    exact_stationary_ratio=(A1+A2)*t*t==A2,
                    exact_frontier_cost=A1*M+A2*N0==inf,
                    smoothed_cost_fits_arbitrarily_close_upper_budget=inf<upper<inf+A2*eps,
                    exact_infimum_boundary_not_accepted=not stress_floor_test(0,A1,A2,inf)['passes_necessary_strict_stress_floor'],
                    every_positive_cost_slack_passes_floor=stress_floor_test(0,A1,A2,inf+eps)['passes_necessary_strict_stress_floor'])))
    return dict(auxiliary_infimum_approximation_not_retained_shape=True,rows=rows,
        checks=dict(all_fifteen_exact_weighted_infimum_controls=all(all(r['checks'].values()) for r in rows)))

def review():
    s,c,g,w=symbolic_review(),capacity_review(),gate_review(),weighted_review()
    return dict(auxiliary_shapes_not_actual_cutoff=True,symbolic=s,capacity=c,gate=g,weighted=w,
        checks=dict(symbolic_capacity_frontier_minimization_and_weighted_floor=all(s['checks'].values()),
            exact_joint_capacity_and_strict_smoothing_controls=all(c['checks'].values()),
            exact_shape_elimination_and_full_retained_affine_composition=all(g['checks'].values()),
            sharp_unattained_weighted_stress_floor_controls=all(w['checks'].values())))
