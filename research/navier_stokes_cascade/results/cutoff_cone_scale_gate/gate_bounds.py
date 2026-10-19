"""Correct retained collar width and affine sufficient gate; no actual scale acceptance."""
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
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def module(name,file):
    p=HERE.parent/name/file;s=importlib.util.spec_from_file_location('combined_gate_'+name,p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def first():return module('first_seminorm_bound','seminorms.py')
def cutoff():return module('cutoff_repair_derivative_budget','derivatives.py')
def ceiling(q):
    q=rat(q);return -(-q.numerator//q.denominator)
def evaluate(row,L1,L2):
    L1,L2=map(rat,(L1,L2));row=tuple(map(rat,row))
    if len(row)!=3 or min(row)<0 or min(L1,L2)<0:raise ValueError('Nonnegative affine data required')
    return row[0]+row[1]*L1+row[2]*L2

def jet_rows(core,coefficient,a,xkeep,tc,rho_root):
    """tc is the ACTUAL t1/(10 sqrt(log C)), never t1. Rows are (constant,L1,L2)."""
    coefficient,a,xkeep,tc,r=map(rat,(coefficient,a,xkeep,tc,rho_root))
    J,Jx,Jxx=(rat(core[k]) for k in ('value','radial','radial2'))
    if min(a,tc,r)<=0 or xkeep<a*a or min(J,Jx,Jxx,coefficient)<0:
        raise ValueError('Positive geometry and nonnegative core/repair bounds required')
    delta=tc/16;p=r*r;B0,B1,B2=cutoff().bump_bounds()
    S=J+coefficient*B0*max(F(1),1/r)
    x0=Jx+coefficient*max(B1/p,(B1+B0/2)/(p*r))
    z0=Jxx+coefficient*max((B2+B1)/p**2,(B2+2*B1+F(3,4)*B0)/(p*p*r))
    return dict(value=(S,F(0),F(0)),
        radial=(x0,J/(delta*xkeep),F(0)),
        radial2=(z0,2*Jx/(delta*xkeep)+J/(delta*xkeep**2),J/(delta**2*xkeep**2)))

def affine_gate(h,nu,a,c,core,coefficient,xkeep,tc,rho_root,leading,tau):
    """Same first operator and compact gate, with every primitive leading input supplied.

    leading=(F3,U3,FX2,UX2,V1norm,w0norm,CF,CU). All are norm bounds,
    not values of a sampled profile. Positive viscosity retained; actual nu=1.
    """
    h,nu,a,c,tau=map(rat,(h,nu,a,c,tau));leading=tuple(map(rat,leading))
    if len(leading)!=8 or min(leading)<0 or min(a,c,tau,nu)<=0 or c<a or not 0<h<F(1,100):
        raise ValueError('Eight nonnegative leading bounds and admissible h, viscosity, geometry/floor required')
    rows=jet_rows(core,coefficient,a,xkeep,tc,rho_root)
    src=[]
    for i in range(3):
        jets={k:v[i] for k,v in rows.items()}
        # The fixed inhomogeneous w0, CF, CU costs occur only in the constant row.
        data=leading if i==0 else (*leading[:5],F(0),F(0),F(0))
        src.append(first().source_bounds(h,nu,a,c,jets,*data))
    S=rows['value'][0];x0,x1,_=rows['radial']
    shape_rows=dict(unit=(F(1),F(0),F(0)),
        value=(5*(2*c*S+F(1,2))/a,F(0),F(0)),
        derivative=(5*(2*c*((2*c+1)*x0+S/a+1)+(2*c*S+F(1,2))/a),10*c*(2*c+1)*x1,F(0)),
        stress=tuple((c*c*src[i]['theta']/2+c*src[i]['axial'])/tau+(F(1,2) if i==0 else 0) for i in range(3)))
    return dict(jets=rows,theta=tuple(s['theta'] for s in src),axial=tuple(s['axial'] for s in src),rows=shape_rows)

def cstar(rows,L1,L2):
    L1,L2=map(rat,(L1,L2))
    if L1<1 or L2<4:raise ValueError('Necessary shape lower bounds L1>=1,L2>=4 required')
    return max(evaluate(row,L1,L2) for row in rows.values())
def shape_region(rows,budget,L1,L2):
    """Exact region for this sufficient bound at B=r*2^(2h(ell-2)); not actual cone iff."""
    budget=rat(budget)
    if budget<=0:raise ValueError('Positive error budget required')
    L1,L2=map(rat,(L1,L2));values={k:evaluate(v,L1,L2) for k,v in rows.items()}
    return dict(row_values=values,necessary_shape_lower_bounds=L1>=1 and L2>=4,
        accepted=L1>=1 and L2>=4 and all(v<=budget for v in values.values()))
def shape_caps(rows,budget,L1):
    budget,L1=map(rat,(budget,L1));d0,d1,d2=rows['derivative'];w0,w1,w2=rows['stress']
    if budget<=0 or L1<1 or d1<=0 or d2!=0 or w2<=0:raise ValueError('Nondegenerate positive-viscosity gate required')
    return dict(L1_upper=(budget-d0)/d1,L2_upper_at_L1=(budget-w0-w1*L1)/w2)
def ceil_log2_positive(q):
    """max(0,ceil(log2 q)) without floating logs, including exact powers of two."""
    q=rat(q)
    if q<=0:raise ValueError('Positive ratio required')
    if q<=1:return 0
    k=max(0,q.numerator.bit_length()-q.denominator.bit_length())
    return k+int(F(2)**k<q)
def dyadic_witness(rows,L1,L2,radius,h):
    """A sufficient later-band witness, NOT verification of an existing ell or other gates."""
    radius,h=map(rat,(radius,h))
    if radius<=0 or not 0<h<F(1,100):raise ValueError('Positive radius and admissible h required')
    C=cstar(rows,L1,L2);m=ceil_log2_positive(C/radius)
    ell=2+ceiling(F(m)/(2*h))
    return dict(Cstar=C,dyadic_exponent=m,sufficient_ell=ell,dyadic_budget=radius*F(2)**m,
        checks=dict(rational_dyadic_budget_suffices=C<=radius*F(2)**m,whole_chart_two_band_shift=2*h*(ell-2)>=m))
def compression_barrier(rows,budget):
    """One admissible compressed shape defeats any finite fixed upper-bound budget.

    This is a barrier for the sufficient certificate, not for the actual stress.
    """
    budget=rat(budget);A2=rat(rows['stress'][2])
    if budget<=0 or A2<=0:raise ValueError('Positive finite budget and quadratic-shape coefficient required')
    n=1+ceiling(budget/(16*A2))
    return dict(compression=n,stress_row_lower=rows['stress'][0]+8*n*rows['stress'][1]+16*n*n*A2,
        exceeds_fixed_budget=16*n*n*A2>budget)

def fixture():
    """Moderate exact auxiliary parameters, not selected enormous-profile samples."""
    h=F(1,200);nu=F(1);a=F(1);c=F(4);tau=F(1,10)
    core=first().collar_jets(2,7,F(1,4),a)
    axial=module('global_axial_source_envelope','envelopes.py').global_bounds(h,nu,25,a*a,7,13,101,303)
    swirl=module('global_swirl_envelope','swirl_bounds.py').swirl_inputs(a*a,c*c,7,101)
    signed=module('signed_pressure_source_moments','moments.py').from_global_axial_norm(h,nu,25,axial['U_C3'])
    mu=first().moment_bounds(2,core['value'],swirl['F_C3'],axial['U_C3'],(signed['I0'],signed['I1']))
    repair=first().repair_bound(mu,2,3,F(1,1000))
    leading=(swirl['F_C3'],axial['U_C3'],swirl['F_X_C2'],axial['U_X_C2'],axial['V_C1'],axial['w0_C1'],swirl['shifted_F_C1'],axial['shifted_U_C1'])
    return dict(h=h,nu=nu,a=a,c=c,core=core,coefficient=repair['coefficient'],xkeep=F(9,8),
        tc=F(1,160),rho_root=F(2),leading=leading,tau=tau)

def symbolic_review():
    import sympy as s
    J,Jx,d,x,c,nu,tau,t1,g,L1,L2=s.symbols('J Jx d x c nu tau t1 g L1 L2',positive=True)
    Sxx=2*Jx*L1/(d*x)+J*(L2/d**2+L1/d)/x**2
    A2=nu*c**3*(c+2)*J/(tau*d**2*x**2)
    stress_shape=(c*c/2+c)*2*nu*c*c*Sxx/tau
    z,lnC=s.symbols('z lnC',positive=True)
    true_delta=t1/(160*s.sqrt(lnC))
    old_alpha=L2/(t1/16)**2/x**2;old_beta=L1/(t1/16)/x**2
    corrected=(L2/(t1/(16*g))**2+L1/(t1/(16*g)))/x**2
    return dict(stress_L2_coefficient=str(A2),checks=dict(
        retained_cutoff_log_width_is_tc_over_sixteen=s.simplify(true_delta-(t1/(10*s.sqrt(lnC)))/16)==0,
        corrected_first_derivative_cost_has_ten_sqrt_logC=s.simplify(1/true_delta/(16/t1)-10*s.sqrt(lnC))==0,
        corrected_second_cost_has_distinct_g_squared_and_g_terms=s.simplify(corrected-g*g*old_alpha-g*old_beta)==0,
        exact_positive_L2_stress_coefficient=s.simplify(s.diff(stress_shape,L2)-A2)==0,
        no_L1_L2_product_or_quadratic_shape_terms=s.diff(stress_shape,L1,L2)==0 and s.diff(stress_shape,L1,2)==0 and s.diff(stress_shape,L2,2)==0,
        missing_viscous_shape_cost_detected=s.simplify(A2.subs(nu,0))==0 and A2.is_positive is True))

def review():
    p=fixture();out=affine_gate(**p);rows=out['rows'];controls=[]
    for L1,L2 in product((F(1),F(8),F(32),F(8000)),(F(4),F(1216),F(16000000))):
        ext=cutoff().propagate(p['core'],p['coefficient'],p['a'],p['xkeep'],p['tc'],L1,L2,p['rho_root'])['extended']
        src=first().source_bounds(p['h'],p['nu'],p['a'],p['c'],ext,*p['leading'])
        direct=first().compact_gate(p['a'],p['c'],ext,src,p['tau'])['Cval']
        controls.append(dict(L1=str(L1),L2=str(L2),Cstar=str(direct),checks=dict(
            exact_jet_affine_identity=ext=={k:evaluate(v,L1,L2) for k,v in out['jets'].items()},
            exact_full_theta_and_axial_affine_identity=src['theta']==evaluate(out['theta'],L1,L2) and src['axial']==evaluate(out['axial'],L1,L2),
            exact_max_affine_compact_gate_identity=direct==cstar(rows,L1,L2))))
    # All stencil inputs here are auxiliary; do not identify these budgets with the old kappa.
    C=cstar(rows,32,1216);budget=2*C;caps=shape_caps(rows,budget,32)
    boundary=(budget-rows['stress'][0]-32*rows['stress'][1])/rows['stress'][2]
    accept=shape_region(rows,budget,32,boundary);reject=shape_region(rows,budget,32,boundary+1)
    corner=cstar(rows,1,4);barrier=compression_barrier(rows,budget)
    leading=p['leading'];erased=dict(p,leading=(*leading[:5],F(0),F(0),F(0)))
    no_fixed=affine_gate(**erased)
    signed=module('signed_pressure_source_moments','moments.py').from_global_axial_norm(p['h'],p['nu'],25,leading[1])
    full_mu=first().moment_bounds(2,p['core']['value'],leading[0],leading[1],(signed['I0'],signed['I1']))
    erased_mu=first().moment_bounds(2,p['core']['value'],leading[0],leading[1],(0,0))
    full_rep=first().repair_bound(full_mu,2,3,F(1,1000));erased_rep=first().repair_bound(erased_mu,2,3,F(1,1000))
    # Keep identical moments/radii here, isolate the missing width cost alone.
    wrong=affine_gate(**dict(p,tc=10*p['tc']))
    width=[]
    for g in (F(2),F(10),F(20),F(40)):
        old=cutoff().cutoff_bounds(p['xkeep'],g*p['tc']/16,32,1216)
        true=cutoff().cutoff_bounds(p['xkeep'],p['tc']/16,32,1216)
        width.append(dict(g=str(g),checks=dict(first_cost_exactly_g=true['radial']==g*old['radial'],
            second_cost_between_g_and_g_squared=g*old['radial2']<true['radial2']<g*g*old['radial2'])))
    tilt=module('leading_cone_margin','cone_bounds.py').tilt_gate(F(1,4),4)
    gap=module('background_gap_certificate','gap_bounds.py').certificate(F(1,2),2,2,5,4,F(1,2),tilt['gamma'],2)
    witness=dyadic_witness(rows,32,1216,gap['radius'],p['h'])
    sym=symbolic_review()
    return dict(auxiliary_rational_not_actual_profile=True,conditional_fixture=serial(p),affine_output=serial(out),
        exact_composition_controls=controls,width_correction_controls=width,symbolic=sym,
        shape_region=dict(budget=str(budget),caps_at_L1_32=serial(caps),stress_boundary_L2=str(boundary),
            boundary_accepted=accept['accepted'],strictly_larger_L2_rejected=not reject['accepted']),
        complete_memory_rows=serial(dict(full_moments=full_mu,erased_moments=erased_mu,
            full_normalized=full_rep['transformed'],erased_normalized=erased_rep['transformed'])),
        fixed_budget_compression_obstruction=serial(barrier),same_tilt=serial(tilt),gap_radius=str(gap['radius']),dyadic_witness=serial(witness),
        checks=dict(all_twelve_exact_complete_operator_compositions=all(all(r['checks'].values()) for r in controls),
            correct_retained_width_chain_and_viscous_shape_coefficient=all(sym['checks'].values()) and all(all(r['checks'].values()) for r in width),
            all_affine_coefficients_nonnegative=all(min(r)>=0 for r in rows.values()),
            strictly_positive_L2_stress_cost=rows['stress'][2]>0,
            complete_fixed_pressure_swirl_and_axial_costs_retained=rows['stress'][0]>no_fixed['rows']['stress'][0] and rows['stress'][1:]==no_fixed['rows']['stress'][1:],
            future_signed_pressure_memory_enters_both_normalized_repair_rows=full_rep['coefficient']==p['coefficient'] and full_rep['transformed'][3]>erased_rep['transformed'][3] and full_rep['transformed'][4]>erased_rep['transformed'][4],
            wider_mislabelled_tc_underestimates_both_shape_slopes=rows['stress'][1]>wrong['rows']['stress'][1] and rows['stress'][2]==100*wrong['rows']['stress'][2],
            exact_closed_shape_region_and_strict_outside=accept['accepted'] and not reject['accepted'],
            lower_shape_corner_is_only_necessary=not shape_region(rows,corner/2,1,4)['accepted'],
            compressed_admissible_shapes_defeat_every_finite_fixed_certificate=barrier['exceeds_fixed_budget'],
            same_tilt_positive_radius_and_two_band_dyadic_witness=tilt['accepted'] and gap['radius']>0 and all(witness['checks'].values())))
