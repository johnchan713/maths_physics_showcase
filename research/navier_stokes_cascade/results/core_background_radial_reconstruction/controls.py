"""Exact reconstruction jets, incoming-memory controls and lower-norm joins."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial
from budgets import HERE,load,reconstruct,fixture,finalized_first,low_comparison,head


@lru_cache(None)
def jets():return load('radial_reconstruction_jets',HERE.parent/'modulation_repair_all_order/map_checks.py').JetClass()


def box_raw(poly,r,k,Xhi):
    out=Q(0)
    for (a,b),c in poly.terms.items():
        if a>=r and b>=k:
            out+=abs(c)*Q(factorial(a),factorial(a-r))*Q(factorial(b),factorial(b-k))*Xhi**(a-r)
    return out


def exact_review():
    J=jets();x,z=[J.variable(i,0,n=2,order=9) for i in range(2)]
    F=1+x*z+x*x;U=1+x+x*x*z*z+x*z**3
    checks={};rows=[];h=Q(1,200);D=Q(1,2)-h;A=Q(1,2)+h
    for m in range(3):
        fr=tuple(sum(box_raw(F,r,k,Q(2))/factorial(k) for k in range(m+2-r)) for r in range(m+1))
        ur=tuple(sum(box_raw(U,r,k,Q(2))/factorial(k) for k in range(m+2-r)) for r in range(m+1))
        for n in (0,1):
            bounds=reconstruct(h,n,m,2,1,2,fr,ur)
            for X0,R0 in ((Q(1,2),Q(1)),(Q(9,8),Q(3,2)),(Q(2),Q(2))):
                for eta0 in (Q(-1),Q(0),Q(3,5),Q(1)):
                    X,eta=[J.variable(i,v,n=2,order=9) for i,v in enumerate((X0,eta0))]
                    R=R0*(X/X0)**Q(1,2);L=1-2*h*eta*eta;d=1-eta*eta
                    f=1+X*eta+X*X;u=1+X+X*X*eta*eta+X*eta**3
                    M=X+X*X/2+X**3*eta*eta/3+X*X*eta**3/2
                    V=(2*eta*X*u-2*(D+2*n*h)*eta*M-d*M.derivative(1))/L
                    fields=dict(E=R*f,U=u,W=V/R,J=2*eta*M/(L*R))
                    key=f'm{m}_n{n}_X{X0}_eta{eta0}'
                    for name,field in fields.items():
                        checks[key+'_'+name+'_all_raw_mixed_jets_bounded']=all(abs(field.partial(a))<=upper for a,upper in bounds['raw'][name].items())
                    vx=(2*(A-2*n*h)*eta*u+2*eta*X*u.derivative(0)-d*u.derivative(1))/L
                    checks[key+'_exact_divergence_derivative']=all((V.derivative(0)-vx).partial(a)==0 for a in ((0,0),(1,0),(0,1)))
            rows.append(dict(order=m,coefficient=n,complete_tuple_bound=str(bounds['through'][m]),angular_input=m+1))
    return dict(checks=checks,rows=rows,scope='exact rational jets of manufactured polynomial F,U with the true axis-to-X primitive, including both eta endpoints; not selected PDE profiles')


def memory_review():
    J=jets();h=Q(1,200);D=Q(1,2)-h;mu=Q(3,7)
    X,eta=[J.variable(i,v,n=2,order=6) for i,v in enumerate((Q(1,2),Q(0)))]
    L=1-2*h*eta*eta;d=1-eta*eta;M=mu*eta;U=0*X
    V0=(-2*D*eta*M-d*M.derivative(1))/L
    V1=(-2*(D+2*h)*eta*M-d*M.derivative(1))/L
    v=V0/X
    # Complete leading source on this plateau, with U=0. Its viscosity
    # term is retained before the exact constant-V cancellation.
    source=(V0.derivative(0)+D*eta*v.derivative(1))/L+v*(V0.derivative(0)-v/2)-2*V0.derivative(0).derivative(0)
    checks={
        'zero_local_U_does_not_erase_incoming_radial_memory':V0.value()==-mu!=0,
        'leading_plateau_V_is_exactly_minus_mu':not (V0+mu).terms,
        'positive_order_radial_memory_is_not_leading_plateau_formula':V1.partial((0,2))!=0,
        'leading_plateau_pressure_source_does_not_vanish':source.value()!=0 and all((source+(mu*mu/2)*X**-2).partial(a)==0 for a in ((0,0),(1,0),(0,1),(2,0),(0,2))),
    }
    for N in (4,16,64):
        c1=Q(1,N)+1;raw=N-Q(2)*(D+2*h)/N
        checks[f'N{N}_low_angular_U_norm_bounded_but_radial_first_jet_grows']=c1<=2 and raw>N/2
        # cos(N eta)/N^2 has raw second derivative one at zero,
        # but factorial Ceta2 norm is <=1/2+1/N+1/N^2<1.
        checks[f'N{N}_factorial_row_is_not_identical_raw_second_jet']=Q(1,2)+Q(1,N)+Q(1,N*N)<1
    return dict(checks=checks,scope='exact primitive-plateau and abstract smooth cosine norm controls; earlier/later disjoint radial bumps can close the final moment while preserving incoming memory, not PDE counterexamples')


def join_review():
    f=fixture();s=f['state'];out=f['result'];checks={}
    checks['first_U_value_agrees_with_frozen_K13']=s['jets']['U'][0]==s['repair']['K13']
    checks['same_full_pressure_Q03_Q13_enter_first_repair']=s['repair']['source']['Q03']>0 and s['repair']['source']['Q13']>0 and s['repair']['discrepancies'][2]>0 and s['repair']['discrepancies'][4]>0
    checks['first_mixed_rows_keep_coefficient_cutoff_derivative']=s['jets']['cut'][1]==f['inner'][1]+f['cutoff'][1]*f['inner'][0]
    checks['repair_radial_bump_derivative_is_retained']=s['jets']['U'][1]>s['jets']['cut'][1]
    checks['first_tuple_first_jet_uses_angular_order_two']=out['first']['angular_input_order']==2
    checks['leading_tuple_second_jet_uses_angular_order_three']=out['leading']['angular_input_order']==3
    checks['Cphi_has_all_seven_parent_groups']=set(out['comparison']['groups'])=={'tilt','radial_box','radial_background','radial_rounding','angular_rounding','angular_box','axial_transport'}
    data=dict(f['first_data']);data['nu']=0
    no_visc=finalized_first(f['h'],f['inner'],f['cutoff'],data,Q(1,2),1)
    checks['nonzero_axis_viscosity_boundary_not_dropped']=s['repair']['source']['Q03']>no_visc['repair']['source']['Q03'] and s['repair']['source']['Q13']==no_visc['repair']['source']['Q13']
    data=dict(f['first_data']);data['leading_u5']=0
    no_pressure=finalized_first(f['h'],f['inner'],f['cutoff'],data,Q(1,2),1)
    checks['removing_full_pressure_memory_changes_the_finalized_first_norm']=s['jets']['mixed_Ceta3']>no_pressure['jets']['mixed_Ceta3']
    old=head.length_compatibility(1,2,out['comparison']['Cphi'],Q(1,6057072),100000,23,Q(1,11))
    checks['manufactured_derived_Cphi_does_not_automatically_certify_selector']=not old['gates']['original_length_budget_dominates_normal_margin']
    # Independently expand the first-order closed formulas (8)--(9).
    h=f['h'];ell=1-2*h;c=Q(1,2)-h+2*h;Xhi=Q(8);Rlo=Q(1,2)
    fr=s['jets']['F'];ur=s['jets']['U'];H0=Xhi*ur[0]*(3+2*c)
    v0=H0/ell;vx=(2*Xhi*ur[1]+(3+2*c)*ur[0])/ell
    veta=Xhi*ur[0]*(8+4*c)/ell+4*h*H0/ell**2
    expected=dict(W={(0,0):v0/Rlo,(1,0):vx/Rlo+v0/Rlo**3,(0,1):veta/Rlo},
                  J={(0,0):2*Xhi*ur[0]/(ell*Rlo),
                     (1,0):2*ur[0]/(ell*Rlo)+2*Xhi*ur[0]/(ell*Rlo**3),
                     (0,1):2*Xhi*ur[0]*(2/ell+4*h/ell**2)/Rlo})
    checks['independent_first_order_W_J_closed_form_equals_majorant']=all(out['first']['raw'][name]==rows for name,rows in expected.items())
    return dict(checks=checks,leading_through_two=str(out['leading']['through'][2]),first_through_one=str(out['first']['through'][1]),Cphi=str(out['comparison']['Cphi']),
                scope='manufactured rational assembly through unchanged historical first coefficient, cutoff, exact linear repair and phase recipes; not actual enormous constant evaluation')


@lru_cache(None)
def review():
    parts=dict(reconstruction=exact_review(),memory=memory_review(),join=join_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
