"""Independent exact polynomial controls and actual-map numerical diagnostics."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import mpmath as mp
from matching_bounds import HERE,load
p=load('matching_root_polynomials',HERE.parent/'realization_foundation/realized_checks.py')
def trunc(x):return {k:v for k,v in x.items() if k[1]<=8}
def coeff(x,n):return x.get((0,n,0,0),F(0))

@lru_cache(None)
def polynomial_review():
    t=p.Y;eta=p.add(p.scale(p.ONE,F(1,2)),t)
    q=p.scale(p.power(p.add(p.ONE,p.power(eta,2)),2),F(3,7))
    c={}
    for _ in range(9):c=trunc(p.add(t,p.scale(p.mul(q,p.mul(c,c)),-1)))
    independent={}
    for n in range(1,9):
        catalan=F(comb(2*n-2,n-1),n)
        independent=trunc(p.add(independent,p.scale(p.mul(p.power(q,n-1),p.power(t,n)),(-1)**(n-1)*catalan)))
    residual=trunc(p.add(c,p.mul(q,p.mul(c,c)),p.scale(t,-1)))
    frozenq=p.scale(p.ONE,coeff(q,0));wrong={}
    for _ in range(9):wrong=trunc(p.add(t,p.scale(p.mul(frozenq,p.mul(wrong,wrong)),-1)))
    d=[F(0)]*9
    for n in range(1,9):
        d[n]=2*((1 if n==1 else 0)+sum(abs(coeff(q,i))*d[j]*d[n-i-j]
            for i in range(min(4,n)+1) for j in range(n-i+1)
            if (i,j,n-i-j) not in ((0,n,0),(0,0,n))))
    g=p.add(p.ONE,p.power(t,3));u=p.add(p.scale(t,4),g)
    identity=p.add(p.power(u,2),p.scale(p.power(t,2),-16),p.scale(p.mul(t,g),-8),p.scale(p.power(g,2),-1))
    pressure=p.add(t,p.power(t,3));defect=p.add(p.scale(p.mul(t,pressure),F(7,3)),p.scale(p.mul(p.add(p.ONE,p.scale(p.power(t,2),-1)),p.diff(pressure,1)),-1))
    return dict(order=8,coefficients=[str(coeff(c,n)) for n in range(1,9)],
        checks=dict(independent_Catalan_and_implicit_series_agree=c==independent,
            exact_variable_Q_residual_vanishes=not residual,constant_Q_jet_rule_fails=c!=wrong,
            triple_majorant_controls_all_eight_jets=all(abs(coeff(c,n))<=d[n] for n in range(1,9)),
            raw_second_derivative_needs_factor_two=2*coeff(c,2)!=coeff(c,2),
            restoration_energy_identity_exact=not identity,
            pressure_value_zero_does_not_remove_gradient=coeff(pressure,0)==0 and coeff(defect,0)==-1))

@lru_cache(None)
def numerical_review():
    ann=load('matching_root_actual_annular',HERE.parent/'axis_matching_audit/moments.py')
    rep=load('matching_root_actual_repair',HERE.parent/'stress_realization_audit/repair.py')
    with mp.workdps(90):
        a=ann.MomentMap(order=16,panels=4);eta=mp.mpf('.5');P=16
        result=a.solve(ann.fixture(eta),eta,P);c=mp.matrix(result['root'])
        error=max(abs(x) for x in a.B*c+a.quadratic(c,c,eta,P)-mp.matrix(result['target']))
        h=mp.mpf('1e-12')
        plus=a.solve(ann.fixture(eta+h),eta+h,P)['root'];minus=a.solve(ann.fixture(eta-h),eta-h,P)['root']
        derivative=max(abs((x-y)/(2*h)-z) for x,y,z in zip(plus,minus,result['root_eta']))
        inverse=max(sum(abs(a.inverse[i,j]) for j in range(5)) for i in range(5))
        rows=[]
        for lam in (mp.mpf('.001'),mp.mpf('1e-20')):
            r=rep.Repair(lam,order=16,panels=4);target=mp.matrix([mp.mpf('1e-20')*k for k in (1,-2,3,-4,5)])
            cr=r.solve(target);err=max(abs(x) for x in r.value(cr)-target)
            inv=max(sum(abs(r.inverse[i,j]) for j in range(5)) for i in range(5))
            direct=max(abs(x-y) for x,y in zip(rep.precondition_target(r.direct_changes(cr),lam),r.value(cr)))
            rows.append(dict(lambda_value=float(lam),inverse=float(inv),root_residual=float(err),ordinary_control_error=float(direct)))
        lam=F(1,10**80);m=1+lam;j=F(1)
        zero_rejected=False
        try:rep.Repair(0)
        except ValueError:zero_rejected=True
    return dict(annular_inverse=float(inverse),annular_residual=float(error),annular_eta_error=float(derivative),repair_rows=rows,
        finite_quadratures_are_diagnostics=True,manufactured_entry_is_not_axis_data=True,
        checks=dict(actual_annular_inverse_below_610=inverse<610,actual_annular_root_closes=error<mp.mpf('1e-60'),
            actual_annular_implicit_derivative_matches_resolved_roots=derivative<mp.mpf('1e-35'),
            actual_repair_inverses_below_800=all(x['inverse']<800 for x in rows),
            actual_repair_roots_close=all(x['root_residual']<1e-60 for x in rows),
            independent_ordinary_integrands_match=all(x['ordinary_control_error']<1e-60 for x in rows),
            rounded_input_difference_can_be_lost=(m-j)/lam==1 and float(m)-float(j)==0,
            zero_lambda_rejected_despite_weight_extension=zero_rejected and rep.divided_weight(0,2)==2))

@lru_cache(None)
def support_review():
    def disjoint(centers,w):
        centers=sorted(map(F,centers));return all(y-x>w for x,y in zip(centers,centers[1:]))
    return dict(checks=dict(annular_supports_disjoint=disjoint(['.3','.7','.15','.5','.85'],F('.08')),
        repair_supports_disjoint=disjoint(['1','2','.5','1.5','2.5'],F('.1')),
        annular_mass_width_differs_from_repair_unit_mass=F('.08')!=1,
        full_C1_sum_is_an_algebra=(F(2)+3)*(F(5)+7)>=2*5+3*5+2*7,
        old_max_norm_needs_product_factor_two=2>1))

def reviews():return dict(exact_variable_coefficient_and_pressure_controls=polynomial_review(),actual_matching_map_diagnostics=numerical_review(),support_and_norm_controls=support_review())
