"""Finite-order radial/physical jet budgets and a bounded-C2 cutoff control.

Only the retained repair bumps have explicit all-order bounds here. Shape
jets of the unchanged coefficient cutoff remain inputs, not guessed values.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import importlib.util
import sys

HERE=Path(__file__).resolve().parent
def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational input required')
    return F(x)
def order(m):
    if type(m) is not int or m<0:raise ValueError('Nonnegative integer order required')
    return m
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def module(name,file):
    p=HERE.parent/name/file;s=importlib.util.spec_from_file_location('higher_cutoff_'+name,p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m

def log_coefficients(m,alpha=0):
    """Coefficients of product_(i=0)^(m-1)(D_y-alpha-i)."""
    m=order(m);alpha=rat(alpha);v=[F(1)]
    for i in range(m):
        w=[F(0)]*(len(v)+1)
        for j,c in enumerate(v):w[j]-=(alpha+i)*c;w[j+1]+=c
        v=w
    return tuple(v)
def norm_list(v):
    v=tuple(map(rat,v))
    if not v or min(v)<0:raise ValueError('Nonempty nonnegative norm list required')
    return v
def radial_cutoff(jets,delta,xkeep):
    """L_j enclose shape jets through m; return radial k_0,...,k_m."""
    L=norm_list(jets);delta,xkeep=map(rat,(delta,xkeep))
    if min(delta,xkeep)<=0 or L[0]!=1:raise ValueError('Positive log width/radius and L0=1 required')
    return tuple(sum(abs(c)*L[j]/delta**j for j,c in enumerate(log_coefficients(m)))/xkeep**m for m in range(len(L)))
def product_jets(k,J):
    k,J=norm_list(k),norm_list(J)
    if len(k)!=len(J):raise ValueError('Every product jet order needs both inputs')
    return tuple(sum(F(comb(m,j))*k[j]*J[m-j] for j in range(m+1)) for m in range(len(k)))
def bump_radial(B,p,rp,kind):
    """U bump beta(y), or F bump beta(y)/sqrt(2X); p lower X."""
    B=norm_list(B);p,rp=map(rat,(p,rp))
    if kind not in ('U','F') or min(p,rp)<=0 or rp*rp>2*p:raise ValueError('Valid bump type and positive radius enclosures required')
    alpha=F(1,2) if kind=='F' else F(0)
    return tuple(sum(abs(c)*B[j] for j,c in enumerate(log_coefficients(m,alpha)))/p**m/(rp if kind=='F' else 1) for m in range(len(B)))
def extended_jets(k,J,B,coefficient,p,rp):
    """Common C_eta^a norms; disjoint additive repair supports, not fivefold."""
    coefficient=rat(coefficient)
    if coefficient<0:raise ValueError('Nonnegative full repair coefficient norm required')
    v=product_jets(k,J);u=bump_radial(B,p,rp,'U');f=bump_radial(B,p,rp,'F')
    if len(v)!=len(u):raise ValueError('Bump and collar jets required through the same order')
    return tuple(v[m]+coefficient*max(u[m],f[m]) for m in range(len(v)))
def physical_coefficients(m):
    """d_r^m f(r^2/(2q)): positive coefficient of r^(2k-m) q^-k f^(k)."""
    m=order(m)
    return {k:factorial(m)//(2**(m-k)*factorial(m-k)*factorial(2*k-m)) for k in range((m+1)//2,m+1)}
def physical_bound(S,rhi,m):
    """Bound after extracting a fixed prefactor times q^(-m/2)."""
    S=norm_list(S);rhi=rat(rhi);m=order(m)
    if len(S)<=m or rhi<=0:raise ValueError('All profile jets and positive normalized upper radius required')
    return sum(F(c)*rhi**(2*k-m)*S[k] for k,c in physical_coefficients(m).items())

def partitions(m,lower=1):
    """Integer partitions as sorted derivative orders; finite Bell expansion."""
    if m==0:yield ()
    for j in range(lower,m+1):
        for tail in partitions(m-j,j):yield (j,)+tail
@lru_cache(None)
def logistic_polynomial(k):
    """P_k(z) for d_t^k logistic(t), exact integer coefficients."""
    k=order(k);v=[0,1]
    for _ in range(k):
        w=[0]*(len(v)+1)
        for j in range(1,len(v)):w[j]+=j*v[j];w[j+1]-=j*v[j]
        v=w
    return tuple(v)
@lru_cache(None)
def generic_sigma_bound(m):
    """Analytical bound for retained squared-flat sigma^(m), no sampling."""
    m=order(m)
    if m==0:return 1
    total=F(0)
    for parts in partitions(m):
        k=len(parts);counts={j:parts.count(j) for j in set(parts)}
        coeff=F(factorial(m))
        for j,n in counts.items():coeff*=F((2*(j+1))**n,factorial(n))
        # z^(m/2+k) <= z^ceil(m/2+k), z>=4.
        p=(m+1)//2+k;power=max(4,p)**p
        D=sum(abs(c) for c in logistic_polynomial(k)[1:])
        total+=coeff*D*power
    if total.denominator!=1:raise ArithmeticError('Bell coefficient sum should be integral')
    return total.numerator
def sigma_bound(m):
    m=order(m)
    known={0:1,1:32,2:1216,3:220000}
    return min(known.get(m,generic_sigma_bound(m)),generic_sigma_bound(m))
def repair_bump_bounds(m):
    """Actual width 1/10, beta=10 sigma'(10(y-c)+1/2)."""
    m=order(m)
    return tuple(10**(j+1)*sigma_bound(j+1) for j in range(m+1))

def control_constants():
    """A concrete auxiliary monotone family, NOT the unchanged cutoff."""
    margin=F(1,2*3**16);psi1=F(256);psi2=F(77824)
    a=margin/(4*(1+psi1+psi2))
    return dict(value_margin=margin,base_slope_lower=16*margin,psi1=psi1,psi2=psi2,a=a)
def control(n):
    if type(n) is not int or n<1:raise ValueError('Positive integer oscillation frequency required')
    c=control_constants();a,p1,p2=(c[k] for k in ('a','psi1','psi2'))
    eps=a/F(n*n);d1=a*(1/F(n)+p1/F(n*n));d2=a*(1+2*p1/F(n)+p2/F(n*n))
    return dict(n=n,value_difference=eps,first_difference=d1,second_difference=d2,
        M=32+d1,N=1216+d2,center_value_difference=F(0),center_first_difference=a/n,
        center_second_difference=F(0),center_third_difference=-a*n,
        center_third_of_base=F(3712),center_third_total=F(3712)-a*n,
        checks=dict(values_in_unit_interval=eps<c['value_margin'],
            strictly_decreasing_on_perturbation_support=d1<c['base_slope_lower'],
            uniform_first_enclosure=32+d1<33,uniform_second_enclosure=1216+d2<1217))

def symbolic_review():
    import sympy as s
    X,y,r,q,delta,nu=s.symbols('X y r q delta nu',positive=True)
    z=s.symbols('z');p=y**9-2*y**6+3*y+1;checks={};rows=[]
    for m in range(9):
        for alpha in (F(0),F(1,2)):
            direct=s.diff(X**(-alpha)*p.subs(y,s.log(X)/delta),X,m)
            expected=X**(-alpha-m)*sum(c*delta**(-j)*s.diff(p,y,j).subs(y,s.log(X)/delta) for j,c in enumerate(log_coefficients(m,alpha)))
            ok=s.simplify(direct-expected)==0;rows.append(dict(order=m,alpha=str(alpha),matches=ok))
        # Independent Bell chain rule on polynomial in X.
        poly=X**9-X**4+2*X
        got=sum(c*r**(2*k-m)*q**(-k)*s.diff(poly,X,k).subs(X,r*r/(2*q)) for k,c in physical_coefficients(m).items())
        checks['physical_chain_order_'+str(m)]=s.expand(s.diff(poly.subs(X,r*r/(2*q)),r,m)-got)==0
    f=s.Function('f')(X);u=s.Function('u')(X)
    Tt=2*nu*X*s.diff(f,X);Tz=2*nu*X*s.diff(u,X)/s.sqrt(2*X)
    checks.update(theta_viscous_primitive=s.simplify(s.diff(X*Tt,X)-2*nu*X*(X*s.diff(f,X,2)+2*s.diff(f,X)))==0,
        axial_viscous_primitive=s.simplify(s.diff(s.sqrt(2*X)*Tz,X)-2*nu*(X*s.diff(u,X,2)+s.diff(u,X)))==0,
        theta_second_stress_principal=s.expand(s.diff(Tt,X,2)).coeff(s.diff(f,X,3))==2*nu*X,
        axial_second_stress_principal=s.simplify(s.expand(s.diff(Tz,X,2)).coeff(s.diff(u,X,3))-nu*s.sqrt(2*X))==0,
        third_log_jacobian_retained=log_coefficients(3)==(0,2,-3,1))
    a,n,t=s.symbols('a n t',positive=True);osc=a/n**2*s.sin(n*t)
    checks['center_third_oscillation']=s.diff(osc,t,3).subs(t,0)==-a*n
    checks['center_second_oscillation_zero']=s.diff(osc,t,2).subs(t,0)==0
    for k in range(1,7):
        poly=sum(c*z**j for j,c in enumerate(logistic_polynomial(k)))
        prev=sum(c*z**j for j,c in enumerate(logistic_polynomial(k-1)))
        checks['logistic_recurrence_'+str(k)]=s.expand(poly-z*(1-z)*s.diff(prev,z))==0
    return dict(log_coordinate_comparisons=rows,checks=checks)

def review():
    symbolic=symbolic_review();B=repair_bump_bounds(5)
    L=tuple(map(F,(1,33,1217,10**9,10**12,10**15)))
    k=radial_cutoff(L,F(1,64),2);J=tuple(map(F,(2,3,5,7,11,13)))
    full=extended_jets(k,J,B,17,9,3)
    old=module('first_seminorm_bound','seminorms.py')
    oldjets=old.extended_jets(dict(value=J[0],radial=J[1],radial2=J[2]),17,1,k[1],k[2],9,3,B[:3])
    # Equality for orders 0..2 because L0=1 and the same explicit bump budgets.
    low=tuple(oldjets[n] for n in ('value','radial','radial2'))
    rows=[control(n) for n in (1,2,1000,10**12,10**24,10**36)]
    # Full w0 memory is present in each absolute ledger. It cancels only in
    # DIFFERENCES at order one, where the leading source is held fixed.
    mom0=old.moment_bounds(2,0,5,7,(11,13))
    pressure=[]
    for row in rows:
        eps=row['value_difference'];fullmom=old.moment_bounds(2,2*eps,5,7,(11,13))
        dm=tuple(x-y for x,y in zip(fullmom,mom0));rep=old.repair_bound(dm,2,3,F(1,1000))
        pressure.append(dict(n=row['n'],full_absolute_moment_bounds=serial(fullmom),
            fixed_pressure_baseline=serial(mom0),difference_bounds=serial(dm),repair_coefficient_bound=str(rep['coefficient']),
            coefficient_bound_over_cutoff_value_error=str(rep['coefficient']/eps)))
    floor=module('sharp_current_scale_floor','cost_floor.py')
    accepted=[]
    for row in rows:
        vals={name:floor.shape().affine_value(v,row['M'],row['N']) for name,v in floor.toy(A1=0).items()}
        accepted.append(dict(n=row['n'],cost=str(max(vals.values())),all_rows_below_common_budget=max(vals.values())<1218))
    checks=dict(all_log_coordinate_identities=all(v['matches'] for v in symbolic['log_coordinate_comparisons']),
        all_symbolic_chain_and_principal_identities=all(symbolic['checks'].values()),
        first_three_actual_bump_bounds_unchanged=B[:3]==(320,121600,220000000),
        all_six_actual_bump_budgets_positive=all(x>0 for x in B),
        low_order_interface_agrees_with_original_first_seminorms=full[:3]==low,
        all_auxiliary_shapes_keep_fixed_low_jet_enclosures=all(all(v['checks'].values()) for v in rows),
        common_affine_budget_accepts_all_auxiliary_controls=all(v['all_rows_below_common_budget'] for v in accepted),
        third_jet_grows_while_first_two_remain_bounded=abs(rows[-1]['center_third_total'])>10**15,
        full_pressure_present_in_each_absolute_ledger=all(v['fixed_pressure_baseline'][2]=='11/2' and v['fixed_pressure_baseline'][4]=='13/2' for v in pressure),
        repair_difference_is_uniformly_value_small=len({v['coefficient_bound_over_cutoff_value_error'] for v in pressure})==1)
    return dict(symbolic=symbolic,actual_repair_bump_log_bounds=list(map(str,B)),
        conditional_finite_jet_example=serial(dict(shape=L,radial_cutoff=k,collar=J,repaired=full,
            physical_normalized_order5_bound=physical_bound(full,8,5))),
        auxiliary_control_constants=serial(control_constants()),auxiliary_controls=serial(rows),
        low_order_affine_controls=accepted,complete_pressure_difference_controls=pressure,
        checks=checks)
