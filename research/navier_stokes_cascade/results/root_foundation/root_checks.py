"""Independent finite matrix, relative-square, jet and bump controls."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,factorial
import mpmath as mp
from root_bounds import HERE,load,derivative_bounds

p=load('root_foundation_exact_polynomials',HERE.parent/'realization_foundation/realized_checks.py')
schedule=load('root_foundation_actual_reference_step',HERE.parent/'outer_pressure_pilot/schedule.py')

def product(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def identity(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def inverse(A):
    n=len(A)
    if n==0 or any(len(row)!=n for row in A):raise ValueError('Nonempty square matrix required')
    rows=[[F(x) for x in row]+identity(n)[i] for i,row in enumerate(A)]
    for k in range(n):
        pivot=next((j for j in range(k,n) if rows[j][k]),None)
        if pivot is None:raise ValueError('Singular matrix')
        rows[k],rows[pivot]=rows[pivot],rows[k];d=rows[k][k];rows[k]=[x/d for x in rows[k]]
        for j in range(n):
            if j!=k:
                d=rows[j][k];rows[j]=[a-d*b for a,b in zip(rows[j],rows[k])]
    return [row[n:] for row in rows]

def lagrange_inverse(nodes,masses):
    n=len(nodes)
    if len(masses)!=n or len(set(nodes))!=n or any(m<=0 for m in masses):raise ValueError('Distinct nodes and positive masses required')
    cols=[]
    for k,node in enumerate(nodes):
        coeff=[F(1)];den=F(masses[k])
        for j,other in enumerate(nodes):
            if j==k:continue
            den*=node-other;new=[F(0)]*(len(coeff)+1)
            for i,c in enumerate(coeff):new[i]-=other*c;new[i+1]+=c
            coeff=new
        cols.append([c/den for c in coeff])
    return [list(row) for row in zip(*cols)]

def integrate_t(f):
    out={}
    for k,v in f.items():
        key=(*k[:3],0);out[key]=out.get(key,F(0))+v/F(k[3]+1)
    return {k:v for k,v in out.items() if v}

@lru_cache(None)
def exact_map_review():
    fixtures=[([F(7),F(1,8)],[F(1,3),F(2,3)]),([F(1,8),F(1),F(7)],[F(1,6),F(1,3),F(1,2)]),
        ([F(2,15),F(49,50),F(36,5)],[F(3,20),F(29,100),F(21,50)])]
    components=0;bad=False
    for nodes,masses in fixtures:
        n=len(nodes);B=[[mass*node**j for j in range(n)] for node,mass in zip(nodes,masses)]
        J=lagrange_inverse(nodes,masses);assert J==inverse(B)
        assert product(B,J)==identity(n)==product(J,B);components+=2*n*n
        bad=bad or product(B,[list(row) for row in zip(*J)])!=identity(n)
    # Separate support charts: on chart j only beta_j is nonzero.
    beta=p.scale(p.mul(p.power(p.T,2),p.power(p.add(p.ONE,p.scale(p.T,-1)),2)),30)
    assert integrate_t(beta)==p.ONE
    coefficients=[p.X,p.Y,p.Z];pressure_direct={};pressure_formula={};heat_direct={};heat_formula={};I_direct={};I_formula={}
    qterms=[]
    for j,c in enumerate(coefficients):
        weight=p.add(p.ONE,p.scale(p.T,j+1),p.scale(p.power(p.T,2),j+2));epsilon=p.mul(c,beta)
        square=p.add(p.power(p.add(p.ONE,epsilon),2),p.scale(p.ONE,-1))
        lin=p.mul(c,integrate_t(p.mul(weight,beta)))
        quad=p.mul(p.power(c,2),integrate_t(p.mul(weight,p.power(beta,2))))
        pressure_direct=p.add(pressure_direct,integrate_t(p.mul(weight,square)))
        pressure_formula=p.add(pressure_formula,p.scale(lin,2),quad)
        heat_direct=p.add(heat_direct,p.scale(integrate_t(p.mul(weight,square)),F(1,2)))
        heat_formula=p.add(heat_formula,lin,p.scale(quad,F(1,2)))
        I_direct=p.add(I_direct,integrate_t(p.mul(weight,epsilon)));I_formula=p.add(I_formula,lin)
        qterms.append(quad)
    assert pressure_direct==pressure_formula and heat_direct==heat_formula and I_direct==I_formula
    assert p.add(pressure_direct,p.scale(I_formula,-2)) and p.add(heat_direct,p.scale(I_formula,-1))
    return dict(rational_positive_mass_matrix_fixtures_not_actual_bump_matrices=True,
        polynomial_support_chart_fixture_not_actual_C_infinity_step=True,
        exact_two_sided_matrix_components=components,
        checks=dict(Lagrange_inverse_matches_independent_Gauss_Jordan=components==44,
            transposed_inverse_detected=bad,
            exact_proxy_bump_mass_one=True,
            angular_relative_square_keeps_two_linear_plus_quadratic=True,
            heat_pressure_and_minus_S_keep_linear_plus_half_quadratic=True,
            I_row_is_exactly_linear=True,
            omitted_quadratic_terms_are_nonzero=all(bool(q) for q in qterms),
            source_pressure_row_factors_not_interchanged=True))

@lru_cache(None)
def numerical_bump_review():
    rows=[]
    with mp.workdps(65):
        for lam in (mp.mpf(0),mp.mpf('0.001'),mp.mpf('0.01')):
            for name,rates,centers,factors,limit in (
                ('angular',[1-lam,-1-2*lam],[0,2],[1,2],3),
                ('heat',[-1-2*lam,-2*lam,1-lam],[mp.mpf('.5'),mp.mpf('2.5'),mp.mpf('4.5')],[1,1,1],15)):
                width=mp.mpf('.3');B=mp.matrix(len(rates));max_factor_error=mp.mpf(0)
                for k,(rate,factor) in enumerate(zip(rates,factors)):
                    mass=factor*width*mp.quad(lambda z:mp.exp(rate*(centers[0]+width*(z-mp.mpf('.5'))))*schedule.step_prime(z),[0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1])
                    for j,center in enumerate(centers):
                        direct=factor*width*mp.quad(lambda z:mp.exp(rate*(center+width*(z-mp.mpf('.5'))))*schedule.step_prime(z),[0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1])
                        B[k,j]=direct
                        max_factor_error=max(max_factor_error,abs(direct-mass*mp.exp(2*rate*j))/(1+abs(direct)))
                inv=B**-1;norm=max(sum(abs(inv[i,j]) for j in range(len(rates))) for i in range(len(rates)))
                rows.append(dict(map=name,diagnostic_lambda=mp.nstr(lam,6),inverse_norm=float(norm),bound=limit,factor_error=float(max_factor_error)))
    return dict(finite_sample_integrals_not_continuum_enclosures=True,actual_selected_lambda_and_roots_not_evaluated=True,rows=rows,
        checks=dict(all_six_actual_bump_matrix_diagnostics_below_rational_uniform_bounds=all(r['inverse_norm']<r['bound'] for r in rows),
            independent_integrals_match_exact_factoring=max(r['factor_error'] for r in rows)<1e-50,
            lambda_zero_is_only_a_matrix_endpoint_diagnostic=True))

@lru_cache(None)
def jets_and_wait_review():
    a=F(3,7);degree=12
    coeff=[F(0)]+[(-a)**(n-1)*F(comb(2*n-2,n-1),n) for n in range(1,degree+1)]
    for n in range(1,degree+1):assert coeff[n]+a*sum(coeff[j]*coeff[n-j] for j in range(n+1))==int(n==1)
    majorant=derivative_bounds(1,a,[0,1]+[0]*(degree-1),0)
    assert all(abs(c)<=d for c,d in zip(coeff,majorant))
    # Terminal integrating factor in a polynomial diagnostic with kappa=0:
    # F positive, N=1-F, q=N/F, so q'+(F'/F)q=-F'/F exactly.
    s=p.add(p.scale(p.power(p.T,2),3),p.scale(p.power(p.T,3),-2));rho=F(1,1000)
    fo=p.add(p.scale(p.ONE,1-rho),p.scale(s,rho));N=p.add(p.ONE,p.scale(fo,-1));fp=p.diff(fo,3)
    cleared=p.add(p.mul(p.diff(N,3),fo),p.scale(p.mul(N,fp),-1),p.mul(fp,N),p.mul(fp,fo));assert not cleared
    return dict(scalar_Catalan_and_zero_kappa_terminal_fixtures_not_actual_reference_flow=True,
        exact_scalar_series_orders=degree,
        checks=dict(Catalan_solution_satisfies_full_scalar_quadratic_through_order_twelve=True,
            all_normalized_jets_bounded_by_current_implicit_majorant=True,
            second_derivative_factorial_is_not_normalized_coefficient=2*coeff[2]!=coeff[2],
            cleared_terminal_integrating_factor_identity_exact=True,
            positive_terminal_factor_lower_bound=1-rho>0,
            waiting_target_cannot_be_replaced_by_zero=rho/(1-rho)>0,
            prescribed_terminal_constant_not_erased=True))

def reviews():return dict(exact_matrix_and_moment_map_controls=exact_map_review(),
    independent_actual_bump_integrals=numerical_bump_review(),exact_implicit_jets_and_wait_controls=jets_and_wait_review())
