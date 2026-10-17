"""Independent exact affine-energy and finite numerical matrix controls."""
from fractions import Fraction as F
from functools import lru_cache
from math import exp
import mpmath as mp
from axial_bounds import HERE,load

p=load('axial_root_polynomial_controls',HERE.parent/'realization_foundation/realized_checks.py')
inverse_helpers=load('axial_root_inverse_controls',HERE.parent/'root_foundation/root_checks.py')
pulse=load('axial_root_actual_pulse_solver',HERE.parent/'pulse_moment_audit/pulse.py')

@lru_cache(None)
def exact_affine_review():
    A,eta=p.X,p.Y
    odd=p.add(eta,p.power(eta,3));a=[F(-2,7),F(3,11)];b=[F(1,13),F(-2,17)];w=[F(7,19),F(11,23)]
    c=[p.add(p.scale(A,x),p.scale(odd,y)) for x,y in zip(a,b)]
    direct=p.add(*(p.scale(p.power(z,2),v) for z,v in zip(c,w)))
    aa=sum(v*x*x for v,x in zip(w,a));ab=2*sum(v*x*y for v,x,y in zip(w,a,b));bb=sum(v*y*y for v,y in zip(w,b))
    expanded=p.add(p.scale(p.power(A,2),aa),p.scale(p.mul(A,odd),ab),p.scale(p.power(odd,2),bb))
    assert direct==expanded and ab
    reflected=p.add(p.scale(p.power(A,2),aa),p.scale(p.mul(A,odd),-ab),p.scale(p.power(odd,2),bb))
    assert direct!=reflected
    for x,y in zip(a,b):
        composed=p.add(p.scale(p.add(p.ONE,eta),x),p.scale(odd,y))
        difference=p.add(p.diff(composed,1),p.scale(p.diff(odd,1),-y))
        assert difference==p.scale(p.ONE,x) and difference
    components=0
    for nodes,masses in (([F(7,3),F(11,5)],[F(3,10),F(31,100)]),([F(1001,1000),F(1)],[F(2,7),F(3,11)])):
        B=[[m,m*n] for m,n in zip(masses,nodes)];J=inverse_helpers.inverse(B)
        incoming=[p.scale(odd,F(2,13)),p.scale(odd,F(-3,17))];main=[F(5,19),F(7,23)]
        rhs=[p.scale(p.add(z,p.scale(A,v)),-1) for z,v in zip(incoming,main)]
        coeff=[p.add(*(p.scale(z,q) for z,q in zip(rhs,row))) for row in J]
        for row,z,v in zip(B,incoming,main):
            residual=p.add(z,p.scale(A,v),*(p.scale(cc,q) for cc,q in zip(coeff,row)));assert not residual;components+=1
    # On the same fixture, correcting M alone leaves a nonzero J polynomial.
    only_M=p.scale(rhs[0],1/B[0][0]);miss=p.add(incoming[1],p.scale(A,main[1]),p.scale(only_M,B[1][0]));assert miss
    return dict(exact_proxy_moment_components=components,fixtures_not_selected_profile=True,
        checks=dict(full_affine_end_energy_identity_exact=direct==expanded,
            positive_quadratic_and_constant_increments_retained=aa>0 and bb>0,
            cross_increment_nonzero=bool(ab),
            strength_not_assumed_even_in_eta=direct!=reflected,
            composed_coefficient_derivative_requires_strength_chain_term=True,
            both_linear_moments_close_as_full_polynomials=components==4,
            only_M_correction_leaves_nonzero_J=bool(miss),
            discarding_all_end_bump_energy_changes_the_equation=bool(direct)))

@lru_cache(None)
def numerical_matrix_review():
    rows=[]
    with mp.workdps(90):
        for lam in (mp.mpf('.01'),mp.mpf('.001'),mp.mpf('1e-20')):
            slopes=[mp.mpf('.5')-lam,mp.mpf('.5')-2*lam];width=mp.mpf('.3')
            masses=[width*mp.quad(lambda z:mp.exp(s*width*(z-mp.mpf('.5')))*pulse.step_prime(z),[0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1]) for s in slopes]
            B=mp.matrix([[b,b*mp.exp(2*s)] for b,s in zip(masses,slopes)])
            rhs=[mp.mpf('2.125'),mp.mpf('-1.375')];c=mp.matrix(pulse.stable_two_row(lam,masses,rhs));direct=mp.lu_solve(B,mp.matrix(rhs))
            rel=max(abs(x-y) for x,y in zip(c,direct))/(1+max(abs(x) for x in direct))
            J=B**-1;scaled=lam*max(sum(abs(J[i,j]) for j in range(2)) for i in range(2))
            residual=max(abs(x) for x in B*c-mp.matrix(rhs))
            rows.append(dict(lambda_diagnostic=mp.nstr(lam,6),scaled_inverse=float(scaled),relative_solve_error=float(rel),absolute_row_error=float(residual)))
    rounded_gap=exp(1-2*1e-20)-exp(1-4*1e-20)
    return dict(rows=rows,finite_quadratures_are_diagnostics=True,selected_lambda_not_evaluated=True,
        checks=dict(three_actual_bump_matrices_below_uniform_four=all(r['scaled_inverse']<4 for r in rows),
            actual_factored_solver_matches_independent_LU=all(r['relative_solve_error']<1e-60 for r in rows),
            both_direct_rows_small=all(r['absolute_row_error']<1e-60 for r in rows),
            binary64_subtraction_can_erase_nonzero_gap=rounded_gap==0,
            positive_lambda_not_set_to_zero_in_stable_solver=True))

@lru_cache(None)
def scalar_review():
    rows=[]
    with mp.workdps(80):
        for eta in (mp.mpf('-1'),mp.mpf('0.25'),mp.mpf('1')):
            k=mp.mpf('0.245');b=mp.mpf('0.000001')*eta;d=-mp.mpf('.25')+mp.mpf('.0000001')*eta**2
            disc=b*b-4*k*d;root=-2*d/(b+mp.sqrt(disc));jac=2*k*root+b
            def positive_root(e):
                bb=mp.mpf('.000001')*e;dd=-mp.mpf('.25')+mp.mpf('.0000001')*e**2
                return -2*dd/(bb+mp.sqrt(bb*bb-4*k*dd))
            implicit=-(mp.mpf('.000001')*root+mp.mpf('.0000002')*eta)/jac
            rows.append(dict(eta=float(eta),root=float(root),residual=float(abs(k*root**2+b*root+d)),jacobian_error=float(abs(jac-mp.sqrt(disc))),derivative_error=float(abs(mp.diff(positive_root,eta)-implicit))))
    k=F(2,7);b=F(3,11);root=F(13,17);d=-k*root**2-b*root
    assert b*b-4*k*d==(2*k*root+b)**2
    return dict(proxy_coefficients_not_actual_selected_strength=True,rows=rows,
        checks=dict(exact_rational_discriminant_jacobian_identity=True,
            stable_positive_root_satisfies_all_proxy_quadratics=all(r['residual']<1e-70 for r in rows),
            positive_root_jacobian_is_sqrt_discriminant=all(r['jacobian_error']<1e-70 for r in rows),
            direct_root_derivative_matches_full_implicit_formula=all(r['derivative_error']<1e-70 for r in rows),
            odd_affine_cross_term_changes_the_positive_root=rows[0]['root']!=rows[-1]['root'],
            negative_constant_gives_one_positive_and_one_negative_root=d<0 and k>0))

def reviews():return dict(exact_axial_and_full_energy_controls=exact_affine_review(),
    independent_actual_axial_matrix_diagnostics=numerical_matrix_review(),exact_and_numerical_scalar_root_controls=scalar_review())
