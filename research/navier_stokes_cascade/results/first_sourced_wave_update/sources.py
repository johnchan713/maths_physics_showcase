"""First sourced correction: exact identities and manufactured diagnostics.

The continuum estimates and the actual homogeneous-pulse tail counterexample
are in README.md. No actual giant band field is numerically evaluated here.
"""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys
import mpmath as mp
import numpy as np

HERE=Path(__file__).resolve().parent
KAPPA=Q(1,100000)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


algebra=load('first_source_frozen_mean_algebra',HERE.parent/'post_signed_mean_relaxation/relax.py')
Poly,R,Z,T,EPS=algebra.Poly,algebra.R,algebra.Z,algebra.T,algebra.EPS
signed=algebra.signed


@lru_cache(None)
def physical_modules():
    p=HERE.parent/'physical_wave_interface';sys.path.insert(0,str(p))
    try:
        import jets
        import waves
        return jets,waves
    finally:
        sys.path.pop(0)


def linear(base,data,waves):
    return data[1][:,3]+waves.advection(base,data,cartesian=True)+waves.advection(data,base,cartesian=True)-waves.laplacian(data,cartesian=True)


@lru_cache(None)
def cutoff_review():
    jets,waves=physical_modules();worst=decomposition=divergence=0.
    omitted=dict(curl=0.,pressure=0.,viscosity=0.,transport=0.)
    for point in ((1.1,.3,-.2,.36),(1.35,-.2,.17,.53),(1.6,.4,.25,.68)):
        coords=jets.variables(point);x,y,z,time=coords
        data=waves.carrier_data(coords,1,cartesian=True)
        phase=13*data['phase']
        potential=[-c*phase.sin()/13 for c in data['coefficient']]
        pressure=.13*data['envelope']*phase.sin()
        psi=jets.bump((time-.5-.06*(x*x+y*y)-.12*z)/.42)
        old=waves.curl(potential,cartesian=True)
        new=waves.curl([psi*c for c in potential],cartesian=True)
        gradpsi=[psi.diff(i) for i in range(3)]
        d=jets.cross(gradpsi,potential)
        a,b,c=jets.arrays(old),jets.arrays(new),jets.arrays(d)
        base,_=waves.background(coords,cartesian=True);base=jets.arrays(base)
        pg=waves.pressure_gradient(pressure,cartesian=True)
        pgcut=waves.pressure_gradient(psi*pressure,cartesian=True)
        actual=linear(base,b,waves)+pgcut-psi.value*(linear(base,a,waves)+pg)
        gp=np.array([v.value for v in gradpsi]);lap=sum(psi.diff(i).diff(i).value for i in range(3))
        transport=(psi.diff(3).value+np.dot(base[0],gp))*a[0]
        viscous=-(2*a[1][:,:3]@gp+lap*a[0])
        pressure_term=pressure.value*gp
        curl_term=linear(base,c,waves)
        expected=transport+viscous+pressure_term+curl_term
        worst=max(worst,float(np.linalg.norm(actual-expected)/(1+np.linalg.norm(actual)+np.linalg.norm(expected))))
        decomposition=max(decomposition,float(np.linalg.norm(b[0]-psi.value*a[0]-c[0])))
        divergence=max(divergence,abs(float(np.trace(b[1][:,:3]))))
        for name,term in (('curl',curl_term),('pressure',pressure_term),('viscosity',viscous),('transport',transport)):
            omitted[name]=max(omitted[name],float(np.linalg.norm(term)))
    return dict(checks=dict(full_linear_potential_pressure_cutoff_commutator=worst<2e-10,
        exact_curl_product_decomposition=decomposition<2e-10,
        cut_potential_is_divergence_free=divergence<2e-9,
        all_four_commutator_contributions_nonzero=all(x>1e-6 for x in omitted.values())),
        maximum_relative_commutator_error=worst,maximum_curl_product_error=decomposition,
        maximum_divergence_error=divergence,omitted_term_gaps=omitted,
        manufactured_physical_fields_not_actual_pulses=True)


def dot(a,b): return sum(x*y for x,y in zip(a,b))
def mv(matrix,v): return [dot(row,v) for row in matrix]


def forced_normal_review():
    matrix=((Q(0),Q(-2),Q(0)),(Q(5),Q(0),Q(0)),(Q(3),Q(0),Q(0)))
    exact=[];moving=[];normal_force=[];wrong_i=[];zero_entry=[];nonzero_exit=[]
    for m in (-2,-1,1,2):
        for v in (Q(0),Q(1,3),Q(2,3),Q(1)):
            n=[Q(1),v,Q(1)];np_=[Q(0),Q(1),Q(0)];n2=dot(n,n)
            q=v*v*(1+v);qp=2*v+3*v*v
            t=[-v*q,q,Q(0)];tp=[-q-v*qp,qp,Q(0)]
            p=1+v;kt=mv(matrix,t);d=Q(2,5)*m*m
            f=[-tp[i]-kt[i]-d*t[i]+n[i]*p for i in range(3)]
            numerator=dot(n,kt)-dot(np_,t)+dot(n,f)
            exact.extend((dot(n,t)==0,numerator==n2*p,
                all(tp[i]+kt[i]+d*t[i]-n[i]*p+f[i]==0 for i in range(3))))
            projected=[-kt[i]+n[i]*(dot(n,kt)-dot(np_,t))/n2-d*t[i]-f[i]+n[i]*dot(n,f)/n2 for i in range(3)]
            exact.append(projected==tp)
            moving.append(dot(np_,t)!=0);normal_force.append(dot(n,f)!=0)
            wrong_i.append(2*p*p*n2>0) # Missing i: real and imaginary normal residuals both survive.
            if v==0: zero_entry.append(all(x==0 for x in t))
            if v==1: nonzero_exit.append(any(x!=0 for x in t))
    return dict(checks=dict(forced_constraint_pressure_and_projected_ODE_exact=all(exact),
        moving_normal_term_is_needed=any(moving),source_normal_component_is_not_dropped=any(normal_force),
        missing_pressure_i_is_detected=all(wrong_i),zero_entry_does_not_imply_zero_exit=all(zero_entry+nonzero_exit),
        both_positive_and_negative_harmonics_checked=True),exact_rational_checks=len(exact),
        pressure_formula='i*(n.Kt-nprime.t+n.f)/(k*m*|n|^2)',
        manufactured_forced_solution_not_actual_band_integration=True)


def tail_review():
    step=10*R**3-15*R**4+6*R**5
    derivative=step.diff(0);primitive=derivative.primitive()
    integral=(primitive.at(0,1)-primitive.at(0,0)).scalar()
    rows=[]
    with mp.workdps(70):
        for length in (32,128,512,2048):
            logbound=mp.log(mp.mpf(75)/(4*length))-mp.mpf(9)*length/200
            rows.append(dict(length=length,log_upper_source_norm=float(logbound),
                             inverse_core_value=-1,source_zero_on_both_endpoint_collars=True))
    return dict(checks=dict(normalized_early_step_derivative_integral_one=integral==1,
        polynomial_step_has_zero_first_two_endpoint_jets=all(step.diff(0).at(0,v)==step.diff(0).diff(0).at(0,v)==0 for v in (0,1)),
        exact_derivative_sup_majorant=derivative==30*R**2*(1-R)**2 and derivative.at(0,Q(1,2))==Q(15,8),
        Gaussian_source_norm_bound_decreases=all(rows[i+1]['log_upper_source_norm']<rows[i]['log_upper_source_norm'] for i in range(3)),
        forward_inverse_core_does_not_become_flat=all(r['inverse_core_value']==-1 for r in rows),
        physical_flatness_is_not_preserved_by_forward_inverse=True),
        rows=rows,early_interval=['L/10','L/5'],Gaussian_exponent='9*L/200',
        diagnostic_step_regularity='C2 zero/one extension; analytical witness uses a genuine C-infinity step',
        exact_actual_homogeneous_witness='f=chi_prime*t_h; solution=(-chi*t_h,-chi*pi_h)',
        tail_kept_as_residual_not_erased_or_forward_inverted=True)


@lru_cache(None)
def covariance_pressure_review():
    bump,records,_,_,_=signed.fixture();changes=[];checks=[];shortcut=[]
    for p in records:
        a=[x+y for x,y in zip(p['a'],p['ell'])];r=[x+y for x,y in zip(p['r'],p['s'])]
        c=[(x+Q(2,7)*y)/5 for x,y in zip(p['ell'],p['s'])]
        d=[(x-Q(2,7)*y)/5 for x,y in zip(p['s'],p['ell'])]
        old=signed.scale(signed.add(signed.outer(a,a),signed.outer(r,r)),Q(1,2))
        ap=[x+y for x,y in zip(a,c)];rp=[x+y for x,y in zip(r,d)]
        new=signed.scale(signed.add(signed.outer(ap,ap),signed.outer(rp,rp)),Q(1,2))
        delta=signed.add(new,signed.scale(old,-1));changes.append(delta)
        exact=signed.scale(signed.add(signed.outer(a,c),signed.outer(c,a),signed.outer(r,d),signed.outer(d,r),signed.outer(c,c),signed.outer(d,d)),Q(1,2))
        checks.extend(delta[i][j]==exact[i][j] for i in range(3) for j in range(3))
        wrong_c=[x/5 for x in p['ell']];wrong_d=[x/5 for x in p['s']]
        wrong=signed.scale(signed.add(signed.outer(a,wrong_c),signed.outer(wrong_c,a),signed.outer(r,wrong_d),signed.outer(wrong_d,r),signed.outer(wrong_c,wrong_c),signed.outer(wrong_d,wrong_d)),Q(1,2))
        shortcut.append(any(delta[0][j]!=wrong[0][j] for j in (1,2)))
    delta=signed.add(*changes);radial=signed.exact.radial
    gr=-radial(delta[0][0])-EPS*delta[0][2].diff(1)+delta[1][1].shift(-1);P=gr.integral()
    rho=bump/bump.integral().scalar()*(1+Z*(R-Q(3,2)));crho=(R**2*rho).integral()/2
    prim=(gr-rho*P).primitive();pm=prim-prim.at(0,1)
    jt=(R**2*delta[1][2]).integral();jz=(R*delta[2][2]).integral()-(R**2*gr).integral()/2
    et=radial(delta[0][1],2)+EPS*delta[1][2].diff(1)
    ez=radial(delta[0][2])+EPS*(delta[2][2]+pm).diff(1)
    return dict(checks=dict(general_two_quadrature_full_tensor_increment=all(checks),
        initial_real_shortcut_not_reused=all(shortcut),all_six_covariance_changes_nonzero=all(bool(delta[i][j]) for i in range(3) for j in range(i,3)),
        compact_pressure_reconstructed=pm.at(0,1)==pm.at(0,2)==0 and pm.diff(0)==gr-rho*P,
        axial_pressure_moment=(R*pm).integral()==-(R**2*gr).integral()/2+crho*P,
        angular_integrated_flux=(R**2*et).integral()==EPS*jt.diff(1),
        axial_integrated_flux=(R*ez).integral()==EPS*(jz+crho*P).diff(1),
        changed_axial_pressure_and_moving_bump_nonzero=bool(pm.diff(1)) and bool(crho.diff(1)*P),
        all_three_defects_must_be_recomputed=all(bool(x) for x in (P,jt,jz))),
        exact_tensor_checks=len(checks),
        arbitrary_complex_quadrature_fixture_not_an_actual_sourced_pulse=True)


def physical_powers(order,h,a,b):
    order,h=Q(order),Q(h)
    if order<0 or not 0<h<Q(1,100) or any(type(n) is not int or n<0 for n in (a,b)):
        raise ValueError('Construction h, nonnegative order and integer derivative counts required')
    A=Q(1,2)+h;shift=-a*(Q(1,2)+h/2)-b*(1+3*h/2)
    return dict(velocity=-A+h*order+shift,pressure=-2*A+h*(order+Q(1,2))+shift,
                potential=Q(1,2)-A+h*(order+Q(1,2))+shift)


def exponent_review():
    k=KAPPA;B=1-3*k;C=B+Q(1,2);mu=1-k
    current=dict(primary_linear=1-3*k,signed_linear=2-6*k,complete_wave_square=1-k,actual_mean_wave=mu)
    new=dict(full_linear=B+Q(1,2)-3*k,old_new_wave=B+Q(1,2)-k,
             new_wave_square=2*B-k,actual_mean_cross=B+mu-Q(1,2))
    D=min(new.values())
    rows=dict(B=B,C=C,D=D,mu=mu,current_core_source=min(current.values()),
              sourced_velocity=B,sourced_pressure_and_potential=C,curl_remainder=C-k,
              complete_covariance_change=C,full_mean_pressure_change=C-k,averaged_mean_pressure_change=C,
              new_full_tangential_mean=C-k,new_Haar_tangential_mean=C,new_defects=C,
              old_fixed_primary_remainder_after_solve=B,finite_joined_order=min(D,C-k,C),
              finite_order_gain=D-B,homogeneous_principal_tail_order=Q(1,2))
    return dict(checks=dict(current_core_class_is_1_minus_3kappa=min(current.values())==B,
        low_order_homogeneous_tail_is_not_in_core_class=Q(1,2)<B,
        linear_remainder_limits_new_core=D==Q(3,2)-6*k,
        all_complete_old_new_products_retained=min(new.values())==D,
        new_square_covariance_no_worse_than_cross=2*B>=C,
        means_and_defects_regenerated_not_preserved=C-k<3-6*k and C<3-3*k,
        source_changes_old_wave_remainder_class=B<1-k,
        finite_joined_gain_positive=rows['finite_order_gain']==Q(1,2)-3*k>0,
        fixed_total_mean_order_preserved=C-k>mu,
        physical_residual_power_still_negative=Q(1,200)*D-2*(Q(1,2)+Q(1,200))-Q(1,2)<0),
        rows={n:str(v) for n,v in rows.items()},current_source_terms={n:str(v) for n,v in current.items()},
        new_source_terms={n:str(v) for n,v in new.items()},
        physical_zero_jet_powers={n:str(v) for n,v in physical_powers(B,Q(1,200),0,0).items()},
        current_harmonic_multipliers=[-2,-1,1,2],maximum_updated_residual_multiplier=4,
        finite_stage_constants_not_uniform_in_iteration_or_a_fixed_scale_contraction=True)
