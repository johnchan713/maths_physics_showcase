"""Native directional transport and a retained-power SAME-band floor."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('native_floor_old',HERE.parent/'core_tangential_normal_closure/budgets.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
exact=old.exact

def nonnegative(x):
    x=exact(x)
    if x<0:raise ValueError('Nonnegative bound required')
    return x

def directional_envelope(a,b):
    """Exact max of 2*a*t+b*(1-t*t), 0<=t<=1."""
    a,b=nonnegative(a),nonnegative(b)
    return 2*a if a>=b else b+a*a/b

def leading_direction(h,Xhi,*,value,radial,angular):
    h,Xhi=exact(h),nonnegative(Xhi)
    if not 0<h<Q(1,100):raise ValueError('Original 0<h<1/100 required')
    v,r,e=map(nonnegative,(value,radial,angular))
    A=Q(1,2)+h
    return directional_envelope(A*v+Xhi*r,e)

def comparison(h,leading_E,leading_U,*,Rmin,Rmax,tilt_upper,lambda_max,gmin,Llow,Lhigh,Dbox):
    h=exact(h)
    if not 0<h<Q(1,100):raise ValueError('Original 0<h<1/100 required')
    lo,hi,U,lam,g,Llo,Lhi,D=map(exact,(Rmin,Rmax,tilt_upper,lambda_max,gmin,Llow,Lhigh,Dbox))
    if min(lo,lam,g,Llo,Lhi)<=0 or hi<lo or Lhi<Llo or min(U,D)<0:raise ValueError('Original positive annular and reference geometry bounds required')
    for data in (leading_E,leading_U):
        if set(data)!={'value','radial','angular'}:raise ValueError('Leading native value, X and eta derivative bounds required')
    Xhi=2*hi*hi
    DE=leading_direction(h,Xhi,**leading_E);DU=leading_direction(h,Xhi,**leading_U)
    carriers=old.head.phase.carrier_bounds(lam,U,g,hi,Llo)
    p,pz=carriers['theta'],carriers['axial'];B=max(Q(1),lam)
    K=4/(1-2*h);weighted_carrier=p/lo+pz
    # q^h <= 2 epsilon, and max_eta [2|eta|(1+A+Xhi)+1-eta^2]
    # is 2(1+A+Xhi). No radial streamfunction component is discarded.
    Hlead=K*(p*DE/lo+pz*DU)
    Hcorr=K*weighted_carrier*4*(1+Q(1,2)+h+Xhi)
    tilt=2*B*U/(Llo*g);rounding=1/lo;box=(p-1)*D/lo**2
    coeff=dict(A=tilt+Lhi*Hlead,B=rounding,C=box+Lhi*Hcorr)
    return dict(coefficients=coeff,HZ_leading=Hlead,HZ_correction_per_epsilon=Hcorr,
        groups=dict(tilt=tilt,rounding=rounding,box=box),carriers=carriers,Lhigh=Lhi,
        native_direction_E=DE,native_direction_U=DU,native_Xlo=lo*lo/8,native_Xhi=Xhi,
        lambda_max=lam,tilt_upper=U,h=h,
        retained_premise='same COMPLETE normalized Euler diagonal through native order1 for EVERY positive coefficient; same original cover/reference gaps',
        actual_leading_norms_diagonal_or_same_band_floor_certified=False)

def normal_gate(result,lambda_min,S,epsilon,k_inverse):
    lam,S,e,ki=map(exact,(lambda_min,S,epsilon,k_inverse))
    if not 0<lam<=result['lambda_max'] or S<1 or min(e,ki)<0:raise ValueError('Consistent original growth floor and nonnegative scale data required')
    U=result['tilt_upper'];beta=lam/(4*(1+U*U)**2)
    g=result['groups'];c=result['coefficients'];slow=S*S*(e+e*e+ki)<=1
    raw=g['tilt']/S+g['rounding']*ki+g['box']/S**3+result['Lhigh']*S*e*(result['HZ_leading']+e*result['HZ_correction_per_epsilon'])
    reduced=c['A']/S+c['B']/S**2+c['C']/S**3
    return dict(beta=beta,raw_error=raw,retained_power_error=reduced,
        gates=dict(original_slow_gate=slow,tangential_margin=4*reduced**2<=beta),
        raw_dominated_if_slow=not slow or raw<=reduced,
        collapsed_same_program_error=sum(c.values(),Q(0))/S,
        inverse_normal_square=4/beta,existing_S=S,
        actual_profile_norms_band_or_field_certified=False,new_scale_or_length_selected=False)

def uniform_slow_bound(h,ell):
    h,ell=exact(h),exact(ell)
    if not 0<h<Q(1,100) or ell<=0:raise ValueError('Original h and supplied positive band required')
    C=Q(3*factorial(16)*4**16)
    bound=C/(h**16*ell**12)
    universal=C*Q(1,100)**8/Q(512)**12
    return dict(original_band_gate=ell>=512/h**2,slow_upper=bound,
        universal_upper_if_original_band_gate=universal,
        uniform_gate_from_original_bound=ell>=512/h**2 and universal<1,
        carrier='ORIGINAL k=ceil(epsilon^-1/2), not a replacement carrier',
        actual_supplied_band_identified=False,new_band_selected=False)

def certificate():
    return dict(parent_commit='9151274fbad956a8e2cabdf99b1e06af919ceabf',
        identity='partial_Z F=s^-1/(R L) D_h E, partial_Z G=s^-1/L D_h U; D_h=2eta(Euler_q-X partial_X-A)+(1-eta^2)partial_eta, A=1/2+h',
        leading='Euler_q E0=Euler_q U0=0; use native E0,U0 value/X/eta jets, no radial primitive M as an independent floor norm',
        correction='same complete native Euler diagonal gives |delta E,delta U|_1<=q^h; D_h correction <=2(1+A+Xhi)q^h',
        heat='E0=C(2X)^(-A) H(2(1-eta^2)/X), U0=0 => D_h E0=D_h U0=0 in the actual retained heat exterior',
        retained_power_gate='4[A/S+B/S^2+C/S^3]^2<=beta; do not collapse distinct powers without need',
        slow='S^2(epsilon+epsilon^2+k^-1)<=3*16!*4^16/(h^16 ell^12); supplied ell>=512/h^2 gives uniform upper <10^-24',
        same_scale=True,actual_leading_input_diagonal_same_band_full_residual_or_blowup_certified=False)

def certify_actual_inputs():
    raise NotImplementedError('Actual leading E0/U0 primitive first-jet enclosures, original complete Euler diagonal and reference normal/root/cone/freezing gates remain to be verified on the unchanged band')
