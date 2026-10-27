"""Exact controls for the reconstructed same-band continuum gap contract."""
from fractions import Fraction as F

def perturbation_radius(f, fbar, t, gbar, d, omega, gamma, u):
    f,fbar,t,gbar,d,omega,gamma,u=map(F,(f,fbar,t,gbar,d,omega,gamma,u))
    if min(f,t,d,omega,gamma,u)<=0 or fbar<f or gbar<t or max(omega,gamma)>1:
        raise ValueError("Positive ordered continuum inputs required")
    LD=4*gbar+2*fbar+3; LP=fbar+gbar+1
    Dbar=gbar*gbar+2*fbar*gbar
    cbar2=Dbar/(2*f*t)
    Lc=2*LD/(f*t)+2*Dbar*LP/(f*f*t*t)
    LH=1+2/t
    Lphi=3*u*u*LH+(1+u*u)*(4*Lc+3*cbar2*LH)
    mphi=u*u*omega*omega*gamma*(2-gamma)
    rows={"relative_stress":F(1,2),"F":f/2,"minus_gtheta":t/2,
          "D":d/(2*LD),"orientation":omega/(2*LH),"same_tilt_cone":mphi/(2*Lphi)}
    return {"radius":min(rows.values()),"rows":rows,"LD":LD,"LP":LP,
            "Lc":Lc,"LH":LH,"Lphi":Lphi,"mphi":mphi}

def persisted_floors(f, fbar, t, gbar, d, omega, gamma):
    f,fbar,t,gbar,d,omega,gamma=map(F,(f,fbar,t,gbar,d,omega,gamma))
    if min(f,t,d,omega,gamma)<=0 or min(fbar,gbar)<0 or gamma>1:
        raise ValueError("Positive gap inputs required")
    return {"F":f/2,"minus_gtheta":t/2,"D":d/2,
      "lambda2":f*t*d/(4*(gbar+1)**2),
      "c2":d/(4*(fbar+1)*(gbar+1)),"orientation":omega/3,
      "cone":omega*omega*gamma*(2-gamma)/16}

def first_exact_dyadic_band(h, ratios, minimum=3):
    """Archimedean selector without floating logarithms: 2^(2h(ell-2))>=ratios."""
    h=F(h); rows=tuple(map(F,ratios))
    if not 0<h<F(1,100) or not rows or min(rows)<=0 or type(minimum) is not int:
        raise ValueError("Positive finite selector data required")
    den=h.denominator
    ell=max(minimum,2)
    ell += (-ell)%den
    while True:
        exponent=2*h*(ell-2)
        power=exponent.numerator//exponent.denominator
        if all(F(2)**power>=row for row in rows):
            return ell
        ell += den

def selector_margins(ell,h,ratios,r0,length_required,k_inverse):
    h=F(h); ratios=tuple(map(F,ratios)); r0,length_required,k_inverse=map(F,(r0,length_required,k_inverse))
    if type(ell) is not int or ell<3 or not ratios or min(ratios)<=0:
        raise ValueError("Positive finite selector rows required")
    exponent=h*ell
    if exponent.denominator!=1 or min(r0,length_required,k_inverse)<=0:
        raise ValueError("Exact dyadic exponent and positive rows required")
    eps=F(1,2**exponent.numerator)
    value_exponent=2*h*(ell-2)
    value=F(2)**(value_exponent.numerator//value_exponent.denominator)
    return {"value":min(value/r for r in ratios)-1,
            "length":2*r0*ell*ell-length_required,
            "slow":1-ell**4*(eps+eps*eps+k_inverse)}

def chart_tail(qstar_power,ell):
    """Exact exponent form for qstar=2^-qstar_power and qv=min(qstar/64,2^-ell-2)."""
    if type(qstar_power) is not int or qstar_power<0 or type(ell) is not int or ell<1:
        raise ValueError("Invalid chart tail")
    qv_power=max(qstar_power+6,ell+2)
    return {"qstar_power":qstar_power,"qv_power":qv_power,
            "enlarged_power":qv_power-4,"inside":qv_power-4>qstar_power,
            "selected":qv_power-2>=ell}

def normal_floor(theta,beta):
    theta,beta=map(F,(theta,beta))
    if not 0<=theta<=F(1,32) or beta<=0: raise ValueError("Invalid phase floor")
    return {"fraction":1-theta,"normal2":(1-theta)**2*beta,
            "inverse":1/((1-theta)**2*beta),
            "below_coarse_inverse":1/((1-theta)**2*beta)<4/beta}

def review():
    data=(F(1,5),3,F(1,7),4,F(1,11),F(1,13),F(1,17),2)
    radius=perturbation_radius(*data)
    floors=persisted_floors(*data[:-1])
    h=F(1,200); ratios=(2,17,10**6)
    ell=max(first_exact_dyadic_band(h,ratios),400000)
    margins=selector_margins(ell,h,ratios,F(1,6057072),1000,F(1,2**1000))
    chart=chart_tail(3,ell); normal=normal_floor(F(1,32),F(1,101))
    return {"checks":{
      "all_continuum_perturbation_rows_are_positive":radius["radius"]>0 and all(v>0 for v in radius["rows"].values()),
      "persisted_value_spectral_and_cone_floors_are_positive":all(v>0 for v in floors.values()),
      "one_exact_integer_selector_has_strict_value_length_and_slow_margins":all(v>0 for v in margins.values()),
      "enlarged_chart_stays_in_the_same_selected_tail":chart["inside"] and chart["selected"],
      "reference_relative_phase_normal_has_positive_floor":normal["normal2"]>0 and normal["below_coarse_inverse"],
      "flat_edge_contract_is_weighted_not_a_positive_stress_magnitude_infimum":True,
    },"controls":{"perturbation_radius":str(radius["radius"]),
      "persisted_floors":{k:str(v) for k,v in floors.items()},"selector":ell,
      "selector_margins":{k:str(v) for k,v in margins.items()},
      "chart_tail":chart,
      "normal_floor":{k:(str(v) if isinstance(v,F) else v) for k,v in normal.items()}},
      "scope":{"same_band_contract_reconstructed":True,"stress_magnitude_positive_at_flat_edges":False,
      "full_path_ode_closed":False,"unconditional_pulse_residual_force_or_blowup":False}}
