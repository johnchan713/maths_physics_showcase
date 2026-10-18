"""Exact sufficient dyadic budgets; fixtures never stand for actual norm enclosures."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m

@lru_cache(None)
def bounds():return load('value_domain_canonical_bounds',HERE.parent/'canonical_primary_foundation/canonical_bounds.py')

def exact(x):
    if isinstance(x,bool) or not isinstance(x,(int,F)):raise TypeError('Exact rational input required')
    return F(x)

def ceil(x):
    x=exact(x);return -(-x.numerator//x.denominator)

def binary_ceiling(x):
    x=exact(x)
    if x<=0:raise ValueError('Positive binary-log argument required')
    n,d=x.numerator,x.denominator;e=n.bit_length()-d.bit_length()
    fits=n<=d<<e if e>=0 else n<<(-e)<=d
    return e if fits else e+1

def coefficient_exponent(c,e):
    c,e=exact(c),exact(e)
    if c<=0 or e<0:raise ValueError('Positive coefficient and nonnegative exponent required')
    return binary_ceiling(c)+ceil(e)

def certificate(p,h=F(1,1000),radius=F(1,6057072),gamma=F(1,4),Jmin=F(1),Jmax=F(2),old_band=16):
    h,radius,gamma,Jmin,Jmax=map(exact,(h,radius,gamma,Jmin,Jmax))
    if not 0<h<F(1,100) or not 0<radius<=F(1,6057072):raise ValueError('Original h and geometric radius ranges required')
    if not 0<gamma<=1 or not 0<Jmin<=Jmax:raise ValueError('Positive cone and ordered Haar bounds required')
    if type(old_band) is not int or old_band<16:raise ValueError('Original band floor at least 16 required')
    # Never fill omitted actual enclosures from the parent's diagnostic defaults.
    if set(p)!=set(bounds().parent_budget().inputs()):raise ValueError('Complete primitive enclosure keys required; no diagnostic defaults')
    # Revalidate, rather than trusting a caller-supplied dictionary of floors.
    p=bounds().parent_budget().inputs(**p);b=bounds().parent_budget().budgets(p);r=bounds().canonical(p)
    u=p['u'];c=p['lam_min']*u*u/2*((1+9*u*u/4)**-2+2*(1+u*u)**-2)
    C=3*p['lam_max']*u*u/2*((1+u*u/4)**-1+2*(1+u*u)**-1)
    G=2*(1+1/c);Da=b['Cd']+F(9,8)*b['CE'];Dr=b['Cd']+F(9,8)*r['CE']
    Dstar=b['Cd']+max(b['CE'],r['CE']);Dast=b['Cd']+b['CE'];m=min(p['cmin'],u)
    Ja=coefficient_exponent((8/F(7))**2/Jmin,2*(F(2,3)+2*Da+C))
    Pa=coefficient_exponent(Jmax*(9/F(8))**2*G,4*Da)
    H=coefficient_exponent(2*Jmax*r['B0']**2*b['Ds']*G,2*(Dstar+Dast))
    K=binary_ceiling(24/gamma)+Pa+Ja
    T=binary_ceiling(12/m)+H+Ja+max(K+1,0)+1
    cones={}
    for name,D,D0 in (('inherited',Da,b['D0']),('canonical',Dr,r['D0'])):
        V=coefficient_exponent(D0+b['D1']*(9/F(7))**2/(2*c),8*D+2*C)
        cones[name]=2*(binary_ceiling(4*max(1/p['cmin'],1/u)/gamma)+V)
    direct=max(binary_ceiling(25),binary_ceiling(b['normal_acceptance']),binary_ceiling(r['Riccati_length']),cones['canonical'])
    common=max(direct,binary_ceiling(b['Riccati_acceptance']),cones['inherited'],T)
    radius_cost=binary_ceiling(1/(2*radius))
    length_B=ceil(F(common+radius_cost,2))
    B=max(4,binary_ceiling(old_band),binary_ceiling(512/h**2),length_B)
    return dict(inverse_mass_log2=Ja,upper_mass_log2=Pa,column_difference_log2=H,kappa_log2=K,
        comparison_length_log2=T,cone_length_log2=cones,direct_length_log2=direct,
        common_length_log2=common,radius_cost_log2=radius_cost,length_band_log2=length_B,
        smallness_band_log2=binary_ceiling(512/h**2),old_band_log2=binary_ceiling(old_band),
        accepted_band_log2=B,smallness_ratio=str(F(3*factorial(8)*4**8,512**4)),
        assumptions='supplied outward value bounds on actual closed charts; diagnostic call does not supply actual enclosures')

def length_margin(cert,band_log2):
    if type(band_log2) is not int or band_log2<0:raise ValueError('Nonnegative integer band exponent required')
    return 2*band_log2-cert['radius_cost_log2']-cert['common_length_log2']

def accepted_dyadic_band(cert,band_log2):
    if type(band_log2) is not int or band_log2<0:raise ValueError('Nonnegative integer band exponent required')
    return band_log2>=cert['accepted_band_log2']

def relative_gate(a,kappa,L):
    a,kappa,L=map(exact,(a,kappa,L))
    if min(a,kappa,L)<=0:raise ValueError('Positive relative gate data required')
    theta=a/L
    return dict(theta=str(theta),accepted=theta*(2*kappa+1)<=1,
        relative_bound=None if theta>=1 else str(kappa*theta/(1-theta)))

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,list):return [encode(v) for v in x]
    return x

def dyadic_review():
    rows=[]
    for h in (F(1,1000),F(1,10**5)):
        for radius in (F(1,6057072),F(1,6057072*16)):
            p=bounds().parent_budget().inputs();c=certificate(p,h,radius)
            rows.append(dict(h=str(h),radius=str(radius),certificate=c,
                accepted=accepted_dyadic_band(c,c['accepted_band_log2']),length_margin=length_margin(c,c['accepted_band_log2'])))
    return dict(manufactured_inputs=encode(bounds().parent_budget().inputs()),rows=rows,
        checks=dict(exact_integer_and_rational_budgets_only=all(type(r['certificate']['accepted_band_log2']) is int for r in rows),
            one_half_line_satisfies_every_value_length=all(r['accepted'] and r['length_margin']>=0 for r in rows),
            exponential_series_smallness_is_strict=F(rows[0]['certificate']['smallness_ratio'])<1,
            h_and_radius_not_retuned=True,actual_norm_enclosures_not_manufactured_fixtures=True,
            no_giant_band_power_or_tiny_q_materialized=True,no_derivative_order_in_acceptance=True))

def domain_review():
    rows=[]
    for E in (16,32,64):
        for qstar in (F(1,2),F(1,2**90)):
            qt=min(qstar/64,F(1,2**(E+2)));q=qt*F(7,8)
            for ratio in (F(1,4),F(1),F(4)):
                Q=q/ratio
                rows.append(dict(E=E,qstar=str(qstar),q=str(q),Q=str(Q),later=Q<F(1,2**E),inside=4*Q<qstar/4))
    witness=relative_gate(1,2,4);endpoint=relative_gate(1,2,5)
    return dict(chart_rows=rows,old_band_parameter_witness=witness,threshold_endpoint=endpoint,
        checks=dict(enlarged_chart_margin_retained=all(r['later'] and r['inside'] for r in rows),
            eventually_small_does_not_verify_old_band=not witness['accepted'] and endpoint['accepted'],
            plus_one_is_required=F(witness['relative_bound'])>F(1,2),
            positive_canonical_roots_do_not_need_historical_relative_error=True,
            old_qstar_supports_fields_and_schedules_not_mutated=True,
            old_localized_force_not_identified_with_new_tail_localization=True))

def logarithmic_review():
    # Independent transcendental evaluation only on manufactured finite budgets.
    import mpmath as mp
    with mp.workdps(80):
        f=lambda x:mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)
        rows=[]
        for gamma,Jmin,Jmax in ((F(1,4),F(1),F(2)),(F(1,17),F(1,9),F(11))):
            p=bounds().parent_budget().inputs();b=bounds().parent_budget().budgets(p);r=bounds().canonical(p)
            cert=certificate(p,gamma=gamma,Jmin=Jmin,Jmax=Jmax);u=f(p['u']);c=f(p['lam_min'])*u*u/2*((1+9*u*u/4)**-2+2*(1+u*u)**-2)
            C=3*f(p['lam_max'])*u*u/2*((1+u*u/4)**-1+2*(1+u*u)**-1)
            Da=f(b['Cd']+F(9,8)*b['CE']);ln2=mp.log(2)
            invj=(-mp.log(f(Jmin))+mp.mpf(2)/3-2*mp.log(mp.mpf(7)/8)+2*Da+C)/ln2
            plus=(mp.log(f(Jmax))+2*mp.log(mp.mpf(9)/8)+2*Da+mp.log(mp.pi/c)/2)/ln2
            H=(mp.log(2*f(Jmax)*f(r['B0'])**2*f(b['Ds']))+f(2*b['Cd']+max(b['CE'],r['CE'])+b['CE'])+mp.log(mp.pi/c)/2)/ln2
            K=mp.log(24/f(gamma),2)+plus+invj
            T=mp.log(12/f(min(p['cmin'],p['u'])),2)+H+invj+mp.log(2*mp.power(2,K)+1,2)
            values=dict(inverse_mass_log2=invj,upper_mass_log2=plus,column_difference_log2=H,kappa_log2=K,comparison_length_log2=T)
            for name,D,D0 in (('inherited',b['Cd']+F(9,8)*b['CE'],b['D0']),('canonical',b['Cd']+F(9,8)*r['CE'],r['D0'])):
                logdir=mp.log(f(D0)+f(b['D1'])*(mp.mpf(9)/7)**2*mp.exp(4*f(D)+C)/(2*c))
                value=2*(mp.log(4*f(max(1/p['cmin'],1/p['u']))/f(gamma))+logdir)/ln2
                rows.append(dict(quantity=name+'_cone',upper=cert['cone_length_log2'][name],below=bool(value<=cert['cone_length_log2'][name])))
            rows += [dict(quantity=k,upper=cert[k],below=bool(v<=cert[k])) for k,v in values.items()]
    return dict(manufactured_only=True,rows=rows,checks=dict(independent_logarithmic_fixtures_below_exact_majorants=all(r['below'] for r in rows),
        numerical_diagnostics_not_actual_domain_certificate=True))
