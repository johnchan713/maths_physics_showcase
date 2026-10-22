"""Original physical representatives and a norm-only tangential floor join."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


native=load('physical_rep_native',HERE.parent/'core_native_chart_cover/budgets.py')
tangent=load('physical_rep_tangent',HERE.parent/'core_tangential_normal_closure/budgets.py')
exact=tangent.exact


def verify_representative(h,ell,Xa,Xb,index,point,root_s):
    """Check SUPPLIED exact chart data; never move or choose a representative.

    s=root_s^(2*denominator(h)) permits exact rational chart checks.
    This checks the physical/native clauses of A_ell. The original qbig
    membership remains a separate defining clause, not evaluated here.
    """
    h,Xa,Xb,t=map(exact,(h,Xa,Xb,root_s))
    if not 0<h<Q(1,100) or type(ell) is not int or ell<16 or not 0<Xa<Xb or t<=0:raise ValueError('Fixed h, existing band, ordered native shell and positive chart root required')
    if len(index)!=3 or any(type(a) is not int for a in index) or len(point)!=3:raise ValueError('Original integer mesh index and supplied three-coordinate representative required')
    R,Z,T=map(exact,point);delta=Q(1,ell**6)
    s=t**(2*h.denominator);sD=t**(h.denominator-2*h.numerator)
    eta=Z/sD;X=R*R/(2*s)
    gates=dict(in_original_closed_box=all(abs(v-delta*a)<=delta for v,a in zip((R,Z,T),index)),
               physical_T_nonnegative=T>=0,positive_cylindrical_radius=R>0,
               original_dyadic_band=Q(1,2)<=s<=2,original_native_active_shell=Xa<=X<=Xb,
               exact_physical_chart=T==s*(1-eta*eta),physical_eta=abs(eta)<=1)
    return dict(gates=gates,all_physical_native_clauses_pass=all(gates.values()),s=s,eta=eta,X=X,
                existing_point=(R,Z,T),existing_index=tuple(index),representative_selected_or_changed=False,
                original_qbig_membership_numerically_checked=False)


def norm_cover(ell,Xa,Xb,rA,rB):
    """Same physical doubled-box norm hull; no value neighborhood is selected."""
    if type(ell) is not int or ell<16:raise ValueError('Existing integer band at least sixteen required')
    Xa,Xb,rA,rB=map(exact,(Xa,Xb,rA,rB))
    if min(rA,rB)<=0 or not 0<rA*rA<=Xa<Xb<=rB*rB:raise ValueError('Original annulus and outward root bounds required')
    d=Q(1,ell**6)
    return dict(existing_ell=ell,Xa=Xa,Xb=Xb,rA=rA,rB=rB,delta=d,
                gates=dict(radial_axis_clearance=6*d<=rA,implicit_s_barrier=100*d<=1),
                band_Rmin=rA/2,band_Rmax=5*rB/2,native_Rlo=rA/4,native_Rhi=5*rB,native_Xhi=25*rB*rB/2,
                physical_T_only=True,value_neighborhood_selected_or_required_for_this_norm_cover=False,
                growth_or_cone_floors_extended_to_norm_hull=False,negative_T_extension_certified=False)


def floor_join(cover,Lead1,lambda_min,geometry):
    """Reference gaps are needed at the ORIGINAL active representative only.

    Leading first tuple and full-background norms must cover the physical
    norm hull. This does not certify the entire wave/cone construction.
    """
    checked=norm_cover(cover['existing_ell'],cover['Xa'],cover['Xb'],cover['rA'],cover['rB'])
    if checked!=cover:raise ValueError('Supplied existing physical norm cover was altered')
    if not all(cover['gates'].values()):raise ValueError('Original physical doubled-box norm coverage fails')
    g=dict(geometry)
    if 'Rmin' in g or 'Rmax' in g:raise ValueError('Band radii come from the same original cover')
    if exact(g['Dbox'])<6:raise ValueError('Doubled-box Euclidean displacement requires the retained safe Dbox>=6')
    r=tangent.comparison(Lead1,cover['band_Rmin'],cover['band_Rmax'],**g)
    lam=exact(lambda_min)
    if not 0<lam<=r['lambda_max']:raise ValueError('Consistent positive growth bound at original active representatives required')
    beta=lam/(4*(1+r['tilt_upper']**2)**2);S=Q(cover['existing_ell']**2)
    return dict(cover=cover,comparison=r,beta=beta,existing_S=S,
                normal_margin=beta*S*S-4*r['Ctan']**2,
                supplied_normal_margin_passes=beta*S*S>=4*r['Ctan']**2,
                independent_native_value_neighborhood_needed_for_this_floor=False,
                negative_T_representative_or_extension_needed_for_this_floor=False,
                original_slow_gate_evaluated=False,actual_closed_annulus_gap_bounds_numerically_verified=False,
                actual_existing_band_floor_certified=False,scales_supports_or_representatives_changed=False)
