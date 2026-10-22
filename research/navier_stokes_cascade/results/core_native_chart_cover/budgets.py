"""Existing doubled-box coverage and support-aware native velocity jet bounds."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m


borel=load('native_cover_borel',HERE.parent/'core_borel_cutoff_certificate/budgets.py')
radial=borel.radial;head=radial.head;exact=head.exact
swirl=load('native_cover_swirl',HERE.parent/'global_swirl_envelope/swirl_bounds.py')
axial=load('native_cover_axial',HERE.parent/'global_axial_source_envelope/envelopes.py')


def cover_existing(ell,Xa,Xb,rA,rB,log_margin):
    """No band choice: delta=ell^-6 for the supplied existing accepted floor.

    rA<=sqrt(Xa), rB>=sqrt(Xb); log_margin is the EXISTING leading
    neighborhood size, not a new support or gap neighborhood selection.
    """
    if type(ell) is not int or ell<16:raise ValueError('Existing integer band floor at least sixteen required')
    Xa,Xb,rA,rB,mu=map(exact,(Xa,Xb,rA,rB,log_margin))
    if min(rA,rB,mu)<=0 or not 0<rA*rA<=Xa<Xb<=rB*rB:raise ValueError('Ordered positive annulus, outward root bounds and existing positive log margin required')
    delta=Q(1,ell**6);width=delta*(12/rA+30)
    gates=dict(radial_axis_clearance=6*delta<=rA,implicit_s_barrier=100*delta<=1,
               existing_leading_native_neighborhood=width<=mu)
    return dict(existing_ell=ell,Xa=Xa,Xb=Xb,rA=rA,rB=rB,log_margin=mu,
                delta=delta,native_log_distance=width,gates=gates,
                norm_cover_gates=gates['radial_axis_clearance'] and gates['implicit_s_barrier'],
                band_Rmin=rA/2,band_Rmax=5*rB/2,native_Rlo=rA/4,native_Rhi=5*rB,
                native_Xhi=25*rB*rB/2,
                band_or_neighborhood_selected_or_changed=False,
                actual_growth_or_cone_gaps_verified=False,
                negative_T_extension_verified=False,
                scope='same selected doubled boxes intersected with physical T>=0; all later bands inherit the same inequalities')


def mesh_gate_compatibility(rA,log_margin,Dbox,dchart):
    rA,mu,D,d=map(exact,(rA,log_margin,Dbox,dchart))
    if min(rA,mu,D,d)<=0:raise ValueError('Positive existing mesh, radius and neighborhood data required')
    return dict(radial_axis_clearance=6*d<=D*rA,implicit_s_barrier=100*d<=D,
                existing_native_neighborhood=d*(12/rA+30)<=D*mu,
                scope='sufficient compatibility of the unchanged Dbox*ell^-6<=dchart gate; no scale reselected')


def leading_global_rows(h,Xa,Xv,Xhi,Kaxis,ME,Ahi,log_first,log_second):
    """Velocity-only extension of the old Ceta3/Ceta2/Ceta1 rows to any Xhi.

    The global swirl proof covers all later radii (including the exact
    heat exterior); U and its fixed jets extend by zero after original Xv.
    This does not extend stress/growth/cone floors to the Cartesian hull.
    """
    h,Xa,Xv,Xhi,Kaxis,ME,Ahi,log_first,log_second=map(exact,(h,Xa,Xv,Xhi,Kaxis,ME,Ahi,log_first,log_second))
    if not 0<h<Q(1,100) or Xv<Xa or min(Xa,Xhi)<=0 or min(Kaxis,ME,Ahi,log_first,log_second)<0:
        raise ValueError('Original support and nonnegative retained regional jet inputs required')
    f=swirl.swirl_inputs(Xa,max(Xa,Xhi),Kaxis,ME)
    u=axial.global_bounds(h,1,Xv,Xa,Kaxis,Ahi,log_first,log_second)
    return dict(F=(f['F_C3'],f['F_X_C2'],f['F_XX_C1']),
                U=(u['U_C3'],u['U_X_C2'],u['U_XX_C1']),
                domain_Xhi=Xhi,U_original_support=Xv,
                growth_cone_or_stress_floors_extended=False)


def supported_reconstruct(h,n,m,Xhi,Rlo,Rhi,Frows,Urows,Xv):
    """Reuse the old positive recipe with V=min(Xhi,Xv) for X*U and M.

    The NEW justification extends its estimate to the larger norm hull:
    every derivative of U is supported in [0,Xv]; M retains its full
    axis-to-X value. R bounds and F bounds still cover the ENTIRE hull.
    Pressure-source/moment integrals are never shortened here.
    """
    Xhi,Xv,Rlo,Rhi=map(exact,(Xhi,Xv,Rlo,Rhi))
    if Xhi<=0 or Xv<Rlo*Rlo/2 or 2*Xhi>Rhi*Rhi:raise ValueError('Whole hull radius bounds and original U support must be supplied')
    volume=min(Xhi,Xv)
    out=radial.reconstruct(h,n,m,volume,Rlo,Rhi,Frows,Urows)
    return dict(out,domain_Xhi=Xhi,primitive_volume=volume,U_original_support=Xv,
                scope='support-aware majorant on the whole hull; only X*U and primitive upper volume use min(Xhi,Xv), never future pressure memory')


def join(h,cover,Xa,Xv,Kaxis,ME,Ahi,log_first,log_second,inner,cutoff,first_data,p_lower,rp_lower,geometry):
    expected=cover_existing(cover['existing_ell'],cover['Xa'],cover['Xb'],cover['rA'],cover['rB'],cover['log_margin'])
    if cover!=expected:raise ValueError('Derived cover or scope fields cannot be replaced by caller assertions')
    if not cover['norm_cover_gates']:raise ValueError('Unchanged band fails doubled-box norm coverage gates')
    if 'Rmin' in geometry or 'Rmax' in geometry:raise ValueError('Band radius enclosures are derived from the same selected boxes')
    Xa,Xv=map(exact,(Xa,Xv))
    if Xa!=cover['Xa']:raise ValueError('Leading row anchor must be the SAME active annulus Xa')
    if exact(first_data['x_velocity'])!=Xv:raise ValueError('First pressure repair must keep the SAME original full source support Xv')
    rows=leading_global_rows(h,Xa,Xv,cover['native_Xhi'],Kaxis,ME,Ahi,log_first,log_second)
    first=radial.finalized_first(h,inner,cutoff,first_data,p_lower,rp_lower)
    args=(cover['native_Xhi'],cover['native_Rlo'],cover['native_Rhi'])
    leading=supported_reconstruct(h,0,2,*args,rows['F'],rows['U'],Xv)
    current=supported_reconstruct(h,1,1,*args,first['jets']['F'][:2],first['jets']['U'][:2],Xv)
    K=borel.cutoff_bounds(2)['K'][2]
    g=dict(geometry,Rmin=cover['band_Rmin'],Rmax=cover['band_Rmax'])
    dep=borel.comparison_dependence(h,leading['through'][1],leading['through'][2],current['through'][1],g)
    return dict(cover=cover,rows=rows,first_state=first,leading=leading,first=current,
                geometry=g,cutoff_K2=K,dependence=dep,Cphi=borel.comparison_value(dep,K),
                existing_value_neighborhood_passes=cover['gates']['existing_leading_native_neighborhood'],
                infinite_Borel_schedule_and_actual_selector_verified=False,
                scope='conditional global velocity-row and unchanged-box cover join; retained regional premises, cutoff identity and all value/support/propagator gates remain')


def fixture():
    selector=load('native_cover_original_selector',HERE.parent/'core_primary_tail_recheck/proof_checks.py')
    p=selector.fixture();selected=selector.select_tail(p)
    cover=cover_existing(selected['ell_tail'],1,64,1,8,Q(1,1000))
    old=radial.fixture();data=dict(old['first_data'],x_velocity=Q(16),leading_u5=Q(1000),leading_u3=Q(1000),leading_f3=Q(512))
    geometry={k:v for k,v in old['geometry'].items() if k not in ('Rmin','Rmax')}
    inputs=dict(h=old['h'],cover=cover,Xa=Q(1),Xv=Q(16),Kaxis=Q(7),ME=Q(128),Ahi=Q(13),
                log_first=Q(101),log_second=Q(303),inner=old['inner'],cutoff=old['cutoff'],first_data=data,
                p_lower=Q(1,2),rp_lower=Q(1),geometry=geometry)
    return dict(selector_inputs=p,selected=selected,inputs=inputs,joined=join(**inputs))
