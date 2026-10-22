"""Exact chart witnesses, native log gates and full-pressure support controls."""
from copy import deepcopy
from fractions import Fraction as Q
from functools import lru_cache
from budgets import cover_existing,mesh_gate_compatibility,leading_global_rows,supported_reconstruct,join,fixture,radial,borel


def geometry_review():
    c=cover_existing(16,Q(1,2),2,Q(1,2),2,Q(1,100));delta=c['delta'];checks={};rows=[]
    for R0 in (Q(1),Q(2)):
        for eta0 in (Q(-1),Q(0),Q(3,5),Q(1)):
            for sign in (-1,1):
                t=1+sign*delta/1000;s=t**200;Z=eta0*t**99;T=(1-eta0*eta0)*s;R=R0+sign*delta
                Z0=eta0;T0=1-eta0*eta0;X=R*R/(2*s)
                log_upper=2*abs(R-R0)/min(R,R0)+200*abs(t-1)/min(t,Q(1))
                key=f'R{R0}_eta{eta0}_sign{sign}'
                checks[key+'_exact_physical_chart_equation']=s-Z*Z*t*t==T
                checks[key+'_within_selected_doubled_box_displacement']=max(abs(R-R0),abs(Z-Z0),abs(T-T0))<=3*delta
                checks[key+'_s_and_native_norm_hull_bounds']=Q(1,4)<=s<=4 and c['band_Rmin']<=R<=c['band_Rmax'] and c['native_Rlo']**2/2<=X<=c['native_Xhi']
                checks[key+'_rational_native_log_distance_upper']=log_upper<=c['native_log_distance']
                rows.append(dict(R0=str(R0),eta0=str(eta0),sign=sign,physical_T_nonnegative=T>=0))
    checks['entire_norm_hull_includes_radii_outside_active_annulus']=c['native_Xhi']>c['Xb'] and c['native_Rlo']**2/2<c['Xa']
    checks['norm_hull_does_not_promote_growth_or_cone_gaps']=not c['actual_growth_or_cone_gaps_verified']
    checks['all_geometry_gates_pass_for_the_manufactured_existing_floor']=all(c['gates'].values())
    return dict(checks=checks,rows=rows,scope='exact rational h=1/200 physical chart witnesses including T=0; continuum cover follows derivative/continuation estimates, not these finite witnesses')


def cover_review():
    f=fixture();c=f['joined']['cover'];ell=c['existing_ell'];delta=c['delta'];checks={}
    checks['all_original_selector_fixture_checks_still_pass']=all(f['selected']['checks'].values())
    for later in (ell,ell+1,ell+17,ell+100000):
        d=cover_existing(later,1,64,1,8,Q(1,1000))
        checks[f'ell{later}_same_all_later_cover_and_native_neighborhood']=all(d['gates'].values()) and d['native_log_distance']<=c['native_log_distance']
    bad=cover_existing(ell,delta*delta,64,delta,8,Q(1,1000))
    narrow=cover_existing(ell,1,64,1,8,c['native_log_distance']/2)
    checks['existing_axis_clearance_failure_is_reported_not_reselected']=not bad['gates']['radial_axis_clearance'] and bad['existing_ell']==ell
    checks['existing_native_neighborhood_failure_is_separate_from_norm_cover']=narrow['norm_cover_gates'] and not narrow['gates']['existing_leading_native_neighborhood']
    checks['supplied_band_and_neighborhood_are_never_changed']=not c['band_or_neighborhood_selected_or_changed']
    compat=mesh_gate_compatibility(1,Q(1,100),6,Q(1,100000))
    checks['old_mesh_gate_can_imply_all_three_cover_inequalities']=all(v for k,v in compat.items() if k!='scope')
    return dict(checks=checks,existing_ell=str(ell),scope='same old selector fixture and supplied native enlargement; no actual historical band constant is numerically certified')


def support_review():
    f=fixture();o=f['joined'];c=o['cover'];checks={};rows=[]
    args=(c['native_Xhi'],c['native_Rlo'],c['native_Rhi'])
    first=o['first_state']['jets']
    for n,m,fr,ur in ((0,2,o['rows']['F'],o['rows']['U']),(1,1,first['F'][:2],first['U'][:2])):
        old=radial.reconstruct(f['inputs']['h'],n,m,*args,fr,ur)
        new=supported_reconstruct(f['inputs']['h'],n,m,*args,fr,ur,f['inputs']['Xv'])
        checks[f'n{n}_all_support_aware_rows_dominated_by_old_whole_hull_bound']=all(new['raw'][name][a]<=v for name,row in old['raw'].items() for a,v in row.items())
        checks[f'n{n}_radial_and_primitive_value_norms_strictly_improve']=new['raw']['W'][(0,0)]<old['raw']['W'][(0,0)] and new['raw']['J'][(0,0)]<old['raw']['J'][(0,0)]
        checks[f'n{n}_full_hull_is_retained_for_swirl_and_axial_jets']=new['raw']['E']==old['raw']['E'] and new['raw']['U']==old['raw']['U']
        checks[f'n{n}_only_axial_primitive_volume_uses_common_support']=new['primitive_volume']==16<c['native_Xhi'] and new['domain_Xhi']==c['native_Xhi']
        rows.append(dict(coefficient=n,old_tuple=str(old['through'][m]),supported_tuple=str(new['through'][m])))
    # When the chart ends before the source support, pressure moments still
    # integrate to the original Xv, even though M(X) only integrates to X.
    inputs=deepcopy(f['inputs']);inputs['Xv']=Q(1600);inputs['first_data']['x_velocity']=Q(1600)
    full=join(**inputs)
    wrong_data=dict(inputs['first_data'],x_velocity=c['native_Xhi'])
    truncated=radial.finalized_first(inputs['h'],inputs['inner'],inputs['cutoff'],wrong_data,inputs['p_lower'],inputs['rp_lower'])
    checks['chart_prefix_primitive_does_not_truncate_future_pressure_source']=full['leading']['primitive_volume']==c['native_Xhi']<1600 and full['first_state']['repair']['source']['Q03']>truncated['repair']['source']['Q03'] and full['first_state']['repair']['source']['Q13']>truncated['repair']['source']['Q13']
    checks['full_future_pressure_changes_the_first_moment_repair']=full['first_state']['repair']['repair']['coefficient_Ceta3']>truncated['repair']['repair']['coefficient_Ceta3']
    return dict(checks=checks,rows=rows,scope='positive exact norm recipes with actual support algebra; not evaluations of selected giant F,U or pressure profiles')


def join_review():
    f=fixture();o=f['joined'];c=o['cover'];g=o['geometry'];checks={}
    checks['derived_native_hull_exactly_meets_parent_coverage_guard']=4*c['native_Rlo']**2==g['Rmin']**2 and c['native_Xhi']==2*g['Rmax']**2 and c['native_Rhi']**2==2*c['native_Xhi']
    checks['all_radius_enclosures_come_from_same_active_shell']=g['Rmin']==c['rA']/2 and g['Rmax']==5*c['rB']/2
    checks['leading_F_rows_are_independent_of_finite_heat_interface_endpoint']=leading_global_rows(Q(1,200),1,16,800,7,128,13,101,303)['F']==leading_global_rows(Q(1,200),1,16,1600,7,128,13,101,303)['F']
    checks['axial_rows_keep_whole_support_primitive_inputs']=o['rows']['U_original_support']==16
    checks['same_implemented_Borel_K2_is_used']=o['cutoff_K2']==4928
    checks['full_phase_budget_remains_positive_and_conditional']=o['Cphi']>0 and not o['infinite_Borel_schedule_and_actual_selector_verified']
    oldL=radial.reconstruct(f['inputs']['h'],0,2,c['native_Xhi'],c['native_Rlo'],c['native_Rhi'],o['rows']['F'],o['rows']['U'])
    oldF=radial.reconstruct(f['inputs']['h'],1,1,c['native_Xhi'],c['native_Rlo'],c['native_Rhi'],o['first_state']['jets']['F'][:2],o['first_state']['jets']['U'][:2])
    dep=borel.comparison_dependence(f['inputs']['h'],oldL['through'][1],oldL['through'][2],oldF['through'][1],g)
    checks['support_aware_norms_improve_manufactured_Cphi_without_shortening_pressure']=o['Cphi']<borel.comparison_value(dep,4928)
    return dict(checks=checks,Cphi=str(o['Cphi']),primitive_volume=str(o['leading']['primitive_volume']),norm_hull_Xhi=str(c['native_Xhi']),
                scope='same formulas and original source support with derived radius enclosures; no new gap floor or accepted band inferred')


@lru_cache(None)
def review():
    parts=dict(geometry=geometry_review(),cover=cover_review(),support=support_review(),join=join_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
