"""Exact representative, physical norm hull and floor-only dependency checks."""
from fractions import Fraction as Q
from functools import lru_cache
import json
from budgets import HERE,verify_representative,norm_cover,floor_join,native,tangent


def nearest(x):return (x+Q(1,2)).numerator//(x+Q(1,2)).denominator


def representative_review():
    h=Q(1,200);ell=16;d=Q(1,ell**6);checks={};rows=[]
    for R in (Q(1),Q(2)):
        for eta in (Q(-1),Q(0),Q(3,5),Q(1)):
            point=(R,eta,1-eta*eta);index=tuple(nearest(v/d) for v in point)
            r=verify_representative(h,ell,Q(1,2),2,index,point,1)
            name=f'R{R}_eta{eta}'
            checks[name+'_original_physical_active_clauses_and_box_pass']=r['all_physical_native_clauses_pass']
            checks[name+'_representative_not_moved']=r['existing_point']==point and r['existing_index']==index and not r['representative_selected_or_changed']
            rows.append(dict(R=str(R),eta=str(eta),native_edge=r['X'] in (Q(1,2),Q(2)),physical_boundary=point[2]==0))
    for t in (Q(99999,100000),Q(100001,100000)):
        s=t**400;eta=Q(3,5);point=(Q(3,2),eta*t**198,s*(1-eta*eta));index=tuple(nearest(v/d) for v in point)
        r=verify_representative(h,ell,Q(1,2),2,index,point,t)
        checks[f'root{t}_nonconstant_s_exact_physical_chart']=r['all_physical_native_clauses_pass'] and r['eta']==eta and r['s']==s
    index=(ell**6,ell**6,-1);point=(Q(1),Q(1),Q(0));rep=verify_representative(h,ell,Q(1,2),2,index,point,1)
    center=(Q(1),Q(1),-d);wrong=verify_representative(h,ell,Q(1,2),2,index,center,1)
    checks['selected_boundary_box_can_have_negative_time_center']=rep['all_physical_native_clauses_pass'] and center[2]<0
    checks['negative_center_is_not_an_original_active_representative']=wrong['gates']['in_original_closed_box'] and not wrong['all_physical_native_clauses_pass'] and not wrong['gates']['physical_T_nonnegative']
    outside=verify_representative(h,ell,1,2,(ell**6,0,ell**6),(1,0,1),1)
    checks['arbitrary_box_point_outside_native_annulus_is_rejected']=outside['gates']['in_original_closed_box'] and not outside['gates']['original_native_active_shell']
    return dict(checks=checks,rows=rows,scope='manufactured supplied rational chart data; original universal membership follows equation (6.8), not these finite witnesses')


def cover_review():
    checks={}
    for ell in (16,32,10000):
        c=norm_cover(ell,Q(1,2),2,Q(1,2),2)
        old=native.cover_existing(ell,Q(1,2),2,Q(1,2),2,Q(1,1000000000000000))
        keys=('delta','band_Rmin','band_Rmax','native_Rlo','native_Rhi','native_Xhi')
        checks[f'ell{ell}_same_original_norm_hull']=all(c[k]==old[k] for k in keys) and all(c['gates'].values())
        checks[f'ell{ell}_no_value_neighborhood_or_negative_time_certificate']=not c['value_neighborhood_selected_or_required_for_this_norm_cover'] and not c['negative_T_extension_certified'] and not c['growth_or_cone_floors_extended_to_norm_hull']
    d=Q(1,16**6);bad=norm_cover(16,d*d,4,d,2)
    checks['failed_original_axis_clearance_is_reported_without_reselection']=not bad['gates']['radial_axis_clearance'] and bad['existing_ell']==16
    # Exact-real sqrt(Xa) clearance: with an exact rational square it is
    # literally the same old geometric gate; no smaller rational root chosen.
    c=norm_cover(16,36*d*d,4,6*d,2)
    checks['original_geometric_clearance_matches_exact_square_root_boundary']=c['gates']['radial_axis_clearance'] and c['rA']**2==c['Xa'] and 6*d==c['rA']
    return dict(checks=checks)


def geometry():return dict(U=2,lambda_max=3,gmin=1,Llow=1,Lhigh=2,Dbox=6)


def join_review():
    checks={};rows=[]
    for ell in (16,10**12):
        c=norm_cover(ell,1,4,1,2);r=floor_join(c,3,1,geometry())
        old=native.cover_existing(ell,1,4,1,2,Q(1,ell**6)*42/2)
        direct=tangent.comparison(3,c['band_Rmin'],c['band_Rmax'],**geometry())
        checks[f'ell{ell}_same_complete_tangential_recipe']=r['comparison']==direct
        checks[f'ell{ell}_arbitrarily_narrow_unrelated_mu_can_fail_without_blocking_floor_recipe']=not old['gates']['existing_leading_native_neighborhood'] and all(c['gates'].values()) and not r['independent_native_value_neighborhood_needed_for_this_floor']
        checks[f'ell{ell}_original_band_and_all_scope_limits_retained']=r['existing_S']==ell**2 and not r['scales_supports_or_representatives_changed'] and not r['actual_existing_band_floor_certified'] and not r['negative_T_representative_or_extension_needed_for_this_floor']
        rows.append(dict(existing_ell=str(ell),Ctan=str(direct['Ctan']),normal_margin=str(r['normal_margin']),supplied_margin_passes=r['supplied_normal_margin_passes']))
    checks['supplied_old_small_scale_margin_failure_is_reported']=not floor_join(norm_cover(16,1,4,1,2),3,1,geometry())['supplied_normal_margin_passes']
    checks['supplied_large_manufactured_scale_passes_without_proving_actual_band']=floor_join(norm_cover(10**12,1,4,1,2),3,1,geometry())['supplied_normal_margin_passes']
    return dict(checks=checks,rows=rows,scope='parameter-space sufficient-margin controls only; neither existing actual band nor profile inputs are evaluated')


def source_review():
    c=json.loads((HERE/'source_contract.json').read_text())
    old=json.loads((HERE.parent/'background_jet_foundation/protocol.json').read_text())
    checks={
        'source_defines_intersection_representative_not_arbitrary_box_point':c['representative_definition']['original_representatives']=='x0_(ell,a) belongs to B_(ell,a) intersect A_ell',
        'source_keeps_physical_native_and_original_q_domain_clauses':c['representative_definition']['active_set_clauses']==['0<q<qbig','1/2<=q/Q<=2','tau>=0','Xa<=r^2/(2q)<=Xb'],
        'source_physical_relative_value_neighborhood_and_one_sided_boundary':c['value_neighborhood_convention']['relative_physical_domain']=='tau>=0' and c['value_neighborhood_convention']['boundary_derivatives']=='one-sided at tau=0',
        'historical_pdf_digest_metadata_matches_without_claiming_recomputed_binary_hash':c['historically_recorded_pdf_sha256']==old['paper_sha256'] and not c['pdf_binary_hash_recomputed_this_checkpoint'],
        'manuscript_theorems_not_assumed_by_using_its_definitions':c['source_definition_is_not_acceptance_of_manuscript_theorems'],
        'global_smooth_extension_remains_separate':'not proved here' in c['separate_extension_obligation'],
    }
    return dict(checks=checks,scope='source-contract integrity and historical metadata checks; not a binary PDF hash recomputation or a proof of manuscript theorems')


@lru_cache(None)
def review():
    parts=dict(representatives=representative_review(),cover=cover_review(),join=join_review(),source=source_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
