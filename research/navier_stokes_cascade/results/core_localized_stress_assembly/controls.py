"""Independent exact finite localization diagnostics, not pulse enclosures."""
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
from assembly import principal_defect,rates,haar,LocalizedStressJoin,prefix


def zero(x):
    if isinstance(x,s.MatrixBase):return all(s.cancel(v)==0 for v in x)
    return s.cancel(x)==0


@lru_cache(None)
def assembly_review():
    R,Z,q,v,x,theta=s.symbols('R Z q v x theta',positive=True)
    # Genuine finite integrals of two cutoffs and angular half. This
    # polynomial bump fixture is not the original flat C-infinity pulse.
    psi=v*(1-v);chi=x*(1-x);J=s.Rational(7,5)
    temporal=s.integrate(psi**2,(v,0,1));transverse=s.integrate(chi**2,(x,0,1))
    half=s.integrate(s.cos(theta)**2,(theta,0,2*s.pi))/(2*s.pi)
    mass=J*half*temporal*transverse
    vectors=[s.Matrix([1,2+Z,1+R]),s.Matrix([1,2+Z,-1-R])]
    H=s.Matrix([[mass*t[0]*t[e] for t in vectors] for e in (1,2)])
    T0=s.Matrix([2+Z,(1+R)/3]);Sigma=T0/q**3
    # h=2, A=h+1/2: different bands, common physical q. A rational
    # squared partition gives independent exact derivatives on a patch.
    z=R*Z/(1+R*Z);etas=[(1-z*z)/(1+z*z),2*z/(1+z*z)]
    labels=[];tensors=[]
    for Qb,eta in zip((s.Rational(1,4),s.Rational(1,8)),etas):
        eps=Qb**2;T=(Qb/q)**3*T0;y=H.inv()*T;a=Qb**-5*eps
        labels.append(dict(eta=eta,H=H,y=y,T=T,a=a))
        tensors.append(eta**2*a*mass*sum((y[i]*(t*t.T) for i,t in enumerate(vectors)),s.zeros(3)))
    r=principal_defect(Sigma,labels);W=sum(tensors,s.zeros(3))
    selected=s.Matrix([W[0,1],W[0,2]])
    checks=dict(
        two_cutoffs_angular_half_and_Haar_J_retained=temporal==transverse==s.Rational(1,30) and half==s.Rational(1,2) and mass==s.Rational(7,9000),
        positive_roots_on_manufactured_patch=all(all(y.subs({R:2,Z:s.Rational(1,3),q:s.Rational(1,6)})>0 for y in lab['y']) for lab in labels),
        cutoff_columns_match_direct_tensor_integrals=zero(selected-r['actual']),
        physical_scale_normalization_exact=all(zero(lab['a']*lab['T']-Sigma) for lab in labels),
        squared_partition_counted_once_per_box=zero(r['cover']-1),
        full_assembled_selected_stress_exact=zero(r['defect']) and zero(selected-Sigma),
        unconditional_three_defect_identity=zero(r['defect']-r['root']-r['normalization']-r['coverage']),
        all_six_covariance_entries_retained=all(not zero(W[i,j]) for i in range(3) for j in range(i,3)),
        unmatched_cutoff_not_needed_for_matched_H=not zero(psi**2-1) and zero(r['defect']),
    )
    for name,der in [('R',lambda a:a.diff(R)),('Z',lambda a:a.diff(Z)),('q',lambda a:a.diff(q)),('RZ',lambda a:a.diff(R,Z))]:
        checks['entire_partition_identity_survives_'+name+'_derivatives']=zero(der(r['defect']))
    checks['partition_derivatives_nonzero_before_exact_sum']=not zero(s.diff(etas[0]**2,R)*Sigma)
    # Negative controls change the actual construction, not just a flag.
    variants={
        'stale_roots_after_recutted_pulse':lambda lab:dict(lab,H=lab['H']*s.Rational(3,4)),
        'missing_temporal_cutoff_in_column':lambda lab:dict(lab,y=lab['y']*temporal),
        'missing_transverse_cutoff_in_column':lambda lab:dict(lab,y=lab['y']*transverse),
        'missing_angular_half_in_column':lambda lab:dict(lab,y=lab['y']*half),
        'missing_Haar_J_in_column':lambda lab:dict(lab,y=lab['y']*J),
        'wrong_physical_velocity_square_scale':lambda lab:dict(lab,a=lab['a']*2),
    }
    for name,change in variants.items():
        bad=principal_defect(Sigma,[change(lab) for lab in labels])
        checks['detect_'+name]=not zero(bad['defect'])
        checks['defect_identity_for_'+name]=zero(bad['defect']-bad['root']-bad['normalization']-bad['coverage'])
    checks['detect_signs_counted_twice_as_boxes']=not zero(principal_defect(Sigma,labels+labels)['defect'])
    checks['detect_incomplete_active_cover']=not zero(principal_defect(Sigma,labels[:1])['defect'])
    # A perturbation can occupy all three terms of the identity at once.
    changed=[dict(labels[0],H=H*2,a=labels[0]['a']*3,eta=labels[0]['eta']/2)]
    bad=principal_defect(Sigma,changed)
    checks['simultaneous_three_nonzero_defects_retained']=all(not zero(bad[k]) for k in ('root','normalization','coverage')) and zero(bad['defect']-bad['root']-bad['normalization']-bad['coverage'])
    return dict(checks=checks,mass=str(mass),scales=['1/4','1/8'],h=2,A='5/2',
                scope='exact two-band local algebra with polynomial bump integrals and rational squared partition; not the selected flat continuum family or an accepted cone evaluation')


@lru_cache(None)
def budget_review():
    # No expensive actual giant-jet construction is materialized.
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k);tag=str(k)
        checks[tag+':full_tangential_postwave_minimum_is_D']=r['postwave_complete_tangential']==r['D']<r['C']
        checks[tag+':final_Haar_and_radial_Haar_preserved']=r['final_Haar']==r['C'] and r['final_radial_Haar']==r['B']+r['C']
        checks[tag+':no_full_radial_tail_or_actual_input_promotion']=not any(r[a] for a in ('full_pointwise_radial_or_tail_bound','independent_actual_input_certification','complete_physical_residual_or_blowup'))
        rows.append({a:str(v) for a,v in r.items()})
    # Trace the stage using an independent stub. Its integer polynomials
    # are finite composition diagnostics, never actual input enclosures.
    P=prefix.P;calls=[]
    class Means:
        reweighting=dict(original_pole='1',L='2',width='1',gap='1')
        def data(self,kind,n):
            calls.append((kind,n))
            if kind=='initial':return {a:P(i+2) for i,a in enumerate(prefix.STATE)}
            if kind=='wave':return dict(U=P(7))
            if kind=='fixed':return {a:P(1) for a in ('b','V','G','r','rho','Kp','Kmom')}
            raise AssertionError(kind)
        def final(self,n):raise AssertionError('Wrong mean stage')
    class Fast:
        def regular_rows(self,n):return dict(rows=dict(Etheta=P(11),Ez=P(13)))
        def postwave_recipe(self,n):
            return dict(exact_A1=dict(haar.fast.flat_recipe(n),original_geometry=dict(Means.reweighting)))
    class Haar:
        embedding=Q(2)
    j=object.__new__(LocalizedStressJoin);j.means=Means();j.fast=Fast();j.haar=Haar();j.kappa=Q(1,100000);j.orders=rates()
    out=j.postwave_envelope(0,principal_defect=P(17),fine_flux=P(19))
    checks['postwave_uses_initial_old_means_not_fresh_final']=('initial',3) in calls and out['stage']=='postwave-before-fresh-fast-and-slow'
    checks['leading_defect_retained_in_both_complete_rows']=all(out['Haar'][a].c for a in ('leading_theta','leading_z'))
    checks['combines_Q1_and_entire_Haar_rows']=all(out['combined_regular'][a].c==(out['Q1_regular']['rows'][a]+out['Haar'][b]+out['Haar'][c]).c for a,b,c in (('Etheta','fine_theta','leading_theta'),('Ez','fine_z','leading_z')))
    checks['full_pressure_source_and_compensator_retained']=all(out['Haar'][a].c for a in ('gr','P','pm'))
    recipe=out['exact_A1']
    core_recipe={k:v for k,v in recipe.items() if k!='original_geometry'}
    checks['exact_original_A1_high_jet_recipe_retained']=core_recipe==haar.fast.flat_recipe(0) and recipe['raw_source_order']==200001
    checks['exact_A1_original_pole_width_gap_geometry_retained']=recipe['original_geometry']==Means.reweighting
    checks['additive_tails_excluded_from_regular_join']='excluded' in out['additive_tail_ledger']
    try:j.actual_input_certificate()
    except NotImplementedError:checks['no_actual_input_certificate_from_recipe']=True
    else:checks['no_actual_input_certificate_from_recipe']=False
    return dict(checks=checks,rates=rows,exact_A1_source_order=recipe['raw_source_order'])


@lru_cache(None)
def review():
    parts=dict(assembly=assembly_review(),budget=budget_review())
    return dict(parts=parts,checks={name+':'+k:v for name,part in parts.items() for k,v in part['checks'].items()})
