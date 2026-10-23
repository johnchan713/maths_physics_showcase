"""Current mean/scalar/radial allocation, independent physical budgets and jets."""
from fractions import Fraction as Q
from functools import lru_cache
import sys
import sympy as s
from budgets import HERE,load,wave,current,mean,physical,P,rates,postwave_rows_envelope,postfast_scalars,current_old_Haar,final_tangential_envelope,full_regular_source,request_trace,compact_physical_orders,regular_physical,SourcedMeanTailJoin


def fixture():
    class Common:reweighting=dict(original_pole='1',L='4',width='1/2',gap='1')
    common=Common()
    class Final:means=common
    class Haar:means=common
    class Signed:
        final=Final();haar=Haar()
        def target_inputs(self,n):return dict(stage='same-current-before-generic-signed')
    m=mean.PostSignedMeanRecipes(Signed(),wave_pole=30,original_pole=1,geometry=common.reweighting)
    return SourcedMeanTailJoin(wave.NextSourcedRecipes(current.CurrentResidualJoin(m,Q(1,200))))


def import_control(name,path,budget):
    before=sys.modules.get('budgets');sys.modules['budgets']=budget
    try:return load(name,path)
    finally:
        if before is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=before


@lru_cache(None)
def ledger_review():
    Ht,Hz,Ft,Fz,A1,rf0,rf1,rs0,rs1,rfc,rsc=s.symbols('Ht Hz Ft Fz A1 rf0 rf1 rs0 rs1 rfc rsc')
    old=s.Matrix([Ht+Ft,Hz+Fz]);principal=s.Matrix([-Ft,-Fz+A1]);rf=s.Matrix([rf0+rfc,rf1]);rs=s.Matrix([rs0+rsc,rs1])
    exact=old+principal+rf+rs;regular=s.Matrix([Ht,Hz])+rf+rs;compact=s.Matrix([0,A1])
    rho,chi,Pold,Pfast,Pslow,A0old,A0change=s.symbols('rho chi Pold Pfast Pslow A0old A0change')
    Pnew=Pold+Pfast+Pslow
    # A0 gr has Haar piece chi*Pnew, cancelled by A0(rho Pnew).
    A0gr=chi*Pnew+A0old+A0change;A0comp=chi*Pnew
    radial=-rho*Pnew-(A0gr-A0comp)
    checks=dict(full_tangential_identity_keeps_entire_fast_and_slow_differences=exact==regular+compact,
        extra_Qcirc_on_differences_would_erase_Haar_terms=exact!=s.Matrix([Ht,Hz])+s.Matrix([rfc+rsc,0])+compact,
        exact_A1_residual_sign_is_positive=compact[1]==A1,
        compensated_A0_zero_Haar_piece_cancels=s.expand(radial+rho*Pnew+A0old+A0change)==0,
        current_radial_compensator_uses_new_pressure=s.diff(radial,Pfast)!=0 and s.diff(radial,Pslow)!=0,
        stale_prewave_A0_residual_is_detected=s.expand(radial-(-rho*Pnew-A0old))==-A0change,
        dropping_compensation_adds_spurious_Haar_A0_pressure=s.expand((-rho*Pnew-A0gr)-radial+chi*Pnew)==0,
        nonzero_radial_A0_does_not_follow_from_Haar_cancellation=s.diff(radial,A0change)==-1,
        complete_pressure_not_declared_scalar_small=rho*Pnew!=Pnew)
    # Whole moving pressure derivative and old/fast correlations survive.
    z=s.symbols('z');cr=1+z;pp=2+z*z;jz=1+z**3
    checks['whole_axial_pressure_derivative_keeps_moving_compensator']=s.expand(s.diff(jz+cr*pp,z)-s.diff(jz,z)-cr*s.diff(pp,z))==pp
    u=s.symbols('u',nonzero=True);oldv=2+u+1/u;newv=(u+1/u)/3
    constant=lambda f:s.expand(f).coeff(u,0)
    checks['FULL_old_fast_correlations_survive_Haar']=constant(oldv*newv)!=constant(oldv)*constant(newv)
    checks['mean_feedback_retains_actual_A1_velocity']=s.expand((1-A1)*z*s.diff(z**2,z)+z**2*s.diff((1-A1)*z,z))!=3*z*z
    return dict(checks=checks,scope='exact formal full-row and compensated-mode algebra, not selected T0 solutions; original compact/difference operators separately replayed')


@lru_cache(None)
def physical_review():
    fc=import_control('fresh_eight_final_physical',HERE.parent/'core_final_mean_row_closure/controls.py',mean.final)
    hc=import_control('fresh_eight_Haar_physical',HERE.parent/'core_preserved_haar_scalar_join/controls.py',mean.preserved)
    checks={};rows=[]
    old={k:P(i+2) for i,k in enumerate(mean.prefix.STATE)}
    fast={k:P(i+5) for i,k in enumerate(mean.prefix.STATE)};slow={k:P(i+8) for i,k in enumerate(mean.prefix.STATE)}
    fixed={k:P(2) for k in ('b','V','G','graph','c','r','rho','Kp','Kmom')}
    for n in (0,1,2):
        out=final_tangential_envelope(n,dict(Etheta=P(3),Ez=P(4)),P(5),old,fast,slow,fixed)
        reset=mean.preserved.haar_slow_reset(n,old,slow,P(2),P(3))
        feedback=full_regular_source(n,P(2),P(3),fast,slow,graph=P(2),reciprocal=P(2),normal=P(2))
        for t in (Q(1,2),Q(1,8)):
            ep=t**10;ev=lambda d:{k:v.evaluate(1) for k,v in d.items()}
            base=ev(fixed);base.update(b=2*ep,graph=2/t)
            o={k:v*t**(19 if k=='x' else 9) for k,v in ev(old).items()}
            f={k:v*t**(27 if k=='x' else 17) for k,v in ev(fast).items()}
            sl={k:v*t**(28 if k=='x' else 18) for k,v in ev(slow).items()}
            rf=fc.physical_difference(n,o,f,base,ep)
            rs=fc.physical_difference(n,{k:o[k]+f[k] for k in o},sl,base,ep)
            for key in ('Etheta','Ez'):
                haar=3 if key=='Etheta' else 4;tail=0 if key=='Etheta' else 5*t**19
                checks[f'n{n}_t{t}:full_current_tangential_C0_'+key]=haar*t**18+tail+rf[key]+rs[key]<=t**18*out[key].evaluate(1)
            checks[f'n{n}_t{t}:increment_pressure_not_promoted_to_full_row_order']=rf['pm']<=t**17*out['fast_difference']['pm'].evaluate(1)
            haarold={k:v*t**(20 if k=='x' else 10) for k,v in ev(old).items()}
            scal=hc.raw_scalar_budget(n,haarold,sl,2*ep,3,ep)
            checks[f'n{n}_t{t}:all_complete_postslow_scalars_C0_plus_one']=all(scal[k]<=t**28*reset[k].evaluate(1) for k in scal)
            d=n+4;Cmw=2*(2*d/t+ep*d+2)+16/t**5
            U=3*t**5;rawfeedback=Cmw*U*(sum(v.evaluate(1) for v in fast.values())*t**17+sum(v.evaluate(1) for v in slow.values())*t**18)
            checks[f'n{n}_t{t}:H4_mean_feedback_has_actual_E_gain']=2*t**17+rawfeedback<=t**17*feedback['F'].evaluate(1)
            # Full final gr includes old AND new fields and complete W.
            final={k:o[k]+f[k]+sl[k] for k in o};xx,uu,zz=(final[k] for k in mean.prefix.STATE)
            d=Q(n+4);dr=d*base['graph'];rr=base['r'];b=base['b'];W=3*t**10
            visc=d*d*(base['graph']**2+1)+rr*dr+rr*rr
            gr=d*(ep+base['c'])*xx+(dr+rr)*(2*b*xx+xx*xx+W)
            gr+=ep*d*(b*zz+base['G']*xx+xx*zz+W)+rr*(2*base['V']*uu+uu*uu+W)+ep*visc*xx
            norm={k:old[k].evaluate(1)+fast[k].evaluate(1)+slow[k].evaluate(1) for k in old}
            normalized=dict(ev(fixed),**norm,sigma=1,w=3)
            checks[f'n{n}_t{t}:full_final_gr_B_with_all_means_and_covariance']=gr<=t**9*mean.prefix.full_mean(n,normalized)['gr']
        oldS={k:P(i+3) for i,k in enumerate(('P','Jtheta','Jz'))}
        pf=postfast_scalars(n,oldS,P(5),P(6),old,fast,P(3))
        direct=postfast_scalars(n,{k:v.evaluate(2) for k,v in oldS.items()},Q(5),Q(6),{k:v.evaluate(2) for k,v in old.items()},{k:v.evaluate(2) for k,v in fast.items()},Q(3))
        checks[f'n{n}:postfast_three_channel_positive_polynomial_substitution']=all(v.evaluate(2)==direct[k] for k,v in pf.items())
        H=current_old_Haar(old,fast,slow)
        checks[f'n{n}:preserved_Haar_is_all_three_historical_slow_summands']=all(H[k].c==(old[k]+fast[k]+slow[k]).c for k in old)
        covfixed={k:P(2) for k in ('graph','r','Kmom','Kp','rho')}
        rows=postwave_rows_envelope(n,dict(Etheta=P(3),Ez=P(4)),P(5),covfixed)
        change=mean.fast.covariance_rows(n,P(5),covfixed)
        checks[f'n{n}:postwave_Qcirc_rows_derived_from_old_C0_and_complete_covariance']=all(rows[k].c==(2*((P(3) if k=='Etheta' else P(4))+change[k])).c for k in rows)
        weak=dict(covfixed,Kp=P())
        checks[f'n{n}:postwave_axial_row_keeps_full_reconstructed_pressure']=rows['Ez'].evaluate(1)>postwave_rows_envelope(n,dict(Etheta=P(3),Ez=P(4)),P(5),weak)['Ez'].evaluate(1)
        calls=[]
        def row(j):calls.append(j);return dict(Etheta=P(j+1),Ez=P(j+2))
        bound=mean.fast_envelope(n,row,P(2),P(3),P(4))
        checks[f'n{n}:actual_fast_recipe_reads_derived_n4_n5_and_keeps_A1']=calls==[n+4,n+5] and bound['z'].evaluate(1)>bound['u'].evaluate(1)
    return dict(checks=checks,scope='independent unreduced rational physical epsilon budgets at kappa=1/10 and exact norm sums; not actual pulse or five-map enclosures')


@lru_cache(None)
def rate_review():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k);tag=str(k);C0=r['C0'];E=r['E']
        checks[tag+':derived_postfast_scalar_and_slow_orders_C0']=r['postfast_scalar']==r['slow']==C0 and r['B']+E>C0 and 2*E>C0
        checks[tag+':postslow_full_scalars_and_axial_pressure_moment']=r['postslow_scalar']==C0+1 and r['whole_axial_pressure_moment']==C0+2
        checks[tag+':full_fast_and_slow_remainders_strictly_above_C0']=r['tangential_terms']['fast_difference']-C0==1-3*k>0 and r['tangential_terms']['slow_difference']-C0==1-2*k>0
        checks[tag+':all_three_regular_mean_rows_C0']=r['final_tangential']==r['final_radial']==r['full_regular_mean']==C0
        checks[tag+':A1_two_integrations_but_not_one_for_strict_rule']=r['A1_integrations']==2 and E+k==C0 and E+2*k>C0
        checks[tag+':A0_strict_count_and_B_source_derived']=r['A0']>C0 and r['B']+(r['A0_integrations']-1)*k<=C0 and r['radial_source']==r['B']
        checks[tag+':whole_regular_E_not_mean_C0']=r['whole_regular']==r['nonzero_fast_feedback']==E<C0==r['nonzero_slow_feedback']
        rows.append({name:str(v) for name,v in r.items() if type(v) in (Q,int)})
    r=rates();checks['normalized_raw_A1_7_and_A0_200003_requests']=r['A1_source_order_zero']==7 and r['A0_source_order_zero']==200003
    return dict(checks=checks,rows=rows,scope='exact sufficient affine order inequalities; original weighted/support inputs not evaluated')


@lru_cache(None)
def recipe_review():
    j=fixture();c=j.certificate(0);nd=c['normalized_domains'];pd=c['physical_domains'];checks={}
    checks['normalized_trace_reads_current_old_source_A0_J_plus14']=nd['trace']['maxima']['old_current_regular_source']==200017 and nd['trace']['maxima']['original_inverse_normal']==200018
    checks['normalized_master_q_traces_older_graph_not_low_output']=nd['older_source_graph_order']==200017 and nd['older_padded_trace']['maxima']['next_signed_target_input']==400034 and nd['common_pole']['q']=='2400216'
    checks['physical_current_A0_and_A1_requests_use_current_E_source']=c['physical_compact']['A0']['source_derivatives']==60200007 and c['physical_compact']['A1']['source_derivatives']==60000011
    checks['physical_trace_reaches_current_source_and_normal_jets']=pd['trace']['maxima']['old_current_regular_source']==60200021 and pd['trace']['maxima']['original_inverse_normal']==60200022
    checks['physical_q_covers_older_target60400038']=pd['older_padded_trace']['maxima']['next_signed_target_input']==60400038 and pd['common_pole']['q']=='362400240'
    for name,d in (('normalized',nd),('physical',pd)):
        q=Q(d['common_pole']['q']);l=d['ledger']
        checks[name+':original_6q_12q_domains_cover_actual_fast_and_final_gr']=l['new_fast']==str(6*q) and l['final_full_gr']==str(12*q) and l['required_operators']['Kp']==[str(6*q),str(12*q)] and d['operator_recipe']['parameter']==str(q)
        checks[name+':full_nonzero_H4_feedback_source_pole9q']=l['full_nonzero_source']==str(9*q)
        checks[name+':all_archived_shifts_and_separate_Ie6q_retained']='K_old_all' in d['common_argument'] and d['oriented_target_recipe']['pole']==str(6*q)
        checks[name+':huge_source_graphs_not_expanded_evaluated_or_geometry_changed']=not d['older_padded_trace']['old_dependency_DAGs_expanded'] and not d['trace']['actual_input_jets_evaluated'] and not d['new_physical_band_roots_width_or_fields_changed'] and not d['operator_recipe']['materialized']
    for key,fac in (('A0',12),('A1',6)):
        q=Q(pd['common_pole']['q']);a=c['physical_compact'][key]
        checks[key+':current_full_source_raw_profiles_factorials_and_poles_retained']=a['raw_kernel']['source_pole']==str(fac*q) and a['raw_kernel']['source_factorial']==a['source_derivatives'] and not a['actual_input_or_factorial_evaluated']
    checks['current_scalar_recipe_uses_entire_Haar_slow_history']='previous generic-block C slow' in c['scalar']['current_Haar'] and 'FULL pre-source means' in c['scalar']['parts'][2]
    checks['whole_pressure_derivative_and_new_radial_compensation_kept']='derivative of c_rho' in c['scalar']['moving_pressure'] and 'gr_final-rho Pnew' in c['radial_identity']
    checks['exact_current_compact_tail_and_assigned_base_ledger']='assigned_once' in c['residual_identity'] and '+A1 Qcirc Epostwave_z' in c['residual_identity'] and '-A0 Qcirc gr_final' in c['residual_identity']
    checks['mean_updates_do_not_create_or_forward_Gaussian_tails']='append no new Gaussian tail' in c['Gaussian'] and 'never inverted' in c['Gaussian']
    checks['conditional_full_regular_gain_E_but_not_actual_complete_theorem']=c['whole_regular_source_order']==str(rates()['E']) and not c['actual_weighted_high_jet_support_base_validity_absolute_decay_or_blowup_certified']
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for m,N in ((0,0),(1,2),(3,7)):
            for e in (0,1):
                a=compact_physical_orders(e,m,N,h)
                checks[f'h{h}_m{m}_N{N}_e{e}:strict_physical_margin_and_count']=a['strict_margin']>0 and a['previous_integer_margin']<=0 and a['raw_source_order']==m+2*a['integration_order']+3
        p=regular_physical(0,h)
        checks[f'h{h}:regular_E_physical_exponent_negative_not_residual_lower_bound']=Q(p['native_physical_exponent'])==-Q(3,2)-3*h*Q(1,100000)<0 and not p['absolute_decay_certified'] and not p['actual_residual_lower_bound']
    try:j.certify_actual_inputs()
    except NotImplementedError:checks['actual_all_scale_input_validity_not_certified_by_recipes']=True
    else:checks['actual_all_scale_input_validity_not_certified_by_recipes']=False
    return dict(checks=checks,certificate=c,scope='finite normalized/physical request graph and deferred raw operator recipes; illustrative h=1/200,p=1 are not selected actual inputs')


@lru_cache(None)
def review():
    parts=dict(ledger=ledger_review(),physical=physical_review(),rates=rate_review(),recipes=recipe_review())
    return dict(parts=parts,checks={name+':'+key:bool(v) for name,p in parts.items() for key,v in p['checks'].items()})
