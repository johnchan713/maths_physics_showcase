"""Exact complete cylindrical Fourier residual controls, not actual pulse solutions."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys

@lru_cache(None)
def algebra():
    p=Path(__file__).resolve().parent.parent/'finite_induction_foundation/finite_checks.py'
    spec=importlib.util.spec_from_file_location('current_residual_finite_algebra',p)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m

def flat_certificate(D,N,b,c):
    if any(isinstance(x,(bool,float)) for x in (D,N,b,c)):raise TypeError('Exact rational controls only')
    D,N,b,c=map(F,(D,N,b,c))
    if min(D,N,b)<0 or c<=0:raise ValueError('Nonnegative losses/power and positive Gaussian constant required')
    x=b/2;n=-(-x.numerator//x.denominator)
    return dict(D=str(D),N=str(N),b=str(b),c=str(c),series_order=n,
        prefactor=str(F(factorial(n))*(2/c)**n),quadratic_log_constant=str((D+N)**2/(2*c)))

def cycle_certificate():
    k=F(1,100000);B0=F(3,2)-2*k;gain=F(1,2)-k
    rows={}
    for step in ('sourced','signed'):
        for term in ('complete_linear','old_wave_both_crosses','actual_mean_both_crosses'):
            rows[step+'_'+term]=(1,gain)
        rows[step+'_complete_square']=(2,-k)
    rows.update(actual_fast_mean_wave=(1,gain),actual_slow_mean_wave=(1,F(1,2)))
    out={n:dict(order_slope=str(a),order_offset=str(b),margin_slope=str(a-1),boundary_margin=str((a-1)*B0+b-gain)) for n,(a,b) in rows.items()}
    return dict(B0=str(B0),gain=str(gain),rows=out,checks=dict(
        whole_half_line_slopes_nonnegative=all(F(r['margin_slope'])>=0 for r in out.values()),
        all_boundary_margins_nonnegative=all(F(r['boundary_margin'])>=0 for r in out.values()),
        exactly_seven_nonzero_groups_tie=sum(r['margin_slope']=='0' and r['boundary_margin']=='0' for r in out.values())==7))

def gaussian_review():
    rows=[]
    for D,N,b,c in ((0,0,0,1),(7,11,3,F(1,8)),(19,23,F(7,2),F(1,10**12)),(100,200,12,F(1,3))):
        r=flat_certificate(D,N,b,c);cc=F(r['c']);d=F(r['D'])+F(r['N'])
        margins=[]
        for ell in (F(1),F(7,3),F(10**3),F(10**9)):
            lhs=-cc*ell**2/2+d*ell
            gap=F(r['quadratic_log_constant'])-lhs
            margins.append(gap==cc*(ell-d/cc)**2/2 and gap>=0)
        r['exact_completed_square_controls']=all(margins);rows.append(r)
    h=F(1,200);delta=F(1,2)-F(1,100000);obstruction=[]
    for j in (8,16,24,32):
        obstruction.append(dict(ell_power_of_two=j,log2_normalized_ratio=str(h*delta*2**j-4*j)))
    spatial=dict(h='1/200',pressure_coefficient_order='1',requested_flat_power='1',
        log_growth_slope=str(1-F(1,200)),nonzero_spatial_derivative_at_temporal_peak=True)
    return dict(rows=rows,normalized_source_obstruction=obstruction,spatial_peak_not_Gaussian_tail=spatial,checks=dict(
        every_exact_square_completion=all(r['exact_completed_square_controls'] for r in rows),
        series_order_dominates_each_polynomial=all(2*r['series_order']>=F(r['b']) for r in rows),
        Gaussian_smallness_is_not_normalized_order_gain=all(F(r['log2_normalized_ratio'])>0 for r in obstruction[1:]),
        generic_spatial_peak_is_not_flat=F(spatial['log_growth_slope'])>0,
        continuum_flatness_argument_not_inferred_from_parameter_scan=True))

@lru_cache(None)
def exact_review():
    f=algebra();a=f.algebra();e=a.exact();R,Z,T,ep=e.R,e.Z,e.T,e.EPS
    z=e.Poly();zero=(z,z);zvec=lambda:[zero for _ in range(3)]
    cp=lambda p:(p,z)
    def cvadd(*vs):return [a.sumc([v[i] for v in vs],z) for i in range(3)]
    def vscale(v,c):return [a.scale(x,c) for x in v]
    def vneg(v):return [a.neg(x) for x in v]
    def diff(x,j):return (x[0].diff(j),x[1].diff(j))
    def shift(x,k):return (x[0].shift(k),x[1].shift(k))
    def nonzero(v):return any(bool(r) or bool(i) for r,i in v)
    def equal(v,w):return all(x==y for x,y in zip(v,w))
    def clean(d):return {m:v for m,v in d.items() if nonzero(v)}
    def join(*ds):
        keys=set().union(*(set(d) for d in ds))
        return clean({m:cvadd(*(d.get(m,zvec()) for d in ds)) for m in keys})
    def negative(d):return {m:vneg(v) for m,v in d.items()}
    def scaled(d,c):return clean({m:vscale(v,c) for m,v in d.items()})
    def perp(d):return {m:v for m,v in d.items() if m}
    def equalfields(d,g):return all(equal(d.get(m,zvec()),g.get(m,zvec())) for m in set(d)|set(g))
    def gradient(p,m):return [diff(p,0),a.scale(a.imag(shift(p,-1)),3*m),a.add(a.scale(diff(p,1),ep),a.scale(a.imag(p),3*m))]
    def visc(v,m):
        out=[];w=3*m
        for x in v:
            out.append(a.sumc([diff(diff(x,0),0),shift(diff(x,0),-1),a.scale(diff(diff(x,1),1),ep**2),
                a.scale(a.add(shift(x,-2),x),-w*w),a.scale(a.imag(diff(x,1)),2*w*ep)],z))
        out[0]=a.add(a.add(out[0],a.neg(shift(v[0],-2))),a.scale(a.imag(shift(v[1],-2)),-2*w))
        out[1]=a.add(a.add(out[1],a.neg(shift(v[1],-2))),a.scale(a.imag(shift(v[0],-2)),2*w))
        return out
    def linear(v,p,m,base):
        return cvadd(vscale([diff(x,2) for x in v],-ep),vscale(visc(v,m),-ep),
            f.transport(base,v,m),f.transport(v,base,0),gradient(p,m))
    def curl(C,m):
        base=[a.scale(diff(C[1],1),-ep),a.add(a.scale(diff(C[0],1),ep),a.neg(diff(C[2],0))),a.add(diff(C[1],0),shift(C[1],-1))]
        cross=[a.add(shift(C[2],-1),a.neg(C[1])),C[0],a.neg(shift(C[0],-1))]
        return cvadd(base,[a.scale(a.imag(x),3*m) for x in cross])
    def divergence(v,m):return a.sumc([diff(v[0],0),shift(v[0],-1),a.scale(diff(v[2],1),ep),a.scale(a.imag(a.add(shift(v[1],-1),v[2])),3*m)],z)
    psiB=R**2*(1+Z/11+T/13)
    base=[cp(-ep*psiB.diff(1)),cp(R*(2+Z+T)),cp(e.radial(psiB))]
    nonlinear_cache={}
    def residual(fields,press):
        # Pressure changes do not alter the already exact velocity convolution.
        key=id(fields)
        if key not in nonlinear_cache:nonlinear_cache[key]=f.convolution(fields,fields)
        out=dict(nonlinear_cache[key])
        for m in set(fields)|set(press):
            v=fields.get(m,zvec());p=press.get(m,zero)
            extra=cvadd(vscale([diff(x,2) for x in v],-ep),vscale(visc(v,m),-ep),gradient(p,m))
            out[m]=cvadd(out.get(m,zvec()),extra)
        return clean(out)
    waves,mean,pots=f.fixture();bump=(R-1)**3*(2-R)**3
    oldpsi=(R-1)**2*(2-R)**2*(1+R/7+Z/11+T/13)
    psi=bump*(1+Z/7+T/11)
    def cutoff(C,p,m,psi):
        u=curl(C,m);cut=curl(vscale(C,psi),m)
        gR,gZ=psi.diff(0),ep*psi.diff(1)
        d=[a.scale(C[1],-gZ),a.add(a.scale(C[0],gZ),a.scale(C[2],-gR)),a.scale(C[1],gR)]
        temporal=vscale(u,-ep*psi.diff(2)+base[0][0]*gR+base[2][0]*gZ)
        mixed=[]
        for x in u:
            dz=a.add(a.scale(diff(x,1),ep),a.scale(a.imag(x),3*m))
            mixed.append(a.scale(a.add(a.scale(diff(x,0),gR),a.scale(dz,gZ)),-2*ep))
        dlap=psi.diff(0).diff(0)+psi.diff(0).shift(-1)+ep**2*psi.diff(1).diff(1)
        lapcut=vscale(u,-ep*dlap);pc=[a.scale(p,gR),zero,a.scale(p,gZ)]
        ld=linear(d,zero,m,base);tail=cvadd(temporal,mixed,lapcut,pc,ld)
        actual=cvadd(linear(cut,a.scale(p,psi),m,base),vneg(vscale(linear(u,p,m,base),psi)))
        return dict(u=u,cut=cut,d=d,tail=tail,actual=actual,groups=dict(temporal=temporal,mixed_viscosity=mixed,
            cutoff_laplacian=lapcut,pressure=pc,cutoff_curl_operator=ld))
    oldp=(R**2*bump*(1+Z/3+T/5),R**3*bump*(1+Z*T/7))
    oldc=cutoff(pots[1],oldp,1,oldpsi)
    oldw={1:oldc['cut'],-1:[a.conj(x) for x in oldc['cut']]}
    oldT={1:oldc['tail'],-1:[a.conj(x) for x in oldc['tail']]}
    pB=R**2+Z**2+ep**2*(Z*T+Z**3)-ep**2*R.shift(-3)/3+ep**4*(Z+T*Z**2)
    pm=R*bump*(1+Z*T);pmnew=pm+ep*R**2*bump*(Z+T)
    pold={0:cp(pB+pm),1:a.scale(oldp,oldpsi),-1:a.conj(a.scale(oldp,oldpsi))}
    oldfields={0:cvadd(base,mean),**oldw};before=residual(oldfields,pold)
    source=join(perp(before),negative(oldT))
    delta={};dpressure={};tails={};remainder={};cutdata={}
    for m in (1,2):
        C=vscale(pots[m],F(1,m+7));p=(R**2*bump*(1+Z/(m+3)+T/(m+4)),R**3*bump*(1+Z*T)/(m+5))
        c=cutoff(C,p,m,psi);cutdata[m]=c
        for n,Cut,pii in ((m,c,a.scale(p,psi)),(-m,None,a.conj(a.scale(p,psi)))):
            delta[n]=c['cut'] if n>0 else [a.conj(x) for x in c['cut']]
            dpressure[n]=pii;tails[n]=c['tail'] if n>0 else [a.conj(x) for x in c['tail']]
            unc=linear(c['u'],p,m,base)
            remainder[n]=cvadd(unc,source.get(n,zvec())) if n>0 else cvadd([a.conj(x) for x in unc],source.get(n,zvec()))
    newfields=join(oldfields,delta);newp={m:a.add(pold.get(m,zero),dpressure.get(m,zero)) for m in set(pold)|set(dpressure)}
    newp[0]=cp(pB+pmnew);after=residual(newfields,newp)
    interacting_old={0:mean,**oldw}
    mixed=f.convolution(interacting_old,delta);reverse=f.convolution(delta,interacting_old);square=f.convolution(delta,delta)
    regular=join(scaled(remainder,psi),perp(mixed),perp(reverse),perp(square))
    ledger=join(oldT,scaled(source,1-psi),tails)
    ids=dict(base_divergence_free=divergence(base,0)==zero,
        old_complete_cutoff_commutator=equal(oldc['actual'],oldc['tail']),
        current_source_covers_each_inverted_mode=set(source)==set(delta),
        complete_current_sourced_residual_ledger=equalfields(perp(after),join(regular,ledger)),
        old_residual_split_is_exact=equalfields(perp(before),join(source,oldT)))
    for m,c in cutdata.items():
        ids[f'full_cutoff_commutator_mode_{m}']=equal(c['actual'],c['tail'])
        ids[f'cut_potential_equals_product_plus_curl_mode_{m}']=equal(c['cut'],cvadd(vscale(c['u'],psi),c['d']))
        ids[f'complete_cut_divergence_zero_mode_{m}']=divergence(c['cut'],m)==zero
        ids[f'cutoff_curl_divergence_cancels_mode_{m}']=a.add(divergence(vscale(c['u'],psi),m),divergence(c['d'],m))==zero
    controls={n:all(nonzero(c['groups'][n]) for c in cutdata.values()) for n in next(iter(cutdata.values()))['groups']}
    controls.update(naive_cut_velocity_is_not_divergence_free=all(divergence(vscale(c['u'],psi),m)!=zero for m,c in cutdata.items()),
        complete_cutoff_curl_not_zero=all(nonzero(c['d']) for c in cutdata.values()),
        both_old_new_cross_orientations_nonzero=bool(perp(mixed)) and bool(perp(reverse)),new_wave_square_nonzero=bool(perp(square)),
        old_tail_retained_nonzero=bool(oldT),uncancelled_source_cutoff_retained_nonzero=bool(scaled(source,1-psi)))
    # A changing mean pressure affects the FULL equation, but not its nonzero projection.
    removed_mean=dict(newp);removed_mean[0]=cp(pB+pm)
    no_mean_change=residual(newfields,removed_mean)
    ids['recomputed_mean_pressure_does_not_change_nonzero_residual']=equalfields(perp(after),perp(no_mean_change))
    ids['full_mean_pressure_difference_is_gradient']=equal(cvadd(after.get(0,zvec()),vneg(no_mean_change.get(0,zvec()))),gradient(cp(pmnew-pm),0))
    controls['actual_mean_pressure_gradient_change_nonzero']=nonzero(gradient(cp(pmnew-pm),0))
    pi2=ep**4*(Z+T*Z**2);without_tail=dict(newp);without_tail[0]=cp(pB+pmnew-pi2)
    trunc=residual(newfields,without_tail)
    ids['same_full_base_pressure_Pi2_cancels_only_in_difference']=equal(cvadd(after.get(0,zvec()),vneg(trunc.get(0,zvec()))),gradient(cp(pi2),0))
    controls['Pi2_axial_gradient_nonzero']=bool(pi2.diff(1))
    # Formula (5) for an actual curl mean increment, with all old/new modes.
    dm=vscale(mean,F(1,13));lastfields=join(newfields,{0:dm});lastp=dict(newp);lastp[0]=a.add(lastp[0],cp(ep*R*bump*(Z**2+T)))
    last=residual(lastfields,lastp);wtotal=join(oldw,delta)
    meanmix=join(f.convolution({0:dm},wtotal),f.convolution(wtotal,{0:dm}))
    ids['actual_mean_increment_divergence_zero']=divergence(dm,0)==zero
    ids['actual_mean_update_nonzero_residual_exact']=equalfields(perp(last),join(perp(after),perp(meanmix)))
    ids['mean_update_keeps_exact_nonzero_tail_ledger']=equalfields(perp(last),join(regular,perp(meanmix),ledger))
    controls['actual_mean_wave_change_nonzero']=bool(perp(meanmix))
    # Pressure-only wave isolates a term that a velocity-only cutoff would miss.
    pc=cutoff(zvec(),oldp,1,psi)
    ids['pressure_only_cutoff_commutator']=equal(pc['actual'],pc['groups']['pressure']) and equal(pc['tail'],pc['groups']['pressure'])
    controls['pressure_only_cutoff_term_nonzero']=nonzero(pc['tail'])
    # Gaussian-ledger cutoff is temporal; a general spatial-cutoff control
    # tests a different identity and cannot certify Gaussian flatness.
    timepsi=(T-1)**3*(2-T)**3;tc=cutoff(pots[1],oldp,1,timepsi)
    ids['temporal_cutoff_has_no_spatial_curl']=not nonzero(tc['d']) and equal(tc['cut'],vscale(tc['u'],timepsi))
    ids['temporal_cutoff_has_no_spatial_pressure_commutator']=not nonzero(tc['groups']['pressure'])
    ids['temporal_cutoff_tail_is_only_time_product']=equal(tc['tail'],vscale(tc['u'],-ep*timepsi.diff(2))) and equal(tc['actual'],tc['tail'])
    controls['temporal_time_product_nonzero']=nonzero(tc['tail'])
    ids['conjugate_complete_residual_reality']=all(equal(after.get(-m,zvec()),[a.conj(x) for x in v]) for m,v in after.items())
    ids['conjugate_exact_tail_reality']=all(equal(ledger.get(-m,zvec()),[a.conj(x) for x in v]) for m,v in ledger.items())
    controls['highest_generated_harmonic_four_present']=4 in after and nonzero(after[4])
    # Six full covariance entries of a cut wave, rather than psi^2 principal stress.
    c=cutdata[1]
    for i in range(3):
        for j in range(i,3):
            cov=lambda v,w:a.add(a.mul(v[i],a.conj(w[j])),a.mul(a.conj(v[i]),w[j]))
            rhs=a.sumc([a.scale(cov(c['u'],c['u']),psi**2),a.scale(a.add(cov(c['u'],c['d']),cov(c['d'],c['u'])),psi),cov(c['d'],c['d'])],z)
            ids[f'complete_cut_covariance_{i}{j}']=cov(c['cut'],c['cut'])==rhs
    return dict(identities=ids,nonzero_omission_controls=controls,source_harmonics=sorted(source),new_harmonics=sorted(perp(after)),
        fixture='exact graph-free normalized cylindrical Fourier polynomials; finite endpoint regularity; no inverse solve or actual giant fields')

def review():
    r=exact_review();g=gaussian_review();b=cycle_certificate()
    return dict(exact_complete_residual_ledger=r,Gaussian_physical_flatness_not_normalized_source_gain=g,
        all_half_line_regular_budget=b,checks=dict(all_complete_current_residual_identities=all(r['identities'].values()),
            all_omitted_groups_detected_nonzero=all(r['nonzero_omission_controls'].values()),
            all_fixed_jet_Gaussian_majorants=all(g['checks'].values()),whole_half_line_regular_cycle=all(b['checks'].values())))
