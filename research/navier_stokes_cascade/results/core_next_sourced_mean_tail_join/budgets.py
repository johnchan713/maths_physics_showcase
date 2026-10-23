"""Fresh means after the next eight-mode sourced wave, current radial repair."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

wave=load('next_mean_sourced_wave',HERE.parent/'core_next_sourced_wave_join/budgets.py')
current=wave.current;mean=current.mean;physical=current.physical
P,exact,order=mean.P,mean.exact,mean.order


def rates(kappa=Q(1,100000)):
    r=wave.rates(kappa);k=exact(kappa);B,C0,E=(r[x] for x in ('B','C0','E'))
    p1=mean.strict_integrations(E,C0,k);p0=mean.strict_integrations(B,C0,k)
    terms=dict(Haar_current=C0,fast_difference=E+1-2*k,slow_difference=C0+1-2*k,A1=E+p1*k)
    radial=dict(fast_time=B+1,wave_covariance_graph=1-k,
        mean_radial_flux=min(B+2-k,2*B+2-k),axial_flux=min(B+2,2*B+2,Q(2)),
        centripetal=min(B,2*B,Q(1)),radial_viscosity=B+2-2*k)
    return dict(B=B,C0=C0,E=E,postwave_Qcirc=E,fast=E,fast_radial=E+1,
        postfast_scalar=min(C0,B+E,2*E),slow=C0,slow_radial=C0+1,
        old_Haar_tangent=Q(1),old_Haar_radial=Q(2),postslow_scalar=min(C0+1,2*C0,C0+2),
        whole_axial_pressure_moment=C0+2,tangential_terms=terms,final_tangential=min(terms.values()),
        radial_source_groups=radial,radial_source=min(radial.values()),radial_compensator=C0+1,
        A0=B+p0*k,final_radial=C0,full_regular_mean=C0,
        A1_integrations=p1,A0_integrations=p0,A1_source_order_zero=2*p1+3,A0_source_order_zero=2*p0+3,
        nonzero_fast_feedback=E,nonzero_slow_feedback=C0,whole_regular=min(C0,E),
        actual_inputs_absolute_decay_or_infinite_field_certified=False)


def postwave_rows_envelope(n,old_rows,cov,fixed):
    order(n);d=mean.fast.covariance_rows(n,cov,fixed)
    return {k:2*(old_rows[k]+d[k]) for k in ('Etheta','Ez')}


def postfast_scalars(n,old_postslow,cov,cov_next,old_full,fast,M):
    """Fresh scalars from old C+1 + current full covariance + FULL old/fast."""
    return mean.postfast_scalar_envelope(n,old_postslow,cov,cov_next,old_full,fast,M)


def current_old_Haar(derived_old_Haar,previous_slow_C,previous_generic_slow_C):
    """PreservedHaarScalarJoin.old_Haar + latest C slow + generic C slow.

    Nonzero waves change no means; complete fast fields have zero Haar.
    This is not the full B norm relabelled as order one.
    """
    a=mean.old_Haar_envelope(derived_old_Haar,previous_slow_C)
    return mean.old_Haar_envelope(a,previous_generic_slow_C)


def final_tangential_envelope(n,Haar_postwave,A1,old_full,fast,slow,fixed):
    return mean.final_row_envelope(n,Haar_postwave,A1,old_full,fast,slow,fixed)


def full_regular_source(n,wave_stage,complete_current_wave,fast,slow,*,graph,reciprocal,normal):
    """H=4: the current wave has input harmonics through four, residual through8."""
    f=mean.mean_feedback(n,complete_current_wave,fast,graph=graph,reciprocal=reciprocal,normal=normal,max_harmonic=4)
    s=mean.mean_feedback(n,complete_current_wave,slow,graph=graph,reciprocal=reciprocal,normal=normal,max_harmonic=4)
    return dict(F=wave_stage+f+s,fast_feedback=f,slow_feedback=s)


def request_trace(n,p1,p0):
    """NEW layer only, then hand exact integer indices to the OLD source DAG.

    No giant polynomial, factorial, kernel or old norm program is evaluated.
    """
    order(n);order(p1);order(p0)
    if min(p1,p0)<1:raise ValueError('Positive raw compact integration counts required')
    calls={}
    def seen(name,j):calls.setdefault(name,set()).add(j)
    def generated(j):
        seen('current_complete_generated_wave',j);seen('old_current_regular_source',j+1);seen('original_inverse_normal',j+2)
    def covariance(j):seen('current_complete_covariance',j);generated(j)
    def row(j):seen('postwave_tangential_row',j);seen('old_pre_source_full_row',j);covariance(j+3)
    def fast(j):seen('new_actual_fast_field',j);row(j+4);row(j+5)
    def scalar(j):
        seen('new_postfast_scalar',j);seen('old_pre_source_postslow_scalar',j)
        covariance(j);covariance(j+1);fast(j+1)
    def slow(j):seen('new_actual_slow_field',j);scalar(j);scalar(j+1)
    def fields(j):seen('old_pre_source_mean_field',j);fast(j);slow(j)
    def gr(j):seen('current_final_full_gr',j);fields(j+3);generated(j+3)
    fields(n+3);generated(n+3)
    J1=n+2*p1+3;J0=n+2*p0+3;row(J1);gr(J0)
    maxima={k:max(v) for k,v in sorted(calls.items())}
    return dict(output_order=n,A1_source_order=J1,A0_source_order=J0,maxima=maxima,
        older_original_root_and_A1_DAGs_expanded=False,actual_input_jets_evaluated=False)


def compact_physical_orders(e,m,N,h,kappa=Q(1,100000)):
    order(m);N=exact(N);h=physical.scales(h,kappa)['h'];r=rates(kappa)
    if e not in (0,1) or N<0:raise ValueError('Original e=0/1 and nonnegative physical power required')
    a=r['B'] if e==0 else r['E'];loss=Q(3,2)+2*h+2*m;target=(N+loss)/h
    p=mean.strict_integrations(a,target,exact(kappa));exponent=h*(a+p*exact(kappa))-loss
    return dict(e=e,physical_jet=m,requested_q_power=N,source_rate=a,normalized_target=target,
        integration_order=p,raw_source_order=m+2*p+3,physical_exponent=exponent,strict_margin=exponent-N,
        previous_integer_margin=h*(a+(p-1)*exact(kappa))-loss-N)


def regular_physical(m,h,kappa=Q(1,100000),reserve=Q(1)):
    order(m);h=physical.scales(h,kappa)['h'];r=rates(kappa);reserve=exact(reserve)
    if reserve<=0:raise ValueError('Positive polynomial reserve required')
    loss=Q(3,2)+2*h+2*m;power=h*r['E']-loss
    return dict(normalized_regular_order=str(r['E']),native_physical_exponent=str(power),absorbed_exponent=str(power-reserve),
        sufficient_later_order=f'h*alpha > N + {loss} + {reserve}',absolute_decay_certified=power-reserve>0,
        actual_residual_lower_bound=False)


class SourcedMeanTailJoin:
    """Same derived next-wave object; no incoming desired row/scalar callbacks."""
    def __init__(self,sourced_recipe):
        if not isinstance(sourced_recipe,wave.NextSourcedRecipes):raise TypeError('Use the SAME next sourced-wave recipe')
        self.wave=sourced_recipe;self.current=sourced_recipe.current;self.before=self.current.mean
        self.kappa=self.before.kappa;self.geometry=self.before.geometry;self.h=self.current.h;self.orders=rates(self.kappa)

    def domain_recipe(self,n,*,physical_power=None):
        r=self.orders
        if physical_power is None:
            p1,p0=r['A1_integrations'],r['A0_integrations']
        else:
            p1=compact_physical_orders(1,n,physical_power,self.h,self.kappa)['integration_order']
            p0=compact_physical_orders(0,n,physical_power,self.h,self.kappa)['integration_order']
        trace=request_trace(n,p1,p0);j=trace['maxima']['old_current_regular_source']
        # Bound the SAME earlier source at its actually requested jet j.
        # The previous mean common-pole rule covers its coefficient calls
        # and retains the older normalized raw A0/A1 graph explicitly.
        c=self.before.common_pole(j);q=Q(c['q']);older=mean.dependency_trace(j,self.before.orders['A1_integrations'],self.before.orders['A0_integrations'])
        ledger=dict(old_current_means=str(2*q),complete_current_wave=str(3*q),postwave_rows=str(6*q),new_fast=str(6*q),
            postfast_full_sources=str(12*q),final_full_gr=str(12*q),new_slow='0',full_nonzero_source=str(9*q),
            required_operators=dict(KT1=[str(6*q)],KA1=[str(6*q)],Kp=[str(6*q),str(12*q)],Kmom=[str(6*q),str(12*q)]))
        return dict(trace=trace,older_source_graph_order=j,older_padded_trace=older,common_pole=c,ledger=ledger,
            original_p_must_cover_ALL_nested_old_and_original_normal_jets=True,
            operator_recipe=dict(method='fresh.shift_recipe',parameter=str(q),geometry=dict(self.geometry),materialized=False),
            oriented_target_recipe=dict(method='edge.primitive_recipe',pole=str(6*q),materialized=False),
            archived_operators='retain all original lower shifts, previous 2q/4q AND previous sourced 6q/12q, plus all distinct oriented primitive constants',
            common_argument='S+K_old_all+K_6q_12q+K_Ie_6q',new_physical_band_roots_width_or_fields_changed=False)

    def postwave_recipe(self,n):
        order(n)
        return dict(method='postwave_rows_envelope',order=n,
            old_rows='SAME previous PostSignedMeanRecipes full tangential C0 rows',
            complete_covariance='FiniteHarmonicWave.covariance through n+3, all six entries and full reconstructed pressure',
            Qcirc_factor=2,order_proved=str(self.orders['E']),actual_input_evaluated=False)

    def fast_recipe(self,n):
        order(n)
        return dict(method='mean.fast_envelope',order=n,tangential_rows=self.postwave_recipe(n+4),
            radial_rows=self.postwave_recipe(n+5),tangential_order=str(self.orders['E']),
            radial_order=str(self.orders['E']+1),actual_gamma='gd-A1 gd, factor 1+KA1 retained',
            raw_factorial_and_original_temporal_reciprocal_kept=True,materialized=False)

    def scalar_recipe(self,n):
        order(n)
        return dict(method='postfast_scalars',order=n,parts=[
            'previous post-generic-slow scalar C+1 from PostSignedMeanRecipes.scalar_recipe',
            'complete next sourced covariance all six entries through n+1',
            'fast_change on FULL pre-source means and ACTUAL next fast fields through n+1'],
            postfast_order=str(self.orders['C0']),slow='SAME response on current scalars n/n+1; whole reserved potential hull G=0',
            current_Haar='PreservedHaarScalarJoin.old_Haar + previous sourced-block C slow + previous generic-block C slow',
            postslow='haar_slow_reset on this exact Haar sum and current slow, inputs through n+2',
            postslow_order=str(self.orders['postslow_scalar']),whole_axial_pressure_moment_order=str(self.orders['whole_axial_pressure_moment']),
            moving_pressure='epsilon partial_Z(Jz+c_rho P), including derivative of c_rho')

    def compact_recipes(self,n,*,physical_power=None):
        d=self.domain_recipe(n,physical_power=physical_power);q=Q(d['common_pole']['q']);out={}
        for e,key in ((0,'A0'),(1,'A1')):
            a=compact_physical_orders(e,n,physical_power,self.h,self.kappa) if physical_power is not None else None
            p=a['integration_order'] if a else self.orders[key+'_integrations']
            src=self.orders['B'] if e==0 else self.orders['E'];target=a['normalized_target'] if a else self.orders['C0']
            r=mean.final.oscillatory_recipe(n,p,e=e,source='Qcirc NEW FINAL full gr' if e==0 else 'Qcirc post-EIGHT-MODE-wave Ez',
                source_rate=src,target_rate=target,kappa=self.kappa,geometry=self.geometry,pole=(12 if e==0 else 6)*q)
            out[key]=dict(raw_kernel=r,orders={k:str(v) for k,v in a.items()} if a else None,
                current_source_recipe='prefix.full_mean with OLD + NEW fast + NEW slow and SAME complete eight-mode covariance' if e==0 else 'next-wave postwave_mean_envelope -> Qcirc tangential rows',
                source_derivatives=r['raw_source_order'],actual_input_or_factorial_evaluated=False)
        return out

    def certificate(self,n,N=0):
        return dict(stage='after-next-eight-mode-sourced-wave-then-actual-fast-and-fresh-slow',
            scalar=self.scalar_recipe(n),normalized_domains=self.domain_recipe(n),normalized_compact=self.compact_recipes(n),
            physical_domains=self.domain_recipe(n,physical_power=N),physical_compact=self.compact_recipes(n,physical_power=N),
            fast=self.fast_recipe(n+3),
            tangential_identity='Efinal=Haar(Epostwave)+(0,A1 Qcirc Epostwave_z)+Rf+Rs, complete unprojected differences',
            radial_identity='Er_final=-rho Pnew-A0(gr_final-rho Pnew); compensated Haar moment zero, A0 source=A0 Qcirc gr_final',
            full_regular_mean_order=str(self.orders['C0']),whole_regular_source_order=str(self.orders['E']),
            nonzero='parent wave-stage E + BOTH new fast/slow mean-wave orientations on COMPLETE H4 wave, residual harmonics through8',
            residual_identity='Rphysical=Rregular+(-A0 Qcirc gr_final,0,+A1 Qcirc Epostwave_z)+TGaussian_current+Tbase_exterior_assigned_once',
            physical_regular=regular_physical(n,self.h,self.kappa),
            Gaussian='same parent Told+(1-psi)CURRENT_Gregular+Tpsi; mean updates append no new Gaussian tail and tails are never inverted',
            physical_tail_flatness_conditional_on_all_specific_actual_inputs=True,
            actual_weighted_high_jet_support_base_validity_absolute_decay_or_blowup_certified=False)

    def certify_actual_inputs(self):
        raise NotImplementedError('Actual all-scale high-jet/root/frame/support/background/base inequalities and infinite diagonal/forcing remain unverified')
