"""Post-generic-signed full mean rows: positive functionals and deferred jets."""
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


signed=load('post_signed_generic',HERE.parent/'core_generic_signed_covariance_join/budgets.py')
final,preserved,scalar,prefix=signed.final,signed.preserved,signed.scalar,signed.prefix
fast=final.fast;fresh=signed.fresh
P,order,exact=prefix.P,prefix.order,prefix.exact


def strict_integrations(source,target,kappa):
    source,target,kappa=map(exact,(source,target,kappa))
    if kappa<=0:raise ValueError('Positive kappa required')
    return max(1,(target-source)//kappa+1)


def rates(kappa=Q(1,100000)):
    a=signed.rates(kappa);k=exact(kappa);B,C,D,C0=(a[x] for x in ('B','C','D','C0'))
    p1=strict_integrations(D,C0,k);p0=strict_integrations(B,C0,k)
    groups=dict(fast_difference=D+1-2*k,slow_difference=C+1-2*k,
                A1=D+p1*k,post_signed_Haar=C0)
    return dict(B=B,C=C,D=D,C0=C0,post_signed_Q1=D,fast=D,fast_radial=D+1,
        postfast_scalar=min(C,B+D,2*D),slow=C,slow_radial=C+1,
        postslow_scalar=min(C+1,2*C,C+2),postslow_axial_moment=C+2,
        tangential_terms=groups,final_tangential=min(groups.values()),
        radial_source=B,radial_compensator=C+1,A0=B+p0*k,final_radial=C0,
        final_regular_mean=C0,A1_integrations=p1,A0_integrations=p0,
        A1_source_order_zero=2*p1+3,A0_source_order_zero=2*p0+3,
        complete_nonzero_residual_or_physical_decay_certified=False)


def fast_envelope(n,rows,KT1,KA1,time_inverse):
    """SAME actual fast map; rows is the DERIVED post-signed Q1 recipe.

    Functional arguments are positive norm programs, not new actual
    fine-row hypotheses. KT1/KA1 MUST be original raw recipes at 2q.
    """
    order(n);a=rows(n+4);b=rows(n+5)
    tang=comb(n+5,5)*24*2**(n+4)*time_inverse*factorial(n+4)
    rad=comb(n+5,5)*24*2**(n+5)*time_inverse*factorial(n+5)*KT1
    return dict(u=tang*a['Etheta'],z=tang*(1+KA1)*a['Ez'],x=rad*b['Ez'])


def postfast_scalar_envelope(n,old_scalars,delta_cov_n,delta_cov_next,old_full,fast_inc,M):
    """Three CURRENT scalar channels: old postslow + wave + actual fast.

    old_full and fast_inc are through n+1; the old field is FULL, not
    Haar-only. No old scalar or covariance entry is reset to zero.
    """
    order(n);w=scalar.wave_change(n,delta_cov_n,delta_cov_next,M)
    f=scalar.fast_change(n,old_full,fast_inc,M)
    return {key:old_scalars[key]+w[key]+f[key] for key in ('P','Jtheta','Jz')}


def old_Haar_envelope(previous_Haar,previous_slow):
    """New signed wave changes no means; new complete fast has zero Haar."""
    if set(previous_Haar)!=set(prefix.STATE) or set(previous_slow)!=set(prefix.STATE):raise ValueError('All three Haar/slow components required')
    return {key:previous_Haar[key]+previous_slow[key] for key in prefix.STATE}


def final_row_envelope(n,post_signed_Haar,A1_tail,old_full,fast_inc,slow_inc,fixed):
    """After exact c Nt cancellation, with full pressure/time/viscosity.

    The complete identity is Efinal=Haar(EpostSigned)+(0,A1 Fz)+Rf+Rs.
    Thus these full-row norms need no extra Q1 factor on Rf or Rs.
    All mean/fixed inputs are through n+3. Actual A1 velocities remain
    in fast_inc in BOTH nonlinear differences.
    """
    order(n)
    f=fast.difference(n,old_full,fast_inc,fixed)
    after={key:old_full[key]+fast_inc[key] for key in prefix.STATE}
    s=fast.difference(n,after,slow_inc,fixed)
    return dict(Etheta=post_signed_Haar['Etheta']+f['Etheta']+s['Etheta'],
                Ez=post_signed_Haar['Ez']+A1_tail+f['Ez']+s['Ez'],fast_difference=f,slow_difference=s)


def mean_feedback(n,wave,increment,*,graph,reciprocal,normal,max_harmonic=2):
    """Mean-change piece ONLY of the existing harmonic-aware residual.

    Complete cut wave through n+3, actual three-component increment
    through n+3. Phase advection costs epsilon^-1/2; graph transport's
    gain 1/2-kappa is dropped, retaining source order of the increment.
    """
    order(n)
    if type(max_harmonic) is not int or max_harmonic<1:raise ValueError('Positive retained maximum input harmonic required')
    if set(increment)!=set(prefix.STATE):raise ValueError('All three actual mean increment components required')
    d=n+4;Cmw=2*(d*graph+d+reciprocal)+2*max_harmonic*normal
    return Cmw*sum(increment.values(),wave*0)*wave


def operator_ledger(wave_pole,original_pole):
    q,p=map(exact,(wave_pole,original_pole))
    if min(q,p)<0:raise ValueError('Nonnegative finite poles required')
    # Current old full rows already have pressure/moment domain 12p;
    # a wave-only q is not automatically sufficient for these mean calls.
    q=max(q,6*p)
    return dict(q=str(q),old_complete_rows=str(12*p),old_current_mean_fields=str(6*p),
        complete_signed_covariance=str(2*q),post_signed_rows=str(2*q),actual_fast_fields=str(2*q),
        postfast_source=str(4*q),final_radial_source=str(4*q),slow_field_pole='0',
        required_domains=dict(KT1=[str(2*q)],KA1=[str(2*q)],Kp=[str(2*q),str(4*q)],Kmom=[str(2*q),str(4*q)]),
        replacement_estimate_parameter=str(q/3),
        lower_pole_embedding='max(1,L)^ceil(4q), before products and moment calls')


def mean_shift(wave_pole,original_pole,geometry):
    ledger=operator_ledger(wave_pole,original_pole)
    shift=fresh.shift_recipe(Q(ledger['replacement_estimate_parameter']),Q(geometry['L']),Q(geometry['width']),Q(geometry['gap']))
    return dict(ledger=ledger,recipe=shift,argument='all nonnegative proof polynomials at S+K_mean; physical S and accepted band unchanged')


def dependency_trace(output_order,A1_integrations,A0_integrations):
    """Trace ONLY newly introduced calls; old DAGs remain opaque references.

    This is an executable call graph, not a guessed total original-root
    jet count. Large integer indices are traversed symbolically, without
    evaluating coefficient polynomials, factorials or original operators.
    """
    order(output_order);order(A1_integrations);order(A0_integrations)
    if min(A1_integrations,A0_integrations)<1:raise ValueError('Positive compact integration counts required')
    calls={}
    def seen(name,n):calls.setdefault(name,set()).add(n)
    def covariance(n):seen('new_covariance',n);seen('next_signed_target_input',n+1)
    def row(n):
        seen('post_signed_row',n);seen('old_final_row_recipe',n)
        covariance(n+3)
    def fast_field(n):seen('new_fast_field',n);row(n+4);row(n+5)
    def scalar(n):
        seen('new_postfast_scalar',n);seen('old_postslow_scalar_recipe',n)
        covariance(n);covariance(n+1);fast_field(n+1)
    def slow_field(n):seen('new_slow_field',n);scalar(n);scalar(n+1)
    def final_fields(n):
        seen('old_final_field_recipe',n);fast_field(n);slow_field(n)
    def gr(n):seen('final_full_gr',n);final_fields(n+3)
    n=output_order
    # Complete tangential difference fields through n+3.
    final_fields(n+3)
    # Exact compact residual estimates at the requested output order.
    J1=n+2*A1_integrations+3;J0=n+2*A0_integrations+3
    row(J1);gr(J0)
    return dict(output_order=n,A1_source_order=J1,A0_source_order=J0,
        maxima={key:max(value) for key,value in sorted(calls.items())},
        old_dependency_DAGs_expanded=False,actual_jet_or_profile_values_evaluated=False)


class PostSignedMeanRecipes:
    """One-way SAME-prefix recipes; do not instantiate giant fine callbacks."""
    def __init__(self,current_signed,*,wave_pole,original_pole,geometry,kappa=Q(1,100000)):
        self.current_signed=current_signed;self.final=current_signed.final;self.haar=current_signed.haar
        if self.final.means is not self.haar.means:raise ValueError('One identical current prefix required')
        self.orders=rates(kappa);self.kappa=exact(kappa);self.geometry=dict(geometry)
        self.wave_pole,self.original_pole=map(exact,(wave_pole,original_pole))
        retained=self.final.means.reweighting
        if any(Q(geometry[key])!=Q(retained[key]) for key in ('L','width','gap')) or self.original_pole!=Q(retained['original_pole']):
            raise ValueError('Retain the original same-prefix pole parameter, width, hull and gap')
        operator_ledger(self.wave_pole,self.original_pole)

    def common_pole(self,n):
        r=self.orders;trace=dependency_trace(n,r['A1_integrations'],r['A0_integrations'])
        N=trace['maxima']['new_covariance'];p=self.original_pole
        # Current rows at pole 12p through ALL requested finite jets give
        # TargetBudget.edge allowance d=12p; choose p for the old DAG
        # before making this inference. This does not evaluate those jets.
        d=12*p
        groups=dict(retained_wave=self.wave_pole,old_full_mean_domain=6*p,
                    original_curl=6*N+9,new_signed_curl=6*N+d+6)
        q=max(groups.values())
        return dict(q=str(q),groups={k:str(v) for k,v in groups.items()},
            maximum_covariance_input=N,maximum_current_target_input=trace['maxima']['next_signed_target_input'],
            target_allowance=str(d),original_pole_covers_ALL_nested_old_jets_required=True,
            actual_old_jet_bounds_verified=False)

    def operators(self,n):
        q=Q(self.common_pole(n)['q']);ledger=operator_ledger(q,self.original_pole)
        return dict(ledger=ledger,recipe=dict(method='fresh.shift_recipe',
            original_parameter=str(q/3),geometry=dict(self.geometry),materialized=False),
            argument='all nonnegative proof polynomials at S+K_mean; physical S and accepted band unchanged')

    def post_signed_rows(self,n):
        order(n)
        return dict(stage='after-next-generic-signed-before-new-fast',order=n,
            expression='Q1(old final rows + covariance_rows on COMPLETE generic covariance_budget)',
            old_rows=dict(method='FinalMeanRowClosure.tangential_recipe',order=n,materialized=False),
            covariance=dict(method='generic covariance_budget -> fast.covariance_rows',input_order=n+3,
                entries='all six, old-error cross + full curl cross + signed square + principal cross'),
            Q1_norm_factor=2,order_proved=str(self.orders['D']),
            original_current_stage=self.current_signed.target_inputs(0)['stage'])

    def fast_recipe(self,n,*,common_request_order=None):
        order(n)
        request=n if common_request_order is None else common_request_order
        return dict(method='fast_envelope',order=n,tangential_rows_order=n+4,radial_rows_order=n+5,
            tangential_order=str(self.orders['D']),radial_order=str(self.orders['D']+1),
            source_method='PostSignedMeanRecipes.post_signed_rows',
            actual_axial_velocity='g_d-A1 g_d; factor 1+KA1 retained',
            operators=self.operators(request)['ledger']['required_domains'],materialized=False)

    def scalar_recipe(self,n):
        order(n)
        return dict(method='postfast_scalar_envelope',order=n,order_proved=str(self.orders['C']),
            parts=['PreservedHaarScalarJoin.postslow(n)',
                'wave_change on COMPLETE generic signed covariance through n+1',
                'fast_change on FULL old current final means and new fast fields through n+1'],
            fixed_five_map='response(n, current postfast scalars n/n+1, SAME original closed-hull profiles)',
            old_Haar='previous old_Haar + previous latest slow; signed wave unchanged and new fast all zero Haar',
            final_postslow='haar_slow_reset through n+2 with fresh current slow response',
            final_postslow_order=str(self.orders['postslow_scalar']),
            whole_axial_pressure_moment='epsilon partial_Z(Jz_new+c_rho Pnew), including derivative of c_rho',
            final_axial_moment_order=str(self.orders['postslow_axial_moment']))

    def compact_recipes(self,n):
        order(n);r=self.orders;q=Q(self.common_pole(n)['q'])
        A1=final.oscillatory_recipe(n,r['A1_integrations'],e=1,source='Q1 post-signed Ez',source_rate=r['D'],target_rate=r['C0'],kappa=self.kappa,geometry=self.geometry,pole=2*q)
        A0=final.oscillatory_recipe(n,r['A0_integrations'],e=0,source='Q1 NEW FINAL full gr',source_rate=r['B'],target_rate=r['C0'],kappa=self.kappa,geometry=self.geometry,pole=4*q)
        A1['source_recipe']=dict(method='PostSignedMeanRecipes.post_signed_rows',order=A1['raw_source_order'],materialized=False)
        A0['source_recipe']=dict(method='prefix.full_mean on old final fields + NEW fast + NEW slow and complete updated wave',order=A0['raw_source_order'],field_order=A0['raw_source_order']+3,materialized=False)
        return dict(A1=A1,A0=A0)

    def certificate(self,n):
        r=self.orders;compact=self.compact_recipes(n)
        return dict(order=n,stage='after-next-generic-signed-then-new-fast-then-fresh-current-slow',
            exact_tangential_identity='Efinal=Haar(EpostSigned)+(0,A1 Q1 EpostSigned_z)+Rf+Rs',
            exact_radial_identity='Er_final=-rho Pnew-A0(gr_final-rho Pnew); A0 source has zero Haar moment',
            full_regular_mean_order=str(r['C0']),postslow_scalar_order=str(r['C']+1),
            scalar_recipe=self.scalar_recipe(n),fast_recipe=self.fast_recipe(n+3,common_request_order=n),compact=compact,
            operator_shift=self.operators(n),common_pole=self.common_pole(n),
            request_trace=dependency_trace(n,r['A1_integrations'],r['A0_integrations']),
            q_scope='ONE finite q must cover ALL displayed and nested coefficient/profile requests, not only the low output order',
            old_frozen_DAGs_and_actual_A1_velocities_retained=True,
            actual_input_validity_complete_nonzero_residual_physical_tail_or_blowup_certified=False)

    def certify_actual_inputs(self):
        raise NotImplementedError('Deferred raw high-jet, common-pole weighted field/profile/support inputs and complete nonzero/physical residual join remain unverified')
