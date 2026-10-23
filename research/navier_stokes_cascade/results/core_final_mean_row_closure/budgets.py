"""Conditional final mean-row closure with explicit deferred A0/A1 recipes."""
from fractions import Fraction as Q
from math import ceil
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


localized=load('final_mean_localized',HERE.parent/'core_localized_stress_assembly/assembly.py')
haar,fast,prefix=localized.haar,localized.fast,localized.prefix
P,order,exact=prefix.P,prefix.order,prefix.exact


def rates(kappa=Q(1,100000)):
    r=localized.rates(kappa);k=exact(kappa);B,C,D=r['B'],r['C'],r['D']
    q0=max(1,(C-B)//k+1)
    groups=dict(fresh_fast_difference=D+1-2*k,fresh_slow_difference=C+1-2*k,
                fresh_A1=D+2*k)
    radial=dict(fast_time_of_beta=B+1,covariance_graph=1-k,
                radial_mean_flux=min(B+2-k,2*B+2-k),
                axial_flux=min(B+2,2*B+2,Q(2)),
                centripetal=min(B,2*B,Q(1)),radial_viscosity=B+2-2*k)
    return dict(B=B,C=C,D=D,final_Q1=min(groups.values()),Q1_terms=groups,
                radial_source_groups=radial,radial_source=min(radial.values()),
                final_tangential=C,A0_integration_order=q0,A0_source_order_zero=2*q0+3,
                A0_order=B+q0*k,radial_compensator=B+C,final_radial=C,
                full_angular_zero_regular=C,joined_regular=min(C,D),
                original_actual_inputs_certified=False,complete_physical_tail_or_blowup_certified=False)


def oscillatory_recipe(n,p,*,e,source,source_rate,target_rate,kappa,geometry,pole):
    """Original raw kernel estimate (8), not an arbitrary operator norm.

    The two factor covers Q1; 6^n converts raw output to five-variable
    factorial through norms. The raw input uses J! times its through norm.
    Giant powers/factorials are finite exact instructions, not evaluated.
    """
    order(n);order(p);e=exact(e);k=exact(kappa);pole=exact(pole)
    if p<1 or e not in (0,1) or pole<0:raise ValueError('Original e=0/1 and positive integration order required')
    if exact(source_rate)+p*k<=exact(target_rate):raise ValueError('A strict retained exponent margin is required')
    J=n+2*p+3;w=Q(geometry['width']);gap=Q(geometry['gap']);L=Q(geometry['L'])
    if L<=0 or not 0<w<=1 or not 0<gap<=L/2:raise ValueError('Original positive hull, width and gap required')
    return dict(operator='A'+str(e),output_order=n,integration_order=p,raw_source_order=J,
                raw_output_scalar=2*6**n*8,two_power=2*n+J+p+3,three_power=p,source_factorial=J,
                factors=['L',f'N_{n}(k_{e} chi_U)',f'N_{J}(f_{e})',source],
                source_pole=str(pole),source_unweighting=f'Eedge({pole},original width)',
                inverse_output_edge_floor_power=ceil((w*w+4)/(gap*gap)),
                original_geometry=dict(geometry),
                floor_factor='(1+means.time_inverse(0))^ceil(original_rho_g*integration_order)',
                source_rate=str(exact(source_rate)),target_rate=str(exact(target_rate)),
                epsilon_order=str(exact(source_rate)+p*k),
                input_jet_derivatives='frozen M,c; whole-line smooth zero extension; original raw kernel profiles',
                zero_mode='zero Haar source for A1; compensated integrated Haar moment for A0',
                materialized_actual_input_bound=False)


def remaining_row_envelope(n,old,fast_inc,slow_inc,fixed):
    """Positive polynomial on actual increment norm arguments.

    The sharp fast increment argument must be the derived deferred recipe,
    not the coarse B norm relabeled D. Supports no temporal principal term:
    c Nt fast cancels the OLD Q1 rows, c Nt slow=0 exactly.
    """
    order(n)
    f=fast.difference(n,old,fast_inc,fixed)
    after={k:old[k]+fast_inc[k] for k in prefix.STATE}
    s=fast.difference(n,after,slow_inc,fixed)
    return dict(fast_difference=f,slow_difference=s,
                Q1={k:2*(f[k]+s[k]) for k in ('Etheta','Ez')})


class FinalMeanRowClosure:
    """Same prefix, actual fresh fast/slow operations, no fine-input callback."""
    def __init__(self,means,kappa=Q(1,100000)):
        self.means=means;self.kappa=exact(kappa);self.orders=rates(kappa)
        self.fast=fast.FastRowGain(means,kappa);self.haar=haar.HaarStressJoin(means,kappa)
        self.embedding=self.haar.embedding

    def A1_recipe(self,n):
        p=Q(self.means.reweighting['original_pole'])
        r=oscillatory_recipe(n,2,e=1,source='Q1 postwave Ez',source_rate=self.orders['D'],
            target_rate=self.orders['C'],kappa=self.kappa,geometry=self.means.reweighting,pole=12*p)
        # Source at only n+7, but it contains the parent's original high
        # A1 recipe. Do not truncate that nested dependency.
        r['source_recipe']=self.fast.postwave_recipe(r['raw_source_order'])
        return r

    def tangential_recipe(self,n):
        order(n);N=n+3;m=self.means;s=self.haar.scalars
        return dict(stage='after-fresh-fast-and-slow',input_order=N,
            old_means={k:(self.embedding*v).record() for k,v in m.data('initial',N).items()},
            sharp_fast_increment=self.fast.fast_recipe(N),
            sharp_slow_increment=dict(method='ScalarGainBudget.slow_increment',order=N,
                components=list(prefix.STATE),embedding=str(self.embedding),
                expression='embedding * scalars.slow_increment(input_order)',materialized=False),
            fixed={k:v.record() for k,v in m.data('fixed',N).items()},
            difference_operator='remaining_row_envelope on the displayed norm arguments',
            temporal_identity='Q1 Efinal=(0,A1 Q1 Ez_postwave)+Q1 Rfast+Q1 Rslow',
            exact_A1=self.A1_recipe(n),Q1_order=str(self.orders['final_Q1']),
            Haar_order=str(self.orders['C']),complete_regular_tangential_order=str(self.orders['C']),
            leading_stress_identity='same canonical localized assembly, unchanged wave during both mean updates',
            additive_tail_ledger='unchanged and excluded from this regular row recipe')

    def radial_recipe(self,n):
        order(n);m=self.means;r=self.orders;p=Q(m.reweighting['original_pole'])
        a=oscillatory_recipe(n,r['A0_integration_order'],e=0,source='Q1 current final gr',
            source_rate=r['radial_source'],target_rate=r['C'],kappa=self.kappa,geometry=m.reweighting,pole=12*p)
        # A specific existing polynomial expression, not a new estimate
        # callback. Evaluate only if the original huge input data are given.
        a['source_recipe']=dict(method='FinalMeanRowClosure.current_radial_source_polynomial',order=a['raw_source_order'],component='gr',
            norm='full graph radial source, all six covariance and actual mean terms',
            field_input_order=a['raw_source_order']+3,
            embedding=str(self.embedding),
            temporal_or_pressure_inverse_added=False)
        return dict(stage='after-fresh-fast-and-slow',exact_identity='Er=-rho Pnew-A0(gr-rho Pnew)',
            zero_mode_identity='int bar(gr-rho Pnew) dR=0; A0 h=A0 Q1 gr',
            compensated_unit_mass_required='int rho dR=1, Pnew=int bar gr dR',
            compensator=dict(method='HaarStressJoin.radial_Haar_compensator',order=n,
                expression='monomial_upper(rho * actual postslow P polynomial)',materialized=False),
            compensator_order=str(r['radial_compensator']),exact_A0=a,
            complete_regular_radial_order=str(r['C']),
            full_pressure_field_not_order_promoted=True)

    def current_radial_source_polynomial(self,n):
        """Same complete gr with lower-pole fields embedded before products.

        This exact existing full_mean polynomial is deferred at order
        100005; it is not an incoming fine norm or an input enclosure.
        """
        order(n);N=n+3;m=self.means;E=self.embedding
        state={k:E*v for k,v in m.final(N).items()}
        Uv=m.data('wave',N)['U']
        return prefix.full_mean(n,dict(m.data('fixed',N),**state,w=E*Uv*Uv))['gr']

    def certificate(self,n):
        return dict(tangential=self.tangential_recipe(n),radial=self.radial_recipe(n),
            angular_zero_regular_order=str(self.orders['C']),retained_nonzero_regular_order=str(self.orders['D']),
            combined_regular_order=str(self.orders['joined_regular']),
            actual_input_validity_or_full_physical_residual=False)

    def actual_input_certificate(self,*args,**kwargs):
        raise NotImplementedError('Deferred polynomial recipes do not independently verify actual high-jet inputs, supports or physical tail bounds')
