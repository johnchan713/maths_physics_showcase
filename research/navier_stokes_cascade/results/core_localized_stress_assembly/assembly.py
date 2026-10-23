"""Exact principal-stress assembly and conditional postwave row recipes.

These functions evaluate supplied fields/majorants. They do not certify
the original continuum column, root, support or profile inputs.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


haar=load('localized_stress_haar',HERE.parent/'core_haar_row_stress_join/budgets.py')
fast,prefix=haar.fast,haar.prefix


def principal_defect(target,labels):
    """Unconditional algebra, including both possible localization defects.

    Each label denotes a slow BOX with both signs inside H, not a sign.
    It supplies eta, H, y, its unpartitioned chart target T and physical
    covariance scale a=Q^(-2A)*epsilon. Matrices may be exact symbolic
    objects. The caller must separately justify the absence of cross-label
    angular covariances. No test result substitutes for that hypothesis.
    """
    if not labels:raise ValueError('A nonempty active cover is required')
    zero=target*0;actual=zero;root=zero;normalization=zero;cover=0
    for lab in labels:
        if set(lab)!=set(('eta','H','y','T','a')):raise ValueError('Complete box data required')
        w=lab['eta']**2;a=lab['a'];cover+=w
        actual+=w*a*(lab['H']*lab['y'])
        root+=w*a*(lab['H']*lab['y']-lab['T'])
        normalization+=w*(a*lab['T']-target)
    coverage=(cover-1)*target
    return dict(actual=actual,defect=actual-target,root=root,
                normalization=normalization,coverage=coverage,cover=cover)


def rates(kappa=Q(1,100000)):
    r=haar.rates(kappa)
    return dict(B=r['B'],C=r['C'],D=r['D'],
                postwave_Haar=r['C'],postwave_Q1=r['D'],
                postwave_complete_tangential=min(r['C'],r['D']),
                final_Haar=r['C'],final_radial_Haar=r['B']+r['C'],
                full_pointwise_radial_or_tail_bound=False,
                independent_actual_input_certification=False,
                complete_physical_residual_or_blowup=False)


class LocalizedStressJoin:
    """Recipe for the SAME canonical prefix, before fresh fast/slow means.

    The exact matched assembly is an analytical premise supplied by the
    displayed construction, not by a Boolean or an incoming fine-row norm.
    The finite original A1 recipe remains explicit and unmaterialized.
    """
    def __init__(self,means,kappa=Q(1,100000)):
        self.means=means;self.kappa=prefix.exact(kappa)
        self.fast=fast.FastRowGain(means,kappa)
        self.haar=haar.HaarStressJoin(means,kappa)
        self.orders=rates(kappa)

    def postwave_envelope(self,n,*,principal_defect,fine_flux):
        """Keep the leading defect until the exact assembly is substituted."""
        prefix.order(n);N=n+3;m=self.means;E=self.haar.embedding
        # FreshMeanBudget.initial is before.final: old means, new wave,
        # BEFORE its own fast/slow increments. Do not substitute m.final.
        state={k:E*v for k,v in m.data('initial',N).items()}
        W=m.data('wave',N)['U'];W=E*W*W
        h=haar.haar_rows(n,state,m.data('fixed',N),W,E*principal_defect,E*fine_flux)
        f=self.fast.regular_rows(n)
        combined={key:f['rows'][key]+h[hkey]+h[lkey]
                  for key,hkey,lkey in (('Etheta','fine_theta','leading_theta'),('Ez','fine_z','leading_z'))}
        return dict(Haar=h,Q1_regular=f,combined_regular=combined,
                    exact_A1=self.fast.postwave_recipe(n)['exact_A1'],
                    leading_defect_order=Q(1),matched_order=self.orders['D'],
                    stage='postwave-before-fresh-fast-and-slow',
                    additive_tail_ledger='unchanged; excluded from this regular row recipe')

    def actual_input_certificate(self,*args,**kwargs):
        raise NotImplementedError('Original actual continuum column/root/profile/support premises still require independent validation')
