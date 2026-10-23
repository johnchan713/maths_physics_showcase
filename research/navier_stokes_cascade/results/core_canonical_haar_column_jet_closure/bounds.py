"""Canonical principal Haar-column jets from the SAME homogeneous pulse recipe."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
entrance=load('haar_column_entrance','core_canonical_entrance_jet_closure/bounds.py')
gaussian=load('haar_column_gaussian','gaussian_mass_foundation/gaussian_bounds.py')
P,exact,order,product=entrance.P,entrance.exact,entrance.order,entrance.product

class HaarColumnBudget:
    """Raw U_j bound |D^j t|/P on the original pulse, uniformly in both signs.

    No Haar-column jet callback, new root or integration domain is supplied.
    Gaussian curvature is a positive fixed constant, never a function of S.
    """
    def __init__(self,homogeneous,curvature,dimension=5):
        self.homogeneous=homogeneous;self.curvature=exact(curvature)
        if self.curvature<=0:raise ValueError('Original strictly positive fixed-tilt curvature required')
        if type(dimension) is not int or dimension<3:raise ValueError('Coefficient chart includes the three slow variables')
        if not callable(getattr(homogeneous,'velocity',None)):raise TypeError('Raw homogeneous pulse jet recipe required')
        self.dimension=dimension
        self.haar_jacobian_upper=Q(2) # original |det(vr,vt)|=1+(sqrt(2)-1)^2 < 2
        self.gaussian_integral_upper=2*max(Q(1),1/self.curvature) # sqrt(pi/c) <= this

    def pulse(self,j):
        order(j);p=self.homogeneous.velocity(j)
        if not isinstance(p,P):raise TypeError('Positive raw homogeneous jet polynomials required')
        return p

    @lru_cache(None)
    def raw_tensor_entry(self,j):
        """One normalized unrooted PRINCIPAL covariance entry; not full curl covariance."""
        order(j)
        return self.haar_jacobian_upper*self.gaussian_integral_upper*product(self.pulse,self.pulse,j)

    def matrix(self,j):
        """J_j(H/s0) in induced matrix infinity norm, with two sign columns."""
        order(j)
        return Q(2*self.dimension**j,factorial(j))*self.raw_tensor_entry(j)

    def wave(self,j,cutoff):
        return self.homogeneous.wave(j,cutoff,self.dimension)

def from_original(p,selector_inputs,ell,*,F,G,btilde,path,tilt,dimension=5):
    """Derive curvature and U_j on the same accepted original tail, no reselection.

    tilt is the already fixed ACTUAL positive tilt. p['u'] is only its upper
    bound: a positive curvature cannot be inferred from an upper bound alone.
    """
    p=entrance.reference.normal.frame_inputs(p);u=exact(tilt)
    if not 0<u<=p['u']:raise ValueError('Actual fixed positive tilt must lie below the retained upper bound')
    phase,h,c,tc=entrance.from_original(p,selector_inputs,ell,F=F,G=G,btilde=btilde,path=path)
    lower,upper=gaussian.curvature(u,p['lam_min'],p['lam_max'])
    column=HaarColumnBudget(h,lower,dimension)
    return phase,column,c,tc,dict(actual_fixed_tilt=u,tilt_upper=p['u'],curvature_lower=lower,curvature_upper=upper,
        actual_curvature_or_primitive_suprema_numerically_certified=False)

def root_join(template,column,cutoff):
    """Rebind H and wave NORM slots together; keep original canonical value gaps."""
    roots=entrance.inverse.roots
    if not isinstance(template,roots.RootBudget) or not isinstance(column,HaarColumnBudget) or not isinstance(column.homogeneous,entrance.CanonicalEntranceBudget):
        raise TypeError('Original root budget and canonical pulse-derived column recipe required')
    return roots.RootBudget(H=column.matrix,target=template.target,wave=lambda j:column.wave(j,cutoff),
        curl=template.curl,loghalf=template.loghalf,K0=template.K0,beta=template.beta,R=template.R,
        mass=template.mass,width=template.width,trace=template.trace)

def certificate():
    return dict(parent_commit='f749704c9ee88e0758640e189a2c024dbe1f579b',
        result='conditional explicit canonical Haar-column mixed-jet closure and joint H/W root-budget substitution',
        integral='Hsigma=(J/2)*(integral chi^2)*ci*integral_0^L psi(v/L)^2*t_r*t_tan dv; s0=r0^2/sqrt(L), ci=2r0/L',
        raw_bound='|D^I(Hsigma/s0)| <= Jmax*sqrt(pi/c)*(U*U)_j; Jmax=2, sqrt(pi/c)<=2 max(1,1/c)',
        matrix='J_j(H/s0)<=2*(dimension^j/j!)*Jmax*G*(U*U)_j; auxiliary derivatives vanish AFTER Haar averaging',
        moving_maps='parameter-dependent torus translations integrate out exactly; normalized finite-cover Haar invariance, original flat cutoff and frozen label/length retained',
        actual_fixed_positive_tilt_required_not_just_upper_bound=True,
        full_curl_covariance_equals_principal_columns=False,
        new_actual_suprema_complete_mean_residual_or_blowup_certified=False)
