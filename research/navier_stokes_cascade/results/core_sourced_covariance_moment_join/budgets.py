"""Generated sourced complete curl -> covariance -> full pressure moments."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


tangent=load('covariance_moment_tangent',HERE.parent/'core_tangent_propagator_closure/budgets.py')
prefix=load('covariance_moment_prefix',HERE.parent/'core_signed_prefix_source_budget/budgets.py')
P,exact,order=tangent.P,tangent.exact,tangent.transfer.order
edge=prefix.haar.edge


def convert(p):return P(p.c)


def through(raw,n,dimension=5):
    order(n)
    if type(dimension) is not int or dimension<1:raise ValueError('Positive coefficient dimension required')
    return sum((Q(dimension**j,factorial(j))*raw(j) for j in range(n+1)),P())


def curl_budget(n,principal,potential,potential_next,graph,reciprocal):
    order(n)
    remainder=2*(n+1)*(graph+1)*potential_next+reciprocal*potential
    return dict(T=principal,R=remainder,U=principal+remainder)


class GeneratedWave:
    """One-way composition from the BEFORE prefix, never its own output.

    All callbacks bound the SAME original fields/common chart. The source
    prefix retains its original-root, coefficient, weighted-map and support
    premises. This constructor does not certify those field enclosures.
    """
    def __init__(self,p,selector_inputs,ell,before,*,F,G,btilde,path):
        self.before=before;self.p=tangent.normal.frame_inputs(p)
        source=before.inverse_source(tangent.transfer.inverse)
        self.modes={}
        for m in (1,2):
            _,f,c=tangent.forced_from_frame(p,selector_inputs,ell,F=F,G=G,btilde=btilde,path=path,source=source,harmonic=m)
            self.modes[m]=f;self.certificate=c
        # A common harmonic-l1 source bound bounds each signed coefficient.
        # Reuse equal majorants, without equating actual mode coefficients.
        self.harmonics=(-2,-1,1,2)

    def cutoff_raw(self,n):
        order(n);return factorial(n)*convert(self.before.linear(n)['cutoff'])

    @lru_cache(None)
    def mode(self,n,m):
        order(n)
        if type(m) is not int or m not in self.harmonics:raise ValueError('First sourced signed harmonics only')
        f=self.modes[abs(m)]
        cut=lambda kind,j:f.cut(kind,self.cutoff_raw,j)
        A=through(lambda j:cut('potential',j),n)
        Anext=through(lambda j:cut('potential',j),n+1)
        fixed=self.before.fixed(n)
        out=curl_budget(n,through(lambda j:cut('velocity',j),n),A,Anext,
                        convert(fixed['graph']),convert(fixed['r']))
        out['Pi']=through(lambda j:cut('pressure',j),n)
        return out

    @lru_cache(None)
    def addition(self,n):
        order(n)
        return {k:sum((self.mode(n,m)[k] for m in self.harmonics),P()) for k in ('T','R','Pi','U')}

    @lru_cache(None)
    def covariance(self,n):
        order(n);old=convert(self.before.wave(n)['U']);new=self.addition(n)['U']
        return 2*old*new+new*new

    def unweighted_covariance(self,n,pole,width):
        pole=exact(pole)
        if pole<0:raise ValueError('Nonnegative incoming common pole required')
        return edge.atlas.edge_constant(6*pole,width)*self.covariance(n)


def moment_budget(n,D,Dnext,a,d,c_rho):
    """Unweighted scalar jets from complete covariance norms through n+1.

    [a,d] contains ALL changed-wave supports, not just the reserved mean
    hull. c_rho is the full unweighted moving second-moment jet norm.
    An extra epsilon in axial-flux terms is conservatively bounded by 1.
    These are CHANGES; the incoming current defects are not reset.
    """
    order(n);a,d=map(exact,(a,d))
    if not 0<a<d:raise ValueError('Full positive radial support hull required')
    values=(D,Dnext,c_rho)
    if not all(isinstance(v,P) for v in values):
        if any(exact(v)<0 for v in values):raise ValueError('Nonnegative norm bounds required')
    M={j:(d-a)*max(a**j,d**j) for j in (-1,0,1,2)}
    axial=(n+1)*Dnext
    Pdef=2*M[-1]*D+M[0]*axial
    Qgr=2*M[1]*D+M[2]*axial
    Jtheta=M[2]*D
    Jz=2*M[1]*D+Q(1,2)*M[2]*axial
    return dict(P=Pdef,Qgr=Qgr,Jtheta=Jtheta,Jz=Jz,
                pressure_moment=Q(1,2)*Qgr+c_rho*Pdef,
                whole_axial_moment=Jz+c_rho*Pdef)


def generated_moments(wave,n,pole,width,c_rho):
    a,d=wave.p['Rmin'],wave.p['Rmax']
    return moment_budget(n,wave.unweighted_covariance(n,pole,width),
                         wave.unweighted_covariance(n+1,pole,width),a,d,c_rho(n))


def generated_residual_moments(wave,n,pole,width,c_rho):
    """After factoring the extra axial epsilon; covariance through n+2."""
    row=generated_moments(wave,n+1,pole,width,c_rho)
    return dict(theta=(n+1)*row['Jtheta'],axial=(n+1)*row['whole_axial_moment'])


def next_pole_ledger(p):
    p=exact(p)
    if p<0:raise ValueError('Nonnegative original finite pole required')
    after=prefix.pole_ledger(3*p)
    return dict(old_wave=p,old_means=2*p,source_and_sourced_wave=3*p,
                covariance_cross=4*p,covariance_square=6*p,covariance=6*p,
                next_common_input=3*p,next_mean_ledger=after,
                required_operator_poles=dict(KT1=[6*p],KA1=[6*p],Kp=[6*p,12*p],Kmom=[6*p,12*p]))


def operator_coverage(p,declared):
    """Match proof domains, not a guessed monotonicity of operator norms."""
    required=next_pole_ledger(p)['required_operator_poles']
    if set(declared)!=set(required):raise ValueError('All weighted operator proof domains required')
    domains={k:{exact(x) for x in v} for k,v in declared.items()}
    if any(x<0 for v in domains.values() for x in v):raise ValueError('Nonnegative pole domains required')
    missing={k:sorted(set(v)-domains[k]) for k,v in required.items()}
    return dict(pass_all=not any(missing.values()),missing=missing)


def orders(kappa):
    k=exact(kappa)
    if not 0<k<Q(1,4):raise ValueError('Original positive kappa range required')
    B=1-k;C=B+Q(1,2)
    return dict(source=B,complete_sourced_velocity=B,curl=B+Q(1,2)-k,
                covariance=C,moments=C,residual_moments=C+1,
                self_square_slack=2*B-C,curl_slack=Q(1,2)-k,
                complete_residual_improvement_certified=False)
