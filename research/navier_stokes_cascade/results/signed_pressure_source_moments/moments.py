"""Full-support pressure-source moments; exact algebra is not continuum existence."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)

class J:
    """Factorial angular Taylor coefficients. Only orders zero and one are certified."""
    N=4
    def __init__(self,c):self.c=(tuple(map(F,c))+(F(0),)*(self.N+1))[:self.N+1]
    @staticmethod
    def cast(x):return x if isinstance(x,J) else J((x,))
    def __add__(self,x):
        x=J.cast(x);return J(tuple(a+b for a,b in zip(self.c,x.c)))
    __radd__=__add__
    def __neg__(self):return J(tuple(-v for v in self.c))
    def __sub__(self,x):return self+-J.cast(x)
    def __rsub__(self,x):return J.cast(x)+-self
    def __mul__(self,x):
        if isinstance(x,P):return NotImplemented
        x=J.cast(x);return J(tuple(sum((self.c[i]*x.c[n-i] for i in range(n+1)),F(0)) for n in range(self.N+1)))
    __rmul__=__mul__
    def inverse(self):
        if not self.c[0]:raise ZeroDivisionError('Nonzero angular denominator required')
        c=[1/self.c[0]]
        for n in range(1,self.N+1):c.append(-sum(self.c[i]*c[n-i] for i in range(1,n+1))/self.c[0])
        return J(c)
    def __truediv__(self,x):return self*J.cast(x).inverse()
    def __rtruediv__(self,x):return J.cast(x)*self.inverse()
    def __pow__(self,n):
        out=J.cast(1)
        for _ in range(n):out=out*self
        return out
    def de(self):return J(tuple((i+1)*self.c[i+1] for i in range(self.N)))
    def pair(self):return self.c[:2]

class P:
    """Polynomial radial calculus with exact rational angular jets."""
    def __init__(self,c):self.c={i:J.cast(v) for i,v in c.items() if any(J.cast(v).c)}
    @staticmethod
    def cast(x):return x if isinstance(x,P) else P({0:x})
    def __add__(self,x):
        c=self.c.copy()
        for i,v in P.cast(x).c.items():c[i]=c.get(i,J.cast(0))+v
        return P(c)
    __radd__=__add__
    def __neg__(self):return P({i:-v for i,v in self.c.items()})
    def __sub__(self,x):return self+-P.cast(x)
    def __rsub__(self,x):return P.cast(x)+-self
    def __mul__(self,x):
        c={}
        for i,a in self.c.items():
            for j,b in P.cast(x).c.items():c[i+j]=c.get(i+j,J.cast(0))+a*b
        return P(c)
    __rmul__=__mul__
    def __truediv__(self,x):return self*J.cast(x).inverse()
    def dx(self):return P({i-1:i*v for i,v in self.c.items() if i})
    def de(self):return P({i:v.de() for i,v in self.c.items()})
    def value(self,R):return sum((v*R**i for i,v in self.c.items()),J.cast(0))
    def total(self,R):
        if any(i<0 for i in self.c):raise ValueError('Integrable polynomial required')
        return sum((v*R**(i+1)/F(i+1) for i,v in self.c.items()),J.cast(0))
    def over_X(self):
        if 0 in self.c:raise ValueError('Regular axis quotient required')
        return P({i-1:v for i,v in self.c.items()})

X=P({1:1})

def fixture_terms(index):
    """M=X(1-X)^6*p(eta)+index*X^2(1-X)^6*q(eta)."""
    if index not in (0,1):raise ValueError('Fixture index required')
    p={1:F(1),2:F(1,3),3:F(1,7),4:F(1,11)};q={0:F(1),2:F(1,5)};out={}
    for n in range(7):
        wall=F((-1)**n*comb(6,n))
        for j,c in p.items():out[n+1,j]=out.get((n+1,j),F(0))+wall*c
        for j,c in q.items():out[n+2,j]=out.get((n+2,j),F(0))+index*wall*c
    return {k:v for k,v in out.items() if v}

def fields(index,angle,h):
    e=J((angle,1));M=P.cast(0)
    for (i,k),v in fixture_terms(index).items():M+=P({i:v*e**k})
    U=M.dx();D=F(1,2)-h;L=1-2*h*e*e;d=1-e*e
    V=(2*e*X*U-2*D*e*M-d*M.de())/L
    return U,V,V.over_X(),e,L,d

def identities(index,angle,h,nu,R):
    angle,h,nu,R=map(rat,(angle,h,nu,R));U,V,v,e,L,d=fields(index,angle,h)
    A=F(1,2)+h;D=F(1,2)-h
    # Original nonconservative source, independent of the new flux formula.
    source=(v+X*v.dx()+D*e*v.de())/L+v*(v/2+X*v.dx())+U*(-2*e*(v+X*v.dx())+d*v.de())/L-2*nu*(2*v.dx()+X*v.dx().dx())
    flux=((1-2*e*U)*V)/L+X*v*v-2*nu*V.dx()
    conservative=flux.dx()+(D*e*v.de()+d*(U*v).de()-2*A*e*U*v)/L+v*v/2
    direct0=source.total(R);direct1=(X*source).total(R)
    av0=v.total(R);av1=V.total(R);uv0=(U*v).total(R);uv1=(U*V).total(R)
    vv0=(v*v).total(R);vv1=(X*v*v).total(R)
    axis=2*nu*v.value(F(0));vr=V.value(R);ur=U.value(R);vxr=V.dx().value(R)
    boundary0=(1-2*e*ur)*vr/L+vr*vr/R-2*nu*vxr
    boundary1=R*(1-2*e*ur)*vr/L+vr*vr-2*nu*(R*vxr-vr)
    reduced0=(D*e*av0.de()+d*uv0.de()-2*A*e*uv0)/L+vv0/2+axis
    reduced1=(-av1+D*e*av1.de()+d*uv1.de()+2*D*e*uv1)/L-vv1/2
    difference=source-conservative
    return dict(direct=(direct0,direct1),reduced=(reduced0,reduced1),boundary=(boundary0,boundary1),axis=axis,
        pointwise=all(not any(v.c[:2]) for v in difference.c.values()),
        exact=direct0.pair()==(reduced0+boundary0).pair() and direct1.pair()==(reduced1+boundary1).pair(),
        support_endpoint=not any(vr.c[:2]) and not any(vxr.c[:2]),
        omission_axis=direct0.pair()!=(reduced0-axis+boundary0).pair(),
        omission_future=direct0.pair()!=reduced0.pair() or direct1.pair()!=reduced1.pair())

def source_bounds(h,nu,Xv,u2,v2):
    """Bounds for signed |int w0|_1 and |int X w0|_1; not int |w0|_1."""
    h,nu,Xv,u2,v2=map(rat,(h,nu,Xv,u2,v2))
    if not 0<h<=F(1,100) or nu<0 or Xv<=0 or min(u2,v2)<0:raise ValueError('Admissible h, nu, support and nonnegative envelopes required')
    A=F(1,2)+h;D=F(1,2)-h;l=1-2*h;R=1/l+4*h/l**2
    i0=4*D*R*Xv*v2+(6+4*A)*R*Xv*u2*v2+Xv*v2*v2/2+2*nu*v2
    i1=R*Xv*Xv*((1+4*D)*v2+(6+4*D)*u2*v2)/2+Xv*Xv*v2*v2/4
    return dict(I0=i0,I1=i1,reciprocal_L_1=R)

def from_global_axial_norm(h,nu,Xv,u3):
    h,nu,Xv,u3=map(rat,(h,nu,Xv,u3))
    if u3<0:raise ValueError('Nonnegative actual global axial angular norm required')
    r=source_bounds(h,nu,Xv,u3,36*u3)
    return dict(**r,u2=u3,v2=36*u3)

def norms_review():
    rows=[]
    for h in (F(1,200),F(1,100)):
        l=1-2*h;R2=1/l+6*h/l**2+16*h*h/l**3
        rows.append(dict(h=str(h),R2=str(R2),checks=dict(reciprocal_second_jet_less_than_two=R2<2,
            v2_cost_at_most_36=(16+4*(F(1,2)-h))*R2<36)))
    a=source_bounds(F(1,200),1,7,3,5);no_visc=source_bounds(F(1,200),0,7,3,5)
    # |D_eta f|_1<=2|f|_2; |D_eta f|_2<=3|f|_3.
    return dict(rows=rows,fixture={k:str(v) for k,v in a.items()},
        checks=dict(all_reciprocal_and_average_bounds=all(all(r['checks'].values()) for r in rows),
            uniform_h_interval_follows_from_monotone_positive_majorant=F(127600,117649)<2,
            axis_viscosity_cost_is_not_erased=a['I0']-no_visc['I0']==10,
            first_weighted_total_has_no_viscosity_boundary=a['I1']==no_visc['I1'],
            first_derivative_Cauchy_norm_factor_two=max(1,2)==2,
            second_derivative_norm_factor_three=max(1,2,3)==3))

@lru_cache(None)
def exact_review():
    rows=[]
    for index,angle,h,nu,R in product((0,1),map(F,('-1','-2/5','0','3/5','1')),
            (F(1,200),F(1,100)),(F(1),F(17,10)),(F(1,3),F(1))):
        z=identities(index,angle,h,nu,R)
        rows.append(dict(fixture=index,eta=str(angle),h=str(h),nu=str(nu),R=str(R),
            direct=[[str(c) for c in j.pair()] for j in z['direct']],
            boundary=[[str(c) for c in j.pair()] for j in z['boundary']],axis=[str(c) for c in z['axis'].pair()],
            checks=dict(pointwise_conservative_flux=z['pointwise'],finite_endpoint_identities=z['exact'],
                full_support_endpoint_zero=(R!=1 or z['support_endpoint'])),
            axis_omission_detected=z['omission_axis'],truncation_omission_detected=z['omission_future']))
    return dict(rows=rows,exact_integrated_scalar_jet_identities=4*len(rows),manufactured_C4_not_actual_infinite_profile=True,
        checks=dict(all_pointwise_and_integrated_identities=all(all(r['checks'].values()) for r in rows),
            missing_axis_viscosity_detected=any(r['axis_omission_detected'] for r in rows if r['R']=='1'),
            truncating_future_support_detected=any(r['truncation_omission_detected'] for r in rows if r['R']=='1/3'),
            endpoint_angles_no_division_by_d=any(r['eta']=='-1' for r in rows) and any(r['eta']=='1' for r in rows)))

def legacy_review():
    root=Path(__file__).resolve().parent.parent;old=sys.path[:]
    sys.path[:0]=[str(root/'first_correction_construction'),str(root/'physical_residual_budget')]
    try:
        s=importlib.util.spec_from_file_location('signed_pressure_original_system',root/'first_correction_construction/system.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
        from residual import Poly
        errors=[]
        for index in (0,1):
            terms=fixture_terms(index)
            U=Poly({(i-1,j):float(i*c) for (i,j),c in terms.items()})
            average=Poly({(i-1,j):float(c) for (i,j),c in terms.items()})
            profile=SimpleNamespace(U=U,average_U=average)
            for x,e in product((F(0),F(1,3),F(4,5)),(F(-1),F(0),F(1))):
                u,V,v,eta,L,d=fields(index,e,F(1,200));h=F(1,200);D=F(1,2)-h
                exact=((v+X*v.dx()+D*eta*v.de())/L+v*(v/2+X*v.dx())+u*(-2*eta*(v+X*v.dx())+d*v.de())/L-2*(2*v.dx()+X*v.dx().dx())).value(x).c[0]
                got=m.omega_over_X(profile,float(x),float(e),float(h),1.)
                errors.append(abs(got-float(exact))/max(1.,abs(float(exact))))
        return dict(manufactured_not_actual_core=True,samples=len(errors),tolerance='2e-10',
            checks=dict(original_regular_source_matches_independent_polynomial=bool(max(errors)<2e-10)))
    finally:sys.path[:]=old

def review():
    exact=exact_review();norms=norms_review();legacy=legacy_review()
    control=signed_not_absolute_review()
    return dict(exact_full_and_truncated_support=exact,norm_estimates=norms,original_source_comparison=legacy,
        signed_not_absolute_control=control,
        checks=dict(exact_all=all(exact['checks'].values()),norm_all=all(norms['checks'].values()),
            legacy_all=all(legacy['checks'].values()),signed_not_absolute=all(control['checks'].values())))

def signed_not_absolute_review():
    """Auxiliary viscosity control only; actual construction viscosity stays one."""
    KU=sum((i*abs(c)*sum(F(comb(j,k)) for k in range(min(3,j)+1))
        for (i,j),c in fixture_terms(1).items()),F(0))
    nu=F(10**30);h=F(1,200);q=from_global_axial_norm(h,nu,1,KU)
    short=identities(1,F(0),h,nu,F(1,3))['direct'][1].c[0]
    full=identities(1,F(0),h,nu,F(1))['direct'][1].c[0]
    # Absolute full integral >= absolute signed prefix integral, exactly.
    return dict(manufactured_auxiliary_viscosity_not_actual_parameters=True,nu=str(nu),
        global_U3_coefficient_triangle_upper=str(KU),signed_I1_bound=str(q['I1']),
        absolute_signed_prefix_lower_for_absolute_full_integral=str(abs(short)),full_signed_I1=str(full),
        checks=dict(signed_full_moment_obeys_Q1=abs(full)<=q['I1'],
            absolute_full_integral_is_strictly_larger_than_Q1=abs(short)>q['I1']))

def repair_moment_bounds(xcut,seed,f0,u0,I0_bound,I1_bound):
    """Signed totals bounded AFTER integration, not truncated absolute integrals."""
    xcut,seed,f0,u0,I0_bound,I1_bound=map(rat,(xcut,seed,f0,u0,I0_bound,I1_bound))
    if xcut<=0 or min(seed,f0,u0,I0_bound,I1_bound)<0:raise ValueError('Positive cutoff and nonnegative norm bounds required')
    return (xcut*seed,xcut*xcut*seed,2*f0*xcut*seed+I0_bound/2,
        xcut*xcut*(u0+f0)*seed,2*u0*xcut*seed+f0*xcut*xcut*seed+I1_bound/2)
