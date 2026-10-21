"""Exact chart jets and geometry certificates; no actual huge field sampling."""
from fractions import Fraction as F
import importlib.util
from math import factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


dynamic=load("reserved_dynamic_budget",HERE.parent/"core_dynamic_mean_recheck/proof_checks.py")
means=dynamic.means


def exact(x):
    if type(x) not in (int,F):
        raise TypeError("Exact rational input required")
    return F(x)


class Jet:
    """Normalized bivariate Taylor coefficients, truncated by total degree."""
    def __init__(self,n,terms=None):
        if type(n) is not int or n<0:
            raise ValueError("Nonnegative integer order required")
        self.n=n
        self.c={k:exact(v) for k,v in (terms or {}).items() if v and sum(k)<=n}

    def __add__(self,other):
        if not isinstance(other,Jet):
            other=Jet(self.n,{(0,0):other})
        if self.n!=other.n:
            raise ValueError("Equal truncation orders required")
        c=dict(self.c)
        for k,v in other.c.items():c[k]=c.get(k,F(0))+v
        return Jet(self.n,c)

    __radd__=__add__

    def __neg__(self):return Jet(self.n,{k:-v for k,v in self.c.items()})
    def __sub__(self,other):return self+-other if isinstance(other,Jet) else self+(-exact(other))
    def __rsub__(self,other):return -self+other

    def __mul__(self,other):
        if not isinstance(other,Jet):other=Jet(self.n,{(0,0):other})
        if self.n!=other.n:raise ValueError("Equal truncation orders required")
        c={}
        for (i,j),v in self.c.items():
            for (k,l),w in other.c.items():
                key=(i+k,j+l)
                if sum(key)<=self.n:c[key]=c.get(key,F(0))+v*w
        return Jet(self.n,c)

    __rmul__=__mul__

    def power_one(self,a):
        """(1+self)^a for a zero-constant jet."""
        a=exact(a)
        if self.c.get((0,0),0):raise ValueError("Zero constant term required")
        out,term,coefficient=Jet(self.n,{(0,0):1}),Jet(self.n,{(0,0):1}),F(1)
        for k in range(1,self.n+1):
            term=term*self;coefficient*=F(a-k+1,k);out=out+coefficient*term
        return out


def implicit_jet(h,eta,n):
    h,eta=exact(h),exact(eta)
    if not 0<h<=F(1,100) or abs(eta)>1:
        raise ValueError("Construction coordinate range required")
    t=Jet(n);p=Jet(n,{(1,0):1});q=Jet(n,{(0,1):1})
    L=1-2*h*eta**2
    for degree in range(1,n+1):
        rhs=q+(eta+p)*(eta+p)*t.power_one(2*h)-eta**2
        terms=dict(t.c)
        for i in range(degree+1):
            key=(i,degree-i);terms[key]=rhs.c.get(key,F(0))/L
        t=Jet(n,terms)
    return t


def c_norm(m,cstar):
    cstar=exact(cstar)
    if type(m) is not int or m<0 or cstar<0:
        raise ValueError("Nonnegative order and moment magnitude required")
    return cstar*(4+sum(32**i*64**(n-i) for n in range(1,m+1) for i in range(n+1)))


def geometry(l,u,Xa,Xb,t1):
    l,u,Xa,Xb,t1=map(exact,(l,u,Xa,Xb,t1))
    if not 0<l<u or not 0<Xa<l*l/2<=u*u/2<Xb or t1<=0:
        raise ValueError("Fixed hull strictly inside active annulus required")
    ra,rb=l*l/(2*Xa),2*Xb/(u*u)
    da,db=2*(ra-1)/(ra+1),2*(rb-1)/(rb+1)
    return dict(Rminus=l/2,Rplus=2*u,left_log_gap_lower=da,right_log_gap_lower=db,
                log_weight_comparison=t1*t1/(da*da)+4/(db*db))


def reserved_budget(m,eps,l,u,cstar,norms):
    """Supply proven shell/c_rho constants only; other field norms are inputs."""
    l,u=exact(l),exact(u)
    if not 0<l<u:raise ValueError("Positive fixed potential hull required")
    if set(norms)!=set(dynamic.KEYS)-{"c_rho"}:
        raise ValueError("Complete remaining dynamic norm ledger required")
    if any(norms[k]!=0 for k in ("Drr","Drt","Drz","Dtt","Dtz","Dzz")):
        raise ValueError("Reserved slow-only budget requires the fresh fixed complete wave")
    inputs=dict(norms);inputs["c_rho"]=c_norm(m+2,cstar)
    return dynamic.budget(m,eps,l/2,2*u,inputs)


def jet_review():
    rows=[];checks={}
    for h,eta in ((F(1,200),F(0)),(F(1,200),F(1)),(F(1,100),F(-1)),(F(1,2**80),F(3,5))):
        n=6;t=implicit_jet(h,eta,n)
        p=Jet(n,{(1,0):1});q=Jet(n,{(0,1):1})
        residual=t-q-(eta+p)*(eta+p)*t.power_one(2*h)+eta**2
        key=f"h{h}_eta{eta}"
        checks[key+":defining_equation"]=not residual.c
        checks[key+":Cauchy_coefficient_bound"]=all(abs(v)<=F(16**sum(k),4) for k,v in t.c.items())
        L=1-2*h*eta**2
        checks[key+":exact_first_jets"]=t.c.get((1,0),0)==2*eta/L and t.c.get((0,1),0)==1/L
        if eta:
            checks[key+":omitted_implicit_divisor_detected"]=t.c[(0,1)]!=1
        chat=F(3,7)
        c=chat*(1+t)
        checks[key+":compensator_norm_bound"]=sum(abs(v) for v in c.c.values())<=c_norm(n,chat)
        P=1+p+q+p*q+q*q
        full=c*P
        checks[key+":whole_pressure_product"]=full.c.get((1,0),0)==c.c.get((1,0),0)+chat
        if eta:
            checks[key+":moving_pressure_term_nonzero"]=full.c[(1,0)]!=chat
        # Exact normalized series powers test the contour estimate for all target factors.
        for a in (1+2*h,-F(1,2)+2*h,2*h):
            power=t.power_one(a)
            checks[key+":power"+str(a)]=all(abs(v)<=F(256,9)*16**sum(k) for k,v in power.c.items())
        rows.append(dict(h=str(h),eta=str(eta),order=n,
                         t_p=str(t.c.get((1,0),F(0))),t_q=str(t.c.get((0,1),F(0))),
                         largest_numerator_bits=max(abs(v.numerator).bit_length() for v in t.c.values()),
                         c_order6_bound=str(c_norm(n,chat))))
    checks["rational_contraction_constant"]=F(6,5)*F(1,50)*F(4,3)==F(4,125)<1
    checks["rational_disk_mapping"]=F(49,256)+F(1,125)==F(6381,32000)<F(1,4)
    return dict(checks=checks,rows=rows,scope="exact order-six diagnostics; all-order argument is the uniform complex contraction")


def radial_review():
    R=means.R;b=(R-1)**4*(2-R)**4
    hat=b/b.integral().scalar();chat=(R**2*hat).integral().scalar()/2
    checks={};rows=[]
    for root in (F(1,2),F(1),F(2)):
        rho=means.Poly({k:v/root**(k[0]+1) for k,v in hat.terms.items()})
        integral=lambda f:f.primitive().at(0,2*root)-f.primitive().at(0,root)
        mass=integral(rho).scalar();moment=integral(R**2*rho).scalar()/2
        checks[f"mass_root{root}"]=mass==1
        checks[f"moment_root{root}"]=moment==root**2*chat
        checks[f"endpoint_jets_root{root}"]=all(rho.diff(0).diff(0).at(0,e)==0 for e in (root,2*root))
        rows.append(dict(s=str(root**2),mass=str(mass),c_rho=str(moment)))
    return dict(checks=checks,chat=str(chat),rows=rows)


def geometry_review():
    l,u,Xa,Xb=F(2),F(3),F(1,2),F(16)
    g=geometry(l,u,Xa,Xb,F(1,2**20))
    checks={"positive_certified_log_gaps":g["left_log_gap_lower"]>0 and g["right_log_gap_lower"]>0}
    for root in (F(1,2),F(3,4),F(1),F(2)):
        checks[f"moving_hull_root{root}"]=g["Rminus"]<=root*l<root*u<=g["Rplus"]
        checks[f"invariant_profile_edges_root{root}"]=(root*l)**2/(2*root**2)==l*l/2 and (root*u)**2/(2*root**2)==u*u/2
    offX=g["Rplus"]**2/(2*F(1,4))
    checks["rectangular_hull_floor_countercontrol"]=offX>Xb
    n={k:F(1) for k in dynamic.KEYS if k!="c_rho"}
    for k in ("Drr","Drt","Drz","Dtt","Dtz","Dzz"):n[k]=F(0)
    budget=reserved_budget(3,F(1,7),l,u,F(3,7),n)
    checks["explicit_dynamic_budget_uses_same_input_order"]=budget["input_order"]==5 and budget["residual_order"]==2
    checks["whole_pressure_budget_positive"]=budget["Mz"]>0
    return dict(checks=checks,geometry={k:str(v) for k,v in g.items()},
                outside_rectangular_point_X=str(offX),
                conditional_budget={k:str(budget[k]) for k in ("P","Qrem","Jtheta","Jz","Mz")})


def review():
    parts={"coordinate":jet_review(),"radial_compensator":radial_review(),
           "geometry_budget":geometry_review(),"pressure_tails":means.pressure_tail_review()}
    return dict(checks={p+":"+k:v for p,r in parts.items() for k,v in r["checks"].items()},parts=parts)
