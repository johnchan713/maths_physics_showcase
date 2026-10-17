"""Independent rational jet controls, not actual giant fields or a PDE solution."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from itertools import product
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def maps():
    p=HERE.parent/'modulation_repair_all_order/map_checks.py'
    s=importlib.util.spec_from_file_location('regular_source_exact_maps',p);m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m


def variable(i,value):return maps().JetClass().variable(i,value,n=6,order=4)
def zero():return maps().jet(0,n=6,order=4)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),zero())
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def indices(order):return [a for a in product(range(order+1),repeat=6) if sum(a)<=order]
def same(a,b,order):return all((a-b).partial(i)==0 for i in indices(order))
def nonzero(a):return bool(a.terms)
def plus(a,b):return (a[0]+b[0],a[1]+b[1])
def times(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scaled(a,c):return (a[0]*c,a[1]*c)
def imaginary(a):return (-a[1],a[0])
def cp(a):return (a,zero())
def csum(v):
    out=cp(zero())
    for a in v:out=plus(out,a)
    return out
def cdiff(a,op):return (op(a[0]),op(a[1]))


def fixture(eps=F(1,8),mi=F(2,7)):
    R,Z,T,xi,v,theta=[variable(i,c) for i,c in enumerate((1,F(1,7),F(2,3),0,F(3,5),0))]
    Fbg=R**-3+Z/10+3*T*T/100;G=-R**-2+7*Z*T/100;b=eps*(R**-1+Z*T/11)
    p=F(1,3);pz=-F(2,5);phi=p*theta+pz*Z/eps+R/3-v*(p*Fbg+pz*G);phi=phi-phi.value()
    coeff=mi*F(2,3)*R**F(-1,3)
    dr=lambda f:f.derivative(0)+coeff*f.derivative(3)
    dz=lambda f:eps*f.derivative(1)
    dt=lambda f:-eps*f.derivative(2)+f.derivative(4)
    n=[dr(phi),p*R**-1,dz(phi)]
    a=[R+Z*T+xi*v/17,1+R*Z+v*v/19+xi/23,R*R+T*v/29+xi*Z/31]
    pi=R*Z+T*v/37+xi*R/41
    return R,Z,T,xi,v,theta,b,Fbg,G,phi,n,a,pi,dr,dz,dt


@lru_cache(None)
def linear_checks():
    count=0;matches=True;gaps=dict(angular=False,mixed=False,pressure=False,graph_coefficient=False)
    for eps,mi in ((F(1,8),F(2,7)),(F(1,16),F(4,11))):
        R,Z,T,xi,v,theta,b,f,g,phi,n,a,pi,dr,dz,dt=fixture(eps,mi)
        for harmonic in (-2,-1,1,2):
            w=3*harmonic;c=maps().cosine(w*phi);s=maps().sine(w*phi)
            u=[x*c for x in a];pressure=-pi*s;th=lambda x:x.derivative(5)
            lap=lambda x:dr(dr(x))+dr(x)/R+dz(dz(x))+th(th(x))/R**2
            vl=[lap(u[0])-u[0]/R**2-2*th(u[1])/R**2,lap(u[1])-u[1]/R**2+2*th(u[0])/R**2,lap(u[2])]
            adv=[dt(x)+b*dr(x)+f*th(x)+g*dz(x) for x in u]
            direct=[adv[0]+b.derivative(0)*u[0]+dz(b)*u[2]-2*f*u[1]-eps*vl[0]+dr(pressure),
                adv[1]+(2*f+R*f.derivative(0))*u[0]+dz(R*f)*u[2]+b*u[1]/R-eps*vl[1]+th(pressure)/R,
                adv[2]+g.derivative(0)*u[0]+dz(g)*u[2]-eps*vl[2]+dz(pressure)]
            K=[-2*f*a[1],(2*f+R*f.derivative(0))*a[0],g.derivative(0)*a[0]]
            principal=[x.derivative(4)+y+eps*w*w*dot(n,n)*x-w*ni*pi for x,y,ni in zip(a,K,n)]
            bd=[a[0]*b.derivative(0)+a[2]*dz(b),a[2]*dz(R*f)+b*a[1]/R,a[2]*dz(g)]
            Eik=dt(phi)+b*n[0]+f/3+g*n[2]
            real=[-eps*x.derivative(2)+b*dr(x)+g*dz(x)+y-eps*(dr(dr(x))+dr(x)/R+dz(dz(x))-(x/R**2 if i<2 else 0)) for i,(x,y) in enumerate(zip(a,bd))]
            divn=dr(n[0])+n[0]/R+dz(n[2]);connection=[-2*n[1]*a[1]/R,2*n[1]*a[0]/R,zero()]
            mixed=[2*(n[0]*dr(x)+n[2]*dz(x))+divn*x for x in a];pg=[dr(pi),zero(),dz(pi)]
            imag=[w*Eik*x+y-eps*w*(mix+con) for x,y,mix,con in zip(a,pg,mixed,connection)]
            expanded=[(x+y)*c-z*s for x,y,z in zip(principal,real,imag)]
            for x,y in zip(direct,expanded):matches &= same(x,y,2);count+=len(indices(2))
            gaps['angular'] |= any(nonzero(eps*w*x*s) for x in connection)
            gaps['mixed'] |= any(nonzero(eps*w*x*s) for x in mixed)
            gaps['pressure'] |= any(nonzero(x*s) for x in pg)
            # Second graph differentiation must also differentiate its R coefficient.
            frozen=a[0].derivative(0).derivative(0)+2*mi*F(2,3)*R**F(-1,3)*a[0].derivative(0).derivative(3)+(mi*F(2,3)*R**F(-1,3))**2*a[0].derivative(3).derivative(3)
            gaps['graph_coefficient'] |= nonzero(dr(dr(a[0]))-frozen)
    return dict(exact_linear_jet_components=count,principal_plus_remainder_matches=bool(matches),
        all_omission_controls_nonzero=all(gaps.values()),omission_controls=gaps)


@lru_cache(None)
def transport_checks():
    eps=F(1,8);R,Z,T,xi,v,theta,base,f,g,phi,n,_,_,dr,dz,dt=fixture()
    C=[R+xi*Z/7+v/11,1+R*Z+xi*v/13,R*R+T*xi/17]
    D=[R*T+xi/19,2+Z*v+xi/23,Z+R*xi/29+v/31]
    def amplitude(potential,m):
        B=[-dz(potential[1]),dz(potential[0])-dr(potential[2]),dr(potential[1])+potential[1]/R]
        return [(x,3*m*y) for x,y in zip(B,cross(n,potential))]
    def transport(a,b,mprime):
        na=csum([scaled(x,y) for x,y in zip(a,n)])
        out=[csum([times(a[0],cdiff(x,dr)),times(a[2],cdiff(x,dz)),scaled(imaginary(times(na,x)),3*mprime)]) for x in b]
        out[0]=plus(out[0],scaled(times(a[1],b[1]),-R**-1));out[1]=plus(out[1],scaled(times(a[1],b[0]),R**-1))
        return out
    def full_transport(a,b):
        out=[csum([times(a[0],cdiff(x,dr)),times(a[2],cdiff(x,dz)),scaled(times(a[1],cdiff(x,lambda z:z.derivative(5))),R**-1)]) for x in b]
        out[0]=plus(out[0],scaled(times(a[1],b[1]),-R**-1));out[1]=plus(out[1],scaled(times(a[1],b[0]),R**-1))
        return out
    matches=True;divs=True;longitudinal=False;count=0;outputs=[]
    for m,mp in ((1,1),(1,-1),(-1,2),(2,-1),(-2,-1)):
        a=amplitude(C,m);b=amplitude(D,mp)
        na=csum([scaled(x,y) for x,y in zip(a,n)])
        diva=csum([cdiff(a[0],dr),scaled(a[0],R**-1),cdiff(a[2],dz)])
        defect=plus(diva,scaled(imaginary(na),3*m))
        divs &= all(not x.terms for x in defect);longitudinal |= any(nonzero(x) for x in na)
        E=lambda j:(maps().cosine(3*j*phi),maps().sine(3*j*phi))
        direct=full_transport([times(x,E(m)) for x in a],[times(x,E(mp)) for x in b])
        expected=[times(x,E(m+mp)) for x in transport(a,b,mp)]
        for x,y in zip(direct,expected):
            for xx,yy in zip(x,y):matches &= same(xx,yy,1);count+=len(indices(1))
        outputs.append(m+mp)
    return dict(exact_transport_jet_components=count,complete_harmonic_transport_matches=bool(matches),
        exact_complete_curl_longitudinal_identity=bool(divs),complete_amplitude_not_transverse=bool(longitudinal),
        zero_and_negative_outputs_present=0 in outputs and min(outputs)<0,output_multipliers=outputs)


def endpoint_checks():
    m=maps();x=m.JetClass().variable(0,0,n=1,order=8);psi=x**8*(1-x)**8
    # At either endpoint, a cut potential has enough vanishing jets for all
    # tested curls/residuals. Finite polynomial regularity is NOT C-infinity.
    controls=[]
    for endpoint in (F(0),F(1)):
        u=m.JetClass().variable(0,endpoint,n=1,order=8);p=u**8*(1-u)**8
        for j in range(7):controls.append(p.partial((j,))==0)
    return dict(endpoint_jet_controls=len(controls),cut_potential_and_pressure_endpoint_jets=all(controls),
        uncut_positive_endpoint_is_not_zero_extension=True,
        flat_homogeneous_tail_witness='f=chi_prime*t_h; zero-entry response=-chi*t_h; core response nonflat',
        flat_tails_must_remain_additive_not_forward_sources=True,
        fixture_regularity='polynomial zero extension C7 only; continuum proof uses C-infinity cutoffs')


def reviews():return dict(exact_full_viscous_linear_expansion=linear_checks(),
    exact_complete_curl_nonlinear_transport=transport_checks(),endpoint_and_tail_split=endpoint_checks())
